program k223d

! To Do : merge new structure for time -> fix surface rupture problem    
!         output rupture time and make it compatible with choice of rupture on the fault   
!         make pdf for slip location part of input vtk file         
!
!

! 1/5/2019
! This code is released on the following Public Domain license :
! The unlicense.
!
! Version 1.1
! This version of the code corresponds to the paper submitted that will be submitted in 2019
! This version has been modified in order to account for surface rupture
! It is written in Fortran90.
!
!  This code builds on k223d which produced stochastic slip distributions with
!   k^-2 power spectra on non-planar faults. In this version it is possible to
!   define the free surface and allow large slip to occur at it.
!
!   Distance is calculated using a trilateration technique on an using an
!   unstructured mesh, the programme for this can be found at: https://github.com/andherit
!   a slightly altered version of this has been used in file lateration.f90
!
!   The program can be compiled using the ./install.sh script. The executable needs to
!   be in the same directory as the input file (called 'input_file') and the mesh file
!   which needs to be in an abacus format (i.e. '*.inp' format). Planar meshes can be
!   generated using the programme make_mesh.x.
!
!    The program can be used in single use producing either a vtk file or a GMT
!      file of the slip distribution. The program can also be implimented in a
!      a job array when used on  high performance computing cluster. In this case
!      the slip distributions are outputed without the grid coordinates with the
!      number at the end of the file indicating the job array number. In this mode
!      3 different slip distributions are generated:
!                    (a) total slip ,
!                    (b) non reflected slip contribution (i.e. called primary)
!                    (c) reflected slip contribution (i.e. called secondary)
!      This is in order to reproduce the results presented in the GJI paper
!      In the "Job Array" output one vtk file is generated using the case where
!      the jobarray is 1 this is done as an example.
!
!      The switch between single and job array use is whether the program finds
!      the environmental variable "PBS_ARRAY_INDEX", if this is not blank the
!      program automatically operates in job array mode.
!
!      For questions or problems : shane.muprhy at ifremer.fr
!
use generic
use LAT_mesh_util
use LAT_mesh
use LAT_source
use distance 
use LAT_time 
! use utils
use typedef
!use lateration
use makepdf

implicit none
  type(mesh) :: amesh
  type(source) :: quake
  type(surface_type) :: surf
  type(reflect) :: surface
  type(model_param) :: model
  type(pdfinputs) :: pdf ! gaussian parameters
  real(pr), parameter :: kilo = 1000.0000
  real(pr) :: a,b       ! Strasser consts, est. of fault area
  integer(pin) :: i
  character(30) :: mesh_file,out_file,surf_file
  character(3) :: out_type
  real(pr), dimension(:), allocatable :: quake_pdf
  character(20) :: en_var
  integer(pin) :: output_level
!  integer(pin) :: nuc_id
  real(pr), dimension(:), allocatable :: dist2

! set random number
  call set_seed()

! ====================  Read input file  ====================
!call getenv("PBS_ARRAYID",en_var)
call getenv("PBS_ARRAY_INDEX",en_var)
if (en_var == '') then
    write(*,*) "In single use mode"
else
    write(*,*) "In Job Array mode"
    write(*,*) "INDEX NUMBER  ",en_var
endif

! call read_input(model, mesh_file,pdf,out_type,out_file,surf,output_level) !TO DO=> simpilfy, put as much into vtk as possible!
call read_input(model, mesh_file,pdf,out_type,out_file,output_level) 

! read mesh file
if (output_level > 2)  write(*,*) "read mesh file  ..."
! call read_mesh_file(amesh,mesh_file,surf) ! reading nodal points and elements from file
call read_mesh_file(amesh,surface,quake%nuc_id)   


! calc mean depth of each element=> TO DO account for this in read in file or read_mesh_file, where to put 'mz'
 call calc_mean_depth(amesh)

! area is in S.I.=> TO DO account for this in read in file or read_mesh_file, where to put 'area'
 call calc_area(amesh)


! read in pdf file if exists =>  !TO DO : remove this and put it input vtk file 
! if (pdf%pdf_type == 'defined') then
!    call read_pdf_file(pdf,amesh)   
! endif

! Define rupture area
if (model%defined_area) then  ! Need to define Earthquake Size on mesh
! relationship based on Strasser et al. 2010  (interface events)
  a = 4.441
  b = 0.846
  model%target_area = 10**((model%mw-a)/b)*(kilo**2.)
  model%moment=10**(1.5*model%mw+9.1)
  a = -0.882
  b = 0.351
  model%target_width = 10**(a+b*model%mw)*kilo
! calculate distance between all points on mesh
  call allvsall2d(amesh,dist)
  ! call select_fault_zone_global(amesh,quake,dist,model,surf%z0,pdf)
   call select_fault_zone_global(amesh,quake,model,pdf) ! dist is sent through LAT_source 
else ! Whole mesh is used for earthquake
  model%target_area  = 0.
  do i = 1,amesh%Ncells
      model%target_area = model%target_area + area(i)
  enddo
  a = 4.441
  b = 0.846
  model%target_area = 10**((model%mw-a)/b)*(kilo**2.)
  model%moment=10**(1.5*model%mw+9.1)
  call select_mesh_bc(amesh,quake,model)  ! select whole mesh for earthquake 
endif

! make a list of the boundary nodes around slipping area
call find_boundary(amesh,quake,surface)

! calculate distance to boundary
call multisource(amesh,quake)

! define the gaussian pdf if required
if (output_level > 2)  write(*,*) 'pdf work '
select case(pdf%pdf_type)
  case('gauss')
	    if (output_level > 2) write(*,*)'set up pdf'
      call set_pdf(pdf,quake,model%target_area) !set the pdf for the distribution of cracks on fault plane
      if (output_level > 2) write(*,*)'calc fault pdf'
      call faultpdf(pdf,quake_pdf,dist,quake,amesh)
  case('defined')    ! pdf has been defined in an ascii file
! need to define quake_pdf here based on global one
       call set_quake_pdf(pdf,quake_pdf,quake,model%defined_area)
!	     call read_pdf_file(pdf%pdf_fname,global_pdf,amesh)
  case('uniform')
       call uniform_pdf(quake_pdf,quake)
  case default
       write(*,*)'ERROR reading defintion of pdf type'
       stop
end select

! place slip on fault
if (output_level > 2)  write(*,*)'run pdftoslip'
call pdftoslip(model,amesh,quake,surface,pdf%pdf_type,quake_pdf,output_level,dist2)

! calculate rupture velocity 
call calc_rupt_front(amesh,quake)

!output
if (output_level > 2) write(*,*)'write out slip distribtion'
call write_out(amesh,quake,surface,out_type,out_file,quake_pdf,en_var,output_level,dist2)
!if (output_level > 2) then
!  call write_model_parameters(model,en_var)
!endif
if (output_level > 2)  write(*,*) 'Program finish okay'

stop
end program k223d
!###############################################################################
