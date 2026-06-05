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
use typedef
use makepdf

implicit none
  type(mesh) :: amesh
  type(surface_type) :: surf
  type(reflect) :: surface
  type(model_param) :: model
  type(mesh_geometry), intent(out) :: geom

  type(pdfinputs) :: pdf ! gaussian parameters
  real(pr), parameter :: kilo = 1000.0000
  ! real(pr) :: a,b       ! Strasser consts, est. of fault area
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
call read_mesh_file(amesh,surface,pdf,quake%nuc_id)   !TO DO => add checks for different elements in mesh e.g. velocity and time 


! calculate geometrical parameters for mesh
call build_mesh_geometry(amesh, surface, geom)

! define the gaussian pdf if required
if (output_level > 2)  write(*,*) 'pdf work '
select case(pdf%pdf_type)
  case('gauss')
	    if (output_level > 2) write(*,*)'set up pdf'
      call set_pdf(pdf,amesh,model%target_area) !set the pdf for the distribution of cracks on fault plane
      if (output_level > 2) write(*,*)'calc fault pdf'
      call faultpdf(pdf,quake_pdf,dist,amesh)   ! in future : what is the difference between this and set_pdf ? 
       if (output_level > 2) write(*,*)'calc quake pdf'
  case('defined')    ! pdf has been defined in an ascii file
! need to define quake_pdf here based on global one
       call set_quake_pdf(pdf,quake_pdf,amesh)
!	     call read_pdf_file(pdf%pdf_fname,global_pdf,amesh)
  case('uniform')
       call uniform_pdf(quake_pdf,amesh)
  case default
       write(*,*)'ERROR reading defintion of pdf type'
       stop
end select

! place slip on fault
if (output_level > 2)  write(*,*)'run pdftoslip'
call pdftoslip(model,amesh,surface,pdf%pdf_type,quake_pdf,output_level,dist2)

! calculate rupture velocity 
call calc_rupt_front(amesh)  ! in reality don't need to pass quake, but it could be used in future to define rupture velocity based on slip or something else

!output
if (output_level > 2) write(*,*)'write out slip distribtion'
call write_out(amesh,surface,out_type,out_file,quake_pdf,en_var,output_level,dist2)
!if (output_level > 2) then
!  call write_model_parameters(model,en_var)
!endif
if (output_level > 2)  write(*,*) 'Program finish okay'

stop
end program k223d
!###############################################################################
