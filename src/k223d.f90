program k223d
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
use utils
use typedef
use lateration
use makepdf
implicit none
  type(mesh) :: amesh
  type(surface_type) :: surf
  type(model_param) :: model
  type(pdfinputs) :: pdf ! gaussian parameters
  real, parameter :: kilo = 1000.0000
  real :: a,b!       ! Strasser consts, est. of fault area
  integer :: i
  character(30) :: mesh_file,out_file,surf_file
  character(3) :: out_type
  real, dimension(:), allocatable :: quake_pdf
  real, dimension(:,:), allocatable :: dist
  character(20) :: en_var
  integer :: output_level

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

call read_input(model, mesh_file,pdf,out_type,out_file,surf,output_level)

! read mesh file
if (output_level > 2)  write(*,*) "read mesh file  ..."
call read_mesh_file(amesh,mesh_file,surf) ! reading nodal points and elements from file

! read in pdf file if exists
if (pdf%pdf_type == 'defined') then
   call read_pdf_file(pdf,amesh)
endif


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
  call select_fault_zone_global(amesh,dist,model,surf%z0,pdf)
else ! Whole mesh is used for earthquake
  model%target_area  = 0.
  do i = 1,amesh%Ncells
      model%target_area = model%target_area + amesh%area(i)
  enddo
  a = 4.441
  b = 0.846
  model%target_area = 10**((model%mw-a)/b)*(kilo**2.)
  model%moment=10**(1.5*model%mw+9.1)
  call select_mesh_bc(amesh,model)
endif

! make a list of the boundary nodes around slipping area
call find_boundary(amesh,surf)

! calculate distance to boundary
call multisource(amesh)

! define the gaussian pdf if required
if (output_level > 2)  write(*,*) 'pdf work '
select case(pdf%pdf_type)
  case('gauss')
	    if (output_level > 2) write(*,*)'set up pdf'
      call set_pdf(pdf,amesh,model%target_area) !set the pdf for the distribution of cracks on fault plane
      if (output_level > 2) write(*,*)'calc fault pdf'
      call faultpdf(pdf,quake_pdf,amesh)
  case('defined')    ! pdf has been defined in an ascii file
! need to define quake_pdf here based on global one
       call set_quake_pdf(pdf,quake_pdf,amesh,model%defined_area)
!	     call read_pdf_file(pdf%pdf_fname,global_pdf,amesh)
  case('uniform')
       call uniform_pdf(quake_pdf,amesh)
  case default
       write(*,*)'ERROR reading defintion of pdf type'
       stop
end select

! place slip on fault
if (output_level > 2)  write(*,*)'run pdftoslip'
call pdftoslip(model,amesh,pdf%pdf_type,quake_pdf,output_level)
!output
if (output_level > 2) write(*,*)'write out slip distribtion'
call write_out(amesh,out_type,out_file,quake_pdf,en_var,output_level)
!if (output_level > 2) then
!  call write_model_parameters(model,en_var)
!endif
if (output_level > 2)  write(*,*) 'Program finish okay'

stop
end program k223d
!###############################################################################
!###############################################################################
subroutine multisource(amesh)
use typedef

  implicit none
  type(mesh) :: amesh
  integer :: i
  logical :: fast

!initialise array
allocate(amesh%Dist2Border(amesh%Nnodes))
amesh%Dist2Border=infinity
! the boundary is set to 1000 and this is finally subtracted from the final values
! for the distance to the border. This is done in order to mimic a plane wave
do i=1,amesh%QuakeBorder_NodesNo
     amesh%Dist2Border(amesh%QuakeBorder_Nodes(i))=10000.
enddo
! calculate distance to boundary
fast = .false.
call onevsall2d(amesh,amesh%Dist2Border,fast)

amesh%Dist2Border = amesh%Dist2Border-10000.

return
end subroutine multisource
!###############################################################################
subroutine write_out(amesh,out_type,out_file,pdf,en_var,output_level)
use typedef
implicit none
type(mesh) :: amesh        ! main structure
character(3) :: out_type   !type of output file
character(20) :: en_var
character(30) :: out_file  ! name of output file
character(34) :: filename  ! name of output file with extension
character (len=20) :: fnamef
integer :: output_level
real, dimension(amesh%QuakeElemNo) :: pdf
real, dimension(:), allocatable :: slipout,pdfout

integer :: i,j

! assigning a slip value to every element in whole mesh
allocate(slipout(amesh%Ncells))
allocate(pdfout(amesh%Ncells))

slipout = 0.
pdfout = 0.
do i = 1,amesh%QuakeElemNo
  slipout(amesh%QuakeElem(i)) = amesh%slip(i)
  pdfout(amesh%QuakeElem(i)) = pdf(i)
enddo

if (en_var == ''.or.en_var=='1') then !   write(*,*) "No pbs_array
   select case(out_type)
    case('gmt')
    write(filename,'(a,a)') trim(adjustl(out_file)),'.gmt'
    open(10,file=trim(adjustl(filename)))

    do i=1,amesh%Ncells
      if (slipout(i) < 10) then  ! if slip is between 0-9
        write(10,'(a4,f9.7)') '> -Z',slipout(i)
       else      ! if slip is between 10 - 99
         write(10,'(a4,f10.7)') '> -Z',slipout(i)
       endif
       do j=1,3
          write(10,*) amesh%px(amesh%cell(i,j)),amesh%py(amesh%cell(i,j)),amesh%pz(amesh%cell(i,j))
       enddo
    enddo
    close(10)

    case('vtk')
         write(filename,'(a,a)') trim(adjustl(out_file)),'.vtk'
         open(11,file=trim(adjustl(filename)),form = 'formatted')
         call dumpmeshvtk(11,amesh)
         call dumpcellattributevtk(11,amesh,slipout,'slip',.true.)
         write (fnamef, "(a,I5.5)") "slip."//trim(en_var)

         if (output_level > 1) then
             slipout=0.
             do i = 1,amesh%QuakeElemNo
                  slipout(amesh%QuakeElem(i))= amesh%slip_prim(i)
             enddo
             call dumpcellattributevtk(11,amesh,slipout,'slip_prim',.false.)
             slipout=0.
             do i=1,amesh%QuakeElemNo
               slipout(amesh%QuakeElem(i))= amesh%slip_sec(i)
             enddo
             call dumpcellattributevtk(11,amesh,slipout,'slip_sec',.false.)
          endif
          if (output_level > 2) then
             call dumpnodeattributevtk(11,amesh,amesh%Dist2Border,'Dist2Border',.true.)
             call dumpcellattributevtk(11,amesh,pdfout,'pdf',.true.)
          endif
          close(11)
    case default
            write(*,*) 'ERROR in creating output file for single event'
            stop
    end select


  else  !  write(*,*) "INDEX NUMBER  ",en_var
    write (fnamef, "(a,I5.5)") "slip."//trim(en_var)
    open (10,file=fnamef, form="formatted",status="unknown")
    do i = 1,amesh%QuakeElemNo
        write(10,*) slipout(i)
    enddo
    close(10)

    if (output_level > 1) then
         slipout=0.
         do i=1,amesh%QuakeElemNo
            slipout(amesh%QuakeElem(i))= amesh%slip_prim(i)
         enddo
         write (fnamef, "(a,I5.5)") "slip_prim."//trim(en_var)
         open (10,file=fnamef, form="formatted",status="unknown")
         do i = 1,amesh%Ncells
            write(10,*) slipout(i)
         enddo
         close(10)
         slipout=0.
         do i=1,amesh%QuakeElemNo
               slipout(amesh%QuakeElem(i))= amesh%slip_sec(i)
         enddo
         open (10,file=fnamef, form="formatted",status="unknown")
         do i = 1,amesh%Ncells
            write(10,*) slipout(i)
         enddo
         close(10)
    endif

    if (output_level > 2) then
         write (fnamef, "(a,I5.5)") "pdf."//trim(en_var)
         open (10,file=fnamef, form="formatted",status="unknown")
         do i = 1,amesh%Ncells
             write(10,*) pdfout(i)
         enddo
         close(10)

    endif
endif
! !--------- ouput boundary nodes of slipping region
!            open(12,file='border.vtk',form = 'formatted')
!            call dumpQuakeBordervtk(12,amesh)
!            close(12)
!          !--------- output all the nodes that make up slipping region
!            open(12,file='quakenodes.vtk',form = 'formatted')
!            call dumpQuakevtk(12,amesh)
!            close(12)

return
end subroutine write_out
!###############################################################################
subroutine write_model_parameters(model,en_var)
use lateration
  implicit none
  character(20) :: en_var  ! name of output file
  character (len=20) :: fnamef

  type(model_param),intent(in) :: model

  if (en_var == '') then
     write (fnamef, "(a,I5.5)") "model_parmaters"
  else
     write (fnamef, "(a,I5.5)") "model_parmaters."//trim(en_var)
  endif
  open(11,file=trim(adjustl(fnamef)),form = 'formatted')
!  open(11,file='model_parmaters.txt',form = 'formatted')
  write(11,*) 'Target Area....',model%target_area
  write(11,*) 'Actual Area....',model%actual_area
  write(11,*) 'Target Width....',model%target_width
  write(11,*) 'Actual Width....',model%actual_width
  write(11,*) 'Length..........',model%actual_length

  close(11)

end subroutine write_model_parameters
!###############################################################################
