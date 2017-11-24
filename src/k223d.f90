program k223d
use utils
use typedef
use lateration
use makepdf
implicit none
  type(mesh) :: amesh
  type(model_param) :: model
  type(pdfinputs) :: gausspar ! gaussian parameters

  real, parameter :: kilo = 1000.0000
  real :: a,b!       ! Strasser consts, est. of fault area
  integer :: i
  character(30) :: mesh_file,out_file
  character(3) :: out_type
  real, dimension(:), allocatable :: pdf


! ====================  Read input file  ====================
!write(*,*) "read input file  ..."
call read_input(model, mesh_file,gausspar,out_type,out_file)


! read mesh file
call read_mesh_file(amesh,mesh_file) ! reading nodal points and elements from file


!set random number seed
call random_seed

! Define rupture area
if (model%defined_area) then  ! Earthquake Size Defined
! relationship based on Strasser et al. 2010  (interface events)
  a = 4.441
  b = 0.846
  model%target_area = 10**((model%mw-a)/b)*(kilo**2.)
  model%moment=10**(1.5*model%mw+9.1)
  call select_fault_zone(amesh,model%target_area)         ! selection of fault area

else ! Whole mesh is used for earthquake
  model%target_area  = 0.
  do i = 1,amesh%Ncells
      model%target_area = model%target_area + amesh%area(i)
  enddo

  call select_mesh_bc(amesh)
endif

! make a list of the boundary nodes
call find_boundary(amesh)

! calculate distance to boundary
call multisource(amesh)


! define the gaussian pdf if required
allocate(pdf(amesh%QuakeElemNo))
select case(gausspar%pdf_type)
  case('gauss')
!  write(*,*)'set up pdf'
  call set_pdf(gausspar,amesh,model%target_area)              !set the pdf for the distribution of cracks on fault plane
!  write(*,*)'calc fault pdf'
  call faultpdf(gausspar,pdf,amesh)
  case('defined')    ! pdf has been defined in an ascii file
  call read_pdf_file(gausspar%pdf_fname,pdf,amesh)
end select

! place slip on fault
!write(*,*)'run pdftoslip'
call pdftoslip(model,amesh,gausspar%pdf_type,pdf)

!output
!write(*,*)'write out slip distribtion'
call  write_out(amesh,out_type,out_file)

write(*,*) 'Program finish okay'
stop
end program k223d
!###############################################################################
subroutine write_out(amesh,out_type,out_file)
use typedef
implicit none
type(mesh) :: amesh        ! main structure
character(3) :: out_type   !type of output file
character(30) :: out_file  ! name of output file
character(34) :: filename  ! name of output file with extension

real, dimension(:), allocatable :: slipout
integer :: i,j

! assigning a slip value to every element in whole mesh
allocate(slipout(amesh%Ncells))
slipout=0.
do i=1,amesh%QuakeElemNo
   slipout(amesh%QuakeElem(i))=amesh%slip(i)
enddo

select case(out_type)
case('gmt')
  write(filename,'(a,a)') trim(adjustl(out_file)),'.out'
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
  call dumpnodeattributevtk(11,amesh,amesh%Dist2Border,'Dist2Border',.true.)

  close(11)

!--------- ouput boundary nodes of slipping region
!  open(12,file='border.vtk',form = 'formatted')
!  call dumpQuakeBordervtk(12,amesh)
!  close(12)
!--------- output all the nodes that make up slipping region
!  open(12,file='quakenodes.vtk',form = 'formatted')
!  call dumpQuakevtk(12,amesh)
!  close(12)

case default
     write(*,*) 'ERROR in creating output file'
     stop
end select

return
end subroutine write_out
