module typedef
! use lateration
use LAT_mesh
use LAT_mesh_util
! use utils
implicit none

contains
!==============================================================================
subroutine read_input(model, mesh_file, pdf,out_type,out_file,output_level)
  use LAT_source 
  type(source) :: quake
  type(pdfinputs) :: pdf
  type(model_param) :: model
  ! type(surface_type) :: surf
real(pr) ::  mu, rmin, rmax, magnitude !param
integer(pin) :: na
integer(pin)  :: io_n
character(30),intent(out) :: mesh_file,out_file
! character(30) :: fname
character(3),intent(out) :: out_type
character(7) :: pdf_type
character(3) :: param_type
integer :: gauss_no,output_level
logical :: defined_area


NAMELIST / GENERAL  / mu, rmin, rmax, na, mesh_file
NAMELIST / SIZE / param_type,magnitude, defined_area
NAMELIST / PDF1   / pdf_type
NAMELIST / PDF2  / gauss_no
NAMELIST / OUTPUT  / out_type, out_file, output_level

io_n = 11
open(io_n,file='input_file',status='old',action='read')

! Values to be read in:
rewind(io_n)
read(io_n,nml=GENERAL,END=200)

rewind(io_n)
read(io_n,nml=SIZE,END=200)

rewind(io_n)
read(io_n,nml=PDF1,END=200)



rewind(io_n)
read(io_n,nml=OUTPUT,END=200)

model%mu = mu
model%rmin = rmin

model%rmax = rmax
model%na = na
model%defined_area = defined_area
model%param_type = param_type


pdf%pdf_type = pdf_type
select case(pdf%pdf_type)
  case('gauss')
    pdf%ng = gauss_no

    rewind(io_n)
    read(io_n,nml=PDF2,END=200)


  case('defined')    ! pdf has been defined in the vtk mesh file, this will be read in later.
    write(*,*) 'Using user defined PDF, this should be included in the input mesh file'
! check that user is using the full the mesh in this case
    if(model%defined_area) then
      write(*,*) 'ERROR: if using a predefined PDF the whole mesh is used'
      write(*,*) '       in input_file either set:                       '
      write(*,*) '                     (a) defined_area=.false.          '
      write(*,*) '                     (b) use a different pdf_type      '
     stop
   endif

end select

select case(model%param_type)
 case('sd')
   model%sd = magnitude
   if ( model%defined_area ) then
        write(0,*) 'ERROR: defining slip based on stress drop requires a predefined fault area'
        stop
   endif
 case('mw')
     model%mw = magnitude
case default
   write(0,*) 'ERROR: param_type in input file not defined correctly'
   stop
end select

close(io_n)
return
! error message on reading input file
200 write(*,*) 'ERROR reading input_file'

stop

end subroutine read_input

!==============================================================================
subroutine  select_mesh_bc(amesh,quake,model)
!  Select the whole mesh as the rupture area
use LAT_source
type(model_param),intent (inout) :: model
type(mesh),intent (in) :: amesh
type(source),intent (inout) :: quake
integer :: i

! assume whole mesh is included in earthquake surface
allocate(quake%QuakeElem(amesh%Ncells),quake%QuakeNodes(amesh%Nnodes))

do i = 1,amesh%Nnodes
  quake%QuakeNodes(i) = i
enddo
do i = 1,amesh%Ncells
  quake%QuakeElem(i) = i
enddo
quake%QuakeElemNo = amesh%Ncells
quake%QuakeNodesNo = amesh%Nnodes

model%actual_width = maxval(abs(amesh%pz(:)))
model%actual_area = maxval(abs(amesh%px(:)))*model%actual_width
model%target_width = model%actual_width

return
end subroutine select_mesh_bc

!###############################################################################

!==============================================================================
subroutine check_face(node,face,candidates,bc)
integer(pin),dimension(3) :: node
integer(pin),dimension(:,:) :: candidates
integer(pin), dimension(2) :: c_face
integer(pin) :: face,i
logical :: bc

if (face == 1 ) then
 c_face = (/ node(1), node(2) /)
elseif (face == 2) then
  c_face = (/ node(2), node(3) /)
else
  c_face = (/ node(3), node(1) /)
endif

 bc = .true.
 do i =1,size(candidates,dim=1)
    if (c_face(1) == candidates(i,1).and.c_face(2) == candidates(i,2)) then
       bc = .false.
    elseif  (c_face(2) == candidates(i,1).and.c_face(1) == candidates(i,2)) then
      bc = .false.
    endif
enddo

return
end subroutine check_face
!==============================================================================
subroutine intersect(A,B,v)
!  Find the common values to both A and B
!
!
integer(pin),dimension(:),intent(in) :: A, B
integer(pin), allocatable,dimension(:),intent(out) :: v
integer(pin), allocatable,dimension(:) :: dummy
integer(pin),allocatable,dimension(:) :: haveit

integer :: i, inx

allocate(dummy(size(A)))
inx = 0
do  i = 1,size(A)
       if (allocated(haveit)) deallocate(haveit)

       call find(B,'==',A(i),haveit)
       if ( haveit(1) > 0 ) then
              inx = inx +1
              dummy(inx) = A(i)
       endif
enddo

allocate(v(inx))
v = dummy(1:inx)
deallocate(dummy,haveit)

return
end subroutine intersect
!==============================================================================
subroutine isamember(A,B,idx)
!  Check if ameshents of A are in B
!      idx(i) = 1 element of A is present in B
!      idx(i) = 0 element of A is not present in B
integer(pin),dimension(:),intent(in) :: A, B
integer(pin), allocatable,dimension(:),intent(out) :: idx
integer(pin),allocatable,dimension(:) :: haveit
integer(pin) :: i

allocate(idx(size(A)))
do  i = 1,size(A)
       call find(B,'==',A(i),haveit)
       if ( haveit(1) > 0 ) then
              idx(i) = 1
       else
              idx(i) = 0
       endif
enddo

return
end subroutine isamember
!==============================================================================
!==============================================================================
! locate returns an array mask based on what points are inside a set an given ameshent
!
FUNCTION locate2d(npts,ptx,pty,amesh) result(arr)
integer(pin),intent(in) ::npts
REAL(pr), INTENT(IN) :: ptx(npts),pty(npts)
REAL(pr), INTENT(IN) :: amesh(4)
LOGICAL :: arr(npts)
real :: ex_min,ex_max,ey_min,ey_max
integer :: i


ex_min = amesh(1)
ex_max = amesh(2)
ey_min = amesh(3)
ey_max = amesh(4)

arr = .false.
do i = 1,npts
       if ((ptx(i) <= ex_max.and.ptx(i) >= ex_min).and. &
            (pty(i) <= ey_max.and.pty(i) >= ey_min)) then

             arr(i) = .true.
       endif
enddo

END FUNCTION locate2d
!!*******************************************************
!*******************************************************
subroutine reorder(a,n)
!######################################################
! Author : André Herrero
! Contact : andherit@gmail.com, andre.herrero@ingv.it
! Public Domain (CC0 1.0 Universal)
!######################################################
implicit none

integer(pin) n
real(pr) a(n)
integer(pin) i
real(pr) mem
logical done

done=.false.
do while (.not.done)
       done=.true.
       do i=1,n-1
              if (a(i).lt.a(i+1)) then
                     mem=a(i)
                     a(i)=a(i+1)
                     a(i+1)=mem
                     done=.false.
              endif
       enddo
enddo

return
end subroutine reorder
!=============================================================================
!########################################################################################
!  ########################################################################################
subroutine interp_lin(x,y,xf,yf)
implicit none
real(pr),dimension(2) :: x,y
real(pr) :: c,m,xf,yf

m = (y(2)-y(1))/(x(2)-x(1))
c = y(1)-m*x(1)

yf = m*xf+c

!end function interp_lin
end subroutine interp_lin
!  ########################################################################################
subroutine locate(n,xx,x,ans)
IMPLICIT NONE
integer(pin),intent(in) ::n
REAL(pr), INTENT(IN) :: xx(n)
REAL(pr), INTENT(IN) :: x
INTEGER(pin) :: ans
INTEGER(pin) :: jl,jm,ju
LOGICAL :: ascnd

ascnd = (xx(n) >= xx(1))
jl=0
ju=n+1
do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
              jl=jm
       else
              ju=jm
       end if
end do
if (abs(x-xx(1)) < epsilon(x)) then
       !locate=1
       ans=1
else if (abs(x-xx(n)) < epsilon(x)) then
       !locate=n-1
       ans=n-1
else
       !locate=jl
       ans=jl
end if
END subroutine locate
!***********************************************************
!***********************************************************
subroutine set_seed()
implicit none
integer(pin) :: jobid,tot_rnos,iseed
integer(pin),dimension(:),allocatable :: a_seed
character (len=30) :: en_var
integer(pin) :: i,s

! initialization of the random generator
open(10, file='/dev/urandom', access='stream', form='UNFORMATTED')
read(10) i
close(10)
!call getenv("PBS_ARRAYID",en_var) ! get environmental variable "test"
call getenv("PBS_ARRAY_INDEX",en_var) ! get environmental variable "test"
if (en_var == '') then
      jobid = 1
else
      read(en_var,'(I5)') jobid
endif
call random_seed(size=iseed)
allocate(a_seed(1:iseed))
call random_seed(get=a_seed)
call system_clock(count=s)
! random seed is choosen based on clock time
if (jobid== 1) then
  a_seed = abs( mod((s*181)*((i-83)*359), 104729) )
else
  a_seed = abs( mod((jobid*181)*((i-83)*359), 104729) )
endif

call random_seed(put=a_seed)     ! generate random seed
deallocate(a_seed)

!call random_seed     ! generate random seed

end subroutine set_seed
!***********************************************************
!----------------------------------------------------------------------
end module typedef
