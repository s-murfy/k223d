program  make_mesh
! program to develop simple mesh that is regularly spaced but with triangular elements
!
!use utils
implicit none
type :: grid
integer :: Nnodes,K
integer :: n_stk,n_dip
integer :: QuakeElemNo,QuakeNodesNo
real :: ds,fwidth,flen
integer,dimension(14) :: mid_ele_pt
       real,allocatable,dimension(:,:) :: dist
       integer,allocatable,dimension(:,:) ::  EtoV,EToE,EToF
       integer,allocatable,dimension(:) :: QuakeElem,QuakeNodes
       integer :: QuakeBorder_NodesNo, QuakeBorder_ElemNo
       integer,allocatable,dimension(:) :: QuakeBorder_Nodes, QuakeBorder_Elem
       real,allocatable,dimension(:) :: VX,VY,VZ,area,slip
end type grid
type(grid) :: elem
NAMELIST / GENERAL  / n_x, n_z, dx, dz

real :: dx,dz
integer :: n_x, n_z

open(11,file='input_file',status='old',action='read')
read(11,nml=GENERAL,END=200)
close(11)
!n_x = 40  ! number of nodes along strike
!n_z = 25   ! number of nodes down dip
!dx  = 1000.
!dz  = 1000.
call create_mesh(n_x,n_z,dx,dz,elem)

!call output_mesh(elem)

call write_inp_file(elem)


return
! error message on reading input file
200 write(*,*) 'ERROR reading input_file'
stop

contains

!==============================================================================
subroutine create_mesh(n_x,n_z,dx,dz,elem)
implicit none
integer, intent(in) :: n_x,n_z
real, intent(in) :: dx,dz
real :: num
type(grid) :: elem

integer :: i,j,inx

elem%Nnodes = n_x*n_z
elem%K = (n_x-1)*(n_z-1)*2
allocate(elem%VX(elem%Nnodes),elem%VY(elem%Nnodes),elem%VZ(elem%Nnodes))
allocate(elem%EToV(elem%K,3))

inx = 0
do j = 1,n_z
  do i = 1,n_x
         inx = inx+1
         elem%VX(inx) = dx*float(i-1)
         elem%VY(inx) = dz*float(j-1)
         elem%VZ(inx) = 0.
    enddo
enddo

call random_seed
inx = 1
do i = 1,(n_x)*(n_z-1)
! there are two options for splitting the element
  if (mod(i,n_x) /= 0 ) then
    call random_number(num)        ! randomly split of cell
    if (num <0.5) then
      elem%EToV(inx,:)   =   (/ i, i+1, n_x+i+1/)
      elem%EToV(inx+1,:) =   (/ i, n_x+i+1, n_x+i/)
    else
        elem%EToV(inx,:)   =   (/ n_x+i, i, i+1 /)
        elem%EToV(inx+1,:) =   (/ n_x+i, i+1, i+1+n_x /)
    endif

    inx = inx+2
    endif
enddo

end subroutine create_mesh
!==============================================================================
subroutine output_mesh(elem)
implicit none
type(grid) :: elem
integer :: i

open(12,file='nodes_test.txt',form = 'formatted', action = 'write')
do i = 1,elem%Nnodes
	 write(12,*) i,elem%VX(i),elem%VY(i),elem%VZ(i)
enddo
close(12)

open(12,file='elements_test.txt',form = 'formatted', action = 'write')
do i = 1,elem%K
	 write(12,*) i,elem%EToV(i,1),elem%EToV(i,2),elem%EToV(i,3)
enddo
close(12)


end subroutine output_mesh
!==============================================================================
!  ########################################################################################
!***********************************************************
subroutine write_inp_file(elem)
implicit none
type(grid) :: elem
integer :: i

open(12,file='test.inp',form = 'formatted', action = 'write')
write(12,*) '* HEADING'
write(12,*) 'homemade'
do i = 1,6
       write(12,*)  '**'
enddo
write(12,*) '* Nodes'

do i = 1,elem%Nnodes
       write(12,*) i,elem%VX(i),elem%VY(i),elem%VZ(i)
enddo
write(12,*) '**'
write(12,*) '** '
write(12,*) '*Elements'
do i = 1,elem%K
	 write(12,*) i,elem%EToV(i,1),elem%EToV(i,2),elem%EToV(i,3)
enddo
write(12,*) '**'
close(12)
end subroutine write_inp_file
!==============================================================================
end program make_mesh
