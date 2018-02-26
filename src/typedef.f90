module typedef
use lateration
use utils
implicit none



contains

!==============================================================================
!=================================================
subroutine read_input(model, mesh_file, gausspar,out_type,out_file)
  type(pdfinputs) :: gausspar
  type(model_param) :: model
real ::  mu, rmin, rmax, param
integer :: na
integer  :: io_n
character(30),intent(out) :: mesh_file,out_file
character(30) :: fname
character(3),intent(out) :: out_type
character(7) :: pdf_type
character(3) :: param_type
integer :: gauss_no
logical :: defined_area

NAMELIST / GENERAL  / mu, rmin, rmax, na, mesh_file
NAMELIST / SIZE / param_type,param, defined_area
NAMELIST / PDF   / pdf_type
NAMELIST / PDF2  / gauss_no, fname
NAMELIST / OUTPUT  / out_type, out_file

io_n = 11
open(io_n,file='input_file',status='old',action='read')

! Values to be read in:
rewind(io_n)
read(io_n,nml=GENERAL,END=200)

rewind(io_n)
read(io_n,nml=SIZE,END=200)

rewind(io_n)
read(io_n,nml=PDF,END=200)

rewind(io_n)
read(io_n,nml=PDF2,END=200)

rewind(io_n)
read(io_n,nml=OUTPUT,END=200)

model%mu = mu
model%rmin = rmin
model%rmax = rmax
model%na = na
model%defined_area = defined_area

model%param_type = param_type


gausspar%pdf_type = pdf_type
select case(gausspar%pdf_type)
  case('gauss')
    gausspar%ng = gauss_no
  case('defined')    ! pdf has been defined in an ascii file
    gausspar%pdf_fname = fname
end select

select case(model%param_type)
 case('sd')
   model%sd = param
   if ( model%defined_area ) then
        write(*,*) 'Error: defining slip base on stress drop requires a predefined fault area'
        stop
   endif
 case('mw')
   model%mw = param
case default
  write(*,*) 'param_type in input file not defined correctly'
  stop
end select


return
! error message on reading input file
200 write(*,*) 'ERROR reading input_file'
stop

end subroutine read_input
!==============================================================================
subroutine read_mesh_file(amesh,mesh_file)
type(mesh) :: amesh
character(30), intent(in) :: mesh_file

! read in quadrangle mesh from file
!call read_quad_mesh(amesh,mesh_file)

! read unstructured mesh from file
call read_tri_mesh(amesh,mesh_file)

! area is in S.I.
call calc_area(amesh)

! calculate distances
call allvsall2d(amesh,amesh%dist)

! identify edge connection between elements
call mesh_EToE(amesh)


end subroutine read_mesh_file
!==============================================================================
subroutine read_pdf_file(pdf_fname,pdf,amesh)
character(30),intent(in) :: pdf_fname
! main structure
  type(mesh) :: amesh
! pdf on cells
real, dimension(amesh%QuakeElemNo) :: pdf
integer :: i

! the pdf is defined for each cell
! the order of the cells is assumed to be the same as the mesh file
open(13,file=trim(adjustl(pdf_fname)),form = 'formatted', action = 'read')
do i =1,amesh%QuakeElemNo
    read(13,*) pdf(i)
enddo
close(13)

return
end subroutine read_pdf_file
!==============================================================================
subroutine read_quad_mesh(amesh,mesh_file)
!!  Reads in a quadrangle mesh and converts it to a triangular
!!
type(mesh) :: amesh
character(30) :: mesh_file
integer :: i,id,quad_cells
character(len=60) :: str
integer :: inx
real :: num
integer :: id_o,a,b,c,d

open(13,file=trim(adjustl(mesh_file)),form = 'formatted', action = 'read')
! skip header
do i = 1,9
       read(13,*)
enddo
! find end of nodal section & element section
i = 0
id = 0
do while ( id == 0 )
    read(13,*) str
    if (trim(str) == '**') then
        id = 1
    endif
    i = i+1
enddo
amesh%Nnodes = i-1
allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))

do i = 1,2
       read(13,*) str
enddo
id = 0
i = 0
do while ( id == 0 )
    read(13,*) str
    if (trim(str) == '**') then
        id = 1
    endif
    i = i+1
enddo
quad_cells = i-1
amesh%Ncells = quad_cells*2
allocate(amesh%cell(amesh%Ncells,3))

!====== now read data in
rewind(13)
! skip header
do i = 1,9
       read(13,*)
enddo

do i = 1,amesh%Nnodes
       read(13,*) id,amesh%px(i),amesh%py(i),amesh%pz(i)
       if (i == 1) id_o = id
enddo

! skip header
do i = 1,3
       read(13,*) str
enddo
! read element to vertices
inx = 1
do i = 1,quad_cells
       read(13,*) id,a,b,c,d
       a = a-id_o+1
       b = b-id_o+1
       c = c-id_o+1
       d = d-id_o+1
       call random_number(num)     ! randomly choose split of cell
         if (num <0.5) then
             amesh%cell(inx,:) = (/ a, b, d/)
             amesh%cell(inx+1,:) = (/ b, c, d/)
         else
             amesh%cell(inx,:) = (/ a, b, c/)
             amesh%cell(inx+1,:) = (/ c, d, a/)
         endif
         inx = inx+2
enddo
close(13)

return
end subroutine read_quad_mesh
!==============================================================================
subroutine read_tri_mesh(amesh,mesh_file)
!  Reads in an unstructured mesh that has .inp format
!
!
type(mesh) :: amesh
integer :: i,id
character(len=60) :: str
real :: x,y,z
character(30),intent(in) :: mesh_file
integer :: a,b,c
integer,allocatable,dimension(:) :: counter,idx


open(13,file=trim(adjustl(mesh_file)),form = 'formatted', action = 'read')

! skip header
do i = 1,9
       read(13,*)
enddo
! find end of nodal section & element section
i = 0
id = 0
do while ( id == 0 )
    read(13,*) str
    if (trim(str) == '**') then
        id = 1
    endif
    i = i+1
enddo
amesh%Nnodes = i-1
print*,'No. of Nodes ....',amesh%Nnodes
allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))
amesh%px(amesh%Nnodes) = 0.
amesh%py(amesh%Nnodes) = 0.
amesh%pz(amesh%Nnodes) = 0.


do i = 1,2
       read(13,*) str
enddo
id = 0
i = 0
do while ( id == 0 )
    read(13,*) str
    if (trim(str) == '**') then
        id = 1
    endif
    i = i+1
enddo
amesh%Ncells = i-1
print*,'No of cells......',amesh%Ncells
allocate(amesh%cell(amesh%Ncells,3))


!====== now read data in =======
rewind(13)
! skip header
do i = 1,9
       read(13,*)
enddo
! read in vertices data
allocate(counter(amesh%Nnodes))
do i = 1,amesh%Nnodes
  read(13,*) id,x,y,z
  counter(i) = id
  amesh%px(i) = x
  amesh%py(i) = y
  amesh%pz(i) = z
enddo

! skip header
do i = 1,3
       read(13,*) str
enddo
! read element to vertices
do i = 1,amesh%Ncells

       read(13,*) id,a,b,c
       if (allocated(idx)) deallocate(idx)
       call  find(counter,'==',a,idx)
       if (idx(1).eq.-1) then
          print*,"problem reading in elements"
          stop
       else
         amesh%cell(i,1) = idx(1)
       endif
       if (allocated(idx)) deallocate(idx)
       call  find(counter,'==',b,idx)

       if (idx(1).eq.-1) then
          print*,"problem reading in elements!"
          stop
       else
         amesh%cell(i,2) = idx(1)
       endif

       if (allocated(idx)) deallocate(idx)

       call  find(counter,'==',c,idx)

       if (idx(1).eq.-1) then
          print*,"problem reading in elements!"
          stop
       else
         amesh%cell(i,3) = idx(1)
       endif
enddo
close(13)
end subroutine
!==============================================================================
subroutine calc_area(amesh)
implicit none
type(mesh) :: amesh
integer :: i
real,dimension(3) :: x,y,z,vec_ab,vec_ac
real :: ab,ac,res,theta
real :: tot_area
allocate(amesh%area(amesh%Ncells))

tot_area = 0
do i = 1,amesh%Ncells
  x = amesh%px(amesh%cell(i,:))
  y = amesh%py(amesh%cell(i,:))
  z = amesh%pz(amesh%cell(i,:))
       vec_ab = (/ x(2)-x(1), y(2)-y(1), z(2)-z(1) /)
       vec_ac = (/ x(3)-x(1), y(3)-y(1), z(3)-z(1) /)
  ac = norm2(vec_ac)
  ab = norm2(vec_ab)
  res = dot_product(vec_ab, vec_ac)/ab/ac
  theta = acos(res)
  amesh%area(i) = ab*ac*sin(theta)/float(2)
  tot_area = tot_area+amesh%area(i)
enddo
print*,'Total area (m^2):......', tot_area
end subroutine calc_area
!==============================================================================
!==============================================================================
subroutine mesh_EToE(amesh)
!  Create a matrix that defines neighbouring cells
!  Dimension :  (No of cells , 3)
!  Row relates to cells
!  Columns give the cell numbers of the 3 cells with adjoining faces to current cell
!  If a cell has less than 3 neighbours the index of the current cell is placed in the column
!
!
type(mesh),intent (inout) :: amesh
integer, dimension(3) :: node_id,neighbour_element
integer,allocatable,dimension(:) :: common_node
integer :: i,k,nofaces,bc_int


allocate(amesh%EToE(amesh%Ncells,3))
do i = 1,amesh%Ncells
       amesh%EToE(i,:) = i
enddo

bc_int = 0
! find boundary of mesh
do i = 1,amesh%Ncells
     node_id = amesh%cell(i,:)
     nofaces = 0
     neighbour_element = 0
! check node id against rest of array
     do k = 1,amesh%Ncells

       if (i /= k) then
           if (any(node_id(1).eq.amesh%cell(k,:)).or.      &
               any(node_id(2).eq.amesh%cell(k,:)).or.      &
               any(node_id(3).eq.amesh%cell(k,:))) then
! compare all the nodes of the two elements
              call intersect(node_id,amesh%cell(k,:),common_node)
              if (size(common_node)==2) then
          ! two common nodes = face
                    nofaces = nofaces+1
                    neighbour_element(nofaces) = k
              endif
            endif
       endif

     enddo

     if (nofaces > 0) then
       do k = 1,nofaces
          amesh%EToE(i,k) = neighbour_element(k)
       enddo
     endif
! 3 faces = internal element
! 2 faces = boundary element
! 1 face = edge element
enddo

end subroutine mesh_EToE
!==============================================================================
!==============================================================================
subroutine  select_mesh_bc(amesh)
!  Select the whole mesh as the rupture area
!
type(mesh),intent (inout) :: amesh
integer :: i

! assume whole mesh is included in earthquake surface
allocate(amesh%QuakeElem(amesh%Ncells),amesh%QuakeNodes(amesh%Nnodes))
do i = 1,amesh%Nnodes
  amesh%QuakeNodes(i) = i
enddo
do i = 1,amesh%Ncells
  amesh%QuakeElem(i) = i
enddo
amesh%QuakeElemNo = amesh%Ncells
amesh%QuakeNodesNo = amesh%Nnodes

end subroutine select_mesh_bc
!==============================================================================
!==============================================================================
subroutine  select_fault_zone(amesh,target_area)
!  Subroutine to pick rupture area
!  Cells are added iteratively until area exceeds target_area value
!  Choice of starting point is randomly choosen but it can be set by
!  prescribing a value to nuc_id
!
type(mesh),intent (inout) :: amesh
integer,allocatable,dimension(:) :: quake_nodes,quake_ameshents,irow,icol,haveit
real, intent(in) :: target_area
integer :: choice_amesh(3)
real :: nuc_x,nuc_y,nuc_z
real :: quake_area
real :: random
integer :: inx,get_out,nuc_id,node_idx,nxt_pt,ele_inx
integer :: kk,i,k


allocate(quake_nodes(amesh%Nnodes),quake_ameshents(amesh%Ncells))


quake_nodes = 0
quake_ameshents = 0
quake_area = 0.   ! initial area for earthquake zone on fault


inx = 1
node_idx = 1
ele_inx = 1
get_out = 0
do while (get_out < 1)
       if (inx == 1) then
! randomly pick starting cell
              call random_number(random)
              nuc_id = ceiling(random*amesh%Nnodes)   !
!            nuc_id = 1071  !centre of Scotia fault
              nuc_x = amesh%px(nuc_id)
              nuc_y = amesh%py(nuc_id);
              nuc_z = amesh%pz(nuc_id);
              nxt_pt = nuc_id
       else
! take next nodal point
              node_idx = node_idx+1
              if (node_idx > size(quake_nodes)) then
                  write(*,*) '***WARNING: requested Mag. poss. too big ***'
                  get_out = 1
                  cycle
              endif
              nxt_pt = quake_nodes(node_idx)
       endif
! find element connected to node
       call find(amesh%cell,'==',nxt_pt,irow,icol)       ! find the index at where amesh%dist = 0
       do kk = 1,size(irow)
              choice_amesh(:) = amesh%EToE(irow(kk),:)   !take element with a joining face
              do i = 1,3
                     if (inx == 1) then
                            quake_nodes(inx) = amesh%cell(choice_amesh(i),1);
                            quake_nodes(inx+1) = amesh%cell(choice_amesh(i),2);
                            quake_nodes(inx+2) = amesh%cell(choice_amesh(i),3);
                            inx = inx+3;
                     else               ! check to see if node is already in list
                         do k = 1,3
                            call find(quake_nodes,'==',amesh%cell(choice_amesh(i),k),haveit)
                            if (haveit(1) == -1) then      ! don't have node
                                   quake_nodes(inx) =  amesh%cell(choice_amesh(i),k);
                                   inx = inx+1;
                            endif
                         enddo
                     endif
! check to see if element is already in quake_ameshents
                    call find(quake_ameshents,'==',choice_amesh(i),haveit)
                    if (haveit(1) == -1) then      ! don't have element
                           quake_ameshents(ele_inx) = choice_amesh(i);
                           quake_area = quake_area + amesh%area(quake_ameshents(ele_inx)) ! assume area is in m^2
                           if (quake_area >= target_area) then
                                    get_out = 1;
                           endif
                           ele_inx = ele_inx+1;
                     endif

              enddo

       enddo

enddo

! Save section of fault used for earthquake
allocate(amesh%QuakeElem(ele_inx-1),amesh%QuakeNodes(inx-1))
amesh%QuakeElem = quake_ameshents(1:ele_inx-1)
amesh%QuakeNodes = quake_nodes(1:inx-1)
amesh%QuakeElemNo = ele_inx-1
amesh%QuakeNodesNo = inx-1
deallocate(quake_nodes,quake_ameshents)


end subroutine select_fault_zone
!==============================================================================
subroutine find_boundary(amesh)
type(mesh),intent (inout) :: amesh

integer,allocatable,dimension(:) :: haveit
integer :: ameshs(3)
integer :: bc_int,int
integer :: inx, kk,i
integer, allocatable,dimension(:) :: id,idx
integer :: i_nodes(3)
integer :: n,nlength
integer,allocatable,dimension(:,:) :: j_nodes
integer,allocatable,dimension(:) :: int_amesh,bc_amesh,bc_nodes,common_node,id_keep



inx = 0
bc_int = 0
int = 0

allocate(int_amesh(amesh%QuakeElemNo),bc_amesh(amesh%QuakeElemNo),bc_nodes(amesh%QuakeNodesNo))
bc_nodes = 0

do i = 1,amesh%QuakeElemNo
  ! ameshs = neighbouring elements to QuakeElem(i)
       ameshs = amesh%EToE(amesh%QuakeElem(i),:)

       call isamember(ameshs,amesh%QuakeElem,id) !check if adjoining elements are in defined fault
       !id = abs(id-1)
       if (id(1) == 0) then
          ameshs(1) = amesh%QuakeElem(i) ! elements is outside defined fault
       elseif (id(2) == 0) then
          ameshs(2) = amesh%QuakeElem(i) ! elements is outside defined fault
       elseif (id(3) == 0) then
          ameshs(3) = amesh%QuakeElem(i) ! elements is outside defined fault
       endif

       ! if ameshs contains an element = QuakeElem(i) it is a boundary point
       call find(ameshs,'==',amesh%QuakeElem(i),haveit)

       if (haveit(1) == -1) then  ! is interior elements
              int = int+1
              int_amesh(int) = amesh%QuakeElem(i)
       else      ! is a boundary element
              bc_int = bc_int+1
              bc_amesh(bc_int) = amesh%QuakeElem(i)

!===== check faces for which nodes make up boundary======
! global nodal id's of current element
              i_nodes = amesh%cell(amesh%QuakeElem(i),:)
              if (allocated(idx)) deallocate(idx)
! find neighboring elements
              call find(ameshs,'\=',amesh%QuakeElem(i),idx)
              allocate(j_nodes(size(idx),3))
! global nodal id's belonging to neightbouring elements within the fault
              j_nodes = amesh%cell(ameshs(idx),:)
              n  = size(idx) ! number of neighbouring elements that are in slipping zone
              if (n == 2) then ! two faces are boundaries, all 3 node pts requried
                     do kk = 1,3
                            inx = inx+1
                            bc_nodes(inx) = i_nodes(kk)
                     enddo
              else  ! one face is a boundary,  2 node pt required
! find the common id node between the two neighbouring elements which are within the fault
                     call intersect(j_nodes(1,:),j_nodes(2,:),common_node)
! keep other two nodes
                     call find(i_nodes,'/=',common_node(1),id_keep)
                     do kk = 1,2
                            inx = inx+1
                            bc_nodes(inx) = i_nodes(id_keep(kk))
                     enddo
              endif
              deallocate(j_nodes)
       endif
enddo

! remove duplicates
call erasedupli(bc_nodes,amesh%QuakeNodesNo,nlength)

amesh%QuakeBorder_NodesNo = nlength-1       ! number of vertices that make up boundary of fault
allocate(amesh%QuakeBorder_Nodes(amesh%QuakeBorder_NodesNo))
amesh%QuakeBorder_Nodes = bc_nodes(1:amesh%QuakeBorder_NodesNo)    ! vercties that make up boundary of fault
amesh%QuakeBorder_elemNo = bc_int     ! number of elements that make up boundary
allocate(amesh%QuakeBorder_elem(amesh%QuakeBorder_elemNo))
amesh%QuakeBorder_elem = bc_amesh(1:bc_int)

end subroutine find_boundary
!==============================================================================



!==============================================================================
subroutine intersect(A,B,v)
!  Find the common values to both A and B
!
!
integer,dimension(:),intent(in) :: A, B
integer, allocatable,dimension(:),intent(out) :: v
integer, allocatable,dimension(:) :: dummy
integer,allocatable,dimension(:) :: haveit

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
!==============================================================================
subroutine isamember(A,B,idx)
!  Check if ameshents of A are in B
!      idx(i) = 1 element of A is present in B
!      idx(i) = 0 element of A is not present in B
integer,dimension(:),intent(in) :: A, B
integer, allocatable,dimension(:),intent(out) :: idx
integer,allocatable,dimension(:) :: haveit
integer :: i

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
integer,intent(in) ::npts
REAL, INTENT(IN) :: ptx(npts),pty(npts)
REAL, INTENT(IN) :: amesh(4)
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

integer n
real a(n)
integer i
real mem
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
!function interp_lin(x,y,xf) result(yf)
implicit none
real,dimension(2) :: x,y
real :: c,m,xf,yf

m = (y(2)-y(1))/(x(2)-x(1))
c = y(1)-m*x(1)

yf = m*xf+c

!end function interp_lin
end subroutine interp_lin
!  ########################################################################################
subroutine locate(n,xx,x,ans)
IMPLICIT NONE
integer,intent(in) ::n
REAL, INTENT(IN) :: xx(n)
REAL, INTENT(IN) :: x
INTEGER :: ans
INTEGER :: jl,jm,ju
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
!----------------------------------------------------------------------
end module typedef
