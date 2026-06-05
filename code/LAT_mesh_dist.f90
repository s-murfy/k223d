module LAT_mesh
use generic
use lists
implicit none

type mesh
! Nnodes : no of nodes ; Ncells : number of cells 
  integer(pin) :: Nnodes,Ncells
! Px,Py,Pz : x,y,z location of nodes (Nnodes)
  real(pr), dimension(:), allocatable :: px,py,pz
! cell : relationship between cell and nodes (Ncells,3)
!        The row gives the cell and the three columns refer 
!        to node indexes
  integer(pin), dimension(:,:), allocatable :: cell 
end type mesh     

contains

!###############################################################################
subroutine compnton(amesh,nton)
! 2024C the whole routine
! compute a node to node array (nton). It is a node array. It is defined by an array
! pointing at the first element of a node list (listd).
! For the ith node in the array, the list represents all the nodes attached to this node.
! In each element listd of the list at the ith node, the listd%idnode represents
! the id node j attached to i, listd%distance the distance between node i and j, 
! list%cellonedge(2) the two Id cell on each side of edge ij
! warning : but fortran does not understand pointer arrays, so I use a container.
! ntoc is not tecnically a pointer array but an array of container containing a
! pointer.
! to avoid redundancy, the edge is always defined only once in the list of the 
! minor Idnode of the edge. Thus the last element of the array nton is a void list.

  type(mesh) :: amesh
  type(containern), dimension(amesh%Nnodes) :: nton

  type(liste), pointer :: pcur
  integer(pin) :: i,j
  integer(pin) :: p1,p2,ponnodarray,othernode
  logical :: foundit
! 2024Cdebug
!  integer(pin), dimension(10) :: hist
!  integer(pin) :: ibor

  if (verbose==1) write(*,*) 'entering in compnton'
! initialisation at null
  do i=1,amesh%Nnodes
     nullify(nton(i)%ptr)
  enddo
! loop on the cells
  do i=1,amesh%Ncells
! for each cell, loop on the nodes forming the cell
     do j=1,3
! define the index in the cell defining the edge
        p1=j
        p2=j+1
        if (p2 > 3) p2=1
! define the node index using the small node index of the edge
        ponnodarray=min(amesh%cell(i,p1),amesh%cell(i,p2))
        othernode=amesh%cell(i,p1)+amesh%cell(i,p2)-ponnodarray
! define a pointer at the beginning of the list of the edges attached to ponnodarray
        if (.not.associated(nton(ponnodarray)%ptr)) then
! if the list is void (null), allocate the space in memory
           allocate(nton(ponnodarray)%ptr)
           pcur=>nton(ponnodarray)%ptr
        else
! looking for the edge inside the list
           pcur=>nton(ponnodarray)%ptr
           foundit=.false.
           do while(.true.)
              foundit=(pcur%idnode==othernode)
              if (foundit) exit
              if (associated(pcur%next)) then
                 pcur=>pcur%next
                 cycle
              endif
              exit
           enddo
! if the edge already exists, it means that we need only add the second cell ID
! along the edge.
           if (foundit) then
              pcur%cellonedge(2)=i
              cycle
           endif
! if the edge is not listed, add an element to the list. NB At this point, pcur
! points to the last element of the list.
           allocate(pcur%next)
           pcur=>pcur%next
        endif
! compiling the attribute of the edge
        nullify(pcur%next)
        pcur%idnode=othernode
        pcur%cellonedge(1)=i
        pcur%cellonedge(2)=0
        pcur%donedge=comp_donedge(amesh,ponnodarray,othernode)
     enddo
  enddo
!  do i=1,10                            ! 2024Cdebug
!     hist(i)=0                         ! 2024Cdebug
!  enddo                                ! 2024Cdebug
!  do i=1,amesh%Nnodes
!     ibor=0
!     pcur=>nton(i)%ptr
!     do while(associated(pcur))
!        ibor=ibor+1  
!        pcur=>pcur%next
!     enddo
!     ibor=ibor+1
!     hist(ibor)=hist(ibor)+1
!  enddo
!  do i=1,10
!     write(*,*) i-1,hist(i)
!  enddo
!  write(*,*) 'sum is : ',sum(hist)
end subroutine compnton
!###############################################################################
subroutine free_nton(nton)
! Deallocates all linked-list nodes inside the nton array.
! The liste nodes are pointer-allocated (not allocatable) so Fortran does NOT
! free them automatically when nton goes out of scope — must be called explicitly.
  use lists
  type(containern), dimension(:), intent(inout) :: nton
  type(liste), pointer :: pcur, pnext
  integer(pin) :: i
  do i = 1, size(nton)
    pcur => nton(i)%ptr
    do while (associated(pcur))
      pnext => pcur%next
      deallocate(pcur)
      pcur => pnext
    enddo
    nullify(nton(i)%ptr)
  enddo
end subroutine free_nton
!###############################################################################
function get_donedge(amesh,nton,p1,p2)

  type(mesh) :: amesh
  type(containern), dimension(amesh%Nnodes) :: nton
  integer(pin) :: p1,p2
  real(pr) :: get_donedge

  type(liste), pointer :: pcur
  integer(pin) :: ponnodarray,othernode

  ponnodarray=min(p1,p2)
  othernode=max(p1,p2)
  pcur=>nton(ponnodarray)%ptr
  do while(associated(pcur))
     if (pcur%idnode==othernode) then
        get_donedge=pcur%donedge
        return
     endif
     pcur=>pcur%next
  enddo
  stop "Why am I here ?"
end function get_donedge

!###############################################################################
subroutine compntoc(amesh,ntoc)
! compute a node to cell array. ntoc is a list array. it is defined by an array
! pointing at the first element of the cell list.
! warning : but fortran does not understand pointer arrays, so I use a container.
! ntoc is not tecnically a pointer array but an array of container containing a
! pointer.
! "last" is a pointer array which, for each node, points at the end of the cell list.
! 2024B computation of the global array donedge

  type(mesh) :: amesh
  type(containerc), dimension(amesh%Nnodes) :: ntoc

  type(containerc), dimension(amesh%Nnodes) :: last
  integer(pin) :: i,j
  logical :: begin

  if (verbose==1) write(*,*) 'entering in compntoc'
! initialisation at null
  do i=1,amesh%Nnodes
     nullify(ntoc(i)%ptr)
     nullify(last(i)%ptr)
  enddo
! loop on the cell number
  do i=1,amesh%Ncells
! for each cell, loop on the nodes forming the cell
     do j=1,3
! if the cell list associated to the point amesh%cell(i,j) is not started, do it
! and point "last" to it
        if (.not.associated(ntoc(amesh%cell(i,j))%ptr)) then
           allocate(ntoc(amesh%cell(i,j))%ptr)
           nullify(ntoc(amesh%cell(i,j))%ptr%next)
           nullify(ntoc(amesh%cell(i,j))%ptr%previous)
           last(amesh%cell(i,j))%ptr=>ntoc(amesh%cell(i,j))%ptr
        else
! else add an element to the list after the last element of the list
           allocate(last(amesh%cell(i,j))%ptr%next)
           nullify(last(amesh%cell(i,j))%ptr%next%next)
           last(amesh%cell(i,j))%ptr%next%previous=>last(amesh%cell(i,j))%ptr
           last(amesh%cell(i,j))%ptr=>last(amesh%cell(i,j))%ptr%next
        endif
! add the cell i to the cell list for the node amesh%cell(i,j)
        last(amesh%cell(i,j))%ptr%idnode=i
!       if (verbose==2) write(*,*) 'last(amesh%cell(i,j))%ptr%idnode',last(amesh%cell(i,j))%ptr%idnode
     enddo
  enddo
  if (verbose==2) write(*,*) 'check on cell list on node :',amesh%Nnodes/2
  if (verbose==2) call printlist(ntoc(amesh%Nnodes/2)%ptr)
  if (verbose==2) write(*,*) 'ntoc(500)%ptr%idnode',ntoc(500)%ptr%idnode
end subroutine compntoc
!###############################################################################
subroutine deallocntoc(ntoc,n)

   integer(pin) :: n
   type(containerc), dimension(n) :: ntoc

   type(listn), pointer :: p1,p2
   integer(pin) :: i

   do i=1,n
      p2=>ntoc(i)%ptr
      if (.not.associated(p2)) then
         write(*,*) "warning... node #",i," is cell orphan"
         cycle
      endif
      do while (associated(p2%next))
         p1=>p2
         p2=>p2%next
         deallocate(p1)
      enddo
   enddo
end subroutine deallocntoc
!###############################################################################
subroutine free_ntoc(ntoc)
! Deallocates all linked-list nodes inside the ntoc array.
! The liste nodes are pointer-allocated (not allocatable) so Fortran does NOT
! free them automatically when ntoc goes out of scope — must be called explicitly.
  use lists
  type(containerc), dimension(:), intent(inout) :: ntoc

   type(listn), pointer :: p1,p2
   integer(pin) :: i

   do i=1,size(ntoc)
      p2=>ntoc(i)%ptr
      if (.not.associated(p2)) then
         write(*,*) "warning... node #",i," is cell orphan"
         cycle
      endif
      do while (associated(p2%next))
         p1=>p2
         p2=>p2%next
         deallocate(p1)
      enddo
   enddo

end subroutine free_ntoc
!###############################################################################
function comp_donedge(amesh,i,j)
! 2024B changing the name of the function. donedge is an array computed in compntoc.

  type(mesh) :: amesh
  integer(pin) :: i,j
  real(pr) :: comp_donedge

  comp_donedge=sqrt((amesh%px(i)-amesh%px(j))**2._pr+&
               (amesh%py(i)-amesh%py(j))**2._pr+&
               (amesh%pz(i)-amesh%pz(j))**2._pr)
  return
end function comp_donedge
!###############################################################################

end module LAT_mesh
