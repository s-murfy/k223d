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
! This version for trilat_time is different from the previous one used in trilat_distance.
! We need here all the neighboring nodes, even those with idnode < idnode of the current node.
!
! The subroutine computes a node to node array (nton). It is a node array. It is
! defined by an array pointing at the first element of a node list (liste).
! For the ith node in the array, the list represents all the nodes attached to this node.
! In each element liste of the list at the ith node, the liste%idnode represents
! the id node j attached to i, listd%distance the distance between node i and j, 
! list%cellonedge(2) the two Id cell on each side of edge ij
! warning : but fortran does not understand pointer arrays, so I use a container.
! ntoc is not tecnically a pointer array but an array of container containing a
! pointer.
! The main loop over the cells is repeated twice: one to get all the edges attached
! to the node idnode with idnode < current node idnode and one to complete the list
! with the edges attached to the node idnode with idnode > current node idnode.
! During the first loop, for each edge, the length of the edge and the neighboring
! cells are computed. During the second loop, we need only the idnode of the other
! node attached to the edge. 

  type(mesh) :: amesh
  type(containern), dimension(amesh%Nnodes) :: nton

  type(liste), pointer :: pcur,pref
  integer(pin) :: i,j
  integer(pin) :: p1,p2,ponnodarray,othernode
  logical :: foundit

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
           nullify(pcur%previous)
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
           pcur%next%previous=>pcur
           pcur=>pcur%next
        endif
! compiling the attribute of the edge
        nullify(pcur%next)
        pcur%idnode=othernode
        pcur%cellonedge(1)=i
        pcur%cellonedge(2)=0
        pcur%donedge=comp_donedge(amesh,ponnodarray,othernode)
        pcur%ref=>pcur
     enddo
  enddo
! second loop on the cells to complete the lists
   do i=1,amesh%Ncells
! for each cell, loop on the nodes forming the cell
      do j=1,3
! define the index in the cell defining the edge
         p1=j
         p2=j+1
         if (p2 > 3) p2=1
! define the node index using the large node index of the edge
         ponnodarray=max(amesh%cell(i,p1),amesh%cell(i,p2))
         othernode=amesh%cell(i,p1)+amesh%cell(i,p2)-ponnodarray
! define a pointer at the beginning of the list of the edges attached to ponnodarray
         pref=>nton(othernode)%ptr
         do while(associated(pref))
            if (pref%idnode == ponnodarray) exit
            pref=>pref%next
         enddo
         if (.not.associated(pref)) stop "compnton: missing canonical edge reference"
         if (.not.associated(nton(ponnodarray)%ptr)) then
! if the list is void (null), allocate the space in memory
            allocate(nton(ponnodarray)%ptr)
            pcur=>nton(ponnodarray)%ptr
            nullify(pcur%previous)
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
! if the edge already exists, we directly go to the next edge
            if (foundit) cycle
! if the edge is not listed, add an element to the list. NB At this point, pcur
! points to the last element of the list.
            allocate(pcur%next)
            pcur%next%previous=>pcur
            pcur=>pcur%next
         endif
! reverse-direction entries are neighbor-only and point to the canonical edge record
         nullify(pcur%next)
         pcur%idnode=othernode
         pcur%ref=>pref
      enddo
   enddo
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
! this function get_donedge is the generalization of the function in trilat_distance.
! this version not only retrieves the distance between two nodes but also retrieves,
! if asked through the optional arguments, the two cells on each side of the edge
! through a dimension(2) array.
! If the edge does not exist, the function stops the program.
!###############################################################################  
function get_donedge(amesh,nton,p1,p2,cellonedge)

  type(mesh) :: amesh
  type(containern), dimension(amesh%Nnodes) :: nton
  integer(pin), intent(in) :: p1,p2
  integer(pin), dimension (2), optional, intent(out) :: cellonedge
  real(pr) :: get_donedge

  type(liste), pointer :: pcur
  integer(pin) :: ponnodarray,othernode

  ponnodarray=min(p1,p2)
  othernode=max(p1,p2)
  pcur=>nton(ponnodarray)%ptr
  do while(associated(pcur))
     if (pcur%idnode==othernode) then
        get_donedge=pcur%donedge
        if (present(cellonedge)) cellonedge=pcur%cellonedge
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
