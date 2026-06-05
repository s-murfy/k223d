module lists

use generic
implicit none

! bidirective elements carrying only one integer (used for a node index)

type listn
  integer(pin) :: idnode
  type(listn), pointer :: previous
  type(listn), pointer :: next
end type listn

type containerc
  type(listn), pointer :: ptr
end type containerc

! 2024C unidirective elements carrying a node index (integer), a real value (distance)
!       and the two (or one) cellId on the side of the edge.
!       This structure is associated thus to an edge

type liste
  integer(pin) :: idnode
  real(pr) :: donedge
  integer(pin), dimension(2) :: cellonedge
  type(liste), pointer :: ref
  type(liste), pointer :: next
  type(liste), pointer :: previous
end type liste

type containern
  type(liste), pointer :: ptr
end type containern

contains
!#########################################################
subroutine wipe(ntodo)
  type(listn), pointer :: ntodo
  type(listn), pointer :: p1,p2

  if (.not.associated(ntodo)) then
     write(*,*) "warning... ntodo is unassociated yet"
     return
  endif
  if (.not.associated(ntodo%next)) then
     deallocate(ntodo)
     return
  endif
  p2=>ntodo
  do while (associated(p2%next))
     p1=>p2
     p2=>p2%next
     deallocate(p1)
  enddo
  deallocate(p2)
end subroutine wipe

!#########################################################
subroutine removepminfromthelist(pmin,ntodo)
! case1 Only one node in the list
! case2 the first node
! case3 the last node
! case4 whatever node inside the list which is not the first or the last

  type(listn), pointer :: ntodo,pmin,pn,pp
  type(listn), pointer :: first,last

  first=>ntodo
  last=>ntodo%previous
! case 1 one node in the list
!  if (pmin%idnode == first%idnode .and. .not.associated(pmin%next)) then
!  the previous test in equivalent because only the last element has a null next pointer
  if (pmin%idnode == first%idnode .and. pmin%idnode ==last%idnode) then
     if (verbose==2) write(*,*) 'remove pmin case1'
     deallocate(ntodo)
     nullify(ntodo)
     return
  endif
! case 2 node is at the beginning of the list
  if (pmin%idnode == first%idnode) then
!  if (.not.associated(pmin%previous)) then
     if (verbose==2) write(*,*) 'remove pmin case2'
     ntodo=>ntodo%next
     deallocate(pmin)
     ntodo%previous=>last
     nullify(last%next)
     return
  endif
  ! case 3 : pmin points on the last element of ntodo list
  if (pmin%idnode == last%idnode) then
     if (verbose==2) write(*,*) 'remove pmin case3'
     pp=>pmin%previous
     deallocate(pp%next)
     nullify(pp%next)
     ntodo%previous=>pp
     return
  endif
  if (pmin%idnode /= first%idnode .and. pmin%idnode /= last%idnode) then
!  if (associated(pmin%previous).and.associated(pmin%next)) then
     if (verbose==2) write(*,*) 'remove pmin generic'
     pn=>pmin%next
     pp=>pmin%previous
     deallocate(pmin)
     pp%next=>pn
     pn%previous=>pp
     return
  endif
  stop "shoudnt be there"
end subroutine removepminfromthelist
!###############################################################################
subroutine printlist(alist)
! Warning ... the output is shifted by 1 for compatibility with paraview

  type(listn), pointer :: alist
  type(listn), pointer :: pcur

  pcur=>alist
  do while (associated(pcur))
!     write(*,'(i5,$)') pcur%idnode
     write(*,'(i5,$)') pcur%idnode-1
     pcur=>pcur%next
  enddo
  write(*,*)
end subroutine printlist
!###############################################################################
subroutine updatelist(anode,alist)
! Check the presence of a node (anode) in a list (alist).
! If the node is not present, it is added to the list at its end.

  type(listn), pointer :: alist
  integer(pin) :: anode

  type(listn), pointer :: pcur,last
  logical :: isinit

! is anode in alist ?
  isinit=.false.
  pcur=>alist
  do while (associated(pcur))
     isinit=(pcur%idnode == anode)
     last=>pcur
     pcur=>pcur%next
!    if yes exit
     if (isinit) return
  enddo
!    if no add anode to alist at the end
  allocate(last%next)
  last%next%previous=>last
  nullify(last%next%next)
  last=>last%next
  last%idnode=anode
  return
end subroutine updatelist
!###############################################################################
function isinit(anode,alist)
! Check the presence of a node (anode) in a list (alist).

  type(listn), pointer :: alist
  integer(pin) :: anode
  logical :: isinit

  type(listn), pointer :: pcur

  isinit=.false.
  pcur=>alist
  do while (associated(pcur))
     isinit=(pcur%idnode == anode)
     if (isinit) return
     pcur=>pcur%next
  enddo
  return
end function isinit
!###############################################################################
end module lists
