module lateration
! work in progress 0710
!   done 071006 - initialization of ntodo inside onevsall
implicit none

  type model_param
    integer :: na
    integer :: gauss_no
    real :: mu, rmin, rmax, param
    real :: sd,mw,target_area,moment
    character(3) :: param_type
    logical :: defined_area
  end type model_param


  type mesh
         integer :: Nnodes,Ncells             ! Number of Nodes ; Number of element
         integer :: QuakeElemNo,QuakeNodesNo ! Number of Elements and nodes for an earthquake
         real,allocatable,dimension(:,:) :: dist ! distance matrix
         integer,allocatable,dimension(:,:) ::  cell
         integer,allocatable,dimension(:,:) ::  EToE  ! Elements to Elemeents
         integer,allocatable,dimension(:) :: QuakeElem,QuakeNodes ! Array of elements and nodes that make up slipping zone
         integer :: QuakeBorder_NodesNo, QuakeBorder_ElemNo ! Number of elements and nodes forming the slipping zone border
         integer,allocatable,dimension(:) :: QuakeBorder_Nodes, QuakeBorder_Elem ! Elements and nodes forming the slipping zone border
         real,allocatable,dimension(:) ::Dist2Border  ! distance matrix to boundary of fault
         real,allocatable,dimension(:) :: px,py,pz,area,slip ! Coordinates of the nodes; area of each element; slip in each element
  end type mesh

  type pdfinputs
     character(30) :: pdf_fname
     character(7) :: pdf_type
     integer :: ng
  ! pointer to a Node number of the Quake
     integer, dimension(:), allocatable :: vertexcenter
     real, dimension(:,:), allocatable :: sizeandhigh
  end type pdfinputs


  type listn
    integer :: idnode
    type(listn), pointer :: previous
    type(listn), pointer :: next
  end type listn

  type container
    type(listn), pointer :: ptr
  end type container

  type point
    real :: x,y,z
  end type point

  type path
    integer :: cellid
    type(point) :: pos
    type(path), pointer :: next
  end type path

  real, parameter :: infinity=1.e32
  logical :: verbose=.false.

contains

!###############################################################################
subroutine allvsall2d(amesh,distarray)

  type(mesh) :: amesh
  real, dimension(:,:), allocatable :: distarray

  type(container), dimension(amesh%Nnodes) :: ntoc
  type(listn), pointer :: ntodo,cellcur,pcur,last
  real :: d01,d02,d12,r1,r2,dface,dedge,dtest
  integer, dimension(2) :: otherninc
  integer :: i,j,k,toggle
  logical :: begin,hasbeenupdated

! intialisation of dist array
!     allocation
  allocate(distarray(amesh%Nnodes,amesh%Nnodes))
!     set to infinity
  distarray=infinity
!     trace to zero
  do i=1,amesh%Nnodes
     distarray(i,i)=0.
  enddo
! computing distance to neighbor nodes cell by cell
  do i=1,amesh%Ncells
!     distance 1-2
      call updclosedist(amesh,distarray,i,1,2)
!     distance 2-3
      call updclosedist(amesh,distarray,i,2,3)
!     distance 3-1
      call updclosedist(amesh,distarray,i,3,1)
  enddo
! computing the node-to-cell array ntoc
  call compntoc(amesh,ntoc)
! main loop on nodes - nodes is a starting point
  if (verbose) write(*,*) 'starting main loop on nodes'
  do k=1,amesh%Nnodes
  if (verbose) write(*,*) '##############################################################'
  if (verbose) write(*,*) 'distance from node #',k,'/',amesh%Nnodes
     call printperc(k,amesh%Nnodes)
! initializing  node todo list
     if (verbose) write(*,*) 'entering vicinode'
     call vicinode(amesh,k,ntoc(k)%ptr,ntodo)
! loop on the to do list
     if (verbose) call printlist(ntodo)
     if (verbose) write(*,*) 'entering in the ntodo list management loop'
     do while (associated(ntodo))
! searching the cells attached to the first element of the to do list
        if (verbose) write(*,*) 'state of ntodo :'
        if (verbose) call printlist(ntodo)
        if (verbose) write(*,*) 'propagating distance from node #',ntodo%idnode
        if (verbose) write(*,*) 'its own distance is :',distarray(k,ntodo%idnode)
        hasbeenupdated=.true.
        toggle=2
        do while(hasbeenupdated)
        hasbeenupdated=.false.
! switch the toggle
        toggle=3-toggle
        if (verbose) then
        if (toggle==1) write(*,*) 'seep forward on cellist'
        if (toggle==2) write(*,*) 'seep backward on cellist'
        endif
        if (toggle==1) then
           pcur=>ntoc(ntodo%idnode)%ptr
        else
           pcur=>last
        endif
        do while (associated(pcur))
!     find the two complemantory nodes in the current cell
            if (verbose) write(*,*) '    working on cell #',pcur%idnode
            call givencomp(amesh,otherninc,pcur,ntodo)
            if (verbose) write(*,*) '    complemantory nodes : ',otherninc
!     loop on the two complemantory nodes
            do i=1,2
               if (verbose) write(*,'(a4,4(i2.2))') 'code',k,ntodo%idnode,pcur%idnode,otherninc(i)
               if (verbose) write(*,*) '       working on complemantory nodes :',otherninc(i)
!      symmetry check
                if (verbose) write(*,*) '       symmetry check:'
                if (verbose) write(*,*) otherninc(i),'<',k,' ???'
                if (otherninc(i) < k) then
                   if (verbose) write(*,*) '       yes'
                   if (verbose) write(*,*) 'is distarray(k,otherninc(i)) > distarray(otherninc(i),k) ?'
                   if (verbose) write(*,*) distarray(k,otherninc(i)),' > ',distarray(otherninc(i),k)
                   if (distarray(k,otherninc(i)) > distarray(otherninc(i),k)) then
                      if (verbose) write(*,*) '       yes then take the complimentary'
                      distarray(k,otherninc(i))=distarray(otherninc(i),k)
                      call updatelist(otherninc(i),ntodo)
                   endif
                   cycle
                endif
!     edge propagation
               dedge=distarray(k,ntodo%idnode)+distarray(ntodo%idnode,otherninc(i))
               if (verbose) write(*,*) 'dedge : ',distarray(k,ntodo%idnode),'+',distarray(ntodo%idnode,otherninc(i))
               if (verbose) write(*,*) 'dedge : ',dedge
!     face propagation
               if (distarray(k,otherninc(3-i)) /= infinity) then
               if (verbose) write(*,*) 'computing dface because ',distarray(k,otherninc(3-i)),' is not infinity'
                  d12=distarray(ntodo%idnode,otherninc(3-i))
                  d01=distarray(otherninc(i),ntodo%idnode)
                  d02=distarray(otherninc(i),otherninc(3-i))
               if (verbose) write(*,*) 'd12,d01,d02 :'
               if (verbose) write(*,*) d12,d01,d02
                  r1=distarray(k,ntodo%idnode)
                  r2=distarray(k,otherninc(3-i))
               if (verbose) write(*,*) 'r1,r2 : ',r1,r2
                  dface=dcircle(d12,d01,d02,r1,r2)
               if (verbose) write(*,*) 'dface : ',dface
               else
                  dface=infinity
               endif
               dtest=min(dedge,dface)
!               dtest=dedge
               if (verbose) write(*,*) 'dtest=min(dedge,dface) : ',dtest
               if (verbose) write(*,*) 'dtest < distance between',k,' and ',otherninc(i)
               if (verbose) write(*,*) dtest,distarray(k,otherninc(i))
               if (dtest < distarray(k,otherninc(i))) then
                   if (verbose) write(*,*) 'better !'
! distance is better : if not in the list, add it
                  distarray(k,otherninc(i))=dtest
                  if (verbose) write(*,*) 'updatelist with :',otherninc(i)
                  call updatelist(otherninc(i),ntodo)
                  if (verbose) write(*,*) 'state of ntodo :'
                  if (verbose) call printlist(ntodo)
               endif
            enddo
! moving on the vicinity list
            if (toggle==1) then
               last=>pcur
               pcur=>pcur%next
            else
               pcur=>pcur%previous
            endif
        enddo
        enddo
! cancelling the cell vicinity list
! removing the first element of the to do list
        if (associated(ntodo%next)) then
           ntodo=>ntodo%next
           deallocate(ntodo%previous)
           nullify(ntodo%previous)
        else
           deallocate(ntodo)
           nullify(ntodo)
        endif
      enddo
  enddo
! deallocating ntoc
  call deallocntoc(ntoc,amesh%Nnodes)
end subroutine allvsall2d
!###############################################################################
subroutine multisource(amesh)
  implicit none
  type(mesh) :: amesh
  integer :: i
  logical :: fast

!initialise array
allocate(amesh%Dist2Border(amesh%Nnodes))
amesh%Dist2Border=infinity
do i=1,amesh%QuakeBorder_NodesNo
     amesh%Dist2Border(amesh%QuakeBorder_Nodes(i))=0.
enddo

! calculate distance to boundary
fast = .false.
call onevsall2d(amesh,amesh%Dist2Border,fast)

return
end subroutine multisource
!###############################################################################
subroutine pre_onevsall2d_offvertex(amesh,distarray,v,d)
! initialisation of onevsall2d in case of a unique source at the
! distance d(1), d(2) and d(3) from vertice v(1), v(2) and v(3) respectively of the SAME face.
! Use this example to create your own routine to initialize distarray
! for a generic multi source problem with n points and not only 3.

  type(mesh) :: amesh
  real, dimension(:), allocatable :: distarray
  integer, dimension(3) :: v
  real, dimension(3) :: d

  integer :: i
  type(listn), pointer :: pcur

!   allocation
  allocate(distarray(amesh%Nnodes))
!   initialisation to infinity
  distarray=infinity
!   initialisation of distarray and ntodo
  do i=1,3
     distarray(v(i))=d(i)
  enddo
  return
end subroutine pre_onevsall2d_offvertex
!###############################################################################
subroutine pre_onevsall2d_onvertex(amesh,k,distarray)
! initialisation of onevsall2d in case of a unique source on a vertex k

  type(mesh) :: amesh
  real, dimension(:), allocatable :: distarray
  integer :: k

!   allocation
  allocate(distarray(amesh%Nnodes))
!   initialisation to infinity
  distarray=infinity
!   initialization of distance at k to zero
  distarray(k)=0.
  return
end subroutine pre_onevsall2d_onvertex
!###############################################################################
subroutine onevsall2d(amesh,dist,fast)
! Compute the distance, in an array format, for all the vertices of the mesh
! amesh which distance is not set to infinity.
! Warning : dist should be initialized before to infinity for all the
! vertices except those representing the sources.

  type(mesh) :: amesh
  real, dimension(amesh%Nnodes) :: dist
  logical :: fast

  type(container), dimension(amesh%Nnodes) :: ntoc
  type(listn), pointer :: cellcur,pcur,last,pmin
  type(listn), pointer :: ntodo
  real :: d13,d23,d12,r1,r2,dface,dedge,dtest,db
  integer, dimension(2) :: otherninc
  integer :: i,j,toggle,idiff,waitfordiff
  logical :: begin,hasbeenupdated,reloop,firstrun,diffoccur
  integer :: nswp,mxswp
  real, dimension(amesh%Nnodes) :: distarray
  logical, dimension(amesh%Nnodes) :: checksecondary
  integer, parameter :: waitdiffthres=25

!  mxswp=0
! computing the node-to-cell array ntoc
  call compntoc(amesh,ntoc)
! initialize the secondary distance array
  distarray=dist
  checksecondary=.false.
! initialization of ntodo from non infinite values inside distarray
  begin=.true.
  do i=1,amesh%Nnodes
     if (distarray(i) > .9*infinity) cycle
! if a vertex has a null distance, it is not a secondary diffraction point
     if (distarray(i) < epsilon(distarray(i))) checksecondary(i)=.true.
     if (begin) then
        allocate(ntodo)
        nullify(ntodo%previous)
        nullify(ntodo%next)
        pcur=>ntodo
        begin = .false.
     else
        allocate(pcur%next)
        nullify(pcur%next%next)
        pcur%next%previous=>pcur
        pcur=>pcur%next
     endif
     pcur%idnode=i
  enddo
! initializing reloop to true
! The condition is changed or not at the end of the loop
! The fast first loop is mandatory with a null distance offset
  reloop=.true.
  firstrun=.true.
  diffoccur=.false.
  db=0.
  do while (reloop)
! loop on the to do list
  if (verbose) call printlist(ntodo)
  if (verbose) write(*,*) 'entering in the ntodo list management loop'
  do while (associated(ntodo))
        if (verbose) write(*,*) 'state of ntodo :'
        if (verbose) call printlist(ntodo)
! searching the cells attached to the first element of the to do list
        call lookformin(ntodo,distarray,amesh%Nnodes,pmin)
        if (verbose) write(*,*) 'propagating distance from node #',pmin%idnode
        if (verbose) write(*,*) 'its own distance is :',distarray(pmin%idnode)
        hasbeenupdated=.true.
        toggle=2
!        nswp=0
        do while (hasbeenupdated)
!        nswp=nswp+1
        hasbeenupdated=.false.
        toggle=3-toggle
        if (verbose) then
        if (toggle==1) write(*,*) 'seep forward on cellist'
        if (toggle==2) write(*,*) 'seep backward on cellist'
        endif
        if (toggle==1) then
           pcur=>ntoc(pmin%idnode)%ptr
        else
           pcur=>last
        endif
        do while (associated(pcur))
!     find the two complemantory nodes in the current cell
            if (verbose) write(*,*) '    working on cell #',pcur%idnode
            call givencomp(amesh,otherninc,pcur,pmin)
            if (verbose) write(*,*) '    complemantory nodes : ',otherninc
!     loop on the two complemantory nodes
            do i=1,2
               if (verbose) write(*,'(a4,3(i4.4))') 'code',pmin%idnode,pcur%idnode,otherninc(i)
               if (verbose) write(*,*) '       working on complemantory nodes :',otherninc(i)
!     edge propagation
               dedge=distarray(pmin%idnode)+donedge(amesh,pmin%idnode,otherninc(i))
               if (verbose) write(*,*) 'dedge : ',distarray(pmin%idnode),'+',donedge(amesh,pmin%idnode,otherninc(i))
               if (verbose) write(*,*) 'dedge : ',dedge
!     face propagation
               if (distarray(otherninc(3-i)) /= infinity) then
               if (verbose) write(*,*) 'computing dface because ',distarray(otherninc(3-i)),' is not infinity'
                  d12=donedge(amesh,pmin%idnode,otherninc(3-i))
                  d13=donedge(amesh,otherninc(i),pmin%idnode)
                  d23=donedge(amesh,otherninc(i),otherninc(3-i))
               if (verbose) write(*,*) 'd12,d13,d23 :'
               if (verbose) write(*,*) d12,d13,d23
                  r1=distarray(pmin%idnode)
                  r2=distarray(otherninc(3-i))
               if (verbose) write(*,*) 'r1,r2 : ',r1,r2
                  dface=dcircle(d12,d13,d23,r1,r2)
               if (verbose) write(*,*) 'dface : ',dface
               else
                  dface=infinity
               endif
               dtest=min(dedge,dface)
!              dtest=dedge
               if (verbose) write(*,*) 'dtest=min(dedge,dface) : ',dtest
               if (verbose) write(*,*) 'dtest < actual minimum on ',otherninc(i)
               if (verbose) write(*,*) dtest,distarray(otherninc(i))
               if (dtest < distarray(otherninc(i))) then
                  hasbeenupdated=.true.
                  if (verbose) write(*,*) 'better !'
! distance is better : if not in the list, add it
                  distarray(otherninc(i))=dtest
                  call updatelist(otherninc(i),ntodo)
                  if (verbose) write(*,*) 'updatelist with :',otherninc(i)
                  if (verbose) write(*,*) 'state of ntodo :'
                  if (verbose) call printlist(ntodo)
               endif
            enddo ! complementary loop on cell
! moving on the vicinity list
            if (toggle==1) then
               last=>pcur
               pcur=>pcur%next
            else
               pcur=>pcur%previous
            endif
        enddo ! neighbor cell to the current node pmin
        enddo ! sweep loop
!        mxswp=max(mxswp,nswp)
! diff case and no diffraction has been detected yet
        if (.not.firstrun.and..not.diffoccur) then
           if (distarray(pmin%idnode)+dist(idiff) < dist(pmin%idnode)) then
!   diffraction detected
              diffoccur=.true.
           else
!   no diffraction
              waitfordiff=waitfordiff-1
              if (waitfordiff == 0) then
!   countdown expired for diffraction. hard exit on ntodo list to pass to another
!                                      potential secondary source
                 call wipe(ntodo)
                 exit
              endif
           endif
        endif
! removing the treated element pmin of the to do list
        call removepminfromthelist(pmin,ntodo)
      enddo ! ntodo loop
! decision to reloop or not
! fast case
      if (fast) then
         dist=distarray
! first case of general exit : only one run is requested
         reloop=.false.
         cycle
      endif
      if (firstrun) then
         dist=distarray
      else
         if (diffoccur) then
! a secondary diff has occured
            do i=1,amesh%Nnodes
               dist(i)=min(distarray(i)+dist(idiff),dist(i))
            enddo
         endif
      endif
      idiff=nextdiffid(dist,amesh%Nnodes,checksecondary)
! second case of general exit : all the potential secondary source have been investigated.
      if (idiff == 0) then
         reloop=.false.
         cycle
      endif
! preparation for the next loop
      checksecondary(idiff)=.true.
      waitfordiff=waitdiffthres
      diffoccur=.false.
! At this point notodo is deallocated. reallocate it with idiff
      allocate(ntodo)
      nullify(ntodo%previous)
      nullify(ntodo%next)
      ntodo%idnode=idiff
! reinitialize distarray
      distarray=infinity
      distarray(idiff)=0.
      firstrun=.false.
  enddo ! end reloop section
! deallocating ntoc
  call deallocntoc(ntoc,amesh%Nnodes)
!  write(0,*) 'maximun sweep numner : ',mxswp
end subroutine onevsall2d
!###############################################################################
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
!###############################################################################
function nextdiffid(dist,n,checksecondary)
! find the minimum value in dist skipping the vertex checked yet and returns
! its position (id)

  integer :: nextdiffid,n
  real, dimension(n) :: dist
  logical, dimension(n) :: checksecondary

  integer :: i
  real :: vmin

  nextdiffid=0
  vmin=infinity
  do i=1,n
     if (checksecondary(i)) cycle
     if (dist(i) < vmin) then
        vmin=dist(i)
        nextdiffid=i
     endif
  enddo
  return
end function nextdiffid
!###############################################################################
subroutine lookformin(ntodo,distarray,n,pmin)

  type(listn), pointer :: pmin,ntodo
  integer :: n
  real, dimension(n) :: distarray

  type(listn), pointer :: pcur
  real :: vmin

  pcur=>ntodo
  vmin=distarray(pcur%idnode)
  pmin=>pcur
  do while (associated(pcur%next))
     pcur=>pcur%next
     if (distarray(pcur%idnode) < vmin) then
        vmin=distarray(pcur%idnode)
        pmin=>pcur
     endif
  enddo
  return
end subroutine lookformin
!###############################################################################
subroutine removepminfromthelist(pmin,ntodo)

  type(listn), pointer :: ntodo,pmin,pn,pp

  if (associated(pmin%previous).and.associated(pmin%next)) then
     if (verbose) write(*,*) 'remove pmin case1'
     pn=>pmin%next
     pp=>pmin%previous
     deallocate(pp%next)
     pp%next=>pn
     pn%previous=>pp
     return
  endif
  if (.not.associated(pmin%previous).and..not.associated(pmin%next)) then
     if (verbose) write(*,*) 'remove pmin case2'
     nullify(ntodo)
     return
  endif
  if (.not.associated(pmin%previous)) then
     if (verbose) write(*,*) 'remove pmin case3'
     ntodo=>ntodo%next
     deallocate(ntodo%previous)
     nullify(ntodo%previous)
     return
  endif
! last case : pmin points on the last element of ntodo list
     if (verbose) write(*,*) 'remove pmin case4'
  pp=>pmin%previous
  deallocate(pp%next)
  nullify(pp%next)
  return
end subroutine removepminfromthelist
!###############################################################################
subroutine timeonevsall2d(amesh,k,time,velocity)

! given a mesh (amesh) and velocities (velocity) defined on
! the mesh cells, timeonevsall2d computes the first arrival
! time from the node number k to all the nodes of the mesh
! (time)

  type(mesh) :: amesh
  real, dimension(:), allocatable :: time
  real, dimension(amesh%Ncells) :: velocity
  integer :: k

  type(container), dimension(amesh%Nnodes) :: ntoc
  type(listn), pointer :: ntodo,cellcur,pcur,last
  real :: d13,d23,d12,t1,t2,tface,tedge,ttest
  integer, dimension(2) :: otherninc
  integer :: i,j,toggle
  logical :: begin,hasbeenupdated
  integer :: nswp,mxswp

  mxswp=0
! intialisation of time array
!     allocation
  allocate(time(amesh%Nnodes))
!     set to infinity
  time=infinity
!     k nodes set to zero
  time(k)=0.
! computing the node-to-cell array ntoc
  call compntoc(amesh,ntoc)
  if (verbose) write(*,*) 'ntoc(500)%ptr%idnode',ntoc(500)%ptr%idnode
  if (verbose) write(*,*) '##############################################################'
  if (verbose) write(*,*) 'distance from node #',k,'/',amesh%Nnodes
! initializing  node todo list and the time associated to them
  if (verbose) write(*,*) 'entering vicinode'
  if (verbose) write(*,*) 'ntoc(k)%ptr%idnode',ntoc(k)%ptr%idnode
  if (verbose) call printlist(ntoc(k)%ptr)
  call startime(amesh,k,ntoc(k)%ptr,ntodo,time,velocity)
! loop on the to do list
  if (verbose) call printlist(ntodo)
  if (verbose) write(*,*) 'entering in the ntodo list management loop'
  do while (associated(ntodo))
! searching the cells attached to the first element of the to do list
        if (verbose) write(*,*) 'state of ntodo :'
        if (verbose) call printlist(ntodo)
        if (verbose) write(*,*) 'propagating time from node #',ntodo%idnode
        if (verbose) write(*,*) 'its own time is :',time(ntodo%idnode)
!    initilisation of the sweep loop
        hasbeenupdated=.true.
        toggle=2
        nswp=0
!    sweep loop up to the stabilisation (hasbeenupdated=.false.)
        do while (hasbeenupdated)
        nswp=nswp+1
        hasbeenupdated=.false.
!    toggle variable swap state between 1 and 2 at each sweep
        toggle=3-toggle
        if (verbose) then
        if (toggle==1) write(*,*) 'seep forward on cellist'
        if (toggle==2) write(*,*) 'seep backward on cellist'
        endif
!    pointer on the beginning or the end of the list of cells containing the node idnode
        if (toggle==1) then
           pcur=>ntoc(ntodo%idnode)%ptr
        else
           pcur=>last
        endif
!    cell list sweeping
        do while (associated(pcur))
!     find the two complemantory nodes in the current cell
            if (verbose) write(*,*) '    working on cell #',pcur%idnode
            call givencomp(amesh,otherninc,pcur,ntodo)
            if (verbose) write(*,*) '    complemantory nodes : ',otherninc
!     loop on the two complemantory nodes
            do i=1,2
               if (verbose) write(*,'(a4,4(i4.4))') 'code',k,ntodo%idnode,pcur%idnode,otherninc(i)
               if (verbose) write(*,*) '       working on complemantory nodes :',otherninc(i)
!     edge propagation
               tedge=time(ntodo%idnode)+tonedge(amesh,ntodo%idnode,otherninc(i),velocity(pcur%idnode))
               if (verbose) write(*,*) 'tedge : ',time(ntodo%idnode),'+',&
                   tonedge(amesh,ntodo%idnode,otherninc(i),velocity(pcur%idnode))
               if (verbose) write(*,*) 'tedge : ',tedge
!     face propagation when time on 2 nodes are available
               if (time(otherninc(3-i)) /= infinity) then
                  if (verbose) write(*,*) 'computing tface because ',time(otherninc(3-i)),&
                     ' is not infinity'
!     edge lengths in the current cell
                  d12=donedge(amesh,ntodo%idnode,otherninc(3-i))
                  d13=donedge(amesh,otherninc(i),ntodo%idnode)
                  d23=donedge(amesh,otherninc(i),otherninc(3-i))
                  if (verbose) write(*,*) 'd12,d13,d23 :',d12,d13,d23
!     time associated to the 2 nodes for the lateration computation
                  t1=time(ntodo%idnode)
                  t2=time(otherninc(3-i))
                  if (verbose) write(*,*) 't1,t2 : ',t1,t2
!     multi lateration kernel
                  tface=tcircle(d12,d13,d23,t1,t2,velocity(pcur%idnode))
                  if (verbose) write(*,*) 'tface : ',tface
               else
                  tface=infinity
               endif
!     result is the minimum between edge propagation and face propagation
               ttest=min(tedge,tface)
!              ttest=tedge
               if (verbose) write(*,*) 'ttest=min(tedge,tface) : ',ttest
               if (verbose) write(*,*) 'ttest < time between',k,' and ',otherninc(i)
               if (verbose) write(*,*) ttest,time(otherninc(i))
!     better result found - updating process
               if (ttest < time(otherninc(i))) then
                  if (verbose) write(*,*) 'better !'
                  time(otherninc(i))=ttest
                  hasbeenupdated=.true.
!     time is better : if not in the list, add it
                  if (verbose) write(*,*) 'updatelist with :',otherninc(i)
                  call updatelist(otherninc(i),ntodo)
                  if (verbose) write(*,*) 'state of ntodo :'
                  if (verbose) call printlist(ntodo)
               endif
            enddo
! moving on the vicinity list
            if (toggle==1) then
               last=>pcur
               pcur=>pcur%next
            else
               pcur=>pcur%previous
            endif
        enddo
        enddo
        mxswp=max(mxswp,nswp)
! move the pointer on the second element of the to do list
! remove the first element of the to do list
        if (associated(ntodo%next)) then
           ntodo=>ntodo%next
           deallocate(ntodo%previous)
           nullify(ntodo%previous)
        else
! no node to do anymore.
           nullify(ntodo)
        endif
      enddo
! deallocating ntoc
  call deallocntoc(ntoc,amesh%Nnodes)
  write(0,*) 'maximun sweep number : ',mxswp
end subroutine timeonevsall2d
!###############################################################################
subroutine printperc(a,b)

  integer :: a,b

!  write(*,*) a,b
  write(*,'(a1,$)') char(8)
  write(*,'(a1,$)') char(8)
  write(*,'(a1,$)') char(8)
!  write(*,*) a,b
  if (a == b) then
     write(*,'(a4)') '100%'
  else
     write(*,'(i2.2,a1,$)') int(100.*(float(a)/float(b))),'%'
  endif
  return
end subroutine printperc
!###############################################################################
subroutine printlist(alist)

  type(listn), pointer :: alist

  type(listn), pointer :: pcur

  pcur=>alist
  do while (associated(pcur))
     write(*,'(i5,$)') pcur%idnode
     pcur=>pcur%next
  enddo
  write(*,*)
end subroutine printlist
!###############################################################################
subroutine updatelist(anode,alist)
! Check the presence of a node (anode) in a list (alist).
! If the node is not present, it is added to the list at its end.

  type(listn), pointer :: alist
  integer :: anode

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
  integer :: anode
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
function donedge(amesh,i,j)

  type(mesh) :: amesh
  integer :: i,j
  real :: donedge

  donedge=sqrt((amesh%px(i)-amesh%px(j))**2+&
               (amesh%py(i)-amesh%py(j))**2+&
               (amesh%pz(i)-amesh%pz(j))**2)
  return
end function donedge
!###############################################################################
function tonedge(amesh,i,j,v)

  type(mesh) :: amesh
  integer :: i,j
  real :: tonedge,v

  tonedge=sqrt((amesh%px(i)-amesh%px(j))**2+&
               (amesh%py(i)-amesh%py(j))**2+&
               (amesh%pz(i)-amesh%pz(j))**2)/v
  return
end function tonedge
!###############################################################################
subroutine updclosedist(amesh,distarray,k,i,j)

   type(mesh) :: amesh
   real, dimension(:,:), allocatable :: distarray
   integer :: i,j,k

   if (distarray(amesh%cell(k,i),amesh%cell(k,j)) == infinity) then
      distarray(amesh%cell(k,i),amesh%cell(k,j))=&
           sqrt((amesh%px(amesh%cell(k,i))-amesh%px(amesh%cell(k,j)))**2+&
                (amesh%py(amesh%cell(k,i))-amesh%py(amesh%cell(k,j)))**2+&
                (amesh%pz(amesh%cell(k,i))-amesh%pz(amesh%cell(k,j)))**2)
! reciprocity
      distarray(amesh%cell(k,j),amesh%cell(k,i))=distarray(amesh%cell(k,i),amesh%cell(k,j))
   endif
end subroutine updclosedist
!###############################################################################
subroutine givencomp(amesh,otherninc,pc,pa)

  type(mesh) :: amesh
  integer, dimension(2) :: otherninc
  type(listn), pointer :: pc,pa

  integer :: i,k

  k=1
  do i=1,3
     if (amesh%cell(pc%idnode,i) /= pa%idnode) then
        otherninc(k)=amesh%cell(pc%idnode,i)
        k=k+1
     endif
  enddo
  return
end subroutine givencomp
!###############################################################################
function dplane(d12,d13,d23,r1,r2)

  real :: dplane,d12,d13,d23,r1,r2

  real :: x,y,xc,yc
  real :: dr,xn,yn,dn,d13n,d23n

! position of v3
  x=(d12**2-d23**2+d13**2)/(2*d12)
  y=sqrt(max(d13**2-x**2,0.))
! difference of distance at the base
  dr=abs(r1-r2)
! computation of the new base length
  dn=sqrt(d12**2.-dr**2.)
  if (r1 < r2) then
! case where v2 remains and V1 is changed
     xn=(d12**2-dn**2+dr**2)/(2*d12)
     yn=(sqrt(max(dr**2-xn**2,0.)))
     d13n=sqrt((x-xn)**2.+(y-yn)**2.)
! reposition of v3 in the new frame
     x=(dn**2-d23**2+d13n**2)/(2*dn)
     if (x > 0 .and. x < dn) then
        y=sqrt(max(d13n**2-x**2,0.))
        dplane=r2+y
     else
        dplane=infinity
     endif
  else
! case where v1 remains and V2 is changed
     xn=(d12**2-dr**2+dn**2)/(2*d12)
     yn=(sqrt(max(dn**2-xn**2,0.)))
     d23n=sqrt((x-xn)**2.+(y-yn)**2.)
! reposition of v3 in the new frame
     x=(dn**2-d23n**2+d13**2)/(2*dn)
     if (x > 0 .and. x < dn) then
        y=sqrt(max(d13**2-x**2,0.))
        dplane=r1+y
     else
        dplane=infinity
     endif
  endif
end function dplane
!###############################################################################
function dcircle(d12,d13,d23,r1,r2)
! computes the coordinates in 2D plane of three vertices V1,V2,V3 making the
! assumptions that V1 is the origin (0,0), V2 is on the x axis (0,d12). The
! first part gives the coordinates x,y of V3.
! Same set of equations are used to estimate the origin of a point distant by
! r1 and r2 from V1 and V2 respectively (xc,+/-yc). The function returns the
! distance between (xc,-yc) and V0 (x,y)

  real :: dcircle,d12,d13,d23,r1,r2

  real :: x,y,xc,yc,a
  logical :: line

! TEST
!  dcircle=dplane(d12,d13,d23,r1,r2)
!  return
! by definition the first lateration computing the V3 coordinates is stable
! if the mesh is not ill formed
  x=(d12**2-d23**2+d13**2)/(2*d12)
  y=sqrt(max(d13**2-x**2,0.))
  if (verbose) write(*,*) 'dcircle - x,y :',x,y
! The second lateration may encounter instability when r1 and r2 are small.
!     First test : is it a line source ?
  line=(r1 < epsilon(r1) .and. r2 < epsilon(r2))
!     if r1 OR r2 are null, no trilateration
  if ((r1 < epsilon(r1) .or. r2 < epsilon(r2)).and..not.line) then
     dcircle=infinity
     if (verbose) write(*,*) 'dcircle - r1 or r2 are null'
     return
  endif
!     in case of a line case, dcircle = y if x is in the gate
  if (line) then
     if (x < 0. .or. x > d12) then
        dcircle=infinity
        if (verbose) write(*,*) 'line case but x is not in the gate'
     else
        dcircle=y
     endif
     return
  endif
! This second test may be useless for distance but to be safe ....
  if (abs(r1-r2) < epsilon(r1) .and. r1 < d12*.5) then
     dcircle=infinity
     return
  endif
  xc=(d12**2-r2**2+r1**2)/(2*d12)
  yc=sqrt(max(r1**2-xc**2,0.))
  if (verbose) write(*,*) 'dcircle - xc,yc :',xc,yc
! The gate test: a is the x value of the virtual path at the abscisse axis crossing.
! It must lay between the two vertices, i.e. passing throuth the gate.
  a=abs(xc-x)/(1+(yc/y))
  if (x < xc) then
     a=x+a
  else
     a=x-a
  endif
  if (a < 0. .or. a > d12) then
     dcircle=infinity
     if (verbose) write(*,*) 'dcircle - a : ',a
     if (verbose) write(*,*) 'dcircle - a outside range'
     return
  endif
! Finally ....
  dcircle=sqrt((x-xc)**2+(y+yc)**2.)  ! remember ... distance to (xc,-yc)
  return
end function dcircle
!###############################################################################
function tcircle(d12,d13,d23,t1,t2,v)
! computes the coordinates in 2D plane of three vertices V1,V2,V3 making the
! assumptions that V1 is the origin (0,0), V2 is on the x axis (0,d12). The
! first part gives the coordinates x,y of V3.
! Same set of equations are used to estimate the origin of a point distant by
! r1 and r2 from V1 and V2 respectively (xc,+/-yc). The function returns the
! distance between (xc,-yc) and V0 (x,y)

  real :: tcircle,d12,d13,d23,t1,t2,v

  real :: x,y,xc,yc,a,r1,r2

  x=(d12**2-d23**2+d13**2)/(2*d12)
  y=sqrt(max(d13**2-x**2,0.))
  if (verbose) write(*,*) 'tcircle - x,y :',x,y
  r1=t1*v
  r2=t2*v
  xc=(d12**2-r2**2+r1**2)/(2*d12)
  yc=sqrt(max(r1**2-xc**2,0.))
  if (verbose) write(*,*) 'tcircle - xc,yc :',xc,yc
  a=abs(xc-x)/(1+(yc/y))
  if (x < xc) then
     a=x+a
  else
     a=x-a
  endif
  if (a < 0. .or. a > d12) then
     tcircle=infinity
     if (verbose) write(*,*) 'dcircle - a : ',a
     if (verbose) write(*,*) 'dcircle - a outside range'
     return
  endif
  tcircle=sqrt((x-xc)**2+(y+yc)**2.)/v  ! remember ... distance to (xc,-yc)
  return
end function tcircle
!###############################################################################
subroutine startdist(amesh,anode,nodelist,distarray)

! Compute the distance at the neightbor nodes list (nodelist) of a given node (anode)
! This routine is needed by onevsall routine after the initial vicinode call.
! The routine allvsall uses the updclosedist routine instead

  type(mesh) :: amesh
  type(listn), pointer :: nodelist
  integer :: anode
  real, dimension(amesh%Nnodes) :: distarray

  type(listn), pointer :: pcur

  pcur=>nodelist
  do while (associated(pcur))
     distarray(pcur%idnode)=donedge(amesh,anode,pcur%idnode)
     pcur=>pcur%next
  enddo
  return
end subroutine startdist
!###############################################################################
subroutine startime(amesh,anode,celllist,nodelist,time,velocity)

  type(mesh) :: amesh
  type(listn), pointer :: nodelist,celllist
  integer :: anode
  real, dimension(amesh%Nnodes) :: time
  real, dimension(amesh%Ncells) :: velocity

  type(listn), pointer :: curcell,curnode
  logical :: begin
  integer :: i,j

!  write(*,*) 'startime : entering vicicell for node ',anode
   if (verbose) write(*,*) 'celllist%idnode',celllist%idnode
   begin=.true.
   curcell=>celllist
   do while (associated(curcell))
      do i=1,3
         if (amesh%cell(curcell%idnode,i).ne.anode) then
            if (begin) then
               allocate(nodelist)
               nullify(nodelist%previous)
               nullify(nodelist%next)
               curnode=>nodelist
               curnode%idnode=amesh%cell(curcell%idnode,i)
               begin=.false.
            else
               if (.not.isinit(amesh%cell(curcell%idnode,i),nodelist)) then
                  allocate(curnode%next)
                  curnode%next%previous=>curnode
                  nullify(curnode%next%next)
                  curnode=>curnode%next
                  curnode%idnode=amesh%cell(curcell%idnode,i)
               endif
           endif
           time(amesh%cell(curcell%idnode,i))=min(tonedge(amesh,anode,&
                amesh%cell(curcell%idnode,i),velocity(curcell%idnode)),&
                time(amesh%cell(curcell%idnode,i)))
        endif
     enddo
     curcell=>curcell%next
  enddo
end subroutine startime
!###############################################################################
subroutine vicinode(amesh,anode,celllist,nodelist)

! Compute the neightbor nodes list (nodelist) of a given node (anode)
! and its given neighbor cell list (celllist).

  type(mesh) :: amesh
  type(listn), pointer :: celllist,nodelist
  integer :: anode

  type(listn), pointer :: curcell,curnode
  logical :: begin
  integer :: i,j

  if (verbose) write(*,*) 'celllist%idnode',celllist%idnode
  begin=.true.
  curcell=>celllist
  do while (associated(curcell))
     do i=1,3
        if (amesh%cell(curcell%idnode,i).ne.anode) then
           if (begin) then
              allocate(nodelist)
              nullify(nodelist%previous)
              nullify(nodelist%next)
              curnode=>nodelist
              curnode%idnode=amesh%cell(curcell%idnode,i)
              begin=.false.
           else
              if (.not.isinit(amesh%cell(curcell%idnode,i),nodelist)) then
                 allocate(curnode%next)
                 curnode%next%previous=>curnode
                 nullify(curnode%next%next)
                 curnode=>curnode%next
                 curnode%idnode=amesh%cell(curcell%idnode,i)
              endif
           endif
        endif
     enddo
     curcell=>curcell%next
  enddo
end subroutine vicinode
!###############################################################################
subroutine compntoc(amesh,ntoc)
! compute a node to cell array. ntoc is a list array. it is defined by an array
! pointing at the first element of the cell list.
! warning : but fortran does not understand pointer arrays, so I use a container.
! ntoc is not tecnically a pointer array but an array of container containing a
! pointer.
! "last" is a pointer array which, for each node, points at the end of the cell list.

  type(mesh) :: amesh
  type(container), dimension(amesh%Nnodes) :: ntoc

  type(container), dimension(amesh%Nnodes) :: last
  integer :: i,j
  logical :: begin

  if (verbose) write(*,*) 'entering in compntoc'
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
!       if (verbose) write(*,*) 'last(amesh%cell(i,j))%ptr%idnode',last(amesh%cell(i,j))%ptr%idnode
     enddo
  enddo
  if (verbose) write(*,*) 'check on cell list on node :',amesh%Nnodes/2
  if (verbose) call printlist(ntoc(amesh%Nnodes/2)%ptr)
  if (verbose) write(*,*) 'ntoc(500)%ptr%idnode',ntoc(500)%ptr%idnode
end subroutine compntoc
!###############################################################################
subroutine deallocntoc(ntoc,n)
   integer :: n
   type(container), dimension(n) :: ntoc

   type(listn), pointer :: p1,p2
   integer :: i

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
subroutine dumpQuakeBordervtk(dev,amesh)
  type(mesh) :: amesh
  integer :: dev
  integer,allocatable, dimension(:) :: ordered
  integer :: i,j
  integer :: id

!  call order_nodes(amesh,ordered)

  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a8)') 'distance'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  write(dev,'(a7,i10,a6)') 'POINTS ',amesh%QuakeBorder_NodesNo,' float'
  do i=1,amesh%QuakeBorder_NodesNo
     id = amesh%QuakeBorder_Nodes(i)
     write(dev,*) amesh%px(id),amesh%py(id),amesh%pz(id)
  enddo
  write(dev,'(a9,2i10)') 'POLYGONS ',1,amesh%QuakeBorder_NodesNo+1
!  do i=1,amesh%Ncells
!     write(dev,'(a2,$)') '3 '
     write(dev,'(i4,$)') amesh%QuakeBorder_NodesNo
     write(dev,*) (j-1,j=1,amesh%QuakeBorder_NodesNo)
!  enddo

  return
end subroutine dumpQuakeBordervtk
!###############################################################################
!###############################################################################
subroutine dumpQuakevtk(dev,amesh)
  type(mesh) :: amesh
  integer :: dev
  integer :: i,j
  integer :: id

!  call order_nodes(amesh,ordered)

  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a8)') 'distance'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  write(dev,'(a7,i10,a6)') 'POINTS ',amesh%QuakeNodesNo,' float'
  do i=1,amesh%QuakeNodesNo
     id = amesh%QuakeNodes(i)
     write(dev,*) amesh%px(id),amesh%py(id),amesh%pz(id)
  enddo
  write(dev,'(a9,2i10)') 'POLYGONS ',1,amesh%QuakeNodesNo+1
!  do i=1,amesh%Ncells
!     write(dev,'(a2,$)') '3 '
     write(dev,'(i4,$)') amesh%QuakeNodesNo
     write(dev,*) (j-1,j=1,amesh%QuakeNodesNo)
!  enddo

  return
end subroutine dumpQuakevtk
!###############################################################################
!###############################################################################
subroutine dumpmeshvtk(dev,amesh)
  type(mesh) :: amesh
  integer :: dev

  integer :: i,j

  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a8)') 'distance'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  write(dev,'(a7,i10,a6)') 'POINTS ',amesh%Nnodes,' float'
  do i=1,amesh%Nnodes
     write(dev,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  write(dev,'(a9,2i10)') 'POLYGONS ',amesh%Ncells,amesh%Ncells*4
  do i=1,amesh%Ncells
     write(dev,'(a2,$)') '3 '
     write(dev,*) (amesh%cell(i,j)-1,j=1,3)
  enddo
  return
end subroutine dumpmeshvtk
!###############################################################################
subroutine dumpnodeattributevtk(dev,amesh,field,attname,init)
 type(mesh) :: amesh
  real, dimension(amesh%Nnodes) :: field
  integer :: dev
  character*(*) :: attname
  logical :: init

  integer :: i

  if (init) write(dev,'(a11,i5)') 'POINT_DATA ',amesh%Nnodes
  write(dev,'(a8,a,a8)') 'SCALARS ',trim(attname),' float 1'
  write(dev,'(a20)') 'LOOKUP_TABLE default'
  do i=1,amesh%Nnodes
     write(dev,*) field(i)
  enddo
  return
end subroutine dumpnodeattributevtk
!###############################################################################
subroutine dumpcellattributevtk(dev,amesh,field,attname,init)
 type(mesh) :: amesh
  real, dimension(amesh%Ncells) :: field
  integer :: dev
  character*(*) :: attname
  logical :: init

  integer :: i

  if (init) write(dev,'(a10,i5)') 'CELL_DATA ',amesh%Ncells
  write(dev,'(a8,a,a8)') 'SCALARS ',trim(attname),' float 1'
  write(dev,'(a20)') 'LOOKUP_TABLE default'
  do i=1,amesh%Ncells
     write(dev,*) field(i)
  enddo
  return
end subroutine dumpcellattributevtk
!###############################################################################
!subroutine backtrilat(amesh,distarray,apath,thres)
! aside the definition of the mesh and matrix array distarray, apath has to be
! initialized (allocated) with the starting point coordinates (apath%x, apath%y, apath%z) and
! the face id containing it
! The threshold thres is the condition to exit
! type(mesh) :: amesh
! real, dimension(amesh%Nnodes) :: distarray
! type(path), pointer :: apath
! real :: thres
!
! real :: dtos
! integer, dimension(2) :: idmin
! type(path), pointer :: pcur
!
!
!! creating a cell local distance array
!  ld=distarray(amesh%cell(apath%cellid,:))
!! finding the minima
!  call getmin(idmin,ld)
!! converting indices on local distance array to vertex indices
!  idmin(1)=amesh%cell(apash%cellid,idmin(1))
!  idmin(2)=amesh%cell(apash%cellid,idmin(2))
!  d12=donedge(amesh,idmin(1),idmin(2))
!  d13=dbpt(amesh%px(idmin(1)),amesh%py(idmin(1)),amesh%pz(idmin(1)),apath%pos%x,apath%pos%y,apath%pos%z)
!  d23=dbpt(amesh%px(idmin(2)),amesh%py(idmin(2)),amesh%pz(idmin(2)),apath%pos%x,apath%pos%y,apath%pos%z)
!  r1=distarray(idmin(1))
!  r2=distarray(idmin(2))
!  call stepback(d12,d13,d23,r1,r2,a,dtos)
!  allocate(apath%next)
!  nullify(apath%next%next)
!  pcur=>apath%next
!  pcur%pos%x=(amesh%px(idmin(2))-amesh%px(idmin(1)))/d12*a+amesh%px(idmin(1))
!  pcur%pos%y=(amesh%py(idmin(2))-amesh%py(idmin(1)))/d12*a+amesh%py(idmin(1))
!  pcur%pos%x=(amesh%pz(idmin(2))-amesh%pz(idmin(1)))/d12*a+amesh%pz(idmin(1))
!  pcur%cellid=nextcell(amesh,apath%cellid,idmin)
!  do while (dtos < thres)
!     ld=distarray(amesh%cell(pcur%cellid,:))
!     call getmin(idmin,ld)
!     idmin(1)=amesh%cell(pcur%cellid,idmin(1))
!     idmin(2)=amesh%cell(pcur%cellid,idmin(2))
!     d12=donedge(amesh,idmin(1),idmin(2))
!     d13=dbpt(amesh%px(idmin(1)),amesh%py(idmin(1)),amesh%pz(idmin(1)),pcur%pos%x,pcur%pos%y,pcur%pos%z)
!     d23=dbpt(amesh%px(idmin(2)),amesh%py(idmin(2)),amesh%pz(idmin(2)),pcur%pos%x,pcur%pos%y,pcur%pos%z)
!     r1=distarray(idmin(1))
!     r2=distarray(idmin(2))
!     call stepback(d12,d13,d23,r1,r2,a,dtos)
!     allocate(pcur%next)
!     nullify(pcur%next%next)
!     remind=pcur%cellid
!     pcur=>pcur%next
!     pcur%pos%x=(amesh%px(idmin(2))-amesh%px(idmin(1)))/d12*a+amesh%px(idmin(1))
!     pcur%pos%y=(amesh%py(idmin(2))-amesh%py(idmin(1)))/d12*a+amesh%py(idmin(1))
!     pcur%pos%x=(amesh%pz(idmin(2))-amesh%pz(idmin(1)))/d12*a+amesh%pz(idmin(1))
!     pcur%cellid=nextcell(amesh,remind,idmin)
!  enddo
!end subroutine backtrilat
!
!!###############################################################################
!subroutine nextcell(amesh,curcell,idmin)
!   integer :: nextcell,curcell
!   type(mesh) :: amesh
!   integer, dimension(2) :: idmin
!
!   logical :: jfound
!   integer :: i,j,k,im
!
!  jfound=.false.
!  i=1
!  do while (.not.jfound)
!     im=0
!     if (i == curcell) cycle
!     do j=1,3
!        do k=1,2
!           if (amesh%cell(i,j) == idmin(k)) im=im+1
!        enddo
!     enddo
!     if (im == 2) then
!        jfound=.true.
!     else
!        i=i+1
!     endif
!     if (i > amesh%Ncells) stop ' something wrong happened'
!  enddo
!  return
!end subroutine nextcell
!###############################################################################
!subroutine getmin(idmin,ld)
!! Warning : in this case, idmin returns the indices containing the two
!! smallest values in the local distance array ld
!  integer, dimension(2) :: idmin
!  real, dimension(3) :: ld
!
!  integer :: ml
!
!  ml=maxloc(ld)
!  if (maxloc == 1) then
!     idmin(1)=2
!     idmin(2)=3
!     return
!  endif
!  if (maxloc == 2) then
!     idmin(1)=1
!     idmin(2)=3
!     return
!  endif
!  idmin(1)=1
!  idmin(2)=2
!  return
!end subroutine getmin
!!###############################################################################
end module lateration
