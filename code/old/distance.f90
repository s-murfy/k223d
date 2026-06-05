module distance
use generic
implicit none

real(pr), allocatable, dimension(:,:) :: dist


contains
!##########################################################
subroutine allvsall2d(amesh,distarray)    ! CHECK: added this from lateration.f90 
   use LAT_mesh
   use lists

   type(mesh) :: amesh
   real(pr), dimension(:,:), allocatable :: distarray
 
   type(container), dimension(amesh%Nnodes) :: ntoc
   type(listn), pointer :: ntodo,cellcur,pcur,last
   real(pr) :: d01,d02,d12,r1,r2,dface,dedge,dtest
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
 !     call printperc(k,amesh%Nnodes)
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
!##########################################################
subroutine pre_onevsall2d_offvertex(amesh,distarray,v,d)
! initialisation of onevsall2d in case of a unique source at the
! distance d(1), d(2) and d(3) from vertice v(1), v(2) and v(3) respectively of the SAME face.
! Use this example to create your own routine to initialize distarray
! for a generic multi source problem with n points and not only 3.
use LAT_mesh

  type(mesh) :: amesh
  real(pr), dimension(:), allocatable :: distarray
  integer(pin), dimension(3) :: v
  real(pr), dimension(3) :: d

  integer(pin) :: i

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
!##########################################################
subroutine pre_onevsall2d_onvertex(amesh,k,distarray)
! initialisation of onevsall2d in case of a unique source on a vertex k
use LAT_mesh

  type(mesh) :: amesh
  real(pr), dimension(:), allocatable :: distarray
  integer(pin) :: k

!   allocation
  allocate(distarray(amesh%Nnodes))
!   initialisation to infinity
  distarray=infinity
!   initialization of distance at k to zero
  distarray(k)=0._pr
  return
end subroutine pre_onevsall2d_onvertex


!###############################################################################
! subroutine onevsall2d(amesh,dist,fast)
!    ! Compute the distance, in an array format, for all the vertices of the mesh
!    ! amesh which distance is not set to infinity.
!    ! Warning : dist should be initialized before to infinity for all the
!    ! vertices except those representing the sources.
!    use LAT_mesh
!    use lists
!      type(mesh) :: amesh
!      real(pr), dimension(amesh%Nnodes) :: dist
!      logical :: fast
   
!      type(container), dimension(amesh%Nnodes) :: ntoc
!      type(listn), pointer :: cellcur,pcur,last,pmin
!      type(listn), pointer :: ntodo
!      real(pr) :: d13,d23,d12,r1,r2,dface,dedge,dtest,db
!      integer(pin), dimension(2) :: otherninc
!      integer(pin) :: i,j,toggle,idiff,waitfordiff
!      logical :: begin,hasbeenupdated,reloop,firstrun,diffoccur
!      integer(pin) :: nswp,mxswp
!      real(pr), dimension(amesh%Nnodes) :: distarray
!      logical, dimension(amesh%Nnodes) :: checksecondary
!      integer(pin), parameter :: waitdiffthres=25
!    !   verbose = .true.
!      write(*,*) 'check 1'
!    !  mxswp=0
!    ! computing the node-to-cell array ntoc
!      call compntoc(amesh,ntoc)
!    ! initialize the secondary distance array
!      distarray=dist
!      checksecondary=.false.
!    ! initialization of ntodo from non infinite values inside distarray
!      begin=.true.
!      do i=1,amesh%Nnodes
!         if (distarray(i) > .9*infinity) cycle
!    ! if a vertex has a null distance, it is not a secondary diffraction point
!         if (distarray(i) < epsilon(distarray(i))) checksecondary(i)=.true.
!         if (begin) then
!            allocate(ntodo)
!            nullify(ntodo%previous)
!            nullify(ntodo%next)
!            pcur=>ntodo
!            begin=.false.
!         else
!            allocate(pcur%next)
!            nullify(pcur%next%next)
!            pcur%next%previous=>pcur
!            pcur=>pcur%next
!         endif
!         pcur%idnode=i
!      enddo
!    ! initializing reloop to true
!    ! The condition is changed or not at the end of the loop
!    ! The fast first loop is mandatory with a null distance offset
!      reloop=.true.
!      firstrun=.true.
!      diffoccur=.false.
!      db=0.
!      do while (reloop)
!    ! loop on the to do list
!      if (verbose) call printlist(ntodo)
!      if (verbose) write(*,*) 'entering in the ntodo list management loop'
!      do while (associated(ntodo))
!            if (verbose) write(*,*) 'state of ntodo :'
!            if (verbose) call printlist(ntodo)
!    ! searching the cells attached to the first element of the to do list
!            call lookformin(ntodo,distarray,amesh%Nnodes,pmin)
!            if (verbose) write(*,*) 'propagating distance from node #',pmin%idnode
!            if (verbose) write(*,*) 'its own distance is :',distarray(pmin%idnode)
!            hasbeenupdated=.true.
!            toggle=2
!    !        nswp=0
!            write(*,*) 'check 2a'

!            do while (hasbeenupdated)
!    !        nswp=nswp+1
!            hasbeenupdated=.false.
!            toggle=3-toggle
!            if (verbose) then
!            if (toggle==1) write(*,*) 'seep forward on cellist'
!            if (toggle==2) write(*,*) 'seep backward on cellist'
!            endif
!            if (toggle==1) then
!               pcur=>ntoc(pmin%idnode)%ptr
!            else
!               pcur=>last
!            endif
!            do while (associated(pcur))
!    !     find the two complemantory nodes in the current cell
!                if (verbose) write(*,*) '    working on cell #',pcur%idnode
!                call givencomp(amesh,otherninc,pcur,pmin)
!                if (verbose) write(*,*) '    complemantory nodes : ',otherninc
!    !     loop on the two complemantory nodes
!                do i=1,2
!                   if (verbose) write(*,'(a4,3(i4.4))') 'code',pmin%idnode,pcur%idnode,otherninc(i)
!                   if (verbose) write(*,*) '       working on complemantory nodes :',otherninc(i)
!    !     edge propagation
!                   dedge=distarray(pmin%idnode)+donedge(amesh,pmin%idnode,otherninc(i))
!                   if (verbose) write(*,*) 'dedge : ',distarray(pmin%idnode),'+',donedge(amesh,pmin%idnode,otherninc(i))
!                   if (verbose) write(*,*) 'dedge : ',dedge
!    !     face propagation
!                   if (distarray(otherninc(3-i)) /= infinity) then
!                   if (verbose) write(*,*) 'computing dface because ',distarray(otherninc(3-i)),' is not infinity'
!                      d12=donedge(amesh,pmin%idnode,otherninc(3-i))
!                      d13=donedge(amesh,otherninc(i),pmin%idnode)
!                      d23=donedge(amesh,otherninc(i),otherninc(3-i))
!                   if (verbose) write(*,*) 'd12,d13,d23 :'
!                   if (verbose) write(*,*) d12,d13,d23
!                      r1=distarray(pmin%idnode)
!                      r2=distarray(otherninc(3-i))
!                   if (verbose) write(*,*) 'r1,r2 : ',r1,r2
!                      dface=dcircle(d12,d13,d23,r1,r2)
!                   if (verbose) write(*,*) 'dface : ',dface
!                   else
!                      dface=infinity
!                   endif
!                   dtest=min(dedge,dface)
!    !              dtest=dedge
!                   if (verbose) write(*,*) 'dtest=min(dedge,dface) : ',dtest
!                   if (verbose) write(*,*) 'dtest < actual minimum on ',otherninc(i)
!                   if (verbose) write(*,*) dtest,distarray(otherninc(i))
!                   if (dtest < distarray(otherninc(i))) then
!                      hasbeenupdated=.true.
!                      if (verbose) write(*,*) 'better !'
!    ! distance is better : if not in the list, add it
!                      distarray(otherninc(i))=dtest
!                      call updatelist(otherninc(i),ntodo)
!                      if (verbose) write(*,*) 'updatelist with :',otherninc(i)
!                      if (verbose) write(*,*) 'state of ntodo :'
!                      if (verbose) call printlist(ntodo)
!                   endif
!                enddo ! complementary loop on cell
!    ! moving on the vicinity list
!                if (toggle==1) then
!                   last=>pcur
!                   pcur=>pcur%next
!                else
!                   pcur=>pcur%previous
!                endif
!            enddo ! neighbor cell to the current node pmin
!            enddo ! sweep loop
!    !        mxswp=max(mxswp,nswp)
!    ! diff case and no diffraction has been detected yet
!            if (.not.firstrun.and..not.diffoccur) then
!               if (distarray(pmin%idnode)+dist(idiff) < dist(pmin%idnode)) then
!    !   diffraction detected
!                  diffoccur=.true.
!               else
!    !   no diffraction
!                  waitfordiff=waitfordiff-1
!                  if (waitfordiff == 0) then
!    !   countdown expired for diffraction. hard exit on ntodo list to pass to another
!    !                                      potential secondary source
!                     call wipe(ntodo)
!                     exit
!                  endif
!               endif
!            endif
!    ! removing the treated element pmin of the to do list
!            write(*,*) 'check 2aa'
!            call printlist(ntodo)

!            call removepminfromthelist(pmin,ntodo)
!            write(*,*) 'check 2aaa'

!          enddo ! ntodo loop
!    ! decision to reloop or not
!    ! fast case
!          write(*,*) 'check 2b'
!          if (fast) then
!             dist=distarray
!    ! first case of general exit : only one run is requested
!             reloop=.false.
!             cycle
!          endif
!          if (firstrun) then
!             dist=distarray
!          else
!             if (diffoccur) then
!    ! a secondary diff has occured
!                do i=1,amesh%Nnodes
!                   dist(i)=min(distarray(i)+dist(idiff),dist(i))
!                enddo
!             endif
!          endif
!          idiff=nextdiffid(dist,amesh%Nnodes,checksecondary)
!    ! second case of general exit : all the potential secondary source have been investigated.
!          if (idiff == 0) then
!             reloop=.false.
!             cycle
!          endif
!    ! preparation for the next loop
!          checksecondary(idiff)=.true.
!          waitfordiff=waitdiffthres
!          diffoccur=.false.
!    ! At this point notodo is deallocated. reallocate it with idiff
!          allocate(ntodo)
!          nullify(ntodo%previous)
!          nullify(ntodo%next)
!          ntodo%idnode=idiff
!    ! reinitialize distarray
!          distarray=infinity
!          distarray(idiff)=0.
!          firstrun=.false.
!      enddo ! end reloop section
!    ! deallocating ntoc
!      write(*,*) 'check end'
!      call deallocntoc(ntoc,amesh%Nnodes)
!    end subroutine onevsall2d
! ########################################################
subroutine onevsall2d(amesh,sdist,fast)
! Compute the distance, in an array format, for all the vertices of the mesh
! amesh which distance is not set to infinity.
! Warning : dist should be initialized before to infinity for all the
! vertices except those representing the sources.
use LAT_mesh
use lists

  type(mesh) :: amesh
  real(pr), dimension(amesh%Nnodes) :: sdist
  logical :: fast

  type(container), dimension(amesh%Nnodes) :: ntoc
  type(listn), pointer :: cellcur,pcur,last,pmin
  type(listn), pointer :: ntodo
  real(pr) :: d13,d23,d12,r1,r2,dface,dedge,dtest,db
  integer(pin), dimension(2) :: otherninc
  integer(pin) :: i,j,toggle,idiff,waitfordiff
  logical :: begin,hasbeenupdated,reloop,firstrun,diffoccur
  integer(pin) :: nswp,mxswp
  real(pr), dimension(amesh%Nnodes) :: distarray
  logical, dimension(amesh%Nnodes) :: checksecondary,inthelist
  integer(pin), parameter :: waitdiffthres=25

!  mxswp=0
! computing the node-to-cell array ntoc
  call compntoc(amesh,ntoc)
! initialize the secondary distance array
  distarray=sdist   ! distarray = secondary sources
  checksecondary=.false.
  inthelist=.false.
! initialization of ntodo from non infinite values inside distarray
  begin=.true.
  ! sound we be using dist rather than distarray (?)
  do i=1,amesh%Nnodes
     if (distarray(i) > .9_pr*infinity) cycle
! if a vertex has a null distance, it is not a secondary diffraction point
   ! if (.not.fast) checksecondary(i)=.true.   ! <sounds better>
   if (distarray(i) < epsilon(distarray(i))) checksecondary(i)=.true.
     if (begin) then
        allocate(ntodo)
        nullify(ntodo%previous)
        nullify(ntodo%next)
        pcur=>ntodo
        begin=.false.
      !   begin=.true.
     else
        allocate(pcur%next)
        nullify(pcur%next%next)
        pcur%next%previous=>pcur
        pcur=>pcur%next
     endif
     pcur%idnode=i
     inthelist(i)=.true.
     ntodo%previous=>pcur
  enddo
! initializing reloop to true
! The condition is changed or not at the end of the loop
! The fast first loop is mandatory with a null distance offset
  reloop=.true.
  firstrun=.true.
  diffoccur=.false.
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
!                  call updatelist(otherninc(i),ntodo)
                  if (.not.inthelist(otherninc(i))) then
                      allocate(ntodo%previous%next)
                      nullify(ntodo%previous%next%next)
                      ntodo%previous%next%previous=>ntodo%previous
                      ntodo%previous%next%idnode=otherninc(i)
                      ntodo%previous=>ntodo%previous%next
                      inthelist(otherninc(i)) = .true.  !added this line, to correct for memory link
                  endif
!                  if (verbose) write(*,*) 'updatelist with :',otherninc(i)
!                  if (verbose) write(*,*) 'state of ntodo :'
!                  if (verbose) call printlist(ntodo)
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
           if (distarray(pmin%idnode)+sdist(idiff) < sdist(pmin%idnode)) then
!   diffraction detected
              diffoccur=.true.
           else
!   no diffraction
              waitfordiff=waitfordiff-1   ! check to see if there is an improvement in cells near diff point
              if (waitfordiff == 0) then
!   countdown expired for diffraction. hard exit on ntodo list to pass to another
!                                      potential secondary source
                 call wipe(ntodo)
                 exit
              endif
           endif
        endif
! removing the treated element pmin of the to do list
        inthelist(pmin%idnode)=.false.
        call removepminfromthelist(pmin,ntodo)
      enddo ! ntodo loop
! decision to reloop or not
! fast case
      if (fast) then
         sdist=distarray
! first case of general exit : only one run is requested
         reloop=.false.
         cycle
      endif
      if (firstrun) then
         sdist=distarray
      else
         if (diffoccur) then
! a secondary diff has occured
            do i=1,amesh%Nnodes
               sdist(i)=min(distarray(i)+sdist(idiff),sdist(i))
            enddo
         endif
      endif
! pick next diffraction point
      idiff=nextdiffid(sdist,amesh%Nnodes,checksecondary)
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
      ntodo%previous=>ntodo
      nullify(ntodo%next)
      ntodo%idnode=idiff
      inthelist=.false.
      inthelist(idiff)=.true.
! reinitialize distarray
      distarray=infinity
      distarray(idiff)=0.
      firstrun=.false.
  enddo ! end reloop section
! deallocating ntoc
  call deallocntoc(ntoc,amesh%Nnodes)
!  write(0,*) 'maximun sweep numner : ',mxswp
end subroutine onevsall2d
!#########################################################
function nextdiffid(dist,n,checksecondary)
! find the minimum value in dist skipping the vertex checked yet and returns
! its position (id)
  integer(pin) :: nextdiffid,n
  real(pr), dimension(n) :: dist
  logical, dimension(n) :: checksecondary

  integer(pin) :: i
  real(pr) :: vmin

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
!########################################################>>>>>>>>>>>>>>>>>>>> module distance
subroutine lookformin(ntodo,anarray,n,pmin)
  use lists

  type(listn), pointer :: pmin,ntodo
  integer(pin) :: n
  real(pr), dimension(n) :: anarray

  type(listn), pointer :: pcur
  real(pr) :: vmin

  pcur=>ntodo
  vmin=anarray(pcur%idnode)
  pmin=>pcur
  do while (associated(pcur%next))
     pcur=>pcur%next
     if (anarray(pcur%idnode) < vmin) then
        vmin=anarray(pcur%idnode)
        pmin=>pcur
     endif
  enddo
  return
end subroutine lookformin

!###############################################################################
function donedge(amesh,i,j)
use LAT_mesh

  type(mesh) :: amesh
  integer(pin) :: i,j
  real(pr) :: donedge

  donedge=sqrt((amesh%px(i)-amesh%px(j))**2._pr+&
               (amesh%py(i)-amesh%py(j))**2._pr+&
               (amesh%pz(i)-amesh%pz(j))**2._pr)
  return
end function donedge
!###############################################################################
subroutine givencomp(amesh,otherninc,pc,pa)
  use LAT_mesh
  use lists

  type(mesh) :: amesh
  integer(pin), dimension(2) :: otherninc
  type(listn), pointer :: pc,pa

  integer(pin) :: i,k

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
function dcircle(d12,d13,d23,r1,r2)
! computes the coordinates in 2D plane of three vertices V1,V2,V3 making the
! assumptions that V1 is the origin (0,0), V2 is on the x axis (0,d12). The
! first part gives the coordinates x,y of V3.
! Same set of equations are used to estimate the origin of a point distant by
! r1 and r2 from V1 and V2 respectively (xc,+/-yc). The function returns the
! distance between (xc,-yc) and V0 (x,y)

  real(pr) :: dcircle,d12,d13,d23,r1,r2
  real(pr) :: x,y,xc,yc,a
  logical :: line

! TEST
!  dcircle=dplane(d12,d13,d23,r1,r2)
!  return
! by definition the first lateration computing the V3 coordinates is stable
! if the mesh is not ill formed
  x=(d12**2._pr-d23**2._pr+d13**2._pr)/(2._pr*d12)
  y=sqrt(max(d13**2._pr-x**2._pr,0._pr))
! if (verbose) write(*,*) 'dcircle - x,y :',x,y
! The second lateration may encounter instability when r1 and r2 are small.
!     First test : is it a line source ?
  line=(r1 < epsilon(r1) .and. r2 < epsilon(r2))
!     if r1 OR r2 are null, no trilateration
  if ((r1 < epsilon(r1) .or. r2 < epsilon(r2)).and..not.line) then
     dcircle=infinity
!    if (verbose) write(*,*) 'dcircle - r1 or r2 are null'
     return
  endif
!     in case of a line case, dcircle = y if x is in the gate
  if (line) then
     if (x < 0._pr .or. x > d12) then
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
  xc=(d12**2._pr-r2**2._pr+r1**2._pr)/(2._pr*d12)
  yc=sqrt(max(r1**2._pr-xc**2._pr,0._pr))
  if (verbose) write(*,*) 'dcircle - xc,yc :',xc,yc
! The gate test. a is the x value of the virtual path at the abscisse axis crossing.
! It must lay between the two vertices, i.e. passing throuth the gate.
  a=abs(xc-x)/(1._pr+(yc/y))
  if (x < xc) then
     a=x+a
  else
     a=x-a
  endif
  if (a < 0._pr .or. a > d12) then
     dcircle=infinity
!    if (verbose) write(*,*) 'dcircle - a : ',a
!    if (verbose) write(*,*) 'dcircle - a outside range'
     return
  endif
! Finally ....
  dcircle=sqrt((x-xc)**2._pr+(y+yc)**2._pr)  ! remember ... distance to (xc,-yc)
  return
end function dcircle
!###############################################################################
subroutine compntoc(amesh,ntoc)
! compute a node to cell array. ntoc is a list array. it is defined by an array
! pointing at the first element of the cell list.
! warning : but fortran does not understand pointer arrays, so I use a container.
! ntoc is not tecnically a pointer array but an array of container containing a
! pointer.
! "last" is a pointer array which, for each node, points at the end of the cell list.
  use LAT_mesh
  use lists

  type(mesh) :: amesh
  type(container), dimension(amesh%Nnodes) :: ntoc

  type(container), dimension(amesh%Nnodes) :: last
  integer(pin) :: i,j
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
   use lists

   integer(pin) :: n
   type(container), dimension(n) :: ntoc

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
!###############################################################################
subroutine updclosedist(amesh,distarray,k,i,j)
   use LAT_mesh
   type(mesh) :: amesh
   real(pr), dimension(:,:), allocatable :: distarray
   integer(pin) :: i,j,k

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
!###############################################################################
subroutine vicinode(amesh,anode,celllist,nodelist)
   use LAT_mesh
   use lists 
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

end module distance
