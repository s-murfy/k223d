module distance
use generic
implicit none
!! main module to compute the geodetic distance through lateration
!! the two main subroutines are
!!  - onevsall
!!  - allvsall

contains
!##########################################################
subroutine allvsall2d(amesh,dist,adiff)
! 2025E new version of allvsall2d based on a loop over an onevsall2d call
   use lists
   use LAT_mesh
   use LAT_time

! input variables

   type(mesh) :: amesh ! the mesh structure
   real(pr), dimension(amesh%Nnodes,amesh%Nnodes) :: dist ! the distance matrix node to node
   type(diff) :: adiff ! the diffraction structure

! local variables
   integer(pin) :: i, j, k ! loop index
   real(pr), dimension(amesh%Nnodes) :: distarray ! the temporary distance array
   type(containerc), dimension(amesh%Nnodes) :: ntoc ! the node to cell array (sparse matrix)
   type(containern), dimension(amesh%Nnodes) :: nton ! the node to node array (sparse matrix)

! computing the node-to-cell array ntoc
   call compntoc(amesh,ntoc)
! computing the node-to-node array nton
   call compnton(amesh,nton)
! initialisation to infinity except for the trace to zero
   dist=infinity
!   trace to zero
   do i=1,amesh%Nnodes
      dist(i,i)=0.
   enddo
! loop on the nodes
   do k=1,amesh%Nnodes-1
! initialisation of the distance array
      distarray(:)=dist(k,:)
      call onevsall2d(amesh,distarray,ntoc,nton,adiff)
! updating the distance matrix
      do i=k+1,amesh%Nnodes
         dist(k,i)=distarray(i)
         dist(i,k)=distarray(i)
      enddo
   enddo
   call free_nton(nton)
   call free_ntoc(ntoc)
end subroutine allvsall2d

!##########################################################
subroutine pre_onevsall2d_offvertex(amesh,distarray,v,d,nv,ntoc,nton)
! initialisation of onevsall2d in case of multiple sources (nv sources)
! defined by the node index v and initialized with the value d
! May be used in case of sources off vertex.
! For example we may insert the distances d(1), d(2) and d(3) from nodes
! v(1), v(2) and v(3) respectively of the SAME cell.
! Use this example to create your own routine to initialize distarray
! for a generic multi source problem.
! 2024A the computation of ntoc is moved here
use LAT_mesh
use lists

  type(mesh) :: amesh
  real(pr), dimension(amesh%Nnodes) :: distarray
  integer(pin) :: nv
  integer(pin), dimension(nv) :: v
  real(pr), dimension(nv) :: d
  type(containerc), dimension(amesh%Nnodes) :: ntoc
  type(containern), dimension(amesh%Nnodes) :: nton

  integer(pin) :: i

! computing the node-to-cell array ntoc
  call compntoc(amesh,ntoc)
! 2024C computing the node-to-node edge array nton
  call compnton(amesh,nton)
!   initialisation to infinity
  distarray=infinity
!   initialisation of distarray and ntodo
  do i=1,nv
     distarray(v(i))=d(i)
  enddo

  return
end subroutine pre_onevsall2d_offvertex
!##########################################################
subroutine pre_onevsall2d_onvertex(amesh,k,distarray,ntoc,nton)
! initialisation of onevsall2d in case of a unique source on a vertex k
! 2024A the computation of ntoc is moved here
! 2024B added computation of nton
use LAT_mesh
use lists

  type(mesh) :: amesh
  real(pr), dimension(amesh%Nnodes) :: distarray
  integer(pin) :: k
  type(containerc), dimension(amesh%Nnodes) :: ntoc
! 2024C
  type(containern), dimension(amesh%Nnodes) :: nton

! computing the node-to-cell array ntoc
  call compntoc(amesh,ntoc)
! 2024C computing the node-to-node edge array nton
  call compnton(amesh,nton)
!   initialisation to infinity
  distarray=infinity
!   initialization of distance at k to zero
  distarray(k)=0._pr
  return
end subroutine pre_onevsall2d_onvertex
!########################################################
subroutine onevsall2d(amesh,dist,ntoc,nton,adiff)
! Compute the distance, in an array format, for all the vertices of the mesh
! amesh which distance is not set to infinity.
! Warning : dist should be initialized before to infinity for all the
! vertices except those representing the sources.
! 2024A the computation of ntoc has been moved in the pre computation phase
! 2024A Adding a logical to track the initial nodes (intheinput)
! 2025A Assessed has been checked and removed
! 2025B Adding a tempering to avoid unecessary trigging of the diffraction
! 2025C commenting the variables  

   use LAT_mesh
   use LAT_time
   use lists

! inout variables

  type(mesh) :: amesh   ! the mesh structure
  real(pr), dimension(amesh%Nnodes) :: dist  ! the distance array for each node
  type(diff) :: adiff ! the diffraction structure
  type(containerc), dimension(amesh%Nnodes) :: ntoc ! the node to cell array (sparse matrix)
  type(containern), dimension(amesh%Nnodes) :: nton ! the node to node array (sparse matrix)

! local variables

  type(listn), pointer :: cellcur,pcur,last,pmin ! pointers to list elements
  type(listn), pointer :: ntodo ! the list of nodes to be treated
  real(pr) :: d13,d23,d12 ! edge distances
  real(pr) :: r1,r2 ! distances between the virtual source and the two vertices on the x axis
  real(pr) :: dface,dedge,dtest ! distances of the different operators
  real(pr) :: dmem ! temporary variable for the swap of d12 and d13
  integer(pin), dimension(2) :: otherninc ! the two complementary nodes off the origin
  integer(pin) :: i,j ! loop index
  integer (pin) :: toggle ! toggle the for the loop
  integer (pin) :: idiff ! the id of the current diffraction point
  integer (pin) :: waitfordiff ! countdown for the diffraction detection
  logical :: begin ! flag to initialize the ntodo list
  logical :: hasbeenupdated ! flag to check if the distance has been updated
  logical :: reloop ! flag to check if the loop should be repeated
  logical :: firstrun ! flag to check if it is the first run
  logical :: diffoccur ! flag to check if a diffraction has occurred
  real(pr), dimension(amesh%Nnodes) :: distarray ! the temporary distance array
  ! 2025C needed work on checksecondary : At a moment, this array may be defined prior
  ! to the call to onevsall2d. It should be an inout variable 
  logical, dimension(amesh%Nnodes) :: checksecondary ! flag to check if the node has to be checked
  logical, dimension(amesh%Nnodes) :: inthelist ! mark the defined nodes at the start
  ! 2025C Is intheinput needed ?
  logical, dimension(amesh%Nnodes) :: intheinput ! for the moment, this is the same as inthelist

  ! 2025A
  integer(pin), parameter :: waitdiffthres=25 ! threshold for the diffraction detection

  verbose=3
! initialize the secondary distance array
  distarray=dist   ! distarray = secondary sources
  checksecondary=.true.
  inthelist=.false.
! initialization of ntodo from non infinite values inside distarray
  begin=.true.

! mark the nodes in checksecondary to false if there are declared diffractors in
! the diff structure. 
  if (.not.adiff%fast) then
      do i = 1,adiff%NdiffNodes
         checksecondary(adiff%nodes(i))=.false.
      enddo
  endif

  ! sound we be using dist rather than distarray (?)
  do i=1,amesh%Nnodes
     if (distarray(i) > .9_pr*infinity) cycle
! if a vertex has a null distance, it is not a secondary diffraction point
!   if (.not.fast) checksecondary(i)=.true.   ! <sounds better>
     if (distarray(i) < epsilon(distarray(i))) checksecondary(i)=.true.
     if (begin) then
        allocate(ntodo)
        nullify(ntodo%previous)
        nullify(ntodo%next)
        pcur=>ntodo
        begin=.false.
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
! 2024A intheinput is the duplicate of inthelist at the beginning
  intheinput=inthelist ! 2024A
! initializing reloop to true
! The condition is changed or not at the end of the loop
! The fast first loop is mandatory with a null distance offset
  reloop=.true.
  firstrun=.true.
  diffoccur=.false.
  do while (reloop)
! loop on the to do list
  if (verbose==2) call printlist(ntodo)
  if (verbose==2) write(*,*) 'entering in the ntodo list management loop'
  do while (associated(ntodo))
        if (verbose==2) write(*,*) 'state of ntodo :'
        if (verbose==2) call printlist(ntodo)
! searching the cells attached to the first element of the to do list
        call lookformin(ntodo,distarray,amesh%Nnodes,pmin)
        if (verbose==2) write(*,*) 'propagating distance from node #',pmin%idnode
        if (verbose==2) write(*,*) 'its own distance is :',distarray(pmin%idnode)
        hasbeenupdated=.true.
        toggle=2

        do while (hasbeenupdated)
        hasbeenupdated=.false.
        toggle=3-toggle
        if (verbose==2) then
        if (toggle==1) write(*,*) 'seep forward on cellist'
        if (toggle==2) write(*,*) 'seep backward on cellist'
        endif
        if (toggle==1) then
           pcur=>ntoc(pmin%idnode)%ptr
        else
           pcur=>last
        endif
        do while (associated(pcur))  ! NB pcur is a pointer on a cell
!     find the two complemantory nodes in the current cell
            if (verbose==2) write(*,*) '    working on cell #',pcur%idnode
            call givencomp(amesh,otherninc,pcur,pmin)
            if (verbose==2) write(*,*) '    complemantory nodes : ',otherninc
!     2024A moving edge computation of the cell pcur before the inner loop
!     2024C the access to the distance needs an small interface to nton
!           idea... it is faster to create a routine to seek the 3 distances alltogether 
            d12=get_donedge(amesh,nton,pmin%idnode,otherninc(2))
            d13=get_donedge(amesh,nton,otherninc(1),pmin%idnode)
            d23=get_donedge(amesh,nton,otherninc(1),otherninc(2))
!     loop on the two complemantory nodes
            do i=1,2
               if (verbose==2) write(*,'(a4,3(i4.4))') 'code',pmin%idnode,pcur%idnode,otherninc(i)
               if (verbose==2) write(*,*) '       working on complemantory nodes :',otherninc(i)
!     2024A cycle if otherninc(i) is intheinput
               if (intheinput(otherninc(i))) cycle
!     2024A swap distance d12 and d13 for the second otherninc
               if (i == 2) then
                  dmem=d12
                  d12=d13
                  d13=dmem
               endif
!     2024A edge propagation
               dedge=distarray(pmin%idnode)+d13
               if (verbose==2) write(*,*) 'dedge : ',distarray(pmin%idnode),'+',get_donedge(amesh,nton,pmin%idnode,otherninc(i))
               if (verbose==2) write(*,*) 'dedge : ',dedge
!     face propagation
               if (distarray(otherninc(3-i)) /= infinity) then
                  if (verbose==2) write(*,*) 'computing dface because ',distarray(otherninc(3-i)),' is not infinity'
                  if (verbose==2) write(*,*) 'd12,d13,d23 :'
                  if (verbose==2) write(*,*) d12,d13,d23
                  r1=distarray(pmin%idnode)
                  r2=distarray(otherninc(3-i))
                  if (verbose==2) write(*,*) 'r1,r2 : ',r1,r2
                  dface=dcircle(d12,d13,d23,r1,r2)
                  if (verbose==2) write(*,*) 'dface : ',dface
               else
                  dface=infinity
               endif
               dtest=min(dedge,dface)
               if (verbose==2) write(*,*) 'dtest=min(dedge,dface) : ',dtest
               if (verbose==2) write(*,*) 'dtest < actual minimum on ',otherninc(i)
               if (verbose==2) write(*,*) dtest,distarray(otherninc(i))
               if (dtest < distarray(otherninc(i))) then
                  hasbeenupdated=.true.
                  if (verbose==2) write(*,*) 'better !'
                  distarray(otherninc(i))=dtest
! distance is better : if not in the list, add it
! but ...
!           2024D Adding a condition to add to the list only if it is better than the first run
                  if (.not.firstrun) then
                     if (diffoccur) then
                        if (distarray(otherninc(i))+dist(idiff) > 1.01*dist(otherninc(i))) then
!                           write(*,*) 'WARNING ! case occuring'
                           cycle 
                        endif
                     endif
                  endif
                  if (.not.inthelist(otherninc(i))) then
                      allocate(ntodo%previous%next)
                      nullify(ntodo%previous%next%next)
                      ntodo%previous%next%previous=>ntodo%previous
                      ntodo%previous%next%idnode=otherninc(i)
                      ntodo%previous=>ntodo%previous%next
                      inthelist(otherninc(i))=.true.
                  endif
!                  if (verbose==2) write(*,*) 'updatelist with :',otherninc(i)
!                  if (verbose==2) write(*,*) 'state of ntodo :'
!                  if (verbose==2) call printlist(ntodo)
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
! 2025B adding a tempering to avoid unecessary trigging of the diffraction
           if (distarray(pmin%idnode)+dist(idiff) < 0.99*dist(pmin%idnode)) then
!   diffraction detected
              diffoccur=.true.
              if (verbose==3) write(*,*) 'diffraction detected at node #',idiff
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
      if (adiff%fast) then
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
! pick next diffraction point
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
! if (verbose==2) write(*,*) 'dcircle - x,y :',x,y
! The second lateration may encounter instability when r1 and r2 are small.
!     First test : is it a line source ?
  line=(r1 < epsilon(r1) .and. r2 < epsilon(r2))
!     if r1 OR r2 are null, no trilateration
  if ((r1 < epsilon(r1) .or. r2 < epsilon(r2)).and..not.line) then
     dcircle=infinity
!    if (verbose==2) write(*,*) 'dcircle - r1 or r2 are null'
     return
  endif
!     in case of a line case, dcircle = y if x is in the gate
  if (line) then
     if (x < 0._pr .or. x > d12) then
        dcircle=infinity
        if (verbose==2) write(*,*) 'line case but x is not in the gate'
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
  if (verbose==2) write(*,*) 'dcircle - xc,yc :',xc,yc
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
!    if (verbose==2) write(*,*) 'dcircle - a : ',a
!    if (verbose==2) write(*,*) 'dcircle - a outside range'
     return
  endif
! Finally ....
  dcircle=sqrt((x-xc)**2._pr+(y+yc)**2._pr)  ! remember ... distance to (xc,-yc)
  return
end function dcircle
!###############################################################################
subroutine updclosedist(amesh,nton,distarray,k,i,j)
   use LAT_mesh

   type(mesh) :: amesh
   type(containern), dimension(amesh%Nnodes) :: nton
   real(pr), dimension(:,:), allocatable :: distarray
   integer(pin) :: i,j,k

   if (distarray(amesh%cell(k,i),amesh%cell(k,j)) == infinity) then
           distarray(amesh%cell(k,i),amesh%cell(k,j))=get_donedge(amesh,nton,amesh%cell(k,i),amesh%cell(k,j))
! reciprocity
      distarray(amesh%cell(k,j),amesh%cell(k,i))=distarray(amesh%cell(k,i),amesh%cell(k,j))
   endif
end subroutine updclosedist
!###############################################################################
end module distance
