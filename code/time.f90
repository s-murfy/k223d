! main module to compute first arrival time on a 2D mesh
! using a fast marching trilateration method
! Version January 2026
! Actions :
! - Homogeneization of the code with trilat-distance module
! - remplacement of FVsNodes and EtoE by nton lists
module time
 use generic
 use LAT_mesh
!  use LAT_time
 implicit none

 ! Things to do:
!     Check the behaviour of the code with null velocity cells  
!     
type localxy   ! used by x_entering in distance ! check if it should stay here
   real(pr) :: x,y
end type localxy

type ops
   real(pr) :: d13,d23,d12       ! current edge lengths 
   real(pr) :: d12a,d23a,t2a     ! edge lengths and time for N2a
   real(pr) :: d12b,d23b,t2b     ! edge lengths and time for N2b
   real(pr) :: t1,t2             ! current times
   real(pr) :: k1,k2             ! distances
   real(pr) :: x3,y3,xc,yc       ! deduced positions in the local reference system
   real(pr) :: a                 ! impact point
   real(pr) :: v1,v2,v3          ! velocity, velocity below, velocity "on the left"
end type ops 

contains

!##########################################################
subroutine pre_timeonevsall2d_onvertex(amesh,k,time,nton)
! initialisation of timeonevsall2d in case of a unique source on a vertex k.
! the time array is initialized to infinity except for the source point k which is set to zero.
use LAT_mesh
use lists

  type(mesh) :: amesh
  real(pr), dimension(:), allocatable :: time
  integer(pin) :: k
! the nton array is an hybrid structure of a list of all nodes connected to
! node nton(i). Each list element contains a pointer to the second node
! connected to the current node nton(i) and an array of size 2 with the cell indexes
! on each side of the edge between the two nodes. It carries also the edge length
! in the cell with the smallest velocity.
  type(containern), dimension(:), allocatable :: nton

! allocation of the time array and nton array
  allocate(time(amesh%Nnodes),nton(amesh%Nnodes))
! computing the node-to-node edge array nton
  call compnton(amesh,nton)
! initialisation to infinity
  time=infinity
!   initialization of time at k to zero
  time(k)=0._pr
  return
end subroutine pre_timeonevsall2d_onvertex

!########################################################
subroutine timeonevsall2d(amesh,velocity,time,nton,adiff)
  use LAT_mesh
  use LAT_time
  use lists
! given a mesh (amesh) and velocities (velocity) defined in
! the mesh cells, timeonevsall2d computes the first arrival
! time on the mesh (time) from all the nodes in time array which
! have a finite time (i.e. not infinite)
! the adiff structure should be also contained inside trilat_distance module
! taking care also of the fast switch for diffraction

! inout variables

  type(mesh) :: amesh ! the mesh structure
  type(diff) :: adiff
  real(pr), dimension(amesh%Nnodes) :: time ! the time array 
  real(pr), dimension(amesh%Ncells) :: velocity ! velocity defined on cells declared in LAT_time
  type(containern), dimension(amesh%Nnodes) :: nton ! the node to node array (sparse matrix)

! local variables

  integer(pin) :: i
  logical, dimension(amesh%Nnodes) :: inthelist
  logical :: begin
  type(listn), pointer :: ntodo
  type(listn), pointer :: pcur
  type(liste), pointer :: pcure
  type(liste), pointer :: last
  type(listn), pointer :: pmin
  integer(pin) :: toggle
  logical :: hasbeenupdated
  integer(pin) :: kn
  integer(pin), dimension(2) :: idcell,idnode,cellonedge12
  type(ops) :: tri_ops
  real(pr) :: tedge
  real(pr) :: kedge
  real(pr), dimension(2) :: tface,tplane,thead
  integer(pin) :: lowercell
  logical :: do_tplane,do_tface,do_thead
  real(pr), dimension(2) :: kface
  real(pr) :: khead
  real(pr) :: ttest
  logical :: skip_neighbor
! diffraction variables 
  integer(pin) :: idiff,waitfordiff,diff_counter
  real(pr), dimension(amesh%Nnodes) :: secondary_time
  logical, dimension(amesh%Nnodes) :: checksecondary
  integer(pin), parameter :: waitdiffthres=25
  logical :: reloop,firstrun,diffoccur

! this version uses only the node-to-node list nton

   verbose=0

! initialize the secondary array
  if (.not.allocated(kappa)) allocate(kappa(amesh%Nnodes))
  if (.not.allocated(mode)) allocate(mode(amesh%Nnodes))
  if (size(kappa) /= amesh%Nnodes) then
     deallocate(kappa)
     allocate(kappa(amesh%Nnodes))
  endif
  if (size(mode) /= amesh%Nnodes) then
     deallocate(mode)
     allocate(mode(amesh%Nnodes))
  endif
  if (verbose == 3) then
     if (.not.allocated(sface)) allocate(sface(amesh%Nnodes))
     if (.not.allocated(shead)) allocate(shead(amesh%Nnodes))
     if (.not.allocated(splane)) allocate(splane(amesh%Nnodes))
     if (.not.allocated(sedge)) allocate(sedge(amesh%Nnodes))
     if (.not.allocated(origidnode)) allocate(origidnode(amesh%Nnodes))
     if (size(sface) /= amesh%Nnodes) then
        deallocate(sface)
        allocate(sface(amesh%Nnodes))
     endif
     if (size(shead) /= amesh%Nnodes) then
        deallocate(shead)
        allocate(shead(amesh%Nnodes))
     endif
     if (size(splane) /= amesh%Nnodes) then
        deallocate(splane)
        allocate(splane(amesh%Nnodes))
     endif
     if (size(sedge) /= amesh%Nnodes) then
        deallocate(sedge)
        allocate(sedge(amesh%Nnodes))
     endif
     if (size(origidnode) /= amesh%Nnodes) then
        deallocate(origidnode)
        allocate(origidnode(amesh%Nnodes))
     endif
     sface = infinity
     shead = infinity
     splane = infinity
     sedge = infinity
     origidnode = -1._pr
  endif

  secondary_time = time   !secondary sources
  checksecondary=.true.   !assume all points have been checked
!  checksecondary=.false.  ! assume no points have been checked

! initialization of ntodo from non infinite values inside time
! same as the previous version
  inthelist=.false. ! check the usefulness of the logical array
  begin=.true.

  ! designate nodal points that will be checked as secondary sources
  if (.not.adiff%fast) then
      do i = 1,adiff%NdiffNodes
         checksecondary(adiff%nodes(i))=.false.
      enddo
  endif

  write(*,*) 'starting trilateration.....'

  do i=1,amesh%Nnodes
   if (secondary_time(i) > .9_pr*infinity) cycle
   ! if (time(i) > .9_pr*infinity) cycle
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
     kappa(i) = 0._pr
     mode(i) = 1
     inthelist(i)=.true.
     ntodo%previous=>pcur ! linking the current list position to the first element ntodo
  enddo

   if (verbose == 2) call printlist(ntodo)
   if (verbose == 1) write(*,*) 'entering in the ntodo list management loop'
   call flush()
! initializing reloop to true
! The condition is changed or not at the end of the loop
! The fast first loop is mandatory with a null distance offset
   reloop=.true.
   firstrun=.true.
   diffoccur=.false.
   do while (reloop)  ! I loop over diffraction points 

! loop on the to do list 
    do while (associated(ntodo)) ! II loop over the nodes in the ntodo list
! searching the node with the minimum time in the ntodo list : pmin
!      if (verbose == 2) call printlist(ntodo)
       if (verbose == 1) call printlist(ntodo)
       call lookformin(ntodo,secondary_time,amesh%Nnodes,pmin)

       if (verbose == 1) write(*,*) 'propagating time from node #',pmin%idnode-1
       if (verbose == 2) write(*,*) 'its own time is :',time(pmin%idnode)

!    initilisation of the sweep loop
      hasbeenupdated=.true.
      toggle=2
      do while (hasbeenupdated) ! III loop on sweeps from pmin
        hasbeenupdated=.false.
        toggle=3-toggle   !    toggle variable swap state between 1 and 2 at each sweep
        if (verbose == 1) then
          if (toggle==1) write(*,*) 'sweep forward'
          if (toggle==2) write(*,*) 'sweep backward'
        endif
! tuning back to a node list reading nton (node to node)
! this time version uses a node-to-node list
        if (toggle==1) then
           pcure=>nton(pmin%idnode)%ptr
        else
           pcure=>last
        endif
! going through all the nodes connected to pmin%idnode - pcure is a pointer to a node
        do while (associated(pcure)) ! IV loop on nodes connected to pmin%idnode  
!
! pmin%idnode is the node called 1
! pcure%idnode is the node called 3
!
! the -1 is to shift the fortran index to the paraview index convention           
           if (verbose == 1) write(*,*) '####################################################'
           if (verbose == 1) write(*,*) 'working to extent node ',pmin%idnode-1,' to ', pcure%idnode-1
           kn=pcure%idnode ! kn is the current node to update
           skip_neighbor = .false.
           if (verbose == 1) write(*,*) 'considering node ',kn-1
! cell index and edge length are stored on the canonical edge record
! integer(pin), dimension(2) :: cellonedge
!          definition of the cell indexes and node indexes on both side of the edge
           idcell=pcure%ref%cellonedge ! array of size 2
! idnode are the third nodes of each cell
           idnode(1)=sum(amesh%cell(idcell(1),:))-kn-pmin%idnode ! sum of the node indexes of the cell minus the 2 known indexes
           if (idcell(2) == 0) then
              idnode(2)=0
           else
              idnode(2)=sum(amesh%cell(idcell(2),:))-kn-pmin%idnode
           endif
!
! idnode(1) is called 2a and idnode(2) is called 2b
!    
!     the axis system is traslated and rotated each time with pmin%idnode at the origin (0,0),
!     idnode(i) on the x-axis and kn in the semispace (y>0)

!     common values before the inner loop on the two side cells
           tri_ops%d13=pcure%ref%donedge
           tri_ops%t1=secondary_time(pmin%idnode)
           tri_ops%k1=kappa(pmin%idnode)

!     first operator : edge propagation between pmin%idnode and kn at max(velocity(idcell(:)))
           if (idcell(2) /= 0) then
              tedge=tri_ops%t1+tri_ops%d13/max(velocity(idcell(1)),velocity(idcell(2)))
           else
              tedge=tri_ops%t1+tri_ops%d13/velocity(idcell(1))
           endif
           kedge=min(tri_ops%k1+tri_ops%d13,infinity)
           if (verbose == 2) write(*,*) 'tedge : ',tedge

!     lists of the operator to consider :
!      - curved 2 tface(2) and two thead(2)
!      - planar 2 tplane(2)
          thead = infinity 
          tface = infinity
          tplane = infinity 

!     loop on the two side cells
           tri_ops%t2a = infinity 
           tri_ops%t2b = infinity

           do i=1,2 ! V loop on side cells
              if (verbose == 2) write(*,*) "cell number : ",i,idcell(i)-1
              if (idcell(i) == 0) cycle ! may only happen for i=2 - The cell is on the border of the mesh
! using the call to get_donedge which also retrieves the two cells on each side of the edge 12
              tri_ops%d12=get_donedge(amesh,nton,pmin%idnode,idnode(i),cellonedge12)
              tri_ops%d23=get_donedge(amesh,nton,idnode(i),kn)
! the lower cell is the other cell respect to idcell(i)
              lowercell=sum(cellonedge12)-idcell(i) ! may be zero if there is no cell beneath - this case is treated below
              tri_ops%v1=velocity(idcell(i))

!             computation of the position of kn in the local system through trilateration
              tri_ops%x3 = (tri_ops%d12**2._pr-tri_ops%d23**2._pr+tri_ops%d13**2._pr)/(2._pr*tri_ops%d12)
              tri_ops%y3 = sqrt(max(tri_ops%d13**2._pr-tri_ops%x3**2._pr,0._pr))

!             do we compute this operator on this cell ?
              do_tplane = .true.
              do_tface = .true.
              do_thead = .true. 
            !   it_is_a_conic=.false.

! cell estimates only if idnode(i) has a finite time
            !   if (time(idnode(i)) == infinity) cycle !!!!! check it
              if (abs(secondary_time(idnode(i))-infinity) <= water_level(infinity)) cycle !!!!! check it
              tri_ops%t2 =secondary_time(idnode(i))
              tri_ops%k2 = kappa(idnode(i))

! thead only if idcell(3-i) /= 0
              if ( idcell(3-i) == 0 ) then ! condizione valida solo per i=1. 
                  do_thead = .false.
               else
                  tri_ops%v3 = velocity(idcell(3-i))
                  if (abs(tri_ops%v1-tri_ops%v3) > epsilon(1.)) then 
                     kappa(kn) = infinity 
                  endif
               endif 

! no cell beneath idcell(i)
               if (lowercell == 0) then 
                  cycle
               else
                  tri_ops%v2 = velocity(lowercell)
               endif 

! One node is attached to a cell with a planar propagation
               if ((abs(tri_ops%k1-infinity) <= water_level(infinity)).or.(abs(tri_ops%k2-infinity) <= water_level(infinity))) then ! check it
                  do_tface = .false. 
                  do_thead = .false. 
               endif
   
! Estimation of thead only if velocity(idcell(3-i))>velocity(idcell(i))
               if (do_thead)  do_thead = (abs(tri_ops%v1-tri_ops%v3) > water_level(tri_ops%v1))  

! Wavefront has crossed from one velocity layer to another ! TO BE CHECKED
               ! if (tri_ops%k1 <= epsilon(1.).or.tri_ops%k2 <= epsilon(1.)) then   ! this use to be used 
               ! write(*,*)  abs(tri_ops%v1-tri_ops%v2)/tri_ops%v2 
               ! if ((abs(tri_ops%v1-tri_ops%v2)/tri_ops%v2 > 0.1_pr)) then 
               if (abs(tri_ops%v1-tri_ops%v2) > epsilon(1.)) then 
                  do_tface = .false.
               endif

                        
! There is a conic propagation in the lower cell
               if (tri_ops%v1 < tri_ops%v2 ) then 
                    if (abs(abs(tri_ops%t2-tri_ops%t1)*tri_ops%v2-tri_ops%d12) <= water_level(tri_ops%d12)) then 
                     do_tface = .false.
                     do_thead = .false. 
                     do_tplane = .true. 
                  endif
               endif

! We can compute a curved operator but the circles do not cross !
               if (do_thead.or.do_tface)  then  
                  if (.not.calc_trilat(tri_ops,verbose) )  then 
                     do_tface = .false.
                     do_thead = .false.
                  endif
               endif

! curved operators                
               if (do_tface) call calc_tcircle(tri_ops,kface(i),tface(i),verbose)
               if (verbose == 2) write(*,*) 'tface : ',do_tface,tface(i)

! head wave operator
               if (do_thead) call calc_thead(tri_ops,khead,thead(i))

! planar operator                
               if (do_tplane) call calc_tplane(tri_ops,tplane(i),kn,pmin%idnode)

               ! if (tplane /= infinity .or. tface /= infinity) tedge = infinity      
               if (i == 1) then    
                  tri_ops%d12a = tri_ops%d12
                  tri_ops%d23a = tri_ops%d23
                  tri_ops%t2a = tri_ops%t2
               else 
                  tri_ops%d12b = tri_ops%d12
                  tri_ops%d23b = tri_ops%d23
                  tri_ops%t2b = tri_ops%t2
               endif    

           enddo ! V loop on side cells

!       best estimate for kn
          ttest=min(tedge,tplane(1),tplane(2),tface(1),tface(2),thead(1),thead(2))


           if (verbose == 2) then 
            write(*,*) 'ttest : ',ttest
            write(*,*) 'edge:   ',tedge
            write(*,*) 'face: ', tface(1),tface(2)
            write(*,*) 'head: ',thead(1),thead(2)
           endif

!       WORSE move to the next node attached to pmin%inode           
           if (ttest > secondary_time(kn)) skip_neighbor = .true.

! save all estimates 
         if ((.not.skip_neighbor).and.(verbose == 3)) then 
           sface(kn) = min(tface(1),tface(2),sface(kn))
           shead(kn) = min(thead(1),thead(2),shead(kn) )
           splane(kn) = min(tplane(1),tplane(2),splane(kn))
           sedge(kn)  = min(tedge,sedge(kn))
           origidnode(kn) = pmin%idnode-1
         endif 
!       order of preferences : tface(),thead,tedge,tplane()
!       modes :                  1       2     4      5     

!    if the result is equal to the time in kn with a smaller mode, do an update of mode and kappa

         if (.not.skip_neighbor) then
         if (abs(ttest-tface(1)) <= water_level(ttest)) then
            if (.not.better_face(ttest,tface(1),kface(1),secondary_time(kn),mode(kn),kappa(kn))) skip_neighbor = .true.
         elseif (abs(ttest-tface(2)) <= water_level(ttest)) then
            if (.not.better_face(ttest,tface(2),kface(2),secondary_time(kn),mode(kn),kappa(kn))) skip_neighbor = .true.
!    thead is better (mode=2) or equal but with a better mode
         elseif (abs(ttest-thead(1)) <= water_level(ttest)) then
            if (.not.better_head(ttest,thead(1),infinity,secondary_time(kn),mode(kn),kappa(kn))) skip_neighbor = .true.
         elseif (abs(ttest-thead(2)) <= water_level(ttest)) then
            if (.not.better_head(ttest,thead(2),infinity,secondary_time(kn),mode(kn),kappa(kn))) skip_neighbor = .true.
!    tedge is better (mode=4) or equal but with a better mode
         elseif (abs(ttest-tedge) <= water_level(ttest)) then
            if (.not.better_edge(ttest,tedge,kedge,secondary_time(kn),mode(kn),kappa(kn))) skip_neighbor = .true.
!    tplane is the best (mode=5)
         elseif (abs(ttest-tplane(1)) <= water_level(ttest)) then
            if (.not.better_plane(ttest,tplane(1),secondary_time(kn),mode(kn),kappa(kn))) skip_neighbor = .true.
         elseif (abs(ttest-tplane(2)) <= water_level(ttest)) then
            if (.not.better_plane(ttest,tplane(2),secondary_time(kn),mode(kn),kappa(kn))) skip_neighbor = .true.
         endif
         endif
         if (skip_neighbor) then
            if (toggle==1) then
               last=>pcure
               pcure=>pcure%next
            else
               pcure=>pcure%previous
            endif
            cycle
         endif
!       if we are here, we have updated kn
        hasbeenupdated=.true.

!       update if needed ntodo list with kn
!       
        if (.not.inthelist(kn)) then
           if (verbose == 2) write(*,*) ' adding ',kn-1,' to the list', mode(kn),secondary_time(kn)
           inthelist(kn)=.true.
           allocate(ntodo%previous%next)
           nullify(ntodo%previous%next%next)
           ntodo%previous%next%previous=>ntodo%previous
           ntodo%previous%next%idnode=kn
           ntodo%previous=>ntodo%previous%next
        endif
! go to the next node connected to pmin%idnode
           if (toggle==1) then
              last=>pcure
              pcure=>pcure%next
           else
              pcure=>pcure%previous
           endif
      enddo ! IV loop on nodes connected to pmin%idnode      
   enddo ! III loop on sweeps from pmin

! diff case and no diffraction has been detected yet
   if (.not.firstrun.and..not.diffoccur) then

      if (secondary_time(pmin%idnode) < time(pmin%idnode)) then
         ! if (secondary_time(pmin%idnode)+time(idiff) < time(pmin%idnode)) then
!   diffraction detected
         write(*,*) 'diffraction detected '
         diffoccur=.true.
       else
!   no diffraction
         waitfordiff=waitfordiff-1   ! check to see if there is an improvement in cells near diff point
         ! write(*,*) waitfordiff
         if (waitfordiff == 0) then
!   countdown expired for diffraction. hard exit on ntodo list to pass to another
!                                      potential secondary source
            call wipe(ntodo)
            exit
       endif
      endif
   endif
! remove pmin element from the to do list
    inthelist(pmin%idnode)=.false.
    if (verbose == 2) write(*,*) "removing ",pmin%idnode-1," from the list"
    call removepminfromthelist(pmin,ntodo)
    if (verbose == 2) read(*,*)

   enddo ! associated(ntodo)

! decision to reloop or not
! fast case
    if (adiff%fast) then
      time = secondary_time
      if (verbose == 1) write(*,*) 'exiting based on fast protocol ....'
! first case of general exit : only one run is requested
      reloop=.false.
      cycle ! to the reloop
    endif
    if (firstrun) then
         time = secondary_time
         ! diff_counter = 0   ! counter for number of diffraction points
         write(*,*) 'first traveltime computation '
    else
         write(*,*) 'next traveltime computation '

        if (diffoccur) then
! a secondary diff has occured
          do i=1,amesh%Nnodes
             time(i)=min(secondary_time(i),time(i))
            !  time(i)=min(secondary_time(i)+time(idiff),time(i))
          enddo
        endif
     endif

   idiff=nextdiffid(time,amesh%Nnodes,checksecondary)
! second case of general exit : all the choosen secondary source have been investigated.
   !  diff_counter = diff_counter+1
   !  if (diff_counter > adiff%NdiffNodes) then
   !    reloop=.false.
   !    cycle
   ! endif
   ! idiff = adiff%nodes(diff_counter)
   write(*,*) 'diffraction points considered    ',idiff 
!  second case of general exit : all the potential secondary source have been investigated.
    if (idiff == 0) then
        reloop=.false.
        cycle
     endif
! preparation for the next loop
    checksecondary(idiff)=.true.
    waitfordiff=waitdiffthres    ! reset counter for diffraction
    diffoccur=.false.            ! reset diffraction check
! At this point notodo is deallocated. reallocate it with idiff
    allocate(ntodo)
    ntodo%previous=>ntodo
    nullify(ntodo%next)
    ntodo%idnode=idiff
    inthelist=.false.
    inthelist(idiff)=.true.
! reinitialize time array
    secondary_time=infinity
    kappa(idiff) = 0._pr
    mode(idiff) = 1
    secondary_time(idiff) = time(idiff) !0._pr
    firstrun=.false.
   enddo ! end reloop section

  return 
end subroutine timeonevsall2d
!###############################################################################
function better_face(ttest,tface,kface,time,mode,kappa)
! check if faces produce fastest time 
real(pr) :: ttest,tface,kface,time,kappa
integer(pin) :: mode
logical :: better_face

better_face=.true.
   !       equal but with a better mode
   if (abs(tface-time) <= water_level(ttest)) then
      if (mode > 1) then  ! better mode
         mode = 1
         kappa= kface
      else 
         better_face=.false.
         ! exit
      endif
   !       better
   else
      mode = 1
      kappa = kface
      time = tface
   endif


return 
end function better_face
!###############################################################################
function better_head(ttest,thead,khead,time,mode,kappa)
   ! check if the head wave produces fastest time 
real(pr) :: ttest,thead,khead,time,kappa
integer(pin) :: mode
logical :: better_head

!       thead is the best (mode=2)
better_head=.true.
!       equal but with a better mode
   if (abs(thead-time) <= water_level(ttest)) then
      if (mode >2) then
         mode = 2
         kappa = khead
      else
         ! cycle
         better_head = .false.
      endif
!       better
   else
      mode = 2
      kappa = khead
      time = thead
   endif
return
end function better_head 
!###############################################################################
function better_edge(ttest,tedge,kedge,time,mode,kappa)
   ! check if edge produces fastest time 
real(pr) :: ttest,tedge,kedge,time,kappa
integer(pin) :: mode
logical :: better_edge
!       tedge is the best (mode=4)
better_edge=.true.
!       equal but with a better mode
   if (abs(tedge-time) <= water_level(ttest)) then
      if (mode > 4) then
         mode = 4
         kappa = kedge
      else
         better_edge=.false.
      endif
!       better
   else
      mode = 4
      kappa = kedge
      time = tedge
   endif

return
end function better_edge
!###############################################################################
function better_plane(ttest,tplane,time,mode,kappa)
   ! check if plane wave produce fastest time 
real(pr) :: tplane,ttest,kappa,time
integer(pin) :: mode
logical :: better_plane

better_plane=.true.
!       equal but with a better mode
      if (abs(tplane-time) <= water_level(ttest)) then
         better_plane=.false.
!       better
       else
         mode = 5
         kappa = infinity
         time = tplane
      endif

end function better_plane
!###############################################################################
function x_entering(po1,po2)
! crossing the y=0 starting from po2 (y < 0) toward po1 (y > 0)
  real(pr) :: x_entering
  type(localxy) :: po1,po2

  x_entering=abs(po2%x-po1%x)/(1._pr+(abs(po2%y)/po1%y))
  if (po1%x < po2%x) then
     x_entering=po1%x+x_entering
  else
     x_entering=po1%x-x_entering
  endif
  return
end function x_entering
!###############################################################################
subroutine calc_tonedge(tri_ops,kappa,t3)
   type(ops) :: tri_ops
   real(pr) :: t3,kappa

      t3 = tri_ops%t1+tri_ops%d13/tri_ops%v1

      kappa = tri_ops%k1+tri_ops%d13
      
   return  
end subroutine calc_tonedge
!###############################################################################
function calc_trilat(tri_ops,verbose)
   type(ops) :: tri_ops
   logical :: calc_trilat
   integer(pin) :: verbose

   ! function calc_trilat(d12,d13,d23,t1,t2,k1,k2,v,x,y,xc,yc,a)
   ! computes the coordinates in 2D plane of three vertices V1,V2,V3 making the
   ! assumptions that V1 is the origin (0,0), V2 is on the x axis (0,d12). The
   ! first part gives the coordinates x,y of V3.
   ! Same set of equations are used to estimate the origin of a point distant by
   ! r1 and r2 from V1 and V2 respectively (xc,+/-yc). The function returns the
   ! distance between (xc,-yc) and V0 (x,y)
   ! The position of V3 and C (the virtual origin) in the local reference system is returned
   ! with po1 and po2 respectively
  
   real(pr) :: yc_sqrt
   type(localxy) :: po1,po2
 
   calc_trilat = .true.


   po1%x=tri_ops%x3
   po1%y=tri_ops%y3
 ! if (verbose==3) write(*,*) 'tcircle - x,y :',x,y

   tri_ops%xc=(tri_ops%d12**2._pr-tri_ops%k2**2._pr+tri_ops%k1**2._pr)/(2._pr*tri_ops%d12)
   yc_sqrt = tri_ops%k1**2._pr-tri_ops%xc**2._pr
   if (yc_sqrt < 0._pr) then 
      calc_trilat = .false.
      return
   endif
   tri_ops%yc = -sqrt(yc_sqrt)  ! minus required as this is the virtual source located below the triangle 

 ! if (verbose==3) write(*,*) 'tcircle - xc,yc :',xc,yc 
   po2%x = tri_ops%xc
   po2%y = tri_ops%yc
 ! if (verbose==3) write(*,*) 'tcircle - xc,yc :',xc,yc
   tri_ops%a = x_entering(po1,po2)

   return 
end function calc_trilat
!###############################################################################
subroutine calc_tcircle(tri_ops,r3,t3,verbose)
   type(ops) :: tri_ops
   real(pr) :: r3,t3
   integer(pin) :: verbose

   !use type ops
   ! d13,d23,d12 ! current (?) edge lengths
   ! d12a,d23a,t2a ! value on the first cell
   ! d12b,d23b,t2b ! value on the second cell
   ! t1,t2 ! current tima at 1 and 2
   ! k1,k2,r1,r2 ! curvature and 
   ! x3,y3,xc,yc ! position of the 3rd point
   ! a1,a2
   ! v1,v2,v3 ! v1 in the cell, v2 in the lowercell
   ! va,vb

   ! subroutine calc_tcircle(d12,d13,d23,t1,t2,x,y,xc,yc,a,k1,k2,v,r3,t3,lmode)
   ! computes the coordinates in 2D plane of three vertices V1,V2,V3 making the
   ! assumptions that V1 is the origin (0,0), V2 is on the x axis (0,d12). The
   ! first part gives the coordinates x,y of V3.
   ! Same set of equations are used to estimate the origin of a point distant by
   ! r1 and r2 from V1 and V2 respectively (xc,+/-yc). The function returns the
   ! distance between (xc,-yc) and V0 (x,y)
   ! The position of V3 and C (the virtual origin) in the local reference system is returned
   ! with po1 and po2 respectively
     real(pr) :: d_c_a,d_a_3,ta

     t3 = infinity 
      if (tri_ops%a < -water_level(tri_ops%a) .or. tri_ops%a-tri_ops%d12 > water_level(tri_ops%a)) then
         t3=infinity 
   !    if (verbose==3) write(*,*) 'dcircle - a : ',a
   !    if (verbose==3) write(*,*) 'dcircle - a outside range'
        return
     endif

     d_c_a = sqrt((tri_ops%a-tri_ops%xc)**2._pr+tri_ops%yc**2._pr)
     d_a_3 = sqrt((tri_ops%a-tri_ops%x3)**2._pr+tri_ops%y3**2._pr)
     r3 = d_c_a + d_a_3

     if (d_c_a-tri_ops%k1 < 0.) then
        ta = tri_ops%t1+(d_c_a-tri_ops%k1)/tri_ops%v1
     else
        ta = tri_ops%t1+(d_c_a-tri_ops%k1)/tri_ops%v2
     endif
     t3 = ta + d_a_3/tri_ops%v1

     if (verbose==3) then
        write(*,*) "d_c_a,tri_ops%k1,d_c_a-tri_ops%k1",d_c_a,tri_ops%k1,d_c_a-tri_ops%k1
        write(*,*) " tri_ops%t1,(r3-tri_ops%k1),t3",tri_ops%t1,(r3-tri_ops%k1),t3
        write(*,*) "tri_ops%t1,tri_ops%t2",tri_ops%t1,tri_ops%t2
     endif

     return
   end subroutine calc_tcircle
!###############################################################################
subroutine calc_thead(tri_ops,r3,t3)
   type(ops) :: tri_ops
   real(pr) :: r3, t3

   ! computes the coordinates in 2D plane of three vertices V1,V2,V3 making the
   ! assumptions that V1 is the origin (0,0), V2 is on the x axis (0,d12). The
   ! first part gives the coordinates x,y of V3.
   ! Same set of equations are used to estimate the origin of a point distant by
   ! r1 and r2 from V1 and V2 respectively (xc,+/-yc). The function returns the
   ! distance between (xc,-yc) and V0 (x,y)
    
   ! e31 is the unit vector pointing from V3 to V1 
   ! e1c is the unit vector pointing from V1 to (xc,-yc)
   ! e3c is the unit vector pointing from V3 to (xc,-yc)

   real(pr) :: r3c,rnc,rn3
   real(pr) :: r_n_a2,r_c1_a2,ta,d_c_a
   real(pr) :: cos_c1m,cos_cnm,cos_c31
    
   t3 = infinity

   if ((abs(tri_ops%xc) <= epsilon(1._pr)).and.(tri_ops%yc <= epsilon(1._pr))) return 
   ! cos_c1m = (-x*xc-y*yc)/sqrt(xc**2._pr+yc**2._pr)/sqrt(x**2._pr+y**2._pr)
   cos_c1m = (-tri_ops%x3*tri_ops%xc-tri_ops%y3*tri_ops%yc)     &
              /sqrt(tri_ops%xc**2._pr+tri_ops%yc**2._pr)      & ! Note: this line produces a division by zero     
              /sqrt(tri_ops%x3**2._pr+tri_ops%y3**2._pr)

   if (cos_c1m < 0._pr) return 
   cos_cnm = tri_ops%v1/tri_ops%v3
 
   if (cos_c1m > cos_cnm) return
   
   r3c = sqrt((tri_ops%xc-tri_ops%x3)**2._pr+(tri_ops%yc-tri_ops%y3)**2._pr)
   ! r3c = sqrt((xc-x)**2._pr+(yc-y)**2._pr)
   cos_c31 = (-tri_ops%x3*(tri_ops%xc-tri_ops%x3)-tri_ops%y3*(tri_ops%yc-tri_ops%y3))   &
                  /r3c/sqrt(tri_ops%x3**2._pr+tri_ops%y3**2._pr)
   ! cos_c31 = (-x*(xc-x)-y*(yc-y))/r3c/sqrt(x**2._pr+y**2._pr)
   
   
   if (cos_cnm > cos_c31) return 
   
   rnc = r3c*sqrt((1._pr-cos_c31**2._pr)/(1._pr-cos_cnm**2._pr))
   rn3 = r3c*cos_c31-rnc*cos_cnm
   ! t3 = rnc/tri_ops%v1+rn3/tri_ops%v3    !

   d_c_a = sqrt((tri_ops%a-tri_ops%xc)**2._pr+tri_ops%yc**2._pr)
   ta = tri_ops%t1+(d_c_a-tri_ops%k1)/tri_ops%v2

   r_c1_a2 = sqrt((tri_ops%a-tri_ops%xc)**2._pr+tri_ops%yc**2._pr)
   r_n_a2 = rnc - r_c1_a2
   t3 = ta + r_n_a2/tri_ops%v1 +rn3/tri_ops%v3    !
   r3 = rnc+rn3 !rn3     ! kappa value 


return 
end subroutine calc_thead   
!###############################################################################
!###############################################################################
function tonedge(amesh,i,j,v)
use LAT_mesh

  type(mesh) :: amesh
  integer(pin) :: i,j
  real(pr) :: tonedge,v

  tonedge=sqrt((amesh%px(i)-amesh%px(j))**2._pr+&
               (amesh%py(i)-amesh%py(j))**2._pr+&
               (amesh%pz(i)-amesh%pz(j))**2._pr)/v
  return
end function tonedge
!###############################################################################
subroutine rotation(x,y,theta,xr,yr)
   real(pr),intent(in) :: x,y,theta
   real(pr), intent(out) :: xr,yr

   xr = x*cos(theta)-y*sin(theta)
   yr = x*sin(theta)+y*cos(theta)

end subroutine rotation
!###############################################################################
subroutine calc_tplane(tri_ops,t3,kn,inode)
   ! function calc_tplane(d12,d13,d23,t1,t2,v)
   type(ops) :: tri_ops
   real(pr) :: t3
   real(pr) :: ta,xd,dt
   integer(pin) :: kn,inode

   t3 = infinity
   !
   ! tri_ops%x3 = (tri_ops%d12**2._pr-tri_ops%d23**2._pr+tri_ops%d13**2._pr)/(2._pr*tri_ops%d12)
   ! tri_ops%y3 = sqrt(max(tri_ops%d13**2._pr-tri_ops%x3**2._pr,0._pr))

   ! xb = d12 
   ! xc = (d12**2._pr-d23**2._pr+d13**2._pr)/(2._pr*d12)
   ! yc = sqrt(d13**2._pr-xc**2._pr)
   ! txb = t1+xc*(t2-t1)/xb 
   ta = tri_ops%t1+tri_ops%x3*(tri_ops%t2-tri_ops%t1)/tri_ops%d12 
   xd = tri_ops%d12-tri_ops%x3
   dt = tri_ops%t2 - tri_ops%t1
   
   ! if ((kn == 6392).and.(inode == 295))  then 
   !    write(*,*) 'do tplane for 6392 using inode 295'
   !    write(*,*) dt-tri_ops%d12*tri_ops%x3/(tri_ops%v1*sqrt(tri_ops%x3**2._pr+tri_ops%y3**2._pr))
   !    write(*,*) dt+tri_ops%d12*xd/(tri_ops%v1*sqrt(xd**2._pr+(tri_ops%y3)**2._pr))
   ! endif

   if ( ( dt-tri_ops%d12*tri_ops%x3/(tri_ops%v1*sqrt(tri_ops%x3**2._pr+tri_ops%y3**2._pr)) <= water_level(1._pr)  ).and.       &
         ( dt+tri_ops%d12*xd/(tri_ops%v1*sqrt(xd**2._pr+(tri_ops%y3)**2._pr)) >= water_level(1._pr)  ) ) then 
   ! if ( ( dt-tri_ops%d12*tri_ops%x3/(tri_ops%v1*sqrt(tri_ops%x3**2._pr+tri_ops%y3**2._pr)) <= epsilon(1._pr)  ).and.       &
   !          ( dt+tri_ops%d12*xd/(tri_ops%v1*sqrt(xd**2._pr+(tri_ops%y3)**2._pr)) >= -epsilon(1._pr) )) then 
               if ((kn == 6392).and.(inode == 295))  then 
                  write(*,*) 'pass the test !'
               endif

              ! wave propagating from left to right
         t3 = ta + tri_ops%y3*sqrt( (1._pr/tri_ops%v1)**2._pr - (dt/tri_ops%d12)**2._pr ) ! 
   endif
   ! if ( ( t2 <= t1+xb*xc/(v*sqrt(xc**2._pr+yc**2._pr))  ).and.       &
   ! (t2 >= t1-xb*(xb-xc)/(v*sqrt((xb-xc)**2._pr+(yc)**2._pr)) ) ) then      ! wave propagating from left to right
   !    calc_tplane = txb + sqrt( (yc/v)**2._pr - (yc*(t1-t2)/xb)**2._pr ) ! 
   ! endif
return 
end subroutine calc_tplane
!###############################################################################
! function tcircle(d12,d13,d23,t1,t2,v)
! ! computes the coordinates in 2D plane of three vertices V1,V2,V3 making the
! ! assumptions that V1 is the origin (0,0), V2 is on the x axis (0,d12). The
! ! first part gives the coordinates x,y of V3.
! ! Same set of equations are used to estimate the origin of a point distant by
! ! r1 and r2 from V1 and V2 respectively (xc,+/-yc). The function returns the
! ! distance between (xc,-yc) and V0 (x,y)
! ! The position of V3 and C (the virtual origin) in the local reference system is returned
! ! with po1 and po2 respectively

!   real(pr) :: tcircle,d12,d13,d23,t1,t2,v
!   real(pr) :: x,y,xc,yc,a,r1,r2
!   type(localxy) :: po1,po2

!   x=(d12**2._pr-d23**2._pr+d13**2._pr)/(2._pr*d12)
!   y=sqrt(max(d13**2._pr-x**2._pr,0._pr))
!   po1%x=x
!   po1%y=y
! ! if (verbose) write(*,*) 'tcircle - x,y :',x,y
!   r1=t1*v
!   r2=t2*v
!   xc=(d12**2._pr-r2**2._pr+r1**2._pr)/(2._pr*d12)
!   yc=sqrt(max(r1**2._pr-xc**2._pr,0._pr))
!   po2%x = xc
!   po2%y = -yc
! ! if (verbose) write(*,*) 'tcircle - xc,yc :',xc,yc
!   a=x_entering(po1,po2)
!   if (a < 0._pr .or. a > d12) then
!      tcircle=infinity
! !    if (verbose) write(*,*) 'dcircle - a : ',a
! !    if (verbose) write(*,*) 'dcircle - a outside range'
!      return
!   endif
!   tcircle=sqrt((x-xc)**2._pr+(y+yc)**2._pr)/v  ! remember ... distance to (xc,-yc)
!   return
! end function tcircle
!###############################################################################
function find_face_cell(cell_no,n1,n2,n,Edge)
   integer(pin) ::find_face_cell,cell_no,n1,n2,n,i
   integer(pin) , dimension(n,4) :: Edge
 
   find_face_cell = -1
   do i = 1,n
       if(Edge(i,3)==n1.and.Edge(i,4)==n2.or.           &
          Edge(i,3)==n2.and.Edge(i,4)==n1) then
              if(Edge(i,1)==cell_no) then
                 find_face_cell = Edge(i,2)
              elseif(Edge(i,2)==cell_no) then
                 find_face_cell = Edge(i,1)
           endif
       endif
   enddo
 
   if (find_face_cell == -1) then
      write(0,*) 'ERROR in find_face_cell: no face found'
      stop
   endif
   return
 end function find_face_cell
!###############################################################################
 function find_face_cell_v2(cell_no,n1,n2,n,FVsNodes)
   integer(pin) ::find_face_cell_v2,cell_no,n1,n2,n,i
   integer(pin) , dimension(n,n,2) :: FVsNodes
   integer(pin), dimension(2) :: candidate_cells
 
   find_face_cell_v2 = -1
   candidate_cells = FVsNodes(n1,n2,:)

   if (candidate_cells(1) == cell_no ) then 
      find_face_cell_v2 = candidate_cells(2)
   elseif (candidate_cells(2) == cell_no ) then 
      find_face_cell_v2 = candidate_cells(1)
   endif
   ! do i = 1,n
   !     if(Edge(i,3)==n1.and.Edge(i,4)==n2.or.           &
   !        Edge(i,3)==n2.and.Edge(i,4)==n1) then
   !            if(Edge(i,1)==cell_no) then
   !               find_face_cell = Edge(i,2)
   !            elseif(Edge(i,2)==cell_no) then
   !               find_face_cell = Edge(i,1)
   !         endif
   !     endif
   ! enddo
 
   if (find_face_cell_v2 == -1) then
      write(0,*) 'ERROR in find_face_cell: no face found'
      stop
   endif
   return
 end function find_face_cell_v2
 !###############################################################################
subroutine calc_theory(amesh,ttheory)
   type(mesh) :: amesh
   real(pr),allocatable, dimension(:) ::ttheory
   real(pr) :: x,y,layer_depth,sx,sy,v1,v2,x1,x2,ddx,theta_c
   real(pr) :: x_cross, tmin,t_test,conic_time
   integer(pin) :: i
   logical :: not_found

   layer_depth = 3000._pr 
   sx = 5000._pr
   sy = 1500._pr
   v1 = 1000._pr
   v2 = 3000._pr
   theta_c = asin(v1/v2)
   ddx = 1._pr

   allocate(ttheory(amesh%Nnodes))
   ttheory = 0._pr
   do i=1,amesh%Nnodes
      x = amesh%px(i)
      y = amesh%py(i)
   ! direct wave 
      if (y <= layer_depth) then
         ttheory(i) = sqrt((x-sx)**2._pr+(y-sy)**2._pr)/v1
         x1 = (layer_depth-sy)*tan(theta_c)
         x2 = (layer_depth-y)*tan(theta_c)   
         if (abs(x-sx) > x1+x2) then 
            conic_time = (2._pr*layer_depth-sy-y)/v1/cos(theta_c)+(abs(x-sx)-x1-x2)/v2
            if (conic_time < ttheory(i)) then 
               ttheory(i) = conic_time 
            endif
         endif
      else
   ! refracted wave 
         x_cross = 0._pr
         tmin = sqrt((x-sx)**2+(y-layer_depth)**2)/v2+(layer_depth-sy)/v1
         not_found = .true.
         do while (not_found)
            x_cross = x_cross+ddx
            t_test = sqrt((abs(x-sx)-x_cross)**2+(y-layer_depth)**2)/v2+ sqrt(x_cross**2._pr+(layer_depth-sy)**2._pr)/v1
            if (t_test < tmin) then 
                 tmin = t_test
            else
                not_found = .false.
            endif
         end do
         ttheory(i) = tmin
      endif
   enddo 
return 
end subroutine calc_theory 
!#########################################################
function nextdiffid(time,n,checksecondary)
! find the minimum value in dist skipping the vertex checked yet and returns
! its position (id)
  integer(pin) :: nextdiffid,n
  real(pr), dimension(n) :: time
  logical, dimension(n) :: checksecondary

  integer(pin) :: i
  real(pr) :: vmin

  nextdiffid=0
  vmin=infinity
  do i=1,n
     if (checksecondary(i)) cycle
     if (time(i) < vmin) then
        vmin=time(i)
        nextdiffid=i
     endif
  enddo
  return
end function nextdiffid
 !###############################################################################
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
end module time
