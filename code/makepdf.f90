module makepdf
use LAT_source
! use distance 
! use lateration
implicit none

contains

!#######################################################################################
subroutine set_quake_pdf(pdf,quake_pdf,quake,defined_area)
  type(pdfinputs) :: pdf
  type(source) :: quake
  ! type(mesh) :: amesh
  logical :: defined_area
  real(pr), allocatable,dimension(:) :: quake_pdf
  integer(pin) :: i
  real(pr) :: ug, pdfint

  allocate(quake_pdf(quake%QuakeElemNo))
  quake_pdf = pdf%g_pdf

return
end subroutine set_quake_pdf
!#######################################################################################
subroutine set_pdf(pdf,quake,target_area)
  use distance 
  ! subroutine set_pdf(pdf,amesh,target_area)
! gaussian parameters
  type(pdfinputs) :: pdf
! main structure
  type(source) :: quake
  ! type(mesh) :: amesh
  integer(pin) :: i
  real(pr), intent(in) :: target_area
  real(pr) :: ran


select case(pdf%pdf_type)
  case('uniform')
! uniform distribution
  case('gauss')
! gaussian distribution
          call random_number(ran)
          i = 1!+int(ran*3.) ! randomly choose between 1-4 gaussian functions

          pdf%ng=i
          allocate(pdf%vertexcenter(pdf%ng),pdf%sizeandhigh(pdf%ng,2))
          do i = 1,pdf%ng
          ! randomly place centres of gaussians with the condtion that they are at least 5km from boundary
              pdf%vertexcenter(i)=findgausscentre(quake,dist,0.e0)
            ! pdf%vertexcenter(i)=findgausscentre(amesh,0.e0)
!              pdf%vertexcenter(i)=findgausscentre(amesh,5.e3)
              pdf%sizeandhigh(i,1)=sqrt(target_area)/3.
              pdf%sizeandhigh(i,2)=1.                           ! relative weight

          enddo
   case default
        write(0,*) 'ERROR in asigning pdf type'
        write(0,*) pdf%pdf_type
        stop
end select

return
end subroutine set_pdf
!#######################################################################################
subroutine faultpdf(pdf,quake_pdf,dist,quake,amesh)
implicit none

! gaussian parameters
  type(pdfinputs) :: pdf
! main structure
  type(mesh) :: amesh
  type(source) :: quake   ! mesh containing earthquake 

! pdf on cells
 real(pr), allocatable, dimension(:) :: quake_pdf
 real(pr), allocatable, dimension(:,:) :: dist

 !  real, dimension(amesh%Ncells) :: pdf

  integer :: i,j,k
  real(pr), dimension(3) :: pdfvertexofcurrentcell
  real(pr) :: d,ug,pdfint
!  real, dimension(amesh%Ncells) :: pdfout

  allocate(quake_pdf(quake%QuakeElemNo))

  quake_pdf = 0.
  do i=1,quake%QuakeElemNo ! loop on all the cells
     do j=1,3
        pdfvertexofcurrentcell(j)=0.
     enddo
     do j=1,pdf%ng ! loop on the gaussian
        do k=1,3 ! loop on the vertex of the cell
           d = dist(quake%QuakeNodes(pdf%vertexcenter(j)),amesh%cell(quake%QuakeElem(i),k))
           ug=exp(-d**2./2./pdf%sizeandhigh(j,1)**2.)
           pdfvertexofcurrentcell(k)=pdfvertexofcurrentcell(k)+ug*pdf%sizeandhigh(j,2)
        enddo
     enddo

     do j=1,3
        quake_pdf(i)=quake_pdf(i)+pdfvertexofcurrentcell(j)
     enddo
     quake_pdf(i)=quake_pdf(i)/3.
  enddo

  pdfint=0.
  do i=1,quake%QuakeElemNo
     pdfint=pdfint+quake_pdf(i)
  enddo
  ug=0.
  do i=1,quake%QuakeElemNo
     quake_pdf(i)=quake_pdf(i)/pdfint
     ug=ug+quake_pdf(i)
  enddo


  !  pdfout=0.
  !  do i=1,amesh%QuakeElemNo
  !     pdfout(amesh%QuakeElem(i))=pdf(i)
  !  enddo
  !open(22,file='pdf.vtk')
  !call dumpmeshvtk(22,amesh)
  !call dumpcellattributevtk(22,amesh,pdfout,'pdf',.true.)
  !close(22)

end subroutine faultpdf
!#######################################################################################
subroutine uniform_pdf(quake_pdf,quake)
implicit none
type(source) :: quake   ! mesh containing earthquake 
! type(mesh) :: amesh
real(pr), allocatable, dimension(:) :: quake_pdf

allocate(quake_pdf(quake%QuakeElemNo))

quake_pdf = 1./float(quake%QuakeElemNo-1)

end subroutine uniform_pdf
!#######################################################################################
subroutine calc_rupt_front(amesh,quake)
  use LAT_time
  use calc_time 
  implicit none
  type(source) :: quake   ! mesh containing earthquake 
  type(mesh) :: amesh,qmesh
  type(diff) :: adiff

   adiff%fast = .true. 

   ! calculate rupture time across whole mesh 
   call timeonevsall2dV2(amesh,adiff)

  return
end subroutine calc_rupt_front
!#######################################################################################
 subroutine pdftoslip(model,amesh,quake,surface,pdf_type,quake_pdf,output_level,dist2)
  use typedef
  use LAT_source
  use LAT_time
  use distance 
! main structure
  type(mesh) :: amesh
  type(source) :: quake   ! mesh containing earthquake 
  type(reflect) :: surface 
  type(model_param) :: model
! pdf on cells
  real(pr), dimension(quake%QuakeElemNo) :: quake_pdf
  real(pr) :: rmin,rmax
  real(pr) :: random
  real(pr) :: cslp,p,length,width,f_area,dx,dmean,mf,db,numrnd,curcum,ave_area
  real(pr), dimension(:), allocatable :: r,dist2
  integer(pin), dimension(:), allocatable :: poss_centres
  real(pr), parameter :: pi=acos(-1.)
!  real :: total_max_slip,prim_max_slip,sec_max_slip
  integer(pin) :: i,j,k,idxmax,cellidx,casp,cnt
  character(7),intent(in) :: pdf_type
  real(pr), dimension(amesh%Ncells) :: pdfout
  integer(pin) :: no_options
  logical :: surf_rupt
!  integer :: nearest_surf_node,nearest_node,jj
  logical :: fast
  integer(pin) :: output_level
  real(pr) :: ave_slip,total_area

! rough estimate of length and width of fault
  call wl(quake,dist,length,width,f_area,dx,output_level)
  f_area = model%actual_area
  length = f_area/width
  model%actual_length = length
  ! set the size of the earthquake
  select case(model%param_type)
   case('sd')
    model%moment = model%sd*f_area*length
    model%Mw = (log10(model%moment)-9.1)*2./3.
    if (output_level > 2) write(0,*) 'given stress drop (MPa)',model%sd/10.**6
    if (output_level > 2) write(0,*) 'estimated magnitude: ',model%mw
   case('mw')
    model%moment=10**(1.5*model%mw+9.1)
    model%sd=model%moment/(length*f_area)
    if (output_level > 2) write(0,*) 'given magnitude',model%mw
    if (output_level > 2) write(0,*) 'estimated stress drop (MPa): ',model%sd/10.**6
  case default
    write(0,*) 'param_type in input file not defined correctly'
    stop
  end select

! Eshelby's constant
  cslp=24./7./pi*model%sd/model%mu
! min/max of the radius distribution
  rmax=width*model%rmax   ! maximum asperity radius set relative to fault width
  rmin=dx*model%rmin      ! minimum fault radius set relative to average elemental area
! fractal parameter
  p=2._pr*7._pr/16._pr*model%moment/model%sd/(rmax-rmin)
! p=7._pr/16._pr*model%moment/model%sd/(rmax-rmin)
! memory allocation
  allocate(r(model%na))
! fractal distribution of the radii (Zeng el al. 1994)
  do i=1,model%na
! D=2
    call random_number(random)
    r(i)=(2._pr*random*float(model%na)/p+rmax**(-2))**(-0.5_pr)
  enddo
! Sort of the distribution by size
  call reorder(r,model%na)

if (abs(maxval(quake_pdf)-minval(quake_pdf)) <= epsilon(1.)) then
  idxmax = ceiling(random*float(quake%QuakeElemNo))  ! random choice
  if (output_level > 2) write(0,*) 'using random location for first subevent '
else
  idxmax = maxloc(quake_pdf,1)      ! locating the maximum of the pdf
  if (output_level > 2) write(0,*)  'using max of pdf fn for location of first subevent '
endif
!print*,'NOTE:: initial index value.......',idxmax
! initializing the slip array
  allocate(slip(quake%QuakeElemNo))
  allocate(slip_prim(quake%QuakeElemNo))
  allocate(slip_sec(quake%QuakeElemNo))
  allocate(poss_centres(quake%QuakeElemNo))
!
  slip = 0._pr
  slip_prim = 0._pr
  slip_sec = 0._pr
! asperity loop
  do i=1,model%na
!    write(0,*) "working on asperity #",i,"/",na
     if (i == 1) then
          cellidx=idxmax
          db=r(i)-1.
          cnt = 0
          do while (r(i) > db)
            cnt=cnt+1
            if (cnt == quake%QuakeElemNo+1 ) then
               write(0,*) "Error: cannot fit largest subevent into designated slipping area"
               write(0,*) "Possible solution is to decrease the Rmax value used in the input file"
               write(0,*) "                 Diameter of largest poss subevent:       ",2.*rmax
               write(0,*) "                 Estimated  earthquake width:             ",width
               stop
            endif
             call random_number(random)
             casp = amesh%cell(quake%QuakeElem(cellidx),min(int(random*3+1),3))
             db = Dist2Border(casp)  ! use this if distance from a line works
             call random_number(random)
             cellidx = ceiling(random*float(quake%QuakeElemNo))  ! random choice
          enddo
     else
        poss_centres = quake%QuakeElem
        no_options = quake%QuakeElemNo
    ! endif
        cnt=0
        db=r(i)-1.
        do while (r(i) > db)
           cnt=cnt+1
           if (cnt == quake%QuakeElemNo+1 ) then
              write(0,*) "Error: cannot fit largest subevent into designated slipping area"
              write(0,*) "Possible solution is to decrease the Rmax value used in the input file"
              write(0,*) "                 Diameter of largest poss subevent:       ",2.*rmax
              write(0,*) "                 Estimated  earthquake width:             ",width
              stop
           endif
           curcum=0.
           call random_number(numrnd)
           cellidx=0
           do while (curcum <= numrnd)
                  cellidx=cellidx+1
!                  if (cellidx == amesh%QuakeElemNo) exit
                  if (cellidx == no_options) exit
                  curcum=curcum+quake_pdf(cellidx)
           enddo
           call random_number(random)
           casp = amesh%cell(poss_centres(cellidx),min(int(random*3+1),3))
           db = Dist2Border(casp)  ! use this if distance from a line works
! remove used choice
           no_options = no_options - 1
           do k = cellidx,no_options
                  poss_centres(k) = poss_centres(k+1)
           enddo
           poss_centres(no_options+1) = 0

        enddo
     endif

! check for secondary pt source
     surf_rupt = .false.
     do j = 1,surface%nnodes
        if ( dist(casp,surface%nodes(j)) <= r(i) ) then
            surf_rupt = .true.
            exit
            ! write(*,*) 'surface rupture....',surface%nodes(j),amesh%py(surface%nodes(j))
        endif
     enddo
    !  do j = 1,quake%NoSurfNodes
    !       if ( dist(casp,quake%SurfNodes(j)) <= r(i) ) then
    !             surf_rupt = .true.
    !             write(*,*) ' we have surface rupture'
    !       endif
    !  enddo
! adding the asperity to the slip distribution
! add loop for secondary source if r(i) > dist2surface
! primary asperity
     do j=1,quake%QuakeElemNo
        dmean=0._pr
        do k=1,3
           dmean=dmean+dist(casp,amesh%cell(quake%QuakeElem(j),k))
        enddo
        dmean=dmean/3._pr

        if (dmean >= r(i) ) then
           cycle
        else
          slip_prim(j) = slip_prim(j)+cslp*sqrt(r(i)**2.-dmean**2.)
        	slip(j) = slip(j)+cslp*sqrt(r(i)**2.-dmean**2.)
         endif
     enddo
! secondary source
     if(surf_rupt) then    ! have surface rupture
       !initialise array: 
      if (.not.allocated(dist2)) allocate(dist2(amesh%Nnodes))
       dist2=infinity
       do j=1,surface%nnodes
           dist2(surface%nodes(j)) = dist(casp,surface%nodes(j)) !surface nodes start with distance from initial subevent
       enddo
       ! calculate distance to from boundary back to nodes in fault
       fast = .true.
       call onevsall2d(amesh,dist2,fast)

       do j=1,quake%QuakeElemNo
          dmean=0._pr
          do k=1,3
             dmean=dmean+dist2(amesh%cell(quake%QuakeElem(j),k))
          enddo
          dmean=dmean/3._pr
          if (dmean >= r(i) ) then
          		cycle
          else
          	slip_sec(j) = slip_sec(j)+cslp*sqrt(r(i)**2.-dmean**2.)
          	slip(j) = slip(j)+cslp*sqrt(r(i)**2.-dmean**2.)
          endif
       enddo
      !  deallocate(dist2)
     endif  ! secondary source conditional 

  enddo

! renormalization moment
  mf=0._pr
! seismic moment = mu*area*average slip
! as average slip is weighted based on the area of the cell the
! total seismic moment is based on the sum of the moments of all the cells
  do i=1,quake%QuakeElemNo
     mf=mf+slip(i)*area(quake%QuakeElem(i))
  enddo
  mf=mf*model%mu

  !write(0,*) "difference to moment target : ",100*(mf-model%moment)/model%moment,"%"
  if (Rmax /rmin < 2.) then
    write(0,*)  "Warning: the ratio of largest to smallest < 2  "
    write(0,*)   "i.e. Rmax/rmin:             " ,Rmax /rmin
    write(0,*)   "a value above 3 is more desirable"
    write(0,*)   "consider increasing Rmax or decreasing rmin"
  endif

  !write(*,*) "comparing actual to moment target : ",mf,model%moment
  ave_slip = 0.
  total_area = 0.
  do i=1,quake%QuakeElemNo
     slip(i)= slip(i)/mf*model%moment
     slip_prim(i) = slip_prim(i)/mf*model%moment
     slip_sec(i) = slip_sec(i)/mf*model%moment
     ave_slip = ave_slip + slip(i)*area(quake%QuakeElem(i))
     total_area = total_area + area(quake%QuakeElem(i))
  enddo
! average slip is weighted based on the area of the elements
!  i.e. ave_slip = Sum (a(i)*slip(i))/ Sum (a(i))

  ave_slip = ave_slip/total_area

   if (output_level > 1)  write(*,*) 'total slipping area      ', total_area
   if (output_level > 1)  write(*,*) 'average slip:            ', ave_slip
   if (output_level > 1)  write(*,*) 'moment of event          ', model%mu*total_area*ave_slip
   if (output_level > 1)  write(*,*) 'Moment magnitude of event', (log10(model%mu*total_area*ave_slip)-9.1)*2./3.

  return
end subroutine pdftoslip
!#######################################################################################
subroutine wl(quake,dist,length,width,calc_area,dx,output_level)
! main structure
  ! type(mesh) :: amesh
  type(source) :: quake   ! mesh containing earthquake 
  real(pr),allocatable,dimension(:,:) :: dist

  real(pr) :: length,width,calc_area,dx
  integer(pin) :: output_level
  integer(pin) :: i,j

! fault length approximation
  length=0.
  do i=1,quake%QuakeBorder_NodesNo
     do j=i,quake%QuakeBorder_NodesNo
        length=max(length,dist(quake%QuakeBorder_Nodes(i),quake%QuakeBorder_Nodes(j)))
     enddo
  enddo
  if (output_level > 2) write(0,*) "estimated length :",length
! area of the fault
  calc_area=0.
  do i=1,quake%QuakeElemNo
    calc_area=calc_area+area(quake%QuakeElem(i))
  enddo

  if (output_level > 2)  write(0,*) " area :",calc_area
! fault width approximation
  width=calc_area/length
  if (output_level > 2)  write(0,*) "estimated width :",width
! dx approximation
  dx=sqrt(calc_area/float(quake%QuakeElemNo))
  return
end subroutine wl
!#######################################################################################
!##########################################################################
function findcentre(quake,dist)
use typedef
  type(source) :: quake   ! mesh containing earthquake 
  real(pr),allocatable,dimension(:,:) :: dist

  integer :: findcentre

  integer :: i,j
  real :: sd,mean,sdmin
  real, parameter :: infini=1.e32

  findcentre=0
  sdmin=infini
  do i=1,quake%QuakeNodesNo
     mean=0.
     do j=1,quake%QuakeBorder_NodesNo
        mean=mean+dist(quake%QuakeNodes(i),quake%QuakeBorder_Nodes(j))
     enddo
     mean=mean/float(quake%QuakeBorder_NodesNo)
     sd=0.
     do j=1,quake%QuakeBorder_NodesNo
        sd=sd+(dist(quake%QuakeNodes(i),quake%QuakeBorder_Nodes(j))-mean)**2.
     enddo
     if (sd < sdmin) then
        sdmin=sd
        findcentre=i
     endif
  enddo
  return
end function findcentre
!##########################################################################
!##########################################################################
function findgausscentre(quake,dist,boundary_dist)
! pick gaussian centre with the condition that it cannot be within distance 'dist'
! from boundary edge
  use typedef
    type(source) :: quake   ! mesh containing earthquake 
    real(pr),allocatable,dimension(:,:) :: dist

    integer :: findgausscentre
    real :: boundary_dist
    real :: x
    integer :: i,j,inx,id
    integer,allocatable,dimension(:) :: keep_id
    logical :: good2use

    allocate(keep_id(quake%QuakeNodesNo))
    inx = 0
    do i=1,quake%QuakeNodesNo
      good2use = .true.
       do j=1,quake%QuakeBorder_NodesNo
         if (dist(quake%QuakeNodes(i),quake%QuakeBorder_Nodes(j)) < boundary_dist) then
              good2use = .false.
              cycle
         endif
       enddo
        if (good2use) then
          inx = inx+1
          keep_id(inx) = i
        endif
    enddo

    call random_number(x)

    if (inx == 0) then
        write(*,*) 'ERROR in findgausscentre, edge is too close relative to fault size'
        stop
    endif
    id = floor(x*float(inx))+1
    findgausscentre = keep_id(id)
!    findgausscentre = 53
!    print*,'WARNING: gaussian centre is set to 53!!!!!!!'

end function findgausscentre
!##########################################################################
end module makepdf
