module makepdf
use lateration
implicit none

contains

!#######################################################################################
subroutine set_quake_pdf(pdf,quake_pdf,amesh,defined_area)
  type(pdfinputs) :: pdf
  type(mesh) :: amesh
  logical :: defined_area
  real, allocatable,dimension(:) :: quake_pdf
  integer :: i
  real :: ug, pdfint

  allocate(quake_pdf(amesh%QuakeElemNo))
  quake_pdf = pdf%g_pdf

return
end subroutine set_quake_pdf
!#######################################################################################
subroutine set_pdf(pdf,amesh,target_area)
implicit none
! gaussian parameters
  type(pdfinputs) :: pdf
! main structure
  type(mesh) :: amesh
  integer :: i
  real, intent(in) :: target_area
  real :: ran


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
              pdf%vertexcenter(i)=findgausscentre(amesh,0.e0)
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
subroutine faultpdf(pdf,quake_pdf,amesh)
implicit none

! gaussian parameters
  type(pdfinputs) :: pdf
! main structure
  type(mesh) :: amesh
! pdf on cells
 real, allocatable, dimension(:) :: quake_pdf
!  real, dimension(amesh%Ncells) :: pdf

  integer :: i,j,k
  real, dimension(3) :: pdfvertexofcurrentcell
  real :: d,ug,pdfint
!  real, dimension(amesh%Ncells) :: pdfout

  allocate(quake_pdf(amesh%QuakeElemNo))

  quake_pdf = 0.
  do i=1,amesh%QuakeElemNo ! loop on all the cells
     do j=1,3
        pdfvertexofcurrentcell(j)=0.
     enddo
     do j=1,pdf%ng ! loop on the gaussian
        do k=1,3 ! loop on the vertex of the cell
           d=amesh%dist(amesh%QuakeNodes(pdf%vertexcenter(j)),amesh%cell(amesh%QuakeElem(i),k))
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
  do i=1,amesh%QuakeElemNo
     pdfint=pdfint+quake_pdf(i)
  enddo
  ug=0.
  do i=1,amesh%QuakeElemNo
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
subroutine uniform_pdf(quake_pdf,amesh)
implicit none

type(mesh) :: amesh
real, allocatable, dimension(:) :: quake_pdf

allocate(quake_pdf(amesh%QuakeElemNo))

quake_pdf = 1./float(amesh%QuakeElemNo-1)

end subroutine uniform_pdf
!#######################################################################################
 subroutine pdftoslip(model,amesh,pdf_type,quake_pdf,output_level)
  use typedef
! main structure
  type(mesh) :: amesh
  type(model_param) :: model
! pdf on cells
  real, dimension(amesh%QuakeElemNo) :: quake_pdf
  real :: rmin,rmax
  real :: random
  real :: cslp,p,length,width,area,dx,dmean,mf,db,numrnd,curcum,ave_area
  real, dimension(:), allocatable :: r,dist2
  integer, dimension(:), allocatable :: poss_centres
  real, parameter :: pi=acos(-1.)
!  real :: total_max_slip,prim_max_slip,sec_max_slip
  integer :: i,j,k,idxmax,cellidx,casp,cnt
  character(7),intent(in) :: pdf_type
  real, dimension(amesh%Ncells) :: pdfout
  integer :: no_options
  logical :: surf_rupt
!  integer :: nearest_surf_node,nearest_node,jj
  logical :: fast
  integer :: output_level
  real :: ave_slip,total_area

! rough estimate of length and width of fault
  call wl(amesh,length,width,area,dx,output_level)
  area = model%actual_area
  length = area/width
  model%actual_length = length
  ! set the size of the earthquake
  select case(model%param_type)
   case('sd')
    model%moment = model%sd*area*length
    model%Mw = (log10(model%moment)-9.1)*2./3.
    if (output_level > 2) write(0,*) 'given stress drop (MPa)',model%sd/10.**6
    if (output_level > 2) write(0,*) 'estimated magnitude: ',model%mw
   case('mw')
    model%moment=10**(1.5*model%mw+9.1)
    model%sd=model%moment/(length*area)
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
  p=2.*7./16.*model%moment/model%sd/(rmax-rmin)
! p=7./16.*model%moment/model%sd/(rmax-rmin)
! memory allocation
  allocate(r(model%na))
! fractal distribution of the radii (Zeng el al. 1994)
  do i=1,model%na
! D=2
    call random_number(random)
    r(i)=(2.*random*float(model%na)/p+rmax**(-2.))**(-.5)
  enddo
! Sort of the distribution by size
  call reorder(r,model%na)

! select location of first and largest asperity
!select case(pdf_type)
!  case('uniform')
!! uniform distribution
!      if (output_level > 2) write(0,*) 'choosen uniform pdf'
!      idxmax = ceiling(random*float(amesh%QuakeElemNo))  ! random choice
!  case('gauss')
!      if (output_level > 2) write(0,*) 'choosen gaussian pdf'
!      idxmax = maxloc(pdf,1)      ! locating the maximum of the pdf
!  case default
!      idxmax = ceiling(random*float(amesh%QuakeElemNo))  ! random choice
!end select

if (abs(maxval(quake_pdf)-minval(quake_pdf)) <= epsilon(1.)) then
  idxmax = ceiling(random*float(amesh%QuakeElemNo))  ! random choice
  if (output_level > 2) write(0,*) 'using random location for first subevent '
else
  idxmax = maxloc(quake_pdf,1)      ! locating the maximum of the pdf
  if (output_level > 2) write(0,*)  'using max of pdf fn for location of first subevent '
endif
!print*,'NOTE:: initial index value.......',idxmax
! initializing the slip array
  allocate(amesh%slip(amesh%QuakeElemNo))
  allocate(amesh%slip_prim(amesh%QuakeElemNo))
  allocate(amesh%slip_sec(amesh%QuakeElemNo))
  allocate(poss_centres(amesh%QuakeElemNo))
!
  amesh%slip = 0.
  amesh%slip_prim = 0.
  amesh%slip_sec = 0.
! asperity loop
  do i=1,model%na
!    write(0,*) "working on asperity #",i,"/",na
     if (i == 1) then
          cellidx=idxmax
          db=r(i)-1.
          cnt = 0
          do while (r(i) > db)
            cnt=cnt+1
            if (cnt == amesh%QuakeElemNo+1 ) then
               write(0,*) "Error: cannot fit largest subevent into designated slipping area"
               write(0,*) "Possible solution is to decrease the Rmax value used in the input file"
               write(0,*) "                 Diameter of largest poss subevent:       ",2.*rmax
               write(0,*) "                 Estimated  earthquake width:             ",width
               stop
            endif
             call random_number(random)
             casp = amesh%cell(amesh%QuakeElem(cellidx),min(int(random*3+1),3))
             db = amesh%Dist2Border(casp)  ! use this if distance from a line works
             call random_number(random)
             cellidx = ceiling(random*float(amesh%QuakeElemNo))  ! random choice
          enddo
     else
        poss_centres = amesh%QuakeElem
        no_options = amesh%QuakeElemNo
    ! endif
        cnt=0
        db=r(i)-1.
        do while (r(i) > db)
           cnt=cnt+1
           if (cnt == amesh%QuakeElemNo+1 ) then
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
           db = amesh%Dist2Border(casp)  ! use this if distance from a line works
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
     do j = 1,amesh%NoSurfNodes
          if ( amesh%dist(casp,amesh%SurfNodes(j)) <= r(i) ) then
                surf_rupt = .true.
          endif
     enddo
! adding the asperity to the slip distribution
! add loop for secondary source if r(i) > dist2surface
! primary asperity
     do j=1,amesh%QuakeElemNo
        dmean=0.
        do k=1,3
           dmean=dmean+amesh%dist(casp,amesh%cell(amesh%QuakeElem(j),k))
        enddo
        dmean=dmean/3.

        if (dmean >= r(i) ) then
           cycle
        else
          amesh%slip_prim(j)=amesh%slip_prim(j)+cslp*sqrt(r(i)**2.-dmean**2.)
        	amesh%slip(j)=amesh%slip(j)+cslp*sqrt(r(i)**2.-dmean**2.)
         endif
     enddo
! secondary source
     if(surf_rupt) then    ! have surface rupture
       !initialise array
       allocate(dist2(amesh%Nnodes))
       dist2=infinity
       do j=1,amesh%NoSurfNodes
            dist2(amesh%SurfNodes(j)) = amesh%dist(casp,amesh%SurfNodes(j))
       enddo
       ! calculate distance to from boundary back to nodes in fault
       fast = .true.
       call onevsall2d(amesh,dist2,fast)

       do j=1,amesh%QuakeElemNo
          dmean=0.
          do k=1,3
             dmean=dmean+dist2(amesh%cell(amesh%QuakeElem(j),k))
          enddo
          dmean=dmean/3.
          if (dmean >= r(i) ) then
          		cycle
          else
          	amesh%slip_sec(j)=amesh%slip_sec(j)+cslp*sqrt(r(i)**2.-dmean**2.)
          	amesh%slip(j)=amesh%slip(j)+cslp*sqrt(r(i)**2.-dmean**2.)
          endif
       enddo
       deallocate(dist2)
     endif

  enddo

! renormalization moment
  mf=0.
! seismic moment = mu*area*average slip
! as average slip is weighted based on the area of the cell the
! total seismic moment is based on the sum of the moments of all the cells
  do i=1,amesh%QuakeElemNo
     mf=mf+amesh%slip(i)*amesh%area(amesh%QuakeElem(i))
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
  do i=1,amesh%QuakeElemNo
     amesh%slip(i)=amesh%slip(i)/mf*model%moment
     amesh%slip_prim(i) = amesh%slip_prim(i)/mf*model%moment
     amesh%slip_sec(i) = amesh%slip_sec(i)/mf*model%moment
     ave_slip = ave_slip+amesh%slip(i)*amesh%area(amesh%QuakeElem(i))
     total_area = total_area + amesh%area(amesh%QuakeElem(i))
  enddo
! average slip is weighted based on the area of the elements
!  i.e. ave_slip = Sum (a(i)*slip(i))/ Sum (a(i)

  ave_slip = ave_slip/total_area

   if (output_level > 1)  write(*,*) 'total slipping area      ', total_area
   if (output_level > 1)  write(*,*) 'average slip:            ', ave_slip
   if (output_level > 1)  write(*,*) 'moment of event          ', model%mu*total_area*ave_slip
   if (output_level > 1)  write(*,*) 'Moment magnitude of event', (log10(model%mu*total_area*ave_slip)-9.1)*2./3.

  return
end subroutine pdftoslip
!#######################################################################################
subroutine wl(amesh,length,width,area,dx,output_level)

! main structure
  type(mesh) :: amesh
  real :: length,width,area,dx
  integer :: output_level
  integer :: i,j

! fault length approximation
  length=0.
  do i=1,amesh%QuakeBorder_NodesNo
     do j=i,amesh%QuakeBorder_NodesNo
        length=max(length,amesh%dist(amesh%QuakeBorder_Nodes(i),amesh%QuakeBorder_Nodes(j)))
     enddo
  enddo
  if (output_level > 2) write(0,*) "estimated length :",length
! area of the fault
  area=0.
  do i=1,amesh%QuakeElemNo
     area=area+amesh%area(amesh%QuakeElem(i))
  enddo

  if (output_level > 2)  write(0,*) " area :",area
! fault width approximation
  width=area/length
  if (output_level > 2)  write(0,*) "estimated width :",width
! dx approximation
  dx=sqrt(area/float(amesh%QuakeElemNo))
  return
end subroutine wl
!#######################################################################################
!##########################################################################
function findcentre(amesh)
use typedef
  type(mesh) :: amesh
  integer :: findcentre

  integer :: i,j
  real :: sd,mean,sdmin
  real, parameter :: infini=1.e32

  findcentre=0
  sdmin=infini
  do i=1,amesh%QuakeNodesNo
     mean=0.
     do j=1,amesh%QuakeBorder_NodesNo
        mean=mean+amesh%dist(amesh%QuakeNodes(i),amesh%QuakeBorder_Nodes(j))
     enddo
     mean=mean/float(amesh%QuakeBorder_NodesNo)
     sd=0.
     do j=1,amesh%QuakeBorder_NodesNo
        sd=sd+(amesh%dist(amesh%QuakeNodes(i),amesh%QuakeBorder_Nodes(j))-mean)**2.
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
function findgausscentre(amesh,dist)
! pick gaussian centre with the condition that it cannot be within distance 'dist'
! from boundary edge
  use typedef
    type(mesh) :: amesh
    integer :: findgausscentre
    real :: dist
    real :: x
    integer :: i,j,inx,id
    integer,allocatable,dimension(:) :: keep_id
    logical :: good2use

    allocate(keep_id(amesh%QuakeNodesNo))
    inx = 0
    do i=1,amesh%QuakeNodesNo
      good2use = .true.
       do j=1,amesh%QuakeBorder_NodesNo
         if (amesh%dist(amesh%QuakeNodes(i),amesh%QuakeBorder_Nodes(j)) < dist) then
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
