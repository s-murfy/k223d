module makepdf
use lateration
implicit none

contains

!#######################################################################################
subroutine set_pdf(gausspar,amesh,target_area)
implicit none
! gaussian parameters
  type(pdfinputs) :: gausspar
! main structure
  type(mesh) :: amesh
  integer :: i
  real, intent(in) :: target_area


select case(gausspar%pdf_type)
  case('uniform')
! uniform distribution
      write(*,*) 'uniform pdf choosen'
  case('gauss')
! gaussian distribution
!          call random_number(ran)
!          i = 1+int(ran*3.) ! randomly choose between 1-4 gaussian functions
!          print*,'number of gaussians.....',i
!          gausspar%ng=i
          allocate(gausspar%vertexcenter(gausspar%ng),gausspar%sizeandhigh(gausspar%ng,2))
          do i = 1,gausspar%ng
          ! randomly place centres of gaussians with the condtion that they are at least 5km from boundary
              gausspar%vertexcenter(i)=findgausscentre(amesh,5.e3)
              gausspar%sizeandhigh(i,1)=sqrt(target_area)/3.
              gausspar%sizeandhigh(i,2)=1.                           ! relative weight
          enddo
   case default
        write(*,*) 'ERROR in asigning pdf type'
        write(*,*) gausspar%pdf_type
        stop
end select

return
end subroutine set_pdf
!#######################################################################################
subroutine faultpdf(gausspar,pdf,amesh)
implicit none

! gaussian parameters
  type(pdfinputs) :: gausspar
! main structure
  type(mesh) :: amesh
! pdf on cells
  real, dimension(amesh%QuakeElemNo) :: pdf

  integer :: i,j,k
  real, dimension(3) :: pdfvertexofcurrentcell
  real :: d,ug,pdfint,totarea

  totarea=0.
  do i=1,amesh%QuakeElemNo ! loop on all the cells
     do j=1,3
        pdfvertexofcurrentcell(j)=0.
     enddo
     do j=1,gausspar%ng ! loop on the gaussian
        do k=1,3 ! loop on the vertex of the cell
           d=amesh%dist(amesh%QuakeNodes(gausspar%vertexcenter(j)),amesh%cell(amesh%QuakeElem(i),k))
           ug=exp(-d**2./2./gausspar%sizeandhigh(j,1)**2.)
           pdfvertexofcurrentcell(k)=pdfvertexofcurrentcell(k)+ug*gausspar%sizeandhigh(j,2)
        enddo
     enddo
     pdf(i)=0.
     do j=1,3
        pdf(i)=pdf(i)+pdfvertexofcurrentcell(j)
     enddo
     pdf(i)=pdf(i)/3.*amesh%area(amesh%QuakeElem(i))
     totarea=totarea+amesh%area(amesh%QuakeElem(i))
  enddo
  pdfint=0.
  do i=1,amesh%QuakeElemNo
     pdf(i)=pdf(i)/totarea
     pdfint=pdfint+pdf(i)
  enddo
  ug=0.
  do i=1,amesh%QuakeElemNo
     pdf(i)=pdf(i)/pdfint
     ug=ug+pdf(i)
  enddo
end subroutine faultpdf
!#######################################################################################
 subroutine pdftoslip(model,amesh,pdf_type,pdf)
  use typedef
! main structure
  type(mesh) :: amesh
  type(model_param) :: model
! pdf on cells
  real, dimension(amesh%QuakeElemNo) :: pdf
  real :: rmin,rmax
  real :: random
  real :: cslp,p,length,width,area,dx,dmean,mf,db,numrnd,curcum
  real, dimension(:), allocatable :: r
  real, parameter :: pi=acos(-1.)
  integer :: i,j,k,idxmax,cellidx,casp,cnt
  character(7),intent(in) :: pdf_type

! rough estimate of length and width of fault
  call wl(amesh,length,width,area,dx)

  ! set the size of the earthquake
  select case(model%param_type)
   case('sd')
    model%moment = model%sd*area*length
    model%Mw = (log10(model%moment)-9.1)*2./3.
    write(*,*) 'given stress drop (MPa)',model%sd/10.**6
    write(*,*) 'estimated magnitude: ',model%mw
   case('mw')
    model%moment=10**(1.5*model%mw+9.1)
    model%sd=model%moment/(length*area)
    write(*,*) 'given magnitude',model%mw
    write(*,*) 'estimated stress drop (MPa): ',model%sd/10.**6
  case default
    write(*,*) 'param_type in input file not defined correctly'
    stop
  end select

! Eshelby's constant
  cslp=24/7/pi*model%sd/model%mu
! min/max of the radius distribution
  rmax=width*model%rmax   ! maximum asperity radius set relative to fault width
  rmin=dx*model%rmin      ! minimum fault radius set relative to average elemental area
! fractal parameter
  p=2.*7./16.*model%moment/model%sd/(rmax-rmin)
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
select case(pdf_type)
case('uniform')
! uniform distribution
  idxmax = ceiling(random*float(amesh%QuakeElemNo))  ! random choice
  case('gauss')
  idxmax=maxloc(pdf,1)      ! locating the maximum of the pdf
end select

! initializing the slip array
  allocate(amesh%slip(amesh%Ncells))
  amesh%slip=0.
! asperity loop
  do i=1,model%na
!    write(0,*) "working on asperity #",i,"/",na
     if (i == 1) then
        cellidx=idxmax
        call random_number(random)
        casp=amesh%cell(amesh%QuakeElem(cellidx),min(int(random*3+1),3))
     else
        cnt=0
        db=r(i)-1.
        do while (r(i) > db)
           cnt=cnt+1
           if (cnt == 100 ) then
              write(0,*) "Error: fault size too small relative to asperity size"
              stop
           endif
           curcum=0.
           call random_number(numrnd)
           cellidx=0

           select case(pdf_type)
             case('uniform')
         ! uniform distribution
                  cellidx =  ceiling(numrnd*float(amesh%QuakeElemNo))
             case('gauss')
         ! gaussian distribution
               do while (curcum <= numrnd)
                  cellidx=cellidx+1
                  if (cellidx == amesh%QuakeElemNo) exit
                  curcum=curcum+pdf(cellidx)
               enddo
           end select

           call random_number(random)
           casp=amesh%cell(amesh%QuakeElem(cellidx),min(int(random*3+1),3))
           db = amesh%Dist2Border(casp)  ! use this if distance from a line works
        enddo
     endif
! adding the asperity to the slip distribution
     do j=1,amesh%QuakeElemNo
        dmean=0.
        do k=1,3
           dmean=dmean+amesh%dist(casp,amesh%cell(amesh%QuakeElem(j),k))
        enddo
        dmean=dmean/3.
        if (dmean > r(i) ) cycle
        amesh%slip(j)=amesh%slip(j)+cslp*sqrt(r(i)**2.-dmean**2.)
     enddo
  enddo

! renormalization
  mf=0.
  do i=1,amesh%QuakeElemNo
     mf=mf+amesh%slip(i)*amesh%area(amesh%QuakeElem(i))
  enddo
  mf=mf*model%mu
  write(*,*) "difference to moment target : ",100*(mf-model%moment)/model%moment,"%"
  do i=1,amesh%QuakeElemNo
     amesh%slip(i)=amesh%slip(i)/mf*model%moment
  enddo
  return
end subroutine pdftoslip

!#######################################################################################
subroutine wl(amesh,length,width,area,dx)

! main structure
  type(mesh) :: amesh
  real :: length,width,area,dx

  integer :: i,j

! fault length approximation
  length=0.
  do i=1,amesh%QuakeBorder_NodesNo
     do j=i,amesh%QuakeBorder_NodesNo
        length=max(length,amesh%dist(amesh%QuakeBorder_Nodes(i),amesh%QuakeBorder_Nodes(j)))
     enddo
  enddo
  write(0,*) "length :",length
! area of the fault
  area=0.
  do i=1,amesh%QuakeElemNo
     area=area+amesh%area(amesh%QuakeElem(i))
  enddo
  write(0,*) "area :",area
! fault width approximation
  width=area/length
  write(0,*) "width :",width
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
              exit
         endif
       enddo
        if (good2use) then
          inx = inx+1
          keep_id(inx) = i
        endif
    enddo

    call random_number(x)

    if (inx == 0) then
        write(*,*) 'ERROR in findgausscentre, edge is too large relative to fault size'
        stop
    endif
    id = floor(x*float(inx))+1
    findgausscentre = keep_id(id)

end function findgausscentre
!##########################################################################
end module makepdf
