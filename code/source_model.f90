module source_model 
use generic
use LAT_mesh
use LAT_source
use mesh_geom
use model_setup
use, intrinsic :: iso_fortran_env, only: error_unit

implicit none

contains 

subroutine pdftoslip(model, amesh, geom, pdf,src)
!==============================================================================
! Generates a stochastic slip distribution on a triangular fault mesh using
! the composite source model of Zeng et al. (1994).
!
! A fractal distribution of circular sub-events with radii between rmin and
! rmax is placed on the fault surface, weighted by the PDF in pdf%distrib.
! Where a sub-event reaches the free surface a secondary reflected source is
! added using a trilateration propagation from the surface boundary back into the
! fault (geom%surface). The total slip is renormalised to match the seismic
! moment implied by the target magnitude model%mw.
!
! Inputs:  amesh  — triangular fault mesh
!          geom   — precomputed mesh geometry (areas, distances, fault width)
!          model  — source parameters (mw, mu, na, rmin, rmax)
!          pdf    — slip probability distribution (uniform or defined)
! Output:  src%slip — slip per cell [m], allocated here
!==============================================================================
  use distance 

  type(mesh),           intent(in)      :: amesh
  type(model_param),    intent(in)      :: model   !
  type(pdfinputs),      intent(in)      :: pdf     ! 
  type(source),         intent(inout)   :: src     ! inout: slip distribution calculated
  type(mesh_geometry),  intent(in)      :: geom
  
! local variables 
  type(diff) :: adiff
  real(pr) :: rmin,rmax
  real(pr) :: random,moment,sd
  real(pr) :: cslp,p,length,width,dx,dmean,mf,db,numrnd,curcum
  real(pr), dimension(:), allocatable :: r,dist2
  integer(pin), dimension(:), allocatable :: poss_centres,cell_options
  real(pr), parameter :: pi=acos(-1._pr)
  integer(pin) :: i,j,k,idxmax,cellidx,casp,cnt
  integer(pin) :: no_options
  logical :: surf_rupt
  real(pr) :: ave_slip,total_area

  adiff%fast = .true. ! used when calculating the distance of nodes to surface

  moment=10._pr**(1.5_pr*model%mw+9.1_pr)
  total_area = sum(geom%area)
  length = total_area/geom%W
  sd = moment/(length*total_area) !estimate expected stress drop 
  dx=sqrt(total_area/real(amesh%Ncells,pr)) ! estimate of cell length
! Eshelby's constant
  cslp=24._pr/7._pr/pi*sd/model%mu
! min/max of the radius distribution
  rmax=geom%W*model%rmax   ! maximum asperity radius set relative to fault width
  rmin=dx*model%rmin      ! minimum fault radius set relative to average elemental area
! fractal parameter
  if (rmax <= rmin) then
    write(error_unit,*) 'ERROR: Rmax must be greater than Rmin'
    stop
  endif
  p=2._pr*7._pr/16._pr*moment/sd/(rmax-rmin)

! memory allocation
  allocate(r(model%na))
  do i=1,model%na
    call random_number(random)
    ! Fractal radius distribution: D=2 case, Zeng et al. (1994) eq. 7
    r(i)=(2._pr*random*real(model%na,pr)/p+rmax**(-2._pr))**(-0.5_pr) ! This is D=2 case
  enddo
! Sort of the distribution by size
  call reorder(r,model%na)

! Uniform PDF: place first sub-event randomly; otherwise use PDF maximum
if (abs(maxval(pdf%distrib)-minval(pdf%distrib)) <= epsilon(1._pr)) then
  idxmax = ceiling(random*real(amesh%Ncells,pr))  ! random choice
  if (verbose > 2) write(error_unit,*) 'using random location for first subevent '
else
  idxmax = maxloc(pdf%distrib,1)      ! locating the maximum of the pdf
  ! write(0,*) 'using max of pdf fn for location of first subevent ',idxmax
  if (verbose > 2) write(error_unit,*)  'using max of pdf fn for location of first subevent '
endif
!print*,'NOTE:: initial index value.......',idxmax
! initializing the slip array
  allocate(src%slip(amesh%Ncells))
  allocate(poss_centres(amesh%Ncells),cell_options(amesh%Ncells))
!
  src%slip = 0._pr
  do i = 1,amesh%Ncells
    cell_options(i) = i 
  enddo 

! asperity loop
  do i=1,model%na
!    write(0,*) "working on asperity #",i,"/",na
     if (i == 1) then
          cellidx=idxmax
          db=r(i)-1._pr ! initialise db < r(i) to enter the while loop
          cnt = 0
          do while (r(i) > db)
            cnt=cnt+1
            if (cnt == amesh%Ncells+1 ) then
               write(error_unit,*) "Error: cannot fit largest subevent into designated slipping area"
               write(error_unit,*) "Possible solution is to decrease the Rmax value used in the input file"
               write(error_unit,*) "                 Diameter of largest poss subevent:       ",2._pr*rmax
               write(error_unit,*) "                 Estimated  earthquake width:             ",geom%W
               stop
            endif
             call random_number(random)
             casp = amesh%cell(cellidx,min(int(random*3+1),3))
             db = geom%dist_to_border(casp)  ! use this if distance from a line works
             call random_number(random)
             cellidx = ceiling(random*real(amesh%Ncells,pr))  ! random choice
          enddo
     else
        poss_centres = cell_options   ! list where used options are removed
        no_options = amesh%Ncells

        cnt=0
        db=r(i)-1._pr
        do while (r(i) > db)
           cnt=cnt+1
           if (cnt == amesh%Ncells+1 ) then
              write(error_unit,*) "Error: cannot fit largest subevent into designated slipping area"
              write(error_unit,*) "Possible solution is to decrease the Rmax value used in the input file"
              write(error_unit,*) "                 Diameter of largest poss subevent:       ",2.*rmax
              write(error_unit,*) "                 Estimated  earthquake width:             ",geom%W
              stop
           endif
           curcum=0._pr
           call random_number(numrnd)
           cellidx=0
           do while (curcum <= numrnd)
                  cellidx=cellidx+1
                  if (cellidx == no_options) exit
                !  if (cellidx == amesh%Ncells) exit 
                  curcum=curcum+pdf%distrib(cellidx)
           enddo
           call random_number(random)
           casp = amesh%cell(poss_centres(cellidx),min(int(random*3+1),3))
           db = geom%dist_to_border(casp)  ! distance from border to interior nodes 
! remove used choice, not sure if this is really needed 
           no_options = no_options - 1
           do k = cellidx,no_options
                  poss_centres(k) = poss_centres(k+1)
           enddo
           poss_centres(no_options+1) = 0

        enddo
     endif

! check for secondary pt source required - i.e. subevent touches surface
     surf_rupt = .false.
     do j = 1,geom%surface%nnodes
        if ( geom%dist(casp,geom%surface%nodes(j)) <= r(i) ) then
            surf_rupt = .true.
            exit
            ! write(*,*) 'surface rupture....',surface%nodes(j),amesh%py(surface%nodes(j))
        endif
     enddo
! adding the asperity to the slip distribution
! add loop for secondary source if r(i) > dist2surface
! primary asperity
     do j=1,amesh%Ncells
        dmean=0._pr
        do k=1,3
           dmean=dmean+geom%dist(casp,amesh%cell(j,k))
        enddo
        dmean=dmean/3._pr

        if (dmean >= r(i) ) then
           cycle
        else
            ! slip_prim(j) = slip_prim(j)+cslp*sqrt(r(i)**2.-dmean**2.)
        	src%slip(j) = src%slip(j)+cslp*sqrt(r(i)**2.-dmean**2.)
        endif
     enddo
! Secondary source: if sub-event reaches the free surface, reflect slip
! back into the fault via distance from surface nodes
     if(surf_rupt) then    ! have surface rupture
       !initialise array: 
      if (.not.allocated(dist2)) allocate(dist2(amesh%Nnodes))
       dist2=infinity
       do j=1,geom%surface%nnodes
           dist2(geom%surface%nodes(j)) = geom%dist(casp,geom%surface%nodes(j)) !surface nodes start with distance from initial subevent
       enddo
       ! calculate distance to from boundary back to nodes in fault
       call onevsall2d(amesh, dist2, geom%ntoc, geom%nton, adiff)
       do j=1,amesh%Ncells
          dmean=0._pr
          do k=1,3
             dmean=dmean+dist2(amesh%cell(j,k))
          enddo
          dmean=dmean/3._pr
          if (dmean >= r(i) ) then
          		cycle
          else
          	src%slip(j) = src%slip(j)+cslp*sqrt(r(i)**2.-dmean**2.)
          endif
       enddo
      !  deallocate(dist2)
     endif  ! secondary source conditional 

  enddo
! Renormalise slip to match target seismic moment
! mf = mu * sum(slip(i) * area(i)) before normalisation
  mf=0._pr
  do i=1,amesh%Ncells
     mf=mf+src%slip(i)*geom%area(i)
  enddo
  mf=mf*model%mu

  !write(0,*) "difference to moment target : ",100*(mf-model%moment)/model%moment,"%"
  if (Rmax /rmin < 2._pr) then
    write(error_unit,*)  "Warning: the ratio of largest to smallest < 2  "
    write(error_unit,*)   "i.e. Rmax/rmin:             " ,Rmax /rmin
    write(error_unit,*)   "a value above 3 is more desirable"
    write(error_unit,*)   "consider increasing Rmax or decreasing rmin"
  endif

  !write(*,*) "comparing actual to moment target : ",mf,model%moment
  ave_slip = 0._pr
  do i=1,amesh%Ncells
     src%slip(i)= src%slip(i)/mf*moment
     ave_slip = ave_slip + src%slip(i)*geom%area(i)
  enddo
! average slip is weighted based on the area of the elements
!  i.e. ave_slip = Sum (a(i)*slip(i))/ Sum (a(i))

  ave_slip = ave_slip/total_area

   write(error_unit,'(a, es12.4)') 'total slipping area      ', total_area
   write(error_unit,'(a, f8.2, a)') 'average slip:            ', ave_slip,' m'
   write(error_unit,'(a, es12.4)') 'moment of event          ', model%mu*total_area*ave_slip
   write(error_unit,'(a, f8.2)') 'Moment magnitude of event', (log10(model%mu*total_area*ave_slip)-9.1_pr)*2._pr/3._pr

  return
end subroutine pdftoslip
!#######################################################################################
!*******************************************************
subroutine reorder(a,n)
!######################################################
! Author : André Herrero
! Contact : andherit@gmail.com, andre.herrero@ingv.it
! Public Domain (CC0 1.0 Universal)
!######################################################
implicit none

integer(pin) n
real(pr) a(n)
integer(pin) i
real(pr) mem
logical done

done=.false.
do while (.not.done)
       done=.true.
       do i=1,n-1
              if (a(i).lt.a(i+1)) then
                     mem=a(i)
                     a(i)=a(i+1)
                     a(i+1)=mem
                     done=.false.
              endif
       enddo
enddo

return
end subroutine reorder
!#######################################################################################
subroutine calc_rupt_front(amesh,geom,src)
  use LAT_time
  use time 

  type(mesh),           intent(in)      :: amesh
  type(source),         intent(inout)   :: src     ! inout: slip distribution calculated
  type(mesh_geometry),  intent(in)      :: geom

  type(diff) :: adiff
  adiff%fast = .true. 

   ! calculate rupture time across whole mesh 
  call timeonevsall2d(amesh,src%velocity,src%rupt_time,geom%nton,adiff)

  return
end subroutine calc_rupt_front
!#######################################################################################

end module source_model 