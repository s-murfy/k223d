module generic

implicit none
integer,parameter :: int4 = selected_int_kind(9)
integer,parameter :: int8 = selected_int_kind(18)
integer,parameter :: real4 = selected_real_kind(6,37) ! 6 significant digits of precision with a exponential range of a least 37
integer,parameter :: real8 = selected_real_kind(15,307)
integer,parameter :: pin = int4
integer,parameter :: pr = real8
integer,parameter :: psav = real4
real(pr), parameter :: infinity=1.e32_pr
! real(pr), parameter :: water_level= 10**-precision(1) ! 10**-12._pr !1000._pr*epsilon(10.**-20._pr)
! real(pr), parameter :: water_level=10._pr*epsilon(1._pr)

integer(pin), parameter :: nimp=100
logical :: verbose=.false.

contains
!###############################################################################
function water_level(x)
real(pr) :: x ,water_level

  water_level =  10._pr**(-1._pr*float(precision(x)-5))

end function water_level
!###############################################################################
subroutine printperc(a,b)
  integer(pin) :: a,b

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


end module generic