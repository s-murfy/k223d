module utils

Integer, Parameter :: kdp = selected_real_kind(15)
public :: find,utm_geo,sort_rows,erasedupli
private :: kdp

private :: find_1dI,find_1dR, find_2dR ,find_2dI

interface find
         module procedure find_1dI,find_1dR, find_2dR,find_2dI
end interface


contains

!=================================================
! arrange array 'a' in ascending order
subroutine sort_rows(a,n,m)
implicit none
integer,intent(in) :: n,m
integer, dimension(n,m),intent(inout) :: a
integer, dimension(n,m) :: a_dummy
integer, dimension(n) :: idx
integer :: i


a_dummy = a

do i = 1,n
! descending order
  idx(i) = minloc(a_dummy(:,1),1)
  a_dummy(idx(i),1) = maxval(a_dummy(:,1),1)+1
enddo

a = a(idx,:)

end subroutine sort_rows

!========================================================
subroutine find_1dR(array,condt,target_value,idx)
implicit none
real,dimension(:),intent(in) :: array
integer,allocatable,dimension(:), intent(out) :: idx
integer,allocatable,dimension(:) :: loc
character(2),intent(in) :: condt
integer :: i,n,target_value,inx

n = size(array)
allocate(loc(n))
loc = 0
inx = 0
if (condt == '==') then
       do i = 1,n
              if(array(i) == target_value) then
                     inx = inx+1
                     loc(inx) = i
              endif
       enddo
else
       do i = 1,n
              if(array(i) /= target_value) then
                     inx = inx+1
                     loc(inx) = i
              endif
       enddo
endif


if (inx > 0) then  ! target_value is in array
	allocate(idx(inx))
	idx = loc(1:inx)
else   ! target_value is not in array
	allocate(idx(1))
	idx = -1

endif
deallocate(loc)

return

end subroutine find_1dR
!========================================================
!========================================================
subroutine find_1dI(array,condt,target_value,idx)
implicit none
integer,dimension(:),intent(in) :: array
integer,allocatable,dimension(:), intent(out) :: idx
integer,allocatable,dimension(:) :: loc
integer :: i,n,target_value,inx
character(2),intent(in) :: condt

n = size(array,1)
allocate(loc(n))
loc = 0
inx = 0

if (condt == '==') then
       do i = 1,n
              if(array(i) == target_value) then
                     inx = inx+1
                     loc(inx) = i
              endif
       enddo
else
       do i = 1,n
              if(array(i) /= target_value) then
                     inx = inx+1
                     loc(inx) = i
              endif
       enddo
endif


if (inx > 0) then  ! target_value is in array
	allocate(idx(inx))
	idx = loc(1:inx)
else   ! target_value is not in array
	allocate(idx(1))
	idx = -1 !loc(1)

endif
deallocate(loc)

return

end subroutine find_1dI
!========================================================
subroutine find_2dR(array,condt,target_value,idx,jdx)
implicit none
real,dimension(:,:),intent(in) :: array
integer,allocatable,dimension(:), intent(out) :: idx,jdx
integer,allocatable,dimension(:) :: loc_i,loc_j
real :: target_value
integer :: i,j,n,m,inx
character(2),intent(in) :: condt

n = size(array,1)
m = size(array,2)

allocate(loc_i(n*m),loc_j(n*m))

loc_i(:) = 0
loc_j(:) = 0

inx = 0
if (condt == '==') then
       do j = 1,m
              do i = 1,n
                     if(array(i,j) == target_value) then
                            inx = inx+1
                            loc_i(inx) = i
                            loc_j(inx) = j
                     endif
              enddo
       enddo
else
       do j = 1,m
              do i = 1,n
                     if(array(i,j) /= target_value) then
                            inx = inx+1
                            loc_i(inx) = i
                            loc_j(inx) = j
                     endif
              enddo
       enddo
endif

if (inx > 0) then  ! target_value is in array
	allocate(idx(inx))
	idx = loc_i(1:inx)

	allocate(jdx(inx))
	jdx = loc_j(1:inx)

else   ! target_value is not in array
	allocate(idx(1))
	idx = -1

	allocate(jdx(1))
	jdx = -1

endif
deallocate(loc_i,loc_j)

return
end subroutine find_2dR
!========================================================
!========================================================
subroutine find_2dI(array,condt,target_value,idx,jdx)
implicit none
integer,dimension(:,:),intent(in) :: array
integer,allocatable,dimension(:), intent(out) :: idx,jdx
integer,allocatable,dimension(:) :: loc_i,loc_j
integer :: i,j,n,m,target_value,inx
character(2),intent(in) :: condt

n = size(array,1)
m = size(array,2)

allocate(loc_i(n*m),loc_j(n*m))

loc_i(:) = 0
loc_j(:) = 0

inx = 0
if (condt == '==') then
       do j = 1,m
              do i = 1,n
                     if(array(i,j) == target_value) then
                            inx = inx+1
                            loc_i(inx) = i
                            loc_j(inx) = j
                     endif
              enddo
       enddo
else
       do j = 1,m
              do i = 1,n
                     if(array(i,j) /= target_value) then
                            inx = inx+1
                            loc_i(inx) = i
                            loc_j(inx) = j
                     endif
              enddo
       enddo
endif


if (inx > 0) then  ! target_value is in array
	allocate(idx(inx))
	idx = loc_i(1:inx)

	allocate(jdx(inx))
	jdx = loc_j(1:inx)

else   ! target_value is not in array
	allocate(idx(1))
	idx = -1 !

	allocate(jdx(1))
	jdx = -1

endif
deallocate(loc_i,loc_j)

return
end subroutine find_2dI
!========================================================
subroutine erasedupli(val,nv,nnv)
! erase duplicates from vector val
   integer :: nv,nnv
   integer, dimension(nv) :: val
   logical, dimension(nv) :: dp
   integer :: i,j

   dp=.false.
   do i=1,nv-1
      if (val(i) == 0) dp(i)=.true.
      if (dp(i)) cycle
      do j=i+1,nv
         if (.not.dp(j)) dp(j)=(val(i) == val(j))
      enddo
   enddo

   i=1
   nnv=nv
   do while (i <= nnv)
      if (dp(i)) then
         do j=i,nnv-1
            val(j)=val(j+1)
            dp(j)=dp(j+1)
         enddo
         nnv=nnv-1
      else
         i=i+1
      endif
   enddo

   return
end subroutine erasedupli
!========================================================

end module utils
