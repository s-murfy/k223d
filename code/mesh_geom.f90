module mesh_geom 

use generic
use LAT_mesh
use LAT_source
use LAT_time
use distance 
use, intrinsic :: iso_fortran_env, only: error_unit

implicit none 

type mesh_geometry
  real(pr)                                    :: W
  type(reflect)                               :: surface       ! free surface boundary
  real(pr),dimension(:), allocatable          :: area,mz
  type(containerc), dimension(:),allocatable  :: ntoc 
  type(containern), dimension(:),allocatable  :: nton
  logical,dimension(:), allocatable           :: is_boundary
  real(pr),dimension(:),allocatable           :: dist_to_border
  real(pr),dimension(:,:),allocatable         :: dist
end type mesh_geometry

contains

!==============================================================================
subroutine build_mesh_geometry(amesh, surface, geom)
    type(mesh),          intent(in)    :: amesh
    type(reflect),       intent(in)    :: surface
    type(mesh_geometry), intent(out)   :: geom

    geom%surface = surface  !assign surface node information to geom
    !allocate geom arrays 
    call allocate_geom(amesh,geom)

    ! Geometry - rough estimate of fault width, area and average depth per cell
    call calc_cell_geometry(amesh, geom)

    ! Connectivity: node to cell and node to node connectivity — needed by everything below
    call calc_mesh_connectivity(amesh, geom%ntoc, geom%nton)

    ! Boundary mask — needs nton and surface
    call find_boundary_nodes(amesh, geom%nton, geom%is_boundary, surface)

    write(error_unit,*) 'calculate distance to border'
    ! Distance to border — needs ntoc, nton, is_boundary; accounting for plane wave from boundary
    call calc_dist_to_border(amesh, geom%ntoc, geom%nton, geom%is_boundary, geom%dist_to_border)

    write(error_unit,*) 'calculate distance between all nodes'
    ! Full distance matrix — most expensive, placed last
    call calc_all_distances(amesh, geom%ntoc, geom%nton, geom%dist)

end subroutine build_mesh_geometry
!==============================================================================
subroutine allocate_geom(amesh,geom)
type(mesh),          intent(in)    :: amesh
type(mesh_geometry), intent(inout) :: geom

allocate(geom%area(amesh%Ncells))
allocate(geom%mz(amesh%Ncells))
allocate(geom%ntoc(amesh%Nnodes))
allocate(geom%nton(amesh%Nnodes))
allocate(geom%is_boundary(amesh%Nnodes))
allocate(geom%dist_to_border(amesh%Nnodes))
allocate(geom%dist(amesh%Nnodes, amesh%Nnodes))

end subroutine allocate_geom
!==============================================================================
subroutine calc_cell_geometry(amesh, geom)
!==============================================================================
! Computes per-cell geometric properties and estimates fault width.
! Results are stored directly in the mesh_geometry derived type.
!
! For each triangular cell, the cross product of two edge vectors gives the
! cell normal vector n. From n, three quantities are derived in a single pass:
!
!   area(i)  = 0.5 * |n|         (cell area in mesh coordinate units^2)
!   mz(i)    = mean vertex depth  (mean depth of cell in mesh coordinate units)
!   sin(dip) = sqrt(nx^2+ny^2)/|n| (avoids acos, sign-independent of pz convention)
!
! The area-weighted mean sin(dip) over all cells is used with the mesh depth
! range to estimate the fault down-dip width W:
!
!   W = (max(pz) - min(pz)) / mean_sin(dip)
!
! This estimate is valid for planar to moderately curved fault surfaces.
! For nearly horizontal faults (mean_sin(dip) < 0.01) the estimate is
! unreliable and a fallback of sqrt(total_area) is used with a warning.
!
! Note: the sign convention of pz (positive or negative downward) does not
! affect the result — sin(dip) depends only on horizontal normal components
! and the depth range is computed as an absolute difference.
!
! W is stored in geom%W for use in source_model when computing the physical
! maximum sub-event radius: r_max_phys = model%rmax * geom%W
!==============================================================================
    type(mesh),          intent(in)    :: amesh
    type(mesh_geometry), intent(inout) :: geom
! local variables 
integer :: i
real(pr),dimension(3) :: x,y,z,vec_ab,vec_ac,n
real(pr)              :: sin_dip,n_mag
real(pr)              :: sum_sin_dip, total_area, mean_sin_dip

sum_sin_dip = 0._pr
total_area  = 0._pr
do i = 1,amesh%Ncells
  x = amesh%px(amesh%cell(i,:))
  y = amesh%py(amesh%cell(i,:))
  z = amesh%pz(amesh%cell(i,:))
  vec_ab = (/ x(2)-x(1), y(2)-y(1), z(2)-z(1) /)
  vec_ac = (/ x(3)-x(1), y(3)-y(1), z(3)-z(1) /)

  n(1) = vec_ab(2)*vec_ac(3) - vec_ab(3)*vec_ac(2)
  n(2) = vec_ab(3)*vec_ac(1) - vec_ab(1)*vec_ac(3)
  n(3) = vec_ab(1)*vec_ac(2) - vec_ab(2)*vec_ac(1)
  n_mag = norm2(n)
  geom%area(i) = 0.5_pr*n_mag
  geom%mz(i) = (z(1) + z(2) + z(3)) / 3._pr

  sin_dip = sqrt(n(1)**2+n(2)**2) /n_mag
  sum_sin_dip = sum_sin_dip + sin_dip * geom%area(i)
  total_area  = total_area  + geom%area(i)
enddo
  ! Fault width from depth range and area-weighted mean sin(dip)
mean_sin_dip = sum_sin_dip / total_area
if (mean_sin_dip < 0.01_pr) then
  write(*,*) 'WARNING: fault nearly horizontal, fault assumed to be square'
    geom%W = sqrt(total_area)  ! fallback: assume fault is square
else
    geom%W = (maxval(amesh%pz) - minval(amesh%pz)) / mean_sin_dip
endif

write(error_unit,'(a, es12.4,a)') 'Total fault area:      ', total_area,' m^2'
write(error_unit,'(a, f10.2,a)')   'Estimated fault width: ', geom%W,' m'

end subroutine calc_cell_geometry
!==============================================================================
subroutine calc_mesh_connectivity(amesh,ntoc,nton)
   use lists
   type(mesh), intent(in) :: amesh ! the mesh structure
   type(containerc), dimension(amesh%Nnodes),intent(out) :: ntoc ! the node to cell array (sparse matrix)
   type(containern), dimension(amesh%Nnodes),intent(out) :: nton ! the node to node array (sparse matrix)

! computing the node-to-cell array ntoc
   call compntoc(amesh,ntoc)
! computing the node-to-node array nton
   call compnton(amesh,nton)

end subroutine calc_mesh_connectivity
!==============================================================================
subroutine find_boundary_nodes(amesh, nton, is_boundary,surface)
  use lists 
  type(mesh), intent(in) :: amesh
  type(containern), dimension(amesh%Nnodes),intent(in) :: nton ! the node to node array (sparse matrix)
  type(reflect), intent(in) :: surface
  ! logical, dimension(amesh%Nnodes), intent(out) :: is_boundary
  logical, dimension(:), allocatable :: is_boundary

  type(liste), pointer :: pcur
  logical, dimension(amesh%Nnodes) :: is_surface
  integer(pin) :: i


  is_boundary = .false.
  is_surface  = .false.

  if (surface%present) then
    do i = 1, surface%Nnodes
      is_surface(surface%nodes(i)) = .true.
    enddo
  endif

  do i = 1, amesh%Nnodes
    pcur => nton(i)%ptr
    do while (associated(pcur))
      if (pcur%ref%cellonedge(2) == 0) then
        if (.not. is_surface(i))            is_boundary(i)            = .true.
        if (.not. is_surface(pcur%idnode))  is_boundary(pcur%idnode)  = .true.
      endif
      pcur => pcur%next
    enddo
  enddo

end subroutine find_boundary_nodes
!==============================================================================
subroutine calc_dist_to_border(amesh, ntoc, nton, is_boundary, dist_to_border)
  use LAT_time
  use distance

! input variables
  type(mesh),intent(in) :: amesh
  logical, dimension(amesh%Nnodes), intent(in) :: is_boundary
  type(containerc), dimension(amesh%Nnodes),intent(in) :: ntoc ! the node to cell array (sparse matrix)
  type(containern), dimension(amesh%Nnodes),intent(in) :: nton ! the node to node array (sparse matrix)
  real(pr), dimension(amesh%Nnodes), intent(out) :: dist_to_border

! local variables 
  integer :: i
  type(diff) :: adiff ! the diffraction structure
  
  ! assume no diffraction points 
  adiff%fast = .true.
  adiff%NdiffNodes = 0

  !initialise array
  dist_to_border=infinity
  ! the boundary is set to 100000 and this is finally subtracted from the final values
  ! for the distance to the border. This is done in order to mimic a plane wave
  where (is_boundary) dist_to_border=100000._pr   

  call onevsall2d(amesh,dist_to_border,ntoc,nton,adiff)
 
  dist_to_border = dist_to_border-100000._pr      !
  
  return
  end subroutine calc_dist_to_border
!==============================================================================
subroutine calc_all_distances(amesh,ntoc,nton,dist)
   use LAT_time 
! input variables
   type(mesh),intent(in) :: amesh ! the mesh structure
   type(containerc), dimension(amesh%Nnodes),intent(in) :: ntoc ! the node to cell array (sparse matrix)
   type(containern), dimension(amesh%Nnodes),intent(in) :: nton ! the node to node array (sparse matrix)
   real(pr), dimension(amesh%Nnodes,amesh%Nnodes),intent(out) :: dist ! the distance matrix node to node

! local variables
   type(diff) :: adiff ! the diffraction structure
   integer(pin) :: i, k ! loop index
   real(pr), dimension(amesh%Nnodes) :: distarray ! the temporary distance array

  ! assume no diffraction points 
  adiff%fast = .true.
  adiff%NdiffNodes = 0

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

end subroutine calc_all_distances
!==============================================================================
end module mesh_geom 