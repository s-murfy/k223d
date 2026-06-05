module LAT_source

use generic
use LAT_mesh
implicit none


type model_param
    integer(pin) :: na
    real(pr) :: mu, rmin, rmax, param
    real(pr) :: sd,mw
end type model_param

type surface_type  ! need to remove reference in typedef.f90 -> read_input 
    real(pr) :: z0
    character(30) :: filename
    logical :: use_file,exists
    integer(pin)  :: NoSurfNodes                         ! no. of surface nodes
    integer(pin),allocatable,dimension(:)  :: SurfNodes  ! nodes that are designated along the surface
end type surface_type

type pdfinputs
    logical :: is_uniform
    real(pr), dimension(:), allocatable :: distrib
end type pdfinputs


! type mesh_geometry
!     real(pr),       allocatable, dimension(:)   :: area  ! area of each cell
!     real(pr),       allocatable, dimension(:)   :: mz    ! mean depth of each cell
!     logical,        allocatable, dimension(:)   :: is_boundary ! mask marking nodes on mesh edge
!     type(containerc), allocatable, dimension(:) :: ntoc   ! sparse list of node to cell 
!     type(containern), allocatable, dimension(:) :: nton   ! sparse list of node to node 
!     real(pr),       allocatable, dimension(:,:) :: dist  ! distance between all nodes 
!     real(pr),       allocatable, dimension(:)   :: dist_to_border ! distance from each node to mesh edge 
! end type mesh_geometry

type source
    integer(pin)                                :: nuc_id  
    real(pr),       allocatable, dimension(:)   :: slip
    ! real(pr),       allocatable, dimension(:)   :: slip_prim
    ! real(pr),       allocatable, dimension(:)   :: slip_sec
    real(pr),       allocatable, dimension(:)   :: velocity
    real(pr),       allocatable, dimension(:)   :: rupt_time
end type source


! type source 
!     integer(pin) :: nuc_id   !cell number for where earthquake will initiate 
!     integer(pin) :: QuakeElemNo,QuakeNodesNo ! Number of Elements and nodes for an earthquake
!     integer(pin),allocatable,dimension(:) :: QuakeElem,QuakeNodes ! Array of elements and nodes that make up slipping zone
!     integer(pin) :: QuakeBorder_NodesNo, QuakeBorder_ElemNo ! Number of elements and nodes forming the slipping zone border
!     integer(pin),allocatable,dimension(:) :: QuakeBorder_Nodes, QuakeBorder_Elem ! Index of elements and nodes forming the slipping zone border
! end type source 



! parameters for slip 
real(pr), allocatable, dimension(:) :: slip,slip_prim,slip_sec
! real(pr), allocatable, dimension(:,:) :: dist
real(pr), allocatable, dimension(:) :: Dist2Border
! real(pr), allocatable, dimension(:) :: velocity
real(pr), allocatable, dimension(:) :: mz  ! array of mean depth of the each cell [whole mesh] 
real(pr), allocatable, dimension(:) :: area ! array of area of each cell [whole mesh] 


end module LAT_source 
