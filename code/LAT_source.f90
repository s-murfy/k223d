module LAT_source

use generic
use LAT_mesh
implicit none


type model_param
    integer(pin) :: na
    integer(pin) :: gauss_no
    real(pr) :: mu, rmin, rmax, param
    real(pr) :: sd,mw,target_area,moment,target_width
    real(pr) :: actual_area,actual_width,actual_length
    character(3) :: param_type
    logical :: defined_area
end type model_param

type surface_type  ! need to remove reference in typedef.f90 -> read_input 
    real(pr) :: z0
    character(30) :: filename
    logical :: use_file,exists
    integer(pin)  :: NoSurfNodes                         ! no. of surface nodes
    integer(pin),allocatable,dimension(:)  :: SurfNodes  ! nodes that are designated along the surface
end type surface_type

type pdfinputs
    character(30) :: pdf_fname
    character(7) :: pdf_type
    integer(pin) :: ng
! pointer to a Node number of the Quake
    integer(pin), dimension(:), allocatable :: vertexcenter
    real(pr), dimension(:,:), allocatable :: sizeandhigh
    real(pr), dimension(:), allocatable ::  g_pdf
end type pdfinputs

type source 
    integer(pin) :: nuc_id   !cell number for where earthquake will initiate 
    integer(pin) :: QuakeElemNo,QuakeNodesNo ! Number of Elements and nodes for an earthquake
    integer(pin),allocatable,dimension(:) :: QuakeElem,QuakeNodes ! Array of elements and nodes that make up slipping zone
    integer(pin) :: QuakeBorder_NodesNo, QuakeBorder_ElemNo ! Number of elements and nodes forming the slipping zone border
    integer(pin),allocatable,dimension(:) :: QuakeBorder_Nodes, QuakeBorder_Elem ! Index of elements and nodes forming the slipping zone border
end type source 


! type mesh
! integer :: Nnodes,Ncells             ! Number of Nodes ; Number of element
! integer :: QuakeElemNo,QuakeNodesNo ! Number of Elements and nodes for an earthquake
! real,allocatable,dimension(:,:) :: dist ! distance matrix
! integer,allocatable,dimension(:,:) ::  cell
! integer,allocatable,dimension(:,:) ::  EToE  ! Elements to Elemeents
! integer,allocatable,dimension(:) :: QuakeElem,QuakeNodes ! Array of elements and nodes that make up slipping zone
! integer :: QuakeBorder_NodesNo, QuakeBorder_ElemNo ! Number of elements and nodes forming the slipping zone border
! integer,allocatable,dimension(:) :: QuakeBorder_Nodes, QuakeBorder_Elem ! Index of elements and nodes forming the slipping zone border
! real,allocatable,dimension(:) ::Dist2Border ! distance matrix to boundary of fault
! real,allocatable,dimension(:) :: px,py,pz,mz,area,slip,slip_prim,slip_sec ! Coordinates of the nodes; area of each element; slip in each element
! integer  :: NoSurfNodes                         ! no. of surface nodes
! integer,allocatable,dimension(:)  :: SurfNodes  ! nodes that are designated along the surface
! end type mesh


! parameters for slip 
real(pr), allocatable, dimension(:) :: slip,slip_prim,slip_sec
! real(pr), allocatable, dimension(:,:) :: dist
real(pr), allocatable, dimension(:) :: Dist2Border
! real(pr), allocatable, dimension(:) :: velocity
real(pr), allocatable, dimension(:) :: mz  ! array of mean depth of the each cell [whole mesh] 
real(pr), allocatable, dimension(:) :: area ! array of area of each cell [whole mesh] 


end module LAT_source 
