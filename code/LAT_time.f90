module LAT_time
  use generic
  use LAT_mesh
  implicit none

  type diff
    logical :: fast
    integer(pin) :: Nnodes
    integer(pin),allocatable,dimension(:) :: nodes
  end type diff

  type interface   ! replace reflect 
    logical :: check = .false.  ! can we do this (preset value?)
    integer(pin) :: nnodes 
    integer(pin),allocatable,dimension(:) :: nodes
  end type interface 

  type reflect
    logical :: present
    integer(pin) :: Nnodes
    integer(pin),allocatable,dimension(:) :: nodes
  end type reflect

  real(pr), allocatable, dimension(:) :: time
  real(pr), allocatable, dimension(:) :: kappa
  integer(pin), allocatable, dimension(:) :: mode

  real(pr), allocatable, dimension(:) :: velocity

! set up for debugging solvers 
  real(pr), allocatable, dimension(:) :: sface
  real(pr), allocatable, dimension(:) :: shead
  real(pr), allocatable, dimension(:) :: splane
  real(pr), allocatable, dimension(:) :: sedge
  real(pr), allocatable, dimension(:) :: scplane
  real(pr), allocatable, dimension(:) :: origidnode

end module LAT_time
