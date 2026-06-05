module LAT_distance
  use generic
  use LAT_mesh
  implicit none

type diff
    logical :: fast
    integer(pin) :: NdiffNodes
    integer(pin),allocatable,dimension(:) :: nodes
  end type diff

end module LAT_distance
