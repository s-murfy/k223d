module LAT_mesh
use generic
implicit none

type mesh
  integer(pin) :: Nnodes,Ncells! Nnodes : no of nodes ; Ncells : number of cells 
  real(pr), dimension(:), allocatable :: px,py,pz  !          Px,Py,Pz : x,y,z location of nodes (Nnodes)
  integer(pin), dimension(:,:), allocatable :: cell !             cell : relationship between cell and nodes (Ncells,3)
                                                    !                    The row gives the cell and the three columns refer 
                                                    !                    to node indexes
  integer(pin),allocatable,dimension(:,:) ::  EToE  !             EToE : element to element connectivity     (Ncells,3)
                                                    !                    This array gives the cell indexes of cells with conjoining faces  
                                                    !                    Row is the query cell, and the three columns are 
                                                    !                    the indexes to the cells that share a face with the query cell
                                                    !                    If column index = row there is no adjointing cell (i.e. EToE(i,3) = i)
  integer(pin), allocatable, dimension(:,:,:) ::   FVsNodes ! FVsNodes : array linking faces (pairs of nodes) with cells (Nnodes,Nodes,2)
                                                    !                    Sparse matrix where the first two indexes give the node index corresponding
                                                    !                    to one face between cells. The two values provided in the 3rd index are the 
                                                    !                    cell indexes corresponding to the face.
                                                    !                    Null values are set to 0 
end type mesh     

contains



end module LAT_mesh
