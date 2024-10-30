module LAT_mesh_util
use generic
use LAT_mesh
implicit none

!To DO: set check for the end of input vtk file properly 

interface find
         module procedure find_1dI,find_1dR,find_2dR,find_2dI
end interface

interface dumpnodeattributevtk
    module procedure idumpnodeattributevtk
    module procedure rdumpnodeattributevtk
end interface 



contains


!###############################################################################
subroutine read_mesh_file(amesh,surface,nuc_id)
  ! subroutine read_mesh_file(amesh,adiff)
  use LAT_time
  use LAT_source 

  integer(pin),allocatable,dimension(:) :: cell_id
  integer(pin),intent(out) :: nuc_id 
  type(mesh) :: amesh
  ! type(diff) :: adiff
  type(reflect) :: surface


  character(30) :: mesh_file

  ! read unstructured mesh from file
  write(*,*) 'enter read_vtk_mesh'
  mesh_file = 'input.vtk'
  ! call read_vtk_mesh_jup(amesh,mesh_file,adiff)
  call read_vtk_mesh_jup_surf(amesh,mesh_file,surface,nuc_id)
  write(*,*) 'exit read_vtk_mesh'

! calculate inter-connectivity of Elemeents 

  call mesh_connectivity(amesh)

  return
end subroutine read_mesh_file
!==============================================================================
subroutine mesh_connectivity(amesh)
  !  Create a matrix that defines neighbouring cells
  !  Dimension :  (No of cells , 3)
  !  Row relates to cells
  !  Columns give the cell numbers of the 3 cells with adjoining faces to current cell
  !  If a cell has less than 3 neighbours the index of the current cell is placed in the column
  type(mesh),intent (inout) :: amesh
  integer(pin), dimension(3) :: node_id,neighbour_element
  integer(pin),allocatable,dimension(:) :: common_node
  integer(pin) :: i,j,k,Edge_no,noedges
  integer(pin), dimension(2) :: test_face
  integer(pin), allocatable, dimension(:,:) ::   Face_array,dummy_Edge
  logical :: new_face

  integer(pin),dimension(3,2) :: cell_side
  integer(pin),dimension(2) :: face
  integer :: Faces,NFaces, id
  integer(pin), allocatable, dimension(:,:) ::  FtoNode,FtoNodeT,FtoF
  integer(pin), allocatable, dimension(:,:,:) ::   FVsNodes
  integer(pin), dimension(3,2):: Edge
  integer(pin) :: cell1,cell2
  integer(pin), dimension(2) :: vertices

  Faces = 3

  write(*,*) 'starting allocation'
  !call flush()

  allocate(amesh%EToE(amesh%Ncells,Faces))
  !allocate(FVsNodes(amesh%Nnodes,amesh%Nnodes,2))  ! interger
  allocate(amesh%FVsNodes(amesh%Nnodes,amesh%Nnodes,2))  ! interger


  Edge(1,:) = (/ 1, 2/)
  Edge(2,:) = (/ 2, 3/)
  Edge(3,:) = (/ 3, 1/)


  write(*,*) 'create FVsNodes'
!  call flush()
!  FVsNodes = 0
  amesh%FVsNodes = 0


  do i = 1, amesh%Ncells
    do j = 1,Faces
       vertices = amesh%cell(i,Edge(j,:))    ! (Ncells, 3)
       do k = 1,2
          if (amesh%FVsNodes(vertices(3-k),vertices(k),1) == 0) then
            amesh%FVsNodes(vertices(3-k),vertices(k),1) = i
          elseif (amesh%FVsNodes(vertices(3-k),vertices(k),2) == 0) then
            amesh%FVsNodes(vertices(3-k),vertices(k),2) = i
          endif
       enddo
    enddo
  enddo

  do i =1,amesh%Ncells
      amesh%EToE(i,:) = i    ! null element for this matrix is identity
  enddo
  write(*,*) 'start EToE'
!   do j = 1,amesh%nnodes-1
!     do k = j+1,amesh%nnodes
!            cell1 = amesh%FVsNodes(j,k,1)
!            cell2 = amesh%FVsNodes(j,k,2)
!            if ((cell1 /= 0).and.(cell2 /= 0 ))then
!              do i = 1,3
!                if (amesh%EToE(cell1,i) == cell1) then
!                    amesh%EToE(cell1,i) = cell2
!                   !  exit
!                endif
!                if (amesh%EToE(cell2,i) == cell2) then
!                    amesh%EToE(cell2,i) = cell1
!                   !  exit
!                endif
!              enddo
!            endif
!      enddo
!  enddo



  do j = 1,amesh%nnodes-1
     do k = j+1,amesh%nnodes
            cell1 = amesh%FVsNodes(j,k,1)
            cell2 = amesh%FVsNodes(j,k,2)
            if ((cell1 /= 0).and.(cell2 /= 0 ))then
              do i = 1,3
                if (amesh%EToE(cell1,i) == cell1) then
                    amesh%EToE(cell1,i) = cell2
                    exit
                endif
              enddo
              do i = 1,3
                if (amesh%EToE(cell2,i) == cell2) then
                    amesh%EToE(cell2,i) = cell1
                    exit
                endif
              enddo
            endif
      enddo
  enddo

  return
end subroutine mesh_connectivity
!###############################################################################
subroutine read_vtk_mesh(amesh,mesh_file)

  type(mesh) :: amesh
  character(30),intent(in) :: mesh_file
  integer(pin),allocatable,dimension(:) :: layer_bc
  integer(pin) :: i, id,nsize,c1,c2,c3,element_type
  integer(pin) :: max_cells,t_cells,t_start,nnodes
  integer(pin) :: ipos_start,ipos_end
  logical :: tri_cells_start
  character(len=60) :: line

  open(13,file=trim(adjustl(mesh_file)),form = 'formatted', action = 'read')
  do i = 1,4
         read(13,*)
  enddo
  read(13,'(a)') line
  ipos_start = scan(line,"S",back=.true.)
  ipos_end = scan(line,"d",back=.false.)
!  write(*,*) ipos_start,ipos_end
  read (line(1+ipos_start:ipos_end),*) nnodes
!  amesh%Nnodes = nnodes-1  !NOTE: REMOVING NODE AT ZERO AS IT IS NOT USED
  amesh%Nnodes = nnodes  !
  write(*,*) 'No. of nodes....',amesh%Nnodes
  allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))
  amesh%px = 0._pr
  amesh%py = 0._pr
  amesh%pz = 0._pr
!  read(13,*)  !NOTE: REMOVING NODE AT ZERO AS IT IS NOT USED
  do i = 1,amesh%Nnodes
     read(13,*) amesh%px(i),amesh%py(i),amesh%pz(i)
!     print*,amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  read(13,*)
  read(13,'(a)')  line
!  write(*,'(a)') line
  ipos_start = scan(line,"CELLS",back=.true.)
!  print*,'ipos ........',ipos_start
  read (line(1+ipos_start:),*) max_cells,nsize
!  read(13,'(5x,i6,i7)')amesh%Ncells,nsize
!  write(*,*) max_cells,nsize
!  write(*,*)'no. of cells',max_cells
!  stop
  t_cells = 0
  tri_cells_start = .false.
  do i = 1,max_cells
      read(13,*) id
      if (id == 3) then
         if (.not.tri_cells_start) then
           t_start = i
            tri_cells_start = .true.
         endif
         t_cells = t_cells+1
       endif
  enddo

!  write(*,*) 'first triangle at.....',t_start
!  write(*,*)'no. of triangular cells.....',t_cells
  rewind(13)
  do i = 1,6+amesh%Nnodes+t_start
!  do i = 1,6+amesh%Nnodes+t_start+1 !add one for the node you jumped
      read(13,*)
  enddo
!  write(*,*) 'current position......',6+amesh%Nnodes+t_start
  amesh%Ncells = t_cells
  allocate(amesh%cell(amesh%Ncells,3))
  write(*,*)  'no. of cells........ ',amesh%Ncells
  do i = 1,amesh%Ncells
   read(13,*) id,c1,c2,c3
   amesh%cell(i,:) = (/ c1+1, c2+1,c3+1 /)
!   amesh%cell(i,:) = (/ c1, c2,c3 /)  !this is because we have removed the first node at the origin

  enddo
  close(13)

return
end subroutine read_vtk_mesh
!###############################################################################
!###############################################################################
subroutine read_vtk_mesh_jup(amesh,mesh_file,surf)
   use LAT_time

  type(mesh) :: amesh
  type(reflect) :: surf
  type(diff) :: adiff


  character(30),intent(in) :: mesh_file
  integer(pin),allocatable,dimension(:) :: layer_bc,list_diff
!  real(pr),allocatable,dimension(:) :: velocity,time
  integer(pin) :: i, id,nsize,c1,c2,c3,element_type
  integer(pin) :: max_cells,t_cells,t_start,nnodes
  integer(pin) :: ipos_start,ipos_end,counter,diff_switch
  real(pr) :: itime
  logical :: tri_cells_start
  character(len=60) :: line

  integer(pin) :: diff_level
  character(100) :: char2int



  open(13,file=trim(adjustl(mesh_file)),form = 'formatted', action = 'read')
  do i = 1,4
         read(13,*)
  enddo
  read(13,'(a)') line
  ipos_start = scan(line,"S",back=.true.)
  ipos_end = scan(line,"d",back=.false.)
!  write(*,*) ipos_start,ipos_end
  read (line(1+ipos_start:ipos_end),*) nnodes
!  amesh%Nnodes = nnodes-1  !NOTE: REMOVING NODE AT ZERO AS IT IS NOT USED
  amesh%Nnodes = nnodes  !
  write(*,*) 'No. of nodes....',amesh%Nnodes
  allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))
  amesh%px = 0._pr
  amesh%py = 0._pr
  amesh%pz = 0._pr
!  read(13,*)  !NOTE: REMOVING NODE AT ZERO AS IT IS NOT USED
  do i = 1,amesh%Nnodes
     read(13,*) amesh%px(i),amesh%py(i),amesh%pz(i)
!     print*,amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  read(13,*)
  read(13,'(a)')  line
!  write(*,'(a)') line
  ipos_start = scan(line,"CELLS",back=.true.)
  read (line(1+ipos_start:),*) max_cells,nsize
!  read(13,'(5x,i6,i7)')amesh%Ncells,nsize
!  write(*,*) max_cells,nsize
!  write(*,*)'no. of cells',max_cells
!  stop
  t_cells = 0
  tri_cells_start = .false.
  do i = 1,max_cells
      read(13,*) id
      if (id == 3) then
         if (.not.tri_cells_start) then
           t_start = i
            tri_cells_start = .true.
         endif
         t_cells = t_cells+1
       endif
  enddo

!  write(*,*) 'first triangle at.....',t_start
!  write(*,*)'no. of triangular cells.....',t_cells
  rewind(13)
  do i = 1,6+amesh%Nnodes+t_start
!  do i = 1,6+amesh%Nnodes+t_start+1 !add one for the node you jumped
      read(13,*)
  enddo
!  write(*,*) 'current position......',6+amesh%Nnodes+t_start
  amesh%Ncells = t_cells
  allocate(amesh%cell(amesh%Ncells,3))
  write(*,*)  'no. of cells........ ',amesh%Ncells
  do i = 1,amesh%Ncells
   read(13,*) id,c1,c2,c3
   amesh%cell(i,:) = (/ c1+1, c2+1,c3+1 /)
!   amesh%cell(i,:) = (/ c1, c2,c3 /)  !this is because we have removed the first node at the origin
  enddo
  read(13,*)     ! blank line


! next read velocity
  do i = 1,3
   read(13,*)     ! skip header informatin on CellEntityIds
  enddo
  allocate(velocity(amesh%Ncells))
  velocity = 0._pr
  do i = 1,amesh%Ncells
    read(13,*) velocity(i)
  enddo
  read(13,*)     ! blank line

! set up time array
   allocate(time(amesh%Nnodes))
   allocate(kappa(amesh%Nnodes))  ! inverse of curvature 
   allocate(mode(amesh%Nnodes))    ! operator type 
  !  kappa = infinity !0._pr
  !  mode = -1 
   time = infinity


! set up for debugging solvers 
  if (verbose) then 
    allocate(sface(amesh%Nnodes))
    allocate(shead(amesh%Nnodes))
    allocate(splane(amesh%Nnodes))
    allocate(sedge(amesh%Nnodes))
    allocate(scplane(amesh%Nnodes))
    allocate(origidnode(amesh%Nnodes))

    sface = infinity
    shead = infinity
    splane = infinity
    sedge = infinity
    scplane = infinity
  endif 
  !------ end of debugging---

   do i = 1,3
      read(13,*)     ! skip header information on initial time
   enddo

   do i = 1,amesh%Nnodes
      read(13,*) itime
      if (itime > -1.0) then
        time(i) = itime
        kappa(i) = 0._pr
        mode(i) = 0
      endif
   enddo

   if(command_argument_count().ne.1) then
     write(*,*)'Assuming no diffraction taking place'
     adiff%fast = .true.  ! turn off diffraction
   else
    call get_command_argument(1,char2int)
    read(char2int,*) diff_level
    if (diff_level== 0) then
        write(*,*) 'No diffraction points considered'
        adiff%fast = .true.  ! turn off diffraction
    elseif (diff_level== 1) then
        read(13,*)
        read(13,*)
        read(13,*)

        write(*,*) 'Defined diffraction points considered'
        adiff%fast = .false.  ! turn on diffraction
        allocate(list_diff(amesh%Nnodes))
        list_diff = 0
        counter = 0
        do i = 1,amesh%Nnodes
          read(13,*) diff_switch
          if (diff_switch > 0) then
              counter=counter+1
              list_diff(counter) = i
          endif
        enddo

        if (counter == 0) then
          write(*,*) 'Error: no diffraction points detected'
          write(*,*) 'Will run assuming no diffraction'
          adiff%fast = .true.  ! turn on diffraction
        endif
        adiff%Nnodes = counter
        allocate(adiff%nodes(adiff%Nnodes))
        adiff%nodes = list_diff(1:counter)
        deallocate(list_diff)
        write(*,*) 'no. of diffraction nodes:     ', adiff%Nnodes
        !HERE we need to
        !(1) read in the number of nodes that are nonzero
        !(2) test if number of nodes ==0  => set adiff%fast = .true. and flag it
        !(3) create an array with all non-zero diffraction nodes
    elseif (diff_level == 2) then
           write(*,*) 'All mesh points considered diffraction points'
           adiff%fast = .false.  ! turn on diffraction
           adiff%Nnodes = amesh%Nnodes
           allocate(adiff%nodes(adiff%Nnodes))
           do i = 1,adiff%Nnodes
             adiff%nodes(i) = i
           enddo
    endif
   endif
   !
   ! adiff%fast = .true.  ! turn off diffraction
   ! if (.not.adiff%fast) then
   !   write(*,*)'number of diffraction nodes.......',ibc
   !   adiff%Nnodes = ibc
   !   allocate(adiff%nodes(adiff%Nnodes))
   !   adiff%nodes = bc_nodes(1:ibc)
   ! else
   !   adiff%Nnodes = 1
   !   allocate(adiff%nodes(adiff%Nnodes))
   !   adiff%nodes = 1
   ! endif

  close(13)

return
end subroutine read_vtk_mesh_jup
!###############################################################################
!###############################################################################
subroutine read_vtk_mesh_jup_surf(amesh,mesh_file,surface,nuc_id)
  use LAT_time
  Use, intrinsic :: iso_fortran_env, Only : iostat_end

 type(mesh) :: amesh
 type(reflect) :: surface
!  type(diff) :: adiff

 integer(pin), intent(out) :: nuc_id
 character(30),intent(in) :: mesh_file
 integer(pin),allocatable,dimension(:) :: layer_bc,list_surf,nuc_nodes
!  real(pr),allocatable,dimension(:) :: velocity,time
 integer(pin) :: i, id,nsize,c1,c2,c3,element_type
 integer(pin) :: max_cells,t_cells,t_start,nnodes
 integer(pin) :: ipos_start,ipos_end,counter,diff_switch
 real(pr) :: itime,random
 logical :: tri_cells_start
 character(len=256) :: line
 integer(pin),dimension(3) :: vec1,vec2,vec3,vec4,vec5,vec6
 integer(pin),allocatable,dimension(:) :: n_cells
 integer(pin) :: diff_level,error,node_count,cell_count
 integer(pin) :: ios,ivar
 logical :: keyword_present
 


 open(13,file=trim(adjustl(mesh_file)),form = 'formatted', action = 'read')
 do i = 1,4
        read(13,*)
 enddo
 read(13,'(a)') line
 ipos_start = scan(line,"S",back=.true.)
 ipos_end = scan(line,"d",back=.false.)
!  write(*,*) ipos_start,ipos_end
 read (line(1+ipos_start:ipos_end),*) nnodes
!  amesh%Nnodes = nnodes-1  !NOTE: REMOVING NODE AT ZERO AS IT IS NOT USED
 amesh%Nnodes = nnodes  !
 write(*,*) 'No. of nodes....',amesh%Nnodes
 allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))
 amesh%px = 0._pr
 amesh%py = 0._pr
 amesh%pz = 0._pr
!  read(13,*)  !NOTE: REMOVING NODE AT ZERO AS IT IS NOT USED
 do i = 1,amesh%Nnodes
    read(13,*) amesh%px(i),amesh%py(i),amesh%pz(i)
 enddo

 read(13,*)
 read(13,'(a)')  line
!  write(*,'(a)') line
 ipos_start = scan(line,"CELLS",back=.true.)
 read (line(1+ipos_start:),*) max_cells,nsize
 t_cells = 0
 tri_cells_start = .false.
 do i = 1,max_cells
     read(13,*) id
     if (id == 3) then
        if (.not.tri_cells_start) then
          t_start = i
           tri_cells_start = .true.
        endif
        t_cells = t_cells+1
      endif
 enddo

!  write(*,*) 'first triangle at.....',t_start
!  write(*,*)'no. of triangular cells.....',t_cells
 rewind(13)
 do i = 1,6+amesh%Nnodes+t_start
!  do i = 1,6+amesh%Nnodes+t_start+1 !add one for the node you jumped
     read(13,*)
 enddo
!  write(*,*) 'current position......',6+amesh%Nnodes+t_start
 amesh%Ncells = t_cells
 allocate(amesh%cell(amesh%Ncells,3))
 write(*,*)  'no. of cells........ ',amesh%Ncells
 do i = 1,amesh%Ncells
  read(13,*) id,c1,c2,c3
  amesh%cell(i,:) = (/ c1+1, c2+1,c3+1 /)
!   amesh%cell(i,:) = (/ c1, c2,c3 /)  !this is because we have removed the first node at the origin
 enddo
 read(13,*)     ! blank line


! next read velocity
 do i = 1,3
  read(13,*)     ! skip header informatin on CellEntityIds
 enddo
 allocate(velocity(amesh%Ncells))
 velocity = 0._pr
 do i = 1,amesh%Ncells
   read(13,*) velocity(i)
   if (velocity(i) <= epsilon(0._pr)) then 
     write(*,*) 'Error: velocity with a value <= 0 detected'
     stop 
   endif 
 enddo
 read(13,*)     ! blank line

! set up time array
  allocate(time(amesh%Nnodes))
  allocate(kappa(amesh%Nnodes))  ! inverse of curvature 
  allocate(mode(amesh%Nnodes))    ! operator type 
 !  kappa = infinity !0._pr
 !  mode = -1 
  time = infinity

! set up for debugging solvers 
if(verbose) then 
  allocate(sface(amesh%Nnodes))
  allocate(shead(amesh%Nnodes))
  allocate(splane(amesh%Nnodes))
  allocate(sedge(amesh%Nnodes))
  allocate(scplane(amesh%Nnodes))
  allocate(origidnode(amesh%Nnodes))

  sface = infinity
  shead = infinity
  splane = infinity
  sedge = infinity
  scplane = infinity
endif 
 !------ end of debugging---

  do i = 1,3
     read(13,*)     ! skip header information on initial time
  enddo
  allocate(nuc_nodes(amesh%nnodes))
! define initial time 
  node_count = 0
  do i = 1,amesh%Nnodes
     read(13,*) itime
     if (itime > -1.0) then
       time(i) = itime
       node_count = node_count+1
       nuc_nodes(node_count) = i 
       kappa(i) = 0._pr
       mode(i) = 0
     endif
  enddo
! find starting cell 
  nuc_id = -1 

    if (node_count == 1) then ! randomly choose from cells connected with node 
       cell_count = 0
       allocate(n_cells(amesh%Ncells))
       n_cells = 0
       do i = 1,amesh%Ncells
        if ( any(amesh%cell(i,:) == nuc_nodes(1)))   then 
              cell_count = cell_count+1
              n_cells(cell_count) = i
        endif
      enddo

      if (cell_count==0) then
        write(*,*) 'ERROR: nucleation node does not correspond to grid'
        stop
      else
        call random_number(random)
        nuc_id = ceiling(random*float(cell_count))   !
        deallocate(n_cells)

        ! write(*,*) 'nuc_nodes(1)   ',nuc_nodes(1)
        ! write(*,*) 'cell_count    ',cell_count
        ! write(*,*)'n_cells    ',n_cells(1:7)
      endif 
      ! deallocate(nuc_nodes)

    elseif (node_count == 3) then! only one choice 
     write(*,*) 'number of nucleation nodes found...',node_count
     write(*,*) nuc_nodes(1:3)
     do i = 1,amesh%Ncells
         vec1 = (/ nuc_nodes(1), nuc_nodes(2), nuc_nodes(3) /)
         vec2 = (/ nuc_nodes(3), nuc_nodes(1), nuc_nodes(2) /)
         vec3 = (/ nuc_nodes(2), nuc_nodes(3), nuc_nodes(1) /)
         vec4 = (/ nuc_nodes(1), nuc_nodes(3), nuc_nodes(2) /)
         vec5 = (/ nuc_nodes(3), nuc_nodes(2), nuc_nodes(1) /)
         vec6 = (/ nuc_nodes(2), nuc_nodes(1), nuc_nodes(3) /)

         if ( all(amesh%cell(i,:) == vec1).or.  &
              all(amesh%cell(i,:) == vec2).or.  & 
              all(amesh%cell(i,:) == vec3).or.  &
              all(amesh%cell(i,:) == vec4).or.  &
              all(amesh%cell(i,:) == vec5).or.  &
              all(amesh%cell(i,:) == vec6))     then 
          nuc_id = i 
          exit
         endif
     enddo 
    else ! randomly decide starting cell (assuming counter == 0 )
      call random_number(random)
      nuc_id = ceiling(random*float(amesh%Ncells))   !
    endif 
  if(nuc_id == -1) write(*,*) 'ERROR: nucleation cell not defined'
  write(*,*) 'nucleation cell id:     ',nuc_id
  deallocate(nuc_nodes)
  

  call find_keyword(13,'surface',keyword_present)
  if (keyword_present) then 
    surface%present = .true.  ! turn on surface 
    write(*,*)'fault  surface found'
    read(13,'(a)',iostat=ios) line 

     !  read in surface nodes from vtk file  
    allocate(list_surf(amesh%Nnodes))
    list_surf = 0._pr
    counter = 0
    do i = 1,amesh%Nnodes
      read(13,*) ivar 
      if (ivar == 1) then 
          counter=counter+1
          list_surf(counter) = i
      endif 
    enddo

    surface%Nnodes = counter
    allocate(surface%nodes(surface%Nnodes))
    surface%nodes = list_surf(1:counter)
    deallocate(list_surf)
    write(*,*) 'no. of surface nodes:     ', surface%Nnodes
  else
    surface%present = .false.  ! turn on surface 
    write(*,*) 'end of file: no surface component detected '
    write(*,*)'Assuming no surface rupture'
  endif 


 close(13)

return
end subroutine read_vtk_mesh_jup_surf
!###############################################################################
subroutine find_keyword(ion,keyword,present)
integer(pin), intent(in) :: ion
character(len=*), intent(in) :: keyword
character(len=256) :: line                
logical, intent(out) :: present 
integer(pin) :: ios

present = .false.
rewind(ion)
do 
  read(ion,'(a)',iostat=ios)  line
  if (ios /=0) then 
    present = .false.  !not found in file 
    exit
  endif
  if (index(line, adjustl(keyword)) > 0) then 
    present = .true.  ! found in file 
    exit
  endif
enddo


return 
end subroutine find_keyword
!###############################################################################
subroutine dumpmeshvtk_jup(dev,amesh)

  type(mesh) :: amesh
  integer(pin) :: dev

  integer(pin) :: i,j

  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a8)') 'distance'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a25:)') 'DATASET UNSTRUCTURED_GRID'
  if( pr == 4) then
     write(dev,'(a7,i10,a6)') 'POINTS ',amesh%Nnodes,' float'
  else
     write(dev,'(a7,i10,a7)') 'POINTS ',amesh%Nnodes,' double'
  endif

  do i=1,amesh%Nnodes
     write(dev,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  write(dev,'(a6,2i10)') 'CELLS ',amesh%Ncells,amesh%Ncells*4
  do i=1,amesh%Ncells
     write(dev,'(a2,$)') '3 '
     write(dev,*) (amesh%cell(i,j)-1,j=1,3)
  enddo

  write(dev,*)
  write(dev,'(a11,i10)') 'CELL_TYPES ',amesh%Ncells
  do i = 1,amesh%Ncells
       write(dev,*) 5
  enddo

  return
end subroutine dumpmeshvtk_jup
!###############################################################################!###############################################################################
subroutine dumpmeshvtk(dev,amesh)

  type(mesh) :: amesh
  integer(pin) :: dev

  integer(pin) :: i,j

  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a8)') 'distance'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  if( pr == 4) then
     write(dev,'(a7,i10,a6)') 'POINTS ',amesh%Nnodes,' float'
  else
     write(dev,'(a7,i10,a7)') 'POINTS ',amesh%Nnodes,' double'
  endif

  do i=1,amesh%Nnodes
     write(dev,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  write(dev,'(a9,2i10)') 'POLYGONS ',amesh%Ncells,amesh%Ncells*4
  do i=1,amesh%Ncells
     write(dev,'(a2,$)') '3 '
     write(dev,*) (amesh%cell(i,j)-1,j=1,3)
  enddo
  return
end subroutine dumpmeshvtk
!###############################################################################
subroutine idumpnodeattributevtk(dev,amesh,field,attname,init)
  type(mesh) :: amesh
  integer(pin), dimension(amesh%Nnodes) :: field
  integer(pin) :: dev
  character*(*) :: attname
  logical :: init

  integer(pin) :: i

  ! if (init) write(dev,'(a11,i10)') 'POINT_DATA ',amesh%Nnodes
  ! if( pr == 4) then
    write(dev,'(a8,a,a8)') 'SCALARS ',trim(attname),' int 1'
  ! else
  !    write(dev,'(a8,a,a10)') 'SCALARS ',trim(attname),' double 1'
  ! endif
  write(dev,'(a20)') 'LOOKUP_TABLE default'
  do i=1,amesh%Nnodes
     write(dev,*) field(i)
  enddo
  return
end subroutine idumpnodeattributevtk
!###############################################################################
subroutine rdumpnodeattributevtk(dev,amesh,field,attname,init)
  type(mesh) :: amesh
  real(pr), dimension(amesh%Nnodes) :: field
  integer(pin) :: dev
  character*(*) :: attname
  logical :: init

  integer(pin) :: i

  if (init) write(dev,'(a11,i10)') 'POINT_DATA ',amesh%Nnodes
  if( pr == 4) then
    write(dev,'(a8,a,a8)') 'SCALARS ',trim(attname),' float 1'
  else
     write(dev,'(a8,a,a10)') 'SCALARS ',trim(attname),' double 1'
  endif
  write(dev,'(a20)') 'LOOKUP_TABLE default'
  do i=1,amesh%Nnodes
     write(dev,*) field(i)
  enddo
  return
end subroutine rdumpnodeattributevtk
!###############################################################################
subroutine dumpcellattributevtk(dev,amesh,field,attname,init)
  type(mesh) :: amesh
  real(pr), dimension(amesh%Ncells) :: field
  integer(pin) :: dev
  character*(*) :: attname
  logical :: init

  integer(pin) :: i

  if (init) write(dev,'(a10,i10)') 'CELL_DATA ',amesh%Ncells
  if( pr == 4) then
     write(dev,'(a8,a,a8)') 'SCALARS ',trim(attname),' float 1'
  else
     write(dev,'(a8,a,a10)') 'SCALARS ',trim(attname),' double 1'
  endif

  write(dev,'(a20)') 'LOOKUP_TABLE default'
  do i=1,amesh%Ncells
     write(dev,*) field(i)
  enddo
  return
end subroutine dumpcellattributevtk
!###############################################################################

!========================================================
subroutine find_1dR(array,condt,target_value,idx)
  implicit none
  real(pr),dimension(:),intent(in) :: array
  integer(pin),allocatable,dimension(:), intent(out) :: idx
  integer(pin),allocatable,dimension(:) :: loc
  character(2),intent(in) :: condt
  integer(pin) :: i,n,target_value,inx
  
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
  integer(pin),dimension(:),intent(in) :: array
  integer(pin),allocatable,dimension(:), intent(out) :: idx
  integer(pin),allocatable,dimension(:) :: loc
  integer(pin) :: i,n,target_value,inx
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
  real(pr),dimension(:,:),intent(in) :: array
  integer(pin),allocatable,dimension(:), intent(out) :: idx,jdx
  integer(pin),allocatable,dimension(:) :: loc_i,loc_j
  real(pr) :: target_value
  integer(pin) :: i,j,n,m,inx
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
  integer(pin),dimension(:,:),intent(in) :: array
  integer(pin),allocatable,dimension(:), intent(out) :: idx,jdx
  integer(pin),allocatable,dimension(:) :: loc_i,loc_j
  integer(pin) :: i,j,n,m,target_value,inx
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
     integer(pin) :: nv,nnv
     integer(pin), dimension(nv) :: val
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
  !###############################################################################
subroutine write_out(amesh,quake,surface,out_type,out_file,pdf,en_var,output_level,dist2)
  ! use typedef
  use LAT_source
  use LAT_time 
  implicit none
  type(mesh) :: amesh        ! complete mesh
  type(source) :: quake   ! mesh containing earthquake 
  type(reflect) :: surface
  
  character(3) :: out_type   !type of output file
  character(20) :: en_var
  character(30) :: out_file  ! name of output file
  character(34) :: filename  ! name of output file with extension
  character (len=20) :: fnamef
  integer :: output_level
  real(pr), dimension(quake%QuakeElemNo) :: pdf
  real(pr), dimension(:), allocatable :: slipout,pdfout,test,rupt_time
  real(pr), dimension(:), allocatable :: dist2
  integer :: i,j
  
  ! assigning a slip value to every element in whole mesh
  allocate(slipout(amesh%Ncells))
  allocate(pdfout(amesh%Ncells))
  allocate(rupt_time(amesh%nnodes))

  slipout = 0.
  pdfout = 0.
  do i = 1,quake%QuakeElemNo
    slipout(quake%QuakeElem(i)) = slip(i)
    pdfout(quake%QuakeElem(i)) = pdf(i)
  enddo
  
  if (en_var == ''.or.en_var=='1') then !   write(*,*) "No pbs_array
     select case(out_type)
      case('gmt')
      write(filename,'(a,a)') trim(adjustl(out_file)),'.gmt'
      open(10,file=trim(adjustl(filename)))
  
      do i=1,amesh%Ncells
        if (slipout(i) < 10) then  ! if slip is between 0-9
          write(10,'(a4,f9.7)') '> -Z',slipout(i)
         else      ! if slip is between 10 - 99
           write(10,'(a4,f10.7)') '> -Z',slipout(i)
         endif
         do j=1,3
            write(10,*) amesh%px(amesh%cell(i,j)),amesh%py(amesh%cell(i,j)),amesh%pz(amesh%cell(i,j))
         enddo
      enddo
      close(10)
  
      case('vtk')
           write(filename,'(a,a)') trim(adjustl(out_file)),'.vtk'
           open(11,file=trim(adjustl(filename)),form = 'formatted')
          !  call dumpmeshvtk(11,amesh) 
           call dumpmeshvtk_jup(11,amesh)  ! changed to this as it is read by meshio (other doesn't recognise POLYDATA)
! save rupture time on the nodes 
           rupt_time = 0._pr
           if (quake%QuakeNodesNo < amesh%nnodes) then 
            ! set rupture front to zero if 
               do i= 1,quake%QuakeNodesNo
                  rupt_time(quake%QuakeNodes(i)) = time(quake%QuakeNodes(i))
               enddo
            else
              rupt_time = time
           endif 
           call dumpnodeattributevtk(11,amesh,rupt_time,'Rupt_time',.true.)

          !  allocate(test(amesh%nnodes))
          !  test = 0._pr 
          !  do i = 1,surface%nnodes
          !   test(surface%nodes(i)) = 1._pr
          !  enddo
          !  call dumpnodeattributevtk(11,amesh,test,'surface_nodes',.true.)
          !  deallocate(test)
          !  if (allocated(dist2)) then 
          !   call dumpnodeattributevtk(11,amesh,dist2,'surface_dist',.false.)
          !  endif 


! save slip on cells 
           call dumpcellattributevtk(11,amesh,slipout,'slip',.true.)
           write (fnamef, "(a,I5.5)") "slip."//trim(en_var)

           call dumpcellattributevtk(11,amesh,velocity,'vel',.false.)
  
           if (output_level > 1) then
               slipout=0._pr
               do i = 1,quake%QuakeElemNo
                    slipout(quake%QuakeElem(i))= slip_prim(i)
               enddo
               call dumpcellattributevtk(11,amesh,slipout,'slip_prim',.false.)
               slipout=0._pr
               do i=1,quake%QuakeElemNo
                 slipout(quake%QuakeElem(i))= slip_sec(i)
               enddo
               call dumpcellattributevtk(11,amesh,slipout,'slip_sec',.false.)
            endif
            if (output_level > 2) then
               call dumpnodeattributevtk(11,amesh,Dist2Border,'Dist2Border',.true.)
               call dumpcellattributevtk(11,amesh,pdfout,'pdf',.true.)
            endif

            close(11)
      case default
              write(*,*) 'ERROR in creating output file for single event'
              stop
      end select
  
  
    else  !  write(*,*) "INDEX NUMBER  ",en_var
      write (fnamef, "(a,I5.5)") "slip."//trim(en_var)
      open (10,file=fnamef, form="formatted",status="unknown")
      do i = 1,quake%QuakeElemNo
          write(10,*) slipout(i)
      enddo
      close(10)
  
      if (output_level > 1) then
           slipout=0.
           do i=1,quake%QuakeElemNo
              slipout(quake%QuakeElem(i))= slip_prim(i)
           enddo
           write (fnamef, "(a,I5.5)") "slip_prim."//trim(en_var)
           open (10,file=fnamef, form="formatted",status="unknown")
           do i = 1,amesh%Ncells
              write(10,*) slipout(i)
           enddo
           close(10)
           slipout=0.
           do i=1,quake%QuakeElemNo
                 slipout(quake%QuakeElem(i))= slip_sec(i)
           enddo
           open (10,file=fnamef, form="formatted",status="unknown")
           do i = 1,amesh%Ncells
              write(10,*) slipout(i)
           enddo
           close(10)
      endif
  
      if (output_level > 2) then
           write (fnamef, "(a,I5.5)") "pdf."//trim(en_var)
           open (10,file=fnamef, form="formatted",status="unknown")
           do i = 1,amesh%Ncells
               write(10,*) pdfout(i)
           enddo
           close(10)
  
      endif
  endif
  ! !--------- ouput boundary nodes of slipping region
  !            open(12,file='border.vtk',form = 'formatted')
  !            call dumpQuakeBordervtk(12,amesh)
  !            close(12)
  !          !--------- output all the nodes that make up slipping region
  !            open(12,file='quakenodes.vtk',form = 'formatted')
  !            call dumpQuakevtk(12,amesh)
  !            close(12)
  
  return
  end subroutine write_out
  !###############################################################################
  subroutine write_model_parameters(model,en_var)
  use LAT_source
    ! use lateration
    implicit none
    character(20) :: en_var  ! name of output file
    character (len=20) :: fnamef
  
    type(model_param),intent(in) :: model
  
    if (en_var == '') then
       write (fnamef, "(a,I5.5)") "model_parmaters"
    else
       write (fnamef, "(a,I5.5)") "model_parmaters."//trim(en_var)
    endif
    open(11,file=trim(adjustl(fnamef)),form = 'formatted')
  !  open(11,file='model_parmaters.txt',form = 'formatted')
    write(11,*) 'Target Area....',model%target_area
    write(11,*) 'Actual Area....',model%actual_area
    write(11,*) 'Target Width....',model%target_width
    write(11,*) 'Actual Width....',model%actual_width
    write(11,*) 'Length..........',model%actual_length
  
    close(11)
  
  end subroutine write_model_parameters
  !###############################################################################

!###############################################################################
end module LAT_mesh_util
