module typedef
use lateration
use utils
implicit none

contains
!==============================================================================
!=================================================
subroutine read_input(model, mesh_file, pdf,out_type,out_file,surf,output_level)
  type(pdfinputs) :: pdf
  type(model_param) :: model
  type(surface_type) :: surf
real ::  mu, rmin, rmax, param
integer :: na
integer  :: io_n
character(30),intent(out) :: mesh_file,out_file
character(30) :: fname,surf_file
character(3),intent(out) :: out_type
character(7) :: pdf_type
character(3) :: param_type
integer :: gauss_no,output_level
logical :: defined_area,use_file,surf_exists


NAMELIST / GENERAL  / mu, rmin, rmax, na, mesh_file
NAMELIST / SIZE / param_type,param, defined_area
NAMELIST / PDF1   / pdf_type
NAMELIST / PDF2  / gauss_no, fname
NAMELIST / SURFACE / surf_exists, use_file ,surf_file
NAMELIST / OUTPUT  / out_type, out_file, output_level

io_n = 11
open(io_n,file='input_file',status='old',action='read')

! Values to be read in:
rewind(io_n)
read(io_n,nml=GENERAL,END=200)

rewind(io_n)
read(io_n,nml=SIZE,END=200)

rewind(io_n)
read(io_n,nml=PDF1,END=200)

rewind(io_n)
read(io_n,nml=PDF2,END=200)

rewind(io_n)
read(io_n,nml=SURFACE,END=200)

surf%exists = surf_exists
if (surf%exists) then
  surf%use_file = use_file
  if(surf%use_file) then
    surf%filename = surf_file
  else
    surf%z0 =  0.0
  endif
endif

rewind(io_n)
read(io_n,nml=OUTPUT,END=200)

model%mu = mu
model%rmin = rmin

model%rmax = rmax
model%na = na
model%defined_area = defined_area
model%param_type = param_type


pdf%pdf_type = pdf_type
select case(pdf%pdf_type)
  case('gauss')
    pdf%ng = gauss_no
  case('defined')    ! pdf has been defined in an ascii file
    pdf%pdf_fname = fname
! check that user is using the full the mesh in this case
    if(model%defined_area) then
      write(*,*) 'ERROR: if using a predefined PDF the whole mesh is used'
      write(*,*) '       in input_file either set:                       '
      write(*,*) '                     (a) defined_area=.false.          '
      write(*,*) '                     (b) use a different pdf_type      '
     stop
   endif

end select

select case(model%param_type)
 case('sd')
   model%sd = param
   if ( model%defined_area ) then
        write(0,*) 'ERROR: defining slip based on stress drop requires a predefined fault area'
        stop
   endif
 case('mw')
     model%mw = param
case default
   write(0,*) 'ERROR: param_type in input file not defined correctly'
   stop
end select


return
! error message on reading input file
200 write(*,*) 'ERROR reading input_file'

stop

end subroutine read_input
!==============================================================================
subroutine read_mesh_file(amesh,mesh_file,surf)
type(mesh) :: amesh
type(surface_type) :: surf
character(30), intent(in) :: mesh_file
!character(30), intent(in) :: surf_file
! read in quadrangle mesh from file
!call read_quad_mesh(amesh,mesh_file)

! read unstructured mesh from file
call read_tri_mesh(amesh,mesh_file)

! identify edge connection between elements
call mesh_EToE(amesh)

! read in surface nodes
if (surf%exists) then
   call surface_nodes(amesh,surf)
endif

! calc mean depth of each element
call calc_mean_depth(amesh)

! area is in S.I.
call calc_area(amesh)

! calculate distances
call allvsall2d(amesh,amesh%dist)



end subroutine read_mesh_file
!==============================================================================
subroutine read_pdf_file(pdf,amesh)
  type(pdfinputs) :: pdf
! main structure
  type(mesh) :: amesh
  integer :: i

allocate(pdf%g_pdf(amesh%Ncells))
! the pdf is defined for each cell
! the order of the cells is assumed to be the same as the mesh file
open(13,file=trim(adjustl(pdf%pdf_fname)),form = 'formatted', action = 'read')
do i =1,amesh%Ncells
    read(13,*) pdf%g_pdf(i)
enddo
close(13)

return
end subroutine read_pdf_file
!==============================================================================
subroutine surface_nodes(amesh,surf)
  type(mesh) :: amesh
  type(surface_type) :: surf
  integer,allocatable,dimension(:) :: surface_mask
  integer :: nlines,i,j
  logical :: not_on_boundary


   if (surf%use_file) then
       open(13,file=trim(adjustl(surf%filename)),form = 'formatted', action = 'read')
!! find out how many data points/lines in the file
       nlines = 0
       do
         read (13,*, END=10)
        nlines = nlines + 1
      enddo
      10 rewind(13)

      amesh%NoSurfNodes = nlines
      allocate(amesh%SurfNodes(amesh%NoSurfNodes))
      do i = 1, amesh%NoSurfNodes
        read(13,*) amesh%SurfNodes(i)
      enddo
      close (13)
   else    ! use a predefined value (i.e. 0)
      allocate(surface_mask(amesh%Nnodes))
      surface_mask = 0
      nlines = 0
      do i = 1,amesh%Nnodes
        if(amesh%pz(i) >= surf%z0) then
!        if(amesh%pz(i) <= surf%z0+epsilon(abs(amesh%pz(i)))) then
         nlines = nlines + 1
         surface_mask(nlines) = i
        endif
      enddo
      amesh%NoSurfNodes = nlines
      allocate(amesh%SurfNodes(amesh%NoSurfNodes))
      amesh%SurfNodes = surface_mask(1:nlines)
      deallocate(surface_mask)
    endif

    if (nlines == 0 ) then
       write(*,*) 'Warning: No surface nodes found'
       write(*,*) '         => will assume buried fault'
    else
    ! check that surface is along boundary of mesh //SHANE
    do i = 1,amesh%NoSurfNodes
      not_on_boundary = .true.
      do j = 1,amesh%Ncells
        if (any(amesh%SurfNodes(i).eq.amesh%cell(j,:))) then
          if (any(j.eq.amesh%EToE(j,:))) then
! element is a boundary element
            not_on_boundary = .false.
          endif
        endif
      enddo
      if (not_on_boundary) then
         write(*,*) 'ERROR: Defined surface node not on edge of mesh'
         write(*,*) 'ERROR: Node not on mesh edge:', amesh%SurfNodes(i)
         stop
      endif
    enddo
  endif
!    do i = 1,amesh%Ncells
!        do j = 1,3
!          ! check if element contains surface node
!            if (any(amesh%cell(i,j).eq.amesh%SurfNodes) ) then
!            ! if it does, check that the element is on boundary of mesh
!              if (any(i.eq.amesh%EToE(i,:))) then
!! element is a boundary element
!              else
!! element is an interior element
!                  write(*,*) 'ERROR: Surface node not on edge of mesh'
!                  print*,'cell:   ',i
!                  print*,'EToE    ',amesh%EToE(i,:)
!                  stop
!                ! if any of k index values in EToE = i value => there are less than 3 faces
!                ! 3 faces = internal element
!                ! 2 faces = boundary element
!                ! 1 face = edge element
!              endif
!            endif
!        enddo
!    enddo



return
end subroutine
!==============================================================================
subroutine read_quad_mesh(amesh,mesh_file)
!!  Reads in a quadrangle mesh and converts it to a triangular
!!
type(mesh) :: amesh
character(30) :: mesh_file
integer :: i,id,quad_cells
character(len=60) :: str
integer :: inx
real :: num
integer :: id_o,a,b,c,d

open(13,file=trim(adjustl(mesh_file)),form = 'formatted', action = 'read')
! skip header
do i = 1,9
       read(13,*)
enddo
! find end of nodal section & element section
i = 0
id = 0
do while ( id == 0 )
    read(13,*) str
    if (trim(str) == '**') then
        id = 1
    endif
    i = i+1
enddo
amesh%Nnodes = i-1
allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))

do i = 1,2
       read(13,*) str
enddo
id = 0
i = 0
do while ( id == 0 )
    read(13,*) str
    if (trim(str) == '**') then
        id = 1
    endif
    i = i+1
enddo
quad_cells = i-1
amesh%Ncells = quad_cells*2
allocate(amesh%cell(amesh%Ncells,3))

!====== now read data in
rewind(13)
! skip header
do i = 1,9
       read(13,*)
enddo

do i = 1,amesh%Nnodes
       read(13,*) id,amesh%px(i),amesh%py(i),amesh%pz(i)
       if (i == 1) id_o = id
enddo

! skip header
do i = 1,3
       read(13,*) str
enddo
! read element to vertices
inx = 1
do i = 1,quad_cells
       read(13,*) id,a,b,c,d
       a = a-id_o+1
       b = b-id_o+1
       c = c-id_o+1
       d = d-id_o+1
       call random_number(num)     ! randomly choose split of cell
         if (num <0.5) then
             amesh%cell(inx,:) = (/ a, b, d/)
             amesh%cell(inx+1,:) = (/ b, c, d/)
         else
             amesh%cell(inx,:) = (/ a, b, c/)
             amesh%cell(inx+1,:) = (/ c, d, a/)
         endif
         inx = inx+2
enddo
close(13)

return
end subroutine read_quad_mesh
!==============================================================================
subroutine read_tri_mesh(amesh,mesh_file)
!  Reads in an unstructured mesh that has .inp format
!
!
type(mesh) :: amesh
integer :: i,id
character(len=60) :: str
real :: x,y,z
character(30),intent(in) :: mesh_file
integer :: a,b,c
integer,allocatable,dimension(:) :: counter,idx


open(13,file=trim(adjustl(mesh_file)),form = 'formatted', action = 'read')

! skip header
do i = 1,9
       read(13,*)
enddo
! find end of nodal section & element section
i = 0
id = 0
do while ( id == 0 )
    read(13,*) str
    if (trim(str) == '**') then
        id = 1
    endif
    i = i+1
enddo
amesh%Nnodes = i-1
!print*,'No. of Nodes ....',amesh%Nnodes
allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))
amesh%px(amesh%Nnodes) = 0.
amesh%py(amesh%Nnodes) = 0.
amesh%pz(amesh%Nnodes) = 0.


do i = 1,2
       read(13,*) str
enddo
id = 0
i = 0
do while ( id == 0 )
    read(13,*) str
    if (trim(str) == '**') then
        id = 1
    endif
    i = i+1
enddo
amesh%Ncells = i-1
!print*,'No of cells......',amesh%Ncells
allocate(amesh%cell(amesh%Ncells,3))


!====== now read data in =======
rewind(13)
! skip header
do i = 1,9
       read(13,*)
enddo
! read in vertices data
allocate(counter(amesh%Nnodes))
do i = 1,amesh%Nnodes
  read(13,*) id,x,y,z
  counter(i) = id
  amesh%px(i) = x
  amesh%py(i) = y
  amesh%pz(i) = z
enddo

! skip header
do i = 1,3
       read(13,*) str
enddo
! read element to vertices
do i = 1,amesh%Ncells

       read(13,*) id,a,b,c
       if (allocated(idx)) deallocate(idx)
       call  find(counter,'==',a,idx)
       if (idx(1).eq.-1) then
          print*,"problem reading in elements"
          stop
       else
         amesh%cell(i,1) = idx(1)
       endif
       if (allocated(idx)) deallocate(idx)
       call  find(counter,'==',b,idx)

       if (idx(1).eq.-1) then
          print*,"problem reading in elements!"
          stop
       else
         amesh%cell(i,2) = idx(1)
       endif

       if (allocated(idx)) deallocate(idx)

       call  find(counter,'==',c,idx)

       if (idx(1).eq.-1) then
          print*,"problem reading in elements!"
          stop
       else
         amesh%cell(i,3) = idx(1)
       endif
enddo
close(13)
end subroutine
!==============================================================================
subroutine calc_mean_depth(amesh)
implicit none
type(mesh) :: amesh
real,dimension(3) :: z
integer :: i

allocate(amesh%mz(amesh%Ncells))
do i = 1,amesh%Ncells
  z = amesh%pz(amesh%cell(i,:))
  amesh%mz(i) = sum(z)/3.
enddo

end subroutine calc_mean_depth
!==============================================================================
subroutine calc_area(amesh)
implicit none
type(mesh) :: amesh
integer :: i
real,dimension(3) :: x,y,z,vec_ab,vec_ac
real :: ab,ac,res,theta
real :: tot_area
allocate(amesh%area(amesh%Ncells))

tot_area = 0
do i = 1,amesh%Ncells
  x = amesh%px(amesh%cell(i,:))
  y = amesh%py(amesh%cell(i,:))
  z = amesh%pz(amesh%cell(i,:))
       vec_ab = (/ x(2)-x(1), y(2)-y(1), z(2)-z(1) /)
       vec_ac = (/ x(3)-x(1), y(3)-y(1), z(3)-z(1) /)
  ac = norm2(vec_ac)
  ab = norm2(vec_ab)
  res = dot_product(vec_ab, vec_ac)/ab/ac
  theta = acos(res)
  amesh%area(i) = ab*ac*sin(theta)/float(2)
  tot_area = tot_area+amesh%area(i)
enddo
!print*,'Total area (m^2):......', tot_area
end subroutine calc_area
!==============================================================================
!==============================================================================
subroutine mesh_EToE(amesh)
!  Create a matrix that defines neighbouring cells
!  Dimension :  (No of cells , 3)
!  Row relates to cells
!  Columns give the cell numbers of the 3 cells with adjoining faces to current cell
!  If a cell has less than 3 neighbours the index of the current cell is placed in the column
!
!
type(mesh),intent (inout) :: amesh
integer, dimension(3) :: node_id,neighbour_element
integer,allocatable,dimension(:) :: common_node
integer :: i,k,nofaces,bc_int


allocate(amesh%EToE(amesh%Ncells,3))
do i = 1,amesh%Ncells
       amesh%EToE(i,:) = i
enddo

bc_int = 0
! find boundary of mesh
do i = 1,amesh%Ncells
     node_id = amesh%cell(i,:)
     nofaces = 0
     neighbour_element = 0
! check node id against rest of array
     do k = 1,amesh%Ncells
       if (i /= k) then
           if (any(node_id(1).eq.amesh%cell(k,:)).or.      &
               any(node_id(2).eq.amesh%cell(k,:)).or.      &
               any(node_id(3).eq.amesh%cell(k,:))) then
! compare all the nodes of the two elements
              call intersect(node_id,amesh%cell(k,:),common_node)
              if (size(common_node)==2) then
          ! two common nodes = face
                    nofaces = nofaces+1
                    neighbour_element(nofaces) = k
              endif
            endif
       endif

     enddo
     if (nofaces > 0) then
       do k = 1,nofaces
          amesh%EToE(i,k) = neighbour_element(k)
       enddo
     endif
! if any of k index values in EToE = i value => there are less than 3 faces
! 3 faces = internal element
! 2 faces = boundary element
! 1 face = edge element

enddo

end subroutine mesh_EToE
!==============================================================================
subroutine  select_mesh_bc(amesh,model)
!  Select the whole mesh as the rupture area
!
type(model_param),intent (inout) :: model
type(mesh),intent (inout) :: amesh
integer :: i

! assume whole mesh is included in earthquake surface
allocate(amesh%QuakeElem(amesh%Ncells),amesh%QuakeNodes(amesh%Nnodes))

do i = 1,amesh%Nnodes
  amesh%QuakeNodes(i) = i
enddo
do i = 1,amesh%Ncells
  amesh%QuakeElem(i) = i
enddo
amesh%QuakeElemNo = amesh%Ncells
amesh%QuakeNodesNo = amesh%Nnodes

model%actual_width = maxval(abs(amesh%pz(:)))
model%actual_area = maxval(abs(amesh%px(:)))*model%actual_width
model%target_width = model%actual_width

return
end subroutine select_mesh_bc
!==============================================================================
subroutine  select_fault_zone_global(amesh,dist,model,z0,pdf)
  type(mesh),intent (inout) :: amesh
  type(model_param),intent (inout) :: model
  type(pdfinputs),intent (in) :: pdf

  integer,allocatable,dimension(:) :: quake_nodes,quake_cells,haveit,potential_cells
  !real,allocatable,dimension(:) :: pdf
  integer :: cell_options(3),c_nodes(3),test_cell(3)
  real :: vec(3), curcum
  real :: random , quake_area,prob_lr, prob_ud
  integer :: get_out, inodes ,icells, c_cell ,i,iu_cells,random_pick
  integer :: lr_count, ud_count, tot_count, lr_elements(3), ud_elements(3)
	integer :: counter,ipcells,ii,k,i_use
	real :: epi_depth,d, ug,sigma,pdfint,min_z,max_z,width_est
  real, dimension(:,:), allocatable :: dist
  real :: z0   ! depth of free surface
  logical :: find_cell,no_holes,reached_width
  real :: minv,maxv,set_width,trial_width,dmax
  real :: numrnd
  real,allocatable,dimension(:) ::dz,dd
  real, allocatable,dimension(:) :: weighting
  integer :: opposite_idx
  integer :: min_loc,max_loc
  integer :: nlength

  get_out = 0
  inodes = 1
  icells = 1
  allocate(quake_nodes(amesh%Nnodes),quake_cells(amesh%Ncells),potential_cells(amesh%Ncells))
  allocate(weighting(amesh%Ncells))
  weighting = 0.

  quake_nodes = 0
  quake_cells = 0
  potential_cells = 0
  quake_area = 0.
  ipcells = 1
  counter = 0
  reached_width = .false.
  do while (get_out < 1)
      counter= counter+1
      if (inodes == 1) then
! randomly pick starting cell
            call random_number(random)
            c_cell = ceiling(random*float(amesh%Ncells))
        quake_nodes(inodes:inodes+2) = amesh%cell(c_cell,:);
        inodes = inodes+3
        epi_depth = amesh%mz(c_cell)
      else
!======================pick next cell ==================================
! check to see if any of the candidate cells are surrounded by a slipping cells
        no_holes = .true.
        do ii = 1,ipcells-1
            test_cell = amesh%EToE(potential_cells(ii),:)
            if (any(test_cell(1).eq.quake_cells(1:icells-1)).and.      &
                any(test_cell(2).eq.quake_cells(1:icells-1)).and.      &
                any(test_cell(3).eq.quake_cells(1:icells-1))) then
! we have a hole!
               no_holes = .false.
			         c_cell = potential_cells(ii)
               i_use = ii
            endif
        enddo
       if (no_holes) then
          do ii = 1,ipcells-1
              test_cell = amesh%EToE(potential_cells(ii),:)
              if (any(test_cell(1).eq.quake_cells(1:icells-1)).and.      &
                  any(test_cell(2).eq.quake_cells(1:icells-1)).or.       &
                  any(test_cell(2).eq.quake_cells(1:icells-1)).and.      &
                  any(test_cell(3).eq.quake_cells(1:icells-1)).or.       &
                  any(test_cell(3).eq.quake_cells(1:icells-1)).and.      &
                  any(test_cell(1).eq.quake_cells(1:icells-1))) then
! this gives twice the weighting to potential cells who have two faces adjoining
! the slipping region twice the likelihood of being picked compared to ones with only one
                  weighting(ii) = 2.
                else
                  weighting(ii) = 1.
              endif
           enddo
!======= randomly choose next cell from list of potental
          call random_number(numrnd)
          ii = 0
          curcum = 0.
          do while (curcum <= numrnd)
               ii=ii+1
               curcum=curcum+weighting(ii)
          enddo
!           call random_number(random)
!            ii = ceiling((random)*float(ipcells-1))
            c_cell = potential_cells(ii)
            i_use = ii
            weighting(ii) = 0.
        endif
     endif
!====================  add cell to slipping list & remove from potential list
     quake_cells(icells) = c_cell
     quake_area = quake_area + amesh%area(quake_cells(icells)) ! assume area is in m^2
     icells = icells+1

! remove used cell from list
    if (counter > 1) then
          do  i = i_use,ipcells
            potential_cells(i) = potential_cells(i+1)
         enddo
         ipcells = ipcells-1
     endif

!===================== check if nodes are already defined as slipping, if not add them
     if (counter > 1) then
		 do i = 1,3
			 call find(quake_nodes(1:inodes-1),'==',amesh%cell(c_cell,i),haveit)
			 if (haveit(1) == -1) then
			 	 quake_nodes(inodes) = amesh%cell(c_cell,i);
				 inodes = inodes+1
			 endif
		  enddo
	  endif
!=====================if ajoining elements are not in potential or quake lists add them to potential
     cell_options = amesh%EToE(c_cell,:)   !find ajoining elements
			do i = 1,3
				    call find(potential_cells(1:ipcells-1),'==',cell_options(i),haveit)
				    if (haveit(1) == -1) then
				       call find(quake_cells(1:icells-1),'==',cell_options(i),haveit)
					     if (haveit(1) == -1) then
! check to see if width constraint is required
!                    if (reached_width) then
!                      if (  (abs(amesh%mz(cell_options(i))) > min_z).and.     &
!                            (abs(amesh%mz(cell_options(i))) < max_z)) then
!                              potential_cells(ipcells) = cell_options(i);
!         		                  ipcells = ipcells+1
!                      endif
!                   else
		                 potential_cells(ipcells) = cell_options(i);
    		             ipcells = ipcells+1
!                   endif   ! checking width
!	               endif   ! checking not above surface
               endif ! haveit check 2
            endif ! haveit check 1
       enddo
!======================= check to see if there are any more elements that can be used
      if(ipcells == 1) then
  	      write(*,*) '***WARNING: run out of cells ***'
	        write(*,*) '***increase mesh ***'
          write(*,*) '*** will keep running but assuming area is smaller than Strasser scaling'
          write(*,*) '*** (this implies a larger stress drop)'
  !                      pause
        	get_out = 1
	        cycle
      endif

! ===================== check to see if fault size has been reached==========
	if (quake_area >= model%target_area) then
		get_out = 1
	elseif (icells > amesh%Ncells) then
		write(*,*) '***WARNING: requested Mag. poss. too big for mesh ***'
!              pause
 		get_out = 1
		cycle
	endif
 enddo

!======================= check for any final holes =========================
no_holes = .true.
do ii = 1,ipcells-1
    test_cell = amesh%EToE(potential_cells(ii),:)
    if (any(test_cell(1).eq.quake_cells(1:icells-1)).and.      &
        any(test_cell(2).eq.quake_cells(1:icells-1)).and.      &
        any(test_cell(3).eq.quake_cells(1:icells-1))) then
! we have a hole!
        no_holes = .false.
        c_cell = potential_cells(ii)
        i_use = ii
        quake_cells(icells) = c_cell
        quake_area = quake_area + amesh%area(quake_cells(icells)) ! assume area is in m^2

        icells = icells+1
    endif
enddo

   model%actual_area = quake_area
 ! model%actual_width = width_est

  ! get rid of duplicates
   call erasedupli(quake_cells,icells-1,nlength)
   allocate(amesh%QuakeElem(nlength))
   amesh%QuakeElem = quake_cells(1:nlength)
   amesh%QuakeElemNo = nlength


   call erasedupli(quake_nodes,inodes-1,nlength)
   allocate(amesh%QuakeNodes(nlength))
   amesh%QuakeNodes = quake_nodes(1:nlength)
   amesh%QuakeNodesNo = nlength


  ! Save section of fault used for earthquake
!  allocate(amesh%QuakeElem(icells-1),amesh%QuakeNodes(inodes-1))
!  amesh%QuakeElem = quake_cells(1:icells-1)
!  amesh%QuakeNodes = quake_nodes(1:inodes-1)
!  amesh%QuakeElemNo = icells-1
!  amesh%QuakeNodesNo = inodes-1
  deallocate(quake_nodes,quake_cells)

  !write(*,*) 'number of elements in slipping zone',amesh%QuakeElemNo
  !write(*,*) 'number of nodes in slipping zone',amesh%QuakeNodesNo


end subroutine select_fault_zone_global
!==============================================================================!==============================================================================
subroutine  select_fault_zone_element(amesh,target_area)
  type(mesh),intent (inout) :: amesh
  integer,allocatable,dimension(:) :: quake_nodes,quake_cells,haveit,unused_cells
  real, intent(in) :: target_area
  integer :: cell_options(3),c_nodes(3)
  real :: vec(3)
  real :: random, alpha , quake_area,prob_lr, prob_ud
  integer :: get_out, inodes ,icells, c_cell ,i,iu_cells,random_pick
  integer :: lr_count, ud_count, tot_count, lr_elements(3), ud_elements(3)
	integer :: counter
  logical :: find_cell
  alpha = 2.
  get_out = 0
  inodes = 1
  icells = 1
	iu_cells = 0
  allocate(quake_nodes(amesh%Nnodes),quake_cells(amesh%Ncells))
  allocate(unused_cells(amesh%Ncells))
  quake_nodes = 0
  quake_cells = 0
  unused_cells = 0

  do while (get_out < 1)
         if (inodes == 1) then
  ! randomly pick starting cell
          call random_number(random)
          c_cell = ceiling(random*float(amesh%Ncells))   !
          quake_nodes(inodes:inodes+2) = amesh%cell(c_cell,:);
          inodes = inodes+3
          quake_cells(icells) = c_cell
          quake_area = quake_area + amesh%area(quake_cells(icells)) ! assume area is in m^2
          icells = icells+1

        else
          lr_count = 0;ud_count=0;
          lr_elements = 0; ud_elements = 0;
          cell_options(:) = amesh%EToE(c_cell,:)   !find ajoining elements
          do i = 1,3
! check to see if elements are already on fault
            call find(quake_cells(1:icells-1),'==',cell_options(i),haveit)

            if (haveit(1) == -1) then      ! don't have element
! work out if element is left-right or up - down of current element
                c_nodes = amesh%cell(c_cell,:)
                vec = amesh%pz(c_nodes)
                if (amesh%mz(cell_options(i)) >= minval(vec).and.     &
                    amesh%mz(cell_options(i)) <= maxval(vec))   then  ! left/right
                   lr_count = lr_count+1
                   lr_elements(lr_count) = cell_options(i)
!                   print*,'l/r '
                else    ! up/down
                   ud_count = ud_count+1
                   ud_elements(ud_count) = cell_options(i)
                endif
            endif
          enddo ! loop over candidate cells

! choose the next element
          tot_count = lr_count+ud_count
  !        print*,tot_count,lr_count,ud_count
          if (tot_count == 1) then ! if there is one candidate cell
            if (lr_count == 1) then
!              print*, '1 candidate cell : l/r'
              c_cell = lr_elements(1)
            else
!              print*, '1 candidate cell : u/d'
              c_cell = ud_elements(1)
            endif
          elseif (tot_count == 2) then ! if there are 2 candidate cells
              call random_number(random)
              if (lr_count == 2) then
!                print*, '2 candidate cells : both l/r'

                  prob_lr = 0.5
                  if (random <= prob_lr) then
                    c_cell = lr_elements(1)
                  else
                    c_cell = lr_elements(2)
                  endif
              elseif (ud_count == 2)  then
!                print*, '2 candidate cells : both u/d'

                  prob_ud = 0.5
                  if (random <= prob_ud) then
                    c_cell = ud_elements(1)
                  else
                    c_cell = ud_elements(2)
                  endif
              else
!                print*, '2 candidate cells : 1 u/d; 1 l/r'

                prob_lr = alpha/(alpha+1.)
                if (random <= prob_lr) then
                  c_cell = lr_elements(1)
                else
                  c_cell = ud_elements(1)
                endif
              endif
          elseif (tot_count == 3) then  ! if there are 3 candidate cells
              if  (lr_count == 2) then ! if there are 2 left/right elements, 1 up/down element
!                print*, '3 candidate cells : 1 u/d; 2 l/r'

                prob_lr = alpha/(2.*alpha+1.)
                call random_number(random)
                if (random <= prob_lr) then ! choose first l/r element
                  c_cell = lr_elements(1)
                elseif (prob_lr <=  random.and.random <= 2*prob_lr) then ! choose 2nd l/r element
                  c_cell = lr_elements(2)
                else ! choose u/d element
                  c_cell = ud_elements(1)
                endif
              else  ! there are 2 up/down elements, 1 left/right element
!                print*, '3 candidate cells : 2 u/d; 1 l/r'

                prob_lr = alpha/(alpha+2.)
                prob_ud = 1./(alpha+2.)
                call random_number(random)
                if (random <= prob_ud) then ! choose first l/r element
                  c_cell = ud_elements(1)
                elseif (prob_ud <=  random.and.random <= 2*prob_ud) then ! choose 2nd l/r element
                  c_cell = ud_elements(2)
                else ! choose u/d element
                  c_cell = lr_elements(1)
                endif
              endif
          else
!              print*,'problem : stuck in a hole !'
		find_cell =.true.
		counter = 0
! randomly choose one of the unused cells
      do while (find_cell)
         call random_number(random)
         random_pick = floor(float(iu_cells)*random)+1
         c_cell = unused_cells(random_pick)
! remove used cell from list
         do  i = random_pick,iu_cells
!	print*, i,unused_cells(i),unused_cells(i+1)
            unused_cells(i) = unused_cells(i+1)
         enddo
         iu_cells = iu_cells-1

! check whether used cell has not already been used
          call find(quake_cells(1:icells),'==',c_cell,haveit)
          if (haveit(1) == -1) then
                 find_cell =.false.
          endif
          counter = counter + 1
        enddo
      endif


    do i = 1,3
        call find(quake_nodes(1:inodes-1),'==',amesh%cell(c_cell,i),haveit)
        if (haveit(1) == -1) then
              quake_nodes(inodes) = amesh%cell(c_cell,i);
              inodes = inodes+1
        endif
    enddo
    quake_cells(icells) = c_cell
    quake_area = quake_area + amesh%area(quake_cells(icells)) ! assume area is in m^2
     icells = icells+1

     if (quake_area >= target_area) then
                get_out = 1
      elseif (icells > amesh%Ncells) then
                write(*,*) '***WARNING: requested Mag. poss. too big ***'
                get_out = 1
                cycle
      endif

! keep cells that haven't been used
         do i =1,3
             call find(quake_cells(1:icells-1),'==',cell_options(i),haveit)
              if (haveit(1) == -1) then
                   iu_cells = iu_cells+1
                   unused_cells(iu_cells) = cell_options(i)
              endif
         enddo
        endif  ! end of conditional on whether first element or not

  end do

  ! Save section of fault used for earthquake
  allocate(amesh%QuakeElem(icells-1),amesh%QuakeNodes(inodes-1))
  amesh%QuakeElem = quake_cells(1:icells-1)
  amesh%QuakeNodes = quake_nodes(1:inodes-1)
  amesh%QuakeElemNo = icells-1
  amesh%QuakeNodesNo = inodes-1

  write(*,*) 'number of elements in slipping zone',amesh%QuakeElemNo
  write(*,*) 'number of nodes in slipping zone',amesh%QuakeNodesNo


  deallocate(quake_nodes,quake_cells)

end subroutine select_fault_zone_element
!==============================================================================
subroutine  select_fault_zone(amesh,target_area)
!  Subroutine to pick rupture area
!  Cells are added iteratively until area exceeds target_area value
!  Choice of starting point is randomly choosen but it can be set by
!  prescribing a value to nuc_id
!
type(mesh),intent (inout) :: amesh
integer,allocatable,dimension(:) :: quake_nodes,quake_ameshents,irow,icol,haveit
real, intent(in) :: target_area
integer :: choice_amesh(3),current_nodes(3)
real :: vec(3)
real :: nuc_x,nuc_y,nuc_z
real :: quake_area
real :: random, alpha
integer :: inx,get_out,nuc_id,node_idx,nxt_pt,ele_inx
integer :: kk,i,k
integer :: current_element
integer :: lr_count,ud_count,ud_elements(3),lr_elements(3),tot_count


allocate(quake_nodes(amesh%Nnodes),quake_ameshents(amesh%Ncells))

alpha = 2. ! length to width ratio
quake_nodes = 0
quake_ameshents = 0
quake_area = 0.   ! initial area for earthquake zone on fault


inx = 1
node_idx = 1
ele_inx = 1
get_out = 0
do while (get_out < 1)
       if (inx == 1) then
! randomly pick starting cell
              call random_number(random)
              nuc_id = ceiling(random*amesh%Nnodes)   !
!            nuc_id = 1071  !centre of Scotia fault
              nuc_x = amesh%px(nuc_id)
              nuc_y = amesh%py(nuc_id);
              nuc_z = amesh%pz(nuc_id);
              nxt_pt = nuc_id
       else
! take next nodal point
              node_idx = node_idx+1
              if (node_idx > size(quake_nodes)) then
                  write(*,*) '***WARNING: requested Mag. poss. too big ***'
                  get_out = 1
                  cycle
              endif
              nxt_pt = quake_nodes(node_idx)
       endif
! find elements connected to node
       call find(amesh%cell,'==',nxt_pt,irow,icol)       ! find the index at where amesh%dist = 0
       do kk = 1,size(irow)   !loop over all possible nodes
         lr_count = 0
         ud_count = 0
         ud_elements = 0
         lr_elements = 0

              choice_amesh(:) = amesh%EToE(irow(kk),:)   !take element with a joining face
              do i = 1,3  ! loop over cells
                     if (inx == 1) then
                            quake_nodes(inx) = amesh%cell(choice_amesh(i),1);
                            quake_nodes(inx+1) = amesh%cell(choice_amesh(i),2);
                            quake_nodes(inx+2) = amesh%cell(choice_amesh(i),3);
                            inx = inx+3;
                     else               ! check to see if node is already in list
                         do k = 1,3
                            call find(quake_nodes,'==',amesh%cell(choice_amesh(i),k),haveit)
                            if (haveit(1) == -1) then      ! don't have node
                                   quake_nodes(inx) =  amesh%cell(choice_amesh(i),k);
                                   inx = inx+1;
                            endif
                         enddo
                     endif
! check to see if element is already in quake_ameshents
                    call find(quake_ameshents,'==',choice_amesh(i),haveit)
                    if (haveit(1) == -1) then      ! don't have element
                           if (inx == 1) then   ! initially include all surrounding elements
                             quake_ameshents(ele_inx) = choice_amesh(i);
                             current_element = choice_amesh(i);
                             quake_area = quake_area + amesh%area(quake_ameshents(ele_inx)) ! assume area is in m^2
                             ele_inx = ele_inx+1;
                           else
! ascertain if element is above/below or left/right of current element : shane
                                 current_nodes = amesh%cell(current_element,:)
                                 vec = amesh%pz(current_nodes)
                               if (amesh%mz(choice_amesh(i)) >= minval(vec).and.     &
                                   amesh%mz(choice_amesh(i)) <= maxval(vec))   then  ! left/right
                                  lr_count = lr_count+1
                                  lr_elements(lr_count) = choice_amesh(i)
                               else    ! up/down
                                  ud_count = ud_count+1
                                  ud_elements(ud_count) = choice_amesh(i)
                               endif
                          endif
                     endif

              enddo
              ! there are 2 cells left-right that can be used
            tot_count = lr_count+ud_count
            print*,tot_count,lr_count,ud_count

       !       prob = alpha/(float(lr_count)*alpha+float(tot_count-1))
            quake_ameshents(ele_inx) = choice_amesh(3);
            current_element = choice_amesh(i);
            quake_area = quake_area + amesh%area(quake_ameshents(ele_inx)) ! assume area is in m^2
            if (quake_area >= target_area) then
                     get_out = 1;
            endif
            ele_inx = ele_inx+1;
       enddo  ! end of loop over all possible nodes for next step
enddo

! Save section of fault used for earthquake
allocate(amesh%QuakeElem(ele_inx-1),amesh%QuakeNodes(inx-1))
amesh%QuakeElem = quake_ameshents(1:ele_inx-1)
amesh%QuakeNodes = quake_nodes(1:inx-1)
amesh%QuakeElemNo = ele_inx-1
amesh%QuakeNodesNo = inx-1
deallocate(quake_nodes,quake_ameshents)


end subroutine select_fault_zone
!==============================================================================
subroutine find_boundary(amesh,surf)
type(mesh),intent (inout) :: amesh
type(surface_type),intent(in) :: surf

integer,allocatable,dimension(:) :: haveit
integer :: ameshs(3)
integer :: bc_int,int
integer :: inx, kk,i,k
integer, allocatable,dimension(:) :: id,idx
integer :: i_nodes(3)
integer :: n,nlength,j
integer,allocatable,dimension(:,:) :: j_nodes,candidate_faces
integer,allocatable,dimension(:) :: int_amesh,bc_amesh,bc_nodes,common_node,id_keep
logical :: not_surface,bc_test
integer :: elem_idx,node_idx,global_node,face_id,face,f_idx
inx = 0
bc_int = 0
int = 0
node_idx = 0
allocate(bc_nodes(amesh%Nnodes),int_amesh(amesh%Ncells),bc_amesh(amesh%QuakeElemNo))

bc_nodes = 0
bc_amesh = 0

do i = 1,amesh%QuakeElemNo
  ! ameshs = neighbouring elements to QuakeElem(i)
       ameshs = amesh%EToE(amesh%QuakeElem(i),:)
       call isamember(ameshs,amesh%QuakeElem,id) !check if adjoining elements are in defined slipping zone
       if (id(1) == 0) then
          ameshs(1) = amesh%QuakeElem(i) ! elements is outside defined fault
       elseif (id(2) == 0) then
          ameshs(2) = amesh%QuakeElem(i) ! elements is outside defined fault
       elseif (id(3) == 0) then
          ameshs(3) = amesh%QuakeElem(i) ! elements is outside defined fault
       endif

       ! if ameshs contains an element = QuakeElem(i) it is a boundary point
       call find(ameshs,'==',amesh%QuakeElem(i),haveit)

       if (haveit(1) == -1) then  ! is interior elements
              int = int+1
              int_amesh(int) = amesh%QuakeElem(i)
       else      ! is a boundary element
         ! check to see if element is part of surface
            not_surface = .true.
            if (surf%exists) then
             do j = 1,3
               if (any(amesh%cell(amesh%QuakeElem(i),j).eq.amesh%SurfNodes) ) then
                     not_surface =.false.
               endif
             enddo
           endif
            if (not_surface) then
              bc_int = bc_int+1
              bc_amesh(bc_int) = amesh%QuakeElem(i)
!===== check faces for which nodes make up boundary======
! global nodal id's of current element
              i_nodes = amesh%cell(amesh%QuakeElem(i),:)
              if (allocated(idx)) deallocate(idx)
! find neighboring elements
              call find(ameshs,'/=',amesh%QuakeElem(i),idx)
              n  = size(idx) ! number of neighbouring elements that are in slipping zone
              if (allocated(candidate_faces)) deallocate(candidate_faces)
              allocate(candidate_faces(n*3,2))
              candidate_faces = 0
              inx = 1
              do kk = 1,n
                 f_idx = ameshs(idx(kk))
                 candidate_faces(inx,:) = amesh%cell(f_idx,1:2)
                 candidate_faces(inx+1,:) = amesh%cell(f_idx,2:3)
                 candidate_faces(inx+2,:) = (/ amesh%cell(f_idx,3), amesh%cell(f_idx,1) /)
                 inx = inx+3
              enddo

              do face = 1,3
               call check_face(i_nodes,face,candidate_faces,bc_test)
               if(bc_test) then
                 ! print*,'current element....',i
                 ! print*,'positive face',face
                 ! if(face == 1)  print*,'positive nodes',i_nodes(1),i_nodes(2)
                 ! if(face == 2)  print*,'positive nodes',i_nodes(2),i_nodes(3)
                 ! if(face == 3)  print*,'positive nodes',i_nodes(3),i_nodes(1)
                 !
                 ! print*,'all nodes',i_nodes
                 ! print*,'no of adjoining elements', n
                 ! do kk = 1,n*3
                 !    print*, 'near faces',candidate_faces(kk,:)
                 ! enddo
                 ! pause
                 if(face == 1) then
                   bc_nodes(node_idx+1) = i_nodes(1)
                   bc_nodes(node_idx+2) = i_nodes(2)
                   node_idx = node_idx+2
                 elseif(face == 2) then
                   bc_nodes(node_idx+1) = i_nodes(2)
                   bc_nodes(node_idx+2) = i_nodes(3)
                   node_idx = node_idx+2
                 elseif (face == 3) then
                   bc_nodes(node_idx+1) = i_nodes(3)
                   bc_nodes(node_idx+2) = i_nodes(1)
                   node_idx = node_idx+2
                 endif
               endif  ! if boundary condition
             enddo  ! loop on faces
        endif ! if not surface
      endif ! if interior / exterior point
  enddo


! remove duplicates
call erasedupli(bc_nodes,amesh%QuakeNodesNo,nlength)

amesh%QuakeBorder_NodesNo = nlength-1       ! number of vertices that make up boundary of fault

allocate(amesh%QuakeBorder_Nodes(amesh%QuakeBorder_NodesNo))
amesh%QuakeBorder_Nodes = bc_nodes(1:amesh%QuakeBorder_NodesNo)    ! vercties that make up boundary of fault
amesh%QuakeBorder_elemNo = bc_int     ! number of elements that make up boundary
!amesh%QuakeBorder_elemNo = elem_idx     ! number of elements that make up boundary

allocate(amesh%QuakeBorder_elem(amesh%QuakeBorder_elemNo))
amesh%QuakeBorder_elem = bc_amesh(1:bc_int)
!amesh%QuakeBorder_elem = bc_amesh(1:elem_idx)

end subroutine find_boundary
!==============================================================================
subroutine check_face(node,face,candidates,bc)
integer,dimension(3) :: node
integer,dimension(:,:) :: candidates
integer, dimension(2) :: c_face
integer :: face,i
logical :: bc

if (face == 1 ) then
 c_face = (/ node(1), node(2) /)
elseif (face == 2) then
  c_face = (/ node(2), node(3) /)
else
  c_face = (/ node(3), node(1) /)
endif

 bc = .true.
 do i =1,size(candidates,dim=1)
    if (c_face(1) == candidates(i,1).and.c_face(2) == candidates(i,2)) then
       bc = .false.
    elseif  (c_face(2) == candidates(i,1).and.c_face(1) == candidates(i,2)) then
      bc = .false.
    endif
enddo

! !allocate(candidate_faces(n*3,2))
!
! if(c_face(1) == candidates(1,1)) then
!     if(c_face(2) /= candidates(1,2).and.c_face(2) /= candidates(3,1)) then
!   ! face 1 is on border
!         bc = .true.
!         bc_nodes = c_face
!     endif
! endif
!
! if(c_face(1) == candidates(2,1)) then
!     if(c_face(2) /= candidates(2,2).and.c_face(2) /= candidates(1,1)) then
!   ! face 2 is on border
!     endif
! endif
!
! if(c_face(3) == candidates(3,1)) then
!     if(c_face(2) /= candidates(3,2).and.c_face(2) /= candidates(2,1)) then
!   ! face 3 is on border
!     endif
! endif
!
! do i =1,size(candidates,dim=1)
!       if (c_face(1) == candidates(i,1)) then
! ! check the other node
!           if (c_face(2) \= candidates(i,2))  then
! ! not sharing a face
!            bc = .true.
!            bc_idx = bc_idx+1
!            bc_test = c_face
!          endif
!       elseif (c_face(1) == candidates(i,2)) then
!         ! check the other node
!           if (c_face(2) \= candidates(i,1))  then
!         ! not sharing a face
!               bc = .true.
!               bc_idx = bc_idx+1
!               bc_test = c_face
!           endif
!       endif
! enddo

return
end subroutine check_face
!==============================================================================
subroutine intersect(A,B,v)
!  Find the common values to both A and B
!
!
integer,dimension(:),intent(in) :: A, B
integer, allocatable,dimension(:),intent(out) :: v
integer, allocatable,dimension(:) :: dummy
integer,allocatable,dimension(:) :: haveit

integer :: i, inx

allocate(dummy(size(A)))
inx = 0
do  i = 1,size(A)
       if (allocated(haveit)) deallocate(haveit)

       call find(B,'==',A(i),haveit)
       if ( haveit(1) > 0 ) then
              inx = inx +1
              dummy(inx) = A(i)
       endif
enddo

allocate(v(inx))
v = dummy(1:inx)
deallocate(dummy,haveit)

return
end subroutine intersect
!==============================================================================
!==============================================================================
subroutine isamember(A,B,idx)
!  Check if ameshents of A are in B
!      idx(i) = 1 element of A is present in B
!      idx(i) = 0 element of A is not present in B
integer,dimension(:),intent(in) :: A, B
integer, allocatable,dimension(:),intent(out) :: idx
integer,allocatable,dimension(:) :: haveit
integer :: i

allocate(idx(size(A)))

do  i = 1,size(A)
       call find(B,'==',A(i),haveit)
       if ( haveit(1) > 0 ) then
              idx(i) = 1
       else
              idx(i) = 0
       endif
enddo

return
end subroutine isamember
!==============================================================================
!==============================================================================
! locate returns an array mask based on what points are inside a set an given ameshent
!
FUNCTION locate2d(npts,ptx,pty,amesh) result(arr)
integer,intent(in) ::npts
REAL, INTENT(IN) :: ptx(npts),pty(npts)
REAL, INTENT(IN) :: amesh(4)
LOGICAL :: arr(npts)
real :: ex_min,ex_max,ey_min,ey_max
integer :: i


ex_min = amesh(1)
ex_max = amesh(2)
ey_min = amesh(3)
ey_max = amesh(4)

arr = .false.
do i = 1,npts
       if ((ptx(i) <= ex_max.and.ptx(i) >= ex_min).and. &
            (pty(i) <= ey_max.and.pty(i) >= ey_min)) then

             arr(i) = .true.
       endif
enddo

END FUNCTION locate2d
!!*******************************************************
!*******************************************************
subroutine reorder(a,n)
!######################################################
! Author : André Herrero
! Contact : andherit@gmail.com, andre.herrero@ingv.it
! Public Domain (CC0 1.0 Universal)
!######################################################
implicit none

integer n
real a(n)
integer i
real mem
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
!=============================================================================
!########################################################################################
!  ########################################################################################
subroutine interp_lin(x,y,xf,yf)
implicit none
real,dimension(2) :: x,y
real :: c,m,xf,yf

m = (y(2)-y(1))/(x(2)-x(1))
c = y(1)-m*x(1)

yf = m*xf+c

!end function interp_lin
end subroutine interp_lin
!  ########################################################################################
subroutine locate(n,xx,x,ans)
IMPLICIT NONE
integer,intent(in) ::n
REAL, INTENT(IN) :: xx(n)
REAL, INTENT(IN) :: x
INTEGER :: ans
INTEGER :: jl,jm,ju
LOGICAL :: ascnd

ascnd = (xx(n) >= xx(1))
jl=0
ju=n+1
do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
              jl=jm
       else
              ju=jm
       end if
end do
if (abs(x-xx(1)) < epsilon(x)) then
       !locate=1
       ans=1
else if (abs(x-xx(n)) < epsilon(x)) then
       !locate=n-1
       ans=n-1
else
       !locate=jl
       ans=jl
end if
END subroutine locate
!***********************************************************
!***********************************************************
subroutine set_seed()
implicit none
integer :: jobid,tot_rnos,iseed
integer,dimension(:),allocatable :: a_seed
character (len=30) :: en_var
integer :: i,s

! initialization of the random generator
open(10, file='/dev/urandom', access='stream', form='UNFORMATTED')
read(10) i
close(10)
!call getenv("PBS_ARRAYID",en_var) ! get environmental variable "test"
call getenv("PBS_ARRAY_INDEX",en_var) ! get environmental variable "test"
if (en_var == '') then
      jobid = 1
else
      read(en_var,'(I5)') jobid
endif
call random_seed(size=iseed)
allocate(a_seed(1:iseed))
call random_seed(get=a_seed)
call system_clock(count=s)
! random seed is choosen based on clock time
if (jobid== 1) then
  a_seed = abs( mod((s*181)*((i-83)*359), 104729) )
else
  a_seed = abs( mod((jobid*181)*((i-83)*359), 104729) )
endif

call random_seed(put=a_seed)     ! generate random seed
deallocate(a_seed)

!call random_seed     ! generate random seed

end subroutine set_seed
!***********************************************************
!----------------------------------------------------------------------
end module typedef
