module mesh_io
!==============================================================================
! Mesh I/O module for k223d
!
! Subroutines:
!   read_vtk_mesh          - main entry point for reading input.vtk
!   read_vtk_header        - validates 4-line VTK header, returns dataset type
!   find_and_read_points   - reads POINTS section (required)
!   find_and_read_cells    - reads CELLS/POLYGONS section (required)
!   find_and_read_cell_types - validates CELL_TYPES section (optional)
!   find_and_read_cell_field - reads a named scalar from CELL_DATA (optional)
!   find_and_read_node_field - reads a named scalar from POINT_DATA (optional)
!   check_mesh             - validates node/cell counts and index bounds
!   find_keyword           - helper: rewind and scan for a VTK section keyword
!   find_field_in_section  - helper: scan forward for a named SCALARS field
!
! Note: VTK legacy format uses 0-based node indexing; this module converts
!       to 1-based on read. If your VTK files are already 1-based (non-standard)
!       remove the +1 in find_and_read_cells.
!
! To be added: read_input, write_vtk_result, write_gmt_result, write_batch_result
!==============================================================================
use generic
use LAT_mesh
use LAT_source
implicit none

integer(pin), parameter, private :: io_vtk = 13  ! reserved file unit for VTK

contains

!==============================================================================
subroutine read_vtk_mesh(amesh, surface, pdf, nuc_id, rupt_time, velocity)
! Reads the VTK mesh file 'input.vtk'.
! Required: POINTS, CELLS (or POLYGONS for POLYDATA)
! Optional cell fields (CELL_DATA): 'velocity', 'pdf'
! Optional node fields (POINT_DATA): 'time', 'surface'
!
! Expected field names (see README):
!   velocity  - rupture velocity per cell (double)
!   pdf       - slip PDF per cell (double)
!   time      - initial rupture time per node (double):
!                 0.0    at the nucleation node
!                 1.e32  at all other nodes (to be computed)
!   surface   - surface rupture flag per node (int): 1 = surface node
!
! Computation mode is determined by field presence:
!   velocity + time both present -> slip and rupture time computed
!   either absent               -> slip only
!
! rupt_time and velocity are only allocated if their fields are found.
! The caller checks allocated(rupt_time) .and. allocated(velocity)
! to decide whether to run calc_rupt_front.
! nuc_id is set to 0 when time field is absent (slip-only mode).
!==============================================================================
  type(mesh),                          intent(out) :: amesh
  type(reflect),                       intent(out) :: surface
  type(pdfinputs),                   intent(inout) :: pdf
  integer(pin),                        intent(out) :: nuc_id
  real(pr), allocatable, dimension(:), intent(out) :: rupt_time
  real(pr), allocatable, dimension(:), intent(out) :: velocity

  character(30) :: dataset_type
  real(pr), allocatable, dimension(:) :: field_r
  logical      :: found
  integer(pin) :: i

  open(io_vtk, file='input.vtk', form='formatted', action='read', status='old')

  ! --- Header validation ---
  call read_vtk_header(io_vtk, dataset_type)

  ! --- Required sections ---
  call find_and_read_points(io_vtk, amesh)
  call find_and_read_cells(io_vtk, amesh, dataset_type)

  ! --- Cell type validation ---
  ! POLYDATA triangle validity is already confirmed by the vertex count check
  ! in find_and_read_cells. CELL_TYPES is only meaningful for UNSTRUCTURED_GRID.
  if (trim(dataset_type) == 'UNSTRUCTURED_GRID') then
    call find_and_read_cell_types(io_vtk, amesh%Ncells, found)
    if (.not. found) write(*,*) 'WARNING: CELL_TYPES absent, assuming all triangles'
  endif

  ! --- Optional cell fields ---
  allocate(field_r(amesh%Ncells))

  ! Velocity: only allocate output if field is present
  call find_and_read_cell_field(io_vtk, 'velocity', amesh%Ncells, found, field_r)
  if (found) then
    allocate(velocity(amesh%Ncells))
    velocity = field_r
  endif

  call find_and_read_cell_field(io_vtk, 'pdf', amesh%Ncells, found, field_r)
  if (found) then
    allocate(pdf%g_pdf(amesh%Ncells))
    pdf%g_pdf    = field_r
    pdf%pdf_type = 'defined'
  endif

  deallocate(field_r)

  ! --- Optional node fields ---
  allocate(field_r(amesh%Nnodes))

  ! Time: initial rupture conditions. Only allocate rupt_time if field present.
  ! Nucleation node has value 0.0; all others 1.e32 (infinity).
  ! nuc_id is set to 0 when absent — slip-only mode, nucleation not required.
  call find_and_read_node_field(io_vtk, 'time', amesh%Nnodes, found, field_r)
  if (found) then
    allocate(rupt_time(amesh%Nnodes))
    rupt_time = field_r
    nuc_id    = minloc(rupt_time, 1)
  else
    nuc_id = 0
  endif

  ! Surface rupture: nodes with field value > 0 are surface nodes
  call find_and_read_node_field(io_vtk, 'surface', amesh%Nnodes, found, field_r)
  if (found) then
    surface%present = .true.
    surface%Nnodes  = count(field_r > 0.5_pr)
    allocate(surface%nodes(surface%Nnodes))
    surface%Nnodes = 0
    do i = 1, amesh%Nnodes
      if (field_r(i) > 0.5_pr) then
        surface%Nnodes = surface%Nnodes + 1
        surface%nodes(surface%Nnodes) = i
      endif
    enddo
  else
    surface%present = .false.
    surface%Nnodes  = 0
  endif

  deallocate(field_r)
  close(io_vtk)

  ! --- Report computation mode ---
  if (allocated(rupt_time) .and. allocated(velocity)) then
    write(*,*) 'MODE: slip and rupture time will be computed'
    write(*,*) '      nucleation node:', nuc_id
    write(*,*) '      surface rupture:', surface%present
  else
    write(*,*) 'MODE: slip only (velocity or time field absent from mesh file)'
  endif

end subroutine read_vtk_mesh

!==============================================================================
subroutine read_vtk_header(unit, dataset_type)
! Reads and validates the mandatory 4-line VTK legacy ASCII header.
! Line 1: # vtk DataFile Version X.X
! Line 2: Title (any string, ignored)
! Line 3: ASCII
! Line 4: DATASET UNSTRUCTURED_GRID or DATASET POLYDATA
!==============================================================================
  integer(pin),  intent(in)  :: unit
  character(30), intent(out) :: dataset_type

  character(100) :: line

  ! Line 1: magic string
  read(unit, '(a)') line
  if (adjustl(line)(1:1) /= '#') then
    write(*,*) 'ERROR: not a valid VTK file - first line must start with #'
    stop
  endif

  ! Line 2: title - skip
  read(unit, '(a)') line

  ! Line 3: data format
  read(unit, '(a)') line
  if (trim(adjustl(line)) /= 'ASCII') then
    write(*,*) 'ERROR: only ASCII VTK files are supported, found: ', trim(adjustl(line))
    stop
  endif

  ! Line 4: dataset type
  read(unit, '(a)') line
  if (index(line, 'UNSTRUCTURED_GRID') > 0) then
    dataset_type = 'UNSTRUCTURED_GRID'
  elseif (index(line, 'POLYDATA') > 0) then
    dataset_type = 'POLYDATA'
  else
    write(*,*) 'ERROR: unsupported VTK dataset type:', trim(adjustl(line))
    stop
  endif

  write(*,*) 'Reading VTK dataset type:', trim(dataset_type)

end subroutine read_vtk_header

!==============================================================================
subroutine find_and_read_points(unit, amesh)
! Locates the POINTS section and reads all node (x,y,z) coordinates.
!==============================================================================
  integer(pin), intent(in)    :: unit
  type(mesh),   intent(inout) :: amesh

  character(100) :: line
  character(20)  :: kw, datatype
  integer(pin)   :: i

  call find_keyword(unit, 'POINTS', line, .true.)
  read(line, *) kw, amesh%Nnodes, datatype
  write(*,*) 'Reading', amesh%Nnodes, 'nodes'

  allocate(amesh%px(amesh%Nnodes), amesh%py(amesh%Nnodes), amesh%pz(amesh%Nnodes))
  do i = 1, amesh%Nnodes
    read(unit, *) amesh%px(i), amesh%py(i), amesh%pz(i)
  enddo

end subroutine find_and_read_points

!==============================================================================
subroutine find_and_read_cells(unit, amesh, dataset_type)
! Locates CELLS (UNSTRUCTURED_GRID) or POLYGONS (POLYDATA) and reads
! triangle connectivity. Cell vertex indices are converted from VTK 0-based
! to Fortran 1-based. Non-triangular cells trigger a fatal error.
!==============================================================================
  integer(pin), intent(in)    :: unit
  type(mesh),   intent(inout) :: amesh
  character(*), intent(in)    :: dataset_type

  character(100) :: line
  character(20)  :: kw
  integer(pin)   :: ncells, nsize, i, nv, n1, n2, n3
  integer        :: ios

  if (trim(dataset_type) == 'POLYDATA') then
    call find_keyword(unit, 'POLYGONS', line, .true.)
  else
    call find_keyword(unit, 'CELLS', line, .true.)
  endif

  read(line, *) kw, ncells, nsize
  write(*,*) 'Reading', ncells, 'cells'

  allocate(amesh%cell(ncells, 3))

  do i = 1, ncells
    read(unit, '(a)') line
    read(line, *, iostat=ios) nv, n1, n2, n3
    if (nv /= 3) then
      write(*,*) 'ERROR: cell', i, 'has', nv, 'vertices (only triangles supported)'
      stop
    endif
    amesh%cell(i, 1) = n1 + 1  ! VTK 0-based -> Fortran 1-based
    amesh%cell(i, 2) = n2 + 1
    amesh%cell(i, 3) = n3 + 1
    ! Validate converted indices immediately — catches 0-vs-1 indexing errors
    if (any(amesh%cell(i,:) < 1) .or. any(amesh%cell(i,:) > amesh%Nnodes)) then
      write(*,*) 'ERROR: cell', i, 'has out-of-range node indices:', amesh%cell(i,:)
      write(*,*) '       Valid range: [1,', amesh%Nnodes, ']'
      write(*,*) '       Check whether mesh uses 0-based or 1-based node indexing'
      stop
    endif
  enddo

  amesh%Ncells = ncells
  write(*,*) 'Mesh read OK:', amesh%Nnodes, 'nodes,', amesh%Ncells, 'cells'

end subroutine find_and_read_cells

!==============================================================================
subroutine find_and_read_cell_types(unit, ncells, found)
! Locates CELL_TYPES and verifies all values are 5 (VTK_TRIANGLE).
! Absence of this section produces a warning rather than an error since
! some VTK writers omit it. Count mismatch with ncells is always fatal.
!==============================================================================
  integer(pin), intent(in)  :: unit
  integer(pin), intent(in)  :: ncells
  logical,      intent(out) :: found

  character(100) :: line
  character(20)  :: kw
  integer(pin)   :: i, cell_type, n_types

  call find_keyword(unit, 'CELL_TYPES', line, .false., found)
  if (.not. found) return

  read(line, *) kw, n_types
  if (n_types /= ncells) then
    write(*,*) 'ERROR: CELL_TYPES count', n_types, 'does not match Ncells', ncells
    stop
  endif

  do i = 1, n_types
    read(unit, *) cell_type
    if (cell_type /= 5) then
      write(*,*) 'ERROR: cell', i, 'has VTK type', cell_type, '(expected 5 = VTK_TRIANGLE)'
      stop
    endif
  enddo

end subroutine find_and_read_cell_types

!==============================================================================
subroutine find_and_read_cell_field(unit, fieldname, ncells, found, field)
! Locates a named SCALARS field within the CELL_DATA section and reads it.
! Validates CELL_DATA count against ncells before reading.
! Returns found=.false. without error if section or field is absent.
!==============================================================================
  integer(pin),                intent(in)  :: unit
  character(*),                intent(in)  :: fieldname
  integer(pin),                intent(in)  :: ncells
  logical,                     intent(out) :: found
  real(pr), dimension(ncells), intent(out) :: field

  character(100) :: line
  character(20)  :: kw
  integer(pin)   :: n, i
  logical        :: section_found

  found = .false.
  call find_keyword(unit, 'CELL_DATA', line, .false., section_found)
  if (.not. section_found) return

  read(line, *) kw, n
  if (n /= ncells) then
    write(*,*) 'ERROR: CELL_DATA count', n, 'does not match Ncells', ncells
    stop
  endif

  call find_field_in_section(unit, fieldname, found)
  if (.not. found) return

  ! Scan forward to LOOKUP_TABLE rather than assuming it is the very next line.
  ! VTK writers sometimes insert blank lines between SCALARS and LOOKUP_TABLE.
  call skip_to_lookup_table(unit)

  do i = 1, ncells
    read(unit, *) field(i)
  enddo

end subroutine find_and_read_cell_field

!==============================================================================
subroutine find_and_read_node_field(unit, fieldname, nnodes, found, field)
! Locates a named SCALARS field within the POINT_DATA section and reads it.
! Validates POINT_DATA count against nnodes before reading.
! Returns found=.false. without error if section or field is absent.
!==============================================================================
  integer(pin),                intent(in)  :: unit
  character(*),                intent(in)  :: fieldname
  integer(pin),                intent(in)  :: nnodes
  logical,                     intent(out) :: found
  real(pr), dimension(nnodes), intent(out) :: field

  character(100) :: line
  character(20)  :: kw
  integer(pin)   :: n, i
  logical        :: section_found

  found = .false.
  call find_keyword(unit, 'POINT_DATA', line, .false., section_found)
  if (.not. section_found) return

  read(line, *) kw, n
  if (n /= nnodes) then
    write(*,*) 'ERROR: POINT_DATA count', n, 'does not match Nnodes', nnodes
    stop
  endif

  call find_field_in_section(unit, fieldname, found)
  if (.not. found) return

  ! Scan forward to LOOKUP_TABLE rather than assuming it is the very next line.
  ! VTK writers sometimes insert blank lines between SCALARS and LOOKUP_TABLE.
  call skip_to_lookup_table(unit)

  do i = 1, nnodes
    read(unit, *) field(i)
  enddo

end subroutine find_and_read_node_field

!==============================================================================
subroutine find_keyword(unit, keyword, line, fatal, found)
! Rewinds the file and scans line by line for a line whose first non-blank
! token matches keyword. Returns the matching line for the caller to parse.
! If fatal=.true. and keyword is not found, stops with an error.
! If fatal=.false., returns result via the optional found argument.
!==============================================================================
  integer(pin),            intent(in)            :: unit
  character(*),            intent(in)            :: keyword
  character(100),          intent(out)           :: line
  logical,                 intent(in)            :: fatal
  logical,       optional, intent(out)           :: found

  integer :: ios
  logical :: hit

  hit = .false.
  rewind(unit)

  do
    read(unit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    if (index(adjustl(line), trim(keyword)) == 1) then
      hit = .true.
      exit
    endif
  enddo

  if (present(found)) found = hit

  if (.not. hit .and. fatal) then
    write(*,*) 'ERROR: required VTK section not found: ', trim(keyword)
    stop
  endif

end subroutine find_keyword

!==============================================================================
subroutine find_field_in_section(unit, fieldname, found)
! Scans forward from the current file position for a SCALARS field with
! the given name. Stops at the next POINT_DATA or CELL_DATA section header
! so that search is contained within the current data section.
! The iostat guard on the internal read handles any malformed SCALARS lines
! that are missing the field name token.
!==============================================================================
  integer(pin), intent(in)  :: unit
  character(*), intent(in)  :: fieldname
  logical,      intent(out) :: found

  character(100) :: line
  character(20)  :: kw, fname
  integer        :: ios

  found = .false.
  do
    read(unit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    line = adjustl(line)
    if (index(line, 'POINT_DATA') == 1) exit  ! left current section
    if (index(line, 'CELL_DATA')  == 1) exit  ! left current section
    if (index(line, 'SCALARS') == 1) then
      read(line, *, iostat=ios) kw, fname      ! iostat guards malformed lines
      if (ios /= 0) cycle                      ! missing name token — skip line
      if (trim(fname) == trim(fieldname)) then
        found = .true.
        return
      endif
    endif
  enddo

end subroutine find_field_in_section

!==============================================================================
subroutine skip_to_lookup_table(unit)
! Scans forward from current position until a LOOKUP_TABLE line is found,
! then returns with the file positioned immediately after it so the next
! read returns the first data value. Handles blank lines between SCALARS
! and LOOKUP_TABLE which some VTK writers insert.
!==============================================================================
  integer(pin), intent(in) :: unit

  character(100) :: line
  integer        :: ios

  do
    read(unit, '(a)', iostat=ios) line
    if (ios /= 0) then
      write(*,*) 'ERROR: LOOKUP_TABLE not found after SCALARS header'
      stop
    endif
    if (index(adjustl(line), 'LOOKUP_TABLE') == 1) return
  enddo

end subroutine skip_to_lookup_table

!==============================================================================
end module mesh_io