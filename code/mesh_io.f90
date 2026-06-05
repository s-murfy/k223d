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
use, intrinsic :: iso_fortran_env, only: error_unit

implicit none

! integer(pin), parameter, private :: io_vtk = 13  ! reserved file unit for VTK

interface dumpnodeattributevtk
    module procedure idumpnodeattributevtk
    module procedure rdumpnodeattributevtk
end interface dumpnodeattributevtk

contains

!==============================================================================
subroutine read_vtk_mesh(filename,amesh, surface, pdf, src)
! Reads the VTK mesh file 'input.vtk'.
! Required: POINTS, CELLS (or POLYGONS for POLYDATA)
! Optional: CELL_TYPES (validation), velocity (cell), pdf (cell),
!           surface (node flag), rupt_time (node flag)
!==============================================================================
use LAT_time
  type(mesh),      intent(out)    :: amesh
  type(reflect),   intent(out)    :: surface
  type(pdfinputs), intent(out)    :: pdf
  type(source),    intent(out)    :: src     !
  character(*), intent(in) :: filename


  character(30) :: dataset_type
  real(pr), allocatable, dimension(:) :: field_r
  logical      :: found
  integer(pin) :: i,nuc_id,ios,io_vtk

  ! open(io_vtk, file='input.vtk', form='formatted', action='read', status='old')

  if (len_trim(filename) == 0 .or. trim(filename) == '-') then
     io_vtk = 6
  else
      io_vtk = 13 
     open(newunit=io_vtk, file=trim(filename), status='old', action='read', iostat=ios)
     if (ios /= 0) then
       write(error_unit,*) 'k223d: cannot open input file: ', trim(filename)
       stop
     endif
     inquire(unit=io_vtk, size=ios)
     if (ios == 0) then
       write(error_unit,*) 'k223d: input file is empty: ', trim(filename)
       stop
     endif
  endif
  ! --- Header validation ---
  call read_vtk_header(io_vtk, dataset_type)

  ! --- Required sections ---
  call find_and_read_points(io_vtk, amesh)
  call find_and_read_cells(io_vtk, amesh, dataset_type)

  ! --- Cell type validation ---
! Cell type validation — only meaningful for UNSTRUCTURED_GRID
! For POLYDATA the vertex count check in find_and_read_cells is sufficient
  if (trim(dataset_type) == 'UNSTRUCTURED_GRID') then
    call find_and_read_cell_types(io_vtk, amesh%Ncells, found)
    if (.not. found) write(error_unit,*) 'WARNING: CELL_TYPES absent, assuming all triangles'
  endif
  ! --- Mesh sanity checks ---
  call check_mesh(amesh)

  ! --- Optional cell fields ---
  allocate(field_r(amesh%Ncells))

  call find_and_read_cell_field(io_vtk, 'velocity', amesh%Ncells, found, field_r)
  if (found) then 
    write(error_unit,*) 'INFO: velocity field found in mesh file'
    allocate(src%velocity(amesh%Ncells))
    src%velocity = field_r 
  endif 

  call find_and_read_cell_field(io_vtk, 'pdf', amesh%Ncells, found, field_r)
  pdf%is_uniform = .true.

  if (found) then
    allocate(pdf%distrib(amesh%Ncells))
    pdf%distrib    = field_r
    pdf%is_uniform = .false.
    write(error_unit,*) 'INFO: PDF field found, pdf_type set to defined'
  endif

  deallocate(field_r)

  ! --- Optional node fields ---
  allocate(field_r(amesh%Nnodes))

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
    write(error_unit,*) 'INFO: surface rupture nodes found:', surface%Nnodes
  else
    surface%present = .false.
    surface%Nnodes  = 0
  endif

  ! Rupture Time: node with lowest field value is taken as nucleation point (i.e. 0) 
  call find_and_read_node_field(io_vtk, 'rupt_time', amesh%Nnodes, found, field_r)
  if (found) then
    nuc_id = minloc(field_r, 1)
    write(error_unit,*) 'INFO: rupt_time field found in mesh file'
    ! write(error_unit,*) minval(field_r)
    allocate(src%rupt_time(amesh%Nnodes))
    src%rupt_time = field_r
  endif

  deallocate(field_r)
  close(io_vtk)

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
  line = adjustl(line)
  if (line(1:1) /= '#') then
    write(error_unit,*) 'ERROR: not a valid VTK file - first line must start with #'
    stop
  endif

  ! Line 2: title - skip
  read(unit, '(a)') line

  ! Line 3: data format
  read(unit, '(a)') line
  if (trim(adjustl(line)) /= 'ASCII') then
    write(error_unit,*) 'ERROR: only ASCII VTK files are supported, found: ', trim(adjustl(line))
    stop
  endif

  ! Line 4: dataset type
  read(unit, '(a)') line
  if (index(line, 'UNSTRUCTURED_GRID') > 0) then
    dataset_type = 'UNSTRUCTURED_GRID'
  elseif (index(line, 'POLYDATA') > 0) then
    dataset_type = 'POLYDATA'
  else
    write(error_unit,*) 'ERROR: unsupported VTK dataset type:', trim(adjustl(line))
    stop
  endif

  write(error_unit,*) 'Reading VTK dataset type:', trim(dataset_type)

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
  ! write(error_unit,'(a,i10,a)') 'Reading  ', amesh%Nnodes, 'nodes'

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
  ! write(error_unit,'(a,i10,a)') 'Reading   ', ncells, 'cells'

  allocate(amesh%cell(ncells, 3))

  do i = 1, ncells
    read(unit, '(a)') line
    read(line, *, iostat=ios) nv, n1, n2, n3
    if (nv /= 3) then
      write(error_unit,'(a,i10,a,i10,a)') 'ERROR: cell', i, 'has', nv, 'vertices (only triangles supported)'
      stop
    endif
    amesh%cell(i, 1) = n1 + 1  ! VTK 0-based -> Fortran 1-based
    amesh%cell(i, 2) = n2 + 1
    amesh%cell(i, 3) = n3 + 1
  enddo

  amesh%Ncells = ncells

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
    write(error_unit,'(a,i10,a,i10)') 'ERROR: CELL_TYPES count', n_types, 'does not match Ncells', ncells
    stop
  endif

  do i = 1, n_types
    read(unit, *) cell_type
    if (cell_type /= 5) then
      write(error_unit,'(a,i10,a,i10,a)') 'ERROR: cell', i, 'has VTK type', cell_type, '(expected 5 = VTK_TRIANGLE)'
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
    write(error_unit,'(a,i10,a,i10)') 'ERROR: CELL_DATA count', n, 'does not match Ncells', ncells
    stop
  endif

  call find_field_in_section(unit, fieldname, found)
  if (.not. found) return

  read(unit, '(a)') line  ! skip LOOKUP_TABLE line
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
    write(error_unit,'(a,i10,a,i10)') 'ERROR: POINT_DATA count', n, 'does not match Nnodes', nnodes
    stop
  endif

  call find_field_in_section(unit, fieldname, found)
  if (.not. found) return

  read(unit, '(a)') line  ! skip LOOKUP_TABLE line
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
    write(error_unit,'(a,a)') 'ERROR: required VTK section not found: ', trim(keyword)
    stop
  endif

end subroutine find_keyword

!==============================================================================
subroutine find_field_in_section(unit, fieldname, found)
! Scans forward from the current file position for a SCALARS field with
! the given name. Stops at the next POINT_DATA or CELL_DATA section header
! so that search is contained within the current data section.
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
    if (ios /= 0) exit    ! keyword not in the file 
    line = adjustl(line)
    if (index(line, 'POINT_DATA') == 1) exit  ! left current section
    if (index(line, 'CELL_DATA')  == 1) exit  ! left current section
    if (index(line, 'SCALARS') == 1) then
      read(line, *, iostat=ios) kw, fname  ! iostat guards against malformed line
      if (ios /= 0) cycle                  ! skip and continue scanning
      if (trim(fname) == trim(fieldname)) then
        found = .true.
        return
      endif
    endif
  enddo

end subroutine find_field_in_section

!==============================================================================
subroutine check_mesh(amesh)
! Validates mesh integrity. Checks:
!   - at least 3 nodes and 1 cell exist
!   - all cell vertex indices are within [1, Nnodes]
! Collects all errors before stopping so all problems are reported at once.
!==============================================================================
  type(mesh), intent(in) :: amesh

  integer(pin) :: i, j
  integer(pin) :: error_count

  error_count = 0

  if (amesh%Nnodes < 3) then
    write(error_unit,'(a,i10,a,i10,a)') 'ERROR: mesh has fewer than 3 nodes (', amesh%Nnodes, ')'
    error_count = error_count + 1
  endif

  if (amesh%Ncells < 1) then
    write(error_unit,'(a)') 'ERROR: mesh has no cells'
    error_count = error_count + 1
  endif

  do i = 1, amesh%Ncells
    do j = 1, 3
      if (amesh%cell(i,j) < 1 .or. amesh%cell(i,j) > amesh%Nnodes) then
        write(error_unit,'(a,i10,a,i10,a,i10,a,i10,a)') 'ERROR: cell', i, 'vertex', j, 'index', amesh%cell(i,j), &
                   'out of range [1,', amesh%Nnodes, ']'
        error_count = error_count + 1
      endif
    enddo
  enddo

  if (error_count > 0) stop

  write(error_unit,'(a,i10,a,i10,a)') 'Mesh OK:  ', amesh%Nnodes, ' nodes,', amesh%Ncells, ' cells'

end subroutine check_mesh
!==============================================================================
subroutine write_vtk_result(filename,amesh, src, pdf)
use mesh_geom
    type(mesh),          intent(in) :: amesh
    type(source),        intent(in) :: src
    type(pdfinputs),     intent(in) :: pdf
    ! character(*),        intent(in) :: outfile
    character(*),        intent(in) :: filename

    integer(pin), allocatable, dimension(:) :: surface_flag
    integer(pin) :: dev, i,ios
    ! dev = 11

    ! open(dev, file=trim(outfile), form='formatted')
    if (len_trim(filename) == 0 .or. trim(filename) == '-') then
      dev = 6
    else
      dev = 11
     open(newunit=dev, file=trim(filename), status='replace', action='write', iostat=ios)
     if (ios /= 0) stop 'k223d: cannot open output file'
    endif


    ! Mesh geometry
    ! call dumpmeshvtk(dev, amesh) ! this outputs in POLYDATA format, can't be read by Meshio
    call dumpmeshvtk_meshio(dev, amesh) ! this version writes out in polydata which can be read by Meshio in python 
    ! CELL_DATA: slip always, velocity and pdf if present
    call dumpcellattributevtk(dev, amesh, src%slip,     'slip',     .true.)
    call dumpcellattributevtk(dev, amesh, pdf%distrib,  'pdf',      .false.)
    if (allocated(src%velocity)) &
        call dumpcellattributevtk(dev, amesh, src%velocity, 'velocity', .false.)

    ! POINT_DATA: rupture time if computed
    if (allocated(src%rupt_time)) &
        call dumpnodeattributevtk(dev, amesh, src%rupt_time, 'rupt_time', .true.)

    close(dev)

end subroutine write_vtk_result
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
!==============================================================================
subroutine dumpmeshvtk_meshio(dev, amesh)
    type(mesh),  intent(in) :: amesh
    integer(pin), intent(in) :: dev
    integer(pin) :: i, j

    write(dev,'(a)') '# vtk DataFile Version 2.0'
    write(dev,'(a)') 'k223d output'
    write(dev,'(a)') 'ASCII'
    write(dev,'(a)') 'DATASET UNSTRUCTURED_GRID'

    if (pr == 4) then
        write(dev,'(a,i10,a)') 'POINTS ', amesh%Nnodes, ' float'
    else
        write(dev,'(a,i10,a)') 'POINTS ', amesh%Nnodes, ' double'
    endif

    do i = 1, amesh%Nnodes
        write(dev,*) amesh%px(i), amesh%py(i), amesh%pz(i)
    enddo

    write(dev,'(a,2i10)') 'CELLS ', amesh%Ncells, amesh%Ncells*4
    do i = 1, amesh%Ncells
        write(dev,'(i2,3i10)') 3, (amesh%cell(i,j)-1, j=1,3)
    enddo

    write(dev,'(a,i10)') 'CELL_TYPES ', amesh%Ncells
    do i = 1, amesh%Ncells
        write(dev,*) 5
    enddo

end subroutine dumpmeshvtk_meshio
!==============================================================================
subroutine dumpcellattributevtk(dev, amesh, field, attname, init)
    type(mesh),                      intent(in) :: amesh
    real(pr), dimension(amesh%Ncells), intent(in) :: field
    integer(pin),                    intent(in) :: dev
    character(*),                    intent(in) :: attname
    logical,                         intent(in) :: init
    integer(pin) :: i

    if (init) write(dev, '(a,i10)') 'CELL_DATA ', amesh%Ncells
    write(dev, '(a)') 'SCALARS '//trim(attname)//' double 1'
    write(dev, '(a)') 'LOOKUP_TABLE default'
    do i = 1, amesh%Ncells
        write(dev, *) field(i)
    enddo
end subroutine dumpcellattributevtk
!==============================================================================
subroutine idumpnodeattributevtk(dev, amesh, field, attname, init)
    type(mesh),                        intent(in) :: amesh
    integer(pin), dimension(amesh%Nnodes), intent(in) :: field
    integer(pin),                      intent(in) :: dev
    character(len=*),                  intent(in) :: attname
    logical,                           intent(in) :: init
    integer(pin) :: i

    if (init) write(dev,'(a,i10)') 'POINT_DATA ', amesh%Nnodes
    write(dev,'(a)') 'SCALARS '//trim(attname)//' int 1'
    write(dev,'(a)') 'LOOKUP_TABLE default'
    do i = 1, amesh%Nnodes
        write(dev,*) field(i)
    enddo
end subroutine idumpnodeattributevtk
!==============================================================================
subroutine rdumpnodeattributevtk(dev, amesh, field, attname, init)
    type(mesh),                      intent(in) :: amesh
    real(pr), dimension(amesh%Nnodes), intent(in) :: field
    integer(pin),                    intent(in) :: dev
    character(len=*),                intent(in) :: attname
    logical,                         intent(in) :: init
    integer(pin) :: i

    if (init) write(dev,'(a,i10)') 'POINT_DATA ', amesh%Nnodes
    write(dev,'(a)') 'SCALARS '//trim(attname)//' double 1'
    write(dev,'(a)') 'LOOKUP_TABLE default'
    do i = 1, amesh%Nnodes
        write(dev,*) field(i)
    enddo
end subroutine rdumpnodeattributevtk
!==============================================================================
end module mesh_io