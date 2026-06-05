module k223d_capi
!==============================================================================
! C API wrapper for k223d.
! Exposes compute_source to Python via ctypes (or any C-compatible caller).
! All arguments use iso_c_binding types — no Fortran derived types at boundary.
!
! Python usage (via k223d.py ctypes wrapper):
!
!   import k223d
!
!   # Slip only (no velocity or time arrays in mesh)
!   slip = k223d.compute_source(nodes, cells, mw=6.5)
!
!   # Slip + rupture time
!   slip, time = k223d.compute_source(nodes, cells, mw=6.5,
!                                      velocity=vel, time=t0)
!
! The Python wrapper (k223d.py) handles numpy array layout and type conversion.
!==============================================================================
use iso_c_binding
use LAT_mesh
use LAT_source
use k223d_core
implicit none

contains

!==============================================================================
subroutine c_compute_source(                                          &
    nnodes, ncells,                                                   &
    px, py, pz,                                                       &
    cells,                                                            &
    mw, mu, na, rmin, rmax,                                           &
    n_surface, surface_nodes,                                         &
    has_velocity, velocity,                                           &
    has_time,     time,                                               &
    has_pdf,      pdf_in,                                             &
    slip,                                                             &
    compute_time_flag)                                                &
    bind(c, name='k223d_compute_source')
!
! Arguments:
!   nnodes, ncells        mesh dimensions
!   px, py, pz            node coordinates [nnodes]
!   cells                 cell-to-node connectivity, 1-based [3, ncells]
!   mw                    moment magnitude (required)
!   mu                    shear modulus, Pa (<=0 -> use default 3.0e10)
!   na                    number of sub-events (<=0 -> use default 5000)
!   rmin                  min sub-event radius, m (<=0 -> derive from mesh)
!   rmax                  max sub-event radius, m (<=0 -> derive from mesh)
!   n_surface             number of surface rupture nodes (0 if none)
!   surface_nodes         surface node indices, 1-based [n_surface]
!   has_velocity          1 if velocity is provided, 0 otherwise
!   velocity              rupture velocity per cell [ncells] (intent in)
!   has_time              1 if time is provided, 0 otherwise
!   time                  dual purpose [nnodes] (intent inout):
!                           on entry: initial conditions (0.0 at nucleation,
!                                     1.e32 elsewhere) when has_time=1
!                           on exit:  computed rupture time when
!                                     compute_time_flag=1
!   has_pdf               1 if pdf_in is provided, 0 for uniform
!   pdf_in                slip probability density per cell [ncells] (intent in)
!   slip                  computed total slip per cell [ncells] (intent out)
!   compute_time_flag     1 if rupture time was computed, 0 if slip only
!==============================================================================
  integer(c_int), value, intent(in)    :: nnodes, ncells
  real(c_double),        intent(in)    :: px(nnodes), py(nnodes), pz(nnodes)
  integer(c_int),        intent(in)    :: cells(3, ncells)  ! column-major for C
  real(c_double), value, intent(in)    :: mw, mu, rmin, rmax
  integer(c_int), value, intent(in)    :: na
  integer(c_int), value, intent(in)    :: n_surface
  integer(c_int),        intent(in)    :: surface_nodes(n_surface)
  integer(c_int), value, intent(in)    :: has_velocity, has_time
  real(c_double),        intent(in)    :: velocity(ncells)
  real(c_double),        intent(inout) :: time(nnodes)
  integer(c_int), value, intent(in)    :: has_pdf
  real(c_double),        intent(in)    :: pdf_in(ncells)
  real(c_double),        intent(out)   :: slip(ncells)
  integer(c_int),        intent(out)   :: compute_time_flag

  type(mesh)        :: amesh
  type(reflect)     :: surface
  type(model_param) :: model
  type(pdfinputs)   :: pdf
  type(source)      :: src
  integer(pin)      :: i

  ! --- Construct amesh from C arrays ---
  amesh%Nnodes = nnodes
  amesh%Ncells = ncells
  allocate(amesh%px(nnodes), amesh%py(nnodes), amesh%pz(nnodes))
  amesh%px = real(px, pr)
  amesh%py = real(py, pr)
  amesh%pz = real(pz, pr)
  allocate(amesh%cell(ncells, 3))
  do i = 1, ncells
    amesh%cell(i, 1) = cells(1, i)  ! transpose: C row-major -> Fortran col-major
    amesh%cell(i, 2) = cells(2, i)
    amesh%cell(i, 3) = cells(3, i)
  enddo

  ! --- Construct surface from optional surface node list ---
  if (n_surface > 0) then
    surface%present = .true.
    surface%Nnodes  = n_surface
    allocate(surface%nodes(n_surface))
    surface%nodes = surface_nodes
  else
    surface%present = .false.
    surface%Nnodes  = 0
  endif

  ! --- Set model parameters (sentinel -1 triggers defaults in compute_source) ---
  model%mw   = real(mw, pr)
  model%mu   = merge(real(mu,   pr), -1._pr, mu   > 0._c_double)
  model%na   = merge(na,              -1,     na   > 0)
  model%rmin = merge(real(rmin, pr), -1._pr, rmin > 0._c_double)
  model%rmax = merge(real(rmax, pr), -1._pr, rmax > 0._c_double)

  ! Apply defaults for any unset parameters
  if (model%mu   < 0._pr) model%mu   = 3.0e10_pr
  if (model%na   < 0)     model%na   = 5000
  if (model%rmin < 0._pr) model%rmin = 2.0_pr
  if (model%rmax < 0._pr) model%rmax = 0.35_pr

  ! --- PDF: user-supplied or uniform ---
  pdf%is_uniform = .true.  ! assume that pdf is uniform unless stated otherwise 
  if (has_pdf == 1) then
    allocate(pdf%distrib(ncells))
    pdf%distrib    = real(pdf_in, pr)
    pdf%is_uniform = .false.
  endif

  ! --- Load optional velocity and initial time arrays ---
  if (has_velocity == 1) then
    allocate(src%velocity(ncells))
    src%velocity = real(velocity, pr)
  endif
  if (has_time == 1) then
    allocate(src%rupt_time(nnodes))
    src%rupt_time = real(time, pr)
  endif

  ! --- Run shared computational core ---
  call compute_source(amesh, surface, model, pdf, src)

  ! --- Return results ---
  ! slip is always computed
  slip = real(src%slip, c_double)

  ! time array is written back if rupture time was computed;
  ! compute_time_flag tells the caller whether to use it
  if (allocated(src%rupt_time)) then
    time              = real(src%rupt_time, c_double)
    compute_time_flag = 1
  else
    compute_time_flag = 0
  endif

end subroutine c_compute_source

!==============================================================================
end module k223d_capi
