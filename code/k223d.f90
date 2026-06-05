program k223d
!==============================================================================
! k223d standalone executable.
! Reads input.vtk, runs compute_source, writes output.vtk.
!
! Usage:
!   k223d mw [mu=value] [na=value] [rmin=value] [rmax=value]
!
! Examples:
!   k223d mw=6.5
!   k223d mw=6.5 mu=3.3e10
!   k223d mw=6.5 na=1000 rmax=25.0
!
! Optional parameters default to:
!   mu   = 3.0e10 Pa
!   na   = 5000
!   rmin = equivalent radius of mean cell area (from mesh)
!   rmax = half fault width estimated from mesh geometry
!   in = input file name
!   out = output file name 
! PDF type is determined from the mesh file:
!   'pdf' field present     -> defined PDF
!   'pdf' field absent      -> uniform PDF
!
! Computation mode is determined from the mesh file:
!   'velocity' + 'time' present -> slip and rupture time computed
!   either absent               -> slip only
!==============================================================================
use mesh_io
use LAT_source
use k223d_core
implicit none
 
type(mesh)        :: amesh
type(reflect)     :: surface
type(model_param) :: model
type(pdfinputs)   :: pdf
type(source)      :: src
character(len=256) :: infile, outfile
integer(pin)      :: nuc_id
 
! Initialise optional model parameters to sentinel (-1 = unset)
call init_model(model)
 
! Parse mw from command line and any optional key=value overrides
call parse_args(model,infile,outfile)
 
! Read mesh, optional PDF, velocity initial conditions and surface flag
call read_vtk_mesh(trim(infile),amesh, surface, pdf, src)

! Run shared computational core
call compute_source(amesh, surface, model, pdf, src)
 
! Write results to outfile
call write_vtk_result(trim(outfile),amesh, src,pdf)


end program k223d 
!==============================================================================
subroutine init_model(model)
! Initialises all optional model parameters to -1 (sentinel = unset).
! mw is left at 0; parse_args will set it or stop with an error.
!==============================================================================
  use LAT_source
  implicit none
  type(model_param), intent(out) :: model
 ! default values are provided 
  model%mw    =  0._pr
  model%mu    =  3.0e10_pr 
  model%na    =  5000
  model%rmin  =  2.0_pr
  model%rmax  =  0.35_pr
 
end subroutine init_model 
!==============================================================================
subroutine parse_args(model,infile,outfile)
! Reads command-line arguments:
!   Argument 1 (required): mw
!   Remaining arguments:   key=value pairs for mu, na, rmin, rmax
!==============================================================================
  use LAT_source
  use generic
  use forparse
  use, intrinsic :: iso_fortran_env, only: error_unit

  implicit none
  type(model_param), intent(inout) :: model
  character(len=256), intent(out) :: infile, outfile

  infile = '-'
  outfile = '-'

if (parse_arg('mw', model%mw) /= PARSE_OK) then
    write(error_unit,*) 'ERROR: mw is required'
    write(error_unit,*) 'Usage: k223d mw=6.5 [mu=3.3e10] [na=5000] [rmin=...] [rmax=...] [in=....] [out=....]'
    stop
endif

if (parse_arg('mu',   model%mu)   == PARSE_TYPE_ERROR) stop
if (parse_arg('na',   model%na)   == PARSE_TYPE_ERROR) stop
if (parse_arg('rmin', model%rmin) == PARSE_TYPE_ERROR) stop
if (parse_arg('rmax', model%rmax) == PARSE_TYPE_ERROR) stop

if (parse_arg('rmax', model%rmax) == PARSE_TYPE_ERROR) stop
if (parse_arg('rmax', model%rmax) == PARSE_TYPE_ERROR) stop

if (parse_arg('in', infile) == PARSE_TYPE_ERROR) stop
if (parse_arg('out', outfile) == PARSE_TYPE_ERROR) stop


if (model%mu   < 0._pr) model%mu   = 3.0e10_pr
if (model%na   < 0)     model%na   = 5000
if (model%rmin < 0._pr) model%rmin = 2.0_pr
if (model%rmax < 0._pr) model%rmax = 0.35_pr
if (model%rmin < 2.0_pr)  write(error_unit,*) 'WARNING: rmin < 2, sub-events may not be well represented'
if (model%rmax >= 0.5_pr) then
    write(error_unit,*) 'ERROR: rmax must be < 0.5'
    stop
endif


end subroutine parse_args
!==============================================================================