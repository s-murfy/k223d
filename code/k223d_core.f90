module k223d_core
!==============================================================================
! Shared computational core for k223d.
! Called identically by the standalone executable (k223d.f90) and the
! C API wrapper (k223d_capi.f90). Contains no file I/O.
!
! Entry point: compute_source
!==============================================================================

use generic
use LAT_mesh
use LAT_source
use mesh_geom
use model_setup
use time
use source_model
implicit none
 
contains

!==============================================================================
subroutine compute_source(amesh, surface, model, pdf, src)
! Runs the full k223d source computation pipeline:
!   1. Build mesh geometry (connectivity, distances, boundary distances)
!   4. Compute slip distribution
!   5. Compute rupture time if velocity and initial time are both present
!
! src%velocity and src%rupt_time must be allocated by the caller if rupture
! time computation is required (see read_vtk_mesh or k223d_capi).
! src%nuc_id is derived internally from the minimum of src%rupt_time (value 0.0).
!==============================================================================
  type(mesh),        intent(in)    :: amesh
  type(reflect),     intent(in)    :: surface
  type(model_param), intent(inout) :: model   ! inout: defaults filled if unset
  type(pdfinputs),   intent(inout) :: pdf     ! inout: pdf_type may be set
  type(source),      intent(inout) :: src     ! inout: velocity/rupt_time may be pre-loaded
 ! internal variables 
  type(mesh_geometry)              :: geom
 
  ! --- Step 1: mesh geometry ---
  write(*,*) 'build_mesh_geometry'
   call build_mesh_geometry(amesh, surface, geom)
 
  ! --- Step 2 : Area-weighted uniform PDF 
    if (pdf%is_uniform) then
        allocate(pdf%distrib(amesh%Ncells))
        pdf%distrib = geom%area / sum(geom%area)
    endif

  ! --- Step 3: slip distribution ---
  write(*,*) 'construct slip distribution'
  call pdftoslip(model, amesh, geom, pdf,src)
 
  ! --- Step 3: rupture time (only if velocity and initial time were provided) ---
  if (allocated(src%rupt_time) .and. allocated(src%velocity)) then
    write(*,*) 'calculate rupture time'
    src%nuc_id = minloc(src%rupt_time, 1)  ! node with value 0.0
    call calc_rupt_front(amesh, geom, src)
  endif
 
end subroutine compute_source 
!==============================================================================







end module k223d_core 