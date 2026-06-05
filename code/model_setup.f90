module model_setup
use LAT_source
implicit none

contains 



!==============================================================================
subroutine set_model_defaults(model)
  type(model_param), intent(inout) :: model   ! inout: defaults filled if unset

if (model%mu   < 0._pr) model%mu  = 3.0e10_pr   ! fixed constant,
if (model%na   < 0)     model%na  = 5000          ! fixed constant
if (model%rmin < 0._pr) model%rmin = 2.0_pr       ! needs geom%area
if (model%rmax < 0._pr) model%rmax = 0.35_pr          ! needs mesh depth + mean dip

! mu=3.0e10, 
! rmin=2.,
! rmax=0.35, 
! na=5000,

end subroutine set_model_defaults
!==============================================================================
subroutine build_default_pdf(pdf,amesh)
implicit none
type(mesh),        intent(in) :: amesh
type(pdfinputs),   intent(inout) :: pdf  

allocate(pdf%distrib(amesh%Ncells))

pdf%distrib = 1._pr/float(amesh%Ncells-1)

end subroutine build_default_pdf
!==============================================================================
end module model_setup