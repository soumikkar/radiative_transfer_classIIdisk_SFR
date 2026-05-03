module common_diffusion
use configure
  implicit none

  ! Constants
  !integer, parameter :: FRSIZE_NDIFF = 500   ! <- replace with actual value
 ! integer, parameter :: FRSIZE_Y_SMALL = 100 ! <- replace with actual value
 ! integer, parameter :: FRSIZE_X = 500       ! <- replace with actual value

  ! Diffusion photon and index variables
  integer :: nphotdiff, ndiff
  integer :: irdiff(FRSIZE_NDIFF)
  integer :: itdiff(FRSIZE_NDIFF)
  integer :: idiff_rleft(FRSIZE_NDIFF)
  integer :: idiff_rright(FRSIZE_NDIFF)
  integer :: idiff_tleft(FRSIZE_NDIFF)
  integer :: idiff_tright(FRSIZE_NDIFF)
  integer :: ipde(FRSIZE_NDIFF)
  integer :: idifcell(FRSIZE_Y_SMALL, FRSIZE_X)

  ! Diffusion matrix and solver arrays
  double precision :: dmatrix(FRSIZE_NDIFF, FRSIZE_NDIFF)
  double precision :: dmat(FRSIZE_NDIFF, FRSIZE_NDIFF)
  double precision :: drhs(FRSIZE_NDIFF)
  double precision :: dsol(FRSIZE_NDIFF)
  double precision :: dsolold(FRSIZE_NDIFF)
  double precision :: diffconst(FRSIZE_NDIFF)

end module common_diffusion

