module common_montecarloarrays
use configure
  implicit none

  ! Constants
  !integer, parameter :: FRSIZE_FREQ = 1000         ! Replace with actual value
  !integer, parameter :: FRSIZE_Y_SMALL = 100       ! Replace with actual value
  !integer, parameter :: FRSIZE_X = 500             ! Replace with actual value
  !integer, parameter :: DUST_SIZE_MAX = 20         ! Replace with actual value
  !integer, parameter :: DUST_SPECIES_MAX = 10       ! Replace with actual value

  ! Monte Carlo radiative transfer arrays
  double precision :: vol(FRSIZE_Y_SMALL, FRSIZE_X)
  double precision :: rrigrid(FRSIZE_X+1)
  double precision :: ttigrid(FRSIZE_Y_SMALL+1)
  double precision :: alpha_a(FRSIZE_FREQ, FRSIZE_Y_SMALL, FRSIZE_X)
  double precision :: alpha_s(FRSIZE_FREQ, FRSIZE_Y_SMALL, FRSIZE_X)
  double precision :: kappa_a(FRSIZE_FREQ, DUST_SIZE_MAX, DUST_SPECIES_MAX)
  double precision :: kappa_s(FRSIZE_FREQ, DUST_SIZE_MAX, DUST_SPECIES_MAX)
  double precision :: cumulener(DUST_SIZE_MAX, DUST_SPECIES_MAX, FRSIZE_Y_SMALL, FRSIZE_X)
  double precision :: cumulener_bk(DUST_SIZE_MAX, DUST_SPECIES_MAX, FRSIZE_Y_SMALL, FRSIZE_X)
  double precision :: freqdistr(FRSIZE_FREQ, FRSIZE_Y_SMALL, FRSIZE_X)

  integer, parameter :: INPUT_MAX_LEN = 5000   ! or your existing limit
  character(len=160), dimension(INPUT_MAX_LEN) :: inpstring
  integer :: nrstring
  integer, dimension(INPUT_MAX_LEN) :: lineok

end module common_montecarloarrays

