module common_montecarlo
use configure
  implicit none
  !integer, parameter :: FRSIZE_FREQ = 1000, FRSIZE_Y = 500
  integer :: mc_integerspec(FRSIZE_FREQ, FRSIZE_Y)
  integer :: itempdecoup, iquantum
  integer :: istarsurf
  !integer, parameter :: FRSIZE_Y_SMALL = 100, FRSIZE_X = 500
  integer :: iphotcurr, nrskips
  integer :: ilastphot(FRSIZE_Y_SMALL, FRSIZE_X)
  integer :: iphotcount(FRSIZE_Y_SMALL, FRSIZE_X)
  !integer, parameter :: DB_NTEMP_MAX = 200
  !integer, parameter :: DUST_SIZE_MAX = 20, DUST_SPECIES_MAX = 10

  integer :: db_ntemp
  double precision :: db_emiss(FRSIZE_FREQ, DB_NTEMP_MAX, DUST_SIZE_MAX, DUST_SPECIES_MAX)
  double precision :: db_cumulnorm(FRSIZE_FREQ+1, DB_NTEMP_MAX, DUST_SIZE_MAX, DUST_SPECIES_MAX)
  double precision :: db_enertemp(DB_NTEMP_MAX, DUST_SIZE_MAX, DUST_SPECIES_MAX)
  double precision :: db_temp(DB_NTEMP_MAX)


end module common_montecarlo

