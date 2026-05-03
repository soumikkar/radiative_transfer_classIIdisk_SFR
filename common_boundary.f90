module common_boundary
use configure
  implicit none

  ! Constants
  !integer, parameter :: FRSIZE_FREQ = 1000  ! <-- set to your actual value

  ! Variables specifying the boundaries
  double precision :: radbnd_rstar
  double precision, dimension(FRSIZE_FREQ) :: radbnd_starspec
  double precision :: radbnd_startemp
  double precision :: radbnd_lstar
  double precision, dimension(FRSIZE_FREQ) :: radbnd_extspec
  double precision :: radbnd_starlumtot
  double precision :: radbnd_extlumtot
  double precision :: radbnd_mstar

  ! Integer control flags for boundaries
  integer :: iradbnd_in_itype
  integer :: iradbnd_out_itype
  integer :: iradbnd_in_spectype
  integer :: iradbnd_read_rstar

  ! Equator disk parameters
  double precision :: eqd_mmstar
  double precision :: eqd_mmdot
  double precision :: eqd_rinn
  double precision :: eqd_rout
  double precision :: eqd_xi
  double precision :: eqd_fracnonevap

  ! Integer flags for equator disk options
  integer :: ieqd_active
  integer :: ieqd_simple
  integer :: ieqd_type

end module common_boundary

