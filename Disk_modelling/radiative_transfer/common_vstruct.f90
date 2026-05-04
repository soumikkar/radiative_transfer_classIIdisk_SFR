module common_vstruct
use configure
  implicit none

  ! Parameters
 ! integer, parameter :: FRSIZE_X=500, DUST_SPECIES_MAX=10

  ! Vertical structure flags
  integer :: ivstr_dostruct(DUST_SPECIES_MAX)

  ! Vertical structure quantities
  double precision :: vstr_sigmadust(FRSIZE_X, DUST_SPECIES_MAX)
  double precision :: vstr_massdust(DUST_SPECIES_MAX)
  double precision :: vstr_massdusttot

end module common_vstruct

