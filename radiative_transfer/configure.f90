module configure
  implicit none

  ! Integer parameters
  integer, parameter :: DB_NTEMP_MAX      = 3000
  integer, parameter :: FRSIZE_NDIFF      = 400
  integer, parameter :: FRSIZE_X          = 200
  integer, parameter :: FRSIZE_Y          = 300
  integer, parameter :: FRSIZE_Y_SMALL     = 300
  integer, parameter :: FRSIZE_FREQ       = 130
  integer, parameter :: FRSIZE_MU         = 36
  integer, parameter :: FRSIZE_MU_HALF = 300
  integer, parameter :: FRSIZE_PHI        = 16
  integer, parameter :: FRSIZE_CHAR       = 24
  integer, parameter :: FRSIZE_DRR        = 11
  integer, parameter :: FRSIZE_INF_PHI    = 40
  integer, parameter :: FRSIZE_INF_R      = (2*FRSIZE_X + 100)
  integer, parameter :: FRSIZE_INF_X      = 200
  integer, parameter :: FRSIZE_INF_Y      = 200
  integer, parameter :: DUST_SPECIES_MAX  = 2
  integer, parameter :: DUST_SIZE_MAX     = 1
  integer, parameter :: DUST_TRANGE_MAX   = 1
  integer, parameter :: FRSIZE_SRC_BK     = 4
  integer, parameter :: MAXRAYNR          = 256000
  integer, parameter :: FRSIZE_CAM_MU     = 100
  integer, parameter :: FRSIZE_CAM_PHI    = 128
  integer, parameter :: FRSIZE_MAX = 5000

  ! Floating-point parameters
  double precision, parameter :: SAFETY_CHECK_TAUMAX = 1.d5
  double precision, parameter :: SHORTCHAR_EPS       = 1.d-10
  double precision, parameter :: TELESC_EPS          = 1.d-10

  ! Compile-time feature flags (logical parameters)
  logical, parameter :: RADGRID_TWODIM         = .true.
  logical, parameter :: INCLUDE_DUST           = .true.
  logical, parameter :: COORD_SPHERICAL        = .true.
  logical, parameter :: ESC_MINIMAL_EXTENSION  = .true.
  logical, parameter :: DUST_OPAC_TEMPDEP      = .true.
  logical, parameter :: SMALL_MEMORY           = .true.
  logical, parameter :: CHECK_NUMBERS          = .true.
  logical, parameter :: MIRROR_THETA           = .true.
  logical, parameter :: MIRROR_PHI             = .true.
  logical, parameter :: INTERPOL_PHI_3         = .true.
  logical, parameter :: INTERPOL_THETA_3       = .true.
  logical, parameter :: INTERPOL_SRC_1         = .true.
  logical, parameter :: INTERPOL_FREQ_1        = .true.
  logical, parameter :: SRCQDR_PNT_2           = .true.
  logical, parameter :: ALI_PNT_2              = .true.
  logical, parameter :: SAFETY_CHECKS_ACTIVE   = .true.
  logical, parameter :: SHORTCHAR_NEWSTYLE     = .true.
  logical, parameter :: CENTRAL_SOURCE         = .true.
  logical, parameter :: DEBUG_FILLINTENS       = .true.
  logical, parameter :: LINUXVERSION           = .true.
  logical, parameter :: INCLUDE_EMISQUANT      = .true.
  logical, parameter :: INCLUDE_QUANTSOURCENEW = .true.

end module configure

