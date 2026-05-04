module common_dust
use configure
  implicit none

  ! Parameters
  !integer, parameter :: FRSIZE_FREQ=1000, FRSIZE_Y_SMALL=100, FRSIZE_X=500
  !integer, parameter :: DUST_TRANGE_MAX=50, DUST_SIZE_MAX=20, DUST_SPECIES_MAX=10

  ! Dust opacity status
  integer :: dust_done_read, dust_opacity_tempdep

  ! Dust opacities and temperature ranges
  double precision :: dust_kappawgt_abs(FRSIZE_FREQ, DUST_TRANGE_MAX, DUST_SIZE_MAX, DUST_SPECIES_MAX)
  double precision :: dust_kappawgt_scat(FRSIZE_FREQ, DUST_TRANGE_MAX, DUST_SIZE_MAX, DUST_SPECIES_MAX)
  double precision :: dust_freq_wgt(FRSIZE_FREQ)
  double precision :: dust_temprange_low(DUST_TRANGE_MAX, DUST_SPECIES_MAX)
  double precision :: dust_temprange_high(DUST_TRANGE_MAX, DUST_SPECIES_MAX)
  double precision :: dust_tmin(DUST_SPECIES_MAX), dust_tmax(DUST_SPECIES_MAX)
  double precision :: dust_dtmin(DUST_SPECIES_MAX), dust_dtmax(DUST_SPECIES_MAX)
  double precision :: dust_n_catom(DUST_SPECIES_MAX)

  ! Dust density and temperatures
  double precision :: dust_rho(DUST_SPECIES_MAX, FRSIZE_Y_SMALL, FRSIZE_X)
  double precision :: dust_temp(DUST_SIZE_MAX, DUST_SPECIES_MAX, FRSIZE_Y_SMALL, FRSIZE_X)

  ! Dust configuration
  integer :: dust_nr_species, dust_warn_zero_temp
  integer :: dust_frwgt_read, dust_warn_few_freqs
  integer :: dust_nr_size(DUST_SPECIES_MAX)
  integer :: dust_nr_temp(DUST_SPECIES_MAX)
  integer :: dust_quantum(DUST_SPECIES_MAX)

  ! Dust setup
  integer :: dust_setup_nrspecies
  integer :: dust_setup_nrsizes(DUST_SPECIES_MAX)

end module common_dust

