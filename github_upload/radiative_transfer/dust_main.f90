module dust_main
  use configure
  use common_grid       
  use common_dust       
  use numerical_receipe
  implicit none

contains

  subroutine read_dustdata()
    implicit none

    integer :: ispec, idum, idustfile, iformat, idum2
    real(8) :: temp0, temp1, dtemp0, dtemp1
    character(len=80) :: comstring
    character(len=256) :: filename, base, ext
    integer :: ios

    ! Check if transfer routine provided expected number of dust species
   ! if (dust_setup_nrspecies == 0) then
    !  write(*,*) 'ERROR: Number of dust species unspecified. Please fix.'
    !  error stop 99
   ! end if

    ! Open master dust file
    open(unit=3, file='dustopac.inp', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Could not open file dustopac.inp'
      error stop 13
    end if

    read(3,*) iformat
    read(3,*) dust_nr_species
    read(3,*) comstring
    print*, iformat, dust_nr_species

    if (dust_nr_species > DUST_SPECIES_MAX) then
      write(*,*) 'ERROR: Too many dust species in dustopac.inp (limit = ', DUST_SPECIES_MAX, ')'
      error stop
    end if

    ! Loop over dust species
    do ispec = 1, dust_nr_species
      read(3,*) idum
      read(3,*) idum2  
      read(3,*) idustfile

      base = 'dustopac_'
      ext  = '.inp'
      call make_indexed_filename(base, idustfile, ext, filename)
      PRINT*, idum, idustfile, dust_nr_species, filename

      select case (idum)
      case (-1)
        ! Temperature-independent opacities
        call read_dustopac_file(ispec, filename)
        print*, 'soumik'
      case (-2)
        ! Opacities with tmin and tmax
        read(3,*) temp0
        read(3,*) temp1
        dust_tmin(ispec)  = temp0
        dust_tmax(ispec)  = temp1
        dust_dtmin(ispec) = 0.d0
        dust_dtmax(ispec) = 0.d0
        call read_dustopac_file(ispec, filename)

      case (-3)
        ! Opacities with tmin, tmax and smooth cutoffs
        read(3,*) temp0
        read(3,*) dtemp0
        read(3,*) temp1
        read(3,*) dtemp1
        dust_tmin(ispec)  = temp0
        dust_tmax(ispec)  = temp1
        dust_dtmin(ispec) = abs(dtemp0)
        dust_dtmax(ispec) = abs(dtemp1)
        call read_dustopac_file(ispec, filename)

      case default
        write(*,*) 'ERROR: Unsupported idum=', idum, ' in dustopac.inp'
        error stop 13
      end select

      ! Read comment string after each species block
      read(3,*) comstring
    end do

    close(3)

    ! Check if number of species matches expected
    if (dust_setup_nrspecies > 0 .and. dust_setup_nrspecies /= dust_nr_species) then
      write(*,*) 'ERROR: Mismatch in dust species count. Expected ', dust_setup_nrspecies, ' found ', dust_nr_species
      error stop 99
    else if (dust_setup_nrspecies < 0) then
      dust_setup_nrspecies = dust_nr_species
    end if

#ifndef EXCLUDE_DUST_ITERATION
    if (dust_frwgt_read == 0) then
      call make_dust_freq_weights()
      dust_frwgt_read = 1
    end if
#endif

    dust_done_read = 1

  end subroutine read_dustdata




subroutine read_dustopac_file(ispec, filename)
   ! use common_dust
  !  use common_grid
    implicit none

    character(len=80), intent(in) :: filename
    integer, intent(in)          :: ispec
    integer                     :: ifr, isize, nsize, ifreq
    double precision            :: dummy
    integer                     :: ifile

    ! Open opacity file
    ifile = 10
    open(unit=ifile, file=filename, status='old', action='read')
   ! if (ifr /= 0) then
    !    write(*,*) 'ERROR: Could not open opacity file: '
    !    stop 31
    !end if

    ! Read header: number of frequencies, number of sizes
    read(ifile, *) ifreq, nsize
    if (ifreq /= freq_nr) then
        write(*,*) 'ERROR: Frequency count mismatch: file=', ifreq, ' expected=', freq_nr
        stop
    end if

    ! Check for size and species limits
    if (nsize > DUST_SIZE_MAX) then
        write(*,*) 'ERROR: nsize exceeds DUST_SIZE_MAX.'
        stop
    end if
    if (ispec > DUST_SPECIES_MAX) then
        write(*,*) 'ERROR: ispec exceeds DUST_SPECIES_MAX.'
        stop
    end if

    ! Store size count for this species
    dust_nr_size(ispec) = nsize
    dust_nr_temp(ispec) = 1  ! No temp dependence

    ! Set temperature range defaults
    dust_temprange_low(1, ispec)  = 0.d0
    dust_temprange_high(1, ispec) = 1.d33

    ! Read absorption opacities
    do ifr = 1, freq_nr
        do isize = 1, nsize
            read(ifile, *) dummy
            if (dummy <= 0.d0) then
                write(*,*) 'ERROR: Absorption opacity invalid: ', dummy
                stop
            end if
            dust_kappawgt_abs(ifr, 1, isize, ispec) = dummy
        end do
    end do

    ! Read scattering opacities
    do ifr = 1, freq_nr
        do isize = 1, nsize
            read(ifile, *) dummy
            if (dummy < 0.d0) then
                write(*,*) 'ERROR: Scattering opacity invalid: ', dummy
                stop
            end if
            dust_kappawgt_scat(ifr, 1, isize, ispec) = dummy
        end do
    end do

    ! Done
    close(ifile)
    return

end subroutine read_dustopac_file




function find_dust_kappa(inu, isize, ispec, temp, iabs, iscat) result(kappa)
    implicit none
    integer, intent(in) :: inu, isize, ispec, iabs, iscat
    real(8), intent(in) :: temp
    real(8) :: kappa

    real(8) :: tmin, tmax, dtmin, dtmax, plindex
    real(8) :: condense, fact, omfact
    integer :: ntemp, itlo, ithi

    condense = 1.0d0
    ntemp    = dust_nr_temp(ispec)
    tmin     = dust_temprange_low(1, ispec)
    tmax     = dust_temprange_high(ntemp, ispec)
    dtmin    = dust_dtmin(ispec)
    dtmax    = dust_dtmax(ispec)
    itlo     = 1
    ithi     = 1
    fact     = 1.0d0

    ! Handle temperature regime limits
    if (temp < tmin) then
      if (dtmin == 0.d0) then
        condense = 0.d0
      else
        plindex  = log(1.d0 - (dtmin / tmin))
        condense = (tmin / temp)**plindex
        if (condense > 1.d0) stop 'Error: condense factor > 1 at low T'
      end if
    else if (temp >= tmax) then
      if (dtmax == 0.d0 .and. tmax > 0.d0) then
        condense = 0.d0
      else
        plindex  = log(1.d0 + (dtmax / tmax))
        condense = (tmin / temp)**plindex
        if (condense > 1.d0) stop 'Error: condense factor > 1 at high T'
      end if
    else
      if (ntemp > 1) then
        call hunt(dust_temprange_low(:, ispec), ntemp, temp, itlo)
        if (temp > dust_temprange_high(itlo, ispec)) then
          if (itlo >= ntemp) then
            stop 'Inconsistency in dust opacity temperature range'
          end if
          ithi = itlo + 1
          fact = (temp - dust_temprange_high(itlo, ispec)) / &
                 (dust_temprange_low(ithi, ispec) - dust_temprange_high(itlo, ispec))
        else
          ithi = itlo
          fact = 1.0d0
        end if
      end if
    end if

    omfact = 1.0d0 - fact

    ! Interpolate opacities
    kappa = 0.d0
    if (iabs /= 0) then
      kappa = kappa + omfact * dust_kappawgt_abs(inu, itlo, isize, ispec) + &
                       fact   * dust_kappawgt_abs(inu, ithi, isize, ispec)
    end if
    if (iscat /= 0) then
      kappa = kappa + omfact * dust_kappawgt_scat(inu, itlo, isize, ispec) + &
                       fact   * dust_kappawgt_scat(inu, ithi, isize, ispec)
    end if

    ! Apply condensation factor
    if (iabs >= 0 .and. iscat >= 0) then
      kappa = kappa * condense
    end if

  end function find_dust_kappa




subroutine make_dust_freq_weights()
  use common_grid
  use common_dust
  implicit none

  integer :: ifreq

  if (freq_nr < 1) then
    write(*,*) 'ERROR: Zero or negative number of frequencies is invalid.'
    stop 99

  elseif (freq_nr == 1) then
    dust_freq_wgt(1) = 1.0d0
    return

  elseif (freq_nr >= 2 .and. freq_nr <= 5) then
    write(*,*) 'WARNING (dust temperature solver):'
    write(*,*) '  Too few frequency bins for consistent dust temperatures.'
    write(*,*) '  Continuing, but results may be inaccurate.'
    dust_warn_few_freqs = 1
  else
    dust_freq_wgt(1)        = 0.5d0 * abs(freq_nu(2) - freq_nu(1))
    dust_freq_wgt(freq_nr)  = 0.5d0 * abs(freq_nu(freq_nr) - freq_nu(freq_nr-1))

    do ifreq = 2, freq_nr-1
      dust_freq_wgt(ifreq) = 0.5d0 * abs(freq_nu(ifreq+1) - freq_nu(ifreq-1))
    end do

    dust_frwgt_read = 1
  endif

end subroutine make_dust_freq_weights


end module dust_main

