module montecarlo
  use configure
  use common_grid
  use common_dust
  use common_boundary
  use common_montecarlo
  use common_montecarloarrays
  use numerical_receipe
  use dust_main
  !use diffusion_array
  !use mod_vertstruct_globals
  implicit none
 

contains

subroutine do_monte_carlo(nphot, spectrum, tdust, scatsrc, &
                          meanint, imethod, ifast, enthres, iseed, &
                          cntdump, irestart, ntemp, temp0, temp1)

  ! Import kinds
  !use iso_fortran_env, only: dp => real64
  implicit none

  ! Input variables
  integer, intent(in)               :: nphot, imethod, ifast, cntdump
  integer, intent(inout)            :: iseed
  integer, intent(in)               :: irestart, ntemp
  real(8), intent(in)              :: enthres, temp0, temp1

  ! Array inputs/outputs (assumed-shape)
  real(8) :: spectrum(FRSIZE_FREQ,FRSIZE_Y_SMALL), tdust(DUST_SIZE_MAX,DUST_SPECIES_MAX,FRSIZE_Y_SMALL,FRSIZE_X)        
  real(8) ::scatsrc(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
  
  real(8), optional :: meanint(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)

  ! Local variables
  real(8) :: energy, theta, mirfact, tempav, rstar
  integer  :: nr, nt, iphot, count, cnt, inu, it, ir, ispec, isize, imirt

  ! Constants
  real(8), parameter :: pi = 3.14159265358979323846d0
  real(8), parameter :: twopi = 2.0d0 * pi
  real(8), parameter :: fourpi = 4.0d0 * pi
  real(8), parameter :: parsec = 3.08572d18

  ! Arrays
  real(8):: starlum(FRSIZE_FREQ), extlum(FRSIZE_FREQ)
  integer  :: integerspec(FRSIZE_FREQ,FRSIZE_Y_SMALL)
  

  ! Initialize
  count = 1; cnt = 1
  

  ! Set star radius
  if (istarsurf == 1) then
    rstar = radbnd_rstar
  else
    rstar = 0.0d0
  end if

  ! Number of grid cells
  nr = irsi_frsizex
  nt    = (irsi_frsizey+1)/2
  imirt = 1
  mirfact = 0.5
!print*, "nt=", nt


!print*,'imirt=', imirt
  ! Prepare simulation
  call prep_monte_carlo(nr, nt, rrigrid, ttigrid, vol, &
                        alpha_a, alpha_s, kappa_a, kappa_s, &
                        tdust, scatsrc, meanint, &
                              cumulener, cumulener_bk, freqdistr, &
                              integerspec, ilastphot, iphotcount, iphotcurr)

  ! Make emissivity database
  call make_temperature_grid(ntemp, temp0, temp1)
  call make_emiss_dbase()


  ! Compute stellar luminosity spectrum
  do inu = 1, freq_nr
    starlum(inu) = pi * radbnd_starspec(inu) * fourpi * radbnd_rstar**2
    !print*, radbnd_starspec(inu)
  end do

 radbnd_starlumtot = 0.d0
      do inu=1,freq_nr
          radbnd_starlumtot = radbnd_starlumtot + &
               starlum(inu) * dust_freq_wgt(inu)
      enddo

  ! Energy per photon packet
  energy = (mirfact * radbnd_starlumtot) / nphot

  ! Start photon emission loop
  do iphot = 1, nphot

     iphotcurr = iphot

    if (mod(iphot, 1000) == 0) print *, 'Photon number: ', iphot
    !print*, "walk full path starts"
    ! Perform photon walk
    call walk_full_path(nr, nt, rrigrid, ttigrid, vol, &
                        alpha_a, alpha_s, kappa_a, rstar, &
                        tdust, scatsrc, meanint, &
                        cumulener, cumulener_bk, &
                        freqdistr,energy, inu, theta, imethod, ifast, enthres, iseed)
   ! print*, "walk full path ends"
    if (theta >= 0.0d0) then
      call hunt(ttigrid, nt, theta, it)
      if(it.gt.nt) it=nt
      if(it.lt.1) it=1
     !print*, it
      integerspec(inu, it) = integerspec(inu, it) + 1
      !print*, "integerspec=",integerspec(inu,it)
    end if
     !print*, "integerspec1=",integerspec(inu,it)
    count = count + 1
    cnt = cnt + 1

  end do

  ! Compute spectrum
  do it = 1, nt
    do inu = 1, freq_nr
      spectrum(inu, it) = integerspec(inu, it) * energy / &
           (twopi * parsec**2 * dust_freq_wgt(inu) * &
           abs(cos(ttigrid(it+1)) - cos(ttigrid(it))))
      !print*, "spectrum calculation=", spectrum(inu,it)
    end do
  end do
   
 
 ! Compute temperatures
    if (itempdecoup == 1) then
        do ir = 1, nr
            do it = 1, nt
                do ispec = 1, dust_nr_species
                    do isize = 1, dust_nr_size(ispec)
                        energy = cumulener(isize, ispec, it, ir) / &
                                 (dust_rho(ispec, it, ir) * vol(it, ir))
                        tdust(isize, ispec, it, ir) = compute_dusttemp_energy_bd(energy, ispec, isize)
                        !print*, "tdust", tdust(isize, ispec, it, ir)
                    end do
                end do
            end do
        end do
    else
        do ir = 1, nr
            do it = 1, nt
                tempav = compute_dusttemp_coupled_bd(cumulener(1,1,it,ir), dust_rho(1,it,ir), vol(it,ir))
                do ispec = 1, dust_nr_species
                    do isize = 1, dust_nr_size(ispec)
                        tdust(isize, ispec, it, ir) = tempav
                    end do
                end do
            end do
        end do
    end if


  print *, 'Monte Carlo simulation completed.'

end subroutine do_monte_carlo




subroutine prep_monte_carlo(nr, nt, rrigrid, ttigrid, vol, &
                              alpha_a, alpha_s, kappa_a, kappa_s, &
                              tdust, scatsrc, meanint, &
                              cumulener, cumulener_bk, freqdistr, &
                              integerspec, ilastphot, iphotcount, iphotcurr)
    
    implicit none

    ! Arguments
    integer, intent(out) :: nr, nt,iphotcurr
    real(8), intent(out) :: rrigrid(FRSIZE_X+1), ttigrid(FRSIZE_Y_SMALL+1)
    real(8), intent(out) :: vol(FRSIZE_Y_SMALL,FRSIZE_X)
    real(8), intent(out) :: alpha_a(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
    real(8):: alpha_s(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
    real(8), intent(out) :: kappa_a(FRSIZE_FREQ,DUST_SIZE_MAX,DUST_SPECIES_MAX), kappa_s(FRSIZE_FREQ,DUST_SIZE_MAX,DUST_SPECIES_MAX)
    real(8), intent(out) :: tdust(DUST_SIZE_MAX,DUST_SPECIES_MAX,FRSIZE_Y_SMALL,FRSIZE_X)
    real(8), intent(out) :: scatsrc(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
    real(8), intent(out) :: cumulener(DUST_SIZE_MAX,DUST_SPECIES_MAX,FRSIZE_Y_SMALL,FRSIZE_X)
    real(8):: cumulener_bk(DUST_SIZE_MAX,DUST_SPECIES_MAX,FRSIZE_Y_SMALL,FRSIZE_X)
    real(8), intent(out) :: freqdistr(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
    integer, intent(out) :: integerspec(FRSIZE_FREQ,FRSIZE_Y_SMALL), ilastphot(FRSIZE_Y_SMALL,FRSIZE_X)
    integer :: iphotcount(FRSIZE_Y_SMALL,FRSIZE_X)
!#ifdef SAVE_MEANINT
    real(8), intent(out) :: meanint(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
!#else
 !   real(8), intent(out) :: meanint
!#endif

    ! Local variables
    integer :: ir, it, isize, ispec, inu, imirt
    real(8) :: temp
    real(8), parameter :: pi = 3.14159265358979323846264338328d0
    real(8), parameter :: pihalf=1.57079632679489661923132169164d0
    logical:: dust_opacity_tempdep

    ! Sanity checks
    if (dust_opacity_tempdep) then
      print *, 'ERROR: Temperature-dependent opacities not supported.'
      stop 13
    end if

    ! Compute opacities (constant temperature assumption)
    temp = 100.0d0
    do inu = 1, freq_nr
      do ispec = 1, dust_nr_species
        do isize = 1, dust_nr_size(ispec)
          kappa_a(inu,isize,ispec) = find_dust_kappa(inu,isize,ispec,temp,1,0)
          kappa_s(inu,isize,ispec) = find_dust_kappa(inu,isize,ispec,temp,0,1)
        !print*, kappa_a(inu,isize,ispec),  kappa_s(inu,isize,ispec)
        end do
      end do
    end do


    ! Copy grid structure
    nr = irsi_frsizex
    nt = (irsi_frsizey + 1) / 2
    imirt = 1
   ! print*, "nt1=", nt

  !  print*, "imirt==", imirt
    do ir = 1, nr+1
      rrigrid(ir) = rsi_x_i(ir,1)
    end do
    do it = 1, nt+1
      ttigrid(it) = rsi_x_i(it,2)
    end do

!if(imirt.eq.1) then
          if(abs(ttigrid(nt+1)-pihalf).gt.1d-3) then
              write(*,*) 'Monte Carlo: ERROR in interface grid!!'
              stop 13
          endif
      !endif
      !if(imirt.eq.0) then
      !    if(abs(ttigrid(nt/2+1)/pihalf-1.d0).gt.1d-6) stop 91477
      !    ttigrid((nt+1)/2)=pihalf
      !else
          if(abs(ttigrid(nt+1)/pihalf-1.d0).gt.1d-6) then
              write(*,*) ttigrid(nt),pihalf
              stop 91478
          endif
          ttigrid(nt+1)=pihalf
     ! endif
      if(abs(ttigrid(1)).gt.1d-6) stop 91477
      ttigrid(1)=0.d0

    ! Compute cell volumes
    do ir = 1, nr
      do it = 1, nt
        vol(it,ir) = (1.0d0/3.0d0) * 2.0d0 * pi * (rrigrid(ir+1)**3 - rrigrid(ir)**3) &
                     * abs(cos(ttigrid(it)) - cos(ttigrid(it+1)))
      end do
    end do

    ! Compute absorption and scattering coefficients
    do ir = 1, nr
      do it = 1, nt
        do inu = 1, freq_nr
          alpha_a(inu,it,ir) = 0.0d0
          alpha_s(inu,it,ir) = 0.0d0
          do ispec = 1, dust_nr_species
            do isize = 1, dust_nr_size(ispec)
              alpha_a(inu,it,ir) = alpha_a(inu,it,ir) + &
                                   dust_rho(ispec,it,ir) * kappa_a(inu,isize,ispec)
              alpha_s(inu,it,ir) = alpha_s(inu,it,ir) + &
                                   dust_rho(ispec,it,ir) * kappa_s(inu,isize,ispec)
            !print*, alpha_a(inu,it,ir), alpha_s(inu,it,ir)
            end do
          end do
        end do
      end do
    end do

    ! Zero source terms and radiation fields
    do ir = 1, nr
      do it = 1, nt
        do inu = 1, freq_nr
          scatsrc(inu,it,ir) = 0.0d0
          freqdistr(inu,it,ir) = 0.0d0
          meanint(inu,it,ir) = 0.0d0

        end do
      end do
    end do


      do ir=1,nr
          do it=1,nt
              do ispec=1,dust_nr_species
                  do isize=1,dust_nr_size(ispec)
                      tdust(isize,ispec,it,ir)        = 0.d0
                      cumulener(isize,ispec,it,ir)    = 0.d0
                      cumulener_bk(isize,ispec,it,ir) = 0.d0
                  enddo
              enddo
          enddo
      enddo

!    Reset the integer spectrum

      do it=1,nt
          do inu=1,freq_nr
              integerspec(inu,it) = 0
          enddo
      enddo

!     Reset the photon counters

      do ir=1,nr
          do it=1,nt
              ilastphot(it,ir) = 0
              iphotcount(it,ir) = 0
          enddo
      enddo
      iphotcurr = 0

  end subroutine prep_monte_carlo




subroutine randomfreq(flux, inupick, iseed)

  implicit none

  ! Arguments
  real(8), intent(in)  :: flux(FRSIZE_FREQ)
  integer,  intent(out) :: inupick
  integer,  intent(inout) :: iseed

  ! Local variables
  real(8) :: cumul(FRSIZE_FREQ+1)
  real(8) :: rn
  integer  :: inu

  ! Check frequency weight initialization
  if (dust_freq_wgt(1) == 0.0d0) then
     stop 'ERROR 8331: dust_freq_wgt not initialized.'
  end if

  ! Build cumulative probability distribution function
  cumul(1) = 0.0d0
  do inu = 1, freq_nr
     cumul(inu+1) = cumul(inu) + flux(inu) * dust_freq_wgt(inu)
  end do

  ! Normalize cumulative distribution to 1.0
   do inu=1,freq_nr+1
          cumul(inu) = cumul(inu) / cumul(freq_nr+1)
      enddo
  ! Draw random number in [0,1)
  rn = ran2(iseed)

  ! Locate corresponding frequency bin in cumulative table
  call hunt(cumul, freq_nr, rn, inupick)

end subroutine randomfreq




subroutine find_position_in_cell(ttigrid, rrigrid, it, ir, r, theta, iseed)
    ! Finds a random (r, theta) position inside the selected grid cell.
    implicit none
    ! Arguments
    real(8), intent(in)    :: ttigrid(FRSIZE_Y_SMALL+1)
    real(8), intent(in)    :: rrigrid(FRSIZE_X + 1)
    integer,  intent(in)    :: it, ir
    real(8), intent(out)   :: r, theta
    integer,  intent(inout) :: iseed

    ! Locals
    real(8) :: rn

    ! Theta position in the selected theta cell
    rn    = ran2(iseed)
    theta = ttigrid(it) + rn * (ttigrid(it+1) - ttigrid(it))

    ! Radial position in the selected radial cell
    rn = ran2(iseed)
    r  = rrigrid(ir) + rn * (rrigrid(ir+1) - rrigrid(ir))

  end subroutine find_position_in_cell





 subroutine period(phi)
    ! Ensures phi lies in [0, 2π)
    real(8), intent(inout) :: phi
    real(8), parameter :: pi = 3.14159265358979323846264338328d0

    phi = mod(phi, 2.0d0 * pi)
    if (phi < 0.0d0) phi = phi + 2.0d0 * pi

  end subroutine period



function bplanck3(temp, nu) result(bval)
    implicit none
    real(8), intent(in) :: temp, nu
    real(8)             :: bval
    real(8), parameter  :: c1 = 1.47455d-47, c2 = 4.7989d-11

    if (temp == 0.d0) then
        bval = 0.d0
    else
        bval = c1 * nu**3 / (exp(c2 * nu / temp) - 1.d0) + 1.d-290
    end if

end function bplanck3





subroutine make_temperature_grid(ntemp, temp0, temp1)
    implicit none
    integer, intent(in)  :: ntemp
    real(8), intent(in) :: temp0, temp1
    !integer, intent(out) :: status

    integer :: itemp
    real(8) :: exponent

    ! Initialize status
    !status = 0

    ! Check grid size limit
    if (ntemp > DB_NTEMP_MAX) then
      write(*,*) "TEMPERATURE DATABASE MISMATCH"
      return
    endif

    ! Store number of temperature points
    db_ntemp = ntemp

    ! Generate logarithmic temperature grid
    do itemp = 1, ntemp
      exponent = (itemp - 1.0d0) / (ntemp - 1.0d0)
      db_temp(itemp) = temp0 * (temp1 / temp0) ** exponent
    end do

  end subroutine make_temperature_grid



function compute_dusttemp_energy_bd(ener,ispec,isize) result(temp)
  implicit none
  integer, intent(in) :: ispec, isize
  double precision, intent(in) :: ener
  double precision :: temp
  integer :: itemp
  double precision :: eps

  ! Search for index in db_enertemp where ener fits
  call hunt(db_enertemp(:,isize,ispec), db_ntemp, ener, itemp)

  ! Check for energy outside tabulated range
  if (itemp >= db_ntemp) then
    write(*,*) 'ERROR: Too high temperature discovered' 
    stop 76823
  end if

  if (itemp <= 0) then
    ! Below lowest energy: extrapolate linearly
    eps = ener / db_enertemp(1,isize,ispec)
    if (eps > 1.0d0) stop 9911
    temp = eps * db_temp(1)
  else
    ! Within range: linear interpolation
    eps = (ener - db_enertemp(itemp,isize,ispec)) / &
          (db_enertemp(itemp+1,isize,ispec) - db_enertemp(itemp,isize,ispec))
    if (eps < 0.0d0 .or. eps > 1.0d0) stop 9912
    temp = (1.0d0 - eps) * db_temp(itemp) + eps * db_temp(itemp+1)
  end if

end function compute_dusttemp_energy_bd




function compute_dusttemp_coupled_bd(energy, rho, vol) result(temp)
  implicit none
  integer :: ispec, isize, itemp
  double precision, intent(in) :: energy(DUST_SIZE_MAX,DUST_SPECIES_MAX)
  double precision, intent(in) :: rho(DUST_SPECIES_MAX)
  double precision, intent(in) :: vol
  double precision :: temp, eps, entot, en1, en2

  ! Compute total energy per gram dust
  entot = 0.d0
  do ispec = 1, dust_nr_species
    do isize = 1, dust_nr_size(ispec)
      entot = entot + energy(isize, ispec)
    end do
  end do
  entot = entot / vol

  ! Find index in temperature-energy grid
  call hunt_temp(rho, entot, itemp)

  ! Check bounds
  if (itemp >= db_ntemp) then
    write(*,*) 'ERROR: Too high temperature discovered'
    stop 76823
  end if

  if (itemp <= 0) then
    ! Below lowest temp in database: linear extrapolation
    en2 = 0.d0
    do ispec = 1, dust_nr_species
      do isize = 1, dust_nr_size(ispec)
        en2 = en2 + rho(ispec) * db_enertemp(1, isize, ispec)
      end do
    end do
    eps = entot / en2
    if (eps > 1.d0) stop 59911
    temp = eps * db_temp(1)
  else
    ! Interpolation between itemp and itemp+1
    en1 = 0.d0
    en2 = 0.d0
    do ispec = 1, dust_nr_species
      do isize = 1, dust_nr_size(ispec)
        en1 = en1 + rho(ispec) * db_enertemp(itemp, isize, ispec)
        en2 = en2 + rho(ispec) * db_enertemp(itemp+1, isize, ispec)
      end do
    end do

    eps = (entot - en1) / (en2 - en1)
    if (eps < 0.d0 .or. eps > 1.d0) stop 99124

    temp = (1.d0 - eps) * db_temp(itemp) + eps * db_temp(itemp+1)
  end if

end function compute_dusttemp_coupled_bd




 function absevfunc(temp, cellalpha,&
                   cellvol, ptclener, fnu_old, fnu_diff) result(demis)
  implicit none
  ! Inputs
  doubleprecision, intent(in) :: temp
 ! doubleprecision, intent(in) :: freq_nu(:)
  doubleprecision, intent(in) :: cellalpha(FRSIZE_FREQ)
 ! doubleprecision, intent(in) :: dust_freq_wgt(:)
  doubleprecision, intent(in) :: cellvol, ptclener
  doubleprecision, intent(in) :: fnu_old(FRSIZE_FREQ)
  ! Output (modified)
  doubleprecision, intent(out) :: fnu_diff(FRSIZE_FREQ)
  doubleprecision :: demis
  ! Constants
  doubleprecision, parameter :: fourpi = 12.5663706143591729538505735331d0
  ! Locals
  integer :: inu
  double precision :: fnu

  demis = 0.d0
  do inu = 1, size(freq_nu)
      fnu           = cellalpha(inu) * bplanck3(temp, freq_nu(inu))
      fnu_diff(inu) = fnu - fnu_old(inu)
      demis         = demis + fnu_diff(inu) * dust_freq_wgt(inu)
  end do

  demis = demis * fourpi * cellvol
  demis = demis - ptclener
end function absevfunc



subroutine make_emiss_dbase()

  !use mod_planck         ! for bplanck3
  implicit none


  ! Locals
  integer :: ispec, isize, itemp, inu
  double precision :: ptclener, cellvol
  double precision:: cellalpha(FRSIZE_FREQ), fnu_old(FRSIZE_FREQ), fnu_diff(FRSIZE_FREQ)
  double precision:: diffemis(FRSIZE_FREQ), cumul(FRSIZE_FREQ + 1)


  ! Loop over species and sizes
  do ispec = 1, dust_nr_species
    do isize = 1, dust_nr_size(ispec)

      ptclener = 0.d0
      cellvol  = 1.d0

      ! Set absorption coefficients and zero old emissivities
      do inu = 1, freq_nr
        cellalpha(inu) = kappa_a(inu, isize, ispec)
        fnu_old(inu)   = 0.d0
      end do

      ! Compute emissivity and cumulative norms for each temperature
      do itemp = 1, db_ntemp

        ! Compute energy emitted at this temperature
        db_enertemp(itemp, isize, ispec) = absevfunc( &
            db_temp(itemp), cellalpha, &
            cellvol, ptclener, fnu_old, fnu_diff)

        ! Store emissivities for this temperature
        do inu = 1, freq_nr
          db_emiss(inu, itemp, isize, ispec) = fnu_diff(inu)
        end do

      end do

      ! First cumulative for itemp=1
      itemp = 1
      do inu = 1, freq_nr
        diffemis(inu) = db_emiss(inu, itemp, isize, ispec)
      end do
      cumul(1) = 0.d0
      do inu = 1, freq_nr
        cumul(inu+1) = cumul(inu) + diffemis(inu) * dust_freq_wgt(inu)
      end do
      do inu = 1, freq_nr+1
        db_cumulnorm(inu, itemp, isize, ispec) = cumul(inu) / cumul(freq_nr+1)
      end do

      ! Remaining cumulatives for itemp ≥ 2
      do itemp = 2, db_ntemp
        do inu = 1, freq_nr
          diffemis(inu) = db_emiss(inu, itemp, isize, ispec) - &
                          db_emiss(inu, itemp-1, isize, ispec)
        end do
        cumul(1) = 0.d0
        do inu = 1, freq_nr
          cumul(inu+1) = cumul(inu) + diffemis(inu) * dust_freq_wgt(inu)
        end do
        do inu = 1, freq_nr+1
          db_cumulnorm(inu, itemp, isize, ispec) = cumul(inu) / cumul(freq_nr+1)
        end do
      end do

    end do
  end do

  !deallocate(cellalpha, fnu_old, fnu_diff, diffemis, cumul)

end subroutine make_emiss_dbase




subroutine do_absorption_event(vol, rhodust, kappa_a, alpha_a, tdust, &
                               cumulener, cumulener_bk, enthres, &
                               freqdistr, energy, ir, it, inucur, &
                               imethod, ifast, iseed, idestroy)

  implicit none

  ! Arguments
  integer, intent(in)    :: ir, it, imethod, ifast
  integer, intent(inout) :: iseed, idestroy,inucur
  double precision, intent(in)    :: vol(FRSIZE_Y_SMALL, FRSIZE_X)
  double precision, intent(in)    :: energy, enthres
  double precision, intent(in)    :: alpha_a(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
  double precision, intent(in)    :: kappa_a(FRSIZE_FREQ, DUST_SIZE_MAX, DUST_SPECIES_MAX)
  double precision, intent(in)    :: rhodust(DUST_SPECIES_MAX, FRSIZE_Y_SMALL, FRSIZE_X)
  double precision, intent(inout) :: tdust(DUST_SIZE_MAX, DUST_SPECIES_MAX, FRSIZE_Y_SMALL, FRSIZE_X)
  double precision, intent(inout) :: cumulener(DUST_SIZE_MAX, DUST_SPECIES_MAX, FRSIZE_Y_SMALL, FRSIZE_X)
  double precision, intent(inout) :: cumulener_bk(DUST_SIZE_MAX, DUST_SPECIES_MAX, FRSIZE_Y_SMALL, FRSIZE_X)
  double precision, intent(inout) :: freqdistr(FRSIZE_FREQ, FRSIZE_Y_SMALL, FRSIZE_X)
 

  ! Local variables
  integer :: inu, isize, ispec
  double precision :: enerpart(DUST_SIZE_MAX, DUST_SPECIES_MAX)
  double precision :: rhodusttot, cumen, tempav, dum
  double precision :: temp(DUST_SIZE_MAX, DUST_SPECIES_MAX)
  logical :: calctemp
 ! external :: absevfunc

  ! External temperature calculation functions
  !double precision :: compute_dusttemp_energy_bd, compute_dusttemp_coupled_bd

  ! Safety check
  if (dust_freq_wgt(1) == 0.0d0) error stop 'dust_freq_wgt(1)=0'

  if ((imethod == 0) .or. (imethod == 1)) then
    print *, 'ERROR: imethod=0 and imethod=1 are obsolete'
    error stop
  end if

  idestroy = 0

  if (inucur >= 1) then
    do ispec = 1, dust_nr_species
      do isize = 1, dust_nr_size(ispec)
        enerpart(isize, ispec) = energy * rhodust(ispec, it, ir) * &
                                kappa_a(inucur, isize, ispec) / alpha_a(inucur, it, ir)
      end do
    end do

  else if (inucur == -1) then
    rhodusttot = sum(rhodust(1:dust_nr_species, it, ir))
    do ispec = 1, dust_nr_species
      if (dust_nr_size(ispec) > 1) error stop 'This qplus mode only allows dust_nr_size=1'
      enerpart(1, ispec) = energy * rhodust(ispec, it, ir) / rhodusttot
    end do

  else
    print *, 'INTERNAL ERROR: inucur = ', inucur
    error stop
  end if

  ! Sanity check: sum should match original energy
  dum = sum(enerpart(1:dust_nr_size(1), 1:dust_nr_species))
  if (abs(dum/energy - 1.0d0) > 1.0d-6) error stop 'Energy partitioning error'


  calctemp = .true.

  if (calctemp) then
    if (itempdecoup == 1) then
      do ispec = 1, dust_nr_species
          do isize = 1, dust_nr_size(ispec)
            cumen = cumulener(isize, ispec, it, ir) / (dust_rho(ispec, it, ir) * vol(it, ir))
            temp(isize, ispec) = compute_dusttemp_energy_bd(cumen, ispec, isize)
            tdust(isize, ispec, it, ir) = temp(isize, ispec)
          end do
      end do
    else
      tempav = compute_dusttemp_coupled_bd(cumulener(1, 1, it, ir), dust_rho(1, it, ir), vol(it, ir))
      do ispec = 1, dust_nr_species
        do isize = 1, dust_nr_size(ispec)
          temp(isize, ispec) = tempav
          tdust(isize, ispec, it, ir) = tempav
        end do
      end do
    end if
  else
    do ispec = 1, dust_nr_species
      do isize = 1, dust_nr_size(ispec)
        temp(isize, ispec) = tdust(isize, ispec, it, ir)
      end do
    end do
  end if

  ! Pick new frequency after absorption
  call pick_randomfreq_db(tdust(1, 1, it, ir), enerpart, inucur, iseed)

  ! Backup cumulative energies if recalculated
  if (calctemp) then
    do ispec = 1, dust_nr_species
      do isize = 1, dust_nr_size(ispec)
        cumulener_bk(isize, ispec, it, ir) = cumulener(isize, ispec, it, ir)
      end do
    end do
  end if


end subroutine do_absorption_event


subroutine walk_full_path(nr, nt, rrigrid, ttigrid, vol, &
                          alpha_a, alpha_s, kappa_a, rstar, &
                          tdust, scatsrc, meanint, &
                          cumulener, cumulener_bk, &
                          freqdistr, energy, inu, theta, &
                          imethod, ifast, enthres, iseed)
  
    implicit none

    ! Arguments
    integer, intent(in) :: nr, nt, imethod, ifast
    integer, intent(inout) :: inu, iseed
    double precision, intent(in) :: rstar, enthres, energy
    double precision, intent(in) :: rrigrid(FRSIZE_X+1), ttigrid(FRSIZE_Y_SMALL+1)
    double precision, intent(in) :: vol(FRSIZE_Y_SMALL,FRSIZE_X)
    double precision, intent(in) :: alpha_a(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X), alpha_s(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
    double precision, intent(in) :: kappa_a(FRSIZE_FREQ,DUST_SIZE_MAX,DUST_SPECIES_MAX)
    double precision, intent(inout) :: tdust(DUST_SIZE_MAX,DUST_SPECIES_MAX,FRSIZE_Y_SMALL,FRSIZE_X)
    double precision, intent(inout) :: scatsrc(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
    !double precision, intent(inout) :: miquant(:,:,:)
!#ifdef SAVE_MEANINT
    double precision, intent(inout) :: meanint(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
!#else
!    double precision, intent(inout) :: meanint
!#endif
    double precision, intent(inout) :: cumulener(DUST_SIZE_MAX,DUST_SPECIES_MAX,FRSIZE_Y_SMALL,FRSIZE_X)
    double precision, intent(inout) :: cumulener_bk(DUST_SIZE_MAX,DUST_SPECIES_MAX,FRSIZE_Y_SMALL,FRSIZE_X)
    double precision, intent(inout) :: freqdistr(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
    double precision, intent(inout) :: theta

    ! Local variables
    integer :: ir, it, isize, ispec, iqactive, idestroy, imirt
    double precision :: mu, phi, r, mu0
    double precision :: taupath, tau, dtau, ds, fr, g, albedo, absorb
    double precision :: rn, mu_out, phi_out, stot
    logical :: trylin
    real(8), parameter :: pi = 3.141592653589793d0

    ! Initialize
    stot = 0.0d0
    
    idestroy = 0
    trylin = .false.
    imirt = 1

    ! Select photon source
    rn = ran2(iseed)
        ! Photon from star
        if (rstar == 0.0d0) then
            mu = 1.0d0
            phi = 0.0d0
            r = 0.0d0
            rn = ran2(iseed)
           
            theta = acos(rn)
           
            call randomfreq(radbnd_starspec, inu, iseed)

            r = 1.0d30
            call huntalt(ttigrid, nt+1, theta, it)
            if (theta == ttigrid(it)) stop "ERROR 32001"

            if ((it < 1) .or. (it > nt)) return
            taupath = -log(1.0d0 - ran2(iseed))
            tau = 0.0d0

            do ir = 1, nr
                ds = rrigrid(ir+1) - rrigrid(ir)
                dtau = (alpha_a(inu,it,ir) + alpha_s(inu,it,ir)) * ds
                albedo = alpha_s(inu,it,ir) / (alpha_a(inu,it,ir) + alpha_s(inu,it,ir))
                absorb = 1.0d0 - albedo

                if (tau + dtau > taupath) then
                    fr = (taupath - tau) / dtau
                    r = rrigrid(ir) + fr * ds
                    stot = stot + fr * ds
                    do ispec = 1, dust_nr_species
                        do isize = 1, dust_nr_size(ispec)
                            cumulener(isize,ispec,it,ir) = cumulener(isize,ispec,it,ir) + &
                                    absorb * (taupath-tau) * energy * dust_rho(ispec,it,ir) * &
                                    kappa_a(inu,isize,ispec) / alpha_a(inu,it,ir)
                        end do
                    end do
                    exit
                end if

                stot = stot + ds
                tau = tau + dtau
                do ispec = 1, dust_nr_species
                        do isize = 1, dust_nr_size(ispec)
                            cumulener(isize,ispec,it,ir) = cumulener(isize,ispec,it,ir) + &
                                absorb * dtau * energy * dust_rho(ispec,it,ir) * &
                                kappa_a(inu,isize,ispec) / alpha_a(inu,it,ir)
                        end do
                end do
                iphotcount(it,ir) = iphotcount(it,ir) + 1
            end do
           print*, albedo
            if (ran2(iseed) < albedo) then
                g = 0.0d0
                mu = 2.0d0 * ran2(iseed) - 1.0d0
                phi = 2 * pi * ran2(iseed)
            else
                call do_absorption_event(vol, dust_rho, kappa_a, alpha_a, tdust, cumulener, cumulener_bk, &
                     enthres, freqdistr, energy, ir, it, inu, imethod, ifast, iseed, idestroy)
                if (idestroy /= 0) then
                    theta = -1.0d30
                    return
                end if
                mu = 2.0d0 * ran2(iseed) - 1.0d0
                phi = 2 * pi * ran2(iseed)
            end if

        else
            ! Extended star source: to be modernized similarly


            rn  = ran2(iseed)
           
            theta = acos(rn)

            r    = rrigrid(1) + 1d-3 * (rrigrid(2)-rrigrid(1))  
            !print*, theta, r

            if(r.gt.1d3*rstar) then
                  mu = 1.d0
              else
                  if(r.le.rstar) stop 28670
                  rn  = ran2(iseed)
                  mu0 = sqrt(rn)
                  mu  = sqrt(1.d0-((rstar/r)**2)*(1.d0-mu0**2))
              endif

              phi   = 2.0d0 * pi*ran2(iseed)
              call randomfreq(radbnd_starspec,inu,iseed)
              !print*, inu
        end if
    

    ! Walk photon through grid
    do
        taupath = -log(1.0d0 - ran2(iseed))
        call walk_cells_twodee(inu, r, theta, mu, phi, nr, nt, rrigrid, ttigrid, vol, &
                               alpha_a, alpha_s, kappa_a, taupath, scatsrc, meanint, &
                               cumulener, energy, imethod, stot, idestroy)
        if (r >= 1.0d29 .or. r < 0.0d0 .or. idestroy /= 0) exit

        call huntalt(rrigrid, nr+1, r, ir)
        call huntalt(ttigrid, nt+1, theta, it)
        if ((ir < 1) .or. (ir > nr) .or. (it < 1) .or. (it > nt)) stop "ERROR 16501"

        albedo = alpha_s(inu,it,ir) / (alpha_a(inu,it,ir) + alpha_s(inu,it,ir))
       
        if (ran2(iseed) < albedo) then
            mu = 2.0d0 * ran2(iseed) - 1.0d0
            phi = 2 * pi * ran2(iseed)
        else
            call do_absorption_event(vol, dust_rho, kappa_a, alpha_a, tdust, cumulener, cumulener_bk, &
                 enthres, freqdistr, energy, ir, it, inu, imethod, ifast, iseed, idestroy)
            if (idestroy /= 0) exit
            mu = 2.0d0 * ran2(iseed) - 1.0d0
            phi = 2 * pi * ran2(iseed)
        end if
    end do

end subroutine walk_full_path



subroutine walk_cells_twodee(inu, r, theta, mu, phi, &
                            nr, nt, rrigrid, ttigrid, vol, &
                            alpha_a, alpha_s, kappa_a, taupath, &
                            scatsrc, meanint, cumulener, &
                            energy, imethod, stot, &
                            idestroy)
  
  implicit none

  ! Arguments
  integer, intent(in) :: inu, nr, nt, imethod
  !logical, intent(in) :: trylin
  double precision, intent(inout) :: r, theta, mu, phi, stot
  double precision, intent(in) :: taupath, energy
  double precision, intent(in) :: rrigrid(FRSIZE_X + 1), ttigrid(FRSIZE_Y_SMALL + 1)
  double precision, intent(in) :: vol(FRSIZE_Y_SMALL,FRSIZE_X)
  double precision, intent(in) :: alpha_a(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X), alpha_s(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
  double precision, intent(in) :: kappa_a(FRSIZE_FREQ,DUST_SIZE_MAX,DUST_SPECIES_MAX)
  double precision, intent(inout) :: scatsrc(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
  double precision, intent(inout) :: meanint(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)

  double precision, intent(inout) :: cumulener(DUST_SIZE_MAX,DUST_SPECIES_MAX,FRSIZE_Y_SMALL,FRSIZE_X)
  integer, intent(out) :: idestroy

  ! Local variables
  integer :: ircell, itcell, ir1, it1, igrid, is, ispec, isize, imirt
  double precision :: rold, told, muold, phiold
  double precision :: r0, r1, t0, t1
  double precision :: tau, dtau, sold, s, albedo, absorb, fr
  double precision :: costinf, ds, b, z0
  double precision :: dummy, sinphi
  real(8):: t_cell, r_cell, r_r, mu_r, s_th, th_th, sdum
  integer:: nteq, ilr_th, ilr_r, ileftright, idir


  real(8), parameter :: pi = 3.141592653589793d0
  real(8), parameter :: pihalf = pi/2.d0
  real(8), parameter :: twopi= 2 * 3.141592653589793d0
  real(8), parameter :: huge_s = 1.d30


imirt = 1
  ! Example checks for bounds
  if (phi < 0.0d0 .or. phi > twopi) stop 'Invalid phi'
  if (r > rrigrid(nr+1)) stop 'R out of bounds'
  if (theta > ttigrid(nt+1) .or. theta < ttigrid(1)) stop 'Theta out of bounds'


  nteq = nt+1
     

  ! Flip theta, phi if needed (mirror theta handling)
  if (theta > pihalf) then
    theta = pi - theta
    phi = twopi - phi
  endif

  ! Locate initial cell indices (use helper module/subroutine)
  call huntalt(rrigrid, nr+1, r, ircell)
  call huntalt(ttigrid, nt+1, theta, itcell)

  ! Boundary safety checks and adjustments
  if (ircell < 1 .or. ircell > nr .or. itcell < 1 .or. itcell > nt) then
      if(itcell.gt.nt) then
         theta = ttigrid(nt+1)-1d-12
         itcell = nt
      else
         r = -1d30
         return
      endif
  endif

  ! Slightly nudge r, theta if on grid boundary
  if (r.eq.rrigrid(ircell)) then
    if (ircell > 1) then
      r = r * (1.0d0 - 1.0d-12)
      ircell = ircell - 1
    else
      r = r * (1.0d0 + 1.0d-12)
    endif
  endif

  if (theta.eq.ttigrid(itcell)) then
    if (itcell > 1) then
      theta = theta - 1.0d-12
      itcell = itcell - 1
    else
      theta = theta + 1.0d-12
    endif
  endif

  igrid = 1
  ! Initialize optical depth
  tau = 0.0d0
  idestroy = 0

 if(iphotcurr.gt.ilastphot(itcell,ircell)) then
          iphotcount(itcell,ircell) = iphotcount(itcell,ircell) + 1
          ilastphot(itcell,ircell)  = iphotcurr
 endif

  ! Save old photon state
  rold = r
  told = theta
  muold = mu
  phiold = phi

  ! Define initial crossing boundaries
  r0 = rrigrid(ircell)
  r1 = rrigrid(ircell+1)
  t0 = ttigrid(itcell)
  t1 = ttigrid(itcell+1)

  ! Determine next crossing using a helper subroutine
  call find_next_crossing(rold, r0, r1, told, t0, t1, muold, phiold, sold, &
                         r, theta, mu, phi, s, idir, ileftright, b, costinf, z0)


  if (r < 0.0d0) then 
      return
  end if

   if(r.gt.rrigrid(nr+1)) then
      write(*,*) 'R out of bounds'
      write(*,*) r,rrigrid(1),rrigrid(nr+1)
      write(*,*) ircell,r0,r1
      stop 7491
   endif
  ! Main loop over cells until scattering event or escape
  do is = 1, 2*nr + 2*nt + 2

    ! Handle out-of-bounds or flipping
    if (theta > pihalf) then
      r = -1.0d0 * 1.0d30
      return
    endif

    if(r.gt.rrigrid(nr+1)) then
       write(*,*) 'A1 = ',is,r,rrigrid(nr+1)
    endif

    !print*, r

    call huntalt(rrigrid, nr+1, r, ir1)
    call huntalt(ttigrid, nt+1, theta, it1)
    !print*,'imirt', imirt, ir1, it1


     if(ir1.gt.nr+1) then
              write(*,*) is,r,rrigrid(nr+1)
              stop 39991
          endif
          if(it1.gt.nt+1) then
              write(*,*) is,theta,ttigrid(nt+1)
              stop 39992
          endif

           if(igrid.eq.1) then
              if((ircell.lt.1).or.(ircell.gt.nr).or. &
                (itcell.lt.1).or.(itcell.gt.nt)) then
                  stop 22371
              endif
          else
              if((ircell.ge.1).and.(ircell.le.nr).and. &
                (itcell.ge.1).and.(itcell.le.nt)) then
                  stop 22372
              endif
          endif

    ! Calculate incremental optical depth and absorption
    ds = s - sold
        if(igrid.eq.1) then
              dtau = ( alpha_a(inu,itcell,ircell) + &
                     alpha_s(inu,itcell,ircell) ) * ds
          else
              dtau = 0.d0
          endif

          if(igrid.eq.1) then
              albedo = alpha_s(inu,itcell,ircell)/ &
                  (alpha_a(inu,itcell,ircell)+ &
                   alpha_s(inu,itcell,ircell))
              absorb = 1.d0-albedo
          endif
     
    ! Check if photon scattering occurs within this step
    if (tau + dtau > taupath) then
       if(igrid.eq.0) stop 94001
      fr = (taupath - tau) / dtau
      
      ! Deposit partial energy into dust cells for absorption
      do ispec = 1, dust_nr_species
        do isize = 1, dust_nr_size(ispec)
          cumulener(isize, ispec, itcell, ircell) = cumulener(isize, ispec, itcell, ircell) + &
                absorb * (taupath - tau) * energy * dust_rho(ispec, itcell, ircell) * &
                kappa_a(inu, isize, ispec) / alpha_a(inu, itcell, ircell)
        end do
      end do

      ! Update photon's position fractionally
      s = sold + fr * ds
      stot = stot + fr * ds

      ! Compute new position and angles
      r = sqrt(b**2 + s**2)
      dummy = (z0 + s * costinf) / r
      theta = acos(dummy)
      mu = s / r

      dummy = b**2 + s**2 - (z0 + s * costinf)**2
      if (dummy <= 0.0d0) stop 'Invalid geometry in phi calculation'

      dummy = b * sqrt(dummy)
      sinphi = (b**2 * costinf - z0 * s) / dummy
      if (phi > 0.0d0) then
        phi = asin(sinphi)
      else
        phi = pi - asin(sinphi)
      endif

      call period(phi)  ! normalize phi angle

      return
    endif

    ! Accumulate optical depth and deposit full energy
    !tau = tau + dtau
    stot = stot + ds

if(igrid.eq.1) then
    do ispec = 1, dust_nr_species
      do isize = 1, dust_nr_size(ispec)
        cumulener(isize, ispec, itcell, ircell) = cumulener(isize, ispec, itcell, ircell) + &
              absorb * dtau * energy * dust_rho(ispec, itcell, ircell) * &
              kappa_a(inu, isize, ispec) / alpha_a(inu, itcell, ircell)
      end do
    end do
end if


if((theta.eq.pihalf)) then
    if(phi.lt.pi) stop 64327
    if(idir.ne.2) stop 64328
    if(ileftright.ne.1) stop 64329
    phi = twopi - phi
    ileftright = -1
endif

   ircell = -10001
   itcell = -10001


if (idir == 1) then

    if (ileftright == 1) then
        if (ir1 == nr+1) then
            costinf = mu * cos(theta) + sqrt(1.d0 - mu*mu) * sin(theta) * sin(phi)
            r       = 1.d30
            theta   = acos(costinf)
            if (theta > pihalf) then 
               theta = pi - theta
            end if
            mu      = 1.d0
            phi     = 0.d0
            return
        endif

        if (ir1 <= 0 .or. ir1 > nr) stop 83761
        if (it1 <= 0 .or. it1 > nt) then
            if(it1.eq.nt+1) then
               theta=pihalf-1d-12
               it1=nt
            else
               theta=ttigrid(1)+1d-12
               it1=1
            endif
        endif

        r0 = rrigrid(ir1)
        r1 = rrigrid(ir1+1)
        t0 = ttigrid(it1)
        t1 = ttigrid(it1+1)
        ircell = ir1
        itcell = it1

    else  ! ileftright ≠ 1

        if (ir1 <= 0 .or. ir1 > nr+1) stop 83771
        if (it1 <= 0 .or. it1 > nt) then
            if(it1.eq.nt+1) then
               theta=pihalf-1d-12
               it1=nt
            else
               theta=ttigrid(1)+1d-12
               it1=1
            endif
        endif
         if(ir1.gt.1) then
            r0 = rrigrid(ir1-1)
         else
            r0 = 1d-10*rrigrid(1)
         endif
       !r0 = merge(rrigrid(ir1-1), 1.d-10 * rrigrid(1), ir1 > 1)
        r1 = rrigrid(ir1)
        t0 = ttigrid(it1)
        t1 = ttigrid(it1+1)
        ircell = ir1-1
        itcell = it1

    endif

else  ! idir ≠ 1

    if (ileftright == 1) then
        if (ir1 <= -1 .or. ir1 > nr) stop 83781
        if (it1 <= 0  .or. it1 > nt) stop 83782
        if (it1 == nteq) stop 7129
        if (ir1 < 0) stop 71291

        r1 = rrigrid(ir1+1)
        if(ir1.eq.0) then
           r0 = 1d-10 * r1
        else
           r0 = rrigrid(ir1)
        endif
       ! r0 = merge(rrigrid(ir1), 1.d-10 * r1, ir1 /= 0)
        t0 = ttigrid(it1)
        t1 = ttigrid(it1+1)
        ircell = ir1
        itcell = it1

    else  ! ileftright ≠ 1

        if (ir1 <= -1 .or. ir1 > nr) stop 83791
        if (it1 <= 0  .or. it1 > nt+1) stop 83792
        if (ir1 < 0) stop 71292

        r1 = rrigrid(ir1+1)
        !r0 = merge(rrigrid(ir1), 1.d-10 * r1, ir1 /= 0)
        !t0 = merge(ttigrid(it1-1), 1.d-10, it1 > 1)
        if(ir1.eq.0) then
           r0 = 1d-10 * r1
        else
           r0 = rrigrid(ir1)
        endif
        if(it1.gt.1) then
           t0 = ttigrid(it1-1)
        else
           t0 = 1.d-10
        endif
        t1 = ttigrid(it1)
        ircell = ir1
        itcell = it1-1

    endif

endif


if (ircell<=-10000 .or. itcell<=-10000) stop 33377
if (r<r0 .or. r>r1 .or. theta<t0 .or. theta>t1) stop 13

! Check grid validity
igrid = 1
if (ircell<1 .or. ircell>nr) igrid = 0
if (itcell<1 .or. itcell>nt) igrid = 0

if (igrid==1) then
    r_cell = 0.5d0 * (rrigrid(ircell) + rrigrid(ircell+1))
    t_cell = 0.5d0 * (ttigrid(itcell) + ttigrid(itcell+1))
    if (r_cell<r0 .or. r_cell>r1 .or. t_cell<t0 .or. t_cell>t1) then
        write(*,*) 'ERROR: Cell wrong...'
        stop 33732
    end if

    if (iphotcurr > ilastphot(itcell, ircell)) then
        iphotcount(itcell, ircell) = iphotcount(itcell, ircell) + 1
        ilastphot(itcell, ircell)  = iphotcurr
    end if
end if

! Propagation
tau    = tau + dtau
rold   = r
told   = theta
muold  = mu
phiold = phi
sold   = s

call find_next_crossing(rold, r0, r1, told, t0, t1, muold, phiold, sdum, &
                        r, theta, mu, phi, s, idir, ileftright, b, costinf, z0)


if (r < 0.d0) then
   return
end if

! Final domain checks
if (r > rrigrid(nr+1)) stop 49001
if (theta > ttigrid(nt+1) .or. theta < ttigrid(1)) stop 49002
if (idir==1 .and. (theta==ttigrid(nt+1) .or. theta==ttigrid(1))) stop 49003
if (idir==2 .and. (r==rrigrid(nr+1) .or. r==rrigrid(1))) stop 49004

if (abs(sold/sdum-1.d0) > 1.d-3) stop 59232
if (s < sold) then
    write(*,*) s, sold, abs((s-sold)/s)
    stop 59233
end if

end do


end subroutine walk_cells_twodee




subroutine find_next_crossing(r, r0, r1, theta, t0, t1, &
                              mu, phi, s, r_next, th_next, mu_next, &
                              phi_next, s_next, idir, ileftright, b, costinf, z0)
  implicit none
  ! Inputs
  double precision, intent(in) :: r, r0, r1, theta, t0, t1
  double precision, intent(in) :: mu, phi
  ! Outputs
  double precision, intent(out) :: r_next, th_next, mu_next, phi_next, s_next
  integer, intent(out) :: idir, ileftright
  double precision, intent(out) :: b, costinf, z0
  ! Locals
  double precision :: s, r_th, th_th, s_th, mu_th, phi_th
  double precision :: r_r, th_r, s_r, mu_r, phi_r
  double precision :: dummy, dum, bdum, cdum, sinphi, costh, qdum
  double precision :: sol1, sol2, ssol1, ssol2
  double precision :: pi, pihalf, twopi
  integer :: ilr_r, ilr_th

  ! Constants
  parameter(pi = 3.14159265358979323846264338328d0)
  parameter(pihalf = 1.57079632679489661923132169164d0)
  parameter(twopi = 6.28318530717958647692528676656d0)

  ! Sanity checks
  if (theta > pihalf .or. t1 > pihalf) stop 'Error: theta or t1 > π/2'
  if (r1 < r0 .or. t1 < t0) stop 'Error: invalid cell boundaries'
  if (r < r0 .or. r > r1) stop 'Error: r out of cell bounds'
  if (theta < t0 .or. theta > t1) stop 'Error: theta out of cell bounds'
  if ((r == r0 .and. mu < 0.d0) .or. (r == r1 .and. mu > 0.d0)) stop 'Error: invalid radial direction'
  if ((theta == t0 .and. phi < pi) .or. (theta == t1 .and. phi > pi)) stop 'Error: invalid theta boundary'

  ! Initialize
  sol1 = 0.d0
  sol2 = 0.d0
  s    = r * mu
  b    = r * sqrt(1.d0 - mu * mu)
  costinf = mu * cos(theta) + sqrt(1.d0 - mu * mu) * sin(theta) * sin(phi)
  z0   = r * ( (1.d0 - mu*mu) * cos(theta) - mu * sqrt(1.d0 - mu*mu) * sin(theta) * sin(phi) )

  ! === Radial crossing ===
  if (mu < 0.d0) then
     if (r0 >= b) then
        r_r = r0
        mu_r = -sqrt(1.d0 - (b / r_r)**2)
        ilr_r = -1
     else
        r_r = r1
        mu_r = sqrt(1.d0 - (b / r_r)**2)
        ilr_r = 1
     end if
  else
     r_r = r1
     mu_r = sqrt(1.d0 - (b / r_r)**2)
     ilr_r = 1
  end if

  s_r = r_r * mu_r
  dummy = b**2 + s_r**2 - (z0 + s_r * costinf)**2
  if (dummy <= 0.d0) stop 'Error: invalid geometry at radial crossing'
  dummy = b * sqrt(dummy)
  sinphi = (b**2 * costinf - z0 * s_r) / dummy
  phi_r = merge(asin(sinphi), pi - asin(sinphi), phi <= 0.d0)
  call period(phi_r)
  th_r = acos((z0 + s_r * costinf) / r_r)

  ! === Theta crossing ===
  call compute_theta_crossing(r, b, z0, costinf, theta, t0, t1, phi, s, &
                              s_th, th_th, ilr_th)
  !print*, s_th
  ! Compute r, mu, phi at theta crossing
  if (s_th < 1d30) then
     r_th  = sqrt(b**2 + s_th**2)
     mu_th = s_th / r_th
     dummy = b**2 + s_th**2 - (z0 + s_th * costinf)**2
     if (dummy <= 0.d0) stop 'Error: invalid geometry at theta crossing'
     dummy = b * sqrt(dummy)
     sinphi = (b**2 * costinf - z0 * s_th) / dummy
     phi_th = merge(asin(sinphi), pi - asin(sinphi), phi <= 0.d0)
     call period(phi_th)
  end if
  !print*, r_th, mu_th
  ! === Choose next crossing ===
  if (s_r < s_th) then
     r_next = r_r; th_next = th_r; mu_next = mu_r; phi_next = phi_r
     s_next = s_r; idir = 1; ileftright = ilr_r
  else
     r_next = r_th; th_next = th_th; mu_next = mu_th; phi_next = phi_th
     s_next = s_th; idir = 2; ileftright = ilr_th
  end if

  ! Consistency check
  if (abs((s_next - s) / s) < 1d-14 .or. s_next < s .or. &
      r_next < r0 .or. r_next > r1 .or. th_next < t0 .or. th_next > t1) then
     r_next = -1d30
     return
  end if

end subroutine find_next_crossing




subroutine compute_theta_crossing(r, b, z0, costinf, theta, t0, t1, phi, s, &
                                  s_th, th_th, ilr_th)
  implicit none

  ! Inputs
  double precision, intent(in) :: r, b, z0, costinf, theta, t0, t1, phi, s
  double precision :: pihalf, costeps

  ! Outputs
  double precision, intent(out) :: s_th, th_th
  integer, intent(out) :: ilr_th

  ! Locals
  double precision :: costh, dum, bdum, cdum, qdum
  double precision :: sol1, sol2, ssol1, ssol2, dummy, sinphi, pi

  parameter(pi = 3.14159265358979323846264338328d0)
  parameter(pihalf = 1.57079632679489661923132169164d0)
  parameter(costeps = 1.0d-6)

  s_th = 1d30
  ilr_th = 0

  if (phi < pi) then
    ! Moving towards the pole
    costh = cos(t0)
    dum = costh**2 - costinf**2
    if (dum == 0.d0) stop 83218

    bdum = -2.d0 * z0 * costinf / dum
    cdum = ( (b*costh)**2 - z0**2 ) / dum
    qdum = bdum**2 - 4.d0 * cdum

    if (qdum >= 0.d0) then
      sol1 = 0.5d0 * ( -bdum - sqrt(qdum) )
      sol2 = 0.5d0 * ( -bdum + sqrt(qdum) )

      if (sol2 > s) then
        if (sol1 < s) then
          s_th = sol2
        else
          s_th = sol1
        end if
        th_th  = t0
        ilr_th = -1

      else
        s_th = 1d30
        if (abs(costinf) > cos(t0) .or. abs(costinf) < cos(t1)) then
          write(*,*) 'Monte Carlo: Internal consistency failure'
          stop 44333
        end if
        return
      end if

    else
      ! No real solution. Check if t1 = pi/2
      if (t1 == pihalf) then
        s_th = -z0 / costinf
        if ((abs((s_th-s)/s) > 1.d-8) .and. (s_th > s)) then
          th_th  = t1
          ilr_th = 1
        else
          s_th = 1d30
          return
        end if

      else
        ! Avoid numerical issue: compute with t1+costeps
        costh = cos(t1+costeps)
        dum = costh**2 - costinf**2
        if (dum == 0.d0) stop 83218
        bdum = -2.d0 * z0 * costinf / dum
        cdum = ( (b*costh)**2 - z0**2 ) / dum
        qdum = bdum**2 - 4.d0 * cdum
        if (qdum < 0.d0) stop 61321

        ssol1 = 0.5d0 * ( -bdum - sqrt(qdum) )
        ssol2 = 0.5d0 * ( -bdum + sqrt(qdum) )

        ! Now at t1
        costh = cos(t1)
        dum = costh**2 - costinf**2
        if (dum == 0.d0) stop 83218
        bdum = -2.d0 * z0 * costinf / dum
        cdum = ( (b*costh)**2 - z0**2 ) / dum
        qdum = bdum**2 - 4.d0 * cdum
        if (qdum < 0.d0) stop 61322

        sol1 = 0.5d0 * ( -bdum - sqrt(qdum) )
        sol2 = 0.5d0 * ( -bdum + sqrt(qdum) )

        if (ssol1 > s) then
          s_th = sol1
          th_th = t1
          ilr_th = 1
        else if (ssol2 > s) then
          s_th = sol2
          th_th = t1
          ilr_th = 1
        else
          s_th = 1d30
          if (abs(costinf) > cos(t0) .or. abs(costinf) < cos(t1+costeps)) then
            write(*,*) 'Monte Carlo: Internal consistency failure'
            write(*,*) cos(t0), cos(t1+costeps), costinf
            stop 44332
          end if
          return
        end if

      end if
    end if

  else
    ! Moving towards equator
    if (t1 == pihalf) then
      s_th = -z0 / costinf
      if (s_th > s) then
        th_th  = t1
        ilr_th = 1
      else
        s_th = 1d30
        return
      end if

    else
      costh = cos(t1)
      dum = costh**2 - costinf**2
      if (dum == 0.d0) stop 83218
      bdum = -2.d0 * z0 * costinf / dum
      cdum = ( (b*costh)**2 - z0**2 ) / dum
      qdum = bdum**2 - 4.d0 * cdum
      if (qdum < 0.d0) then
        write(*,*) 'Discriminant negative: ', r, theta, phi, z0, costh, costinf
        stop 99538
      end if

      sol1 = 0.5d0 * ( -bdum - sqrt(qdum) )
      sol2 = 0.5d0 * ( -bdum + sqrt(qdum) )

      if (sol2 > s) then
        if (sol1 < s) then
          s_th = sol2
        else
          s_th = sol1
        end if
        th_th  = t1
        ilr_th = 1
      else
        s_th = 1d30
        if (abs(costinf) > cos(t0) .or. abs(costinf) < cos(t1)) then
          write(*,*) 'Monte Carlo: Internal consistency failure'
          stop 44331
        end if
        return
      end if
    end if
  end if


end subroutine compute_theta_crossing



SUBROUTINE pick_randomfreq_db(temp, enerpart, inupick, iseed)

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: inupick
  INTEGER, INTENT(INOUT) :: iseed
  DOUBLE PRECISION, INTENT(IN) :: temp(DUST_SIZE_MAX, DUST_SPECIES_MAX)
  DOUBLE PRECISION, INTENT(IN) :: enerpart(DUST_SIZE_MAX, DUST_SPECIES_MAX)

  DOUBLE PRECISION :: enercum(DUST_SIZE_MAX+1, DUST_SPECIES_MAX)
  DOUBLE PRECISION :: enercumm(DUST_SPECIES_MAX+1)
  DOUBLE PRECISION :: cumul(FRSIZE_FREQ+1)
  DOUBLE PRECISION :: rn, eps
  INTEGER :: inu, ispec, isize, itemp

  ! Build cumulative energy distributions for sizes within each species
  enercumm(1) = 0.0D0
  DO ispec = 1, dust_nr_species
    enercum(1, ispec) = 0.0D0
    DO isize = 1, dust_nr_size(ispec)
      enercum(isize+1, ispec) = enercum(isize, ispec) + enerpart(isize, ispec)
    END DO
    enercumm(ispec+1) = enercumm(ispec) + enercum(dust_nr_size(ispec)+1, ispec)
  END DO

  ! Select random species
  IF (dust_nr_species > 1) THEN
    rn = ran2(iseed) * enercumm(dust_nr_species+1)
    CALL hunt(enercumm, dust_nr_species, rn, ispec)
    IF (ispec < 1 .OR. ispec > dust_nr_species) STOP 50209
  ELSE
    ispec = 1
  END IF

  ! Select random size within species
  IF (dust_nr_size(ispec) > 1) THEN
    rn = ran2(iseed) * enercum(dust_nr_size(ispec)+1, ispec)
    CALL hunt(enercum(1, ispec), dust_nr_size(ispec), rn, isize)
    IF (isize < 1 .OR. isize > dust_nr_size(ispec)) STOP 50210
  ELSE
    isize = 1
  END IF

  ! Find temperature index in database
  CALL hunt(db_temp, db_ntemp, temp(isize, ispec), itemp)

  ! Check temperature bounds
  IF (itemp >= db_ntemp) THEN
    WRITE(*,*) 'ERROR: Too high temperature discovered'
    STOP 77823
  ELSE IF (itemp <= 0) THEN
    itemp = 1
    eps = 0.0D0
  ELSE
    eps = (temp(isize, ispec) - db_temp(itemp)) / (db_temp(itemp+1) - db_temp(itemp))
    IF (eps > 1.0D0 .OR. eps < 0.0D0) STOP 94941
  END IF

  ! Build cumulative probability distribution for frequencies
  DO inu = 1, freq_nr+1
    cumul(inu) = (1.0D0 - eps) * db_cumulnorm(inu, itemp, isize, ispec) + &
                 eps * db_cumulnorm(inu, itemp+1, isize, ispec)
  END DO

  ! Final sanity check on cumulative normalization
  IF (ABS(cumul(freq_nr+1) - 1.0D0) > 1.0D-11) STOP 43465

  ! Pick random frequency bin
  rn = ran2(iseed)
  CALL hunt(cumul, freq_nr, rn, inupick)
  IF (inupick < 1 .OR. inupick > freq_nr) STOP 8189

END SUBROUTINE pick_randomfreq_db




SUBROUTINE hunt_temp(rho, energy, jlo)
 
  IMPLICIT NONE

  INTEGER, INTENT(INOUT) :: jlo
  DOUBLE PRECISION, INTENT(IN) :: energy
  DOUBLE PRECISION, INTENT(IN) :: rho(DUST_SPECIES_MAX)

  INTEGER :: n, inc, jhi, jm, ispec, isize
  LOGICAL :: ascnd
  DOUBLE PRECISION :: x, xx

  ! Initialize
  n = db_ntemp
  x = energy
  ascnd = .TRUE.

  IF (jlo <= 0 .OR. jlo > n) THEN
    jlo = 0
    jhi = n + 1
    GOTO 100
  END IF

  inc = 1

  ! Compute xx at initial jlo
  xx = 0.0D0
  DO ispec = 1, dust_nr_species
    DO isize = 1, dust_nr_size(ispec)
      xx = xx + rho(ispec) * db_enertemp(jlo, isize, ispec)
    END DO
  END DO

  ! Upward hunt if x ≥ xx
  IF ((x >= xx) .EQV. ascnd) THEN
    DO
      jhi = jlo + inc
      IF (jhi > n) THEN
        jhi = n + 1
        EXIT
      END IF

      xx = 0.0D0
      DO ispec = 1, dust_nr_species
        DO isize = 1, dust_nr_size(ispec)
          xx = xx + rho(ispec) * db_enertemp(jhi, isize, ispec)
        END DO
      END DO

      IF ((x >= xx) .EQV. ascnd) THEN
        jlo = jhi
        inc = inc * 2
      ELSE
        EXIT
      END IF
    END DO

  ELSE
    ! Downward hunt if x < xx
    jhi = jlo
    DO
      jlo = jhi - inc
      IF (jlo < 1) THEN
        jlo = 0
        EXIT
      END IF

      xx = 0.0D0
      DO ispec = 1, dust_nr_species
        DO isize = 1, dust_nr_size(ispec)
          xx = xx + rho(ispec) * db_enertemp(jlo, isize, ispec)
        END DO
      END DO

      IF ((x < xx) .EQV. ascnd) THEN
        jhi = jlo
        inc = inc * 2
      ELSE
        EXIT
      END IF
    END DO
  END IF

100 CONTINUE

  ! Bisection loop
  DO WHILE (jhi - jlo > 1)
    jm = (jhi + jlo) / 2

    xx = 0.0D0
    DO ispec = 1, dust_nr_species
      DO isize = 1, dust_nr_size(ispec)
        xx = xx + rho(ispec) * db_enertemp(jm, isize, ispec)
      END DO
    END DO

    IF ((x > xx) .EQV. ascnd) THEN
      jlo = jm
    ELSE
      jhi = jm
    END IF
  END DO

END SUBROUTINE hunt_temp






end module montecarlo
