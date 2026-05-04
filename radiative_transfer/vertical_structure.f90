module vertical_structure
use common_boundary
use common_grid
use common_dust
use common_montecarlo
use common_vstruct
!use mc_arrays
  implicit none
  !private
  !public :: calc_sigmadust

  double precision, parameter :: pi = 3.1415926535897932385d0

contains



subroutine read_vstructinp()
 ! use mod_grid
 ! use mod_dust_globals
 ! use mod_vertstruct_globals
  implicit none

  integer :: iformat, nspec, ispec, idum
  logical :: fex

  inquire(file='vstruct.inp', exist=fex)

  if (fex) then
    open(unit=5, file='vstruct.inp', status='old', action='read')
    read(5, *) iformat
    read(5, *) nspec
    if (nspec /= dust_nr_species) then
      write(*,*) 'ERROR: Number of dust species in vstruct.inp (', nspec, &
                 ') does not match dustopac.inp (', dust_nr_species, ')'
      close(5)
      stop
    endif

    do ispec = 1, nspec
      read(5, *) idum
      if (idum < 0) then
        write(*,*) 'WARNING: Negative iteration flag for species ', ispec, '. Setting to zero.'
        idum = 0
      endif

      ivstr_dostruct(ispec) = idum

      if (idum /= 0) then
        write(*,*) '...For dust species ', ispec, ' vertical structure iteration enabled.'
      endif
    end do

    close(5)

  else
    write(*,*) 'vstruct.inp not found. Doing vertical structure iteration for all dust species.'
    do ispec = 1, dust_nr_species
      ivstr_dostruct(ispec) = 1
    end do
  endif

end subroutine read_vstructinp




  subroutine calc_sigmadust()
!use boundaries
!use mod_grid
!use mod_dust_globals
!use mc_spec_module
!use mod_vertstruct_globals
    implicit none
   
    ! Arguments
    !integer, intent(in) :: irsi_frsizex, irsi_frsizey, dust_nr_species
   ! double precision, intent(in)  :: dust_rho(:, :, :)
   ! double precision, intent(in)  :: rsi_x_c(:, :)
   ! double precision, intent(in)  :: rsi_x_i(:, :)
   ! double precision, intent(out) :: vstr_sigmadust(:, :)
   ! double precision, intent(out) :: vstr_massdust(:)
   ! double precision, intent(out) :: vstr_massdusttot

    ! Local variables
    integer :: ir, it, ispec, nt
    double precision :: dum
    real(8), parameter:: pi=3.1415926535897932385d0

    ! Number of vertical cells on one side
    nt = irsi_frsizey / 2
    !print*, nt

    ! Compute surface density for each radial cell and dust species
    do ir = 1, irsi_frsizex
       do ispec = 1, dust_nr_species
          dum = 0.d0
          do it = 1, nt
             dum = dum + dust_rho(ispec, it, ir) * rsi_x_c(ir, 1) * &
                    abs(rsi_x_i(it, 2) - rsi_x_i(it + 1, 2))
          end do
          vstr_sigmadust(ir, ispec) = 2.d0 * dum
       end do
    end do
    
    ! Compute total dust mass for each dust species and total
    vstr_massdusttot = 0.d0
  !  vstr_massdust(:) = 0.d0

    do ispec = 1, dust_nr_species
        vstr_massdust(ispec) = 0.d0
       do ir = 1, irsi_frsizex
          vstr_massdust(ispec) = vstr_massdust(ispec) + &
               pi * (rsi_x_i(ir + 1, 1)**2 - rsi_x_i(ir, 1)**2) * &
               vstr_sigmadust(ir, ispec)
       end do
       vstr_massdusttot = vstr_massdusttot + vstr_massdust(ispec)
    end do

  end subroutine calc_sigmadust






!module constants_mod
!  implicit none
!  integer, parameter :: dp = selected_real_kind(15, 300)
!  real(dp), parameter :: pi = 3.1415926535897932385d0
!  real(dp), parameter :: GG = 6.672d-8
!  real(dp), parameter :: mp = 1.6726d-24
!  real(dp), parameter :: kk = 1.3807d-16
!end module constants_mod


subroutine vertstruct_integrate(tdust, ivstrt, error, iwarn)
use common_boundary
use common_grid
use common_dust
use common_vstruct
  implicit none

  ! Arguments
  integer, intent(in) :: ivstrt
  real(8), intent(inout) :: tdust(:,:,:,:)
  real(8), intent(out) :: error
  integer, intent(out) :: iwarn

  ! Local variables
  integer :: ispec, ir, it, nr, nt
  real(8) :: mugas, gravc, r, dz, dlgt, grv, errloc
  real(8) :: rho0(FRSIZE_Y), z(FRSIZE_Y), zi(FRSIZE_Y+1)
  real(8) :: sig0

  
  real(8), parameter :: pi = 3.1415926535897932385d0
  real(8), parameter :: GG = 6.672d-8
  real(8), parameter :: mp = 1.6726d-24
  real(8), parameter :: kk = 1.3807d-16

  ! Initialization
  error = 0.d0
  mugas = 2.3d0
  nr = irsi_frsizex
  nt = irsi_frsizey / 2

  ! Loop over species
  do ispec = 1, dust_nr_species
    if (ivstr_dostruct(ispec) /= 0) then
      do ir = 1, nr
        r = rsi_x_c(ir,1)

        ! Compute z and zi grids
        do it = 1, nt
          z(it)  = r * (0.5d0 * pi - rsi_x_c(it,2))
        end do
        do it = 1, nt+1
          zi(it) = r * (0.5d0 * pi - rsi_x_i(it,2))
        end do

        gravc = (GG * radbnd_mstar / r**3) * mugas * mp / kk

        ! Integrate density structure
        rho0(nt) = 1.d0
        do it = nt, 2, -1
          dz   = z(it-1) - z(it)
          dlgt = (log(tdust(1,ivstrt,it-1,ir)) - log(tdust(1,ivstrt,it,ir))) / dz
          grv  = 0.5d0 * gravc * (z(it-1)/tdust(1,ivstrt,it-1,ir) + z(it)/tdust(1,ivstrt,it,ir))
          rho0(it-1) = rho0(it) * exp( -dz * (grv + dlgt) )
        end do

        ! Compute Sigma
        sig0 = 0.d0
        do it = nt, 1, -1
          dz   = abs(zi(it) - zi(it+1))
          sig0 = sig0 + rho0(it) * dz
        end do
        sig0 = 2.d0 * sig0

        ! Normalize density
        do it = nt, 1, -1
          rho0(it) = rho0(it) * vstr_sigmadust(ir,ispec) / sig0
          if (rho0(it) < 1.d-90) rho0(it) = 1.d-90
        end do

        ! Midplane error
        errloc = 2.d0 * abs(rho0(nt)-dust_rho(ispec,nt,ir)) / (rho0(nt) + dust_rho(ispec,nt,ir))
        if (errloc > error) error = errloc

        ! Update rho array
        do it = nt, 1, -1
          dust_rho(ispec,it,ir) = rho0(it)
        end do

        ! Warning if midplane discontinuity too large
        if (abs(rho0(nt)-rho0(nt-1))/(rho0(nt)+rho0(nt-1)) > 0.3d0) then
          iwarn = 1
        else
          iwarn = 0
        end if

      end do
    end if
  end do

end subroutine vertstruct_integrate


end module vertical_structure


