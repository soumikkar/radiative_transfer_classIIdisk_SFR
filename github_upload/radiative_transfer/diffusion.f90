module diffusion
  use configure
  use common_grid
  use common_boundary
  use common_dust
  use common_montecarlo
  use common_diffusion
  use numerical_receipe
implicit none

contains 
subroutine diffusion_main(nr, nt, nf, rc, tc, freq, alpha_a, alpha_s, tdust, errtol, ierror)
 
  !use mod_grid
  !use boundaries
  ! use mod_dust_globals
  ! use mc_spec_module
  ! use diffusion_array
  implicit none

  ! Arguments
  integer, intent(in) :: nr, nt, nf
  real(8), intent(in) :: rc(FRSIZE_X), tc(FRSIZE_Y_SMALL)
  real(8), intent(in) :: freq(FRSIZE_FREQ)
  real(8), intent(in) :: alpha_a(FRSIZE_FREQ, FRSIZE_Y_SMALL, FRSIZE_X)
  real(8), intent(in) :: alpha_s(FRSIZE_FREQ, FRSIZE_Y_SMALL, FRSIZE_X)
  real(8), intent(inout) :: tdust(DUST_SIZE_MAX, DUST_SPECIES_MAX, FRSIZE_Y_SMALL, FRSIZE_X)
  real(8), intent(in) :: errtol
  integer, intent(out) :: ierror

  ! Local variables
  integer :: ir, it, idiff, ndiff0, idiff_r, idiff_l,iter,isize,ispec
  real(8) :: tguess
  logical :: gotit
  real(8):: dxc,dxm,dxp,dcp,dcm,solmax, dummy, error,thghost
  integer :: indx(FRSIZE_NDIFF)

  real(8), parameter :: pi = 3.1415926535897932385d0
  integer, parameter :: niter = 20
  !integer :: number_invalid

  ! Reset error and initialize
  ierror = 0
  tguess = 100.0d0

#ifdef CHECK_NUMBERS
  do ispec = 1, dust_nr_species
    do ir = 1, nr
      do it = 1, nt
        if (number_invalid(tdust(1,ispec,it,ir)) /= 0) then
          write(*,*) 'ERROR: Invalid temperature found before diffusion'
          write(*,*) ir, it, ispec, tdust(1,ispec,it,ir)
          stop 5727
        end if
      end do
    end do
  end do
#endif

  ! Reset matrix and RHS
 do idiff=1,FRSIZE_NDIFF
     do it=1,FRSIZE_NDIFF
         dmatrix(idiff,it) = 0.d0
     enddo
     drhs(idiff) = 0.d0
     ipde(idiff) = 0
 enddo

  ! Start counting diffusion cells
  ndiff = 0

   do ir=1+1,nr-1
          do it=1+1,nt    ! In case of theta: go all the way to equator
              if(iphotcount(it,ir).lt.nphotdiff) then
                  gotit = .false.
                  do idiff=1,ndiff
                      if((irdiff(idiff).eq.ir).and. &
                        (itdiff(idiff).eq.it)) then
                          gotit=.true.
                          ipde(idiff) = 1
                          ndiff0 = idiff
                      endif
                  enddo
                  if(.not.gotit) then
                      ndiff = ndiff + 1
                      if(ndiff.gt.FRSIZE_NDIFF) then
                          write(*,*) 'ERROR: DIFFUSION ALGORITHM'
                          write(*,*) '  CRASHED BECAUSE OF TOO FEW '
                          write(*,*) '  diffusion cell places.'
                          write(*,*) ndiff,FRSIZE_NDIFF
                          stop
                      endif
                      irdiff(ndiff) = ir
                      itdiff(ndiff) = it
                      ipde(ndiff)   = 1
                      ndiff0 = ndiff
                  endif
                  gotit = .false.
                  do idiff=1,ndiff
                      if((irdiff(idiff).eq.ir-1).and. &
                        (itdiff(idiff).eq.it)) then
                          gotit=.true.
                          idiff_rleft(ndiff0) = idiff
                          idiff_rright(idiff) = ndiff0
                      endif
                  enddo
                  if(.not.gotit) then
                      ndiff = ndiff + 1
                      idiff = ndiff
                      if(ndiff.gt.FRSIZE_NDIFF) then
                          write(*,*) 'ERROR: DIFFUSION ALGORITHM'
                          write(*,*) '  CRASHED BECAUSE OF TOO FEW '
                          write(*,*) '  diffusion cell places.'
                          write(*,*) ndiff,FRSIZE_NDIFF
                          stop
                      endif
                      irdiff(ndiff) = ir-1
                      itdiff(ndiff) = it
                      dmatrix(ndiff,ndiff) = 1.d0
                      drhs(ndiff) = tdust(1,1,it,ir-1)**4
                      idiff_rleft(ndiff0) = idiff
                      idiff_rright(idiff) = ndiff0
                  endif
                  gotit = .false.
                  do idiff=1,ndiff
                      if((irdiff(idiff).eq.ir+1).and. &
                        (itdiff(idiff).eq.it)) then
                          gotit=.true.
                          idiff_rright(ndiff0) = idiff
                          idiff_rleft(idiff)   = ndiff0
                      endif
                  enddo
                  if(.not.gotit) then
                      ndiff = ndiff + 1
                      idiff = ndiff
                      if(ndiff.gt.FRSIZE_NDIFF) then
                          write(*,*) 'ERROR: DIFFUSION ALGORITHM'
                          write(*,*) '  CRASHED BECAUSE OF TOO FEW '
                          write(*,*) '  diffusion cell places.'
                          write(*,*) ndiff,FRSIZE_NDIFF
                          stop
                      endif
                      irdiff(ndiff) = ir+1
                      itdiff(ndiff) = it
                      dmatrix(ndiff,ndiff) = 1.d0
                      drhs(ndiff) = tdust(1,1,it,ir+1)**4
                      idiff_rright(ndiff0) = idiff
                      idiff_rleft(idiff)   = ndiff0
                  endif
                  gotit = .false.
                  do idiff=1,ndiff
                      if((itdiff(idiff).eq.it-1).and. &
                        (irdiff(idiff).eq.ir)) then
                          gotit=.true.
                          idiff_tleft(ndiff0) = idiff
                          idiff_tright(idiff) = ndiff0
                      endif
                  enddo
                  if(.not.gotit) then
                      ndiff = ndiff + 1
                      idiff = ndiff
                      if(ndiff.gt.FRSIZE_NDIFF) then
                          write(*,*) 'ERROR: DIFFUSION ALGORITHM'
                          write(*,*) '  CRASHED BECAUSE OF TOO FEW '
                          write(*,*) '  diffusion cell places.'
                          write(*,*) ndiff,FRSIZE_NDIFF
                          stop
                      endif
                      irdiff(ndiff) = ir
                      itdiff(ndiff) = it-1
                      dmatrix(ndiff,ndiff) = 1.d0
                      drhs(ndiff) = tdust(1,1,it-1,ir)**4
                      idiff_tleft(ndiff0) = idiff
                      idiff_tright(idiff) = ndiff0
                  endif
                  if(it.lt.nt) then 
                      gotit = .false.
                      do idiff=1,ndiff
                          if((itdiff(idiff).eq.it+1).and. &
                            (irdiff(idiff).eq.ir)) then 
                              gotit=.true.
                              idiff_tright(ndiff0) = idiff
                              idiff_tleft(idiff)   = ndiff0
                          endif
                      enddo
                      if(.not.gotit) then
                          ndiff = ndiff + 1
                          idiff = ndiff
                          if(ndiff.gt.FRSIZE_NDIFF) then
                              write(*,*) 'ERROR: DIFFUSION ALGORITHM'
                              write(*,*) '  CRASHED BECAUSE OF TOO FEW'
                              write(*,*) '  diffusion cell places.'
                              write(*,*) ndiff,FRSIZE_NDIFF
                              stop
                          endif
                          irdiff(ndiff) = ir
                          itdiff(ndiff) = it+1
                          dmatrix(ndiff,ndiff) = 1.d0
                          drhs(ndiff) = tdust(1,1,it+1,ir)**4
                          idiff_tright(ndiff0) = idiff
                          idiff_tleft(idiff)   = ndiff0
                      endif
                      ipde(ndiff0) = 1
                  endif
              endif
          enddo
      enddo

  if (ndiff == 0) return

  write(*,*) 'Now fixing low-statistics cells near midplane with a diffusion method'



  do idiff = 1, ndiff
    ir = irdiff(idiff)
    it = itdiff(idiff)
    diffconst(idiff) = 1.0d0 / rossmeanalp(tguess, nf, freq, &
                          alpha_a(1, it, ir), alpha_s(1, it, ir))
  end do


do iter = 1, niter
    write(*,*) 'Diffusion iteration ', iter

    do idiff = 1, ndiff
        if (ipde(idiff) /= 1) cycle

        ir = irdiff(idiff)
        it = itdiff(idiff)

        ! R-direction
        dxc = 0.5d0 * (rc(ir+1) - rc(ir-1))
        dxp = rc(ir+1) - rc(ir)
        dxm = rc(ir)   - rc(ir-1)
        dcp = 0.5d0 * (diffconst(idiff) + diffconst(idiff_rright(idiff)))
        dcm = 0.5d0 * (diffconst(idiff) + diffconst(idiff_rleft(idiff)))

        dmatrix(idiff, idiff_rleft(idiff)) = (1.0d0 / dxc) * (dcm / dxm)
        dmatrix(idiff, idiff_rright(idiff)) = (1.0d0 / dxc) * (dcp / dxp)
        dmatrix(idiff, idiff) = - (1.0d0 / dxc) * (dcm / dxm + dcp / dxp)

        ! Theta-direction
        if (it < nt) then
            dxc = 0.5d0 * rc(ir) * (tc(it+1) - tc(it-1))
            dxp = rc(ir) * (tc(it+1) - tc(it))
            dxm = rc(ir) * (tc(it) - tc(it-1))
            dcp = 0.5d0 * (diffconst(idiff) + diffconst(idiff_tright(idiff)))
            dcm = 0.5d0 * (diffconst(idiff) + diffconst(idiff_tleft(idiff)))

            dmatrix(idiff, idiff_tleft(idiff))  = (1.0d0 / dxc) * (dcm / dxm)
            dmatrix(idiff, idiff_tright(idiff)) = (1.0d0 / dxc) * (dcp / dxp)
            dmatrix(idiff, idiff) = dmatrix(idiff, idiff) - (1.0d0 / dxc) * (dcm / dxm + dcp / dxp)
        else
            ! Equator: special ghost cell handling
            thghost = 0.5d0 * pi + abs(0.5d0 * pi - tc(it))
            dxc = 0.5d0 * rc(ir) * (thghost - tc(it-1))
            dxp = 1.0d99
            dxm = rc(ir) * (tc(it) - tc(it-1))
            dcp = 1.0d99
            dcm = 0.5d0 * (diffconst(idiff) + diffconst(idiff_tleft(idiff)))

            dmatrix(idiff, idiff_tleft(idiff)) = (1.0d0 / dxc) * (dcm / dxm)
            dmatrix(idiff, idiff) = dmatrix(idiff, idiff) - (1.0d0 / dxc) * (dcm / dxm)
        end if

        ! Reset RHS
        drhs(idiff) = 0.0d0
    end do



! Initialize
dummy = 1.0d0
dmat(1:ndiff, 1:ndiff) = dmatrix(1:ndiff, 1:ndiff)
dsolold(1:ndiff) = dsol(1:ndiff)
dsol(1:ndiff)    = drhs(1:ndiff)

! LU Decomposition & Solution
call ludcmp(dmat, ndiff, FRSIZE_NDIFF, indx, dummy)
call lubksb(dmat, ndiff, FRSIZE_NDIFF, indx, dsol)

! Check for convergence
error = 0.0d0
do idiff = 1, ndiff
    dummy = abs((dsol(idiff) / (dsolold(idiff) + 1.0d-30)) - 1.0d0)
    if (dummy > error) error = dummy
end do
write(*,*) 'Nr of cells = ', ndiff
write(*,*) 'The error   = ', error

! Exit if converged
if (error < errtol) then
    write(*,*) 'Convergence in diffusion reached in ', iter, ' iterations.'
    exit
end if

! Check for invalid temperatures and correct if needed
solmax = maxval(dsol(1:ndiff))
do idiff = 1, ndiff
    ir = irdiff(idiff)
    it = itdiff(idiff)
    if (dsol(idiff) <= 0.0d0) then
        write(*,*) 'In diffusion: zero or negative temperature detected'
        write(*,*) 'Correcting cell:', idiff, ir, it, dsol(idiff)
        dsol(idiff) = 1.0d-2 * solmax
        ierror = 1
    end if
    tdust(1,1,it,ir) = dsol(idiff)**0.25d0
end do

! Update Rosseland mean inverse opacity
do idiff = 1, ndiff
    ir = irdiff(idiff)
    it = itdiff(idiff)
    diffconst(idiff) = 1.0d0 / rossmeanalp(tdust(1,1,it,ir), nf, freq, &
                          alpha_a(1,it,ir), alpha_s(1,it,ir))
end do

end do

do idiff = 1, ndiff
    ir = irdiff(idiff)
    it = itdiff(idiff)
    dummy = dsol(idiff)**0.25d0
    do ispec = 1, dust_nr_species
        do isize = 1, dust_nr_size(ispec)
            tdust(isize, ispec, it, ir) = dummy
        end do
    end do
end do


end subroutine diffusion_main


  subroutine add_neighbor(irn, itn, parent_idiff, ndiff, tdust, gotit, link_from, link_to)
   ! use diffusion_array
    implicit none
    integer, intent(in) :: irn, itn, parent_idiff
    integer, intent(inout) :: ndiff
    real(8), intent(in) :: tdust(DUST_SIZE_MAX, DUST_SPECIES_MAX, FRSIZE_Y_SMALL, FRSIZE_X)
    logical, intent(out) :: gotit
    integer, intent(inout) :: link_from(:), link_to(:)
    integer :: idiff, ierror

    gotit = .false.
    do idiff = 1, ndiff
      if (irdiff(idiff) == irn .and. itdiff(idiff) == itn) then
        gotit = .true.
        link_from(parent_idiff) = idiff
        link_to(idiff)          = parent_idiff
        exit
      end if
    end do

    if (.not. gotit) then
      ndiff = ndiff + 1
      if (ndiff > FRSIZE_NDIFF) then
        write(*,*) 'ERROR: Too few diffusion cell places.', ndiff, FRSIZE_NDIFF
        ierror = 1
        return
      end if
      idiff = ndiff
      irdiff(idiff) = irn
      itdiff(idiff) = itn
      dmatrix(idiff,idiff) = 1.0d0
      drhs(idiff) = tdust(1,1,itn,irn)**4
      link_from(parent_idiff) = idiff
      link_to(idiff)          = parent_idiff
    end if

  end subroutine add_neighbor


subroutine add_neighbor_theta(irn, itn, parent_idiff, ndiff, tdust, gotit, link_from, link_to)
 ! use diffusion_array
  implicit none 
  integer, intent(in) :: irn, itn, parent_idiff
  integer, intent(inout) :: ndiff
  real(8), intent(in) :: tdust(DUST_SIZE_MAX, DUST_SPECIES_MAX, FRSIZE_Y_SMALL, FRSIZE_X)
  logical, intent(out) :: gotit
  integer, intent(inout) :: link_from(:), link_to(:)
  integer :: idiff, ierror

  gotit = .false.
  do idiff = 1, ndiff
    if (irdiff(idiff) == irn .and. itdiff(idiff) == itn) then
      gotit = .true.
      link_from(parent_idiff) = idiff
      link_to(idiff)          = parent_idiff
      exit
    end if
  end do

  if (.not. gotit) then
    ndiff = ndiff + 1
    if (ndiff > FRSIZE_NDIFF) then
      write(*,*) 'ERROR: Too few diffusion cell places.', ndiff, FRSIZE_NDIFF
      ierror = 1
      return
    end if
    idiff = ndiff
    irdiff(idiff) = irn
    itdiff(idiff) = itn
    dmatrix(idiff, idiff) = 1.0d0
    drhs(idiff) = tdust(1,1,itn,irn)**4
    link_from(parent_idiff) = idiff
    link_to(idiff)          = parent_idiff
  end if
end subroutine add_neighbor_theta



function rossmeanalp(tdust, nf, freq, alpha_a, alpha_s) result(rmeanalp)
    implicit none
    ! Arguments
    integer, intent(in)             :: nf
    real(8), intent(in)             :: tdust
    real(8), intent(in)             :: freq(FRSIZE_FREQ)
    real(8), intent(in)             :: alpha_a(FRSIZE_FREQ)
    real(8), intent(in)             :: alpha_s(FRSIZE_FREQ)
    ! Result
    real(8)                         :: rmeanalp
    ! Local variables
    real(8)                         :: dum1, dum2
    real(8), dimension(size(freq))  :: dumarr1, dumarr2
    integer                         :: inu

    ! External functions
    !real(8) :: bplanckdt, integrate

    ! Compute Planck function and ratios
    do inu = 1, nf
        dumarr1(inu) = bplanckdt(tdust, freq(inu))
        dumarr2(inu) = dumarr1(inu) / (alpha_a(inu) + alpha_s(inu))
    end do

    ! Integrate
    dum1 = integrate(nf, freq, dumarr1)
    dum2 = integrate(nf, freq, dumarr2)

    ! Rosseland mean opacity
    rmeanalp = dum1 / dum2

end function rossmeanalp




function bplanckdt(temp, nu) result(bdt)
    implicit none
    ! Arguments
    real(8), intent(in) :: temp, nu
    ! Result
    real(8)             :: bdt
    ! Local variables
    real(8)             :: theexp

    ! Compute exponent
    theexp = exp(4.7989d-11 * nu / temp)

    ! Compute Planck derivative
    if (theexp < 1.d33) then
        bdt = 7.07661334104d-58 * nu**4 * theexp / &
              ( (theexp - 1.d0)**2 * temp**2 ) + 1.d-290
    else
        bdt = 7.07661334104d-58 * nu**4 / &
              ( theexp * temp**2 ) + 1.d-290
    end if

end function bplanckdt




function integrate(n, x, f) result(intg)
    implicit none
    ! Arguments
    integer, intent(in)        :: n
    real(8), intent(in)        :: x(n), f(n)
    ! Result
    real(8)                    :: intg
    ! Local variable
    integer                    :: i
    real(8)                    :: integral

    ! Initialize integral
    integral = 0.d0

    ! Trapezoidal integration
    do i = 2, n
        integral = integral + 0.5d0 * (f(i) + f(i-1)) * (x(i) - x(i-1))
    end do

    ! Adjust sign if necessary
    if (x(n) > x(1)) then
        intg = integral
    else
        intg = -integral
    end if

end function integrate




end module diffusion


































