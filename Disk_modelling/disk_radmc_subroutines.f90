module disk_radmc_subroutines
    implicit none

    type opacity_data_type
        integer :: nf, ns, nrt
        real(kind=8), allocatable, dimension(:) :: freq
        real(kind=8), allocatable, dimension(:) :: wave
        real(kind=8), allocatable, dimension(:,:,:) :: cabs, csca
        real(kind=8), allocatable, dimension(:) :: trange
    end type opacity_data_type

    !A derived type to represent the disk data structure
    type disk_data_type
        real(kind=8), allocatable, dimension(:,:,:) :: rho
        real(kind=8), allocatable, dimension(:) :: r
        real(kind=8), allocatable, dimension(:) :: theta
    end type disk_data_type


    type optical_depth_type
        real(kind=8), allocatable, dimension(:,:) :: taur
        real(kind=8), allocatable, dimension(:,:) :: taut
    end type optical_depth_type


   
contains

    function find_peak_starspec(nu, nulnu) result(inu)
   
     !Finds the index of the peak value in a spectrum
     real(kind=8), intent(in), dimension(:) :: nu, nulnu
        integer :: inu, i, nf
        real(kind=8) :: val
        nf = size(nu)
        inu = -1
        val = -1.0d0
        do i = 1, nf
            if (nulnu(i) > val) then
               
               inu = i 
                val = nulnu(i)
            end if
        end do
        return
    end function find_peak_starspec

    function findchi(kappa, alpha, sigma) result(chi_result)
 
        !Solves for chi using an iterative method
        real(kind=8), intent(in) :: kappa, alpha, sigma
        real(kind=8) :: chi_result, rhs, chi, chiold, dum
        integer :: iloop
        
        rhs = 2.0d0 * alpha / (kappa * sigma)
        chi = 4.0d0
        iloop = 1
        
        if (0.5d0 * sigma &
 * kappa / alpha < 3.0d0) then 
            chi_result = 0.0d0
            return
        end if
        
        do
            chiold = chi
            dum = -2.0d0 * log(rhs * exp(-chi**2 / 2.0d0) / (1.0d0 - erf(chi * sqrt(0.5d0))))
   
            chi = sqrt(abs(dum)) 
            iloop = iloop + 1
            if (iloop > 20) then
                write(*,*) 'No convergence in Hs solver'
                stop
            end if
       
            if (abs(chi - chiold) < 1.0d-4 * abs(chi + chiold)) exit 
        end do
        
        chi_result = chi
    end function findchi

    function bplanck(nnu, TT) result(bpl)
        !
        ! [cite_start]Calculates the Planck function. [cite: 14]
        real(kind=8), intent(in), dimension(:) :: nnu
        real(kind=8), intent(in) :: TT
        real(kind=8), allocatable, dimension(:) :: bpl
        real(kind=8) :: nu_val, T_val, x
        real(kind=8), parameter :: cc=2.9979d10, hh=6.6262d-27, kk=1.3807d-16
        integer :: n, i
        n = size(nnu)
        allocate(bpl(n))
        
       
        if (TT == 0.0d0) then 
            bpl = 0.0d0
            return
        end if

        do i = 1, n
            x = hh * nnu(i) / (kk * TT)
            if (x > 100.0d0) then
            
                bpl(i) = 0.0d0 
            else if (x > 1.0d-3) then
                bpl(i) = (2.0d0*hh*nnu(i)**3/cc**2) / (exp(x)-1.0d0)
            else
                bpl(i) = 2.0d0*nnu(i)**2*kk*TT/cc**2
            end if
        end do
     
        return 
    end function bplanck

    function kappaplanck(opac, temp) result(kappa)
        
        !Calculates the Planck-averaged opacity
        type(opacity_data_type), intent(in) :: opac
        real(kind=8), intent(in) :: temp
        real(kind=8) :: kappa, aa, bb, dnu
        integer :: i, nf
        real(kind=8), allocatable, dimension(:) :: bnu
        
        nf = size(opac%freq)
        bnu = bplanck(opac%freq, temp)
        
        aa &
 = 0.0d0 
        bb = 0.0d0
        do i = 2, nf - 1
            dnu = 0.5d0 * opac%freq(i) * (log(opac%freq(i+1)) - log(opac%freq(i-1)))
            aa = aa + bnu(i) * opac%cabs(i, 1, 1) * dnu ! Assuming a single species
            bb = bb + bnu(i) * dnu
        end do
        kappa = aa / bb
        return
    end function kappaplanck
    
    function solverinselfirr(rin0, Sigma0, rsig, plsigma, kap, Hp0, &
                            fixhs, chi, hsback, chback) result(rin_result)
        !Solves for the inner disk radius with self-irradiation
        real(kind=8), intent(in) :: rin0, Sigma0, rsig, plsigma, kap, Hp0 
        real(kind=8), intent(in), optional :: fixhs, chi
        real(kind=8), intent(out), optional :: hsback, chback
        real(kind=8) :: rin_result, rin, Sig, Hpin, Hsin, rinold, chi_local
        integer :: i
        
        rin = rin0
        do i = 1, 11
          
            if (present(chi)) then 
                chi_local = chi
            else
                Sig = Sigma0 * (rin / rsig)**plsigma
                chi_local = findchi(kap, 1.0d0, Sig)
            end if

           
            if (present(fixhs)) then 
                Hsin = fixhs
                Hpin = Hsin / chi_local
            else
                Hpin = Hp0 * (rin / rin0)**1.5d0
                Hsin = chi_local * Hpin
     
            end if

            rinold = rin
            rin = rin0 * sqrt(1.0d0 + Hsin/rin)
            
            if (abs(rin - rinold) < 1.0d-3 * (rin + rinold)) then
                if (present(hsback)) hsback = Hsin
      
               if (present(chback)) chback = chi_local 
                rin_result = rin
                return
            end if
        end do
        write(*,*) "MAJOR PROBLEM: solve Rin with self-irradiation!"
        stop 
    end function solverinselfirr

    function mass(rhodust, r, t, gastodust) result(total_mass)
        !
        !Calculates the total disk mass from a 2D density grid
        real(kind=8), intent(in), dimension(:,:,:) :: rhodust 
        real(kind=8), intent(in), dimension(:) :: r, t
        real(kind=8), intent(in), optional :: gastodust
        
        real(kind=8) :: total_mass, m_val, vol
        integer :: nr, nt, ir, it, imirt
        real(kind=8), allocatable, dimension(:) :: ri, ti
        real(kind=8) :: gastodust_val
        real(kind=8), parameter :: pi = 3.141592653589793d0

 
       if (present(gastodust)) then 
            gastodust_val = gastodust
        else
            gastodust_val = 100.0d0
        end if
        
        nr = size(r, dim=1)
        nt = size(t, dim=1)
        
        
        if (maxval(t) <= pi/2.0d0) then 
            imirt = 1
        else
            imirt = 0
            write(*,*) 'Warning: mirror theta not active'
        end if
        
        allocate(ri(nr + 1), ti(nt + 1))
        ri(1) = r(1)
   
        ri(nr + 1) = r(nr) 
        ti(1) = 0.0d0
        if (imirt == 1) then
            ti(nt + 1) = pi / 2.0d0
        else
            ti(nt + 1) = pi
        end if
        
        do ir &
 = 2, nr 
            ri(ir) = 0.5d0 * (r(ir) + r(ir-1))
        end do
        do it = 2, nt
            ti(it) = 0.5d0 * (t(it) + t(it-1))
        end do
        
        m_val = 0.0d0
        do ir = 1, nr
  
            do it = 1, nt 
                vol = (2.0d0 * pi / 3.0d0) * (ri(ir+1)**3 - ri(ir)**3) * &
                      abs(cos(ti(it)) - cos(ti(it+1)))
                m_val = m_val + vol * rhodust(ir, it, 1) !Assuming a single species
            end do
        end do
        
        m_val = m_val * gastodust_val
        if (imirt == 1) m_val = m_val * 2.0d0
        total_mass = m_val
    end function mass

    function readopac(nr) result(opac)
        !
        ! Reads a dust opacity file
        integer, intent(in), optional :: nr
        type(opacity_data_type) :: opac
        
        ! Local variables
        character(len=20) :: filename
        integer :: nr_val, nf_val, ns_val, nrt_val
        integer :: i, is, k
        real(kind=8) :: dum, a, b
        integer :: iunit, io_stat
        logical :: file_exists
        character(len=256) :: freq_file
        
        !Determine file number from optional argument
        if (present(nr)) then
            nr_val = nr
        else
            nr_val = 1
        end if
        
        !Construct the filename
        write(filename, '(A, I0, A)') 'dustopac_', nr_val, '.inp'
        write(*,*) 'Reading ', trim(filename)
        
        !Open the opacity file
        iunit = 10
        open(unit=iunit, file=trim(filename), status='old', iostat=io_stat)
        if (io_stat /= 0) then
            write(*,*) 'Error: Could not open file ', trim(filename)
            stop
        end if
        
        !
        !Read format and number of species/frequencies
        read(iunit, *) nf_val, ns_val
        
        !
        !Check the file format
        if (nf_val >= 1) then
            !
            ! Format 1: Simple absorption/scattering cross-section
            opac%nf = nf_val
            opac%ns = ns_val
            opac%nrt = 1
            
            allocate(opac%cabs(nf_val, ns_val, 1), opac%csca(nf_val, ns_val, 1))
            allocate(opac%trange(1))
            
            
            ! Read absorption data.Loop order is now correct
            do k = 1, nf_val
                do is = 1, ns_val
                    read(iunit, *) dum
                    opac%cabs(k, is, 1) = dum
                end do
    
            end do 
            
            !
            ! Read scattering data.Loop order is now correct
            do k = 1, nf_val
                do is = 1, ns_val
                    read(iunit, *) dum
                    opac%csca(k, is, 1) = dum
                end do
     
            end do 
            opac%trange = 0.0d0
            
        else
            !
            !Format 2: More complex, includes temperature ranges
            opac%nf = ns_val
            opac%ns = nf_val
            
            read(iunit, *) nf_val, a, nrt_val
            opac%nrt = nrt_val
            
            !
            !Check for smoothing, which is not supported
            if (nint(a) /= 0) then
                write(*,*) "Smoothing is not supported. Stopping."
                stop 
            end if
            
            allocate(opac%cabs(opac%nf, opac%ns, opac%nrt), &
                     opac%csca(opac%nf, opac%ns, opac%nrt), &
                     opac%trange(opac%nrt + 1))
            
 
            !Read temperature ranges and data
            do i = 1, opac%nrt 
                read(iunit, *) a, opac%trange(i+1)
                
                !
                ! Read absorption data.Loop order is now correct
                do k = 1, opac%nf
                    do is = 1, opac%ns
                        read(iunit, *) dum
                        opac%cabs(k, is, i) = dum
          
                    end do 
                end do
                
                !
                ! Read scattering data.Loop order is now correct
                do k = 1, opac%nf
                    do is = 1, opac%ns
                        read(iunit, *) dum
                        opac%csca(k, is, i) = dum
          
                   end do 
                end do
            end do
            
        end if
        close(iunit)

        !
        ! Read the frequency file
        freq_file = 'frequency.inp'
        inquire(file=trim(freq_file), exist=file_exists)
        if (.not. file_exists) then
            write(*,*) "Could not find frequency.inp. Taking frequency.dat"
            freq_file = 'frequency.dat'
            inquire(file=trim(freq_file), exist=file_exists)
            if (.not. file_exists) then
              
                write(*,*) "Could not find frequency.dat either" 
                stop
            end if
        end if
        
        iunit = 11
        open(unit=iunit, file=trim(freq_file), status='old')
        read(iunit, *) nf_val
        if (nf_val /= opac%nf) then
     
            write(*,*) "ERROR: Frequency file has a different number of points than the dustopac file." 
            stop 
        end if
        
        allocate(opac%freq(nf_val), opac%wave(nf_val))
        
        do k = 1, nf_val
            read(iunit, *) opac%freq(k)
            opac%wave(k) = 2.9979d14 / opac%freq(k)
        end do
        close(iunit)
       
        return 
    end function readopac
    


    function maketau(a, kappa) result(tau_result)
       
        !This function calculates the optical depth in the radial and vertical directions.
        type(disk_data_type), intent(in) :: a 
        real(kind=8), intent(in), dimension(:) :: kappa
        type(optical_depth_type) :: tau_result
        
       
        !Local variables
        integer :: nr, nt, nspec, nspecc
        real(kind=8) :: dummy
        integer :: ir, it, ispec
        
       
        !Get array sizes
        nr = size(a%rho, dim=1)
        nt = size(a%rho, dim=2)
        nspecc = size(a%rho, dim=3)
        nspec = size(kappa)
        
        if (nspec /= nspecc) then
            write(*,*) 'ERROR: Number of opacities must equal the number of dust components.'
            stop
        end if
        
        !
        !Allocate the output arrays
        allocate(tau_result%taur(nr, nt), tau_result%taut(nr, nt))
        
        !
        ! --- Calculate radial optical depth (taur) ---
        do it = 1, nt
            tau_result%taur(1, it) = 0.0d0
            do ir = 2, nr
                dummy = 0.0d0
                do ispec = 1, nspec
             
                    dummy = dummy + 0.5d0 * (a%rho(ir, it, ispec) + a%rho(ir-1, it, ispec)) * &
                            kappa(ispec) * (a%r(ir) - a%r(ir-1)) 
                end do
                tau_result%taur(ir, it) = tau_result%taur(ir-1, it) + dummy
        
            end do 
        end do
        
        !
        !--- Calculate vertical optical depth (taut) --- 
        do ir = 1, nr
            tau_result%taut(ir, 1) = 0.0d0
            do it = 2, nt
                dummy = 0.0d0
                do ispec = 1, nspec
             
                    dummy = dummy + 0.5d0 * (a%rho(ir, it, ispec) + a%rho(ir, it-1, ispec)) * &
                            kappa(ispec) * (a%theta(it) - a%theta(it-1)) * a%r(ir) 
                end do
                tau_result%taut(ir, it) = tau_result%taut(ir, it-1) + dummy
      
            end do 
        end do
        
        return
    end function maketau



    subroutine smooth_rim(r, theta, rhodust, kappa, sigdust, drsm, tautol)
        !
        !This subroutine smoothes the inner rim of a disk if the optical depth is too high. 
        real(kind=8), intent(in), dimension(:) :: r, theta, kappa 
        real(kind=8), intent(inout), dimension(:,:,:) :: rhodust
        real(kind=8), intent(inout), dimension(:), optional :: sigdust
        real(kind=8), intent(in), optional :: drsm, tautol

        !
        !Local variables
        integer :: nr, nt, nspec
        integer :: ispec, ir, i
        real(kind=8) :: tau, drsm_val, tautol_val, rdxsm
        real(kind=8) :: r1, x, q1, q2, z_val, redux
        real(kind=8), allocatable, dimension(:) :: x_arr, z_arr, sigdust_flat
        real(kind=8), parameter :: ten = 10.0d0

        !
        !Get array sizes
        nr = size(r)
        nt = size(theta)
        nspec = size(rhodust, dim=3)
        
        !
        !Set default values for optional arguments
        if (present(drsm)) then
            drsm_val = drsm
        else
            drsm_val = 0.1d0
        end if
        
        if (present(tautol)) then
            tautol_val = tautol
        else
       
            tautol_val = 0.5d0 
        end if

        !
        ! --- Find the optical depth of the first cell at the midplane
        
        tau = 0.0d0
        do ispec = 1, nspec
            tau = tau + 0.5d0 * (rhodust(1, nt, ispec) + rhodust(2, nt, ispec)) * &
                  kappa(ispec) * (r(2) - r(1))
        end do
        
        !
        !--- If this is too high, then smooth the inner rim a bit
        if (tau > tautol_val) then
            !
            !Calculate the required reduction factor
            rdxsm = tautol_val / tau
            
            write(*,*) 'Smoothing inner rim with a factor ', rdxsm
            
            !
            !Make a smooth curve from 0 to 1
            r1 = r(1) * (1.0d0 + drsm_val)
            
            allocate(x_arr(nr), z_arr(nr))
            do i = 1, nr
                x_arr(i) = (r(i) - r(1)) / (r1 - r(1))
          
                end do 
            
            q1 = 2.0d1
            q2 = 1.0d0
            
            !
            !This is a simplified translation of the IDL `alog10`.
            do i = 1, nr
                if (x_arr(i) < 0.0d0) then
                    x_arr(i) = 0.0d0
                end if
                z_arr(i) = -log10(exp(-x_arr(i)*q1) + ten**(-q2)) / q2
      
            end do 
            
            !
            !Ensure `z_arr(1)` is 0.0d0.
            z_val = z_arr(1)
            do i = 1, nr
                z_arr(i) = 1.0d0 - (1.0d0 - z_arr(i)) / (1.0d0 - z_val)
            end do
            
            !
            !Now reduce the density.
            do ir = 1, nr
                x = 1.0d0 - z_arr(ir)
                redux = rdxsm**x
                
                rhodust(ir, :, :) = rhodust(ir, :, :) * redux
              
            
                if (present(sigdust)) then
                    sigdust(ir) = sigdust(ir) * redux
                end if
            end do

            deallocate(x_arr, z_arr)
        end if
  
    end subroutine smooth_rim 



    function integrate(x, f, cons, prim) result(integral_result)
        !
        ! This function performs a numerical integration of a 1D array.
       real(kind=8), intent(in), dimension(:) :: x, f 
        logical, intent(in), optional :: cons, prim
        real(kind=8), allocatable :: integral_result(:)
        
        !
        !Local variables 
        integer :: n, i
        real(kind=8) :: total_integral, sign_val
        real(kind=8), allocatable, dimension(:) :: ff
        logical :: prim_present

        n = size(x)

        if (n /= size(f)) then
            write(*,*) "Error in function integrate: x and f are not equally long."
            stop 
        end if
        
        prim_present = present(prim)
        if (prim_present) then
            allocate(ff(n))
            ff = 0.0d0
        end if
        
        sign_val = sign(1.0d0, x(n) - x(1))
       
      total_integral = 0.0d0
        
        if (present(cons)) then
            !
            !Constant integration method 
            do i = 2, n
                total_integral = total_integral + f(i-1) * (x(i) - x(i-1))
                if (prim_present) then
                    ff(i) = total_integral
                end if
 
            end do 
        else
            !
            ! Trapezoidal integration method (default)
            do i = 2, n
                total_integral = total_integral + 0.5d0 * (f(i) + f(i-1)) * (x(i) - x(i-1))
                if (prim_present) then
                    ff(i) = total_integral
             
            end if 
            end do
        end if
        
        if (.not. prim_present) then
            allocate(integral_result(1))
            integral_result(1) = total_integral * sign_val
        else
            if (x(1) < x(n)) then
   
                integral_result = ff 
            else
                integral_result = ff - ff(n)
            end if
        end if
        
        return
    end function integrate

    
end module disk_radmc_subroutines
