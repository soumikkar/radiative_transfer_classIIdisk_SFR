module disk_opacity
use disk_global_constants, only: cc,mp, AU, RS, MS, LS, pc, kk, GG, ss

    implicit none

    
   

contains

    subroutine useopac(inputfiles, pllongs, fresmd, scat, onlynf, &
                       nf, pnc, pnh, pz, nspec, npah, ntherm)
        ! Arguments
        character(len=*), intent(in), dimension(:) :: inputfiles
        real(kind=8), intent(in), dimension(:) :: pllongs

        ! Keywords / Optional Arguments
        integer, intent(in), optional :: fresmd
        logical, intent(in), optional :: scat, onlynf
        integer, intent(out), optional :: nf, nspec, npah, ntherm
        integer, intent(in), optional, dimension(:) :: pnc, pnh, pz

        ! Local variables
        integer :: fresmd_val
        logical :: scat_val, onlynf_val, l_exists
        !real(kind=8) :: au, cc, rs, ms, ls, pc, ss, kk, mp, gg
        real(kind=8) :: lambda0, lambda1, lambda2, lambda3
        integer :: nfrnew01, nfrnew12, nfrnew23
        real(kind=8) :: freq0, freq1, freq2, freq3
        real(kind=8) :: dlfreq01, dlfreq12, dlfreq23
        real(kind=8), allocatable, dimension(:) :: freqnw01, freqnw12, freqnw23
        real(kind=8), allocatable, dimension(:) :: freq_hardcode
        real(kind=8), allocatable, dimension(:) :: freqnew
        integer :: nf_val, nfiles, ifile, i, ilam, inu, iopac
        character(len=256) :: filename
        character(len=4) :: file_ext
        integer :: iformat, nlam, n_elements_pnc, n_elements_pnh, n_elements_pz
        real(kind=8), allocatable, dimension(:,:) :: data
        real(kind=8), allocatable, dimension(:) :: lam, kabs, ksca, kap_a, kap_s
        integer, allocatable, dimension(:) :: ii
        real(kind=8), allocatable, dimension(:) :: dum
        integer :: count, ntherm_val, nspec_val
        
        ! PAH variables
        integer :: npah_val
        real(kind=8) :: hh
        real(kind=8), allocatable, dimension(:) :: E_ref, C_ref, cross, kappa
        character(len=20) :: pnc_string
        real(kind=8), allocatable :: freq(:), kabs_rev(:),ksca_rev(:)
        real(kind=8), allocatable :: freq_tmp(:), kabs_tmp(:),ksca_tmp(:)
        integer :: n
        
      !  au = 1.496d13
      !  cc = 2.9979d10
      !  rs = 6.96d10
      !  ms = 1.99d33
      !  ls = 3.8525d33
      !  pc = 3.08572d18
      !  ss = 5.6703d-5
      !  kk = 1.3807d-16
      !  mp = 1.6726d-24
      !  gg = 6.672d-8
      !  hh = 6.62607004d-27

        scat_val = .false.
        onlynf_val = .false.
        fresmd_val = 0
        if (present(scat)) scat_val = scat
        if (present(onlynf)) onlynf_val = onlynf
        if (present(fresmd)) fresmd_val = fresmd

        select case (fresmd_val)
            case (0)
                lambda0  = 1.0d4; lambda1  = 50.0d0; lambda2  = 6.0d0; lambda3  = 0.1d0
                nfrnew01 = 24; nfrnew12 = 22; nfrnew23 = 20
            case (1)
                lambda0  = 1.0d4; lambda1  = 50.0d0; lambda2  = 6.0d0; lambda3  = 0.1d0
                nfrnew01 = 20; nfrnew12 = 80; nfrnew23 = 20
            case (2)
                lambda0  = 1.0d4; lambda1  = 30.0d0; lambda2  = 6.0d0; lambda3  = 0.1d0
                nfrnew01 = 20; nfrnew12 = 26; nfrnew23 = 19
            case (3)
                lambda0  = 1.0d4; lambda1  = 30.0d0; lambda2  = 5.5d0; lambda3  = 0.1d0
                nfrnew01 = 20; nfrnew12 = 30; nfrnew23 = 30
            case (4)
                lambda0  = 1.0d4; lambda1  = 30.0d0; lambda2  = 7.0d0; lambda3  = 0.1d0
                nfrnew01 = 20; nfrnew12 = 30; nfrnew23 = 30
                allocate(freq_hardcode(38))
                freq_hardcode = (/4.28e13, 4.54e13, 4.84e13, 5.0e13, 5.17e13, &
                                 5.450723e13, 6.258446e13, 7.185862e13, &
                                 7.89e13, 8.250724e13, 8.8e13, &
                                 9.473368e13, 9.677e13, 1.087719e14, &
                                 1.248904e14, 1.38888e14, 1.646470e14, &
                                 1.8181818e14, 2.170594e14, 2.400000e14, &
                                 2.861575e14, 3.285621e14, 3.772505e14, &
                                 4.331539e14, 4.973413e14, 5.710405e14, &
                                 6.556609e14, 7.528209e14, 8.643787e14, &
                                 9.924678e14, 1.139538e15, 1.308407e15, &
                                 1.502295e15, 1.724914e15, 1.980523e15, &
                                 2.274010e15, 2.610987e15, 2.997899e15/)
            case (30)
                lambda0  = 1.0d4; lambda1  = 40.0d0; lambda2  = 5.5d0; lambda3  = 0.1d0
                nfrnew01 = 20; nfrnew12 = 80; nfrnew23 = 30
            case (31)
                lambda0  = 1.0d4; lambda1  = 75.0d0; lambda2  = 8.0d0; lambda3  = 0.1d0
                nfrnew01 = 20; nfrnew12 = 200; nfrnew23 = 30
            case default
                write(*,*) "Error: Unknown frequency resolution mode."
                stop
        end select
        print*, cc
        freq0 = 1.0d4 * cc / lambda0
        freq1 = 1.0d4 * cc / lambda1
        freq2 = 1.0d4 * cc / lambda2
        freq3 = 1.0d4 * cc / lambda3

        dlfreq01 = log(freq1) - log(freq0)
        dlfreq12 = log(freq2) - log(freq1)
        dlfreq23 = log(freq3) - log(freq2)

        allocate(freqnw01(nfrnew01))
        do i = 1, nfrnew01
            freqnw01(i) = exp((real(i-1, kind=8) / real(nfrnew01, kind=8)) * dlfreq01 + log(freq0))
        end do
        
        allocate(freqnw12(nfrnew12))
        do i = 1, nfrnew12
            freqnw12(i) = exp((real(i-1, kind=8) / real(nfrnew12, kind=8)) * dlfreq12 + log(freq1))
        end do

        allocate(freqnw23(nfrnew23))
        do i = 1, nfrnew23
            freqnw23(i) = exp((real(i-1, kind=8) / (real(nfrnew23, kind=8)-1.0d0)) * dlfreq23 + log(freq2))
        end do

        if (fresmd_val == 4) then
            deallocate(freqnw23)
            allocate(freqnw23(size(freq_hardcode)))
            freqnw23 = freq_hardcode
        end if
        
        nf_val = size(freqnw01) + size(freqnw12) + size(freqnw23)
        allocate(freqnew(nf_val))
        freqnew(1:size(freqnw01)) = freqnw01
        freqnew(size(freqnw01)+1:size(freqnw01)+size(freqnw12)) = freqnw12
        freqnew(size(freqnw01)+size(freqnw12)+1:nf_val) = freqnw23
        
        if (present(nf)) nf = nf_val

        deallocate(freqnw01, freqnw12, freqnw23)
        if (fresmd_val == 4) deallocate(freq_hardcode)

        if (onlynf_val) then
            if (present(pnc)) then
                npah_val = size(pnc)
            else
                npah_val = 0
            end if
            
            nfiles = size(inputfiles)
            if (present(nspec)) nspec = nfiles + npah_val
            if (present(ntherm)) ntherm = nfiles
            if (present(npah)) npah = npah_val
            return
        end if
        
        open(unit=10, file='frequency.inp', status='replace')
        write(10, '(I5)') nf_val
        write(10, *) ''
        do i = 1, nf_val
            write(10, '(E13.6)') freqnew(i)
        end do
        close(unit=10)

        nfiles = size(inputfiles)
        if (size(pllongs) .ne. nfiles) then
            write(*,*) 'ERROR: Nr of elements of inputfiles not equal to'
            write(*,*) '       nr of elements of pllongs.'
            stop
        end if

        do ifile = 1, nfiles
            inquire(file=trim(inputfiles(ifile)), exist=l_exists)
            if (.not. l_exists) then
                write(*,*) 'ERROR: Could not find file ', trim(inputfiles(ifile))
                stop
            end if

            open(unit=11, file=trim(inputfiles(ifile)), status='old')
            read(11, *) iformat
            read(11, *) nlam
            
            select case (iformat)
                case (1)
                    allocate(data(2, nlam))
                case (2)
                    allocate(data(3, nlam))
                case (3)
                    allocate(data(4,nlam))
                case default
                    write(*,*) 'ERROR: Could not interpret iformat = ', iformat
                    stop
            end select
            read(11, *) data
            close(unit=11)

            select case (iformat)
                case (1)
                    allocate(lam(nlam), kabs(nlam), ksca(nlam))
                    lam = data(1, :)
                    kabs = data(2, :)
                    ksca = 0.0d0
                case (2)
                    allocate(lam(nlam), kabs(nlam), ksca(nlam))
                    lam = data(1, :)
                    kabs = data(2, :)
                    ksca = data(3, :)
                case (3)
                    allocate(lam(nlam), kabs(nlam), ksca(nlam))
                    lam = data(1, :)
                    kabs = data(2, :)
                    ksca = data(3, :)
                case default
                    write(*,*) 'ERROR: Could not interpret iformat = ', iformat
                    stop
            end select
            deallocate(data)

            do ilam = 2, nlam
                if (lam(ilam) <= lam(ilam-1)) then
                    write(*,*) 'ERROR: Input files must have monotonically'
                    write(*,*) '       increasing wavelength grid'
                    stop
                end if
            end do
            n = size(lam)
            if (allocated(freq)) deallocate(freq)
            if (allocated(kabs_rev)) deallocate(kabs_rev)
            if (allocated(ksca_rev)) deallocate(ksca_rev)
            allocate(freq(n), kabs_rev(n),ksca_rev(n))

           
           do i = 1, n
              freq(i) = 3.0d14 / lam(i)
           end do

! Use temporaries to reverse arrays so 'freq' used for interpolation is strictly increasing
           
            allocate(freq_tmp(n), kabs_tmp(n),ksca_tmp(n))
            do i = 1, n
                freq_tmp(i) = freq(n - i + 1)
                kabs_tmp(i) = kabs(n - i + 1)
                ksca_tmp(i) = ksca(n - i + 1)
            end do
            freq = freq_tmp
            kabs_rev = kabs_tmp
            ksca_rev = ksca_tmp
            deallocate(freq_tmp, kabs_tmp, ksca_tmp)

            ! Debug prints (comment out later) to diagnose the issue:
             print *, 'DEBUG: file=', trim(inputfiles(ifile)), ' nlam=', n, ' kabs(1)=', kabs(1), ' kabs(nlam)=', kabs(nlam)
             print *, 'DEBUG: freq new range: ', minval(freq), maxval(freq), ' freqnew range: ', minval(freqnew), maxval(freqnew)

            kap_a = interpol(freq, kabs_rev, freqnew)
            kap_s = interpol(freq, ksca_rev, freqnew)   ! reverse ksca on the fly

            ! handle extrapolation at long wavelengths (lower freq)
            if (maxval(lam) < maxval(3.0d14/freqnew)) then
                count = 0
                do i = 1, nf_val
                    if (freqnew(i) <= 3.0d14/maxval(lam)) then
                        count = count + 1
                    end if
                end do

                if (count > 0) then
                    allocate(ii(count))
                    allocate(dum(count))
                    count = 0
                    do i = 1, nf_val
                        if (freqnew(i) <= 3.0d14/maxval(lam)) then
                            count = count + 1
                            ii(count) = i
                        end if
                    end do

                    do i = 1, count
                        dum(i) = (3.0d14/freqnew(ii(i)))**pllongs(ifile)
                    end do
                    ! normalize by the last element (safeguard: check not zero)
                    if (abs(dum(size(dum))) < 1.0d-300) then
                        write(*,*) 'WARNING: dum(last) is zero or extremely small => skipping longwave scaling'
                    else
                        dum = dum / dum(size(dum))
                        dum = dum * kabs(nlam)
                        kap_a(ii(1:count)) = dum
                        kap_s(ii(1:count)) = 0.0d0
                    end if

                    deallocate(dum)
                    deallocate(ii)
                end if
            end if

            ! handle extrapolation at short wavelengths (higher freq)
            if (minval(lam) > minval(3.0d14/freqnew)) then
                count = 0
                do i = 1, nf_val
                    if (freqnew(i) >= 3.0d14/minval(lam)) then
                        count = count + 1
                    end if
                end do
                if (count > 0) then
                    allocate(ii(count))
                    count = 0
                    do i = 1, nf_val
                        if (freqnew(i) >= 3.0d14/minval(lam)) then
                            count = count + 1
                            ii(count) = i
                        end if
                    end do
                    kap_a(ii(1:count)) = kabs(1)
                    kap_s(ii(1:count)) = ksca(1)
                    deallocate(ii)
                end if
            end if

            ! enforce small floor (avoid zeros)
            do i = 1, nf_val
                kap_a(i) = max(kap_a(i), 1.0d-10)
                kap_s(i) = max(kap_s(i), 1.0d-10)
            end do
            
            write(file_ext, '(I4)') ifile
            filename = 'dustopac_' // trim(adjustl(file_ext)) // '.inp'
            open(unit=12, file=trim(filename), status='replace')
            write(12, '(I4, I4)') nf_val, 1
            write(12, *) ''
            do inu = 1, nf_val
                write(12, '(E13.6)') kap_a(inu)
            end do
            write(12, *) ''
            do inu = 1, nf_val
                write(12, '(E13.6)') kap_s(inu)
            end do
            close(unit=12)

            deallocate(lam, kabs, ksca, kap_a, kap_s)
        end do
        ntherm_val = nfiles
        
        
        nspec_val = ntherm_val
        if (present(nspec)) nspec = nspec_val

        open(unit=14, file='dustopac.inp', status='replace')
        write(14, '(A)') '2               Format number of this file'
        write(file_ext, '(I4)') nspec_val
        write(14, '(A)') trim(adjustl(file_ext)) // '               Nr of dust species'
        write(14, '(A)') '============================================================================'
        
        do iopac = 1, ntherm_val
            write(14, '(A)') '-1              Way in which this dust species is read (-1=file)'
            write(14, '(A)') '0               0=Thermal grain'
            write(file_ext, '(I4)') iopac
            write(14, '(A)') trim(adjustl(file_ext)) // '               Extension of name of dustopac_***.inp file'
            write(14, '(A)') '----------------------------------------------------------------------------'
        end do

        if (npah_val > 0) then
            do iopac = ntherm_val + 1, nspec_val
                write(14, '(A)') '-1              Way in which this dust species is read (-1=file)'
                write(14, '(A)') '2               2=Quantum-heated grain'
                write(pnc_string, '(E13.6)') real(pnc(iopac-ntherm_val))
                write(14, '(A)') trim(adjustl(pnc_string)) // '    Number of C-atoms of PAH/smallcarbgrain'
                write(file_ext, '(I4)') iopac
                write(14, '(A)') trim(adjustl(file_ext)) // '            Extension of name of dustopac_***.inp file'
                write(14, '(A)') '----------------------------------------------------------------------------'
            end do
        end if
        close(unit=14)

        deallocate(freqnew)
    print*, nspec_val, ntherm, ntherm_val
    end subroutine useopac


FUNCTION interpol(x, y, x0) RESULT(y0)
    REAL(8), INTENT(IN) :: x(:), y(:), x0(:)
    REAL(8), DIMENSION(SIZE(x0)) :: y0
    INTEGER :: i, j

    DO j = 1, SIZE(x0)
        IF (x0(j) <= x(1)) THEN
            y0(j) = y(1)
        ELSE IF (x0(j) >= x(SIZE(x))) THEN
            y0(j) = y(SIZE(x))
        ELSE
            DO i = 1, SIZE(x) - 1
                IF (x0(j) >= x(i) .AND. x0(j) <= x(i + 1)) THEN
                    y0(j) = y(i) + (y(i + 1) - y(i)) * (x0(j) - x(i)) / (x(i + 1) - x(i))
                    EXIT
                END IF
            END DO
        END IF
    END DO
END FUNCTION interpol


end module disk_opacity

