module disk_mixopacity
use disk_global_constants, only: CC
    implicit none

type kappa_data_type
        real(kind=8), allocatable, dimension(:) :: freq
        real(kind=8), allocatable, dimension(:) :: lambda
        real(kind=8), allocatable, dimension(:) :: kappa_abs
        real(kind=8), allocatable, dimension(:) :: kappa_sca
        real(kind=8), allocatable, dimension(:) :: g_sca
    end type kappa_data_type



contains

    subroutine mixopacities(files, fileo, abun)
        ! This subroutine mixes multiple opacity files according to given abundances.
        
        character(len=*), intent(in), dimension(:) :: files
        character(len=*), intent(in) :: fileo
        real(kind=8), intent(in), dimension(:) :: abun

        ! Local variables
        integer :: nfiles0, nfiles, ifile, i
        real(kind=8), allocatable, dimension(:) :: ab
        integer :: iformat, nf, nf2, idum
        real(kind=8), allocatable, dimension(:,:) :: data1, data2
        real(kind=8), allocatable, dimension(:) :: lambda, kappa_abs, kappa_sca
        logical :: file_exists

        nfiles0 = size(files)
        nfiles = nfiles0
        do ifile = nfiles0, 1, -1
            if (trim(files(ifile)) == '') then
                nfiles = ifile - 1
            end if
        end do
        
        if (nfiles < 2) return
        if (size(abun) < nfiles) stop 'Error: Abun array is too small.'

        ! Normalize abundances
        allocate(ab(nfiles))
        ab = abun(1:nfiles)
        ab = ab / sum(ab)

        ! Read the first file
        inquire(file=trim(files(1)), exist=file_exists)
        if (.not. file_exists) stop 'Error: File not found: ' // trim(files(1))
        open(unit=10, file=trim(files(1)), status='old')
        read(10, *) iformat
        select case (iformat)
            case (1)
                idum = 2
            case (2)
                idum = 3
            case (3)
                idum = 4
            case default
                stop 'Error: Unknown file format.'
        end select
        read(10, *) nf
        allocate(data1(idum, nf))
        read(10, *) data1
        close(10)

        ! Set the wavelength grid and the main kappa file
        allocate(lambda(nf), kappa_abs(nf))
        lambda = data1(1, :)
        kappa_abs = ab(1) * data1(2, :)
        
        if (iformat > 1) then
            allocate(kappa_sca(nf))
            kappa_sca = ab(1) * data1(3, :)
        else
            allocate(kappa_sca(nf))
            kappa_sca = 0.0d0
        end if
        deallocate(data1)

        ! Read the other files
        do ifile = 2, nfiles
            inquire(file=trim(files(ifile)), exist=file_exists)
            if (.not. file_exists) stop 'Error: File not found: ' // trim(files(ifile))
            open(unit=10, file=trim(files(ifile)), status='old')
            read(10, *) iformat
            select case (iformat)
                case (1)
                    idum = 2
                case (2, 3)
                    idum = 3
                case default
                    stop 'Error: Unknown file format.'
            end select
            read(10, *) nf2
            allocate(data2(idum, nf2))
            read(10, *) data2
            close(10)

            kappa_abs = kappa_abs + ab(ifile) * exp(interpol(log(data2(2, :)), log(data2(1, :)), log(lambda)))
            
            if (iformat > 1) then
                kappa_sca = kappa_sca + ab(ifile) * exp(interpol(log(data2(3, :)), log(data2(1, :)), log(lambda)))
            end if
            deallocate(data2)
        end do
        
        ! Write the opacity mixture
        open(unit=10, file=trim(fileo), status='replace')
        write(10, '(I5)') 3
        write(10, '(I5)') nf
        write(10, *) ''
        do i = 1, nf
            write(10, '(4(E13.6,1X))') lambda(i), kappa_abs(i), kappa_sca(i), 0.0d0
        end do
        close(10)
        
        deallocate(ab, lambda, kappa_abs, kappa_sca)

    end subroutine mixopacities



    function readkappa(file) result(kappa_data)
        ! This function reads a kappa opacity file and returns a derived type
        ! with the extracted data.
        
        character(len=*), intent(in) :: file
        type(kappa_data_type) :: kappa_data
        
        ! Local variables
        integer :: iunit, io_stat, iformat, nf, idum, i
        real(kind=8), allocatable, dimension(:,:) :: data
        
        ! Open the file
        iunit = 10
        open(unit=iunit, file=trim(file), status='old', iostat=io_stat)
        if (io_stat /= 0) then
            write(*,*) 'Error opening file: ', trim(file)
            stop
        end if
        
        read(iunit, *) iformat
        read(iunit, *) nf
        
        select case (iformat)
            case (1)
                idum = 2
            case (2)
                idum = 3
            case (3)
                idum = 4
            case default
                write(*,*) 'Error: Unknown file format in ', trim(file)
                stop
        end select
        
        allocate(data(idum, nf))
        read(iunit, *) data
        close(iunit)
        
        allocate(kappa_data%lambda(nf), &
                 kappa_data%kappa_abs(nf), &
                 kappa_data%kappa_sca(nf), &
                 kappa_data%g_sca(nf), &
                 kappa_data%freq(nf))
                 
        ! Extract data from the read array and transpose
        kappa_data%lambda = data(1, :)
        kappa_data%kappa_abs = data(2, :)
        
        select case (iformat)
            case (1)
                kappa_data%kappa_sca = 0.0d0
                kappa_data%g_sca = 0.0d0
            case (2)
                kappa_data%kappa_sca = data(3, :)
                kappa_data%g_sca = 0.0d0
            case (3)
                kappa_data%kappa_sca = data(3, :)
                kappa_data%g_sca = data(4, :)
        end select
        
        deallocate(data)
        
        ! Calculate frequency from lambda
        do i = 1, nf
            kappa_data%freq(i) = CC / kappa_data%lambda(i)
        end do
        
        return
    
    end function readkappa

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



end module disk_mixopacity

