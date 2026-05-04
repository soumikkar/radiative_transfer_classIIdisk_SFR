module disk_params
use disk_global_constants
    
    implicit none

     type r_refine_type
        integer :: nspanr
        integer :: nlevr
        integer :: nstepr
    end type r_refine_type
    
    ! Main type to hold all simulation parameters.
    type problem_params
        ! Main switch
        integer :: imakedisk
        
        ! Simulation names
        character(len=10) :: simser
        character(len=10) :: simnr
        character(len=10) :: simweb
        
        ! Central star
        integer :: istar
        real(kind=8) :: rstar, mstar, tstar
        integer :: kurucz
        integer :: ifinstar
        
        ! Opacity and frequency grid
        real(kind=8) :: abunc
        character(len=256), allocatable, dimension(:) :: mixnames
        character(len=256), allocatable, dimension(:,:) :: mixspecs
        real(kind=8), allocatable, dimension(:,:) :: mixabun
        integer :: fresmd
        character(len=256), allocatable, dimension(:) :: infile
        real(kind=8), allocatable, dimension(:) :: pll
        integer :: scat
        
        ! General parameters
        real(kind=8) :: gastodust
        real(kind=8) :: rhofloor
        integer :: run
        
        ! Spatial grid
        integer :: nr, nt, ntex
        real(kind=8) :: rin, rout, hrgrid, hrgmax
        real(kind=8) :: thintin
        
        ! Grid refinement
        type(r_refine_type) :: rrefine
        real(kind=8) :: drsm
        
        ! Sigma(R) setup
        real(kind=8) :: rdisk, sig0, mdmstr, mdisk
        real(kind=8) :: plsig1, plsig2
        real(kind=8) :: bgdens
        
        ! Radial mixing
        real(kind=8) :: schm
        real(kind=8) :: ab_r0
        real(kind=8), allocatable, dimension(:) :: ab_ab0, ab_pl, ab_min
        
        ! H(R) settings
        real(kind=8) :: hrdisk, hrmin, plh
        
        ! Artificial inner rim
        real(kind=8) :: rpfrin, hrpuff
        
        ! RADMC radiative transfer
        integer :: nphot, npdiff
        real(kind=8) :: errtol
        integer :: ifast
        
        ! RADMC vertical structure
        integer :: nvstr
        real(kind=8) :: vserrt
        integer :: ivstrt
        integer, allocatable, dimension(:) :: dostr
        
        ! Paths (assuming they will be set at runtime)
        character(len=256) :: bindir, kuruczdir, src_dust, dustdata
    end type problem_params

contains

    subroutine get_params(params, disk_id)
        ! This subroutine initializes the problem parameters.
        ! If disk_id is provided and > 0, reads from parameter file
        
        type(problem_params), intent(out) :: params
        integer, intent(in), optional :: disk_id
        
        ! Local variables for reading file
        character(len=256) :: param_file, line, full_path
        integer :: ios, file_unit, i, id_read
        real(kind=8) :: mstar_read, rstar_read, tstar_read, mdisk_read
        real(kind=8) :: rdisk_read, rout_read, rin_read
        logical :: found, file_exists
        
        ! Set default parameters first
        params%imakedisk = 1
        params%simser = '1'
        params%simnr  = '1'
        params%simweb = '1'
        
        ! Default star parameters (will be overwritten if reading from file)
        params%istar = 10
        params%rstar = 2.0d0 * RS
        params%mstar = 2.5d0 * MS
        params%tstar = 10000.0d0
        params%kurucz = 0
        params%ifinstar = 1

        ! Default disk parameters (will be overwritten if reading from file)
        params%rin    = 0.5d0 * AU
        params%rout   = 250.0d0 * AU
        params%rdisk  = 200.0d0 * AU
        params%mdmstr = 1.0d-2
        params%mdisk  = params%mdmstr * MS
        
        ! If disk_id is provided, read from parameter file
        if (present(disk_id)) then
            if (disk_id > 0) then
                param_file = 'disk_parameters.tsv'
                file_unit = 20
                found = .false.
                
                ! Check if file exists
                inquire(file=trim(param_file), exist=file_exists)
                if (.not. file_exists) then
                    write(*,*) "Error: Cannot open parameter file: ", trim(param_file)
                    write(*,*) "Current directory contents:"
                    call system("ls -la disk_parameters.tsv 2>/dev/null || echo 'File not found'")
                    stop
                end if
                
                open(unit=file_unit, file=trim(param_file), status='old', iostat=ios)
                if (ios /= 0) then
                    write(*,*) "Error: Cannot open parameter file: ", trim(param_file)
                    write(*,*) "I/O status:", ios
                    stop
                end if
                
                ! Skip header line (contains column names with quotes)
                read(file_unit, '(A)', iostat=ios) line
                
                ! Read through file to find matching disk_id
                ! Format: disk_id mdisk rdisk rin rstar mstar tstar rout
                do
                    read(file_unit, *, iostat=ios) id_read, mdisk_read, rdisk_read, &
                                                    rin_read, rstar_read, mstar_read, &
                                                    tstar_read, rout_read
                    if (ios /= 0) exit
                    
                    if (id_read == disk_id) then
                        ! Found the matching disk
                        params%mstar = mstar_read * MS
                        params%rstar = rstar_read * RS
                        params%tstar = tstar_read
                        params%mdisk = mdisk_read * MS
                        params%mdmstr = mdisk_read
                        params%rdisk = rdisk_read * AU
                        params%rout = rout_read * AU
                        params%rin = rin_read * AU
                        found = .true.
                        exit
                    end if
                end do
                
                close(file_unit)
                
                if (.not. found) then
                    write(*,*) "Error: disk_id ", disk_id, " not found in parameter file"
                    stop
                end if
                
                ! Update simulation identifiers
                write(params%simnr, '(I0)') disk_id
            end if
        end if

        ! Set opacity parameters
        params%abunc = 0.50d0
        allocate(params%mixnames(2))
        params%mixnames(1) = 'silicate_0.1.Kappa'
        params%mixnames(2) = 'silicate_10.Kappa'
        
        allocate(params%mixspecs(2,2))
        params%mixspecs(1,1) = 'silicate_0.1.Kappa'
        params%mixspecs(1,2) = 'silicate_0.1.Kappa'
        params%mixspecs(2,1) = 'silicate_10.Kappa'
        params%mixspecs(2,2) = 'silicate_10.Kappa'

        allocate(params%mixabun(2,2))
        params%mixabun(1,1) = params%abunc
        params%mixabun(1,2) = 1.0d0 - params%abunc
        params%mixabun(2,1) = params%abunc
        params%mixabun(2,2) = 1.0d0 - params%abunc

        params%fresmd = 30
        allocate(params%infile(2))
        params%infile(1) = 'silicate_0.1.Kappa'
        params%infile(2) = 'silicate_10.Kappa'
        
        allocate(params%pll(2))
        params%pll(1) = -1.0d0
        params%pll(2) = -1.0d0
        
        params%scat = 1
        params%gastodust = 100.0d0
        params%rhofloor  = 1.0d-26
        params%run       = 0

        params%nr     = 130
        params%nt     = 50
        params%ntex   = 10
        params%hrgrid = 1.5d0
        params%hrgmax = PI / 2.0d0

        params%rrefine%nlevr  = 3
        params%rrefine%nspanr = 3
        params%rrefine%nstepr = 3
        params%drsm  = 0.03d0

        params%sig0   = 0.0d0
        params%plsig1 = -1.0d0
        params%plsig2 = -12.0d0
        params%bgdens = 1.0d-40

        params%schm = 1.0d0/3.0d0
        params%ab_r0 = 1.0d0 * AU
        allocate(params%ab_ab0(1), params%ab_pl(1), params%ab_min(1))
        params%ab_ab0(1) = 1.0d0
        params%ab_pl(1)  = -3.0d0 * params%schm / 2.0d0
        params%ab_min(1) = 1.0d-6

        params%hrdisk = 0.15d0
        params%hrmin  = 0.01d0
        params%plh    = 2.0d0 / 7.0d0

        params%rpfrin = 0.0d0
        params%hrpuff = 0.0d0

        params%nphot  = 10000000
        params%npdiff = 30
        params%errtol = 1.0d-6
        params%ifast  = 1

        params%nvstr = 10
        params%vserrt = 0.01d0
        params%ivstrt = 1
        allocate(params%dostr(2))
        params%dostr = [1, 1]

        params%bindir = ''
        params%kuruczdir = ''
        params%src_dust = ''
        params%dustdata = ''

    end subroutine get_params
end module disk_params
