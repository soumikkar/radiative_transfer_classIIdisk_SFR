! This is the main program to run the disk model setup.
! It uses the modules provided to generate all the necessary input files
! for a radiative transfer code.

program main_program
    
    use disk_global_constants
    use disk_params, only: problem_params, get_params
    use disk_grid, only: make_rgrid, make_tgrid
    use disk_opacity, only: useopac
    use disk_structure, only: simpledisk_vertstruct
    use disk_radmc_inpfiles, only: write_star, write_radmc, write_vstruct
    use disk_mixopacity, only: mixopacities

    implicit none
    
    ! Declare variables for problem parameters.
    type(problem_params) :: params
    
    ! Variables for disk_id input
    integer :: disk_id
    character(len=20) :: arg
    integer :: iostat
    
    ! Other local variables
    integer :: nf, npah, nspec, ntherm, nfiles, nspec_val
    integer :: iunit
    character(len=256), allocatable, dimension(:) :: infile, pll, opac_files
    character(len=256), allocatable, dimension(:,:) :: mixspecs
    character(len=256) :: file_out
    real(kind=8), allocatable, dimension(:) :: abun
    real(kind=8) :: lampk
    integer :: inupk
    
    ! Variables for the main subroutine call
    real(kind=8), allocatable, dimension(:) :: r, theta, sigdust1, hrstore
    real(kind=8), allocatable, dimension(:,:,:) :: rhodust
    integer, allocatable, dimension(:) :: dostr
    
    ! ============================================
    ! Get disk_id from command line or stdin
    ! ============================================
    disk_id = 0  ! Default to 0 (use default parameters)
    
    if (command_argument_count() > 0) then
        ! Read from command line argument
        call get_command_argument(1, arg)
        read(arg, *, iostat=iostat) disk_id
        if (iostat /= 0) then
            write(*,*) "Error: Invalid disk_id argument"
            stop
        end if
    else
        ! Read from stdin
        read(*,*, iostat=iostat) disk_id
        if (iostat /= 0) then
            write(*,*) "Error reading disk_id from stdin"
            write(*,*) "Using default parameters (disk_id = 0)"
            disk_id = 0
        end if
    end if
    
    ! ============================================
    ! Get parameters for this disk
    ! ============================================
    if (disk_id > 0) then
        write(*,*) "=============================================="
        write(*,*) "Loading parameters for disk_id:", disk_id
        call get_params(params, disk_id)
    else
        write(*,*) "=============================================="
        write(*,*) "Using default parameters"
        call get_params(params)
    end if
    
    ! Print loaded parameters for verification
    write(*,*) "=============================================="
    write(*,*) "LOADED PARAMETERS:"
    write(*,*) "=============================================="
    write(*,*) "Stellar parameters:"
    write(*,'(A,F10.4,A)') "  M* = ", params%mstar/MS, " Msun"
    write(*,'(A,F10.4,A)') "  R* = ", params%rstar/RS, " Rsun"
    write(*,'(A,F10.1,A)') "  T* = ", params%tstar, " K"
    write(*,*) "Disk parameters:"
    write(*,'(A,E12.4,A)') "  Mdisk = ", params%mdisk/MS, " Msun"
    write(*,'(A,F10.2,A)') "  Rdisk = ", params%rdisk/AU, " AU"
    write(*,'(A,F10.2,A)') "  Rin   = ", params%rin/AU, " AU"
    write(*,'(A,F10.2,A)') "  Rout  = ", params%rout/AU, " AU"
    write(*,*) "Grid parameters:"
    write(*,'(A,I5)') "  Nr    = ", params%nr
    write(*,'(A,I5)') "  Nt    = ", params%nt
    write(*,*) "=============================================="
    write(*,*)
    
    write(*,*) 'Starting disk model setup...'
    
    ! ============================================
    ! 1. Mix opacities if specified
    ! ============================================
    nfiles = size(params%infile)
    allocate(abun(nfiles), opac_files(nfiles))
    
    abun = 1.0d0 / real(nfiles) ! Assuming equal mixing for simplicity
    opac_files = params%infile

    write(*,*) 'Mixing opacities...'
    file_out = 'silicate_mix.Kappa'
    call mixopacities(opac_files, file_out, abun)

    ! ============================================
    ! 2. Generate frequency grid and opacity files
    ! ============================================
    write(*,*) 'Generating opacity and frequency files...'
    call useopac(params%infile, params%pll, nf=nf, nspec=nspec, &
                npah=npah, ntherm=ntherm, fresmd=params%fresmd)
     
    ! ============================================
    ! 3. Call the main disk structure routine
    ! ============================================
    write(*,*) 'Generating disk structure and density files...'
    
    call simpledisk_vertstruct( &
        rstar_in=params%rstar, tstar_in=params%tstar, mstar_in=params%mstar, &
        ifinstar_in=params%ifinstar, fresmd_in=params%fresmd, &
        scat_in=params%scat, nrr_in=params%nr, ntt_in=params%nt, &
        rin_in=params%rin, rout_in=params%rout, &
        hrgrid_in=params%hrgrid, hrgmax_in=params%hrgmax, ntex_in=params%ntex, &
        rrefine_in=params%rrefine, drsm_in=params%drsm, &
        rdisk_in=params%rdisk, mdisk_in=params%mdisk, &
        plsig1_in=params%plsig1, plsig2_in=params%plsig2, &
        opacnames_in=params%infile, pllongs_in=params%pll, &
        schmidt_in=params%schm, ab_r0_in=params%ab_r0, ab_ab0_in=params%ab_ab0, ab_pl_in=params%ab_pl, &
        gastodust_in=params%gastodust, ab_min_in=params%ab_min, hrdisk_in=params%hrdisk, &
        hrmin_in=params%hrmin, plh_in=params%plh, rpfrin_in=params%rpfrin, &
        hrpuff_in=params%hrpuff, nvstr_in=params%nvstr, ivstrt_in=params%ivstrt, &
        vserrtol_in=params%vserrt, nphot_in=params%nphot, npdiff_in=params%npdiff, &
        errtol_in=params%errtol, rhofloor_in=params%rhofloor, &
        imakedisk_in=params%imakedisk, run_in=params%run, &
        hrstore_in=hrstore, radius_in=r, theta_in=theta, &
        rhodust_in=rhodust, sigdust1_in=sigdust1, &
        kurucz_in=(params%kurucz == 1),  &
        ifast_in=(params%ifast == 1), dostr_in=dostr)
    
    ! ============================================
    ! 4. Generate other RADMC input files
    ! ============================================
    write(*,*) 'Generating star input files...'
    call write_star(params%mstar, params%rstar, params%tstar, &
                    (params%kurucz == 1), inupk, lampk)

    
    write(*,*) 'Generating radmc.inp...'
    call write_radmc(params%nphot, 12345, params%ifast, 1.0d-4, 100, 1, &
                     params%ifinstar, params%npdiff, params%nvstr, &
                     params%ivstrt, params%vserrt)

    write(*,*) 'Generating vstruct.inp...'
    call write_vstruct(params%dostr)

    write(*,*)
    write(*,*) '=============================================='
    write(*,*) 'Disk model setup complete!'
    write(*,*) 'Input files generated successfully.'
    write(*,*) '=============================================='
    
end program main_program
