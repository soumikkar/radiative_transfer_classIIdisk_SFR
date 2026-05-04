module disk_grid
use disk_params
    implicit none

    ! This derived type contains parameters for radial grid refinement.
  !  type r_refine_type
   !     integer :: nspanr
   !     integer :: nlevr
   !     integer :: nstepr
   ! end type r_refine_type

    ! This derived type contains parameters for polar angle grid refinement.
    type z_refine_type
        real(kind=8) :: zrref
        integer      :: nzref
    end type z_refine_type

contains

    function make_rgrid(rin, rout, nr, rrefine) result(r_grid)
        real(kind=8), intent(in) :: rin, rout
        integer, intent(in) :: nr
        type(r_refine_type), intent(in), optional :: rrefine
        
        real(kind=8), allocatable, dimension(:) :: r_grid

        ! Local variables
        integer :: nextra, nnr, i, ilev, ispan, i_val
        real(kind=8) :: drr, fr
        real(kind=8), allocatable, dimension(:) :: rins
        real(kind=8), allocatable, dimension(:) :: r_temp

        ! Check if the optional rrefine argument is present.
        if (.not. present(rrefine)) then
            ! Case 1: No refinement. Create a simple log-linear grid.
            allocate(r_grid(nr))
            do i = 1, nr
                r_grid(i) = rin * (rout/rin)**((real(i-1, kind=8)) / real(nr-1, kind=8))
            end do
        else
            ! Case 2: With refinement.
            nextra = rrefine%nspanr * rrefine%nlevr * (2**rrefine%nstepr - 1)
           ! print*, nextra, nr, rrefine%nspanr, rrefine%nlevr, rrefine%nstepr
            if (nextra > nr - 5) then
                write(*,*) "Sorry not enough nr for this refinement level"
                stop
            end if

            ! Allocate the initial base grid of size (nr - nextra).
            allocate(r_grid(nr - nextra))
            do i = 1, nr - nextra
                r_grid(i) = rin * (rout/rin)**((real(i-1, kind=8)) / real(nr - 1 - nextra, kind=8))
            end do

            ! Refinement loops, inserting new points.
            do ilev = 1, rrefine%nlevr
                do ispan = rrefine%nspanr, 1, -1
                    allocate(rins(2**rrefine%nstepr - 1))
                    fr = (r_grid(ispan+1)/r_grid(ispan))**(0.5d0**rrefine%nstepr)
                    
                    do i = 1, 2**rrefine%nstepr - 1
                        rins(i) = r_grid(ispan) * fr**i
                    end do

                    ! Fortran 90 array concatenation.
                    nnr = size(r_grid)
                    allocate(r_temp(nnr + size(rins)))
                    r_temp(1:ispan) = r_grid(1:ispan)
                    r_temp(ispan+1:ispan+size(rins)) = rins
                    r_temp(ispan+size(rins)+1:size(r_temp)) = r_grid(ispan+1:nnr)
                    
                    deallocate(r_grid)
                    allocate(r_grid(size(r_temp)))
                    r_grid = r_temp
                    
                    deallocate(r_temp, rins)
                end do
            end do
        end if

        ! Check the final grid size
        nnr = size(r_grid)
        if (nnr .ne. nr) then
            write(*,*) "ERROR: Final grid size mismatch. Something went wrong."
            stop
        end if

        ! Check the DeltaR/R ratio.
        drr = (r_grid(nr) - r_grid(nr-1)) / r_grid(nr-1)
        
        if (drr > 0.15d0) then
            write(*,*) "ERROR: Radial grid too coarse..."
            stop
        end if
        
        if (drr < 0.05d0) then
            write(*,*) 'Radial grid too finely spaced... This will take too much'
            write(*,*) 'computational time... Are you sure this is okay? (Type 1)'
            read(*,*) i_val
            if (i_val .ne. 1) then
                stop
            end if
        end if
    end function make_rgrid


    function make_tgrid(hrgrid, nt, hrgmax, ntex, hrlg, zrefine) result(theta)
        ! This function generates a polar angle grid (in radians).
        
        real(kind=8), intent(in) :: hrgrid
        integer, intent(in) :: nt
        real(kind=8), intent(in), optional :: hrgmax, hrlg
        integer, intent(in), optional :: ntex
        type(z_refine_type), intent(in), optional :: zrefine

        real(kind=8), allocatable, dimension(:) :: theta

        ! Local variables
        integer :: ntt, i, iz, old_size
        real(kind=8) :: thmax, thmin, A, B, pi_half
        real(kind=8), allocatable, dimension(:) :: thex, r_temp, new_points, indices

        pi_half = acos(-1.0d0) / 2.0d0

        if (present(ntex)) then
            ntt = nt - ntex
        else
            ntt = nt
        end if

        if (ntt < 4) then
            write(*,*) "Error: ntt must be at least 4."
            stop
        end if

        if (.not. present(hrlg)) then
            thmax = pi_half
            thmin = pi_half - hrgrid
            
            allocate(theta(ntt))
            do i = 1, ntt
                theta(i) = (thmax - thmin) * (real(i-1, kind=8) / (real(ntt-1, kind=8) + 0.5d0)) + thmin
            end do
            
            if (present(hrgmax)) then
                if (.not. present(ntex)) then
                    write(*,*) 'PROBLEM: If you define hrgmax, must also define ntex'
                    stop
                end if
                
                if (hrgmax <= hrgrid) then
                    write(*,*) 'PROBLEM: hrgmax must be larger than hrgrid'
                    stop
                end if
                
                thmax = pi_half - hrgrid
                thmin = pi_half - hrgmax
                
                allocate(thex(ntex))
                do i = 1, ntex
                    thex(i) = (thmax - thmin) * (real(i-1, kind=8) / real(ntex, kind=8)) + thmin
                end do
                
                old_size = size(theta)
                allocate(r_temp(size(thex) + old_size))
                r_temp(1:ntex) = thex
                r_temp(ntex+1:size(r_temp)) = theta
                
                deallocate(theta)
                allocate(theta(size(r_temp)))
                theta = r_temp
                deallocate(r_temp, thex)
            end if
        else
            B = log(hrgrid)
            A = (B - log(hrlg)) / real(nt - 1, kind=8)
            
            allocate(theta(nt))
            do i = 1, nt
                theta(i) = pi_half - exp(-A * real(i-1, kind=8) + B)
            end do
        end if

        if (present(zrefine)) then
            iz = 0
            do i = 1, size(theta)
                if (theta(i) > pi_half - zrefine%zrref) then
                    iz = i
                    exit
                end if
            end do
            if (iz == 0) then
                iz = size(theta) + 1
            end if

            allocate(new_points(zrefine%nzref))
            allocate(indices(zrefine%nzref))
            do i = 1, zrefine%nzref
                indices(i) = real(zrefine%nzref - i + 1, kind=8) - 0.5d0
            end do
            
            do i = 1, zrefine%nzref
                new_points(i) = pi_half - indices(i) * &
                                (pi_half - theta(iz)) / (real(zrefine%nzref, kind=8) + 1.0d0)
            end do
            
            old_size = size(theta)
            allocate(r_temp(iz - 1 + size(new_points)))
            if (iz > 1) then
                r_temp(1:iz-1) = theta(1:iz-1)
            end if
            r_temp(iz:size(r_temp)) = new_points
            
            deallocate(theta)
            allocate(theta(size(r_temp)))
            theta = r_temp
            deallocate(r_temp, new_points, indices)
        end if
        
    end function make_tgrid

end module disk_grid

