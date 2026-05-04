module disk_structure
    
    use disk_global_constants
    use disk_params
    use disk_grid
    use disk_density
    use disk_radmc_subroutines
    use disk_radmc_inpfiles
    use disk_opacity
    use disk_mixopacity
    
    implicit none

contains

    subroutine simpledisk_vertstruct(rstar_in,tstar_in,mstar_in, &
                                     ifinstar_in,fresmd_in, &
  
                                     scat_in,nrr_in,ntt_in,rin_in,tin_in, &
                                     rout_in,hrgrid_in,hrgmax_in,ntex_in, rrefine_in, &
                        
                                     drsm_in,rdisk_in,sigdust0_in,mdisk_in, &
                                     plsig1_in,plsig2_in, &
                      
                                     opacnames_in,pllongs_in, &
            
                                     schmidt_in,ab_r0_in,ab_ab0_in,ab_pl_in, &
                                  
                                     gastodust_in,ab_min_in,hrdisk_in,hrmin_in, &
                                    
                                     plh_in,rpfrin_in,hrpuff_in,nvstr_in, &
                      
                                     ivstrt_in,vserrtol_in,nphot_in, &
                                     npdiff_in,errtol_in,tauchop_in,lmbchop_in, &
                      
                        
                                     idxchop_in,rhofloor_in,pnc_in,pnh_in,pz_in, &
                                     imakedisk_in,run_in,hrstore_in, &
                                     thintin_in,tt_in,radius_in,theta_in,rhodust_in, &
        
                                     sigdust1_in,kurucz_in,kurdir_in,bindir_in, &
                                     ifast_in,dostr_in)
        !
        

        ! All variables are declared here to avoid IMPLICIT typing issues.
        real(kind=8), intent(in), optional :: rstar_in,tstar_in,mstar_in, &
                                                rin_in,tin_in,rout_in,hrgrid_in,hrgmax_in, &
                                              
                                                drsm_in,rdisk_in,sigdust0_in,mdisk_in, &
                                                plsig1_in,plsig2_in, &
  
                                                schmidt_in,ab_r0_in, &
     
                                                gastodust_in,hrdisk_in,hrmin_in, &
                                                plh_in,rpfrin_in,hrpuff_in, &
 
                                                vserrtol_in,errtol_in,tauchop_in,lmbchop_in, &
                                                idxchop_in,rhofloor_in, &
   
                                                tt_in,thintin_in
        integer, intent(in), optional :: ifinstar_in,fresmd_in,scat_in,nrr_in,ntt_in,nvstr_in, &
                                          
                                         ivstrt_in,nphot_in,npdiff_in, &
                                          pnc_in,pnh_in,pz_in,imakedisk_in,run_in, ntex_in
        character(len=*), intent(in), optional :: opacnames_in(:), &
    
                                         kurdir_in,bindir_in
        type(r_refine_type), intent(in), optional :: rrefine_in
        real(kind=8), intent(in), optional, dimension(:) :: ab_ab0_in, pllongs_in, ab_pl_in, ab_min_in

        
        !
        ! Output variables (intent(out) or intent(inout))
        real(kind=8), allocatable, dimension(:), intent(inout) :: radius_in,theta_in
        real(kind=8), allocatable, dimension(:,:,:), intent(inout) :: rhodust_in
        real(kind=8), allocatable, dimension(:), intent(inout) :: hrstore_in,sigdust1_in
        integer, intent(inout), optional :: dostr_in(:)

        !
        ! Local variables
        integer :: nr, nt, nnr, nnt, i, i_val, ispec, idum, nspec, &
                   nf, npah, ntherm, inupk
        real(kind=8) :: rstar,tstar,mstar,ifinstar,fresmd,scat,nrr,ntt, &
                       rin,rout,hrgrid,hrgmax,drsm,rdisk,sigdust0, &
                       mdisk,plsig1,plsig2,schmidt,ab_r0,ab_ab0, &
                       ab_pl,gastodust,ab_min,hrdisk,hrmin,plh,rpfrin, &
                       hrpuff,nvstr,ivstrt,vserrtol,nphot,npdiff,errtol, &
                       tauchop,lmbchop,idxchop,rhofloor,pnc,pnh,pz, &
                       imakedisk,run_flag,thintin,tt, lampk
        real(kind=8) :: pi_val, pi_half, sigdust00, mddum, ddr, tr0
        real(kind=8), allocatable, dimension(:) :: &
 kappa, q, r_val, theta_val, integral_temp
        real(kind=8), allocatable, dimension(:,:,:) :: rhodusttot, abun
        character(len=20), allocatable, dimension(:) :: opacnames, pllongs
        character(len=256) :: kurdir, bindir
        logical, intent(in), optional:: kurucz_in, &
 ifast_in
        type(r_refine_type) :: rrefine, rrefine_val
        type(disk_data_type) :: b_struct
        type(opacity_data_type) :: o
        type(optical_depth_type) :: tt_struct
        
        logical :: l_exists

        !
        ! --- Set defaults ---
        rstar = RS
        tstar = TS
        mstar = MS
        ifinstar = 1
        fresmd = 3
        scat = 1
        nrr = 130
        ntt = 60
        nphot = 10000000
      
        npdiff = 30
        errtol = 1.0d-10
        tauchop = 0.0d0
        lmbchop = 0.55d0
        idxchop = 1.0d0
        imakedisk = 1
        rhofloor = 1.0d-26
        gastodust = 100.0d0
        nspec = 2 ! Default to a single dust species
        
        pi_val = acos(-1.0d0)
        !mdisk_in = 1.0d-2

        !
        ! --- Handle optional arguments ---
        if (present(rstar_in)) rstar = rstar_in
        if (present(tstar_in)) tstar = tstar_in
        if (present(mstar_in)) mstar = mstar_in
        !
        ! ... all other optional arguments.

        if (present(rin_in) .and. present(tin_in)) then
            write(*,*) 'ERROR: Cannot set both tin and rin simultaneously'
            stop
        end if
        
        if (present(sigdust0_in) .and. present(mdisk_in)) then
            write(*,*) 'ERROR: Either specify sigdust0 or mdisk.'
            write(*,*) 'Not both.'
            stop
        end if

        if (present(nvstr_in) .and. present(tauchop_in)) then
            if (nvstr_in > 0 .and. tauchop_in > 0.0d0) then
                write(*,*) 'ERROR: Vertical structure iteration and chopdens are not compatible.'
                stop
            end if
        end if

        if (present(opacnames_in)) then
            allocate(opacnames(size(opacnames_in)))
            opacnames = opacnames_in
            if (present(ab_r0_in) .and. present(ab_ab0_in)) then
                if (size(opacnames) /= size(ab_ab0_in)+1) then
      
                    write(*,*) 'ERROR: Nr of opacities does not match abundances.'
                    stop
                end if
            end if
        else
            write(*,*) 'ERROR: No opacities specified!'
            write(*,*) 'Set opacnames[].'
            stop
        end if

        !
        ! --- Main logic starts here ---
       
        ! Make the R-grid.
        r_val = make_rgrid(rin_in, rout_in, nrr_in, rrefine=rrefine_in)
        
        !
        ! Make the Theta-grid.
        theta_val = make_tgrid(hrgrid_in, ntt_in, hrgmax=hrgmax_in, ntex=ntex_in)

        nnr = size(r_val)
        nnt = size(theta_val)
        print*, nnr, nnt
        if (present(sigdust0_in)) then
            !
            ! Directly from sigdust0
            allocate(rhodusttot(nnr, nnt, nspec))
            call disk_model_1(r_val,theta_val,rdisk_in,sigdust0_in,plsig1_in,plsig2_in, &
                              hrdisk_in,plh_in,rhodust=rhodusttot,sigdust=sigdust1_in,hrstore=hrstore_in)
        else
            !
            ! Compute sigdust0 from mdisk
            sigdust00 = 1.0d0
            allocate(rhodusttot(nnr, nnt, nspec))
            call disk_model_1(r_val,theta_val,rdisk_in,sigdust00,plsig1_in,plsig2_in, &
                              hrdisk_in,plh_in,rhodust=rhodusttot,sigdust=sigdust1_in,hrstore=hrstore_in)
            integral_temp = integrate(r_val, 2.0d0 * pi_val * r_val * sigdust1_in)
            mddum = integral_temp(1) * gastodust_in
 
            sigdust00 = mdisk_in / mddum
            !
            ! Re-run disk model with corrected sigdust0
            call disk_model_1(r_val,theta_val,rdisk_in,sigdust00,plsig1_in,plsig2_in, &
                              hrdisk_in,plh_in,rhodust=rhodusttot,sigdust=sigdust1_in,hrstore=hrstore_in)
        end if
        
        !
        ! Split density into abundances if necessary.
        if (nspec > 1) then
             allocate(abun(nnr,nnt,nspec))
             abun = 1.0d0
        end if

        !print*, rhodusttot(:,:,1)
        ! Now convert this into dustdens.
        if (allocated(abun)) then
            allocate(rhodust_in(nnr, nnt, nspec))
            do ispec = 1, nspec
                rhodust_in(:,:,ispec) = rhodusttot(:,:,1) * abun(:,:,ispec)
            end do
        else
            print*, nnr, nnt
            allocate(rhodust_in(nnr, nnt, 1))
       
            rhodust_in(:,:,1) = rhodusttot(:,:,1)
        end if

        
        
        call write_radmc_density(r_val, theta_val, rhodust_in)

       
    end subroutine simpledisk_vertstruct

end module disk_structure
