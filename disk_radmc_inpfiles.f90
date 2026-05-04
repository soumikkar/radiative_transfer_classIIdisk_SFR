module disk_radmc_inpfiles
use disk_radmc_subroutines
    implicit none

    
    ! [cite_start]The derived type for opacity data must be defined before any interfaces that use it. [cite: 152]
   
   ! type opacity_data_type [cite: 153]
    !    integer :: nf, ns, nrt
   
   ! real(kind=8), allocatable, dimension(:) :: freq [cite: 154]
   !     real(kind=8), allocatable, dimension(:) :: wave
   
   ! real(kind=8), allocatable, dimension(:,:,:) :: cabs, csca [cite: 155]
    
   ! real(kind=8), allocatable, dimension(:) :: trange [cite: 156]
    !end type opacity_data_type

    
    ! [cite_start]Interfaces for external functions (assumed to be in other modules). [cite: 157]
  !  interface
  !
   ! function readopac(nr) result(opac) [cite: 158]
   !         integer, intent(in), optional :: nr
    !        type(opacity_data_type) :: opac !
  ! end function readopac [cite: 159]
  !  end interface
    
  !  interface
  ! function find_peak_starspec(nu, nulnu) result(inu) [cite: 160]
  !          real(kind=8), intent(in), dimension(:) :: nu, nulnu
  !          integer :: inu
   ! end function find_peak_starspec [cite: 161]
  !  end interface

  !  interface
   !     function bplanck(nnu, TT) result(bpl)
   ! real(kind=8), intent(in), dimension(:) :: nnu [cite: 162]
   !         real(kind=8), intent(in) :: TT
   ! real(kind=8), allocatable, dimension(:) :: bpl [cite: 163]
 !       end function bplanck
 ! end interface [cite: 164]

contains

    subroutine write_radmc_density(r, theta, rhodust)
        !
        ! [cite_start]Writes the radius, theta, and dust density files for RADMC. [cite: 165]
        real(kind=8), intent(in), dimension(:) :: r, theta 
        real(kind=8), intent(in), dimension(:,:,:) :: rhodust
        
        integer :: nr, nt, nspec, ir, it, ispec
        
        nr = size(r)
        nt = size(theta)
        nspec = size(rhodust, dim=3)
        
        open(unit=1, file='radius.inp', status='replace')
   
        write(1, '(I5)') nr 
        write(1, *) ''
        do ir = 1, nr
            write(1, '(E13.6)') r(ir)
        end do
        close(1)

        open(unit=2, file='theta.inp', status='replace')
        write(2, '(I5, I5)') nt, 1
        write(2, *) ''
        
       do it = 1, nt 
            write(2, '(E13.6)') theta(it)
        end do
        close(2)

        open(unit=3, file='dustdens.inp', status='replace')
        write(3, '(I5, I5, I5, I5)') nspec, nr, nt, 1
        write(3, *) ''
        do ispec = 1, nspec
            do ir = 1, nr
 
                do it = 1, nt 
                    write(3, '(E13.6)') rhodust(ir, it, ispec) + 1.0d-90
                end do
            end do
            write(3, *) ''
        end do
  
        close(3) 

    end subroutine write_radmc_density


    subroutine write_star(mstar, rstar, tstar, kurucz, kurdir, &
                          inupk, lampk)
        !
        ! [cite_start]Writes the star info and stellar spectrum files. [cite: 171]

        ! Natural constants (from problem_natconst.pro)
        real(kind=8), parameter :: pc = 3.08572d18
        real(kind=8), parameter :: ls = 3.8525d33
        real(kind=8), parameter :: rs = 6.96d10
        real(kind=8), parameter :: ts = 5780.0d0
        real(kind=8), parameter :: pi_val = 3.14159265359d0
        real(kind=8), parameter :: cc = 2.9979d10
        
  
        real(kind=8), intent(in) :: mstar, rstar, tstar 
        logical, intent(in) :: kurucz
        character(len=*), intent(in) :: kurdir
        integer, intent(out) :: inupk
        real(kind=8), intent(out) :: lampk
        
        integer :: nf, nnf, i
        real(kind=8), allocatable, dimension(:,:) :: data
        real(kind=8), allocatable, dimension(:) :: nu, &
 nulnu, temp_freq_arr 
        real(kind=8), allocatable, dimension(:) :: bpl_arr
        real(kind=8) :: lstar, planck_val
        type(opacity_data_type) :: o
       
         
        o = readopac(nr = 1)
        
        nf = size(o%freq)
         
        
        ! [cite_start]Create stellar spectrum (Kurucz or Blackbody) [cite: 174]
        if (kurucz) then
            lstar = ls * (rstar/rs)**2 * (tstar/ts)**4
            
            ! [cite_start]This is a placeholder for a Kurucz model routine [cite: 175]
            write(*,*) 'Kurucz model not implemented, using a placeholder.'
        else 
            open(unit=4, file='starspectrum.inp', status='replace')
            write(4, '(I5)') nf
            write(4, *) ''
            do i = 1, nf
               
                ! [cite_start]Corrected syntax for calling bplanck and getting the value [cite: 177]
                allocate(temp_freq_arr(1))
                temp_freq_arr(1) = o%freq(i)
                bpl_arr = bplanck(temp_freq_arr, tstar)
                planck_val = bpl_arr(1)
                deallocate(temp_freq_arr, bpl_arr)
    
                write(4, '(E13.6, 1X, E13.6)') o%freq(i), pi_val * (rstar**2/pc**2) * planck_val 
            end do
            close(4)
        end if
        
      
        ! [cite_start]Make starinfo.inp [cite: 179]
        open(unit=5, file='starinfo.inp', status='replace')
        write(5, '(I5)') 1
        write(5, '(E13.6)') rstar
        write(5, '(E13.6)') mstar
        write(5, '(E13.6)') tstar
        close(5)
        
    
        ! [cite_start]Determine peak of spectrum [cite: 180]
        open(unit=6, file='starspectrum.inp', status='old')
        read(6, *) nnf
        if (nnf /= nf) stop 'Error: Frequency grid mismatch in starspectrum.inp'
        
        allocate(data(2, nf))
        do i = 1, nf
            read(6, *) data(1, i), data(2, i)
        end do
    
        close(6) 

        allocate(nu(nf), nulnu(nf))
        nu = data(1, :)
        nulnu = nu * data(2, :) * 4.0d0 * pi_val * pc**2
        deallocate(data)

        inupk = find_peak_starspec(nu, nulnu)
        lampk = 1.0d4*cc/nu(inupk)
        deallocate(nu, nulnu)

    end subroutine write_star


    subroutine write_chopdens(tauchop, lmbchop, idxchop)
      
        ! [cite_start]Writes the chopdens.inp file. [cite: 182]

        real(kind=8), intent(in) :: tauchop, lmbchop, idxchop
        
        open(unit=7, file='chopdens.inp', status='replace')
        write(7, '(A, E13.6)') 'taumax  = ', tauchop
        write(7, '(A, E13.6)') 'lambda  = ', lmbchop
        write(7, '(A, E13.6)') 'smooth  = ', idxchop
        close(7)

    end subroutine write_chopdens


    subroutine write_radmc(nphot, iseed, ifast, enthres, cntdump, &
 
                           iquant, ifinstar, npdiff, nvstr, ivstrt, vserrtol)
        !
        ! [cite_start]Writes the main radmc.inp file. [cite: 184]
        
        integer, intent(in) :: nphot, iseed, ifast, cntdump, iquant, ifinstar
        integer, intent(in) :: npdiff, nvstr, ivstrt
        real(kind=8), intent(in) :: enthres, vserrtol

        open(unit=8, file='radmc.inp', status='replace')
        write(8, '(A, I8)')    'nphot       = ', nphot        
        write(8, '(A, I5)')    'iseed       &
 = ', iseed        
        write(8, '(A, I5)')    'imethod     = ', 2            
        write(8, '(A, I5)')    'ifast       = ', ifast        
        write(8, '(A, E13.6)') 'enthres     = ', enthres      
  
        write(8, '(A, I5)')    'cntdump     = ', cntdump     
        write(8, '(A, I5)')    'irestart    = ', 0            
        write(8, '(A, I5)')    &
 'itempdecoup = ', 1            
        write(8, '(A, I5)')    'iquantum    = ', iquant       
        write(8, '(A, I5)')    'istarsurf   = ', ifinstar     
        write(8, '(A, I5)')    'nphotdiff   = ', npdiff       
        write(8, '(A, E13.6)') 'errtol      = ', 1.0d-10      
        write(8, '(A, I5)')  &
   'nvstr       = ', nvstr       
        write(8, '(A, E13.6)') 'vserrtol    = ', vserrtol     
        write(8, '(A, I5)')    'ivstrt      = ', ivstrt       
        write(8, '(A, I5)')    'ntemp       = ', 3000       
        write(8, '(A, E13.6)') 'temp0       = ', 1.0d-2       
        write(8, '(A, E13.6)') 'temp1       = ', 1.0d6        
        close(8)

    end subroutine write_radmc


    subroutine write_raytrace(nextra, nrref)
        !
        ! [cite_start]Writes the raytrace.inp file. [cite: 190]

        integer, intent(in) :: nextra, nrref
        
        open(unit=9, file='raytrace.inp', status='replace')
        write(9, '(A, I5)')    'nrphiinf    = ', 32
        write(9, '(A, I5)')    'nrrayextra  = ', -nextra
        
        if (nrref > 0) then
            write(9, '(A, I5)') 'imethod &
     = ', 1 
            write(9, '(A, I5)') 'nrref       = ', nrref
        else
            write(9, '(A, I5)') 'imethod     = ', 0
        end if
        
        write(9, '(A, I5)')    'dbdr        = &
 ', 1 
        write(9, '(A, E13.6)') 'inclination = ', 45.0d0
        close(9)

    end subroutine write_raytrace


    subroutine write_vstruct(dostruct)
        !
        ! [cite_start]Writes the vstruct.inp file. [cite: 193]

        integer, intent(in), dimension(:) :: dostruct
        integer :: nspec, i
        
        nspec = size(dostruct)

        open(unit=10, file='vstruct.inp', status='replace')
        write(10, '(I5)') 1
        write(10, '(I5)') nspec
        do i = 1, nspec
            write(10, '(I5)') dostruct(i)
       
        end do 
        close(10)

    end subroutine write_vstruct
    
end module disk_radmc_inpfiles
