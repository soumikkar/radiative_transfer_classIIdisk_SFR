module main
  use configure
  use common_grid
  use common_dust
  use common_boundary
  use common_montecarlo
  use common_montecarloarrays
  use numerical_receipe
  use dust_main
  use common_diffusion
  use common_vstruct
  use numerical_receipe
  use diffusion
  use vertical_structure
  use montecarlo
  implicit none

  

contains


subroutine radmc_main_program

implicit none

      double precision dummy,enthres
      double precision spectrum(FRSIZE_FREQ,FRSIZE_Y_SMALL)
      double precision tdust(DUST_SIZE_MAX,DUST_SPECIES_MAX, FRSIZE_Y_SMALL,FRSIZE_X)
      double precision scatsrc(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
      double precision miquant(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X)
      double precision errtol
!#ifdef SAVE_MEANINT
      double precision, dimension(FRSIZE_FREQ,FRSIZE_Y_SMALL,FRSIZE_X) :: meanint
!#else
     ! double precision :: meanint
!#endif

      integer inu,nphot,iseed,iformat,imethod,irestart,ir,it,idiff
      integer ispec,cntdump,ifast,itmax,ifile,iredo,i,l
      integer ntemp,nvstr,ivstrt,iterstr,iwarn,ierror_diff
      doubleprecision temp0,temp1,vserrtol,error
      logical fex1,fex2,writescat,iqy,writespec, file_exists,ierror
      character*160 string


inquire(file='radmc.inp', exist=file_exists)

  if (.not. file_exists) then
     write(*,*) 'File radmc.inp not found. Using default parameters.'
     ! Set some default values if you want
     nphot = 1000000
     iseed = -1
     imethod = 2
     ifast = 0
     enthres = 0.01d0
     cntdump = 10
     irestart = 0
     itempdecoup = 0
     iquantum = 0
     istarsurf = 0
     nphotdiff = 100000
     errtol = 0.01d0
     nvstr = 3
     vserrtol = 0.01d0
     ivstrt = 0
     ntemp = 100
     temp0 = 10.d0
     temp1 = 3000.d0
     return
  endif

  open(unit=10, file='radmc.inp', status='old', action='read')
  
  call read_input_file(10)  ! your existing subroutine to read lines into buffer

  call parse_input_integer('nphot@', nphot)
  call parse_input_integer('iseed@', iseed)
  call parse_input_integer('imethod@', imethod)
  call parse_input_integer('ifast@', ifast)
  call parse_input_double('enthres@', enthres)
  call parse_input_integer('cntdump@', cntdump)
  call parse_input_integer('irestart@', irestart)
  call parse_input_integer('itempdecoup@', itempdecoup)
  call parse_input_integer('iquantum@', iquantum)
  call parse_input_integer('istarsurf@', istarsurf)
  call parse_input_integer('nphotdiff@', nphotdiff)
  call parse_input_double('errtol@', errtol)
  call parse_input_integer('nvstr@', nvstr)
  call parse_input_double('vserrtol@', vserrtol)
  call parse_input_integer('ivstrt@', ivstrt)
  call parse_input_integer('ntemp@', ntemp)
  call parse_input_double('temp0@', temp0)
  call parse_input_double('temp1@', temp1)

  call check_all_lines_ok(ierror)
  if (ierror) then
     write(*,*) 'Error parsing radmc.inp'
     stop
  endif

  close(10)

  if (iseed > 0) iseed = -iseed  ! Ensure negative seed as per your original logic



  ! Read input files and parameters
  call read_all()

  ! Compute dust surface densities
  call calc_sigmadust()

  ! Write diagnostics
  if (dust_nr_species > 1) then
     do ispec = 1, dust_nr_species
        write(*,*) 'Dust mass in species ', ispec, ' = ', &
                   vstr_massdust(ispec) / 1.99d33, ' Msun'
     end do
  end if

  write(*,*) 'Dust mass total = ', vstr_massdusttot / 1.99d33, ' Msun'

  ! Vertical structure iteration block
  if (nvstr > 0) then

     ! Validate diffusion requirement
     if (nphotdiff == 0) then
        write(*,*) 'ERROR: Vertical structure iteration requires diffusion enabled.'
        write(*,*) 'Set nphotdiff to a value like 30.'
        stop 1
     end if

     ! Read vertical structure input
     call read_vstructinp()

     ! Backup dust density input file
     write(*,*) '...Copying dustdens.inp to dustdens_orig.inp'

     ! Modern safer system call (F2008+)
     call execute_command_line('cp -f dustdens.inp dustdens_orig.inp')

  end if

do iterstr = 0, nvstr

  ! Message
  if (iterstr >= 1) then
    write(*,*) '=====> Vertical structure iteration ', iterstr
  endif

  ! Emergency retry control for Monte Carlo
  iredo = 0

  do
    ! Run the Monte Carlo simulation
    call do_monte_carlo(nphot, spectrum, tdust, scatsrc, &
                        meanint, imethod, ifast, enthres, iseed,   &
                        cntdump, irestart, ntemp, temp0, temp1)

    ! Message for external luminosity contribution
    if (radbnd_extlumtot > 0.d0) then
      write(*,*) 'L_ext/L_* = ', radbnd_extlumtot / (radbnd_starlumtot + 1.d-99)
    endif

    ! Report skipped photon packages
    if (nrskips /= 0) then
      write(*,*) 'Number of photon packages skipped due to algorithmic problems = ', nrskips
      write(*,*) '  amounting to ', 100.d0 * nrskips / real(nphot) , '% of photon packages lost.'
    endif

    ! Do diffusion if requested
    if (nphotdiff > 0) then

      ! Write unmodified Monte Carlo result
      call write_results(1, spectrum, tdust, scatsrc, miquant, meanint, .false., .false.)

      ! Call the diffusion routine
      itmax = irsi_frsizey / 2
      call diffusion_main(irsi_frsizex, itmax, freq_nr, rsi_x_c(1,1), rsi_x_c(1,2), &
                     freq_nu, alpha_a, alpha_s, tdust, errtol, ierror_diff)

      ! If problem in diffusion, retry once
      if (ierror_diff /= 0) then
        if (iredo == 0) then
          iredo = 1
          write(*,*) 'Problem: The diffusion algorithm caught a problem.'
          ierror_diff = 0
          cycle  ! retry the Monte Carlo and diffusion
        else
          write(*,*) '##############################################'
          write(*,*) 'Repeated diffusion failure. Stopping.'
          stop
        endif
      endif

      ! Write diffusion info file
      call write_diffusion_info(irsi_frsizex, itmax, ndiff, irdiff, itdiff, ipde)

    else
      ! Clean up diffusion info files if no diffusion performed
      call delete_file('diffusion.info')
      call delete_file('dusttemp_mc.dat')
    endif

    exit  ! successful iteration
  end do


 ! Write final results after MC and (if any) diffusion
    writescat = (imethod == 2)
    writespec = .true.
    ifile = 0
    call write_results(ifile, spectrum, tdust, scatsrc, miquant, meanint, writescat, writespec)

  

    ! Vertical structure iteration
    if (iterstr < nvstr) then
      call vertstruct_integrate(tdust, ivstrt, error, iwarn)

      if (iwarn /= 0) then
        write(*,*) 'WARNING: Coarse grid near midplane.'
      endif

      call write_dustdens()

      write(*,*) 'Error in vertical structure = ', error
      if (error < errtol) then
        write(*,*) 'DONE...'
        exit
      endif
    endif

end do


 if (nvstr > 0) then
    write(*,*) 'Final vertical structure error = ', error
    if (nvstr >= 4) then
      write(*,*) 'Note: Likely due to Monte Carlo photon noise.'
    endif
    write(*,*) '==============================================================='
  else
    write(*,*) 'DONE...'
  endif


end subroutine radmc_main_program




  subroutine read_all()
    implicit none
    integer :: i, inu, ir, it, nr, nt, k, j, ispec, itmax, imirt
    integer :: iformat
    double precision :: dummy, frq, dum
    logical :: fex
    character(len=256) :: msg

    ! Read R grid
    call read_radius()

    ! Read Theta grid
    call read_theta()

    ! Read frequency grid
    call read_frequency()

    ! Read dust opacity data
    call read_dustdata()

    ! Read dust density
    call read_dust_density()

    ! Read stellar information
    call read_starinfo()

  end subroutine read_all


  subroutine read_radius()
    use common_grid
    implicit none
    integer :: ix
    real(8):: pi
    open(unit=10, file='radius.inp', status='old')
    read(10,*) irsi_frsizex
    if (irsi_frsizex > FRSIZE_X) then
      write(*,*) 'ERROR: radius.inp exceeds FRSIZE_X'
      stop 13
    endif

    do ix = 1, irsi_frsizex
      read(10,*) rsi_x_c(ix,1)
    enddo
    close(10)

    rsi_x_c(0,1)  = rsi_x_c(1,1)**2 / rsi_x_c(2,1)
    rsi_x_c(-1,1) = rsi_x_c(0,1)**2 / rsi_x_c(1,1)
    rsi_x_c(irsi_frsizex+1,1) = rsi_x_c(irsi_frsizex,1)**2 / rsi_x_c(irsi_frsizex-1,1)
    rsi_x_c(irsi_frsizex+2,1) = rsi_x_c(irsi_frsizex+1,1)**2 / rsi_x_c(irsi_frsizex,1)

    do ix = 0, irsi_frsizex+2
      rsi_x_i(ix,1) = sqrt(rsi_x_c(ix,1)*rsi_x_c(ix-1,1))
    enddo
  end subroutine read_radius


  subroutine read_theta()
    use common_grid
    implicit none
    integer :: iy, itmax, imirt
    real(8) :: pi
    open(unit=11, file='theta.inp', status='old')
    read(11,*) irsi_frsizey, imirt
    itmax = irsi_frsizey

    do iy=1,itmax
      read(11,*) rsi_x_c(iy,2)
      if(rsi_x_c(iy,2) < 0.d0) then
        write(*,*) 'ERROR: theta.inp has negative values.'
        stop 13
      endif
    enddo

    if(imirt /= 0) then
      irsi_frsizey = 2 * irsi_frsizey
    endif

    if(irsi_frsizey > FRSIZE_Y) then
      write(*,*) 'ERROR: theta.inp exceeds FRSIZE_Y'
      stop 13
    endif

    if(imirt /= 0) then
      do iy=1,itmax
        rsi_x_c(2*itmax+1-iy,2) = 3.1415926d0 - rsi_x_c(iy,2)
      enddo
    endif
    close(11)
    
    rsi_x_c(0,2)  = -rsi_x_c(1,2)
    rsi_x_c(-1,2) = -rsi_x_c(2,2)
    rsi_x_c(irsi_frsizey+1,2) = 2*3.1415926d0 - rsi_x_c(irsi_frsizey,2)
    rsi_x_c(irsi_frsizey+2,2) = 2*3.1415926d0 - rsi_x_c(irsi_frsizey-1,2)

    do iy = 0, irsi_frsizey+2
      rsi_x_i(iy,2) = 0.5d0 * (rsi_x_c(iy-1,2) + rsi_x_c(iy,2))
    enddo
  end subroutine read_theta


  subroutine read_frequency()
    use common_grid
    implicit none
    integer :: inu
    open(unit=12, file='frequency.inp', status='old')
    read(12,*) freq_nr
    if(freq_nr > FRSIZE_FREQ) then
      write(*,*) 'ERROR: frequency.inp exceeds FRSIZE_FREQ'
      stop 13
    endif
    do inu=1,freq_nr
      read(12,*) freq_nu(inu)
    enddo
    close(12)
  end subroutine read_frequency


  subroutine read_dust_density()
    use common_dust
    implicit none
    integer :: ir, it, ispec, nr, nt, k, j
    double precision :: dum

    open(unit=13, file='dustdens.inp', status='old')
    read(13,*) k, nr, nt, j
    if ((nr /= irsi_frsizex) .or. (nt /= irsi_frsizey/2)) stop 721
    if (k /= dust_nr_species) then
      write(*,*) 'ERROR: Dust species count mismatch.'
      stop 13
    endif

    do ispec=1,dust_nr_species
      do ir=1,nr
        do it=1,nt
          read(13,*) dum
          dust_rho(ispec,it,ir) = max(dum, 1d-90)
        enddo
      enddo
    enddo
    close(13)
  end subroutine read_dust_density


  subroutine read_starinfo()
    implicit none
    integer :: i, iformat
    double precision :: frq, dummy

    open(unit=14, file='starinfo.inp', status='old')
    read(14,*) iformat
    read(14,*) radbnd_rstar
    read(14,*) radbnd_mstar
    close(14)

    open(unit=15, file='starspectrum.inp', status='old')
    read(15,*) i
    if(i /= freq_nr) then
      write(*,*) 'ERROR: starspectrum.inp freq count mismatch.'
      stop 13
    endif

    if(i == 1) then
      read(15,*) frq, dummy
      radbnd_starspec(1) = 3.0308410d36 * dummy / radbnd_rstar**2
    else
      do i=1,freq_nr
        read(15,*) frq, dummy
        if(abs(frq-freq_nu(i))/(frq+freq_nu(i)) > 1.d-3) then
          write(*,*) 'ERROR: Frequency grid mismatch in stellar spectrum.'
          stop 13
        endif
        radbnd_starspec(i) = 3.0308410d36 * dummy / radbnd_rstar**2
      enddo
    endif
    close(15)
  end subroutine read_starinfo

!end module mod_input



subroutine write_results(ifile, spectrum, tdust, scatsrc, miquant, meanint, writescat, writespec)
  use configure
  use common_grid
  use common_dust
  use common_montecarlo
  use common_boundary
  implicit none
  ! Arguments
  integer, intent(in) :: ifile
  logical, intent(in) :: writescat, writespec
  double precision, intent(in) :: spectrum(FRSIZE_FREQ, FRSIZE_Y_SMALL)
  double precision, intent(in) :: tdust(DUST_SIZE_MAX, DUST_SPECIES_MAX, FRSIZE_Y_SMALL, FRSIZE_X)
  double precision, intent(in) :: scatsrc(FRSIZE_FREQ, FRSIZE_Y_SMALL, FRSIZE_X)
  double precision, intent(in) :: miquant(FRSIZE_FREQ, FRSIZE_Y_SMALL, FRSIZE_X)
  double precision, intent(in) :: meanint(FRSIZE_FREQ, FRSIZE_Y_SMALL, FRSIZE_X)

  ! Local variables
  integer :: inu, it, ir, itmax, imirt, isize, ispec, dummy, unit
  double precision :: dum

  ! Named constants
  double precision, parameter :: rad_to_deg = 57.2957795d0
  double precision, parameter :: four_pi_inv = 0.0795774715457d0


  ! Determine itmax based on MIRROR_THETA setting
!#ifdef MIRROR_THETA
  itmax = (irsi_frsizey+1)/2
  imirt = 1
!#else
!  itmax = irsi_frsizey
 ! imirt = 0
!#endif

  ! Find maximum number of dust sizes
  dummy = 0
  do ispec = 1, dust_nr_species
    if (dust_nr_size(ispec) > dummy) dummy = dust_nr_size(ispec)
  end do

  ! Write dusttemp.info
  open(newunit=unit, file='dusttemp.info', status='replace')
  write(unit, *) -2, 1, dummy, 1, 1
  close(unit)

  ! Write dust temperature
  if (ifile == 0) then
    write(*,*) '...Writing dust temperature (dusttemp_final.dat)'
    open(newunit=unit, file='dusttemp_final.dat', status='replace')
  else
    write(*,*) '...Writing dust temperature (dusttemp_mc.dat)'
    open(newunit=unit, file='dusttemp_mc.dat', status='replace')
  end if

  write(unit, *) dust_nr_species, irsi_frsizex, itmax, imirt
  write(unit, *)
  do ispec = 1, dust_nr_species
    write(unit, *) dust_nr_size(ispec)
    do isize = 1, dust_nr_size(ispec)
      do ir = 1, irsi_frsizex
        do it = 1, itmax
          write(unit, *) tdust(isize, ispec, it, ir)
        end do
      end do
    end do
  end do
  close(unit)

  ! Write spectrum
  if (writespec) then
    write(*,*) '...Writing Monte Carlo spectra (spectrum_all.dat)'
    open(newunit=unit, file='spectrum_all.dat', status='replace')
    write(unit, *) freq_nr, itmax
    write(unit, *)
    do it = 1, itmax
      write(unit, *) rad_to_deg * rsi_x_i(it, 2), rad_to_deg * rsi_x_i(it+1, 2)
      do inu = 1, freq_nr
        if (spectrum(inu, it) > 1.d-97) then
          write(unit, '(E21.13,1X,E21.13)') freq_nu(inu), spectrum(inu, it)
        else
          write(unit, '(E21.13,1X,E21.13)') freq_nu(inu), 0.d0
        end if
      end do
    end do
    close(unit)
  end if

  ! Write scattering source function
  if (writescat) then
    write(*,*) '...Writing isotropic scattering sources (scatsource.dat)'
    open(newunit=unit, file='scatsource.dat', status='replace')
    write(unit, *) freq_nr, irsi_frsizex, itmax, imirt
    write(unit, *)
    do inu = 1, freq_nr
      do it = 1, itmax
        do ir = 1, irsi_frsizex
          write(unit, *) scatsrc(inu, it, ir) * four_pi_inv
        end do
      end do
    end do
    close(unit)
  end if

end subroutine write_results




subroutine write_dustdens()
    use common_dust
    use common_grid
    implicit none
    integer :: it, ir, ispec
    integer :: nr, nt

    ! Set grid sizes
    nr = irsi_frsizex
    nt = irsi_frsizey / 2

    write(*,*) '...Writing dust density (dustdens.inp)'

    open(unit=10, file='dustdens.inp', status='replace', action='write', iostat=ir)
    if (ir /= 0) then
      write(*,*) 'Error opening file dustdens.inp'
      return
    end if

    ! Write header: number of species, nr, nt, 1
    write(10, *) dust_nr_species, nr, nt, 1

    ! Write dust densities: species, radius index, theta index
    do ispec = 1, dust_nr_species
      do ir = 1, nr
        do it = 1, nt
          write(10, '(E15.6)') dust_rho(ispec, it, ir)
        end do
      end do
    end do

    close(10)

  end subroutine write_dustdens




subroutine write_diffusion_info(irsize, itsize, ndiff, irdiff, itdiff, ipde)
  implicit none
  integer, intent(in) :: irsize, itsize, ndiff
  integer, intent(in) :: irdiff(ndiff), itdiff(ndiff), ipde(ndiff)
  integer :: ir, it, idiff
  integer, allocatable :: idifcell(:,:)

  allocate(idifcell(itsize, irsize))
  idifcell = 0

  do idiff = 1, ndiff
    ir = irdiff(idiff)
    it = itdiff(idiff)
    if (ipde(idiff) == 0) then
      idifcell(it, ir) = 1
    else
      idifcell(it, ir) = 2
    endif
  end do

  open(unit=10, file='diffusion.info', status='replace', action='write')
  write(10,*) 1
  write(10,*) irsize, itsize
  write(10,*)
  do it = 1, itsize
    write(10,*) (idifcell(it, ir), ir=1, irsize)
  end do
  close(10)

  deallocate(idifcell)
end subroutine write_diffusion_info



subroutine delete_file(filename)
  implicit none
  character(len=*), intent(in) :: filename
  logical :: exists
  inquire(file=filename, exist=exists)
 

  if (exists) then
    open(unit=20, file=filename, status='old')
    close(20, status='delete')
  endif
end subroutine delete_file




function stringcompare(string1, string2, len) result(is_equal)
  implicit none
  character(len=*), intent(in)  :: string1, string2
  integer, intent(in)      :: len
  logical                  :: is_equal
  integer                  :: i

  ! Error check on length
  if (len > 80) then
    write(*,*) 'ERROR: stringcompare length exceeds 80'
    stop
  endif

  if (len <= 0) then
    is_equal = .false.
    return
  endif

  ! Compare character-by-character
  is_equal = .true.
  do i = 1, len
    if (string1(i:i) /= string2(i:i)) then
      is_equal = .false.
    endif
  end do

end function stringcompare



subroutine parse_input_double(name, value)
  use common_montecarloarrays
  implicit none

  ! Arguments
  character(len=*), intent(in)  :: name
  double precision, intent(out) :: value

  ! Local variables
  integer :: iline, ichar, lenname, lenvalue, lenn
  character(len=80) :: strname, strvalue
  logical :: found, icomm

  ! Initialize
  found = .false.

  ! Check if name has a terminating '@'
  do ichar = 2, len_trim(name)
    if (name(ichar:ichar) == '@') exit
  end do

  if (ichar > len_trim(name)) then
    write(*,*) 'INTERNAL ERROR: Keywords must end with @!!!'
    stop
  end if

  lenn = ichar - 1

  ! Loop over all input lines to find the matching parameter
  do iline = 1, nrstring
    call parse_name_value(inpstring(iline), strname, lenname, strvalue, lenvalue, icomm)

    if (icomm) then
      lineok(iline) = 1

    elseif (stringcompare(strname, name, lenn) .and. (lenn == lenname)) then
      if (found) then
        write(*,*) 'Found one parameter more than once: ', strname(1:lenname)
        write(*,*) '*** ABORTING ***'
        stop
      end if

      read(strvalue(1:lenvalue), *) value
      found = .true.
      lineok(iline) = 1

    end if
  end do

end subroutine parse_input_double




subroutine parse_input_integer(name, value)
  use common_montecarloarrays
  implicit none

  ! Arguments
  character(len=*), intent(in)  :: name
  integer, intent(out)          :: value

  ! Local variables
  integer :: iline, ichar, lenname, lenvalue, lenn
  character(len=80) :: strname, strvalue
  logical :: found, icomm

  ! Initialize
  found = .false.

  ! Check if name has a terminating '@'
  do ichar = 2, len_trim(name)
    if (name(ichar:ichar) == '@') exit
  end do

  if (ichar > len_trim(name)) then
    write(*,*) 'INTERNAL ERROR: Keywords must end with @!!!'
    stop
  end if

  lenn = ichar - 1

  ! Loop over all input lines to find the matching parameter
  do iline = 1, nrstring
    call parse_name_value(inpstring(iline), strname, lenname, strvalue, lenvalue, icomm)

    if (icomm) then
      lineok(iline) = 1

    elseif (stringcompare(strname, name, lenn) .and. (lenn == lenname)) then
      if (found) then
        write(*,*) 'Found one parameter more than once: ', strname(1:lenname)
        write(*,*) '*** ABORTING ***'
        stop
      end if

      read(strvalue(1:lenvalue), *) value
      found = .true.
      lineok(iline) = 1

    end if
  end do

end subroutine parse_input_integer



subroutine check_all_lines_ok(ierror)
  use common_montecarloarrays
  implicit none

  logical, intent(out) :: ierror

  integer :: iline, ilen
  character(len=160) :: string

  ierror = .false.

  do iline = 1, nrstring
    if (lineok(iline) == 0) then
      write(*,*) 'ERROR in input file: variable unknown:'
      ilen = len_trim(inpstring(iline))
      string = inpstring(iline)
      write(*,*) string(1:ilen)
      ierror = .true.
    end if
  end do

end subroutine check_all_lines_ok



subroutine read_input_file(unit)
  use common_montecarloarrays
  implicit none

  integer, intent(in) :: unit
  integer :: iline
  character(len=160) :: string

  do iline = 1, INPUT_MAX_LEN
    read(unit, '(A)', end=10) string
    inpstring(iline) = string
    lineok(iline) = 0
  end do

  write(*,*) 'Input file contains too many lines (limit =', INPUT_MAX_LEN, ')'
  stop

10 continue
  nrstring = iline - 1

end subroutine read_input_file



subroutine parse_name_value(string, strname, lenname, strvalue, lenvalue, icomm)
  implicit none
  character(len=160), intent(in)  :: string
  character(len=80),  intent(out) :: strname, strvalue
  integer,             intent(out) :: lenname, lenvalue
  logical,             intent(out) :: icomm

  integer :: ilen, i, istart, istop

  ! Initialize
  icomm    = .false.
  strname  = ''
  strvalue = ''
  lenname  = 0
  lenvalue = 0

  ilen = len_trim(string)
  istart = 1

  ! Skip leading spaces
  do i = 1, ilen-1
    if (string(i:i) == ' ') then
      istart = i + 1
    else
      exit
    end if
  end do

  ! If line is comment or empty
  if (istart > ilen) then
    icomm = .true.
    return
  end if
  if (string(istart:istart) == ';') then
    icomm = .true.
    return
  end if

  ! Search for end of name: space or '='
  do i = istart+1, ilen-1
    if (string(i:i) == ' ' .or. string(i:i) == '=') then
      exit
    end if
  end do
  istop = i - 1
  lenname = istop - istart + 1

  if (lenname > 80) then
    write(*,*) 'Variable name too long'
    stop
  end if

  strname = string(istart:istop)

  ! Look for '=' after name
  do i = istop+1, ilen-1
    if (string(i:i) == '=') exit
  end do

  if (i >= ilen) then
    lenname  = 0
    strname  = ''
    lenvalue = 0
    strvalue = ''
    return
  end if

  ! Skip spaces after '='
  istart = i + 1
  if (istart > ilen) stop

  do i = istart, ilen
    if (string(i:i) /= ' ') then
      istart = i
      exit
    end if
  end do

  ! Find end of value: ';' or space
  do i = istart+1, ilen
    if (string(i:i) == ';' .or. string(i:i) == ' ') then
      exit
    end if
  end do
  istop = i - 1
  lenvalue = istop - istart + 1

  if (lenvalue > 80) then
    write(*,*) 'Value length too long'
    stop
  end if

  strvalue = string(istart:istop)

end subroutine parse_name_value


end module main


