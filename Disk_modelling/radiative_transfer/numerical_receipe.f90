module numerical_receipe
implicit none 

contains 

subroutine ludcmp(a, n, np, indx, d)
  implicit none

  ! Arguments
  integer, intent(in)              :: n, np
  real(8), intent(inout)           :: a(np, np)
  integer, intent(out)             :: indx(n)
  real(8), intent(out)             :: d

  ! Local variables
  real(8), parameter :: TINY = 1.0d-20
  integer            :: i, j, k, imax
  real(8)            :: aamax, dum, sum
  real(8), allocatable :: vv(:)

  ! Check np is large enough
  if (np < n) error stop 'np must be >= n in ludcmp'

  ! Allocate scaling vector
  allocate(vv(n))

  ! Initialize d
  d = 1.0d0

  ! Set up scaling factors
  do i = 1, n
    aamax = maxval(abs(a(i,1:n)))
    if (aamax == 0.0d0) error stop 'Singular matrix in ludcmp'
    vv(i) = 1.0d0 / aamax
  end do

  ! Crout's method
  do j = 1, n

    ! Compute upper part
    do i = 1, j-1
      sum = a(i,j)
      do k = 1, i-1
        sum = sum - a(i,k) * a(k,j)
      end do
      a(i,j) = sum
    end do

    ! Find pivot element
    aamax = 0.0d0
    do i = j, n
      sum = a(i,j)
      do k = 1, j-1
        sum = sum - a(i,k) * a(k,j)
      end do
      a(i,j) = sum
      dum = vv(i) * abs(sum)
      if (dum >= aamax) then
        imax = i
        aamax = dum
      end if
    end do

    ! Interchange rows if necessary
    if (j /= imax) then
      do k = 1, n
        dum         = a(imax,k)
        a(imax,k)   = a(j,k)
        a(j,k)      = dum
      end do
      d       = -d
      vv(imax)= vv(j)
    end if

    indx(j) = imax

    ! Avoid division by zero
    if (a(j,j) == 0.0d0) a(j,j) = TINY

    ! Divide by pivot element
    if (j /= n) then
      dum = 1.0d0 / a(j,j)
      do i = j+1, n
        a(i,j) = a(i,j) * dum
      end do
    end if

  end do

  ! Deallocate scaling vector
  deallocate(vv)

end subroutine ludcmp



SUBROUTINE lubksb(a,n,np,indx,b)
  INTEGER n,np,indx(n)
  DOUBLEPRECISION a(np,np),b(n)
  INTEGER i,ii,j,ll
  DOUBLEPRECISION sum

  ii=0
  do i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     if (ii.ne.0) then
        do j=ii,i-1
           sum=sum-a(i,j)*b(j)
        end do
     else if (sum.ne.0.d0) then
        ii=i
     end if
     b(i)=sum
  end do

  do i=n,1,-1
     sum=b(i)
     do j=i+1,n
        sum=sum-a(i,j)*b(j)
     end do
     b(i)=sum/a(i,i)
  end do
  RETURN
END


 subroutine hunt(xx, n, x, jlo)
  implicit none

  ! Arguments
  integer, intent(in)          :: n
  real(8), intent(in)          :: x
  real(8), intent(in)          :: xx(n)
  integer, intent(inout)       :: jlo

  ! Local variables
  integer                      :: inc, jhi, jm
  logical                      :: ascnd

  ! Check monotonicity
  ascnd = (xx(n) > xx(1))

  ! Initialize search bounds if needed
  if (jlo <= 0 .or. jlo > n) then
    jlo = 0
    jhi = n + 1
  else
    inc = 1
    if ((x >= xx(jlo)) .eqv. ascnd) then
      do
        jhi = jlo + inc
        if (jhi > n) then
          jhi = n + 1
          exit
        else if ((x >= xx(jhi)) .eqv. ascnd) then
          jlo = jhi
          inc = 2 * inc
        else
          exit
        end if
      end do
    else
      jhi = jlo
      do
        jlo = jhi - inc
        if (jlo < 1) then
          jlo = 0
          exit
        else if ((x < xx(jlo)) .eqv. ascnd) then
          jhi = jlo
          inc = 2 * inc
        else
          exit
        end if
      end do
    end if
  end if

  ! Bisection search
  do while (jhi - jlo > 1)
    jm = (jhi + jlo) / 2
    if ((x > xx(jm)) .eqv. ascnd) then
      jlo = jm
    else
      jhi = jm
    end if
  end do

end subroutine hunt




  subroutine huntalt(xx, n, x, i)
  implicit none

  ! Arguments
  integer, intent(in)           :: n
  real(8), intent(in)           :: x
  real(8), intent(in)           :: xx(n)
  integer, intent(inout)        :: i

  ! Handle exact matches at boundaries
  if (x == xx(1)) then
    i = 1
    return
  end if

  if (x == xx(n)) then
    i = n
    return
  end if

  ! Call hunt to bracket x in xx
  call hunt(xx, n, x, i)

  ! Check for out-of-bounds results
  if (i >= n) then
    i = n + 1
    return
  else if (i <= 0) then
    i = 0
    return
  end if

  ! If exact match to xx(i+1), adjust index
  if (x == xx(i+1)) then
    i = i + 1
  end if

end subroutine huntalt



  function ran2(idum) result(randnum)
    implicit none
    integer, intent(inout) :: idum
    real(8) :: randnum

    ! Parameters
    integer, parameter :: IM1=2147483563, IM2=2147483399
    integer, parameter :: IA1=40014, IA2=40692
    integer, parameter :: IQ1=53668, IQ2=52774
    integer, parameter :: IR1=12211, IR2=3791
    integer, parameter :: NTAB=32, IMM1=IM1-1
    integer, parameter :: NDIV=1 + IMM1/NTAB
    real(8), parameter  :: AM=1.0d0/IM1
    real(8), parameter  :: EPS=1.2d-7
    real(8), parameter  :: RNMX=1.0d0 - EPS

    ! Local variables
    integer :: idum2, j, k
    integer, save :: iy = 0
    integer, save :: iv(NTAB) = 0
    integer, save :: idum2_save = 123456789

    if (idum <= 0) then
      idum = max(-idum, 1)
      idum2_save = idum
      do j = NTAB+8, 1, -1
        k = idum / IQ1
        idum = IA1 * (idum - k*IQ1) - k*IR1
        if (idum < 0) idum = idum + IM1
        if (j <= NTAB) iv(j) = idum
      end do
      iy = iv(1)
    end if

    k = idum / IQ1
    idum = IA1 * (idum - k*IQ1) - k*IR1
    if (idum < 0) idum = idum + IM1

    k = idum2_save / IQ2
    idum2_save = IA2 * (idum2_save - k*IQ2) - k*IR2
    if (idum2_save < 0) idum2_save = idum2_save + IM2

    j = 1 + iy / NDIV
    iy = iv(j) - idum2_save
    iv(j) = idum
    if (iy < 1) iy = iy + IMM1

    randnum = min(AM * dble(iy), RNMX)

  end function ran2




subroutine make_indexed_filename(base, index, ext, filename)
    implicit none
    ! Arguments
    character(len=*), intent(in)  :: base
    integer,          intent(in)  :: index
    character(len=*), intent(in)  :: ext
    character(len=*), intent(out) :: filename

    ! Local variables
    character(len=12) :: ch

    ! Check for valid index range
    if ((index < 0) .or. (index >= 1000)) then
        write(*,*) 'ERROR in make_indexed_filename()'
        stop 729
    end if

    ! Format index into character string according to value
    if (index < 10) then
        write(ch, 11) index
11      format(I1)
    elseif (index < 100) then
        write(ch, 12) index
12      format(I2)
    elseif (index < 1000) then
        write(ch, 13) index
13      format(I3)
    end if

    ! Construct the final filename
    filename = base(1:len_trim(base)) // ch(1:len_trim(ch)) // ext(1:len_trim(ext))

end subroutine make_indexed_filename





function number_invalid(a) result(flag)
    implicit none
    ! Arguments
    real(8), intent(in) :: a
    ! Result
    integer             :: flag
    ! Local variables
    real(8)             :: b, c
    logical             :: div, sub

    ! Perform operations to check numerical validity
    b = a * 2.d0
    b = b / 2.d0
    c = a - 1.d100

    div = (b == a)
    sub = (c < a)

    if (div .and. sub) then
        flag = 0    ! valid number
    elseif (div) then
        flag = 1    ! NaN
    else
        flag = 2    ! Inf
    end if

end function number_invalid

end module numerical_receipe
