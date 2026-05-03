module disk_density
    implicit none

contains

    subroutine disk_model_1(r, theta, rdisk, sigmadust0, p1, p2, hr0, plh, &
                             rhodust, sigdust, hrstore, hrmin, hrpuff, rpuff)

        real(kind=8), intent(in), dimension(:) :: r, theta
        real(kind=8), intent(in) :: rdisk, sigmadust0, p1, p2, hr0, plh
        real(kind=8), intent(out), allocatable, dimension(:,:,:) :: rhodust
    
        real(kind=8), intent(out), allocatable, dimension(:) :: sigdust, hrstore
        real(kind=8), intent(in), optional :: hrmin, hrpuff, rpuff

        integer :: nr, nt, ir, it, iwarn_negpuff
        real(kind=8) :: thmax, hr_local, hrmin_val, hrpuff0, eps
        real(kind=8) :: pi_val, normal_const, pl_val

        nr = size(r)
        nt = size(theta)
     
        pi_val = acos(-1.0d0) 
        thmax = pi_val / 2.0d0
        normal_const = 1.0d0 / sqrt(2.0d0 * pi_val)

        allocate(rhodust(nr, nt, 1), sigdust(nr), hrstore(nr))

        if (present(hrmin)) then
            hrmin_val = hrmin
        else
            hrmin_val = 0.0d0
   
        end if 

        hrpuff0 = 0.0d0
        if (present(rpuff)) then
            hrpuff0 = hr0 * (rpuff / rdisk)**plh
            if (.not.&
 present(hrpuff)) then 
                write(*,*) "ERROR: If you specify rpuff, you must also specify hrpuff"
                stop
            end if
        end if
        
        iwarn_negpuff = 0

        do ir = 1, nr
     
            hr_local = hr0 * (r(ir) / rdisk)**plh 
            hr_local = sqrt(hr_local**2 + hrmin_val**2)

            if (present(rpuff) .and. r(ir) < rpuff) then
                eps = (log(r(ir)) - log(r(1))) / (log(rpuff) - log(r(1)))
                hr_local = (1.0d0 - eps) * hrpuff + eps * hrpuff0
   
                if (hr_local < (hr0 * (r(ir)/rdisk)**plh)) then 
                    iwarn_negpuff = 1
                end if
            end if
            hrstore(ir) = hr_local
        end do

      
        if (iwarn_negpuff /= 0) then 
            write(*,*) 'Warning: The parameterized puffed-up inner rim has thickness'
            write(*,*) '          smaller than the normal flaring part.'
            write(*,*) 'Are you sure' 
            write(*,*) '          that you want this?'
            end if 

        do ir = 1, nr
            if (r(ir) <= rdisk) then
                pl_val = p1
            else
                pl_val = p2
            end if
          
            hr_local = hrstore(ir) 
            sigdust(ir) = sigmadust0 * (r(ir) / rdisk)**pl_val
            
            do it = 1, nt
                rhodust(ir, it, 1) = normal_const * sigdust(ir) * &
              
                        exp(-0.5d0 * ((thmax - theta(it)) / hr_local)**2) / &
                                     (hr_local * r(ir))
            end do
        end do

    end subroutine disk_model_1

end module disk_density 
