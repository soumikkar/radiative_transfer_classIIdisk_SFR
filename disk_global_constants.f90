module disk_global_constants
    implicit none

    real(kind=8), parameter :: GG   = 6.672d-8       ! Gravitational constant [cm^3/g/s^2]
    real(kind=8), parameter :: mp   = 1.6726d-24     ! Mass of proton [g]
    real(kind=8), parameter :: me   = 9.1095d-28     ! Mass of electron [g]
    real(kind=8), parameter :: kk   = 1.3807d-16     ! Boltzmann's constant [erg/K]
    real(kind=8), parameter :: hh   = 6.6262d-27     ! Planck's constant [erg.s]
    real(kind=8), parameter :: ee   = 4.8032d-10     ! Unit charge [esu]
    real(kind=8), parameter :: cc   = 2.9979d10      ! Light speed [cm/s]
    real(kind=8), parameter :: st   = 6.6524d-25     ! Thompson cross-section [cm^2]
    real(kind=8), parameter :: ss   = 5.6703d-5      ! Stefan-Boltzmann const [erg/cm^2/K^4/s]
    real(kind=8), parameter :: aa   = 7.5657d-15     ! 4 * SS / CC [erg/cm^3/K^4]
    
    real(kind=8), parameter :: muh2 = 2.3000d0       ! Mean molec weight H2+He+Metals
    
    real(kind=8), parameter :: ev   = 1.6022d-12     ! Electronvolt [erg]
    real(kind=8), parameter :: kev  = 1.6022d-9      ! Kilo electronvolt [erg]
    real(kind=8), parameter :: micr = 1.d-4          ! Micron [cm]
    real(kind=8), parameter :: km   = 1.d5           ! Kilometer [cm]
    real(kind=8), parameter :: angs = 1.d-8          ! Angstroem [cm]
    
    real(kind=8), parameter :: LS   = 3.8525d33      ! Solar luminosity [erg/s]
    real(kind=8), parameter :: RS   = 6.96d10        ! Solar radius [cm]
    real(kind=8), parameter :: MS   = 1.99d33        ! Solar mass [g]
    real(kind=8), parameter :: TS   = 5.78d3         ! Solar temperature [K]
    real(kind=8), parameter :: AU   = 1.496d13       ! Astronomical Unit [cm]
    real(kind=8), parameter :: pc   = 3.08572d18     ! Parsec [cm]
    
    real(kind=8), parameter :: pi   = 3.1415926535897932385d0 ! Mathematical constant pi

end module disk_global_constants

