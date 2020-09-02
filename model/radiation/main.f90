program main
!==================================================================================
! Main program for testing the radiative transfer module
!==================================================================================

use radiation_mod

implicit none

real :: pg0     = 1.0
real :: fh2     = 0.0
real :: fco2    = 4.0e-4
real :: tg0     = 288.0
real :: zy      = 60.0
real :: surfalb = 0.20
real :: olr, palb

call radiation_init

call getOLR( pg0, fh2, fco2, tg0, olr )
print *, "OLR = ", olr

call getPALB( pg0, fh2, fco2, tg0, zy, surfalb, palb )
print *, "PALB = ", palb

call radiation_end

end program main
