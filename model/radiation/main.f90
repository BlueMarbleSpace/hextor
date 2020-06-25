program main
!==================================================================================
! Main program for testing the radiative transfer module
!==================================================================================

use radiation_mod

implicit none

real :: pg0     = 3.0
real :: fh2     = 0.0
real :: fch4    = 0.0
real :: tg0     = 300.0
real :: zy      = 60.0
real :: surfalb = 0.20
real :: olr, palb

call radiation_init

call getOLR( pg0, fh2, fch4, tg0, olr )
print *, "OLR = ", olr

call getPALB( pg0, fh2, fch4, tg0, zy, surfalb, palb )
print *, "PALB = ", palb

call radiation_end

end program main
