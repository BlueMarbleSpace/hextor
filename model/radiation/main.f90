program main
!==================================================================================
! Main program for testing the radiative transfer module
!==================================================================================

use radiation_mod

implicit none

real                :: pg0     = 1.0
real                :: zy      = 60.0
real                :: surfalb = 0.30
real                :: fco2, pco2, ptot, olr, palb, tg0
real, dimension(28) :: pco2levels

integer             :: temp, tstart=190, tend=370, tstep=1
integer             :: i

pco2levels = (/ 1.e-6,3.e-6,6.e-6,1.e-5, 3.e-5, 6.e-5,1.e-4, &
3.e-4, 6.e-4,1.e-3, 3.e-3, 6.e-3,1.e-2, 3.e-2, 6.e-2,1.e-1,  &
3.e-1, 6.e-1,1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)

call radiation_init

do temp = tstart, tend, tstep
  do i = 1, 28

    pco2 = pco2levels(i)

    ptot = pg0 + pco2
    fco2 = pco2 / ptot

    tg0 = real( temp )

    call getOLR( fco2, tg0, olr )
    !print *, olr / 1000.

    call getPALB( fco2, tg0, zy, surfalb, palb )
    print *, palb
  
  end do

end do

call radiation_end

end program main
