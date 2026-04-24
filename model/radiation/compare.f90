program compare
!==================================================================================
! Sweep OLR and PALB across (fco2, T, zenith, surfalb) and print results.
! Run this against both radiation.f90 versions, then compare outputs.
! Output columns: fco2  tg0  zy  surfalb  OLR  PALB
!==================================================================================

use radiation_mod
implicit none

real :: pg0 = 1.0
real :: fco2, pco2, ptot, olr, palb, tg0, zy, surfalb

! CO2 partial pressures spanning many decades (same set as main.f90)
real, dimension(28) :: pco2levels = (/ &
  1.e-6, 3.e-6, 6.e-6, 1.e-5, 3.e-5, 6.e-5, 1.e-4,  &
  3.e-4, 6.e-4, 1.e-3, 3.e-3, 6.e-3, 1.e-2, 3.e-2,  &
  6.e-2, 1.e-1, 3.e-1, 6.e-1, 1.,    2.,    3.,      &
  4.,    5.,    6.,    7.,    8.,    9.,    10. /)

! Temperatures: on-grid values (round hundreds) and off-grid midpoints
real, dimension(8) :: tmptest = (/ 200., 220., 250., 265., 288., 300., 325., 350. /)

! Zenith angles: 0 (edge case), midpoints between grid nodes, and grid nodes
real, dimension(7) :: zentest = (/ 0., 15., 30., 45., 60., 75., 90. /)

! Surface albedos: grid nodes and midpoints
real, dimension(7) :: albtest = (/ 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0 /)

integer :: i, j, k, l

call radiation_init( "./radiation_N2_CO2_Sun.h5" )

write(*,'(A)') '# fco2  tg0  zy  surfalb  OLR  PALB'

do l = 1, 7
  surfalb = albtest(l)
  do k = 1, 7
    zy = zentest(k)
    do j = 1, 8
      tg0 = tmptest(j)
      do i = 1, 28
        pco2 = pco2levels(i)
        ptot = pg0 + pco2
        fco2 = pco2 / ptot

        call getOLR(  fco2, tg0, olr )
        call getPALB( fco2, tg0, zy, surfalb, palb )

        write(*,'(6ES16.8)') fco2, tg0, zy, surfalb, olr, palb
      end do
    end do
  end do
end do

call radiation_end

end program compare
