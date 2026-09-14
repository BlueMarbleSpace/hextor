program tabletest
!==================================================================================
! Standalone probe for a radiation lookup table.
!
!   ./tabletest <table.h5>  <  queries.txt
!
! Each input line is  pdry[bar]  fco2  fch4  tg0[K]  zenith[deg]  surfalb
! and produces one output line  pdry fco2 fch4 tg0 zy surfalb OLR PALB
! with OLR in W/m^2 (the table's mW/m^2 divided by 1000, as driver.f does).
! Against a table with no CH4 axis any fch4 is clamped away, so the same
! query file works for v1, v2 and v3 tables.
!
! Used by tools/check_table_reader.py to check the Fortran interpolation
! against an independent implementation reading the same HDF5 file.
!==================================================================================

use radiation_mod
implicit none

character(len=512) :: radfile
real :: pdry, fco2, fch4, tg0, zy, surfalb, olr, palb
integer :: ios

if ( command_argument_count() .lt. 1 ) then
  write(*,*) 'usage: tabletest <table.h5>  <  queries.txt'
  stop 1
end if
call get_command_argument( 1, radfile )

call radiation_init( trim(radfile) )

do
  read(*,*,iostat=ios) pdry, fco2, fch4, tg0, zy, surfalb
  if ( ios .ne. 0 ) exit
  call getOLR(  pdry, fco2, fch4, tg0, olr )
  call getPALB( pdry, fco2, fch4, tg0, zy, surfalb, palb )
  write(*,'(8ES20.10)') pdry, fco2, fch4, tg0, zy, surfalb, olr/1000., palb
end do

call radiation_end

end program tabletest
