module radiation_mod

!==================================================================================

  use hdf5

!==================================================================================
implicit none

public :: radiation_init, radiation_end, getOLR, getPALB
!==================================================================================

  integer, parameter :: nolr = 4
  integer, parameter :: nalb = 6
  integer, parameter :: dolr = 1748
  integer, parameter :: dalb = 34960

  character*100 :: file_name = "./radiation/radiation_N2_CO2.h5"
  !character*100 :: file_name = "./radiation/radiation_N2_CO2_2600K.h5"
  !character*100 :: file_name = "./radiation_N2_CO2.h5"   !uncomment for standalone test with main.f90
  character*100 :: sds_olr   = "/olr"
  character*100 :: sds_alb   = "/palb"
  integer(HID_T):: file_id, olr_id, alb_id
  integer       :: status, error 

  real, dimension(nolr,dolr)        :: database_olr
  real, dimension(nalb,dalb)        :: database_alb
  real, dimension(nolr)             :: weights
  integer(HSIZE_T), dimension(nolr) :: dims_olr 
  integer(HSIZE_T), dimension(nalb) :: dims_alb

contains

!==================================================================================

subroutine radiation_init

  call h5open_f( status )
  call h5fopen_f( file_name, H5F_ACC_RDONLY_F, file_id, status )

  call h5dopen_f( file_id, sds_olr, olr_id, error )
  call h5dread_f( olr_id, H5T_NATIVE_DOUBLE, database_olr, dims_olr, error )

  call h5dopen_f( file_id, sds_alb, alb_id, error )
  call h5dread_f( alb_id, H5T_NATIVE_DOUBLE, database_alb, dims_alb, error )

  return

end subroutine radiation_init

!==================================================================================

subroutine getOLR( pg0, fh2, fco2, tg0, olr )

  real, intent(in)  :: pg0
  real, intent(in)  :: fh2
  real, intent(in)  :: fco2
  real, intent(in)  :: tg0
  real, intent(out) :: olr

  real, dimension(nolr-1,dolr) :: difference
  real, dimension(nolr-1)      :: diff
  real, dimension(dolr)        :: distance
  real    :: sigma, weight1, weight2, rbf1, rbf2, rbf, dist
  integer :: neighbor1, neighbor2
  integer :: i, j

  difference(1,:) = ( database_olr(1,:) - pg0  )**2
  difference(2,:) = ( database_olr(2,:) - fco2  )**2
  !difference(3,:) = ( database_olr(3,:) - fh2 )**2
  !difference(4,:) = ( database_olr(4,:) - tg0  )**2
  difference(3,:) = ( database_olr(3,:) - tg0  )**2
  distance        = sum( difference, dim=1 )
  distance        = sqrt( distance )

  neighbor1 = 1
  neighbor2 = 2
  do i = 1, dolr
    if ( distance(i) .lt. distance( neighbor1 ) ) then
      neighbor1 = i
    end if    
    if ( ( distance(i) .lt. distance( neighbor2 ) ) .and. ( distance(i) .ne. distance( neighbor1 ) ) ) then
      neighbor2 = i
    end if
  end do
  !print *, neighbor1, neighbor2

  !if ( database_olr(5,neighbor1) .le. -1 ) then
  if ( database_olr(4,neighbor1) .eq. -1 ) then
    olr = -1
    return
  end if


  diff(1) = ( database_olr(1,neighbor1) - database_olr(1,neighbor2) )**2
  diff(2) = ( database_olr(2,neighbor1) - database_olr(2,neighbor2) )**2
  diff(3) = ( database_olr(3,neighbor1) - database_olr(3,neighbor2) )**2
  !diff(4) = ( database_olr(4,neighbor1) - database_olr(4,neighbor2) )**2
  dist      = sqrt( sum( diff ) )
  sigma     = ( distance( neighbor1 ) + distance( neighbor2 ) ) / 2  
  sigma     = 100000.0
  rbf       = exp( -dist / ( sigma**2 ) )

  weight1 = ( database_olr(4,neighbor2)*rbf - database_olr(4,neighbor1) ) / ( rbf*rbf - 1 )
  weight2 = ( database_olr(4,neighbor1)*rbf - database_olr(4,neighbor2) ) / ( rbf*rbf - 1 )
  !weight1 = ( database_olr(5,neighbor2)*rbf - database_olr(5,neighbor1) ) / ( rbf*rbf - 1 )
  !weight2 = ( database_olr(5,neighbor1)*rbf - database_olr(5,neighbor2) ) / ( rbf*rbf - 1 )
  rbf1    = exp( -distance( neighbor1) / ( sigma**2 ) )
  rbf2    = exp( -distance( neighbor2) / ( sigma**2 ) )

  olr = weight1*rbf1 + weight2*rbf2
  olr = -olr

  return

end subroutine getOLR

!==================================================================================

subroutine getPALB( pg0, fh2, fco2, tg0, zy, surfalb, palb )

  real, intent(in)  :: pg0
  real, intent(in)  :: fh2
  real, intent(in)  :: fco2
  real, intent(in)  :: tg0
  real, intent(in)  :: zy
  real, intent(in)  :: surfalb
  real, intent(out) :: palb

  real, dimension(nalb-1,dalb) :: difference
  real, dimension(nalb-1)      :: diff
  real, dimension(dalb)        :: distance
  real    :: sigma, weight1, weight2, rbf1, rbf2, rbf, dist
  integer :: neighbor1, neighbor2
  integer :: i, j

  difference(1,:) = ( database_alb(1,:) - pg0  )**2
  difference(2,:) = ( database_alb(2,:) - fco2  )**2
  !difference(3,:) = ( database_alb(3,:) - fh2 )**2
  !difference(4,:) = ( database_alb(4,:) - tg0  )**2
  !difference(5,:) = ( database_alb(5,:) - zy  )**2
  !difference(6,:) = ( database_alb(6,:) - surfalb  )**2
  difference(3,:) = ( database_alb(3,:) - tg0  )**2
  difference(4,:) = ( database_alb(4,:) - zy  )**2
  difference(5,:) = ( database_alb(5,:) - surfalb  )**2
  distance        = sum( difference, dim=1 )
  distance        = sqrt( distance )

  neighbor1 = 1
  neighbor2 = 2
  do i = 1, dalb
    if ( distance(i) .lt. distance( neighbor1 ) ) then
      neighbor1 = i
    end if    
    if ( ( distance(i) .lt. distance( neighbor2 ) ) .and. ( distance(i) .ne. distance( neighbor1 ) ) ) then
      neighbor2 = i
    end if
  end do
  !print *, neighbor1, neighbor2

  !if ( database_alb(5,neighbor1) .le. -1 ) then
  if ( database_alb(4,neighbor1) .le. -1 ) then
    palb = -1
    return
  end if

  diff(1) = ( database_alb(1,neighbor1) - database_alb(1,neighbor2) )**2
  diff(2) = ( database_alb(2,neighbor1) - database_alb(2,neighbor2) )**2
  diff(3) = ( database_alb(3,neighbor1) - database_alb(3,neighbor2) )**2
  diff(4) = ( database_alb(4,neighbor1) - database_alb(4,neighbor2) )**2
  diff(5) = ( database_alb(5,neighbor1) - database_alb(5,neighbor2) )**2
  !diff(6) = ( database_alb(6,neighbor1) - database_alb(6,neighbor2) )**2
  dist      = sqrt( sum( diff ) )
  sigma     = ( distance( neighbor1 ) + distance( neighbor2 ) ) / 2  
  sigma     = 100000.0
  rbf       = exp( -dist / ( sigma**2 ) )

  weight1 = ( database_alb(6,neighbor2)*rbf - database_alb(6,neighbor1) ) / ( rbf*rbf - 1 )
  weight2 = ( database_alb(6,neighbor1)*rbf - database_alb(6,neighbor2) ) / ( rbf*rbf - 1 )
  !weight1 = ( database_alb(7,neighbor2)*rbf - database_alb(7,neighbor1) ) / ( rbf*rbf - 1 )
  !weight2 = ( database_alb(7,neighbor1)*rbf - database_alb(7,neighbor2) ) / ( rbf*rbf - 1 )
  rbf1    = exp( -distance( neighbor1) / ( sigma**2 ) )
  rbf2    = exp( -distance( neighbor2) / ( sigma**2 ) )

  palb = weight1*rbf1 + weight2*rbf2

  return

end subroutine getPALB

!==================================================================================

subroutine radiation_end

  ! Terminate access to the datasets
  call h5dclose_f( olr_id, error )
  call h5dclose_f( alb_id, error )

  ! Terminate access to the file
  call h5fclose_f( file_id, error )

  ! Close FORTRAN interface.
  call h5close_f( status )

end subroutine radiation_end

!==================================================================================

end module radiation_mod
