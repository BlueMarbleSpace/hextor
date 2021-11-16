module radiation_mod
!==================================================================================

  use hdf5
  use m_ffhash

!==================================================================================
implicit none

public :: radiation_init, radiation_end, getOLR, getPALB
!==================================================================================

  integer, parameter :: nolr = 4
  integer, parameter :: nalb = 6
  integer, parameter :: dolr = 1748
  integer, parameter :: dalb = 34960
  integer, parameter :: nco2 = 92
  integer, parameter :: ntmp = 19
  integer, parameter :: nzen = 4
  integer, parameter :: nsab = 5


  character*100 :: file_name = "./radiation/radiation_N2_CO2_Sun.h5"
  !character*100 :: file_name = "./radiation/radiation_N2_CO2_2600K.h5"
  !character*100 :: file_name = "./radiation_N2_CO2_2600K.h5"   !uncomment for standalone test with main.f90
  !character*100 :: file_name = "./radiation_N2_CO2_Sun.h5"   !uncomment for standalone test with main.f90
  character*100 :: sds_olr   = "/olr"
  character*100 :: sds_alb   = "/palb"
  integer(HID_T):: file_id, olr_id, alb_id
  integer       :: status, error 

  real, dimension(nolr,dolr)        :: database_olr
  real, dimension(nalb,dalb)        :: database_alb
  real, dimension(nolr)             :: weights
  integer(HSIZE_T), dimension(nolr) :: dims_olr 
  integer(HSIZE_T), dimension(nalb) :: dims_alb
  real, dimension(nco2)             :: fco2levels
  real, dimension(ntmp)             :: templevels
  real, dimension(nzen)             :: zenilevels
  real, dimension(nsab)             :: albelevels

  type(ffh_t) :: h_olr, h_alb

contains

!==================================================================================

subroutine radiation_init

  character*10 :: string1, string2, string3, string4
  character*40 :: stringkey
  integer      :: i

  call h5open_f( status )
  call h5fopen_f( file_name, H5F_ACC_RDONLY_F, file_id, status )

  call h5dopen_f( file_id, sds_olr, olr_id, error )
  call h5dread_f( olr_id, H5T_NATIVE_DOUBLE, database_olr, dims_olr, error )

  call h5dopen_f( file_id, sds_alb, alb_id, error )
  call h5dread_f( alb_id, H5T_NATIVE_DOUBLE, database_alb, dims_alb, error )

  do i = 1, dolr
    write( string1, '(F10.8)' ) database_olr(2,i)
    write( string2, '(F5.1)' ) database_olr(3,i)
    stringkey = string1 // string2
    call h_olr%ustore_value( stringkey, database_olr(4,i) )
  end do

  do i = 1, dalb
    write( string1, '(F10.8)' ) database_alb(2,i)
    write( string2, '(F5.1)' ) database_alb(3,i)
    write( string3, '(F4.1)' ) database_alb(4,i)
    write( string4, '(F3.1)' ) database_alb(5,i)
    stringkey = string1 // string2 // string3 // string4
    call h_alb%ustore_value( stringkey, database_alb(6,i) )
  end do

  fco2levels = database_olr(2,:nco2)
  templevels = database_olr(3,::nco2)
  zenilevels = (/ 0., 30., 60., 90. /)
  albelevels = (/ 0.2, 0.4, 0.6, 0.8, 1.0 /)

  return

end subroutine radiation_init

!==================================================================================

subroutine getOLR( fco2, tg0, olr )

  real, intent(in)  :: fco2
  real, intent(in)  :: tg0
  real, intent(out) :: olr

  real, dimension(nco2) :: fco2diff
  real, dimension(ntmp) :: tempdiff
  real, dimension(2)    :: co2pts, tmppts, olrpts
  real, dimension(2)    :: frac
  integer               :: co2ind, tmpind, nextind

  character*10 :: string1co2, string2co2, string1tmp, string2tmp
  character*20 :: stringkey1, stringkey2

  fco2diff = abs( fco2levels - fco2 ) / fco2
  tempdiff = abs( templevels - tg0 ) / tg0
  co2ind   = minloc( fco2diff, 1 )
  tmpind   = minloc( tempdiff, 1 )

  if ( ( fco2diff( co2ind ) .eq. 0.0 ) .or. ( fco2 .lt. fco2levels(1) ) .or. ( fco2 .gt. fco2levels(nco2) ) ) then
    nextind = co2ind
  else 
    if ( co2ind .eq. 1 ) then 
      nextind = co2ind + 1
    else if ( co2ind .eq. nco2 ) then 
      nextind = co2ind - 1
    else
      if ( fco2diff( co2ind - 1 ) .le. fco2diff( co2ind + 1 ) ) then
        nextind = co2ind - 1
      else 
        nextind = co2ind + 1
      end if
    end if
  end if
  co2pts(1) = fco2levels( co2ind )
  co2pts(2) = fco2levels( nextind )
  write( string1co2, '(F10.8)' ) co2pts(1)
  write( string2co2, '(F10.8)' ) co2pts(2)

  if ( ( tempdiff( tmpind ) .eq. 0.0 ) .or. ( tg0 .lt. templevels(1) ) .or. ( tg0 .gt. templevels(ntmp) ) ) then
    nextind = tmpind
  else 
    if ( tmpind .eq. 1 ) then 
      nextind = tmpind + 1
    else if ( tmpind .eq. ntmp ) then 
      nextind = tmpind - 1
    else
      if ( tempdiff( tmpind - 1 ) .le. tempdiff( tmpind + 1 ) ) then
        nextind = tmpind - 1
      else
        nextind = tmpind + 1
      end if
    end if
  end if
  tmppts(1) = templevels( tmpind )
  tmppts(2) = templevels( nextind )
  write( string1tmp, '(F5.1)' ) tmppts(1)
  write( string2tmp, '(F5.1)' ) tmppts(2)

  stringkey1 = string1co2 // string1tmp
  stringkey2 = string2co2 // string2tmp

  if ( stringkey1 .eq. stringkey2 ) then
    call h_olr%uget_value( stringkey1, olr )   
  else 
    call h_olr%uget_value( stringkey1, olrpts(1) )   
    call h_olr%uget_value( stringkey2, olrpts(2) )   
    if ( ( olrpts(1) .eq. -1.0 ) .or. ( olrpts(2) .eq. -1.0 ) ) then
      olr = 1.0
    else
      if ( co2pts(1) .eq. co2pts(2) ) then
        frac(1) = 0.0
      else
        frac(1) = abs( log10( fco2 ) - log10( co2pts(1) ) ) / abs( log10( co2pts(2) ) - log10( co2pts(1) ) )
      end if
      if ( tmppts(1) .eq. tmppts(2) ) then
        frac(2) = 0.0
      else
        frac(2) = abs( tg0 - tmppts(1) ) / abs( tmppts(2) - tmppts(1) )
      end if

      olr = olrpts(1) + maxval( frac ) * ( olrpts(2) - olrpts(1) )
      !olr = olrpts(1) + ( frac(2) ) * ( olrpts(2) - olrpts(1) )
      !olr = ( olrpts(2) + olrpts(1) ) / 2

    end if
  end if

  olr = -olr

  return

end subroutine getOLR

!==================================================================================

subroutine getPALB( fco2, tg0, zy, surfalb, palb )

  real, intent(in)  :: fco2
  real, intent(in)  :: tg0
  real, intent(in)  :: zy
  real, intent(in)  :: surfalb
  real, intent(out) :: palb

  real, dimension(nco2)             :: fco2diff
  real, dimension(ntmp)             :: tempdiff
  real, dimension(nzen)             :: zenidiff
  real, dimension(nsab)             :: albediff
  real, dimension(2)                :: co2pts, tmppts, zenpts, albpts, palbpts
  real, dimension(4)                :: frac
  integer                           :: co2ind, tmpind, zenind, albind, nextind

  character*10 :: string1co2, string2co2, string1tmp, string2tmp, string1zen, string2zen, string1alb, string2alb
  character*40 :: stringkey1, stringkey2

  fco2diff = abs( fco2levels - fco2 ) / fco2
  tempdiff = abs( templevels - tg0 ) / tg0
  zenidiff = abs( zenilevels - zy ) / zy
  albediff = abs( albelevels - surfalb ) / surfalb
  co2ind   = minloc( fco2diff, 1 )
  tmpind   = minloc( tempdiff, 1 )
  zenind   = minloc( zenidiff, 1 )
  albind   = minloc( albediff, 1 )

  if ( ( fco2diff( co2ind ) .eq. 0.0 ) .or. ( fco2 .lt. fco2levels(1) ) .or. ( fco2 .gt. fco2levels(nco2) ) ) then
    nextind = co2ind
  else
    if ( co2ind .eq. 1 ) then
      nextind = co2ind + 1
    else if ( co2ind .eq. nco2 ) then
      nextind = co2ind - 1
    else
      if ( fco2diff( co2ind - 1 ) .le. fco2diff( co2ind + 1 ) ) then
        nextind = co2ind - 1
      else
        nextind = co2ind + 1
      end if
    end if
  end if
  co2pts(1) = fco2levels( co2ind )
  co2pts(2) = fco2levels( nextind )
  write( string1co2, '(F10.8)' ) fco2levels( co2ind )
  write( string2co2, '(F10.8)' ) fco2levels( nextind )

  if ( ( tempdiff( tmpind ) .eq. 0.0 ) .or. ( tg0 .lt. templevels(1) ) .or. ( tg0 .gt. templevels(ntmp) ) ) then
    nextind = tmpind
  else
    if ( tmpind .eq. 1 ) then
      nextind = tmpind + 1
    else if ( tmpind .eq. ntmp ) then
      nextind = tmpind - 1
    else
      if ( tempdiff( tmpind - 1 ) .le. tempdiff( tmpind + 1 ) ) then
        nextind = tmpind - 1
      else
        nextind = tmpind + 1
      end if
    end if
  end if
  tmppts(1) = templevels( tmpind )
  tmppts(2) = templevels( nextind )
  write( string1tmp, '(F5.1)' ) templevels( tmpind )
  write( string2tmp, '(F5.1)' ) templevels( nextind )

  if ( ( zenidiff( zenind ) .eq. 0.0 ) .or. ( zy .lt. zenilevels(1) ) .or. ( zy .gt. zenilevels(nzen) ) ) then
    nextind = zenind
  else
    if ( zenind .eq. 1 ) then
      nextind = zenind + 1
    else if ( zenind .eq. nzen ) then
      nextind = zenind - 1
    else
      if ( zenidiff( zenind - 1 ) .le. zenidiff( zenind + 1 ) ) then
        nextind = zenind - 1
      else
        nextind = zenind + 1
      end if
    end if
  end if
  zenpts(1) = zenilevels( zenind )
  zenpts(2) = zenilevels( nextind )
  write( string1zen, '(F4.1)' ) zenilevels( zenind )
  write( string2zen, '(F4.1)' ) zenilevels( nextind )

  if ( ( albediff( albind ) .eq. 0.0 ) .or. ( surfalb .lt. albelevels(1) ) .or. ( surfalb .gt. albelevels(nsab) ) ) then
    nextind = albind
  else
    if ( albind .eq. 1 ) then
      nextind = albind + 1
    else if ( albind .eq. nsab ) then
      nextind = albind - 1
    else
      if ( albediff( albind - 1 ) .le. albediff( albind + 1 ) ) then
        nextind = albind - 1
      else
        nextind = albind + 1
      end if
    end if
  end if
  albpts(1) = albelevels( albind )
  albpts(2) = albelevels( nextind )
  write( string1alb, '(F3.1)' ) albelevels( albind )
  write( string2alb, '(F3.1)' ) albelevels( nextind )

  stringkey1 = string1co2 // string1tmp // string1zen // string1alb
  stringkey2 = string2co2 // string2tmp // string2zen // string2alb

  if ( stringkey1 .eq. stringkey2 ) then
    call h_alb%uget_value( stringkey1, palb )
  else
    call h_alb%uget_value( stringkey1, palbpts(1) )
    call h_alb%uget_value( stringkey2, palbpts(2) )
    if ( ( palbpts(1) .eq. -1.0 ) .or. ( palbpts(2) .eq. -1.0 ) ) then
      palb = -1.0
    else
      if ( co2pts(1) .eq. co2pts(2) ) then
        frac(1) = 0.0
      else
        frac(1) = abs( log10( fco2 ) - log10( co2pts(1) ) ) / abs( log10( co2pts(2) ) - log10( co2pts(1) ) )
      end if
      if ( tmppts(1) .eq. tmppts(2) ) then
        frac(2) = 0.0
      else
        frac(2) = abs( tg0 - tmppts(1) ) / abs( tmppts(2) - tmppts(1) )
      end if
      if ( zenpts(1) .eq. zenpts(2) ) then
        frac(3) = 0.0
      else
        frac(3) = abs( zy - zenpts(1) ) / abs( zenpts(2) - zenpts(1) )
      end if
      if ( albpts(1) .eq. albpts(2) ) then
        frac(4) = 0.0
      else
        frac(4) = abs( surfalb - albpts(1) ) / abs( albpts(2) - albpts(1) )
      end if

      palb = palbpts(1) + maxval( frac ) * ( palbpts(2) - palbpts(1) )
      !palb = ( palbpts(2) + palbpts(1) ) / 2
      !palb = palbpts(1) + ( frac(2) ) * ( palbpts(2) - palbpts(1) )
    end if
  end if

  return

end subroutine getPALB

!==================================================================================

subroutine radiation_end

  call h5dclose_f( olr_id, error )
  call h5dclose_f( alb_id, error )
  call h5fclose_f( file_id, error )
  call h5close_f( status )

  call h_olr%reset()
  call h_alb%reset()

end subroutine radiation_end

!==================================================================================

end module radiation_mod
