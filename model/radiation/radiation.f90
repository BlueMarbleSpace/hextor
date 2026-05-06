module radiation_mod
!==================================================================================

  use hdf5

!==================================================================================
implicit none

public  :: radiation_init, radiation_end, getOLR, getPALB
private :: bracket
!==================================================================================

  integer, parameter :: nolr = 4
  integer, parameter :: nalb = 6
  integer, parameter :: dolr = 1748
  integer, parameter :: dalb = 34960
  integer, parameter :: nco2 = 92
  integer, parameter :: ntmp = 19
  integer, parameter :: nzen = 4
  integer, parameter :: nsab = 5


  character*100 :: sds_olr = "/olr"
  character*100 :: sds_alb = "/palb"
  integer(HID_T):: file_id, olr_id, alb_id
  integer       :: status, error

  real, dimension(nco2)             :: fco2levels
  real, dimension(ntmp)             :: templevels
  real, dimension(nzen)             :: zenilevels
  real, dimension(nsab)             :: albelevels

  real, dimension(nco2, ntmp)             :: olr_table
  real, dimension(nco2, ntmp, nzen, nsab) :: palb_table

  ! Effective OLR power-law exponents at upper and lower table boundaries,
  ! fitted from the last/first two temperature levels so extrapolation
  ! matches the actual table gradient rather than assuming T^4.
  real, dimension(nco2) :: n_eff_upper
  real, dimension(nco2) :: n_eff_lower

contains

!==================================================================================

subroutine radiation_init( radfile )

  character*(*), intent(in) :: radfile
  real, allocatable :: database_olr(:,:)
  real, allocatable :: database_alb(:,:)
  integer(HSIZE_T)  :: dims_olr(2)
  integer(HSIZE_T)  :: dims_alb(2)
  integer :: i, j, ico2, itmp, izen, ialb

  allocate( database_olr(nolr, dolr) )
  allocate( database_alb(nalb, dalb) )
  dims_olr = (/ nolr, dolr /)
  dims_alb = (/ nalb, dalb /)

  call h5open_f( status )
  call h5fopen_f( radfile, H5F_ACC_RDONLY_F, file_id, status )
  if ( status .ne. 0 ) then
    write(*,*) "radiation_init: cannot open HDF5 file: ", trim(radfile)
    stop
  end if

  call h5dopen_f( file_id, sds_olr, olr_id, error )
  call h5dread_f( olr_id, H5T_NATIVE_DOUBLE, database_olr, dims_olr, error )

  call h5dopen_f( file_id, sds_alb, alb_id, error )
  call h5dread_f( alb_id, H5T_NATIVE_DOUBLE, database_alb, dims_alb, error )

  ! Extract grid axis levels
  fco2levels = database_olr(2, :nco2)
  templevels = database_olr(3, ::nco2)
  zenilevels = (/ 0., 30., 60., 90. /)
  albelevels = (/ 0.2, 0.4, 0.6, 0.8, 1.0 /)

  ! OLR table: CO2 is inner loop, T is outer in the HDF5 data
  do j = 1, ntmp
    do i = 1, nco2
      olr_table(i, j) = database_olr(4, (j-1)*nco2 + i)
    end do
  end do

  ! PALB table: ordering in HDF5 is unknown, so map each row via its axis values
  do i = 1, dalb
    ico2 = minloc( abs( fco2levels - database_alb(2,i) ), 1 )
    itmp = minloc( abs( templevels - database_alb(3,i) ), 1 )
    izen = minloc( abs( zenilevels - database_alb(4,i) ), 1 )
    ialb = minloc( abs( albelevels - database_alb(5,i) ), 1 )
    palb_table(ico2, itmp, izen, ialb) = database_alb(6,i)
  end do

  ! Fit effective extrapolation exponents from boundary gradients in the table.
  ! olr_table stores negative values, so ratios are positive and log is valid.
  do i = 1, nco2
    n_eff_upper(i) = log( olr_table(i, ntmp)   / olr_table(i, ntmp-1) ) &
                   / log( templevels(ntmp)       / templevels(ntmp-1)   )
    n_eff_lower(i) = log( olr_table(i, 2)       / olr_table(i, 1)      ) &
                   / log( templevels(2)           / templevels(1)        )
  end do

  deallocate( database_olr )
  deallocate( database_alb )

  return

end subroutine radiation_init

!==================================================================================

! Find bracketing indices in a sorted array such that levels(lo) <= val < levels(hi).
! Clamps: returns lo=hi=1 if val <= levels(1), lo=hi=n if val >= levels(n).
subroutine bracket( levels, n, val, lo, hi )

  integer, intent(in)  :: n
  real,    intent(in)  :: levels(n), val
  integer, intent(out) :: lo, hi
  integer :: i

  if ( val .le. levels(1) ) then
    lo = 1
    hi = 1
  else if ( val .ge. levels(n) ) then
    lo = n
    hi = n
  else
    lo = 1
    do i = 1, n - 1
      if ( levels(i) .le. val .and. val .lt. levels(i+1) ) then
        lo = i
        exit
      end if
    end do
    hi = lo + 1
  end if

end subroutine bracket

!==================================================================================

subroutine getOLR( fco2, tg0, olr )

  real, intent(in)  :: fco2
  real, intent(in)  :: tg0
  real, intent(out) :: olr

  integer :: co2lo, co2hi, tmplo, tmphi
  real    :: co2frac, tmpfrac, olr0, olr1
  real    :: tg0_eval, n_eff

  ! Clamp T to table range for interpolation; extrapolate outside using the
  ! effective exponent fitted from the boundary gradient in radiation_init.
  tg0_eval = max( templevels(1), min( templevels(ntmp), tg0 ) )

  call bracket( fco2levels, nco2, fco2,     co2lo, co2hi )
  call bracket( templevels, ntmp, tg0_eval, tmplo, tmphi )

  if ( co2lo .eq. co2hi ) then
    co2frac = 0.0
  else
    co2frac = ( log10(fco2)              - log10(fco2levels(co2lo)) ) &
            / ( log10(fco2levels(co2hi)) - log10(fco2levels(co2lo)) )
  end if

  if ( tmplo .eq. tmphi ) then
    tmpfrac = 0.0
  else
    tmpfrac = ( tg0_eval - templevels(tmplo) ) / ( templevels(tmphi) - templevels(tmplo) )
  end if

  olr0 = olr_table(co2lo, tmplo) + co2frac * ( olr_table(co2hi, tmplo) - olr_table(co2lo, tmplo) )
  olr1 = olr_table(co2lo, tmphi) + co2frac * ( olr_table(co2hi, tmphi) - olr_table(co2lo, tmphi) )
  olr  = olr0 + tmpfrac * ( olr1 - olr0 )

  olr = -olr

  if ( tg0 .lt. templevels(1) .or. tg0 .gt. templevels(ntmp) ) then
    if ( tg0 .gt. templevels(ntmp) ) then
      n_eff = n_eff_upper(co2lo) + co2frac * ( n_eff_upper(co2hi) - n_eff_upper(co2lo) )
    else
      n_eff = n_eff_lower(co2lo) + co2frac * ( n_eff_lower(co2hi) - n_eff_lower(co2lo) )
    end if
    olr = olr * ( tg0 / tg0_eval )**n_eff
  end if

  return

end subroutine getOLR

!==================================================================================

subroutine getPALB( fco2, tg0, zy, surfalb, palb )

  real, intent(in)  :: fco2
  real, intent(in)  :: tg0
  real, intent(in)  :: zy
  real, intent(in)  :: surfalb
  real, intent(out) :: palb

  integer :: co2lo, co2hi, tmplo, tmphi, zenlo, zenhi, alblo, albhi
  real    :: co2frac, tmpfrac, zenfrac, albfrac

  ! pv(ico2, itmp, izen, ialb), each index 0=lo 1=hi
  real :: pv(0:1, 0:1, 0:1, 0:1)
  real :: p3d(0:1, 0:1, 0:1)
  real :: p2d(0:1, 0:1)
  real :: p1d(0:1)

  integer :: ic, it, iz, ia
  integer :: idx_co2(0:1), idx_tmp(0:1), idx_zen(0:1), idx_alb(0:1)

  call bracket( fco2levels, nco2, fco2,    co2lo, co2hi )
  call bracket( templevels, ntmp, tg0,     tmplo, tmphi )
  call bracket( zenilevels, nzen, zy,      zenlo, zenhi )
  call bracket( albelevels, nsab, surfalb, alblo, albhi )

  if ( co2lo .eq. co2hi ) then
    co2frac = 0.0
  else
    co2frac = ( log10(fco2)              - log10(fco2levels(co2lo)) ) &
            / ( log10(fco2levels(co2hi)) - log10(fco2levels(co2lo)) )
  end if

  if ( tmplo .eq. tmphi ) then
    tmpfrac = 0.0
  else
    tmpfrac = ( tg0 - templevels(tmplo) ) / ( templevels(tmphi) - templevels(tmplo) )
  end if

  if ( zenlo .eq. zenhi ) then
    zenfrac = 0.0
  else
    zenfrac = ( zy - zenilevels(zenlo) ) / ( zenilevels(zenhi) - zenilevels(zenlo) )
  end if

  if ( alblo .eq. albhi ) then
    albfrac = 0.0
  else
    albfrac = ( surfalb - albelevels(alblo) ) / ( albelevels(albhi) - albelevels(alblo) )
  end if

  idx_co2 = (/ co2lo, co2hi /)
  idx_tmp = (/ tmplo, tmphi /)
  idx_zen = (/ zenlo, zenhi /)
  idx_alb = (/ alblo, albhi /)

  ! Retrieve all 16 corners directly from the lookup table
  do ia = 0, 1
    do iz = 0, 1
      do it = 0, 1
        do ic = 0, 1
          pv(ic,it,iz,ia) = palb_table( idx_co2(ic), idx_tmp(it), idx_zen(iz), idx_alb(ia) )
        end do
      end do
    end do
  end do

  ! Quadrilinear interpolation: reduce one axis at a time
  do ia = 0, 1
    do iz = 0, 1
      do it = 0, 1
        p3d(it,iz,ia) = pv(0,it,iz,ia) + co2frac * ( pv(1,it,iz,ia) - pv(0,it,iz,ia) )
      end do
    end do
  end do

  do ia = 0, 1
    do iz = 0, 1
      p2d(iz,ia) = p3d(0,iz,ia) + tmpfrac * ( p3d(1,iz,ia) - p3d(0,iz,ia) )
    end do
  end do

  do ia = 0, 1
    p1d(ia) = p2d(0,ia) + zenfrac * ( p2d(1,ia) - p2d(0,ia) )
  end do

  palb = p1d(0) + albfrac * ( p1d(1) - p1d(0) )

  return

end subroutine getPALB

!==================================================================================

subroutine radiation_end

  call h5dclose_f( olr_id, error )
  call h5dclose_f( alb_id, error )
  call h5fclose_f( file_id, error )
  call h5close_f( status )

end subroutine radiation_end

!==================================================================================

end module radiation_mod
