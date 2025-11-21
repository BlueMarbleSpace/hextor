      program  energy_balance_climate_model
c--------------------------------------------------------------------c
c  This program calculates seasonal latitudinal temperature 
c  gradients for Earth-like planets with N2,O2,CO2, and H20 
c  atmospheres (e.g. Earth and Mars). Energy balance is treated as in 
c  Caldeira and Kasting (1992). Input files: 'oceans.dat', for
c  modeling Earth with present geography, and 'fresnel_reflct.dat'.
c
c  Implemented lookup table for radiative transfer 01--24-18 JDH
c
c  Updated carbonate-silicate cycle 03-10-15 JDH
c
c  Implemented a stochastic mode and time-dependent orbital properties 02-13-14 JDH
c
c  Added namelist input: 05-26-09 JDH
c
c  Changes made and commented: 06-30-04 JDH
c
c  Original version: 11-3-95 DMW
c
c
c
c
c          lat => latitude (radians)
c            x => SIN(latitude)
c     latangle => latitude (degrees) 
c       nbelts => number of latitudinal belts  
c            t => time (seconds)
c         temp => surface temperature (K) 
c         pco2 => atmospheric partial pressure of CO2
c            c => thermal heat capacity (cl,ci,cw = land,ice,water)
c       focean => oceanic surface fraction within latitude belt
c        ocean => fraction of planet covered by water
c         fice => frozen water surface fraction in latitude belt  
c            h => solar hour-angle
c            z => solar zenith-angle (degrees)
c           mu => COS(solar zenith-angle)
c            s => diurnal-average solar insolation (Wm^-2)
c           q0 => solar constant at 1 AU (Wm^-2)
c            q => solar irradiance at arbitrary orbital distance (Wm^-2)
c           ir => outgoing infrared flux (Wm^-2)
c            d => thermal diffusion coefficient
c     coastlat => coastal latitude for polar and equatorial geographies (deg)
c      surfalb => surface albedo
c         atoa => top-of-atmosphere (planetary) albedo
c          dec => solar declination (radians)
c     decangle => solar declination (degrees)
c          obl => obliquity (degrees)
c          ecc => orbit eccentricity
c            a => orbit semi-major axis (cm)
c         peri => longitude of perihelion wrt vernal equinox (degrees)
c      cloudir => reduction of outgoing infrared by clouds (Wm^-2)
c         msun => solar mass (g)
c            m => mean anomaly (radians)
c            e => eccentric anomaly (radians)
c            w => orbital angular velocity (radians per second)
c            r => orbital distance from sun (cm)
c
c----------------------------------------------------------------------c

      use radiation_mod
    
      real*4  lat,latangle,mu,ir,msun,landalb,icealb,irsum,irave,mp,
     &  iravesum,iceline,icelat,icesht
      real*8  ecc,m,e

      parameter (nbelts=18)
      parameter (pi=3.14159265,grav=6.6732e-8)
      parameter (mp=1.67e-24,cnvg=1.e-1)
      parameter (sbc=5.67e-8,emis=0.64)
      parameter (twopi=2*pi)
      !parameter (niter=1)
      !parameter (niter=300000)
      !parameter (niter=1200000)  ! for do_futuresol and yrstep = 1.e3 (future Earth to 1.2 Gyr from present)
      !parameter (niter=120000)  ! for do_futuresol and yrstep = 1.e4 (future Earth to 1.2 Gyr from present)
      !parameter (niter=24000)  ! for do_futuresol and yrstep = 5.e4 (future Earth to 1.2 Gyr from present)
      !parameter (niter=12000)  ! for do_futuresol and yrstep = 1.e5 (future Earth to 1.2 Gyr from present)
      !parameter (niter=1200)  ! for do_futuresol and yrstep = 1.e6 (future Earth to 1.2 Gyr from present)
      !parameter (niter=120)  ! for do_futuresol and yrstep = 1.e7 (future Earth to 1.2 Gyr from present)
      !parameter (niter=4570)  ! for do_futuresol and yrstep = 1.e6 (Mars)
      !parameter (niter=45700)  ! for do_futuresol and yrstep = 1.e5 (Mars)
      !parameter (niter=91400)  ! for do_futuresol and yrstep = 5.e4 (Mars)
      !parameter (niter=457000)  ! for do_futuresol and yrstep = 1.e4 (Mars)
      parameter (niter=50)
      !parameter (niterhalf=1001)
      !parameter (niterquarter=1501)

      parameter (ndays=1627)
   
      common/geogblk/cntrlat(nbelts),cntrlatangle(nbelts),coastlat,
     &   focean(nbelts),lat(0:nbelts+1),latangle(0:nbelts+1),ocean

      dimension  temp(0:nbelts+1),c(nbelts),dx(0:nbelts),
     &  tprime(0:nbelts),tprimeave(nbelts),t2prime(nbelts),
     &  x(0:nbelts+1),h20alb(0:90),ir(nbelts),
     &  mu(nbelts),s(nbelts),atoa(nbelts),surfalb(nbelts),
     &  acloud(nbelts),zntempmin(nbelts),zntempmax(nbelts),area(nbelts),
     &  zntempsum(nbelts),zntempave(0:nbelts+1),zndecmax(nbelts),
     &  zndecmin(nbelts),obstemp(nbelts),iceline(0:5),fice(nbelts),
     &  wthrate(nbelts),warea(nbelts),imco2(nbelts), diff(nbelts),
     &  stab(0:nbelts), znalbsum(nbelts), znalbave(nbelts),
     &  znolrsum(nbelts), znolrave(nbelts), 
     &  znsurfalbsum(nbelts), znsurfalbave(nbelts)

      character  header*80,file(0:3)*8
      logical seasons, last, linrad, linalb, cloudalb
      logical do_stochastic, soladj, constheatcap, diffadj, iterhalt
      logical do_cs_cycle, do_h2_cycle, oceanalbconst 
      logical do_longitudinal, do_manualseasons, do_gough, do_marshist
      logical fillet, do_dailyoutput, do_futuresol, do_bioprod
      real landsnowfrac, RAND, boxmuller, noisevar, heatcap, ocnalb
      real outgassing, weathering, betaexp, kact, krun, q0
      real pg0, ir2, fh2, co2sat, h2escape, ph2, ncolh2, h2outgas
      real icelineN, icelineS
      real icelineNMax, icelineNMin, icelineSMax, icelineSMin
      real cl, cw, ci, pco20, pco2soil0, pco2soil, gamma, gammaout
      real pmin, phalf, gamma0
      integer yrcnt, yrstep, radparam, co2flag
      integer ISEED, resfile, nt, daynum
      integer*4 now(3)
      real total, snowalb, tempinit, solarcon, fco2, icetemp
      dimension solcon(niter),prec(niter),ecce(niter),
     &  yrlabel(niter),obliq(niter)
      dimension solconD(ndays),precD(ndays),ecceD(ndays),
     &  yrlabelD(ndays),obliqD(ndays)

      data  temp/20*273./       !**use some temperature > 265K, otherwise iceball
      data  file/ 'spng.out', 'summ.out', 'fall.out', 'wint.out' /


      NAMELIST /ebm/ seasons, tend, dt, rot, a, ecc, peri, 
     &               obl, ocean, igeog, yrstep, resfile, d0,
     &               constheatcap, heatcap, diffadj,
     &               iterhalt, fco2, fh2, pg0, tempinit, msun,
     &               do_longitudinal, do_manualseasons,
     &               cl, cw, ci, do_dailyoutput

      NAMELIST /radiation/ relsolcon, radparam, groundalb, snowalb,
     &               landsnowfrac, cloudir, fcloud, cloudalb, soladj,
     &               linrad, linalb, solarcon, oceanalbconst, ocnalb,
     &               do_gough, do_futuresol

      NAMELIST /co2cycle/ do_cs_cycle, outgassing, weathering,
     &               betaexp, kact, krun, pco20, pco2soil0, 
     &               do_bioprod

      NAMELIST /h2cycle/ do_h2_cycle, h2outgas

      NAMELIST /stochastic/ do_stochastic, noisevar


c  NAMELIST PARAMETERS
      INTEGER namelistid
      OPEN( namelistid, FILE='./input.nml', DELIM='APOSTROPHE' )

c  INITIALIZE VARIABLES
      total = 0
      ann_tempave = 10.         !just some constant greater than zero
      seasons = .true.     !if .true. do seasons: otherwise do mean-annual
      last = .false.       !should always be .false.
      do_stochastic = .false. !if .true. then include stochastic perturbation
      soladj = .false.     !if .true. then include solar forcing adjustment
      do_gough = .false.    !if .true. then use Gough (1981) solar forcing
      do_futuresol = .false.    !if .true. then use solar forcing from Caldeira and Kasting (1992)
      resfile = 0          !start from previous (1), present Earth (2) hothouse (3), large ice (4)
      q0 = 1360.           !solar constant
      tend = 7.e11          !calculation length (sec)
      dt = 8.64e4            !time step (sec)
      rot0 = 7.27e-5       !Earth's present rotation rate (rad/sec)
      rot = rot0           
      a0 = 1.49597892E13   !Earth's orbit semi-major axis (cm)
      a = a0               !semi-major axis for this run
      ecc = 0.016708617   !orbit eccentricity
      peri = 76.25     !longitude of perihelion wrt vernal equinox (degrees)
      obl = 23.45               !obliquity (degrees)
      cloudir = -9.5   !correction to outgoing infrared for clouds (W/m^-2)
      pco20 = 3.e-4   !present atmospheric CO2 concentration (bars)
      pco2 = pco20              !initial co2 concentration for this run (bars)
      pco2soil0 = 3.e-3   !present soil CO2 concentration (bars)
      pco2soil = pco2soil0 !initial co2 soil concentration for this run (bars)
      ocean = 0.7      !planet fraction covered by water
      igeog = 1        !geography ('1'=present,'2'=polar,'3'=equatorial,
     &                 ! '4'=100% water,'5'=100% land, '6'=equal lat frac)
      groundalb = 0.20  !base land surface albedo
      !groundalb = 0.25
      relsolcon = 1.0  !relative solar constant (1.0 for present Earth) 
      alpha = -0.078            !alpha: cloud albedo = alpha + beta*zrad
      beta = 0.65      !beta
      cw = 2.1e8       !heat capacity over ocean
      cl = 5.25e6      !heat capacity over land
      ci = 1.05e7      !heat capacity over lice
      d0 = 0.58        !thermal diffusion coefficient
      v = 3.3e14  !volcanic outgassing (g/year) [from Holland (1978)]
      wco2 = 9.49e14  !carbonate-silicate weathering constant
      avemol0 = 28.89*mp  !mean-molecular weight of present atmosphere (g)
      hcp0 = 0.2401    !heat capacity of present atmosphere (cal/g K)
      hcpco2 = 0.2105  !heat capacity of co2 (cal/g K)
      hcpn2 = .2484     ! heat capacity of n2
      hcph2 = 3.420      !heat capacity of h2 
      landsnowfrac = 1.0 !snowfall fraction on land
      fcloud = 0.5     !fractional cloud cover
      noisevar = 0.15  ! noise variance (K^2 per year)
      co2forcing = 3.7  ! radiative forcing from doubling co2 (w/m^2)
      linrad = .false.   ! set .true. for A+BT absorption
      linalb = .false.   ! set .true. for constant TOA albedo
      constheatcap = .false. ! set .true. for constant heat capacity
      diffadj = .true.  ! set .false. to turn off diffusion parameter adjustment
      cloudalb = .true. ! set .false. to disable cloud albedo
      iterhalt = .false.  ! set .true. to enable halt based on iterations
      do_cs_cycle = .false. !set .true. to enable carbonate-silicate cycle
      do_h2_cycle = .false. !set .true. to enable H2 cycle
      do_bioprod = .false.  !set .true. to enable the biological productivity function from Caldeira and Kasting (1992)
      outgassing = 7.0  ! volcanic outgassing rate (bars/Gyr)
      weathering = 7.0  ! weathering rate (bars/Gyr)
      betaexp = 0.50       ! weathering exponent
      kact = 0.09       ! activation energy factor
      krun = 0.045      ! runoff efficiency factor
      pg0 = 1.0         ! surface pressure (bars)
      fh2 = 0.0
      h2outgas = 2.67e12
      radparam = 4      ! OLR/albedo parameterization: (0) Williams & Kasting, (1) CO2/N2, (2) CO2/H2, (3) lookup table, (4) New CO2/H2
      nfile = 0
      yrstep = 1
      yrcnt = 0
      yricnt = 1
      msun = 1.9891e33
      heatcap = cl
      snowalb = 0.663 ! changed this from 0.7 as per Caldiera & Kasting (JDH), added to namelist
      tempinit = 273.16
      do_longitudinal = .false.
      nt = 1          !counter for number of timesteps 
      daynum = 1          !counter for number of days per orbit
      icetemp = 263.15    !threshold for iceline
      gammaout = 1.0      !initial value of biological productivity function relative to present Earth

      do_dailyoutput = .false.
      fillet = .true.
      do_marshist = .false.

      CALL itime( now )
      CALL SRAND( now(3) )

      READ( namelistid, NML=ebm )
      READ( namelistid, NML=radiation )
      READ( namelistid, NML=co2cycle )
      READ( namelistid, NML=h2cycle )
      READ( namelistid, NML=stochastic )
      CLOSE( namelistid )

 
c  OPEN FILES 
      open (unit=2,file='data/oceans.dat',status='old')
      open (unit=3,file='data/fresnel_reflct.dat',status='old')
      open (unit=4,file='data/presentEarth.dat',status='old')
      open (unit=5,file='data/hothouse.dat',status='old')
      open (unit=6,file='data/restart.dat',status='unknown')
      open (unit=7,file='data/bigice.dat',status='unknown')
      open (unit=15,file='out/model.out',status='unknown')
      open (unit=16,file='out/seasons.out',status='unknown')
      open (unit=17,file='out/co2clouds.out',status='unknown')
      open (unit=18,file='out/tempseries.out',status='unknown')
      open (unit=19,file='out/icelines.out',status='unknown')
      open (unit=20,file='out/dailytempseries.out',status='unknown')   
      open (unit=98,file='data/Altair_inc90_a3.3.txt',status='old')
      open (unit=99,file='data/insolaout.dat',status='old')
      if ( fillet ) then
        open (unit=50,file='out/fillet.out',status='unknown')
        open (unit=51,file='out/fillet_global.out',status='unknown')
      end if

c  WRITE OBLIQUITY TO OUTPUT
      write(15,2) 'OBLIQUITY: ', obl, ' degrees'
 2    format(/ a,f5.2,a)

c  WRITE CO2-CLOUD FILE HEADER
      write(17,3)
 3    format(58x,'IMCO2')
      write(17,4)
 4    format(2x,'solar dec.(deg)',2x,'-85',1x,'-75',1x,'-65',
     &  1x,'-55',1x,'-45',1x,'-35',1x,'-25',1x,'-15',2x,'-5',3x,'5',2x,
     &  '15',2x,'25',2x,'35',2x,'45',2x,'55',2x,'65',2x,'75',2x,'85',
     &  2x,'global cloud coverage' /)

c----------------------------------------------------------------------c
c SET UP INITIAL TEMPERATURE PROFILE

      ! set partial pressures
      if( radparam .eq. 1 ) then
        pco2 = fco2*pg0
        pn2 = pg0 - pco2
        ph2 = 0.0
        fh2 = 0.0
      else if ( radparam .eq. 2 ) then
        pco2 = pg0*fco2
        ph2 = pg0*(1-fco2)
        pn2 = 0.0
        fh2 = 1-fco2
      else if ( radparam .eq. 3 ) then

        if ( do_h2_cycle ) then
          pco2 = pg0*fco2
          ph2 = pg0*(1.-fco2)
          pn2 = 0.0
          fh2 = 1.-fco2
        else 
          pn2 = 1.0
          pco2 = pg0*fco2
          pg0 = pn2 + pco2
          !pco2 = pg0*fco2
          !ph2 = pg0*fh2
          !pn2 = pg0
        end if
        call radiation_init

      else if (radparam .eq. 4 ) then
        pco2 = pg0*fco2
        ph2 = pg0*(1.-fco2)
        pn2 = 0.0
      end if

      q0 = solarcon

      if( resfile .eq. 0 ) then
        do k = 1, nbelts, 1
          temp(k) = tempinit
        end do   
      else if( resfile .eq. 1 ) then !restart from last run
        rewind(6)
        do k = 1, nbelts, 1
           read(6,*) temp(k)
        end do
      else if( resfile .eq. 2 ) then !restart from present Earth
        rewind(4)
        do k = 1, nbelts, 1
           read(4,*) temp(k)
        end do
      else if( resfile .eq. 3 ) then !restart from hothouse
        rewind(5)
        do k = 1, nbelts, 1
           read(5,*) temp(k)
        end do
      else if( resfile .eq. 4 ) then !restart from big ice cap
        rewind(7)
        do k = 1, nbelts, 1
           read(7,*) temp(k)
        end do
      end if

c----------------------------------------------------------------------c
c  SET UP LATITUDINAL GRID (BELT BOUNDARIES)
      lat(0) = -pi/2
      latangle(0) = lat(0)*180./pi
      lat(nbelts+1) = pi/2
      latangle(nbelts+1) = lat(nbelts+1)*180./pi
      x(0) = sin(lat(0))
      x(nbelts+1) = sin(lat(nbelts+1))
      do 15 k = 1,nbelts,1
         lat(k) = lat(0) + k*pi/nbelts
         latangle(k) = lat(k)*180/pi
         cntrlat(k) = lat(0) + (1+2*(k-1))*pi/(2*nbelts)  !**belt centers
         cntrlatangle(k) = cntrlat(k)*180./pi             !**belt centers
 15   continue

c SET UP GEOGRAPHY
      call geog(igeog)


c Set diffusion coefficient
      diff(:) = d0
      !do k = 1, nbelts, 1
      !  diff(k) = d0
      !end do

c----------------------------------------------------------------------c
c  SET SOLAR CONSTANT, ECCENTRICITY, AND OBLIQUITY

      if( soladj ) then   
        rewind(99)
        do p = 1, niter, 1
           read(99,*) yrlabel(p),ecce(p),prec(p),obliq(p),solcon(p)
           solcon(p) = 4*solcon(p) * relsolcon    !for insolaout.dat
           prec(p)   = ASIN(prec(p)/ecce(p))*180/pi
        end do
      else if ( do_gough ) then
        do p = 1, niter, 1
          solcon(p) = relsolcon *
     &                  q0/(1 + 0.4*(1 - (p*yrstep/1.e9/4.57)))
          ecce(p)   = ecc
          prec(p)   = peri
          if ( do_marshist ) then ! do_oblvar
            !obliq(p)  = 25*sin(pi*p*yrstep/1.e6/250.) + 25.
            !obliq(p)  = 25*sin(pi*p*yrstep/1.e6/500.) + 25.
            !obliq(p)  = 25*sin(pi*p*yrstep/1.e6/400.) + 25.
            obliq(p)  = 25*sin(pi*p*yrstep/1.e6/450.) + 25.
          else
            obliq(p)  = obl
          end if
        end do
      else if ( do_futuresol ) then
        do p = 1, niter, 1
          solcon(p) = relsolcon *
     &                  q0/(1 - 0.38*((p-1)*yrstep/1.e9/4.55))
          ecce(p)   = ecc
          prec(p)   = peri
          obliq(p)  = obl
        end do
      else
        do p = 1, niter, 1
          solcon(p) = q0 * relsolcon
          ecce(p)   = ecc
          obliq(p)  = obl
          prec(p)   = peri
        end do
      end if

      if ( do_manualseasons ) then
        rewind(98)
        do p = 1, ndays, 1
           read(98,*) yrlabelD(p),ecceD(p),precD(p),
     &                obliqD(p),solconD(p)
           solconD(p) = solconD(p) * relsolcon       !for Altair_inc90_a3.3.txt
        end do
      end if

c----------------------------------------------------------------------c
c  SET UP TEMPERATURE GRID - TEMPERATURE IS CALCULATED AT BOTH POLES
c  AND BELT CENTERS. POLE TEMPERATURES ARE HELD EQUAL TO ADJACENT BELT
c  TEMPERATURES AS BOUNDARY CONDITIONS.
c
c **adjust latitude grid to belt centers   
      do 115 k = 1,nbelts,1
         lat(k) = cntrlat(k)
         latangle(k) = cntrlatangle(k)
         x(k) = sin(lat(k))
 115  continue

c **calculate grid spacing
      do 125 k = 0,nbelts,1
         dx(k) = abs(x(k+1) - x(k))
         stab(k)=((dx(k)**2.))/(2.*7.07*(510e12)*.2*50.*(1./.8368))
 125  continue

c  WRITE GRID DATA TO OUTPUT FILE

c      write(15,140) 
 140  format(/ / '**LATITUDE GRID DATA**')
c      write(15,142)
 142  format(10x,'x',6x,'dx',4x,'latitude(rad)',2x,'latitude(deg)')
      do 145 k = 0,nbelts+1,1
         !write(15,143) k,x(k),dx(k),lat(k),latangle(k),stab(k)
 143     format(3x,i2,2x,f6.3,2x,f6.3,5x,f6.3,8x,f5.1,8x,f6.3)
 145  continue
c
c----------------------------------------------------------------------c
c  READ FRESNEL REFLECTANCE TABLE
 146  format(2x,10(f5.3,2x))
      read(3,147) header
 147  format(/ a)
c      write(15,147) 'FRESNEL REFLECTANCE TABLE'
      do 149 i=1,9,1
         n1 = 10*(i-1) 
         n2 = n1+9
         read(3,148) (h20alb(n),n=n1,n2)
 148     format(17x,2p,10(f4.1,1x))
c         write(15,146) (h20alb(n),n=n1,n2)
 149  continue
      h20alb(90) = 1.00

c----------------------------------------------------------------------c
c  CALCULATE BELT AREAS - NORMALIZED TO PLANET SURFACE AREA
      do 180 k = 1,nbelts,1
         area(k) = abs(sin(lat(k) + pi/(2*nbelts)) -
     &     sin(lat(k) - pi/(2*nbelts)))/2.
 180  continue

c----------------------------------------------------------------------c
c  BEGIN INTEGRATION AT VERNAL EQUINOX. FIRST CALCULATE TIME SINCE
c  PERIHELION.

      d = d0  !**initialize diffusion coefficient
      w = (grav*msun)**(0.5)*a**(-1.5)  !2*pi/(orbital period)


c **write to 'co2clouds.out', at most, 1000 times.
      nwrite = amin0(1000,int((2*pi/w)/dt))
      twrite = (2*pi/w)/nwrite
 200  if(seasons) goto 230
      
      !if not seasons, then do mean annual insolation
      q = solcon(yricnt)
      goto 260

c  SEASONS
 230  e = acos((cos(prec(yricnt)*pi/180.)
     &    +ecce(yricnt))/(1.+cos(prec(yricnt)*pi/180.)*ecce(yricnt)))
      m = e - ecce(yricnt)*sin(e)   !**Kepler's Equation
      tperi = -m/w    !**tperi is the time of perihelion (vern.eqnx. is 0)

c      write(15,232)
 232  format(/ 'ORBITAL DATA')
c      write(15,233) 
 233  format(2x,'time (sec)',2x,'true anomaly (deg)',2x,
     & 'solar declination (deg)',2x,'orbital distance (cm)',2x,
     & 'irradiance (Wm^-2)')

 240  if (tcalc.lt.tend) then
         t = t + dt
         tcalc = tcalc + dt
         tlast = tlast + dt
         nt = nt + 1
      else
!         write(*,*) 'Calculation time has elapsed.'
!         goto 1000
      end if

      m = w*(t-tperi)
      if (m.gt.2*pi) m = m - 2*pi

      call keplereqn(m,ecce(yricnt),e)
      trueanom = acos((cos(e) - ecce(yricnt))/(1 - 
     &   ecce(yricnt)*cos(e)))

c  correct trueanomaly for pi < m < 2pi
      if (m .gt. pi) trueanom = 2*pi-trueanom

      r = a*(1-ecce(yricnt)*cos(e))
      q = solcon(yricnt) * (a/r)**2
      !q = q0 * relsolcon * (a/r)**2
      !q = q0 * relsolcon * (a0/r)**2
      thetawin = prec(yricnt)*pi/180 + 3*pi/2
      dec = asin(-sin(obliq(yricnt)*pi/180.)*cos(trueanom-thetawin))
      decangle = dec*180./pi
c      write(15,255) t,trueanom*180./pi,decangle,r,q
 255  format(2x,e9.3,6x,f7.3,14x,f7.3,16x,e11.6,13x,f9.3)

      goto 270      !**skip mean-annual stuff if doing seasonal model

c----------------------------------------------------------------------c
c MEAN-ANNUAL CALCULATION

 260  if (tcalc.lt.tend) then
         t = t + dt
         tcalc = tcalc + dt
         tlast = tlast + dt
         nt = nt + 1
      else
         write(*,*) 'Calculation time has elapsed.'
         goto 1000
      end if

 270  tempsum = 0.
      fluxsum = 0.
      albsum = 0.
      irsum = 0.
      globwthrate = 0.
      co2cldsum = 0.


c----------------------------------------------------------------------c
c  **THE BIG LOOP**
c----------------------------------------------------------------------c
c  FINITE DIFFERENCING - ATMOSPHERIC and OCEANIC ADVECTIVE HEATING 

      do 300 k=0,nbelts,1   !**first derivatives between grid points
         tprime(k) = (temp(k+1) - temp(k))/dx(k)
 300  continue

      do 310 k=1,nbelts,1   !**start belt loop 
                            !**first derivatives at grid points
         tprimeave(k) = (tprime(k)*dx(k) + tprime(k-1)*
     &      dx(k-1))/(dx(k) + dx(k-1))
                            !**second derivatives at grid points
         t2prime(k) = ((((dx(k)/2)**2)-((dx(k-1)/2)**2))*(1-x(k)**2)*
     &      tprimeave(k) + ((dx(k-1)/2)**2)*(1-(x(k)+(dx(k)/2))**2)*
     &      tprime(k) - (((dx(k)/2)**2)*(1-(x(k)-(dx(k-1)/2))**2))*
     &      tprime(k-1))/((dx(k)/2)*(dx(k-1)/2)*(dx(k)/2 + dx(k-1)/2)) 

c----------------------------------------------------------------------c
c  OUTGOING INFRARED (Obtained from fits to Jim Kasting's radiative
c    -convective model earthtem_ir. The fit can be applied for
c    190 < T < 370 K and 10^-5 < PCO2 < 10 bars with an rms error
c    of 4.59 Wm^-2.)                                 

      if ( linrad ) then

      !use a linear relationship for outgoing IR following North et al. 1981
      ir(k) = 203.3 + 2.09*(temp(k) - 273.15)

      else

      !use an interpolated relationship for outgoing IR
      if ( radparam .eq. 0 ) then 

        phi = log(pco2/3.3e-4)
     
      ! interpolation from Williams & Kasting (1997)
      ir(k) = -3.102011e-2 - 7.714727e-5*phi -
     &  3.547406e-4*phi**2 - 3.442973e-3*phi**3 -
     &  3.367567e-2*phi**4 - 2.794778*temp(k) - 
     &  3.244753e-3*phi*temp(k) + 2.229142e-3*phi**2*temp(k) + 
     &  9.173169e-3*phi**3*temp(k) - 1.631909e-4*phi**4*temp(k) + 
     &  2.212108e-2*temp(k)**2 + 3.088497e-5*phi*temp(k)**2 - 
     &  2.789815e-5*phi**2*temp(k)**2 - 7.775195e-5*phi**3*temp(k)**2 + 
     &  3.663871e-6*phi**4*temp(k)**2 - 3.361939e-5*temp(k)**3 - 
     &  1.679112e-7*phi*temp(k)**3 + 6.590999e-8*phi**2*temp(k)**3 + 
     &  1.528125e-7*phi**3*temp(k)**3 - 9.255646e-9*phi**4*temp(k)**3

      else if ( radparam .eq. 1 ) then
    
      phi = log10( pco2 )
      tmpk = log10( temp(k) )   

        ! Updated interpolation for CO2-dominated atmospheres (2015)
        term1=8.31425009785720803279e+00*tmpk**4+
     &  8.40111812576741634473e+00*tmpk**3*phi-
     &  7.20889378557585018825e+01*tmpk**3+
     &  9.27250707901223991669e-01*tmpk**2*phi**2-
     &  5.45101996969859712294e+01*tmpk**2*phi+
     &  2.34312071020079656591e+02*tmpk**2-
     &  1.72985690812803238892e-01*tmpk*phi**3-
     &  5.31819979906857742691e+00*tmpk*phi**2+
     &  1.15145551793545493524e+02*tmpk*phi

        term2=-3.36759529484533402410e+02*tmpk+
     &  1.03341023690065331869e-02*phi**4+
     &  4.62978730291089213278e-01*phi**3+
     &  7.43991998573408519491e+00*phi**2-
     &  7.87033076292617721492e+01*phi+
     &  1.84241516203608512114e+02

      ir(k) = 10**(term1 + term2) / 1000.

      else if ( radparam .eq. 2 ) then

      phi = log10( pg0 )
      tmpk = log10( temp(k) )   
        term1=4.69677278520229535275e+01*tmpk**4-
     &  8.52299915167291111118e+00*tmpk**3*(1.0-fh2)-
     &  3.65809180220086638258e+00*tmpk**3*phi-
     &  4.70385559455918155436e+02*tmpk**3-
     &  2.16216471567635526441e+00*tmpk**2*(1.0-fh2)**2-
     &  8.94315663928571424890e-01*tmpk**2*(1.0-fh2)*phi+
     &  6.18255931368686049154e+01*tmpk**2*(1.0-fh2)+
     &  2.80284049882792629660e+00*tmpk**2*phi**2+
     &  4.08137171809756296170e+01*tmpk**2*phi

        term2=+1.76947559145794366486e+03*tmpk**2-
     &  4.06932623054890107994e+00*tmpk*(1.0-fh2)**3+
     &  4.79116412732104701711e-01*tmpk*(1.0-fh2)**2*phi+
     &  2.01383744413982483934e+01*tmpk*(1.0-fh2)**2+
     &  2.61051880414713693979e-02*tmpk*(1.0-fh2)*phi**2+
     &  3.54072864661716391055e+00*tmpk*(1.0-fh2)*phi-
     &  1.57011107863713078814e+02*tmpk*(1.0-fh2)-
     &  3.18837733944702772515e-01*tmpk*phi**3-
     &  1.50232919353366050075e+01*tmpk*phi**2-
     &  1.34198870508522446698e+02*tmpk*phi

        term3=-2.95796874141556781979e+03*tmpk-
     &  7.75539632021813529761e-01*(1.0-fh2)**4+
     &  5.86922379610697264596e-01*(1.0-fh2)**3*phi+
     &  1.28924953407765237046e+01*(1.0-fh2)**3-
     &  2.43107734436077034534e-02*(1.0-fh2)**2*phi**2-
     &  2.52019655592552638268e+00*(1.0-fh2)**2*phi-
     &  4.02446385182495376398e+01*(1.0-fh2)**2-
     &  1.85054943621701661893e-02*(1.0-fh2)*phi**3-
     &  5.12824782188033845287e-02*(1.0-fh2)*phi**2-
     &  2.13570330172963895876e+00*(1.0-fh2)*phi

        term4=+1.41375290564461124632e+02*(1.0-fh2)+
     &  1.95159594394843996512e-02*phi**4+
     &  8.69351946778762929569e-01*phi**3+
     &  1.99432868566345007366e+01*phi**2+
     &  1.36829651234430457407e+02*phi+
     &  1.85540240841698141594e+03


      ir(k) = 10**(term1 + term2 + term3 + term4) / 1000.

      else if ( radparam .eq. 3 ) then

        olrval = 0.0
        call getOLR( fco2, temp(k), olrval ) 
        ir(k) = olrval / 1000.
        if ( ir(k) .le. -1. ) then
          print *, "radiation_mod: OLR solution unstable"
          print *, ir(k)
          stop
        end if

      else if ( radparam .eq. 4 ) then !bph new fits for CO2/H2

      X1 = Log10(temp(k))
      X2 = fh2
      X3 = Log10(pg0)      

        term1=1.01864613117570970644e+04*X1-
     &  1.11554429632794676763e+03*X2-
     &  4.06853660590399726971e+02*X3+
     &  1.51340773593768244609e+03*X1*X2+
     &  5.21925003735954646800e+02*X1*X3+
     &  7.14602011428987964337e+01*X2*X3+
     &  7.21839374617950511492e+02*X1*X2**2-
     &  6.80644630896995295188e+02*X1**2*X2+
     &  7.14769438781747510347e+01*X1*X2**3

        term2=-4.45170108948733544985e-01*X1*X3**2-
     &  2.22815355269315318765e+02*X1**2*X3+
     &  1.01412044533366866972e+02*X1**3*X2-
     &  4.39026854285687440083e-01*X1*X3**3+
     &  5.31784355134182895597e+00*X2*X3**2+
     &  3.16354045100431839899e+01*X1**3*X3-
     &  4.83090340477328794577e+00*X2**2*X3+
     &  1.14060709016931771664e-01*X2*X3**3-
     &  8.22622947819270677883e+00*X2**3*X3-
     &  6.19179994089639239974e+03*X1**2

        term3=+1.67115388935491182565e+03*X1**3-
     &  8.39924935826250589344e+02*X2**2-
     &  1.68929051690792618956e+02*X1**4-
     &  2.19433071347752303382e+02*X2**3+
     &  2.15128054099794452370e+00*X3**2+
     &  7.62663407275418876452e+01*X2**4+
     &  1.12522348042480802022e+00*X3**3+
     &  1.93311922625575838275e-02*X3**4-
     &  1.52902480259684608654e+02*X1**2*X2**2-
     &  1.91680176022660714308e-01*X1**2*X3**2

        term4=+2.92840400214665175227e-01*X2**2*X3**2-
     &  5.44825633272778162564e+01*X1*X2*X3-
     &  2.09292381681761252565e+00*X1*X2*X3**2+
     &  3.69021874454336895610e+00*X1*X2**2*X3+
     &  1.01745095615026599489e+01*X1**2*X2*X3-
     &  6.27515164298742547544e+03




      ir(k) = 10.**(term1 + term2 + term3 + term4) / 1000.      
c      ir(k) = log10(term1 + term2 + term3 + term4) / 1000.
      end if

        ir(k) = ir(k) - cloudir   !**reduction of outgoing-infrared by clouds

      end if
c      print *, X1,X2,X3,temp,ir
c      STOP

c----------------------------------------------------------------------c
c  DIURNALLY-AVERAGED ZENITH ANGLE                            

      if ((tan(dec)*tan(asin(x(k)))).ge.1.) then 
         h = pi
         mu(k) = x(k)*sin(dec) + cos(asin(x(k)))*cos(dec)*sin(h)/h
      else if ((tan(dec)*tan(asin(x(k)))).le.-1.) then 
         h = 0.
         mu(k) = 0.
      else
         h = acos(-x(k)*tan(dec)*((1-x(k)**2)**(-0.5)))
         mu(k) = x(k)*sin(dec) + cos(asin(x(k)))*cos(dec)*sin(h)/h
      end if

      !if ( do_longitudinal ) then  
      !  mu(k) = 1.0
      !end if

      z = acos(mu(k))*180./pi
      zrad = acos(mu(k))
                 
c----------------------------------------------------------------------c
c  HEAT CAPACITY and SURFACE ALBEDO                    
 
      if ( oceanalbconst ) then
        oceanalb = ocnalb
      else
        oceanalb = h20alb(int(z))
      end if

      ! Added this check so fractional ice cover would only 
      ! occur for 263 < T < 273.  (JDH)
      if (temp(k).ge.273.15) then
         fice(k) = 0
      else if (temp(k).lt.icetemp) then
         fice(k) = 1
      else
         fice(k) = 1. - exp((temp(k)-273.15)/10.)
      end if

      if ( cloudalb ) then
        acloud(k) = alpha + beta*zrad
      else
        acloud(k) = 0.0 
      end if

      if (temp(k).le.273.15) goto 420

      landalb = groundalb
      fice(k) = 0.
      fwthr = 1.       !**100% if T > 273K, 0% otherwise
      c(k) = focean(k)*cw + (1 - focean(k))*cl
      goto 430

c CO2 ICE ALBEDO IF CO2 IS CONDENSING AT SURFACE
 415  CALL SATCO2( temp(k), co2sat )
      if ( pco2 .gt. co2sat ) then
      snowalb = 0.33
      goto 420
      end if

c LAND WITH STABLE SNOW COVER; SEA-ICE WITH SNOW COVER
 !420  snowalb = 0.663 ! changed this from 0.7 as per Caldiera & Kasting (JDH) 
 !420  snowalb = 0.45 ! changed this from 0.7 as per Caldiera & Kasting (JDH) 
 420  landalb = snowalb*landsnowfrac + groundalb*(1 - landsnowfrac)
      icealb = snowalb

      c(k) = (1-fice(k))*focean(k)*cw + fice(k)*focean(k)*ci + 
     &  (1 - focean(k))*cl
      fwthr = 0.

      ! Added this conditional: if T < 273 then choose the max
      ! between the ice albedo and the cloud albedo (JDH)
 430  if (temp(k).ge.273.15) then         
         surfalb(k) = (1-fcloud)*((1-focean(k))*landalb +   
     &        focean(k)*((1-fice(k))*oceanalb + fice(k)*icealb)) + 
     &        fcloud*acloud(k)      
      else
         surfalb(k) = max(
     &        ((1-focean(k))*landalb +
     &        focean(k)*((1-fice(k))*oceanalb + fice(k)*icealb)),
     &        fcloud*acloud(k))
      end if      

c      if(.not.last) goto 450
c----------------------------------------------------------------------c
c  CARBONATE-SILICATE WEATHERING
 450  warea(k) = area(k)*(1-focean(k))*fwthr

      !carbonate-silicate cycle on long time scale, following Menou (2014)
      if ( temp(k) .ge. 273.15 ) then
!        wthrate(k) = warea(k)*weathering*((pco2/pco2soil)**betaexp)
!     &   *exp(kact*(temp(k)-288.))*(1+(krun*(temp(k)-288.)))**0.65
        wthrate(k) = warea(k)*weathering*((pco2/pco20)**betaexp)
     &   *exp(kact*(temp(k)-288.))*(1+(krun*(temp(k)-288.)))**0.65
      else
        wthrate(k) = 0.0
      end if

c----------------------------------------------------------------------c
c  TOP-OF-ATMOSPHERE ALBEDO (Obtained from fits to Jim Kasting's radiative
c    -convective model 'earthtem_alb.' Both fits can be applied for
c    10^-5 < pco2 < 10 bars, 0 < surfalb < 1, and 0 < z < 90 degrees
c    with r.m.s. errors (given in planetary-average incident solar
c    flux {340 W/m^-2}) of 7.58 and 4.66 Watts/m^2 for the low and
c    high temperature regimes respectively.

      as = surfalb(k)
c      as = .216   

      if ( linalb ) then 
        !set TOA albedo to surface albedo (JDH)
        atoa(k) = surfalb(k)
      else 

      if ( radparam .eq. 0 ) then
     
      ! albedo fits from Williams & Kasting (1997)
      if(temp(k).ge.280.) goto 520   !**goto high-temp fit

      atoa(k) = -6.891041e-1 + 1.046004*as + 7.323946e-2*as**2 - 
     &   2.889937e-1*mu(k)+2.012211e-1*as*mu(k)+8.121773e-2*mu(k)**2 -  
     &   2.837280e-3*pco2 - 3.741199e-2*as*pco2 - 
     &   6.349855e-3*mu(k)*pco2 + 6.581704e-4*pco2**2 + 
     &   7.805398e-3*temp(k) - 1.850840e-3*as*temp(k) + 
     &   1.364872e-4*mu(k)*temp(k) + 9.858050e-5*pco2*temp(k) - 
     &   1.655457e-5*temp(k)**2
      goto 530

 520  atoa(k) = 1.108210 + 1.517222*as + 7.588651e-2*as**2 - 
     &   1.867039e-1*mu(k)+2.098557e-1*as*mu(k)+6.329810e-2*mu(k)**2 +  
     &   1.970523e-2*pco2 - 3.135482e-2*as*pco2 - 
     &   1.021418e-2*mu(k)*pco2 - 4.132671e-4*pco2**2 - 
     &   5.79932e-3*temp(k) - 3.709826e-3*as*temp(k) - 
     &   1.133523e-4*mu(k)*temp(k) + 5.371405e-5*pco2*temp(k) + 
     &   9.269027e-6*temp(k)**2
 530  continue

      else if ( radparam .eq. 1 ) then

        phi = log10( pco2 )
        tmpk = log10( temp(k) )

        ! updated albedo fits for CO2/N2 atmospheres (2015)
        if ( temp(k) .le. 250 ) then

        term1 =-3.80922339639159390767e-01*mu(k)**3-
     &  6.86181617120208420246e-01*mu(k)**2*as+
     &  2.31793230316548221071e-01*mu(k)**2*tmpk-
     &  2.36223436336026049176e-02*mu(k)**2*phi+
     &  7.67301400772310682186e-01*mu(k)**2+
     &  6.28824306230756496783e-02*mu(k)*as**2-
     &  1.78348383597633913800e-01*mu(k)*as*tmpk-
     &  2.37093058356872832260e-02*mu(k)*as*phi+
     &  1.42993398589078402061e+00*mu(k)*as

        term2=-9.42874140148848183252e-01*mu(k)*tmpk**2+
     &  1.58597859799572450668e-02*mu(k)*tmpk*phi+
     &  4.17871254486884513568e+00*mu(k)*tmpk+
     &  2.93837265687994829769e-03*mu(k)*phi**2+
     &  9.72952667477233745785e-03*mu(k)*phi-
     &  5.97332227647180502572e+00*mu(k)+
     &  4.43312679710318527371e-02*as**3-
     &  5.01823197516588874467e-03*as**2*tmpk+
     &  2.66700143595150242215e-02*as**2*phi+
     &  1.87706863898554551784e-02*as**2

        term3=-3.17215890172350567511e+00*as*tmpk**2+
     &  1.55888179766147973171e-02*as*tmpk*phi+
     &  1.44250066694850413995e+01*as*tmpk-
     &  3.08351412489488475865e-02*as*phi**2-
     &  1.85426061631201261060e-01*as*phi-
     &  1.60883135091017592799e+01*as-
     &  1.17673787301924637205e+01*tmpk**3+
     &  2.22840757353930218887e-01*tmpk**2*phi+
     &  8.21723621055410262670e+01*tmpk**2+
     &  1.59651944857305064240e-02*tmpk*phi**2

        term4=-9.81711981544405865030e-01*tmpk*phi-
     &  1.91137300771063621596e+02*tmpk+
     &  2.31632261445296507019e-03*phi**3-
     &  1.05856783039563161208e-02*phi**2+
     &  1.14830209690318607585e+00*phi+
     &  1.48662609849158229736e+02

        else

        term1 =-5.04232716745127151903e-01*mu(k)**3-
     &  3.33831975323527430088e-01*mu(k)**2*as+
     &  1.10989674004723704037e+00*mu(k)**2*tmpk-
     &  5.12985634933501161159e-02*mu(k)**2*phi-
     &  1.37002953772578139890e+00*mu(k)**2+
     &  7.42133599320679293587e-02*mu(k)*as**2-
     &  1.54926696559962251420e+00*mu(k)*as*tmpk-
     &  3.65656241697920325606e-02*mu(k)*as*phi+
     &  4.34163847358376830954e+00*mu(k)*as

        term2=+2.66525684665264384066e+00*mu(k)*tmpk**2-
     &  4.76404344529663126284e-02*mu(k)*tmpk*phi-
     &  1.39651863088575627359e+01*mu(k)*tmpk+
     &  6.66267692133716647046e-03*mu(k)*phi**2+
     &  2.08349173453556052449e-01*mu(k)*phi+
     &  1.69542841293293449212e+01*mu(k)+
     &  7.46932461848770906654e-02*as**3-
     &  1.83677092125004437495e-01*as**2*tmpk+
     &  2.57348239818812260515e-02*as**2*phi+
     &  3.93925768617288063478e-01*as**2

        term3=-3.95929907396172842127e+00*as*tmpk**2+
     &  3.29522425670107643736e-01*as*tmpk*phi+
     &  1.84045452486614244947e+01*as*tmpk-
     &  3.13359092463711685905e-02*as*phi**2-
     &  9.42011930045157264146e-01*as*phi-
     &  2.10503344700073924400e+01*as+
     &  3.29676405480184655516e+01*tmpk**3+
     &  7.90783265312129834967e-01*tmpk**2*phi-
     &  2.46889915210692521441e+02*tmpk**2+
     &  8.30108618133403586281e-02*tmpk*phi**2

        term4=-3.72584714043850429022e+00*tmpk*phi+
     &  6.15857421042532109823e+02*tmpk+
     &  4.09373012948141025424e-03*phi**3-
     &  1.69079209243857525591e-01*phi**2+
     &  4.45447135967024365755e+00*phi-
     &  5.11168467178811397389e+02

        end if

        atoa(k) = term1 + term2 + term3 + term4 

      else if ( radparam .eq. 2 ) then

        phi = log10( pg0 )
        tmpk = log10( temp(k) )

        ! updated albedo fits for CO2/H2 atmospheres (2015)
        if ( temp(k) .le. 250 ) then

        term1 =-3.41888251549252286665e-01*mu(k)**3-
     &  7.90177446947207928751e-01*mu(k)**2*as-
     &  2.72262737687278599807e-01*mu(k)**2*tmpk+
     &  1.03466725555122224939e-02*mu(k)**2*(1.0-fh2)+
     &  8.50029147631192055767e-03*mu(k)**2*phi+
     &  1.93212037131758607167e+00*mu(k)**2+
     &  4.14893171933844936983e-02*mu(k)*as**2-
     &  3.07137890718480709162e-02*mu(k)*as*tmpk+
     &  1.14792944390371322999e-02*mu(k)*as*(1.0-fh2)

        term2=-5.54314356783364436954e-02*mu(k)*as*phi+
     &  1.17466585092971143034e+00*mu(k)*as+
     &  2.88632104728009819539e+00*mu(k)*tmpk**2-
     &  3.19654435462561198333e-02*mu(k)*tmpk*(1.0-fh2)-
     &  3.92811662257368909845e-01*mu(k)*tmpk*phi-
     &  1.35814776466665882992e+01*mu(k)*tmpk-
     &  3.58803689291265458239e-03*mu(k)*(1.0-fh2)**2+
     &  1.15604667327326630466e-04*mu(k)*(1.0-fh2)*phi+
     &  6.58513327478918025770e-02*mu(k)*(1.0-fh2)+
     &  1.54683743964890690198e-02*mu(k)*phi**2

        term3=+9.58075324811706985351e-01*mu(k)*phi+
     &  1.46116232993565215992e+01*mu(k)+
     &  3.59418044120485766224e-02*as**3+
     &  6.11902687210621396008e-02*as**2*tmpk+
     &  1.07831449663317829706e-03*as**2*(1.0-fh2)+
     &  5.41411219832318574285e-02*as**2*phi-
     &  1.02186234477593285153e-01*as**2-
     &  3.19054737071549521232e+00*as*tmpk**2+
     &  8.11189374199748208794e-02*as*tmpk*(1.0-fh2)-
     &  8.41274612039302413513e-03*as*tmpk*phi

        term4=+1.42668358155923016284e+01*as*tmpk-
     &  1.52636953777691378870e-02*as*(1.0-fh2)**2-
     &  2.00486768715006907077e-03*as*(1.0-fh2)*phi-
     &  1.73661142630655435104e-01*as*(1.0-fh2)-
     &  5.15509475686814835904e-02*as*phi**2-
     &  2.33460669190694719566e-01*as*phi-
     &  1.56829735069136404491e+01*as-
     &  2.01079212886263860582e+01*tmpk**3+
     &  7.67597139570774911199e-01*tmpk**2*(1.0-fh2)+
     &  1.17365804685178165556e+00*tmpk**2*phi

        term5=+1.38294517497544603657e+02*tmpk**2-
     &  8.58907463110849672683e-02*tmpk*(1.0-fh2)**2-
     &  6.28204459275867321821e-02*tmpk*(1.0-fh2)*phi-
     &  3.52960104567058063907e+00*tmpk*(1.0-fh2)-
     &  5.73195603682496515607e-02*tmpk*phi**2-
     &  5.29024591649556086281e+00*tmpk*phi-
     &  3.16603696488367347683e+02*tmpk+
     &  1.87184850687503318012e-02*(1.0-fh2)**3-
     &  4.99139908059316201455e-03*(1.0-fh2)**2*phi+
     &  1.58858270898283582273e-01*(1.0-fh2)**2

        term6=+2.03127920236113651206e-03*(1.0-fh2)*phi**2+
     &  1.59315444513521503600e-01*(1.0-fh2)*phi+
     &  4.08859633431839153417e+00*(1.0-fh2)+
     &  4.92942902461866814395e-03*phi**3+
     &  1.76050957604708957494e-01*phi**2+
     &  6.08822715747315967860e+00*phi+
     &  2.41835068101150113762e+02


        else

        term1 =-3.70851269169879815824e-01*mu(k)**3-
     &  3.45445707097643828209e-01*mu(k)**2*as+
     &  7.75813106851914513484e-01*mu(k)**2*tmpk-
     &  4.37926598821664850997e-03*mu(k)**2*(1.0-fh2)+
     &  8.30907012476914466625e-03*mu(k)**2*phi-
     &  8.08121898461139132053e-01*mu(k)**2+
     &  6.19434328882089643709e-02*mu(k)*as**2-
     &  1.86910174613248014630e+00*mu(k)*as*tmpk+
     &  2.97484543244541030371e-02*mu(k)*as*(1.0-fh2)

        term2=-6.47159644790028937278e-02*mu(k)*as*phi+
     &  5.09685365127381206918e+00*mu(k)*as-
     &  2.22139207043435993327e+00*mu(k)*tmpk**2+
     &  8.04155552223042024984e-02*mu(k)*tmpk*(1.0-fh2)+
     &  1.62894211373932590314e-01*mu(k)*tmpk*phi+
     &  1.08147183735960812356e+01*mu(k)*tmpk-
     &  5.47050254503252987581e-03*mu(k)*(1.0-fh2)**2-
     &  6.47303913823305855646e-03*mu(k)*(1.0-fh2)*phi-
     &  2.01172121790839208977e-01*mu(k)*(1.0-fh2)+
     &  1.68709677231692820043e-02*mu(k)*phi**2

        term3=-3.73870273052445190043e-01*mu(k)*phi-
     &  1.42527376906675655732e+01*mu(k)+
     &  8.02788220193820512005e-02*as**3-
     &  5.72715994416816420731e-02*as**2*tmpk+
     &  2.26901395282027001227e-03*as**2*(1.0-fh2)+
     &  2.71046463949741710253e-02*as**2*phi+
     &  5.70829781212874656782e-02*as**2-
     &  6.10825133251704155413e+00*as*tmpk**2+
     &  4.07976543429671101304e-01*as*tmpk*(1.0-fh2)+
     &  6.95094459838085376724e-01*as*tmpk*phi

        term4=+2.87052467136008431225e+01*as*tmpk-
     &  4.65629570929901731580e-02*as*(1.0-fh2)**2-
     &  1.55652000112216845618e-02*as*(1.0-fh2)*phi-
     &  9.29906450324752942294e-01*as*(1.0-fh2)-
     &  3.30578735873280696311e-02*as*phi**2-
     &  1.85654436072478024045e+00*as*phi-
     &  3.34271433377927493780e+01*as+
     &  4.19075909323097306469e+01*tmpk**3-
     &  2.41437621421361292562e+00*tmpk**2*(1.0-fh2)-
     &  3.65931088974145612980e-01*tmpk**2*phi

        term5=-3.05463565653780563025e+02*tmpk**2+
     &  2.36004773312411507413e-02*tmpk*(1.0-fh2)**2-
     &  6.45897400020475354054e-02*tmpk*(1.0-fh2)*phi+
     &  1.16961025823086650632e+01*tmpk*(1.0-fh2)+
     &  1.46636590367962948989e-01*tmpk*phi**2+
     &  1.81985447344254880342e+00*tmpk*phi+
     &  7.41386912374803728198e+02*tmpk+
     &  4.92589810954733656190e-02*(1.0-fh2)**3+
     &  5.17074814273357774574e-03*(1.0-fh2)**2*phi-
     &  1.43439261962234887449e-01*(1.0-fh2)**2

        term6=-3.42234657730983984555e-03*(1.0-fh2)*phi**2+
     &  1.48010570194314516890e-01*(1.0-fh2)*phi-
     &  1.41040031286142664158e+01*(1.0-fh2)+
     &  3.78842516310583807618e-03*phi**3-
     &  3.24961066540855103568e-01*phi**2-
     &  2.11874640435505945391e+00*phi-
     &  5.98624493625917807549e+02


        end if

        atoa(k) = term1 + term2 + term3 + term4 + term5 + term6

      else if ( radparam .eq. 3 ) then
        zendeg = mu(k)*180/pi

        call getPALB( fco2, temp(k), zendeg,        
     &                surfalb(k), atoa(k) )

        if ( atoa(k) .le. -1 ) then
          print *, "radiation_mod: PALB solution unstable"
          stop
        end if

      else if ( radparam .eq. 4 ) then !bph new fits for CO2/H2

      X1 = mu(k)   ! [Don’t forget that the zenith angle should be inradians when using trig functions]
      X2 = surfalb(k)
      X3 = Log10(temp(k))
      X4 = fh2
      X5 = Log10(pg0)

      if ( temp(k) .le. 270. ) then

           term1 =0.80296306939259698420841004917747*X4
     &    -24.150741094789406560039424221031*X2
     &    -61.065201590327610858821572037414*X3
     &    -5.3782272788744664993032529309858*X1
     &    +1.6437580835326983663691180481692*X5
     &    +1.5831092500738721540187725622673*X1*X2
     &    +3.6197193177596478363966525648721*X1*X3
     &    +0.13963248046409618141616704178887*X1*X4
     &    +21.193937158480956384210003307089*X2*X3

        term2=+0.76768756182805653054401773260906*X1*X5  -
     &  0.25077715461460736712240304768784*X2*X4
     &    -0.095959108873828663499594426866679*X2*X5
     &    -0.46529584482769553721581701211107*X3*X4
     &    -1.6450982916493717134187591000227*X3*X5
     &    +0.026848546460768573512778445433469*X4*X5
     &    +0.12137525453631232974505138599852*X1*X2**2
     &    -0.49309418199260346415968569999677*X1**2*X2
     &    -0.75632099741337199816371139604598*X1*X3**2
     &    +0.094428872883463313425700391690043*X1**2*X3

        term3=-0.014855881492708888913512055296451*X1*X4**2  -
     &  4.5827299423053489135782001540065*X2*X3**2
     &    -0.0012115354318920068792814204172714*X1**2*X4
     &    -0.052965960792519910238507208077863*X2**2*X3
     &    +0.01804858048660996222567476365839*X1*X5**2
     &    +0.10603464766816737829824290884062*X2*X4**2
     &    -0.040387612533265754288791526960267*X1**2*X5
     &    +0.0013613095000860285770788449255519*X2**2*X4
     &    -0.0627993104179862954294932819721*X2*X5**2
     &    +0.31535992165419451938035422244866*X3*X4**2

        term4=+0.084658641118897079436145247655077*X2**2*X5  +
     &  0.057695384295692299370106326250607*X3**2*X4
     &    +0.050186317445398531678080189522007*X3*X5**2
     &    +0.43478919151567818612846849646303*X3**2*X5
     &    +0.0032878546504428078980919458729204*X4*X5**2
     &    -0.018375797562613351010529783025049*X4**2*X5
     &    +0.70833794354776935531248227562173*X1**2
     &    -0.26379782951504221477634359871445*X1**3
     &    -0.020566650420521317976785269365791*X2**2
     &    +0.14704301075048747948770255788986*X2**3

        term5=+27.513796527502233146833532373421*X3**2  -
     &  4.1100064964847300430506038537715*X3**3
     &    -0.83485730627398302416963815630879*X4**2
     &    -0.010262784778034033522642332059149*X4**3
     &    -0.070200084614304783525717823522427*X5**2
     &    +0.003016187703942105596127687405783*X5**3
     &    -0.37054493943301297642634040130361*X1*X2*X3
     &    -0.0060270300344690318564144781987579*X1*X2*X4
     &    -0.074000559250107830755460724958539*X1*X2*X5
     &    -0.054771945860188242516386480929214*X1*X3*X4

        term6=-0.27897964906680727814602960279444*X1*X3*X5  +
     &  0.093187713827676016409284898145415*X2*X3*X4
     &    +0.0033342553127215395804283648573119*X1*X4*X5
     &    -0.092234013887297169875800761928986*X2*X3*X5
     &    -0.0024142092841804457283017004698422*X2*X4*X5
     &    -0.0054167410669458786187657750588187*X3*X4*X5
     &    +45.513033481175355632331047672778

      else

           term1 =5.0267861572060361652347637573257*X1
     &    -49.779706023277014992345357313752*X2
     &    +598.85651989270729700365336611867*X3
     &    +3.3926239759534917439509627001826*X4
     &    -4.2322551607841862875147853628732*X5
     &    +4.574581396833711899319041549461*X1*X2
     &    -4.463113834055816298018726229202*X1*X3
     &    +0.22624891432245175515447499492439*X1*X4
     &    +41.57155382677778732158913044259*X2*X3

        term2=-0.20526532347849119064520095889748*X1*X5  +
     &  0.97538836508197179764323436756968*X2*X4
     &    -2.0234707561887592675020641763695*X2*X5
     &    -2.8715303418845814498183699470246*X3*X4
     &    +3.7954083281352213319337352004368*X3*X5
     &    -0.0070211817023723326405826306029212*X4*X5
     &    +0.11589347467574391548073009516884*X1*X2**2
     &    -0.26363811728246194787672607162676*X1**2*X2
     &    +0.81950724181666811940516481627128*X1*X3**2
     &    +0.86113814605079042063806582518737*X1**2*X3

        term3=+0.015793743209526141940690635578903*X1*X4**2  -
     &  8.6257903394618402614923979854211*X2*X3**2
     &    +0.0064596687689501171292238979049216*X1**2*X4
     &    -0.08653781973037230623724269662489*X2**2*X3
     &    +0.016503398844075947010656690849828*X1*X5**2
     &    +0.040530059545890124994560466120674*X2*X4**2
     &    -0.0057978371881241066998491362483037*X1**2*X5
     &    -0.0040854762053125675061449051383988*X2**2*X4
     &    -0.051633464212411959992632404237156*X2*X5**2
     &    -0.10328010320951024059166201141124*X3*X4**2

        term4=+0.047197403614696556428675933148043*X2**2*X5  +
     &  0.60788329293614473680662513288553*X3**2*X4
     &    +0.10859389477605767282408777418823*X3*X5**2
     &    -0.81056731709101870197997641298571*X3**2*X5
     &    +0.0054025109908369299896446769082559*X4*X5**2
     &    -0.014388580103157301917637944654871*X4**2*X5
     &    -1.1619007049553844446165840054164*X1**2
     &    -0.33604612466494437894581892578572*X1**3
     &    +0.051431131875346092208189929806395*X2**2
     &    +0.14545992943627544802254192291002*X2**3

        term5=-243.79506982643917467612482141703*X3**2  +
     &  33.06608079382326792483581812121*X3**3
     &    +0.19629279658223547366802108626871*X4**2
     &    +0.0053382547639579999015913003290734*X4**3
     &    -0.21064302982024196708721319737379*X5**2
     &    +0.0062048147912155751648732859848678*X5**3
     &    -1.6874099652687377659532330653747*X1*X2*X3
     &    -0.023197237403772792002643043929311*X1*X2*X4
     &    -0.067444400354935873931339074260904*X1*X2*X5
     &    -0.092197958264046600551644417009811*X1*X3*X4

        term6=+0.10243131958860483832207677323822*X1*X3*X5  -
     &  0.4002046194796735067988890932611*X2*X3*X4
     &    +0.0011372380330880419059547126181542*X1*X4*X5
     &    +0.72202066774105577717790538372356*X2*X3*X5
     &    +0.011860560688605140816642702361605*X2*X4*X5
     &    +0.0080176660219647769012318860859523*X3*X4*X5
     &    -489.48662878273108844950911588967

      

      end if
       
      atoa(k) = term1 + term2 + term3 + term4 + term5 + term6
c      print *, ann_albave
c      print *, term1, term2, term3
c      STOP
      end if

      end if
c----------------------------------------------------------------------c
c  DIURNALLY-AVERAGED INSOLATION 
      
      if ( do_manualseasons ) then
        q = solconD( daynum )
        if ( k .eq. nbelts ) then
          daynum = daynum + 1
          if ( mod( ndays, nt ) .eq. 0 ) daynum = 1
          !print *, "day = ", daynum
        end if
      end if

      if ( do_longitudinal ) then
        if ( x(k) .ge. 0.0 ) then
          s(k) = (q/pi) * x(k)
        else
           s(k) = 0.0
           ir(k) = ir(k) + cloudir
        end if 

        !if ((x(k) .le. sin(pi/4)) .and. (x(k) .ge. sin(-pi/4))) then 
        !  s(k) = (q/pi) * cos( 2*asin(x(k)) )
        !else
        !  s(k) = 0.0
        !  ir(k) = ir(k) + cloudir
        !end if
      else if (seasons) then
        s(k) = (q/pi)*(x(k)*sin(dec)*h + 
     &        cos(asin(x(k)))*cos(dec)*sin(h))
      else 
        ! for mean annual use the second Legendre polynomial
        p2   = (3*sin(lat(k))*sin(lat(k)) - 1) / 2
        s(k) = (q/pi)*(1 - 0.5*p2)
      end if

c----------------------------------------------------------------------c
c  SURFACE TEMPERATURE - SOLVE ENERGY-BALANCE EQUATION 

      if ( constheatcap ) then
         c(k) = heatcap
         !c(k) = 5.25e6 * 40
      end if

      if( do_stochastic ) then
        temp(k) = (diff(k)*t2prime(k)-ir(k)+s(k)*(1-atoa(k)))*dt/c(k) 
     &            + sqrt(noisevar)*boxmuller()
     &            + temp(k)
      else
        temp(k) = (diff(k)*t2prime(k)-ir(k)+s(k)*(1-atoa(k)))*dt/c(k) 
     &            + temp(k)
      end if

c      if ( temp(k) .le. 195.0 ) then    !changed to work with database
c        temp(k) = 195.0
c      end if


      imco2(k) = 0.
      co2flag  = 0

c-nb----------------------------------------------------
c-nb   CO2 Saturation Temperature 
c-nb   Set minimum temperature based on saturation temperature of CO2 
 
      CALL SATCO2( temp(k), co2sat )
      
      if ( pco2 .gt. co2sat ) then
  
        imco2(k) = 1.
        co2flag  = 1
        temp(k) = 3.1823*log10(pco2)**3 + 10.5165*log10(pco2)**2 + 
     &           28.5760*log10(pco2) + 192.2084 
c      print *, "CO2 is condensing!"
c      print *, temp, pco2, co2sat
      end if

!      if ((yrcnt.gt.134.6e6).and.(yrcnt.lt.135e6)) then
!      write(*,1235)(temp(iii),iii=1,18)
!      print *,'second time'
!      pause
!1235  format(f8.2)
!      endif

c  SUM FOR GLOBAL AVERAGING
      irsum = irsum + area(k)*ir(k)
      fluxsum = fluxsum + area(k)*s(k)
      albsum = albsum + area(k)*s(k)*atoa(k)
      tempsum = tempsum + area(k)*temp(k)
      globwthrate = globwthrate + wthrate(k)
      co2cldsum = co2cldsum + area(k)*imco2(k)*s(k)

      if(.not.last) goto 310
c  ZONAL STATISTICS - if last loop 
      zntempmin(k) = amin1(zntempmin(k),temp(k))
      if (zntempmin(k).eq.temp(k)) zndecmin(k) = decangle
      zntempmax(k) = amax1(zntempmax(k),temp(k))
      if (zntempmax(k).eq.temp(k)) zndecmax(k) = decangle
      zntempsum(k) = zntempsum(k) + temp(k)
      znalbsum(k)  = znalbsum(k) + atoa(k)
      znolrsum(k)  = znolrsum(k) + ir(k)
      znsurfalbsum(k) = znsurfalbsum(k) + surfalb(k)

 310  continue                                     !**end of belt loop


c  **set pole temps equal to adjacent belt temps
      temp(0) = temp(1)
      temp(nbelts+1) = temp(nbelts)

      ! daily output could get large, set do_dailyoutput = .false. for long integrations
      if ( do_dailyoutput ) then
      if( .not. last) then
        write(20,609) t/dt+366*yrcnt,q,temp(1),temp(5),temp(9),pg0,
     &    pco2,fh2,co2flag,d
 609    format(2x,f12.2,f12.2,f12.2,f12.2,2x,f8.3,2x,e12.3,2x,e12.6,
     &         2x,f8.5,2x,i12,2x,f8.5)
! 609    format(2x,f12.2,f12.2,f12.2,f12.2,2x,f8.3,2x,e12.3,2x,e12.6,
!     &         2x,i12,2x,f8.5,2x,f8.5)
        !write(20,609) t/dt+366*yrcnt,tempsum,pg0,pco2,co2flag,fh2,d
! 609    format(2x,f12.2,2x,f8.3,2x,e12.3,2x,e12.6,2x,i12,
!     &    2x,f8.5,2x,f8.5)
      end if
      end if


c  WRITE OUTPUT FILES - ONLY ON LAST LOOP      
      if(.not.last ) goto 710
      if(tlast.lt.twrite) goto 710
      tlast = 0.

c  EQUINOXES AND SOLSTICES

      if (seasons) then
         decold = decmid
         decmid = decnew
         decnew = decangle

         if (((abs(decmid).gt.(abs(decold))).and.
     &    (abs(decmid).gt.(abs(decnew)))).or.(decnew*decmid.le.0.)) then 
     
c     write(15,610) file(nfile),decangle
 610        format('data written to file ',a8, ' at declination ',f6.2,
     &           ' degrees')
 615           format(f6.2,5(3x,f7.3))
 620        continue
            nfile = nfile + 1
         end if         
      end if

c  CO2-CLOUD DATA
      write(17,630) decangle,imco2,co2cldsum/fluxsum
 630  format(3x,f6.2,9x,18(3x,i1),9x,f5.3)
     
c----------------------------------------------------------------------c
c  GLOBAL AVERAGING
     
 710  irave = irsum
      fluxave = fluxsum
      co2cldave = co2cldsum/fluxave
      albave = albsum/fluxave
      tempave = tempsum

      iravesum = iravesum + irave
      fluxavesum = fluxavesum + fluxave
      albavesum = albavesum + albave
      tempavesum = tempavesum + tempave
      wthratesum = wthratesum + globwthrate
      co2cldavesum = co2cldavesum + co2cldave
      nstep = nstep + 1

c  SEASONAL AVERAGING
      if(t.lt.2*pi/w) goto 800    !**one orbit since last averaging
      ann_tempave = tempavesum/nstep
      ann_albave = albavesum/nstep
      ann_fluxave = fluxavesum/nstep
      ann_irave = iravesum/nstep
      ann_wthrateave = wthratesum/nstep
      ann_co2cldave = co2cldavesum/nstep

      if(.not.last) then
!        write(18,711) yrcnt,maxval(temp),minval(temp),pg0,pco2,co2flag,
!     &    fh2,d
!        write(18,711) yrcnt,maxval(temp),minval(temp),pg0,pco2,co2flag,
!     &    q,d
        write(18,711) yrcnt,ann_tempave,pg0,pco2,pco2soil,
     &    gammaout,q,d
 711    format(2x,i12,2x,f8.3,2x,f8.3,2x,e12.3,2x,e12.3,2x,f5.2,
     &    2x,f8.2,2x,f8.5)
! 711    format(2x,i12,2x,f8.3,2x,f8.3,2x,e12.3,2x,e12.3,2x,i12,
!     &    2x,f8.5,2x,f8.5)

c      write(19,712) pco2,zntempave(1:nbelts)
c 712  format(2x,f12.3,20(2x,f8.3))

c        write(19,712) pco2,temp
 712    format(2x,f12.3,20(2x,f8.3))
      
      if ((yrcnt.gt.126e6).and.(yrcnt.lt.192e6)) then
!      if ((yrcnt.gt.134.6e6).and.(yrcnt.lt.135e6)) then
c print out latitude dependent tempreatures 
c            open(unit=765,file= 'out/LATdependence.out')

c            do 996 k =1,nbelts,1 
c                  write(765,9997) temp(k), pco2
c996         continue
9997        format(3(3x,1P6E14.6))
      end if

c

      end if

      !H2 cycle: calculate ncolh2 and update based on outgassing rate
      if ( do_h2_cycle ) then
        h2escape = 1.6e13 * fh2 / ( 1 + fh2 )  !cgs units: 1/(cm^2 s)
        ph2 = (2./44.*pco2 + 2./28.*pn2) * (fh2 / (1.-fh2))
        ncolh2 = ph2 / ( 3.3e-30 * 373. )  !cgs units
        ncolh2 = ncolh2 + (yrstep*3.1557e7)*(h2outgas-h2escape)
      end if

c     !Carbonate-Silicate pCO2 update
      if ( do_cs_cycle ) then
        pco2 = pco2 + (yrstep/1.e9)*(outgassing-sum(wthrate))
        if ( pco2 .lt. 0 ) then
          pco2 = 1.e-12
        end if

        ! adjust pg0
        if ( radparam .eq. 1 ) then

          pg0 = pco2 + pn2

        else if ( radparam .eq. 3 ) then

          pg0  = pco2 + pn2
          fco2 = pco2 / pg0

        else if ( radparam .eq. 4 ) then

          !H2 cycle: calculate new value of fh2
          if ( do_h2_cycle ) then
            ph2 = ncolh2 * ( 3.3e-30 * 373. )  !cgs units
            fh2 = ph2 / ( ph2 + 2./44.*pco2 + 2./28.*pn2 )
            if ( fh2 .lt. 0 ) then
              fh2 = 0.0
            else if ( fh2 .gt. 1 ) then
              fh2 = 1.0
            end if
          end if

c          pg0 = ( pco2 + pn2 ) / ( 1 - fh2 )
          pg0 = pco2 + pco2*2.0/44.0*fh2/(1-fh2)
        end if

        !print *, pg0, pco2, fh2, ann_tempave
        
      end if

c     !Biological Productivity update to soil pCO2
      if ( do_bioprod ) then

        pmin = 1.e-5  !10ppm co2 limit
        !pmin = 1.e-6  !1ppm co2 limit

        phalf = 2.108e-4  !assumes 10ppm limit, see Caldeira & Kasting (1992) for calculation method
        !phalf = 2.1692e-4  !assumes 1ppm limit       

        gamma = 2*(1-((ann_tempave-tempinit)-25)/25)**2 *
     &           ((pco2-pmin)/(phalf+(pco2-pmin)))
        if ( yricnt .eq. 1 ) then
          gamma0 = gamma
        end if
        gammaout = gamma / gamma0

        pco2soil = pco2soil0*gamma*(1-(pco20/pco2soil0)) + pco2

      end if


c  ADJUST DIFFUSION COEFFICIENT
      avemol = mp*(28.0*pn2+44.*pco2 + 2.0*ph2)/(pg0) 
      hcp = (hcpn2*pn2 + pco2*hcpco2 + hcph2*ph2)/(pg0)

      if ( diffadj ) then
      d = d0*(pg0)*((avemol0/avemol)**2)*(hcp/hcp0)*  
     &   (rot0/rot)**2
!        print *, hcp 
c-nb  
      do k = 1, nbelts, 1
          diff(k) = d
      end do 

      end if
          
      if(.not.last) goto 790  
c ZONAL SEASONAL AVERAGES 
      shtempave = 0
      nhtempave = 0
      do 732 k = 1, nbelts, 1
         zntempave(k) = zntempsum(k) / nstep
         znalbave(k) = znalbsum(k) / nstep
         znolrave(k) = znolrsum(k) / nstep
         znsurfalbave(k) = znsurfalbsum(k) / nstep
         if ( k .le. nbelts/2 ) then
           shtempave = shtempave + zntempave(k)
         else
           nhtempave = nhtempave + zntempave(k)
         end if     
 732  continue
      shtempave    = shtempave / (nbelts/2)
      nhtempave    = nhtempave / (nbelts/2)
      print *, "SH/NH temperature difference = ", shtempave - nhtempave
      zntempave(nbelts+1) = zntempave(nbelts)  !**for ice-line calculation

c  FIND ICE-LINES (ANNUAL-AVERAGE TEMP < icetemp [=263.15K])
      nedge = 0
      do 740 k=1,nbelts,1
      if ((zntempave(k+1)-icetemp)*(zntempave(k)-icetemp) .lt. 0.) then
         icelat = latangle(k) + ((latangle(k+1)-latangle(k))/
     &   (zntempave(k+1)-zntempave(k)))*(icetemp-zntempave(k))
         nedge = nedge + 1
         iceline(nedge) = icelat
      end if
 740  continue

      do 750 k=1,nbelts,1
         write(15,751) latangle(k),zntempave(k),zntempmin(k),
     &      zndecmin(k),zntempmax(k),zndecmax(k),
     &      znalbave(k)
c-nb     &      (zntempmax(k)-zntempmin(k))/2.
 751     format(4x,f4.0,9x,f8.3,5x,f8.3,5x,f6.2,5x,f8.3,5x,f6.2,5x,f8.3)
         write(16,752) latangle(k),(zntempmax(k)-zntempmin(k))/2.,
     &     10.8*abs(x(k)),15.5*abs(x(k)),6.2*abs(x(k))
 752     format(4x,f4.0,4x,4(3x,f8.3))
         write(6,753) zntempave(k)
 753     format(2x,f8.3)
         if ( fillet ) then
           write(50,754) latangle(k), zntempave(k), znsurfalbave(k), 
     &         znalbave(k), znolrave(k)
 754       format(f5.1,1x,f6.2,1x,f4.2,1x,f4.2,1x,f6.2)
         end if
 750  continue
      
 755  format(/ 'SURFACE DATA')
 756  format(2x,'latitude(deg)',2x,'temp(k)',2x,'belt area',
     &  2x,'weathering area',2x,'zonal weathering rate (g/yr)')
       write(15,760)
 760  format(/ 'ICE LINES (Tave = 263.15K)')
      if((nedge.eq.0).and.(zntempave(nbelts/2).le.icetemp)) then
        icelineN = 0.0
        icelineS = 0.0
        icelineNMax = 90.0
        icelineNMin = 0.0
        icelineSMax = 0.0
        icelineSMin = -90.0
      	write(15,*) '  planet is an ice-ball.' 
      else if((nedge.eq.0).and.(zntempave(nbelts/2).gt.icetemp)) then
        icelineN = 90.0
        icelineS = -90.0
        icelineNMax = 90.0
        icelineNMin = 90.0
        icelineSMax = -90.0
        icelineSMin = -90.0
        write(15,*) '  planet is ice-free.'
      else
        do 765 k=1,nedge,1
           write(15,762) iceline(k)
 762       format(2x,'ice-line latitude = ',f5.1,' degrees.')
 765    continue
        if((nedge.eq.2)) then
          if(zntempave(nbelts).le.icetemp) then
            icelineN = iceline(2)
            icelineS = iceline(1)
            icelineNMax = 90.0
            icelineNMin = iceline(2)
            icelineSMax = iceline(1)
            icelineSMin = -90.0
          else
            icelineN = iceline(2)
            icelineS = iceline(1)
            icelineNMax = iceline(2)
            icelineNMin = 0.0
            icelineSMax = 0.0
            icelineSMin = iceline(1)
          end if
        else
          if(iceline(1) .gt. 0.0) then
             icelineN = iceline(1)
             icelineS = -90.0
             icelineNMax = 90.0
             icelineNMin = iceline(1)
             icelineSMax = -0.0
             icelineSMin = -90.0
          else
             icelineN = 90.0
             icelineS = iceline(1)
             icelineNMax = 90.0
             icelineNMin = 90.0
             icelineSMax = iceline(1)
             icelineSMin = -90.0
          end if
        end if
      end if
      write(19,766) relsolcon, icelineS, icelineN
 766  format(f5.3,2x,f8.3,2x,f8.3)

      if ( fillet ) then
        write(51,768) 0, relsolcon*solarcon/1361.0, obl, 
     &                   fco2*1.e6, 
     &                   ann_tempave, icelineNMax, icelineNMin, 
     &                   icelineSMax, icelineSMin, d0, ann_irave
 768    format(i1,1x,f4.2,1x,f4.1,1x,f8.1,1x,f6.2,1x,f4.1,1x,f4.1,
     &         1x,f6.1,1x,f6.2,1x,f4.2,1x,f6.2)
      end if

c  CO2 CLOUDS
      write(15,770)
 770  format(/ 'CO2 CLOUDS')
      write(15,772) ann_co2cldave
 772  format(2x,'planet average co2 cloud coverage = ',f5.3)

c  GEOGRAPHY INFORMATION 

 779  format (/ 'GEOGRAPHIC DATA')
      write(15,779)
      write(15,780) ocean
 780  format('planet ocean fraction: ',f3.1)
      write(15,782) coastlat
 782  format('coast latitude: ',f5.1)
      write(15,784) igeog
 784  format('geography: ',i1)
      write(15,785) nbelts
 785  format('number of belts: ',i2)
      write(15,786) 
 786  format(2x,'latitude of belt center',2x,'belt ocean fraction')
      do 788 k = 1,nbelts,1
         write(15,787) cntrlatangle(k),focean(k)
 787     format(11x,f5.1,19x,f5.3)
 788  continue

 790  fluxavesum = 0.
      albavesum = 0.
      tempavesum = 0.
      iravesum = 0.
      wthratesum = 0.
      t = 0
      nstep = 0
      if(last) stop

c      write(19,712) pco2,zntempave(1:nbelts)
c 712  format(2x,f12.3,20(2x,f8.3))

c      open(unit=765,file= 'out_600_2.1c3/LATdependence.out')
c      !write(765,9997) zntempmave(1:nbelts)
c      write(765,9997) zntempave(1:nbelts)
c9997  format(20(2x,f8.3))


      !check for convergence
      if ( iterhalt ) then
        !print *, "YEAR = ", yrcnt      

        if( yricnt .ge. niter ) goto 1000
        yricnt = yricnt + 1
        yrcnt  = yrcnt + yrstep

        !print *, "year = ", yrcnt
        print *, "iteration = ", yrcnt / yrstep

      else

        if(abs(prevtempave-ann_tempave).lt.cnvg) goto 1000

      end if

      prevtempave = ann_tempave      

 800  if(seasons) goto 240
      goto 260

c----------------------------------------------------------------------c
c  WRAP UP

 1000 write(*,*) 'The energy-balance calculation is complete.'
      write(15,1010)
 1010 format(/ 'SUMMARY')
      write(15,1015) 'planet average temperature = ', ann_tempave, 
     & ' Kelvin'
      write(15,1015) 'planet average albedo = ', ann_albave
      write(15,1015) 'planet average insolation = ', ann_fluxave, 
     & ' Watts/m^2'
      write(15,1015) 'planet average outgoing infrared = ',ann_irave, 
     & ' Watts/m^2'
      write(15,1016) 'total pressure pg0 = ', pg0, ' bars'
      write(15,1016) 'co2 partial pressure = ', pco2, ' bars'
      write(15,1016) 'n2 partial pressure = ', pn2, ' bars'
      write(15,1016) 'h2 partial pressure = ', pf2, ' bars' 
 1015 format(3x,a,f7.3,a)
      write(15,1016) 'thermal diffusion coefficient (D) = ', d, 
     & ' Watts/m^2 K'
 1016 format(3x,a,e8.3,a)
      write(15,1020) 'convergence to final temperature profile in ',
     & tcalc, ' seconds'
 1020 format(3x,a,e9.3,a)

c  LOOP ONE MORE TIME, AND WRITE OUTPUT TO FILES
      last = .TRUE.

c  initialize zntempmin matrix

      do 1125 k = 1,nbelts,1
         zntempmin(k) = 500. !**some large temperature
 1125 continue
c
      write(15,1130)
 1130 format(/ 'ZONAL STATISTICS')
      write(15,1135)
 1135 format(1x,'latitude(deg)',2x,'ave temp(K)',2x,'min temp(K)',
     &  2x,'@dec(deg)',2x,'max temp(K)',2x,'@dec(deg)',2x,
     &  'planet albe')

      if ( fillet ) then
        write(50,1137)
 1137   format(/ '# Lat Tsurf Asurf ATOA OLR')
        write(51,1138)
 1138   format(/ '# Case Inst Obl XCO2 Tglob IceLineNMax IceLineNMin 
     &IceLineSMax IceLineSMin Diff OLRglob')
      end if

 1100 format(/ 'OUTPUT FILES')
      tcalc = 0.
      if(seasons) goto 240
      goto 260

      if ( radparam .eq. 3 ) then
        call radiation_end
      end if

      end

c----------------------------------------------------------------------c
      subroutine  geog(igeog)

c  This subroutine calculates continental/oceanic zone fractions for five
c  geography types (1) present geography (2) polar supercontinent (3) equa-
c  torial supercontinent (4) 100% oceanic and (5) 100% land (10-27-95)

      parameter (pi=3.141592653589793,nbelts=18)

      common/geogblk/cntrlat(nbelts),cntrlatangle(nbelts),coastlat,
     &   focean(nbelts),lat(0:nbelts+1),latangle(0:nbelts+1),ocean

      real*4  lat,latangle
      logical coast

      goto (10,20,30,40,50,60) igeog
      
c  PRESENT GEOGRAPHY from Sellers (1965)
 10   rewind(2)
      n = 1
      coastlat = 0.
 12   read(2,*,end=15) focean(n)
      focean(n) = focean(n)*ocean/0.7
      n = n + 1
      goto 12
c **make sure there is as many ocean data as nbelts
 15   if (n .lt. nbelts+1) then
         write(*,18) nbelts
         stop
      end if
 18   format('Need more ocean data to run model for ',i2,' lat belts!'
     & 'Energy-balance calculation aborted.')

      return


c  SOUTH POLAR SUPERCONTINENT 
 20   coast = .TRUE.    !** still on land
      coastlat = asin(1-2*ocean)*180/pi
      do 25 k = 1,nbelts,1
         if (latangle(k).le.coastlat) then
            focean(k) = 0.

c **see whether you are on the coast
            if (latangle(k) .eq. coastlat) coast = .FALSE.
  
         else if ((latangle(k).gt.coastlat).and.(coast)) then
            focean(k) = (sin(lat(k)) - sin(coastlat*pi/180.))/
     &         (sin(lat(k)) - sin(lat(k)-pi/nbelts))          
            coast = .FALSE.     !** now in water
         else
            focean(k) = 1.
         end if
 25   continue

      return

c  EQUATORIAL SUPERCONTINENT
 30   coast = .TRUE.     !**still on land
      coastlat = asin(1-ocean)*180/pi
      do 35 k = nbelts/2+1,nbelts,1
         if (latangle(k) .le. coastlat) then
            focean(k) = 0.
            focean(nbelts-k+1) = 0.

c **see whether you are on the coast
            if (latangle(k) .eq. coastlat) coast = .FALSE. 
             
         else if ((latangle(k) .gt. coastlat) .and. (coast)) then
            focean(k) = (sin(lat(k)) - sin(coastlat*pi/180.))/
     &         (sin(lat(k)) - sin(lat(k)-pi/nbelts))          
            focean(nbelts-k+1) = focean(k)
            coast = .FALSE.    !**now in water
         else
            focean(k) = 1.
            focean(nbelts-k+1) = 1.
         end if
 35   continue

      return

c  100% WATER COVERED***
 40   do 45 k = 1,nbelts,1
         focean(k) = 1.
 45   continue
      return

c  100% LAND SURFACE***
 50   do 55 k = 1,nbelts,1
         focean(k) = 0.
 55   continue
      return

c  EQUAL LAND/OCEAN FRACTION AT ALL LATITUDE BANDS
 60   do 65 k = 1,nbelts,1
         focean(k) = ocean
 65   continue
      return

      end

c---------------------------------------------------------------------c
      subroutine  keplereqn(m,e,x)
c---------------------------------------------------------------------c
c  This subroutine is used to find a solution to Keplers' equation,
c  namely M = E - eSin(E), where M is the mean anomaly, E is the
c  eccentric anomaly, and e is the eccentricity of the orbit.
c  Input is M and e, and output is E. The chosen convergence is quartic 
c  as outlined in Danby (1988), pp. 149-54.  coded from Danby (5-11-95)

      implicit  integer(n), real*8(a-m,o-z)
      parameter (pi=3.141592653589793)

c  initial guess is x = M + k*e with k = 0.85
      k = 0.85
c  bound on the permissible error for the answer
      del = 1.e-13

      ms = m - int(m/(2*pi))*(2*pi)
      sigma = sign(1.,sin(ms))
      x = ms + sigma*k*e
      nc = 0
 100  f = x - e*sin(x) - ms
c  check for convergenge
      if (abs(f) .lt. del) goto 200
      fp = 1. - e*cos(x)
      fpp = e*sin(x)
      fppp = e*cos(x)
      dx1 = -f/fp                 ! Newton's Method (quadratic convergence)
      dx2 = -f/(fp + dx1*fpp/2.)  ! Halley's Method (cubic convergence)
      dx3 = -f/(fp + dx2*fpp/2. + dx2**2*fppp/6.) ! (quartic convergence)
      x = x + dx3
c      write(*,*)  nc,x,f,dx3  !**uncomment this to check for convergence
      nc = nc + 1
c  stop calculation if number of iterations is too large
      if (nc .gt. 10) then
         write(*,*) 'no solution to keplers equation'
         stop
      end if
      goto 100
 200  return
      end

c------------------------------------------------------------------
c Subroutines for random number generation
c------------------------------------------------------------------

      SUBROUTINE SRAND(ISEED)

C  This subroutine sets the integer seed to be used with the
C  companion RAND function to the value of ISEED.  A flag is 
C  set to indicate that the sequence of pseudo-random numbers 
C  for the specified seed should start from the beginning.

      COMMON /SEED/JSEED,IFRST

      JSEED = ISEED
      IFRST = 0

      RETURN
      END
      REAL FUNCTION RAND()

C  This function returns a pseudo-random number for each invocation.
C  It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal 
C  standard number generator whose Pascal code appears in the article:
C
C     Park, Steven K. and Miller, Keith W., "Random Number Generators: 
C     Good Ones are Hard to Find", Communications of the ACM, 
C     October, 1988.

      PARAMETER (MPLIER=16807,MODLUS=2147483647,MOBYMP=127773,
     +           MOMDMP=2836)

      COMMON  /SEED/JSEED,IFRST
      INTEGER HVLUE, LVLUE, TESTV, NEXTN
      SAVE    NEXTN

      IF (IFRST .EQ. 0) THEN
        NEXTN = JSEED
        IFRST = 1
      ENDIF

      HVLUE = NEXTN / MOBYMP
      LVLUE = MOD(NEXTN, MOBYMP)
      TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
      IF (TESTV .GT. 0) THEN
        NEXTN = TESTV
      ELSE
        NEXTN = TESTV + MODLUS
      ENDIF
      RAND = REAL(NEXTN)/REAL(MODLUS)

      RETURN
      END
      BLOCKDATA RANDBD
      COMMON /SEED/JSEED,IFRST

      DATA JSEED,IFRST/123456789,0/

      END

c------------------------------------------------------------------
c Polar form of the Box-Muller transformation 
c------------------------------------------------------------------

      real function boxmuller()
        real x1, x2, w, y1, y2

 1800   continue  
          x1 = 2.0 * RAND() - 1.0;
          x2 = 2.0 * RAND() - 1.0;
          w = x1 * x1 + x2 * x2;
        if ( w .ge. 1.0 ) goto 1800

        w = sqrt( (-2.0 * log( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;

        boxmuller = y1
        return
      end

c------------------------------------------------------------------
c Subroutines to check for CO2 condensation
c------------------------------------------------------------------

      SUBROUTINE SATCO2(T,PSCO2)
c This subroutine calculates the vapor saturation pressure for CO2

      IF (T.LT.216.56) GO TO 1

c   VAPOR PRESSURE OVER LIQUID
      PSL = 3.128082 - 867.2124/T + 1.865612E-2*T
     2  - 7.248820E-5*T*T + 9.3E-8*T*T*T
      PATM = 10.**PSL
      PSCO2 = 1.013*PATM

      RETURN

c   VAPOR PRESSURE OVER SOLID
C   (Altered to match vapor pressure over liquid at triple point)
   1  PSL = 6.760956 - 1284.07/(T - 4.718) + 1.256E-4*(T - 143.15)
      PATM = 10.**PSL
      PSCO2 = 1.013*PATM 
      RETURN
      END
