      program  energy_balance_climate_model
c--------------------------------------------------------------------c
c  This program calculates seasonal latitudinal temperature 
c  gradients for Earth-like planets with N2,O2,CO2, and H20 
c  atmospheres (e.g. Earth and Mars). Energy balance is treated as in 
c  Caldeira and Kasting (1992). Input files: 'oceans.dat', for
c  modeling Earth with present geography, and 'fresnel_reflct.dat'.
c
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
    
      real*4  lat,latangle,mu,ir,msun,landalb,icealb,irsum,irave,mp,
     &  iravesum,iceline,icelat,icesht
      real*8  ecc,m,e

      parameter (nbelts=18)
      parameter (pi=3.14159265,grav=6.6732e-8,msun=1.9891e33)
      parameter (gacc=9.81)
      parameter (mp=1.67e-24,q0=1360.,cnvg=1.e-2)
      parameter (sbc=5.67e-8,emis=0.64)
      parameter (twopi=2*pi)
      parameter (niter=10)
c     parameter (niter=1000)
c     parameter (niter=3000)
c     parameter (niter=5000)
c     parameter (niter=50000)
c     parameter (niter=2001)
      parameter (niterhalf=1001)
      parameter (niterquarter=1501)
   
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
     &  co2ice(nbelts), co2melt(nbelts), zthick(nbelts),
     &  znco2icemax(nbelts)

      character  header*80,file(0:3)*8
      logical seasons, last, snowball, linrad, linalb, cloudalb
      logical stochastic, soladj, constheatcap, diffadj, iterhalt
      logical do_cs_cycle, do_h2_cycle, doco2cond, doco2ice 
      logical doco2albedo, dobandalbedo, doseaweather
      real landsnowfrac, RAND, boxmuller, noisevar
      real seaweather, seawthrate, betaexpsea
      real outgassing, weathering, betaexp, kact, krun, fn2
      real pg0, ir2, fh2, co2sat, h2escape, ph2, ncolh2, h2outgas
      real co2albedo, visalbedo, iralbedo, fvisible, fnearir
      real kco2ice, rhoco2ice, zmax, geoflux, newfreeze, pco20s
      integer(kind=8) yrcnt, yrstep, radparam, co2flag
      integer ISEED, resfile
      integer*4 now(3)

      dimension solcon(niter),prec(niter),ecce(niter),
     &  yrlabel(niter),obliq(niter)
       
      data  temp/20*273./       !**use some temperature > 265K, otherwise iceball
      data  file/ 'spng.out', 'summ.out', 'fall.out', 'wint.out' /

c  NAMELIST PARAMETERS
      INTEGER namelistid
      OPEN( namelistid, FILE='./input.nml', DELIM='APOSTROPHE' )

c  INITIALIZE VARIABLES
      ann_tempave = 10.         !just some constant greater than zero
      seasons = .true.     !if .true. do seasons: otherwise do mean-annual
      last = .false.       !should always be .false.
      snowball = .false.   !if .true. then start as an ice-ball
      stochastic = .false. !if .true. then include stochastic perturbation
      soladj = .false.     !if .true. then include solar forcing adjustment
      resfile = 0          !start from previous (1), present Earth (2) hothouse (3), large ice (4)
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
      pco20 = 3.3e-4   !present atmospheric CO2 concentration (bars)
      pco2 = pco20              !initial co2 concentration for this run (bars)
      ocean = 0.7      !planet fraction covered by water
      igeog = 1        !geography ('1'=present,'2'=polar,'3'=equatorial,
     &                 ! '4'=100% water,'5'=100% land)
      groundalb = 0.20  !base land surface albedo
      !groundalb = 0.25
      relsolcon = 1.0  !relative solar constant (1.0 for present Earth) 
      alpha = -0.078            !alpha: cloud albedo = alpha + beta*zrad
      beta = 0.65      !beta
      cw = 2.1e8       !heat capacity over ocean
      cl = 5.25e6      !heat capacity over land
      d0 = 0.58        !thermal diffusion coefficient
      v = 3.3e14  !volcanic outgassing (g/year) [from Holland (1978)]
      wco2 = 9.49e14  !carbonate-silicate weathering constant
      avemol0 = 28.89*mp  !mean-molecular weight of present atmosphere (g)
      hcp0 = 0.2401    !heat capacity of present atmosphere (cal/g K)
      hcpco2 = 0.2105  !heat capacity of co2 (cal/g K)
      hcpn2 = .2484     ! heat capacity of n2
      hcph2 = 3.420      !heat capacity of h2 
      kco2ice = 0.6      !thermal conductivity of CO2 ice (W/m/K)
      rhoco2ice = 1600.  !density of CO2 ice (kg/m^3)
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
      doco2cond = .true. ! set .false. to disable CO2 condensation
      doco2ice = .false. ! set .true. to enable CO2 ice mass balance
      doco2albedo = .false. ! set .true. to enable CO2 ice albedo
      co2albedo = 0.30  ! albedo of CO2 ice
      dobandalbedo = .false. !set .true. to enable 2-band ice albedo
      visalbedo = 0.90 ! visible albedo for 2-band ice albedo
      iralbedo = 0.70  ! near-IR albedo for 2-band ice albedo
      fvisible = 0.60  ! fraction of visible radiation for 2-band ice albedo
      fnearir  = 0.40  ! fraction of near-IR radiation for 2-band ice albedo
      geoflux = 0.0    ! geothermal heat flux (W/m^2)
      outgassing = 7.0  ! volcanic outgassing rate (bars/Gyr)
      weathering = 7.0  ! weathering rate (bars/Gyr)
      betaexp = 0.50       ! weathering exponent
      kact = 0.09       ! activation energy factor
      krun = 0.045      ! runoff efficiency factor
      pg0 = 1.0         ! surface pressure (bars)
      fh2 = 0.0
      h2outgas = 2.67e12
      radparam = 0      ! OLR/albedo parameterization: (0) Williams & Kasting, (1) CO2/N2, (2) CO2/H2
                        ! (3) F-star (4) G-star (5) K-star (6) M-star
      doseaweather = .false.  ! seafloor weathering
      seaweather   = 7.0      ! seafloor weathering rate (bars/Gyr)
      betaexpsea   = 0.0      ! seafloor weathering exponent
      pco20s = 0.01   !present soil CO2 concentration (bars)
      !pco20s = 3.3e-4   !present soil CO2 concentration (bars)
      nfile = 0
      yrstep = 1
      yrcnt = 0
      yricnt = 1

      CALL itime( now )
      CALL SRAND( now(3) )


      NAMELIST /ebm/ seasons, snowball, tend, dt, rot, a, ecc, peri, 
     &               obl, cloudir, ocean, igeog, groundalb, 
     &               relsolcon, landsnowfrac, fcloud, yrstep,
     &               resfile, noisevar, stochastic, soladj, d0,
     &               linrad, linalb, constheatcap, diffadj,
     &               cloudalb, iterhalt, do_cs_cycle, outgassing, 
     &               weathering, betaexp, kact, krun, pco2, fh2,
     &               radparam, do_h2_cycle, h2outgas, doco2cond,
     &               doco2ice, doco2albedo, co2albedo, dobandalbedo,
     &               visalbedo, iralbedo, geoflux, doseaweather, 
     &               seaweather, betaexpsea, pco20s
      READ( namelistid, NML=ebm )
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
      open (unit=99,file='data/insolaout.dat',status='old')

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
        pn2 = 1.0
        pg0 = pco2 + pn2
        ph2 = 0.0
        fh2 = 0.0
      else if ( radparam .eq. 2 ) then
        fn2 = 0.05
        pg0 = pco2 / ( 1 - fn2 - fh2 )
        pn2 = fn2 * pg0
        !pn2 = 1.0
        !pg0 = pn2+pco2
        !fn2 =  pn2/pg0
         ph2 = fh2*pg0
      else if ( radparam .ge. 3 ) then
        pn2 = 1.0
        pg0 = pco2 + pn2
        ph2 = 0.0
        fh2 = 0.0
      end if


      if( snowball .and. ( resfile .ne. 0 ) ) then
        print *, 'Namelist parameters SNOWBALL and RESFILE cannot
     & both be used.' 
        print *, 'EBM Calculation aborted.'
        stop
      end if
  
      if( snowball ) then
        do k = 0, nbelts+1, 1
          temp(k) = 233.
        end do
      else if( resfile .eq. 0 ) then
        do k = 1, nbelts, 1
          temp(k) = 273.15
        end do   
      else if( resfile .eq. 1 ) then	!restart from last run
        rewind(6)
        do k = 1, nbelts, 1
           read(6,*) temp(k)
        end do
      else if( resfile .eq. 2 ) then	!restart from present Earth
        rewind(4)
        do k = 1, nbelts, 1
           read(4,*) temp(k)
        end do
      else if( resfile .eq. 3 ) then	!restart from hothouse
        rewind(5)
        do k = 1, nbelts, 1
           read(5,*) temp(k)
        end do
      else if( resfile .eq. 4 ) then	!restart from big ice cap
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
      do k = 1, nbelts, 1

        diff(k) = d0

      end do

c----------------------------------------------------------------------c
c  SET SOLAR CONSTANT, ECCENTRICITY, AND OBLIQUITY

      if( soladj ) then   
        rewind(99)
        do p = 1, niter, 1
           read(99,*) yrlabel(p),ecce(p),prec(p),obliq(p),solcon(p)
           solcon(p) = 4*solcon(p) * relsolcon
           prec(p)   = ASIN(prec(p)/ecce(p))*180/pi
        end do
      else
        do p = 1, niter, 1
          solcon(p) = q0 * relsolcon
          ecce(p)   = ecc
          obliq(p)  = obl
          prec(p)   = peri
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
 125  continue

c  WRITE GRID DATA TO OUTPUT FILE

c      write(15,140) 
 140  format(/ / '**LATITUDE GRID DATA**')
c      write(15,142)
 142  format(10x,'x',6x,'dx',4x,'latitude(rad)',2x,'latitude(deg)')
      do 145 k = 0,nbelts+1,1
c         write(15,143) k,x(k),dx(k),lat(k),latangle(k)
 143     format(3x,i2,2x,f6.3,2x,f6.3,5x,f6.3,8x,f5.1)
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
 148     format(17x,2p10(f4.1,1x))
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
      else
         write(*,*) 'Calculation time has elapsed.'
         goto 1000
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
      else
         write(*,*) 'Calculation time has elapsed.'
         goto 1000
      end if

 270  tempsum = 0.
      tempmaxsum = 0.
      tempminsum = 0.
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
   
        if ( doco2ice ) then
          phi = log10( pco2 - co2ice(k) )
        else
          phi = log10( pco2 )
        end if
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

        if ( doco2ice ) then
          phi = log10( pg0 - co2ice(k) )
        else
          phi = log10( pg0 )
        end if
        tmpk = log10( temp(k) )   

        term1=4.76808656814761633314e+01*tmpk**4+
     &  5.49935481767358158578e+00*tmpk**3*fh2-
     &  1.09429061839906305309e+00*tmpk**3*phi-
     &  4.78173291533181611612e+02*tmpk**3+
     &  8.49190549496696478471e-01*tmpk**2*fh2**2+
     &  1.37765937761845114196e+00*tmpk**2*fh2*phi-
     &  3.66017123367183145888e+01*tmpk**2*fh2+
     &  2.39017525956834697709e+00*tmpk**2*phi**2+
     &  1.94856225493008743399e+01*tmpk**2*phi

        term2=+1.79716986886717290872e+03*tmpk**2-
     &  5.08081407778601223946e-01*tmpk*fh2**3+
     &  3.00410810274085282590e-01*tmpk*fh2**2*phi-
     &  3.10430848409695459011e+00*tmpk*fh2**2-
     &  4.70658106780907359301e-02*tmpk*fh2*phi**2-
     &  6.75476435103275996141e+00*tmpk*fh2*phi+
     &  8.02983658642542934558e+01*tmpk*fh2-
     &  2.78693163775443009111e-01*tmpk*phi**3-
     &  1.28002847689644276841e+01*tmpk*phi**2-
     &  7.61706327338787190229e+01*tmpk*phi

        term3=-2.99763768092468944815e+03*tmpk-
     &  4.64112012085224553970e-02*fh2**4+
     &  1.80172294709195594808e-02*fh2**3*phi+
     &  1.31650299007109294891e+00*fh2**3-
     &  1.53163365643916478398e-02*fh2**2*phi**2-
     &  8.55693979833223372644e-01*fh2**2*phi+
     &  2.28114230486580993329e+00*fh2**2+
     &  1.79311735513174204393e-02*fh2*phi**3+
     &  1.34635973577179490768e-01*fh2*phi**2+
     &  8.13967439897796829484e+00*fh2*phi

        term4=-5.80918681555363392022e+01*fh2+
     &  1.83043738623349175332e-02*phi**4+
     &  7.45155230077319874482e-01*phi**3+
     &  1.69576563201193479813e+01*phi**2+
     &  8.54071288991185468831e+01*phi+
     &  1.87594411270478985898e+03

      ir(k) = 10**(term1 + term2 + term3 + term4) / 1000.

      else if ( radparam .eq. 3 ) then   

        fvisible = 0.67
        fnearir  = 0.33
        if ( doco2ice ) then
          phi = log10( pco2 - co2ice(k) )
        else
          phi = log10( pco2 )
        end if
        tmpk = log10( temp(k) )   
        INCLUDE './data/write_Fortran_OLR7200K.dat'
        ir(k) = 10**(term1 + term2) / 1000.

      else if ( radparam .eq. 4 ) then   

        fvisible = 0.52
        fnearir  = 0.48
        if ( doco2ice ) then
          phi = log10( pco2 - co2ice(k) )
        else
          phi = log10( pco2 )
        end if
        tmpk = log10( temp(k) )  
        INCLUDE './data/write_Fortran_OLR5800K.dat'
        ir(k) = 10**(term1 + term2) / 1000.

      else if ( radparam .eq. 5 ) then   

        fvisible = 0.32
        fnearir  = 0.68
        if ( doco2ice ) then
          phi = log10( pco2 - co2ice(k) )
        else
          phi = log10( pco2 )
        end if
        tmpk = log10( temp(k) )   
        INCLUDE './data/write_Fortran_OLR4600K.dat'
        ir(k) = 10**(term1 + term2) / 1000.

      else if ( radparam .eq. 6 ) then   

        fvisible = 0.10
        fnearir  = 0.90
        if ( doco2ice ) then
          phi = log10( pco2 - co2ice(k) )
        else
          phi = log10( pco2 )
        end if
        tmpk = log10( temp(k) )   
        INCLUDE './data/write_Fortran_OLR3400K.dat'
        ir(k) = 10**(term1 + term2) / 1000.

      end if

        ir(k) = ir(k) - cloudir   !**reduction of outgoing-infrared by clouds

      end if

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
  
      z = acos(mu(k))*180./pi
      zrad = acos(mu(k))
                 
c----------------------------------------------------------------------c
c  HEAT CAPACITY and SURFACE ALBEDO                    

      oceanalb = h20alb(int(z))

      ! Added this check so fractional ice cover would only 
      ! occur for 263 < T < 273.  (JDH)
      if (temp(k).ge.273.15) then
         fice(k) = 0
      else if (temp(k).lt.263.15) then
         fice(k) = 1
      else
         fice(k) = 1. - exp((temp(k)-273.15)/10.)
      end if

      if ( cloudalb ) then
        acloud(k) = alpha + beta*zrad
      else
        acloud(k) = 0.0 
      end if

      if (temp(k).le.273.15) goto 415

      landalb = groundalb
      fice(k) = 0.
      fwthr = 1.       !**100% if T > 273K, 0% otherwise
      c(k) = focean(k)*cw + (1 - focean(k))*cl
      goto 430

c CO2 ICE ALBEDO IF CO2 IS CONDENSING AT SURFACE
 415  if ( doco2albedo ) then
        CALL SATCO2( temp(k), co2sat )
        if ( pco2 .gt. co2sat ) then
          snowalb = co2albedo
          goto 421
        end if
      end if

      if ( dobandalbedo ) then
        snowalb = fvisible*visalbedo + fnearir*iralbedo
        goto 421
      end if

c LAND WITH STABLE SNOW COVER; SEA-ICE WITH SNOW COVER
 420  snowalb = 0.663  ! changed this from 0.7 as per Caldiera & Kasting (JDH)
 421  landalb = snowalb*landsnowfrac + groundalb*(1 - landsnowfrac)
      icealb = snowalb

      ci = 1.05e7
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
c        wthrate(k) = warea(k)*weathering*((pco2/pco20)**betaexp)
c     &   *exp(kact*(temp(k)-288.))*(1+(krun*(temp(k)-288.)))**0.65
        wthrate(k) = warea(k)*weathering*((pco2/pco20s)**betaexp)
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

        if ( doco2ice ) then
          phi = log10( pco2 - co2ice(k) )
        else
          phi = log10( pco2 )
        end if
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

        if ( doco2ice ) then
          phi = log10( pg0 - co2ice(k) )
        else
          phi = log10( pg0 )
        end if
        tmpk = log10( temp(k) )

        ! updated albedo fits for CO2/H2 atmospheres (2015)
        if ( temp(k) .le. 250 ) then

        term1 =-3.41881382273888756451e-01*mu(k)**3-
     &  7.82426187892112934286e-01*mu(k)**2*as-
     &  2.98234592912095253237e-01*mu(k)**2*tmpk-
     &  1.04812948671497343373e-02*mu(k)**2*fh2+
     &  8.78831353460830799751e-03*mu(k)**2*phi+
     &  1.99493626489599984453e+00*mu(k)**2+
     &  4.00889848532216958032e-02*mu(k)*as**2-
     &  2.63781358942446579607e-02*mu(k)*as*tmpk-
     &  1.13738660707575038600e-02*mu(k)*as*fh2

        term2=-5.92075899402270330873e-02*mu(k)*as*phi+
     &  1.16592410386607769901e+00*mu(k)*as+
     &  3.05817932525109226205e+00*mu(k)*tmpk**2+
     &  4.55404679827312378060e-02*mu(k)*tmpk*fh2-
     &  4.09807970694542444967e-01*mu(k)*tmpk*phi-
     &  1.44076413354266730238e+01*mu(k)*tmpk-
     &  5.90864025206625214509e-03*mu(k)*fh2**2-
     &  1.37541407363503485201e-03*mu(k)*fh2*phi-
     &  9.09945469993787836582e-02*mu(k)*fh2+
     &  1.64338897369647428393e-02*mu(k)*phi**2

        term3=+1.00187017644873965772e+00*mu(k)*phi+
     &  1.56017563623685315832e+01*mu(k)+
     &  4.37284669736712483523e-02*as**3+
     &  1.43728397075257047916e-02*as**2*tmpk+
     &  1.08542075534173851001e-03*as**2*fh2+
     &  5.62908389505368811356e-02*as**2*phi-
     &  5.96205015905240578306e-03*as**2-
     &  3.37641094945521613724e+00*as*tmpk**2-
     &  8.30553802255556555822e-02*as*tmpk*fh2+
     &  3.93819274461570598134e-03*as*tmpk*phi

        term4=+1.52729776360046312078e+01*as*tmpk-
     &  6.64563123009531860064e-03*as*fh2**2+
     &  4.09531648961306438822e-03*as*fh2*phi+
     &  2.06557952087644880468e-01*as*fh2-
     &  5.19553098335195778779e-02*as*phi**2-
     &  2.68917955568870714611e-01*as*phi-
     &  1.70134666260915210501e+01*as-
     &  2.05383049531839283475e+01*tmpk**3-
     &  8.51825973656873070006e-01*tmpk**2*fh2+
     &  1.17558766493159172484e+00*tmpk**2*phi

        term5=+1.41908353208508884791e+02*tmpk**2-
     &  6.20670129333104891867e-02*tmpk*fh2**2+
     &  6.82602368680743920581e-02*tmpk*fh2*phi+
     &  4.06748498001836988891e+00*tmpk*fh2-
     &  5.91040203477738376736e-02*tmpk*phi**2-
     &  5.35717146939420985063e+00*tmpk*phi-
     &  3.26479492305605447200e+02*tmpk+
     &  7.32643081307610393588e-03*fh2**3-
     &  1.79603740710287279082e-03*fh2**2*phi+
     &  1.41628420988978820372e-01*fh2**2

        term6=-2.36279193983561118172e-03*fh2*phi**2-
     &  1.65259997145445591826e-01*fh2*phi-
     &  4.85728489815857855660e+00*fh2+
     &  4.95864284669994548338e-03*phi**3+
     &  1.82075298055801054753e-01*phi**2+
     &  6.24336067599479171975e+00*phi+
     &  2.50680778713592985696e+02

        else

        term1 =-3.64953635588242231158e-01*mu(k)**3-
     &  3.36184863117452492620e-01*mu(k)**2*as+
     &  7.42855566182375470774e-01*mu(k)**2*tmpk+
     &  2.35250045071474329569e-03*mu(k)**2*fh2+
     &  1.03140984013693603333e-02*mu(k)**2*phi-
     &  7.50158525635715056623e-01*mu(k)**2+
     &  6.19678741078302092182e-02*mu(k)*as**2-
     &  1.88843608069845148023e+00*mu(k)*as*tmpk-
     &  3.51228279746902474767e-02*mu(k)*as*fh2

        term2=-6.59746479346320358061e-02*mu(k)*as*phi+
     &  5.16167787639029462810e+00*mu(k)*as-
     &  2.12155537739597477298e+00*mu(k)*tmpk**2-
     &  7.85503169940939272031e-02*mu(k)*tmpk*fh2+
     &  1.67083962229760435436e-01*mu(k)*tmpk*phi+
     &  1.04507082880023816074e+01*mu(k)*tmpk-
     &  6.71874522052537132694e-03*mu(k)*fh2**2+
     &  4.54865891120825524552e-03*mu(k)*fh2*phi+
     &  2.12861526648225163338e-01*mu(k)*fh2+
     &  1.72902842337580193999e-02*mu(k)*phi**2

        term3=-3.91468482134287643071e-01*mu(k)*phi-
     &  1.39550191937482175319e+01*mu(k)+
     &  7.66012173959446096561e-02*as**3-
     &  5.33258519884163770253e-02*as**2*tmpk-
     &  3.24457891729390632968e-03*as**2*fh2+
     &  2.59435678221003190869e-02*as**2*phi+
     &  5.62019786112206570783e-02*as**2-
     &  5.97292118888921041986e+00*as*tmpk**2-
     &  3.87009770489884163958e-01*as*tmpk*fh2+
     &  7.19742438017603736178e-01*as*tmpk*phi

        term4=+2.84282223773413313950e+01*as*tmpk-
     &  3.66795093171868660797e-02*as*fh2**2+
     &  1.73632089463807148810e-02*as*fh2*phi+
     &  9.70804413891795059399e-01*as*fh2-
     &  3.25015124546070566236e-02*as*phi**2-
     &  1.93009348053676799140e+00*as*phi-
     &  3.35460424111451231965e+01*as+
     &  4.25837679307632441805e+01*tmpk**3+
     &  2.33414006785506211727e+00*tmpk**2*fh2-
     &  4.63279637872022642675e-01*tmpk**2*phi

        term5=-3.12848837120424946079e+02*tmpk**2-
     &  1.10661505683727779542e-02*tmpk*fh2**2+
     &  6.43891660759138700909e-02*tmpk*fh2*phi-
     &  1.13498362686029636848e+01*tmpk*fh2+
     &  1.50766229910903215572e-01*tmpk*phi**2+
     &  2.22835151136342046740e+00*tmpk*phi+
     &  7.65279777491220670527e+02*tmpk-
     &  2.97440540952105363093e-02*fh2**3+
     &  5.87836752705261922358e-03*fh2**2*phi+
     &  6.57729484225070132331e-02*fh2**2

        term6=+3.35025642146789467968e-03*fh2*phi**2-
     &  1.58732057910619589469e-01*fh2*phi+
     &  1.37637199533833634035e+01*fh2+
     &  3.67386094751287863372e-03*phi**3-
     &  3.39029543003353861508e-01*phi**2-
     &  2.53878600425094846926e+00*phi-
     &  6.22751304763008533882e+02

        end if

        atoa(k) = term1 + term2 + term3 + term4 + term5 + term6

      else if ( radparam .eq. 3 ) then

        if ( doco2ice ) then
          phi = log10( pco2 - co2ice(k) )
        else
          phi = log10( pco2 )
        end if
        tmpk = log10( temp(k) )
        if ( temp(k) .le. 250 ) then
          INCLUDE './data/write_Fortran_alb7200K_250.dat'
        else
          INCLUDE './data/write_Fortran_alb7200K_350.dat'
        end if
        atoa(k) = term1 + term2 + term3 + term4

      else if ( radparam .eq. 4 ) then

        if ( doco2ice ) then
          phi = log10( pco2 - co2ice(k) )
        else
          phi = log10( pco2 )
        end if
        tmpk = log10( temp(k) )

        if ( temp(k) .le. 250 ) then
          INCLUDE './data/write_Fortran_alb5800K_250.dat'
        else
          INCLUDE './data/write_Fortran_alb5800K_350.dat'
        end if
        atoa(k) = term1 + term2 + term3 + term4

      else if ( radparam .eq. 5 ) then

        if ( doco2ice ) then
          phi = log10( pco2 - co2ice(k) )
        else
          phi = log10( pco2 )
        end if
        tmpk = log10( temp(k) )
        if ( temp(k) .le. 250 ) then
          INCLUDE './data/write_Fortran_alb4600K_250.dat'
        else
          INCLUDE './data/write_Fortran_alb4600K_350.dat'
        end if
        atoa(k) = term1 + term2 + term3 + term4

      else if ( radparam .eq. 6 ) then

        if ( doco2ice ) then
          phi = log10( pco2 - co2ice(k) )
        else
          phi = log10( pco2 )
        end if
        tmpk = log10( temp(k) )
        if ( temp(k) .le. 250 ) then
          INCLUDE './data/write_Fortran_alb3400K_250.dat'
        else
          INCLUDE './data/write_Fortran_alb3400K_350.dat'
        end if
        atoa(k) = term1 + term2 + term3 + term4

      end if

      end if

c----------------------------------------------------------------------c
c  DIURNALLY-AVERAGED INSOLATION 
      
      if (seasons) then
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
         c(k) = 5.25e6 * 40
      end if

! Ravi's stuff
!      atoa(k) = 0.45

      if( stochastic ) then
        temp(k) = (diff(k)*t2prime(k)-ir(k)+s(k)*(1-atoa(k)))*dt/c(k) 
     &            + sqrt(noisevar)*boxmuller()
     &            + temp(k)
      else
        temp(k) = (diff(k)*t2prime(k)-ir(k)+s(k)*(1-atoa(k)))*dt/c(k) 
     &            + temp(k)
      end if

      if ( temp(k) .le. 175.0 ) then
        temp(k) = 175.0
      end if


      imco2(k) = 0.
      co2flag  = 0

c-nb----------------------------------------------------
c-nb   CO2 Saturation Temperature 
c-nb   Set minimum temperature based on saturation temperature of CO2 

      if ( doco2cond ) then
 
      CALL SATCO2( temp(k), co2sat )
      
      if ( pco2 .gt. co2sat ) then
  
        imco2(k) = 1.
        co2flag  = 1
        temp(k) = 3.1823*log10(pco2)**3 + 10.5165*log10(pco2)**2 + 
     &           28.5760*log10(pco2) + 192.2084 

        newfreeze = pco2 - co2sat - co2ice(k)
        co2ice(k) = pco2 - co2sat
        zthick(k) = zthick(k) + newfreeze*1.e5 / ( gacc * rhoco2ice )
        zmax      = kco2ice * ( 273.16 - temp(k) ) / geoflux
        if ( zthick(k) .gt. zmax ) then
          print *, "CO2 ice is too thick!"
          stop
        end if
        if ( newfreeze .lt. 0.0 ) then
          co2melt(k) = area(k) * newfreeze * gacc * rhoco2ice / 1.e5
        else
          co2melt(k) = 0.0
        end if

      else 

        if ( co2ice(k) .ne. 0.0 ) then         
          co2melt(k) = area(k) * zthick(k) * gacc * rhoco2ice / 1.e5
        else
          co2melt(k) = 0.0
        end if
        co2ice(k) = 0.0
        zthick(k) = 0.0

      end if

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

c      if(.not.last) goto 310
c  ZONAL STATISTICS - if last loop 
      zntempmin(k) = amin1(zntempmin(k),temp(k))
      if (zntempmin(k).eq.temp(k)) zndecmin(k) = decangle
      zntempmax(k) = amax1(zntempmax(k),temp(k))
      if (zntempmax(k).eq.temp(k)) zndecmax(k) = decangle
      zntempsum(k) = zntempsum(k) + temp(k)

      tempmaxsum = tempmaxsum + area(k)*zntempmax(k)
      tempminsum = tempminsum + area(k)*zntempmin(k)

      znco2icemax(k) = amax1(znco2icemax(k),zthick(k))

 310  continue                                     !**end of belt loop

c  **set pole temps equal to adjacent belt temps
      temp(0) = temp(1)
      temp(nbelts+1) = temp(nbelts)

      ! daily output could get large--comment this block if needed for long integrations
      if( .not. last) then
        write(20,609) t/dt+366*yrcnt,tempsum,pg0,pco2,co2flag,fh2,d
 609    format(2x,f12.2,2x,f8.3,2x,e12.3,2x,e12.6,2x,i12,
     &    2x,f8.5,2x,f8.5)
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

        yrwrite = yrcnt / 1.e6


        write(18,711) yrwrite,tempmaxsum,pg0,pco2,co2flag,fh2,d
 711    format(2x,f12.6,2x,f8.3,2x,e12.3,2x,e12.3,2x,i12,
     &    2x,f8.5,2x,f8.5)
c        write(18,711) yrcnt,ann_tempave,pg0,pco2,co2flag,fh2,d
c 711    format(2x,i12,2x,f8.3,2x,e12.3,2x,e12.3,2x,i12,
c     &    2x,f8.5,2x,f8.5)
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

        if ( doco2ice ) then
          pco2 = pco2 + sum( co2melt )
        end if

        if ( doseaweather ) then
          seawthrate = seaweather*((pco2/pco20s)**betaexpsea)
          pco2 = pco2 + (yrstep/1.e9)*(outgassing-sum(wthrate)
     &      - seawthrate)
c          seawthrate = seaweather*((pco2/pco20)**betaexpsea)
c          pco2 = pco2 + (yrstep/1.e9)*(outgassing-sum(wthrate)
c     &      - seawthrate)
        else
          pco2 = pco2 + (yrstep/1.e9)*(outgassing-sum(wthrate))
        end if

        if( stochastic ) then
          pco2 = pco2 + sqrt(noisevar)*boxmuller()*1.e-2
        end if        

        if ( pco2 .lt. 0 ) then
          !pco2 = 3.e-4
          pco2 = 1.e-5
          !pco2 = 1.e-12
        end if

        ! adjust pg0
        if ( radparam .eq. 1 ) then

          pg0 = pco2 + pn2

        else if ( radparam .eq. 2 ) then

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

          pg0 = ( pco2 + pn2 ) / ( 1 - fh2 )

        else if ( radparam .ge. 3 ) then

          pg0 = pco2 + pn2

        end if

        !print *, pg0, pco2, sum(wthrate)

      end if


c  ADJUST DIFFUSION COEFFICIENT
      avemol = mp*(28.0*pn2+44.*pco2 + 2.0*ph2)/(pg0) 
      hcp = (hcpn2*pn2 + pco2*hcpco2 + hcph2*ph2)/(pg0)

      if ( diffadj ) then
      d = d0*(pg0)*((avemol0/avemol)**2)*(hcp/hcp0)*  
     &   (rot0/rot)**2
 
c-nb  
      do k = 1, nbelts, 1
          diff(k) = d
      end do 

      end if
          
c      if(.not.last) goto 790  
c ZONAL SEASONAL AVERAGES 
      shtempave = 0
      nhtempave = 0
      do 732 k = 1, nbelts, 1
         zntempave(k) = zntempsum(k) / nstep
         if ( k .le. nbelts/2 ) then
           shtempave = shtempave + zntempave(k)
         else
           nhtempave = nhtempave + zntempave(k)
         end if     

 732  continue
      shtempave    = shtempave / (nbelts/2)
      nhtempave    = nhtempave / (nbelts/2)
      !print *, "SH/NH temperature difference = ", shtempave - nhtempave
      zntempave(nbelts+1) = zntempave(nbelts)  !**for ice-line calculation

c  FIND ICE-LINES (ANNUAL-AVERAGE TEMP < 263.15K)
      nedge = 0
      do 740 k=1,nbelts,1
      if ((zntempave(k+1)-263.15)*(zntempave(k)-263.15) .lt. 0.) then
         icelat = latangle(k) + ((latangle(k+1)-latangle(k))/
     &   (zntempave(k+1)-zntempave(k)))*(263.15-zntempave(k))
         nedge = nedge + 1
         iceline(nedge) = icelat
      end if
 740  continue

      if(.not.last) goto 790  
      do 750 k=1,nbelts,1
         write(15,751) latangle(k),zntempave(k),zntempmin(k),
     &      zndecmin(k),zntempmax(k),zndecmax(k),
     &      atoa(k)
c-nb     &      (zntempmax(k)-zntempmin(k))/2.
 751     format(4x,f4.0,9x,f8.3,5x,f8.3,5x,f6.2,5x,f8.3,5x,f6.2,5x,f8.3)
         write(16,752) latangle(k),(zntempmax(k)-zntempmin(k))/2.,
     &     10.8*abs(x(k)),15.5*abs(x(k)),6.2*abs(x(k))
 752     format(4x,f4.0,4x,4(3x,f8.3))
         write(6,753) zntempave(k)
 753     format(2x,f8.3)
 750  continue
      
 755  format(/ 'SURFACE DATA')
 756  format(2x,'latitude(deg)',2x,'temp(k)',2x,'belt area',
     &  2x,'weathering area',2x,'zonal weathering rate (g/yr)')
       write(15,760)
 760  format(/ 'ICE LINES (Tave = 263K)')
      if((nedge.eq.0).and.(zntempave(nbelts/2).le.263.)) then
      	write(15,*) '  planet is an ice-ball.' 
        write(19,766) relsolcon, 0.0, 0.0
      else if((nedge.eq.0).and.(zntempave(nbelts/2).gt.263.)) then
     	write(15,*) '  planet is ice-free.'
        write(19,766) relsolcon, -90.0, 90.0
      else
        do 765 k=1,nedge,1
           write(15,762) iceline(k)
 762       format(2x,'ice-line latitude = ',f5.1,' degrees.')
 765    continue
        if((nedge.eq.2)) then
          write(19,766) relsolcon, iceline(1), iceline(2)
        else
          if(iceline(1) .gt. 0.0) then
            write(19,766) relsolcon, -90.0, iceline(1)
          else
            write(19,766) relsolcon, iceline(1), 90.0
          end if
        end if
 766    format(f5.3,2x,f8.3,2x,f8.3)
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

c print out latitude dependent temperatures 
      open(unit=765,file= 'out/LATdependence.out')
      !write(765,9997) zntempmave(1:nbelts)
      write(765,9997) zntempmax(1:nbelts)
9997  format(20(2x,f8.3))

c print out co2 ice thickness
      open(unit=766,file= 'out/CO2ice.out')
      write(766,9998) znco2icemax
9998  format(20(2x,f8.3))

      znco2icemax(:) = 0
      zntempmin(:) = 0
      zntempmax(:) = 0
      zntempsum(:) = 0

      !check for convergence
      if ( iterhalt ) then

        if( yricnt .ge. niter ) goto 1000
        yricnt = yricnt + 1
        yrcnt  = yrcnt + yrstep

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

 1100 format(/ 'OUTPUT FILES')
      tcalc = 0.
      if(seasons) goto 240
      goto 260

      end

c----------------------------------------------------------------------c
      subroutine  geog(igeog)

c  This subroutine calculates continental/oceanic zone fractions for five
c  geography types (1) present geography (2) polar supercontinent (3) equa-
c  torial supercontinent (4) 100% oceanic (5) 100% land and 
c  (6) equal land/ocean fraction at all latitude bands

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
 18   format('Need more ocean data to run model for ',i2,' lat belts!' /
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

C   VAPOR PRESSURE OVER LIQUID
      PSL = 3.128082 - 867.2124/T + 1.865612E-2*T
     2  - 7.248820E-5*T*T + 9.3E-8*T*T*T
      PATM = 10.**PSL
      PSCO2 = 1.013*PATM

      RETURN

C   VAPOR PRESSURE OVER SOLID
C   (Altered to match vapor pressure over liquid at triple point)
   1  PSL = 6.760956 - 1284.07/(T - 4.718) + 1.256E-4*(T - 143.15)
      PATM = 10.**PSL
      PSCO2 = 1.013*PATM 
      RETURN
      END
