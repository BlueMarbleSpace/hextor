      program  energy_balance_climate_model
c--------------------------------------------------------------------c
c  This program calculates seasonal latitudinal temperature 
c  gradients for Earth-like planets with N2,O2,CO2, and H20 
c  atmospheres (e.g. Earth and Mars). Energy balance is treated as in 
c  Caldeira and Kasting (1992). Input files: 'oceans.dat', for
c  modeling Earth with present geography, and 'fresnel_reflct.dat'.
c
c  Original version: 11-3-95 DMW 
c
c  Changes made and commented: 06-30-04 JDH
c
c  Added namelist input: 05-26-09 JDH
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
c
      parameter (nbelts=18)
      parameter (pi=3.14159265,grav=6.6732e-8,msun=1.9891e33)
      parameter (mp=1.67e-24,q0=1360.,cnvg=1.e-1)
      parameter (twopi=2*pi)
c   
      common/geogblk/cntrlat(nbelts),cntrlatangle(nbelts),coastlat,
     &   focean(nbelts),lat(0:nbelts+1),latangle(0:nbelts+1),ocean
c
      dimension  temp(0:nbelts+1),c(nbelts),dx(0:nbelts),
     &  tprime(0:nbelts),tprimeave(nbelts),t2prime(nbelts),
     &  x(0:nbelts+1),h20alb(0:90),ir(nbelts),
     &  mu(nbelts),s(nbelts),atoa(nbelts),surfalb(nbelts),
     &  acloud(nbelts),zntempmin(nbelts),zntempmax(nbelts),area(nbelts),
     &  zntempsum(nbelts),zntempave(0:nbelts+1),zndecmax(nbelts),
     &  zndecmin(nbelts),obstemp(nbelts),iceline(0:5),fice(nbelts),
     &  wthrate(nbelts),warea(nbelts),imco2(nbelts)
c
      character  header*80,file(0:3)*8
      logical  seasons, last, snowball
      real landsnowfrac
      data  temp/20*288./       !**use some temperature > 265K, otherwise iceball
      !data temp/20*233./ 
      data  file/ 'spng.out', 'summ.out', 'fall.out', 'wint.out' /

c  NAMELIST PARAMETERS
      INTEGER namelistid
      OPEN( namelistid, FILE='../input.nml', DELIM='APOSTROPHE' )
c

c  INITIALIZE VARIABLES
      ann_tempave = 10.         !just some constant greater than zero
      seasons = .true.     !if .true. do seasons: otherwise do mean-annual
      last = .false.       !should always be .false.
      snowball = .false.   !if .true. then start as an ice-ball
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
      avemol0 = 28.965*mp  !mean-molecular weight of present atmosphere (g)
      hcp0 = 0.2401    !heat capacity of present atmosphere (cal/g K)
      hcpco2 = 0.2105  !heat capacity of co2 (cal/g K)
      landsnowfrac = 1.0 !snowfall fraction on land
      fcloud = 0.5     !fractional cloud cover
      nfile = 0


      NAMELIST /ebm/ seasons, snowball, tend, dt, rot, a, ecc, peri, 
     &               obl, cloudir, pco2, ocean, igeog, groundalb, 
     &               relsolcon, landsnowfrac, fcloud
      READ( namelistid, NML=ebm )
      CLOSE( namelistid )


      oblrad = obl*pi/180
      if( snowball ) then
        do k = 1, nbelts, 1
          temp(k) = 233.
        end do
      end if
c 
c  OPEN FILES 
      open (unit=2,file='data/oceans.dat',status='old')
      open (unit=3,file='data/fresnel_reflct.dat',status='old')
      open (unit=15,file='out/model.out',status='unknown')
      open (unit=16,file='out/seasons.out',status='unknown')
      open (unit=17,file='out/co2clouds.out',status='unknown')
c
c  WRITE OBLIQUITY TO OUTPUT
      write(15,2) 'OBLIQUITY: ', obl, ' degrees'
 2    format(/ a,f5.2,a)
c
c  WRITE CO2-CLOUD FILE HEADER
      write(17,3)
 3    format(58x,'IMCO2')
      write(17,4)
 4    format(2x,'solar dec.(deg)',2x,'-85',1x,'-75',1x,'-65',
     &  1x,'-55',1x,'-45',1x,'-35',1x,'-25',1x,'-15',2x,'-5',3x,'5',2x,
     &  '15',2x,'25',2x,'35',2x,'45',2x,'55',2x,'65',2x,'75',2x,'85',
     &  2x,'global cloud coverage' /)
c
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
c
c SET UP GEOGRAPHY
      call geog(igeog)
c
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
c
c **calculate grid spacing
      do 125 k = 0,nbelts,1
         dx(k) = abs(x(k+1) - x(k))
 125  continue
c
c  WRITE GRID DATA TO OUTPUT FILE
c
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
c
c----------------------------------------------------------------------c
c  CALCULATE BELT AREAS - NORMALIZED TO PLANET SURFACE AREA
      do 180 k = 1,nbelts,1
         area(k) = abs(sin(lat(k) + pi/(2*nbelts)) -
     &     sin(lat(k) - pi/(2*nbelts)))/2.
 180  continue

c
c----------------------------------------------------------------------c
c  BEGIN INTEGRATION AT VERNAL EQUINOX. FIRST CALCULATE TIME SINCE
c  PERIHELION.
c
      d = d0  !**initialize diffusion coefficient
      w = (grav*msun)**(0.5)*a**(-1.5)  !2*pi/(orbital period)

c **write to 'co2clouds.out', at most, 1000 times.
      nwrite = amin0(1000,int((2*pi/w)/dt))
      twrite = (2*pi/w)/nwrite
 200  if(seasons) goto 230
      
      ! Do integration for mean-annual case with variable obliquity
      ! taken from Ward.  Integration is done in the insolation 
      ! section and uses a function definied at the bottom of this
      ! code (JDH)
      !obl = 0.
      !dec = 0.
      q = q0 * relsolcon
      goto 260
c
c  SEASONS
 230  e = acos((cos(peri*pi/180.)+ecc)/(1.+cos(peri*pi/180.)*ecc))
      m = e - ecc*sin(e)   !**Kepler's Equation
      tperi = -m/w    !**tperi is the time of perihelion (vern.eqnx. is 0)
c
c      write(15,232)
 232  format(/ 'ORBITAL DATA')
c      write(15,233) 
 233  format(2x,'time (sec)',2x,'true anomaly (deg)',2x,
     & 'solar declination (deg)',2x,'orbital distance (cm)',2x,
     & 'irradiance (Wm^-2)')
c
 240  if (tcalc.lt.tend) then
         t = t + dt
         tcalc = tcalc + dt
         tlast = tlast + dt
      else
         write(*,*) 'Calculation time has elapsed.'
         goto 1000
      end if
c
      m = w*(t-tperi)
      if (m.gt.2*pi) m = m - 2*pi
c
      call keplereqn(m,ecc,e)
      trueanom = acos((cos(e) - ecc)/(1 - ecc*cos(e)))
c
c  correct trueanomaly for pi < m < 2pi
      if (m .gt. pi) trueanom = 2*pi-trueanom
c
      r = a*(1-ecc*cos(e))
      q = q0 * relsolcon * (a/r)**2
      !q = q0 * relsolcon * (a0/r)**2
      thetawin = peri*pi/180 + 3*pi/2
      dec = asin(-sin(obl*pi/180.)*cos(trueanom-thetawin))
      decangle = dec*180./pi
c
c      write(15,255) t,trueanom*180./pi,decangle,r,q
 255  format(2x,e9.3,6x,f7.3,14x,f7.3,16x,e11.6,13x,f9.3)
c
      goto 270      !**skip mean-annual stuff if doing seasonal model
c
c----------------------------------------------------------------------c
c MEAN-ANNUAL CALCULATION
c
 260  if (tcalc.lt.tend) then
         t = t + dt
         tcalc = tcalc + dt
         tlast = tlast + dt
      else
         write(*,*) 'Calculation time has elapsed.'
         goto 1000
      end if
c
 270  tempsum = 0.
      fluxsum = 0.
      albsum = 0.
      irsum = 0.
      globwthrate = 0.
      co2cldsum = 0.
c
c----------------------------------------------------------------------c
c  **THE BIG LOOP**
c----------------------------------------------------------------------c
c  FINITE DIFFERENCING - ATMOSPHERIC and OCEANIC ADVECTIVE HEATING 
c
      do 300 k=0,nbelts,1   !**first derivatives between grid points
         tprime(k) = (temp(k+1) - temp(k))/dx(k)
 300  continue
c
      do 310 k=1,nbelts,1   !**start belt loop 
c                           !**first derivatives at grid points
         tprimeave(k) = (tprime(k)*dx(k) + tprime(k-1)*
     &      dx(k-1))/(dx(k) + dx(k-1))
c                           !**second derivatives at grid points
         t2prime(k) = ((((dx(k)/2)**2)-((dx(k-1)/2)**2))*(1-x(k)**2)*
     &      tprimeave(k) + ((dx(k-1)/2)**2)*(1-(x(k)+(dx(k)/2))**2)*
     &      tprime(k) - (((dx(k)/2)**2)*(1-(x(k)-(dx(k-1)/2))**2))*
     &      tprime(k-1))/((dx(k)/2)*(dx(k-1)/2)*(dx(k)/2 + dx(k-1)/2)) 
c
c----------------------------------------------------------------------c
c  OUTGOING INFRARED (Obtained from fits to Jim Kasting's radiative
c    -convective model earthtem_ir. The fit can be applied for
c    190 < T < 370 K and 10^-5 < PCO2 < 10 bars with an rms error
c    of 4.59 Wm^-2.)                                 
c
c      if ( temp(k) .le. 190.0 ) then
c        print *, k, ": ", temp(k), " Temp too low"
c      else if ( temp(k) .ge. 370.0 ) then
c        print *, k, ": ", temp(k), " Temp too high"
c      end if

      phi = log(pco2/3.3e-4)
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
c
      ir(k) = ir(k) - cloudir   !**reduction of outgoing-infrared by clouds
c
c----------------------------------------------------------------------c
c  DIURNALLY-AVERAGED ZENITH ANGLE                            
c
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
c  
      z = acos(mu(k))*180./pi
      zrad = acos(mu(k))
c                 
c----------------------------------------------------------------------c
c  HEAT CAPACITY and SURFACE ALBEDO                    
c
      oceanalb = h20alb(int(z))


      ! Added this check so fractional ice cover would only 
      ! occur for 263 < T < 273.  (JDH)
c      if(focean(k).ge.0.5) then        
c         if (temp(k).ge.238.15) then
c            fice(k) = 0
c         else if (temp(k).le.238.15) then
c            fice(k) = 1
c         else
c            !fice(k) = 1. - exp((temp(k)-273.15)/10.)
c            fice(k) = 0.1*(248.15 - temp(k)) ! uncomment for a linear decrease
c         end if 
c      else
         if (temp(k).ge.273.15) then
            fice(k) = 0
         else if (temp(k).le.263.15) then
            fice(k) = 1
         else
            fice(k) = 1. - exp((temp(k)-273.15)/10.)
         end if
c      end if   
 

      acloud(k) = alpha + beta*zrad
      if (temp(k).le.263.15) goto 420
      if (temp(k).le.273.15) goto 420 !removed variable snow albedo (JDH)
      !if (temp(k).le.273.15) goto 410
c
      landalb = groundalb
      fice(k) = 0.
      fwthr = 1.       !**100% if T > 273K, 0% otherwise
      c(k) = focean(k)*cw + (1 - focean(k))*cl
      goto 430
c
c LAND WITH UNSTABLE SNOW COVER; SEA-ICE FORMS
 410  snowalb = 0.45
      landalb = snowalb
      icealb = 0.55
      ci = 4.846e7
      c(k) = (1-fice(k))*focean(k)*cw + fice(k)*focean(k)*ci + 
     &  (1 - focean(k))*cl
      fwthr = 0.
      goto 430
c
c LAND WITH STABLE SNOW COVER; SEA-ICE WITH SNOW COVER
 420  snowalb = 0.663  ! changed this from 0.7 as per Caldiera & Kasting (JDH)
      landalb = snowalb*landsnowfrac + groundalb*(1 - landsnowfrac)
      icealb = snowalb

      ci = 1.05e7
      c(k) = (1-fice(k))*focean(k)*cw + fice(k)*focean(k)*ci + 
     &  (1 - focean(k))*cl
      fwthr = 0.
c
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

c      surfalb(k) = (1-fcloud)*((1-focean(k))*landalb +   
c     &     focean(k)*((1-fice(k))*oceanalb + fice(k)*icealb)) + 
c     &     fcloud*acloud(k)
      
c
      if(.not.last) goto 450
c
c----------------------------------------------------------------------c
c  CO2 CLOUDS, IF LAST LOOP
      tstrat = -188.0968076496237 - 1.954500243043531*phi +
     &  3.809523755282467*phi**2 + 2.327695124784982*temp(k) + 
     &  0.0003732556589321478*phi*temp(k) -
     &  0.02855747090737606*phi**2*temp(k) -
     &  0.00332927322412119*temp(k)**2 + 
     &  0.00002213983847712599*phi*temp(k)**2 + 
     &  0.0000460529090383333*phi**2*temp(k)**2
c
      nh = 0
      tlapse = 6.5    !**degrees Kelvin per kilometer
 442  htemp = temp(k) - tlapse*nh
      if (htemp.lt.tstrat) goto 450
c
      hpco2 = pco2*((temp(k)-tlapse*nh)/temp(k))**
     &  (avemol*981./(1.38e-16*tlapse*1.e-5))
c
c  SUBROUTINE SATCO2 from Jim Kasting "earthtem_0.f"
      if (htemp.lt.216.56) goto 445
c
c   VAPOR PRESSURE OVER LIQUID
      psl = 3.128082 - 867.2124/htemp + 1.865612e-2*htemp -
     &  7.248820e-5*htemp**2 + 9.3e-8*htemp**3
c
c   VAPOR PRESSURE OVER SOLID
c   (Altered to match vapor pressure over liquid at triple point)
 445  psl = 6.760956-1284.07/(htemp-4.718) + 1.256e-4*(htemp-143.15)
c
      patm = 10.**psl
      psat = 1.013*patm
c
      imco2(k) = 0.
      if (hpco2.ge.psat) then      !**co2 clouds form
         imco2(k) = 1.
         goto 450
      end if
      nh = nh + 1
      goto 442
c
c----------------------------------------------------------------------c
c  CARBONATE-SILICATE WEATHERING
 450  warea(k) = area(k)*(1-focean(k))*fwthr
      wthrate(k) = warea(k)*(1+0.087*(temp(k)-288.) + 
     &   1.862e-3*(temp(k)-288.)**2)
c
c----------------------------------------------------------------------c
c  TOP-OF-ATMOSPHERE ALBEDO (Obtained from fits to Jim Kasting's radiative
c    -convective model 'earthtem_alb.' Both fits can be applied for
c    10^-5 < pco2 < 10 bars, 0 < surfalb < 1, and 0 < z < 90 degrees
c    with r.m.s. errors (given in planetary-average incident solar
c    flux {340 W/m^-2}) of 7.58 and 4.66 Watts/m^2 for the low and
c    high temperature regimes respectively.
c
      as = surfalb(k)
      ! Added this check for out of range temperatures (JDH)
c      if(temp(k).ge.370) then
c         print *, 'Temp too high!:',k, temp(k)
c         goto 530
c      end if
      if(temp(k).ge.280.) goto 520   !**goto high-temp fit
c
      atoa(k) = -6.891041e-1 + 1.046004*as + 7.323946e-2*as**2 - 
     &   2.889937e-1*mu(k)+2.012211e-1*as*mu(k)+8.121773e-2*mu(k)**2 -  
     &   2.837280e-3*pco2 - 3.741199e-2*as*pco2 - 
     &   6.349855e-3*mu(k)*pco2 + 6.581704e-4*pco2**2 + 
     &   7.805398e-3*temp(k) - 1.850840e-3*as*temp(k) + 
     &   1.364872e-4*mu(k)*temp(k) + 9.858050e-5*pco2*temp(k) - 
     &   1.655457e-5*temp(k)**2
      goto 530
c
 520  atoa(k) = 1.108210 + 1.517222*as + 7.588651e-2*as**2 - 
     &   1.867039e-1*mu(k)+2.098557e-1*as*mu(k)+6.329810e-2*mu(k)**2 +  
     &   1.970523e-2*pco2 - 3.135482e-2*as*pco2 - 
     &   1.021418e-2*mu(k)*pco2 - 4.132671e-4*pco2**2 - 
     &   5.79932e-3*temp(k) - 3.709826e-3*as*temp(k) - 
     &   1.133523e-4*mu(k)*temp(k) + 5.371405e-5*pco2*temp(k) + 
     &   9.269027e-6*temp(k)**2
 530  continue
c
c----------------------------------------------------------------------c
c  DIURNALLY-AVERAGED INSOLATION 
c      
      if (seasons) then
         s(k) = (q/pi)*(x(k)*sin(dec)*h + 
     &        cos(asin(x(k)))*cos(dec)*sin(h))
      else
         ! for mean annual do the Ward isolation integration (JDH)
         call qtrap(0, twopi, s(k), lat(k), oblrad)
         s(k) = (q/(2*(pi**2))) * 1/(sqrt(1 - ecc**2)) * s(k)
         !print *, 'insolation at', latangle(k), ':', s(k)
      end if
c
c----------------------------------------------------------------------c
c  SURFACE TEMPERATURE - SOLVE ENERGY-BALANCE EQUATION 
c
      temp(k) = (d*t2prime(k)-ir(k)+s(k)*(1-atoa(k)))*dt/c(k) + temp(k)
c
c  SUM FOR GLOBAL AVERAGING
      irsum = irsum + area(k)*ir(k)
      fluxsum = fluxsum + area(k)*s(k)
      albsum = albsum + area(k)*s(k)*atoa(k)
      tempsum = tempsum + area(k)*temp(k)
      globwthrate = globwthrate + wthrate(k)
      co2cldsum = co2cldsum + area(k)*imco2(k)*s(k)
c
      if(.not.last) goto 310
c  ZONAL STATISTICS - if last loop 
      zntempmin(k) = amin1(zntempmin(k),temp(k))
      if (zntempmin(k).eq.temp(k)) zndecmin(k) = decangle
      zntempmax(k) = amax1(zntempmax(k),temp(k))
      if (zntempmax(k).eq.temp(k)) zndecmax(k) = decangle
      zntempsum(k) = zntempsum(k) + temp(k)
c
 310  continue                                     !**end of belt loop
c  **set pole temps equal to adjacent belt temps
      temp(0) = temp(1)
      temp(nbelts+1) = temp(nbelts)
c
c  WRITE OUTPUT FILES - ONLY ON LAST LOOP      
      if(.not.last ) goto 710
      if(tlast.lt.twrite) goto 710
      tlast = 0.
c
c  EQUINOXES AND SOLSTICES
c
      ! added seasons check here (JDH)
      if (seasons) then
         decold = decmid
         decmid = decnew
         decnew = decangle

         if (((abs(decmid).gt.(abs(decold))).and.
     &    (abs(decmid).gt.(abs(decnew)))).or.(decnew*decmid.le.0.)) then 
c     
c     write(15,610) file(nfile),decangle
 610        format('data written to file ',a8, ' at declination ',f6.2,
     &           ' degrees')
c            print *, 'nfile=', nfile
c            open(unit=nfile+7,file=file(nfile),status='unknown')
c            do 620 k = 1,nbelts,1
c               write(nfile+7,615) latangle(k),temp(k),atoa(k),ir(k),
c     &              d*t2prime(k),s(k)*(1-atoa(k))
 615           format(f6.2,5(3x,f7.3))
 620        continue
            nfile = nfile + 1
         end if         
      end if
c
c  CO2-CLOUD DATA
      write(17,630) decangle,imco2,co2cldsum/fluxsum
 630  format(3x,f6.2,9x,18(3x,i1),9x,f5.3)
c     
c----------------------------------------------------------------------c
c  GLOBAL AVERAGING
c     
 710  irave = irsum
      fluxave = fluxsum
      co2cldave = co2cldsum/fluxave
      albave = albsum/fluxave
      tempave = tempsum
c
      iravesum = iravesum + irave
      fluxavesum = fluxavesum + fluxave
      albavesum = albavesum + albave
      tempavesum = tempavesum + tempave
      wthratesum = wthratesum + globwthrate
      co2cldavesum = co2cldavesum + co2cldave
      nstep = nstep + 1
c
c  SEASONAL AVERAGING
      if(t.lt.2*pi/w) goto 800    !**one orbit since last averaging
      ann_tempave = tempavesum/nstep
      ann_albave = albavesum/nstep
      ann_fluxave = fluxavesum/nstep
      ann_irave = iravesum/nstep
      ann_wthrateave = wthratesum/nstep
      ann_co2cldave = co2cldavesum/nstep
c
c      write(*,*) ann_tempave,pco2,d


c  ADJUST PCO2 EVERY 5 SEASONAL CYCLES
      if(nwthr.lt.5) goto 730
      nwthr = 0
      wthr = wco2*ann_wthrateave
      if (wthr.lt.v) goto 725
c       pco2 = pco2*(1+exp((v-wthr)/v))/2.  ! removed co2 adjusting (JDH)
      goto 730
 725  if(wthr.lt.1.e-2) wthr = 1.
c       pco2 = 2*pco2*(1-0.5*exp((wthr-v)/wthr))
c
c  ADJUST DIFFUSION COEFFICIENT
      avemol = mp*(28.965+44.*pco2)/(1.+pco2) 
      hcp = (hcp0 + pco2*hcpco2)/(1.+pco2)
      write(*,*) pco2,d
      d = d0*(1. + pco2)*((avemol0/avemol)**2)*(hcp/hcp0)*
     &   (rot0/rot)**2
 730  nwthr = nwthr + 1
c
      if(.not.last) goto 790  
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
      print *, "SH/NH temperature difference = ", shtempave - nhtempave
      zntempave(nbelts+1) = zntempave(nbelts)  !**for ice-line calculation
c
c  FIND ICE-LINES (ANNUAL-AVERAGE TEMP < 263.15K)
      nedge = 0
      do 740 k=1,nbelts,1
      if ((zntempave(k+1)-263.15)*(zntempave(k)-263.15) .lt. 0.) then
         icelat = latangle(k) + ((latangle(k+1)-latangle(k))/
     &   (zntempave(k+1)-zntempave(k)))*(263.15-zntempave(k))
         nedge = nedge + 1
         iceline(nedge) = icelat
      end if
c      write(*,*) nedge,iceline
 740  continue
c
      do 750 k=1,nbelts,1
         write(15,751) latangle(k),zntempave(k),zntempmin(k),
     &      zndecmin(k),zntempmax(k),zndecmax(k),
     &      (zntempmax(k)-zntempmin(k))/2.
 751     format(4x,f4.0,9x,f8.3,5x,f8.3,5x,f6.2,5x,f8.3,5x,f6.2,5x,f8.3)
         write(16,752) latangle(k),(zntempmax(k)-zntempmin(k))/2.,
     &     10.8*abs(x(k)),15.5*abs(x(k)),6.2*abs(x(k))
 752     format(4x,f4.0,4x,4(3x,f8.3))
 750  continue
c      
c      write(15,755)
 755  format(/ 'SURFACE DATA')
c      write(15,756)
 756  format(2x,'latitude(deg)',2x,'temp(k)',2x,'belt area',
     &  2x,'weathering area',2x,'zonal weathering rate (g/yr)')
      do 758 k = 1,nbelts,1
c         write(15,757) latangle(k),temp(k),area(k),warea(k),
c     &     wthrate(k)
 757     format(6x,f4.0,8x,f8.3,10x,f5.3,10x,f5.3,15x,e8.2)
 758  continue
c
       write(15,760)
 760  format(/ 'ICE LINES (Tave = 263K)')
      if((nedge.eq.0).and.(zntempave(nbelts/2).le.263.)) then
      	write(15,*) '  planet is an ice-ball.' 
      else if((nedge.eq.0).and.(zntempave(nbelts/2).gt.263.)) then
     	write(15,*) '  planet is ice-free.'
      else
        do 765 k=1,nedge,1
           write(15,762) iceline(k)
 762       format(2x,'ice-line latitude = ',f5.1,' degrees.')
 765    continue
      end if
c
c  CO2 CLOUDS
      write(15,770)
 770  format(/ 'CO2 CLOUDS')
      write(15,772) ann_co2cldave
 772  format(2x,'planet average co2 cloud coverage = ',f5.3)
c
c  GEOGRAPHY INFORMATION 
c
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
c
 790  fluxavesum = 0.
      albavesum = 0.
      tempavesum = 0.
      iravesum = 0.
      wthratesum = 0.
      t = 0
      nstep = 0
      if(last) stop
c
c check for convergence
c      if((nwthr.ge.4).and.(abs(wthr-v)/v.lt.cnvg)) goto 1000
      if(abs(prevtempave-ann_tempave).lt.cnvg) goto 1000
      prevtempave = ann_tempave      
      !print *, "tempave:", ann_tempave
c
 800  if(seasons) goto 240
      goto 260
c
c----------------------------------------------------------------------c
c  WRAP UP
c
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
      write(15,1016) 'seasonal-average weathering rate = ',
     & ann_wthrateave*wco2,' g/year'
      write(15,1016) 'co2 partial pressure = ', pco2, ' bars'
 1015 format(3x,a,f7.3,a)
      write(15,1016) 'thermal diffusion coefficient (D) = ', d, 
     & ' Watts/m^2 K'
 1016 format(3x,a,e8.3,a)
      write(15,1020) 'convergence to final temperature profile in ',
     & tcalc, ' seconds'
 1020 format(3x,a,e9.3,a)
c
c  LOOP ONE MORE TIME, AND WRITE OUTPUT TO FILES
      last = .TRUE.
c
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
     &  'seas.amp(K)')
c
c      write(15,1100) 
 1100 format(/ 'OUTPUT FILES')
      tcalc = 0.
      if(seasons) goto 240
      goto 260
c
      end
c
c----------------------------------------------------------------------c
      subroutine  geog(igeog)
c
c  This subroutine calculates continental/oceanic zone fractions for five
c  geography types (1) present geography (2) polar supercontinent (3) equa-
c  torial supercontinent (4) 100% oceanic and (5) 100% land (10-27-95)
c
      parameter (pi=3.141592653589793,nbelts=18)
c
      common/geogblk/cntrlat(nbelts),cntrlatangle(nbelts),coastlat,
     &   focean(nbelts),lat(0:nbelts+1),latangle(0:nbelts+1),ocean
c
      real*4  lat,latangle
      logical coast
c
      goto (10,20,30,40,50) igeog
c      
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
c
      return
c
c
c  SOUTH POLAR SUPERCONTINENT 
 20   coast = .TRUE.    !** still on land
      coastlat = asin(1-2*ocean)*180/pi
      do 25 k = 1,nbelts,1
         if (latangle(k).le.coastlat) then
            focean(k) = 0.
c
c **see whether you are on the coast
            if (latangle(k) .eq. coastlat) coast = .FALSE.
c  
         else if ((latangle(k).gt.coastlat).and.(coast)) then
            focean(k) = (sin(lat(k)) - sin(coastlat*pi/180.))/
     &         (sin(lat(k)) - sin(lat(k)-pi/nbelts))          
            coast = .FALSE.     !** now in water
         else
            focean(k) = 1.
         end if
 25   continue
c
      return
c
c  EQUATORIAL SUPERCONTINENT
 30   coast = .TRUE.     !**still on land
      coastlat = asin(1-ocean)*180/pi
      do 35 k = nbelts/2+1,nbelts,1
         if (latangle(k) .le. coastlat) then
            focean(k) = 0.
            focean(nbelts-k+1) = 0.
c
c **see whether you are on the coast
            if (latangle(k) .eq. coastlat) coast = .FALSE. 
c             
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
c
      return
c
c  100% WATER COVERED***
 40   do 45 k = 1,nbelts,1
         focean(k) = 1.
 45   continue
      return
c
c  100% LAND SURFACE***
 50   do 55 k = 1,nbelts,1
         focean(k) = 0.
 55   continue
      return
c
      end
c
c---------------------------------------------------------------------c
      subroutine  keplereqn(m,e,x)
c---------------------------------------------------------------------c
c  This subroutine is used to find a solution to Keplers' equation,
c  namely M = E - eSin(E), where M is the mean anomaly, E is the
c  eccentric anomaly, and e is the eccentricity of the orbit.
c  Input is M and e, and output is E. The chosen convergence is quartic 
c  as outlined in Danby (1988), pp. 149-54.  coded from Danby (5-11-95)
c
      implicit  integer(n), real*8(a-m,o-z)
      parameter (pi=3.141592653589793)
c
c  initial guess is x = M + k*e with k = 0.85
      k = 0.85
c  bound on the permissible error for the answer
      del = 1.e-13
c
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


c---------------------------------------------------------------
c     Integrand from Ward, 1974
c---------------------------------------------------------------
      FUNCTION insolation(x,d,t)
      REAL :: insolation
      REAL, INTENT(IN) :: x,d,t
      insolation = SQRT(1 - (sin(d)*cos(t) - cos(d)*sin(t)*sin(x))**2)
      END FUNCTION


c------------------------------------------------------------------
c     Subroutines for numerical integration from Numerical Recipes
c------------------------------------------------------------------
      
      SUBROUTINE trapzd(a,b,s,n,lat,ob) 
      INTEGER n 
      REAL a,b,s,insolation,lat,ob
      EXTERNAL insolation
      INTEGER it,j 
      REAL del,sum,tnm,x 
      if (n.eq.1) then 
         s=0.5*(b-a)*(insolation(a,lat,ob)+insolation(b,lat,ob)) 
      else 
         it=2**(n-2) 
         tnm=it 
         del=(b-a)/tnm 
         x=a+0.5*del 
         sum=0. 
         do 1701 j=1,it 
            sum=sum+insolation(x,lat,ob) 
            x=x+del 
 1701    continue
         s=0.5*(s+(b-a)*sum/tnm)         
      endif 
      return 
      END


      SUBROUTINE qtrap(a,b,s,lat,ob) 
      INTEGER JMAX 
      REAL a,b,s,EPS,lat,ob 
      PARAMETER (EPS=1.e-6, JMAX=20) 
      INTEGER j 
      REAL olds olds=-1.e30 
      do 1711 j=1,JMAX 
         call trapzd(a,b,s,j,lat,ob) 
         if (j.gt.5) then 
            if (abs(s-olds).lt.EPS*abs(olds).or. 
     &           (s.eq.0..and.olds.eq.0.)) return 
         endif 
         olds=s 
 1711 continue
      pause 'too many steps in qtrap'
      END
