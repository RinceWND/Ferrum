!
!
      subroutine bulk(im,jm,tbias,fsm,tsurf,alon,alat,iyr,imo,iday,
     $    ihour,imin,wusurf,wvsurf,wtsurf,swrad,pme,
     $                                             uair,vair,
     $                                    tair,rhum,rain,cloud,pres)
!*******************************************************************************
! THIS SUBROUTINE PROVIDES SURFACE BOUNDARY CONDITIONS FOR MOMENTUM,
! HEAT AND SALT EQUATIONS SOLVED BY THE HYDRODYNAMIC MODEL IN CASES WHEN
! THE ATMOSPHERIC MODEL PROVIDES CLOUD COVER DATA INSTEAD OF THE NET SOLAR
! RADIATION FLUX AND THE DOWNWARD LONGWAVE RADIATION FLUX AS IT IS ASSUMED
! IN THE FIRST VERSION OF BULK CODE. THE NET SOLAR RADIATION IS CALCULATED
! ACCORDING TO THE REED FORMULA WHILE THE NET LONGWAVE RADIATION CAN BE
! CALCULATED ACCORDING TO BIGNAMI OR MAY FORMULA (SEE LOGICAL VARIABLES
! BIGNAMI_FORMULA & MAY_FORMULA BELOW)
!
! MOMENTUM, HEAT AND FRESHWATER FLUXES ARE CALCULATED FROM THE ATMOSPHERIC
! PARAMETERS (WIND VELOCITY, AIR TEMPERATURE, RELATIVE HUMIDITY, PRECIPITATION,
! CLOUD COVERAGE) PROVIDED BY THE
! WEATHER PREDICTION MODEL AND THE MODEL'S SST USING PROPER AIR-SEA BULK
! FORMULAE (Castellari et al., 1998, Korres and Lascaratos 2003).
! ( Castellari et al., 1998. Journal of Marine Systems, 18, 89-114 ;
!   Korres and Lascaratos, 2003. Annales Geophysicae, 21, 205-220.)
!
!  ALL UNITS ARE S.I. (M.K.S.)
!
! THE USER HAS THE OPTION TO CALCULATE THE NET LONGWAVE RADIATION FLUX
! ACCORDING TO BIGNAMI OR MAY FORMULA. THIS IS DONE THROUGH THE LOGICAL
! VARIABLES BIGNAMI_FORMULA & MAY_FORMULA
!
! SUBROUTINE BULK PROVIDES ITS OUTPUT INTO THE ARRAYS:
!
! 1. WUSURF(): X-COMPONENT OF WIND STRESS DIVIDED BY (-1)*RHO
! 2. WVSURF(): Y-COMPONENT OF WIND STRESS DIVIDED BY (-1)*RHO
! 3. SWRAD() : NET SOLAR RADIATION HEAT FLUX DIVIDED BY (-1)*RHO*Cpw
! 4. WTSURF(): NET UPWARD HEAT FLUX DIVIDED BY RHO*Cpw
! 5. PME()   : PRECIPITATION MINUS EVAPORATION RATE
! ( RHO: sea water density, Cpw: specific heat capacity of seawater )
!
!
! SUBROUTINE BULK NEEDS THE FOLLOWING INPUT:
!
!
! MODEL RELATED DATA
! 1. IM     : NUMBER OF GRID POINTS IN X
! 2. JM     : NUMBER OF GRID POINTS IN Y
! 3. TBIAS  : CONSTANT VALUE SUBTRACTED FROM MODEL'S TEMPERATURE FIELD
! 4. FSM()  : THE MODEL MASK (1.:SEA GRID POINT, 0.:LAND GRID POINT)
! 5. TSURF(): MODEL'S TEMPERATURE FIELD AT THE TOP VERTICAL LEVEL AND AT THE
!             CENTRAL TIME LEVEL
! 6. ALON() : LOGNITUDE OF GRID POINTS
! 7. ALAT() : LATITUDE OF GRID POINTS
! 9. IYR    : INTEGER VALUE CORRESPONDING TO CURRENT YEAR(for example iyr=2002)
! 10. IMO   : INTEGER VALUE CORRESPONDING TO CURRENT MONTH (1-->12)
! 11. IDAY  : INTEGER VALUE CORRESPONDING TO CURRENT DAY (1-->31)
! 12. IHOUR : INTEGER VALUE CORRESPONDING TO CURRENT HOUR (0-->23)
! 13: IMIN  : INEGER VALUE CORRESPONDING TO CURRENT MINUTES (0--59)
!
! ATMOSPHERIC DATA:
! 6.  UAIR() : X-COMPONENT OF AIR VELOCITY (in m/s) AT 10m ABOVE SEA SURFACE
! 7.  VAIR() : Y-COMPONENT OF AIR VELOCITY (in m/s) AT 10m ABOVE SEA SURFACE
! 8.  TAIR() : AIR TEMPERATURE (in deg Kelvin) AT 2m ABOVE SEA SURFACE
! 9.  RHUM() : RELATIVE HUMIDITY (%) AT 2m ABOVE SEA SURFACE
! 10. RAIN() : PRECIPITATION RATE (in m/s)
! 11. CLOUD()  : CLOUD COVERAGE IN TENTHS (0.-->1.)
! 12. PRES() : ATMOSPHERIC PRESSURE AT SURFACE (hPa)
!
!
! Important:
! SUBROUTINE BULK REQUIRES THAT THE ATMOSPHERIC DATA ARE ALREADY
! INTERPOLATED IN TIME AND SPACE (i.e. mapped onto the model grid and
! interpolated to the current time step. The user has to write his/her own
! subroutines in order to map the raw atmospheric data onto the model grid
! and interpolate them to the current time step)
!
!*******************************************************************************
!      implicit none

      include 'realkind'

      logical BIGNAMI_FORMULA, MAY_FORMULA, herzfeld_formula
      parameter (BIGNAMI_FORMULA=.true., MAY_FORMULA=.false.,
     $           herzfeld_formula=.false.)

      real(kind=rk), dimension(im, jm) ::
     $ fsm(im,jm),pme(im,jm),swrad(im,jm),wusurf(im,jm),
     $ wvsurf(im,jm),wtsurf(im,jm),tsurf(im,jm),alon(im,jm),
     $ alat(im,jm)

      real(kind=rk), dimension(im, jm) ::
     $  uair(im,jm),vair(im,jm),tair(im,jm),rhum(im,jm),
     $  rain(im,jm),cloud(im,jm),pres(im,jm)

      real(kind=rk) unow, vnow, tnow, pnow, e, es, precip, cld
     &   , sst_model, QBW

      real(kind=rk), external :: cd, heatlat, esk

      real(kind=rk) arho
     &            , cd1, ce2, ch2
     &            , deltemp
     &            , ea12, emiss, esatair, esatoce, evap
     &            , fe, fh
     &            , Qe, Qh, Qu
     &            , rhnow, rhom
     &            , sigma, sol_net, sp, ss, sstk, stp
     &            , taux, tauy, tnowk
     &            , wair, wflux, wsatair, wsatoce
!
!--------------------------------------------------------------------
!       coefficients ( in MKS )  :
!-----------------------------------------------------------------
!
! --- Sea water density

      data rho/1023./

! --- Sea water density times the specific heat of seawater

      data rho_cpw/4.082793d6/

! --- surface air pressure, expsi, dry air gas constant
!
      data ps,expsi,rd / 1013., 0.622, 287./
!
! --- turbulent exchange coefficients ( from Budyko 1963 )
!
      data  ce2,ch2  / 2.1d-3, 2.1d-3/
!
! --- air density, Stefan-Boltzmann constant , ocean emissivity
!
      data arho,sigma,emiss   /1.2,   5.67d-8, .97/ ! 5.669e-8
!
! --- Solar constant , specific heat capacity of air
!
      data solar ,cp   /1350., 1005./
!
!
      data ckelv /273.16/
!
      const = 0.622/1013.
!
!
      do j = 1,jm
        do i = 1,im
          if (fsm(i,j) /= 0.) then
            unow      = uair(i,j)
            vnow      = vair(i,j)
            tnow      = tair(i,j)+ckelv
            pnow      = pres(i,j)
            e         = rhum(i,j)*pnow / ( 0.378*rhum(i,j) + 0.622 )
            es        = 6.112 *
     &           exp( (( 17.67*tair(i,j) )/( tair(i,j) + 243.25 )) )
            rhnow     = 100.*e/es ! rwnd: specific to relative humidity
            if (rhnow>100.) rhnow=100.
            if (rhnow<0.) rhnow=0.
            precip    = rain(i,j)/1000. ! rwnd: precipitation rate from kg/(m2*s) to m/s
            cld       = cloud(i,j)/100. ! rwnd: total cloud cover from % to tenths
            sst_model = tsurf(i,j)+tbias

! SST_MODEL IS THE MODEL'S TEMPERATURE (in deg Celsius) AT THE TOP LEVEL
!
! --- compute wind speed magnitude for bulk coeff.
!
            SP = sqrt(unow*unow+vnow*vnow)

!
! --- SST data converted in Kelvin degrees
! --- TNOW is already in Kelvin
!
            sstk = sst_model + ckelv
            tnowk = tnow
!
!
! ---calculates the Saturation Vapor Pressure at air temp. and at sea temp.
! ---esat(Ta) , esat(Ts)
!
            esatair = esk(tnowk)
            esatoce = esk(sstk)
!
! --- calculates the saturation mixing ratios at air temp. and sea temp.
! --- wsat(Ta) , wsat(Ts)
!
            wsatair = (expsi/pnow) * esatair
            wsatoce = (expsi/pnow) * esatoce
!
! --- calculates the mixing ratio of the air
! --- w(Ta)
!
            wair = 0.01 * rhnow * wsatair
!
! --- calculates the density of  moist air
!
            rhom = 100.*(pnow/rd)*(expsi*(1.+wair)/(tnowk*(expsi+wair)))
!
!---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
! Calculate the net longwave radiation flux at the sea surface (QBW)
! according to Bignami (Bignami et al., 1995) or May formula (May,1986)
!
! Bignami et al., 1995: Longwave radiation budget in the Mediterranean
! Sea, J.Geophys Res., 100, 2501-2514.
!---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
!

            ea12 = 0.01*rhnow*esatair

            if (bignami_formula) then
              QBW=0.98*sigma*sstk**4-sigma*tnowk**4*(0.653+0.00535*ea12)
     $                 *(1.+0.1762*cld*cld)
            end if

            if (may_formula) then
              QBW = (1.-0.75*(cld**3.4))
     $             * (sigma*(tnowk**4.)*(0.4 -0.05*sqrt(ea12))
     $           + 4.*sigma*(tnowk**3.)*(sstk-tnowk))
            end if

            if (herzfeld_formula) then
              QBW = (sigma*0.96*(1-(0.92e-5*tnowk*tnowk))*tnowk**4 +
     $             4*sigma*0.96*(ckelv+tair(i,j)**3)*(sstk-tnowk)) *
     $              (1-cos(alat(i,j))*cld) ! cos(phi) here is an improvised `beta` coefficient as a function of latitude from Herzfeld
            end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! --- calculate the term : ( Ts - Ta )
            deltemp = sstk - tnowk
!
! Calculate turbulent exchange coefficients according to Kondo scheme
! ( Kondo, J., 1975. Boundary Layer Meteorology, 9, 91-112)

! ---- Kondo scheme
            ss = 0.
            fe = 0.
            fh = 0.
            if (sp > 0.0) ss = deltemp/(sp**2.)
!
! --- calculate the Stability Parameter :
!
            stp = ss*abs(ss)/(abs(ss)+0.01)
!
! --- for stable condition :
            if (ss < 0.) then
              if ((stp > -3.3).and.(stp < 0.)) then
                fh = 0.1+0.03*stp+0.9*exp(4.8*stp)
                fe = fh
              else
                if (stp <= -3.3) then
                  fh = 0.
                  fe = fh
                end if
              end if
!
! --- for unstable condition :
            else
              fh = 1.0+0.63*sqrt(stp)
              fe = fh
            end if
!
            ch2 = 1.3e-03*fh
            ce2 = 1.5e-03*fe

!---- --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
! Calculate the sensible (QH) latent (QE) heat flux
!                    and the evaporation rate (EVAP)
!---- --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
!
            QH = rhom*cp*ch2*sp*deltemp
!
! --- calculates the term : esat(Ts)-r*esat(Ta)
!
            EVAP = esatoce - rhnow*0.01*esatair
!
! --- calculate the term : Ce*|V|*[esat(Ts)-r*esat(Ta)]0.622/1013
! --- Evaporation rate [kg/(m2*sec)]
!
            EVAP = rhom*ce2*sp*evap*const
!
! --- calculate the LATENT HEAT FLUX  QE in MKS ( watt/m*m )
! --- QE = L*rhom*Ce*|V|*[esat(Ts)-r*esat(Ta)]0.622/1013
!
            QE = evap*heatlat(sst_model)
!
!
!---- --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
! Calculate the water flux (WFLUX) in m/sec
!---- -- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
!
            WFLUX =  evap/rho - precip

            pme(i,j)= -wflux
!
! Important note for Princeton Ocean Model users:
! THE SALT FLUX ( WSSURF() ) IN POM MAIN CODE (REQUIRED FOR PROFT) SHOULD
! BE CALCULATED AS:
!       do 3072 j = 1, jm
!       do 3072 i = 1, im
!  3072 WSSURF(I,J)=PME(I,J)*(VF(I,J,1)+SBIAS)

!---- --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
! Calculate  the net upward flux (QU) at the sea surface
!---- --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
!
! --- calculates : Qu = Qb + QH + QE
!
            QU = qbw + qh + qe
!
! --- 1. Devide upward heat flux by rho*Cpw
!
            wtsurf(i,j) = QU/rho_cpw

! --- Calculate net solar radiation flux according to Reed (Reed,1977) formula
!
! Reed, R.K.,1977: On estimating insolation over the ocean, J.Phys.
! Oceanogr. 17, 854-871.

            call sol_rad(sol_net,cld,alon(i,j),alat(i,j)
     $                  ,iyr,imo,iday,ihour,imin)

! --- 1. Divide net solar radiation flux by rho*Cpw and reverse sign
            swrad(i,j) = -sol_net/rho_cpw
!

!---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
! Calculate the wind stress components (TAUX, TAUY)
!---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
!
! --- Calculate  the Drag Coefficient Cd accoriding to
!                         Hellerman & Rosenstein (1983)
!
            cd1 = cd(sp,deltemp)

! --- Calculate  the wind stresses in MKS ( newton/m*m )
! --- taux= rhom*Cd*|V|u     tauy= rhom*Cd*|V|v
!
            TauX = rhom*cd1*sp*unow
            TauY = rhom*cd1*sp*vnow

! --- Reverse Sign and divide by sea water density
            wusurf(i,j) = -taux/rho
            wvsurf(i,j) = -tauy/rho

          end if

!---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
! Multiply all fluxes by model mask
!---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
!
          wusurf(i,j) = wusurf(i,j)*fsm(i,j)
          wvsurf(i,j) = wvsurf(i,j)*fsm(i,j)
          wtsurf(i,j) = wtsurf(i,j)*fsm(i,j)
          pme(i,j)    =    pme(i,j)*fsm(i,j)

        end do
      end do
!
!
! --------------------------------------------------
      return

      end
!
      function HEATLAT(t)
!
! --- calculates the Latent Heat of Vaporization ( J/kg ) as function of
! --- the temperature ( Celsius degrees )
! --- ( from A. Gill  pag. 607 )
!
! --- Constant Latent Heat of Vaporization
!     L = 2.501e+6  (MKS)
!
        implicit none

        include 'realkind'

        real(kind=rk) heatlat
        real(kind=rk), intent(in)  :: t

        heatlat = 2.5008d+6 -2.3d+3*t

        return
      end
!
      function cd(sp,delt)
!
! --- calculates the Drag Coefficient as a function of the abs. value of the
! --- wind velocity
! --- ( Hellermann and  Rosenstein, 1983. JPO, 13,1093-1104)
        implicit none

        include 'realkind'

        real(kind=rk)                  cd
        real(kind=rk), intent(in)   :: sp, delt
        real(kind=rk), dimension(6) :: a
!
        data a / 0.934d-3,0.788d-4,0.868d-4
     $         ,-0.616d-6,-.120d-5,-.214d-5/
!
        cd = a(1) + a(2)*sp + a(3)*delt + a(4)*sp*sp
     $       + a(5)*delt*delt  + a(6)*sp*delt
!
        return
      end
!
!==============================================================================
!
      function ESK(t)
!
! --- compute the saturation water vapor pressure from
! --- temperature in kelvin  (Lowe,1977; JAM,16,100,1977)
        implicit none

        include 'realkind'

        real(kind=rk) esk
        real(kind=rk), intent(in) :: t

        real(kind=rk), dimension(7) :: a
!
        data  a /6984.505294,-188.9039310,2.133357675,-1.288580973d-2,
     $           4.393587233d-5,-8.023923082d-8,6.136820929d-11  /
!
        esk = a(1) +t*(a(2)+t*(a(3)+t*(a(4)+t*(a(5)+t*(a(6)+t*a(7))))))
!
        return
      end
!
      subroutine sol_rad( sola,cld,alonP,alatP,iyr,imt,idy,ihr,ime )
!--------------------------------------------------------------------

        include 'realkind'

        real(kind=rk), intent(out) :: sola
        real(kind=rk), intent(in)  :: cld, alonp, alatp
        integer      , intent(in)  :: iyr, imt, idy, ihr, ime
        real(kind=rk) qsw, blat, blon, acl

        real(kind=rk), external :: qshort

        data PI/3.141593/,RAD/.01745329/
!
!
! Panos scale solar radiation
        scal=1.0
!
! --- Rosati and Miyakoda :
!
        blon = alonp*pi/180.
        blat = alatp*pi/180.
        acl = cld
        qsw = qshort(iyr,imt,idy,ihr,ime,blat,blon,acl)
        sola = qsw
!
        return
      end

      function qshort(iyr,imt,idy,ihr,ime,alat,alon,acl)

        implicit none
!
        include 'realkind'

        real(kind=rk), intent(in)  :: alat, alon, acl
        integer      , intent(in)  :: iyr, imt, idy, ihr, ime

        real(kind=rk) degrad, eclips, raddeg, pi

        parameter(pi=3.1415927,degrad=pi/180.,raddeg=180./pi,
     $            eclips=23.439*degrad)
!
        dimension alpham(12),alb1(20),za(20),dza(19)

        real(kind=rk) qshort

        integer imt1, iyr1, intT1, intT2, jab
        real(kind=rk)
     &   albedo, alb1, alpha, alpham, aozone
     & , bb1, bb2
     & , capC, capG, capL, cosZen, DEC, DTOR, dza, dZen
     & , epsiln, g360, gha, gha360
     & , solar, SolAlt, SunBet, SunDec
     & , tau, ThSun, TRM111, TRM112, TRM11, UT, SHA, SMLT
     & , qatten, qdiff, qdir, qtot, qzer
     & , xl360, XLCT, yrdays, za, zen
!
! ---   alat,alon - (lat, lon)  in radians !!
!
        data solar/1350./
        data tau /0.7/
        data aozone /0.09/
        data yrdays /365./
        data alb1/.719, .656, .603, .480, .385, .300, .250, .193, .164
     $  ,.131 , .103, .084, .071, .061, .054, .039, .036, .032, .031
     $  ,.030 /
!
        data za/ 90., 88., 86., 84., 82., 80., 78., 76., 74., 70., 66.
     $  ,62., 58., 54., 50., 40., 30., 20., 10., 0.0 /
!
        data dza/8*2.0, 6*4.0, 5*10.0/
!
! --- albedo monthly values from Payne (1972) as means of the values
! --- at 40N and 30N for the Atlantic Ocean ( hence the same latitudinal
! --- band of the Mediterranean Sea ) :
!
        data alpham /0.09,0.08,0.06,0.06,0.06,0.06,0.06,0.06,
     $               0.06,0.07,0.09,0.10/
!
!--------------------- calculations start -----------------------------
!
! --- sun hour angle :
!
        DTOR = DEGRAD
        XLCT = ( 1.*ihr ) + ( 1.*ime / 60. )
        UT   = XLCT

        if ( imt > 2 ) then
          iyr1 = iyr
          imt1 = imt-3
        else
          iyr1 = iyr-1
          imt1 = imt+9
        end if

        intT1 = int(  30.6 * imt1      + 0.5 )
        intT2 = int( 365.25*(iyr1-1976)      )
        SMLT  = ( (UT/24.) + idy + intT1 + intT2 - 8707.5) / 36525.
        epsiln=  23.4393 -     0.013*SMLT
        capG  = 357.528  + 35999.050*SMLT

        if ( capG > 360. ) then
          g360 = capG - int( capG/360. )*360.
        else
          g360 = capG
        end if

        capC = 1.915*sin( g360*DTOR ) + .020*sin( 2.*g360*DTOR )
        capL = 280.46 + 36000.770*SMLT + capC

        if ( capL > 360. ) then
          xl360 = capL - int( capL/360. ) *360.
        else
          xl360 = capL
        end if

        alpha = xl360 -
     &     2.466*sin( 2.*xl360*DTOR ) + .053*sin( 4.*xl360*DTOR )
        gha = 15.*UT - 180.- capC + xl360 - alpha

        if ( gha > 360. ) then
          gha360 = gha - int( gha/360. )*360.
        else
          gha360 = gha
        end if

        DEC = atan( tan( epsiln*DTOR ) * sin( alpha*DTOR ) )/DTOR

!     Calculate Solar Hour Angle
        ThSun = ( GHA360 + alon*RADDEG )*degrad
        SHA = GHA360 + ( alon*RADDEG )


! --- sun declination :
        SUNDEC = DEC * DEGRAD

        TRM111 = sin( alat ) * sin( DEC*DTOR )
        TRM112 =-cos( alat ) * cos( DEC*DTOR )
        TRM11  = TRM111 - TRM112

! --- solar noon altitude in degrees :
        SolAlt = asin( TRM11 ) / DTOR
        SunBet = SolAlt

!
! --- cosine of the solar zenith angle :
!
        coszen =sin(alat)*sin(sundec)+cos(alat)*cos(sundec)*cos(ThSun)
!
        if ( coszen <= 0.+1.d-4 ) then
          coszen = 0.0
          qatten = 0.0
        else
          qatten = tau**(1./coszen)
        end if
        qzer  = coszen * solar
        qdir  = qzer * qatten
        qdiff = ((1.-aozone)*qzer - qdir) * 0.5
        qtot  =  qdir + qdiff
!
! --- ( radiation as from Reed(1977), Simpson and Paulson(1979) )
!
!
!-----------------------------------------------------------------------
! --- calculates the albedo as a function of the solar zenith angle :
! --- ( after Payne jas 1972 )
!-----------------------------------------------------------------------
!
! --- solar zenith angle in degrees :
!
        zen=raddeg*acos(coszen)
!
        if ( zen >= 74. ) then
          jab = int( .5 *(90.-zen) +  1. )
        elseif ( zen >= 50. ) then
          jab = int( .23*(74.-zen) +  9. )
        else
          jab = int( .10*(50.-zen) + 15. )
        end if
!
        dzen = ( za(jab) - zen ) / dza(jab)
!
        albedo=alb1(jab)+dzen*(alb1(jab+1)-alb1(jab))
!
! --- calculates SHORT WAVE FLUX ( watt/m*m )
! --- ( Rosati,Miyakoda 1988 ; eq. 3.8 )
!
!
        bb1=0.62      ! Reed
        bb2=0.636     ! Isemer et al. (1989)
!
        qshort = qtot*(1. - bb1*acl + 0.0019*sunbet)*(1. - albedo)
        if (qshort.gt.qtot) qshort=qtot
!
        return
      end
