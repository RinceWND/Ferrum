!
!
      subroutine bulk(im,jm,tbias, fsm, tsurf, alon,alat
     &               ,iyr,imo,iday,ihour,imin
     &               ,wusurf,wvsurf,wtsurf,swrad,pme
     &               ,uair,vair,usurf,vsurf,rsea
     &               ,tair,rhum,rain,cloud,pres,lwrd,my_task)
!*******************************************************************************
! THIS SUBROUTINE PROVIDES SURFACE BOUNDARY CONDITIONS FOR MOMENTUM,
! HEAT AND SALT EQUATIONS SOLVED BY THE HYDRODYNAMIC MODEL IN CASES WHEN
! THE ATMOSPHERIC MODEL PROVIDES CLOUD COVER DATA INSTEAD OF THE NET SOLAR
! RADIATION FLUX AND THE DOWNWARD LONGWAVE RADIATION FLUX AS IT IS ASSUMED
! IN THE FIRST VERSION OF BULK CODE. THE NET SOLAR RADIATION IS CALCULATED
! ACCORDING TO THE REED FORMULA WHILE THE NET LONGWAVE RADIATION CAN BE
! CALCULATED ACCORDING TO BIGNAMI, MAY, HERZFELD OR BERLIAND FORMULA
! (SEE `LWRAD_FORMULA` FLAG BELOW)
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
! ACCORDING TO BIGNAMI, MAY, (pseudo)HERZFELD, and BERLIAND FORMULA.
! THIS IS DONE THROUGH THE `LWRAD_FORMULA` PARAMETER (0 - 3)
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

      integer*1 LWRAD_FORMULA,
     &          lwBERLIAND, lwBIGNAMI, lwHERZFELD, lwMAY
!-------------------------------------------------------------------------------
! Longwave radiation formula (defaulted to lwBERLIAND)
      parameter ( lwBIGNAMI  = 0,
     &            lwMAY      = 1,
     &            lwHERZFELD = 2,
     &            lwBERLIAND = 3 )
      logical, parameter :: use_coare = .true., calc_swr = .false.

      integer, intent(in) :: my_task

      real(kind=rk), dimension(im, jm) ::
     $ fsm(im,jm),pme(im,jm),swrad(im,jm),wusurf(im,jm),
     $ wvsurf(im,jm),wtsurf(im,jm),tsurf(im,jm),alon(im,jm),
     $ alat(im,jm),usurf(im,jm),vsurf(im,jm),rsea(im,jm)

      real(kind=rk), dimension(im, jm) ::
     $  uair(im,jm),vair(im,jm),tair(im,jm),rhum(im,jm),
     $  rain(im,jm),cloud(im,jm),pres(im,jm)

      real(kind=rk) unow, vnow, tnow, pnow, precip, cld
     &   , sst_model, QBW, lwrd

      real(kind=rk), external :: cd, heatlat, esk

      real(kind=rk) arho
     &            , cd1, ce2, ch2
     &            , deltemp
     &            , ea12, emiss, esatair, esatoce, evap
     &            , fe, fh
     &            , Qe, Qh, Qu
     &            , rhnow, rho, rhom, rnow
     &            , sigma, sol_net, sp, ss, sstk, stp
     &            , taux, tauy, tnowk
     &            , usrf, vsrf
     &            , wair, wflux, wsatair, wsatoce
!
!--------------------------------------------------------------------
!       coefficients ( in MKS )  :
!-----------------------------------------------------------------
!
! --- Sea water density

      data rho/1025./

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
      LWRAD_FORMULA = -1  ! Default to Berliand formula (same as 3)
!
      do j = 1,jm
        do i = 1,im
          if (fsm(i,j) /= 0.) then
            unow      = uair(i,j)
            vnow      = vair(i,j)
            tnow      = tair(i,j)+ckelv
            pnow      = pres(i,j)
            rnow      = rsea(i,j)*rho + 1000.
!            e         = rhum(i,j)*pnow / ( 0.378*rhum(i,j) + 0.622 )
!            es        = 6.112 *
!     &           exp( (( 17.67*tair(i,j) )/( tair(i,j) + 243.25 )) )
!            rhnow     = 100.*e/es ! rwnd: specific to relative humidity
!            if (rhnow>100.) rhnow=100.
!            if (rhnow<0.) rhnow=0.
            rhnow     = rhum(i,j)
            precip    = rain(i,j)/1000. ! rwnd: precipitation rate from kg/(m2*s) to m/s
            cld       = cloud(i,j)/100. ! rwnd: total cloud cover from % to tenths
            sst_model = tsurf(i,j)+tbias

            if (i<im) then
              usrf = .5*(usurf(i+1,j)+usurf(i,j))
            else
              usrf = .5*(usurf(i-1,j)+usurf(i,j))
            end if

            if (j<jm) then
              vsrf = .5*(vsurf(i,j+1)+vsurf(i,j))
            else
              vsrf = .5*(vsurf(i,j-1)+vsurf(i,j))
            end if

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

            select case ( LWRAD_FORMULA )

              case ( lwBIGNAMI )
                QBW = 0.98*sigma*sstk**4 - sigma*tnowk**4
     $                 *(0.653+0.00535*ea12)
     $                 *(1.+0.1762*cld*cld)

              case ( lwMAY )
                QBW = (1.-0.75*(cld**3.4))
     $             * (sigma*(tnowk**4.)*(0.4 -0.05*sqrt(ea12))
     $           + 4.*sigma*(tnowk**3.)*(sstk-tnowk))

              case ( lwHERZFELD )
                QBW = (sigma*0.96*(1-(0.92e-5*tnowk*tnowk))*tnowk**4 +
     $             4*sigma*0.96*(ckelv+tair(i,j)**3)*(sstk-tnowk)) *
     $              (1-.75*cos(alat(i,j)*3.14/180.)*cld) ! cos(phi) here is an improvised `beta` coefficient as a function of latitude from Herzfeld

              case default  ! lwBERLIAND - Berliand (1952)
                QBW = 0.97*sigma*
     &               (     tnowk**4 * (.39-.05*sqrt(ea12))
     &                *(1.-.6823*cld*cld)
     &                + 4.*tnowk**3 * (sstk-tnowk))

            end select

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! --- calculate the term : ( Ts - Ta )
            deltemp = sstk - tnowk

            if ( calc_swr ) then
! --- Calculate net solar radiation flux according to Reed (Reed,1977) formula
!
! Reed, R.K.,1977: On estimating insolation over the ocean, J.Phys.
! Oceanogr. 17, 854-871.

              call sol_rad(sol_net,cld,alon(i,j),alat(i,j)
     $                  ,iyr,imo,iday,ihour,imin)

            else
              sol_net = swrad(i,j)
            end if
! --- 1. Divide net solar radiation flux by rho*Cpw and reverse sign
            swrad(i,j) = -sol_net/rho_cpw
!

          if ( use_coare ) then

            call coare35vn( (sqrt((unow-usrf)**2+(vnow-vsrf)**2)),
     &                      10._rk, tair(i,j), 2._rk, rhnow, 2._rk,
     &                      pnow,sst_model,sol_net,lwrd,alat(i,j),
     &                      600._rk, (3.6e6*precip), 0._rk, 0._rk,
     &                      QH, QE, Evap )
!            Evap = Evap*rho/3.6e6
!            if (my_task==1.and.i==50.and.j==50) then
!!            print *, "3.5 EVAP: ", Evap
!            write(61, '(6(f12.7,x),f12.7)') QE,QH,QBW,sol_net,
!     &                           (qbw+qh+qe-sol_net),lwrd,3.6e6*precip
!            end if

!            call coare30((unow-usurf(i,j)),(vnow-vsurf(i,j)),
!     &                    10._rk, tair(i,j), 2._rk, rhnow, 2._rk,
!     &                    pnow,sst_model,rnow,cld,precip*1000.,sol_net,
!     &                   -QBW, QH, QE, Evap )

!            if (my_task==1.and.i==50.and.j==50) then
!!            print *, "3.5 EVAP: ", Evap
!            write(62, '(6(f12.7,x),f12.7)') QE,QH,QBW,sol_net,
!     &                           (qbw+qh+qe-sol_net),lwrd,precip*1000.
!            end if
!            if (i==50.and.j==50) then
!            print *, "3.0 EVAP: ", Evap
!            print *, "3.0 QE:   ", QE
!            print *, "3.0 QH:   ", QH
!            print *, "RAIN===   ", precip
!            end if

          else
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
          end if
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
            if (isnan(QU)) then
            print *, ":", i,j,qbw, qh, qe
            end if
!
! --- 1. Devide upward heat flux by rho*Cpw
!
            wtsurf(i,j) = QU/rho_cpw
!

!---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
! Calculate the wind stress components (TAUX, TAUY)
!---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
!
! --- Calculate  the Drag Coefficient Cd accoriding to
!                         Hellerman & Rosenstein (1983)
!
!            cd1 = cd(sp,deltemp)

! --- Calculate  the wind stresses in MKS ( newton/m*m )
! --- taux= rhom*Cd*|V|u     tauy= rhom*Cd*|V|v
!
!            TauX = rhom*cd1*sp*unow
!            TauY = rhom*cd1*sp*vnow

! --- Reverse Sign and divide by sea water density
!            wusurf(i,j) = -taux/rho
!            wvsurf(i,j) = -tauy/rho

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
      pure function HEATLAT(t)
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

        heatlat = 2.5008e+6_rk - 2.3e+3_rk * t

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

        real(kind=rk), external :: qshort, qshort_simple

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
        coszen = sin(alat)*sin(sundec)+cos(alat)*cos(sundec)*cos(ThSun)
!
        if ( coszen <= 1.d-4 ) then
          coszen = 0.
          qatten = 0.
        else
          qatten = tau**(1./coszen)
        end if
        qzer  = coszen * solar
        qdir  = qzer * qatten
        qdiff = ((1.-aozone)*qzer - qdir) * .5
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

      function qshort_simple(iyr,imt,idy,ihr,ime,alat,alon,acl)

        implicit none
!
        include 'realkind'

        real(kind=rk), intent(in)  :: alat, alon, acl
        integer      , intent(in)  :: iyr, imt, idy, ihr, ime

        real(kind=rk) degrad, eclips, raddeg, pi

        parameter(pi=3.1415927,degrad=pi/180.,raddeg=180./pi,
     $            eclips=23.439*degrad)
!
        dimension alpham(12)

        real(kind=rk) qshort_simple

        integer imt1, iyr1, intT1, intT2
        real(kind=rk)
     &   albedo, alpha, alpham
     & , capC, capG, capL, coseca, DEC, DTOR, delta
     & , epsiln, g360, gha, gha360, phi
     & , SC, sinalt
     & , theta, UT, SHA, SMLT
     & , xl360, XLCT
!
! ---   alat,alon - (lat, lon)  in radians !!
!
        data SC/1370./ ! Solar constant [W/m2]
        data delta/0.85/ ! Atmospheric transmission coefficient
        data albedo/0.06/
!
! --- albedo monthly values from Payne (1972) as means of the values
! --- at 40N and 30N for the Atlantic Ocean ( hence the same latitudinal
! --- band of the Mediterranean Sea ) :
!
        data alpham /0.09,0.08,0.06,0.06,0.06,0.06,0.06,0.06,
     $               0.06,0.07,0.09,0.10/
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
        phi = ( GHA360 + alon*RADDEG )*degrad
        SHA = GHA360 + ( alon*RADDEG )

! --- sun declination :
        theta = DEC * DEGRAD

        sinalt = cos(theta)*cos(alat)*cos(phi)+sin(theta)*sin(alat)
        coseca = 1. / sinalt
!
        if ( coseca < -40. ) coseca = -40.

        rewind( 40 )
        write( 40, * ) delta, coseca, sinalt
        qshort_simple = SC*( delta**coseca * (sinalt-1.) + .91 )
     &                    *( 1 - .71*acl )*( 1. - albedo )
        if ( qshort_simple > SC ) qshort_simple = SC
!
        return
      end


      subroutine coare30( uw, vw, blk_ZW, Tair, blk_ZT, Hair, blk_ZQ
     &                  , Pair, Tsea, Rsea, cloud, rain, srflx, LRad
     &                  , shflx, lhflx, evap )

        implicit none

        include 'realkind'

        real(kind=rk), intent(in) :: uw, vw, Pair, Tair, Hair, Tsea
     &                             , LRad, Rsea, cloud, rain
        real(kind=rk) srflx, lrflx, evap

        real(kind=rk) Bf, blk_beta, blk_Cpa, blk_Cpw, blk_dter
     &              , blk_Rgas, blk_tcw, blk_visw, blk_Zabl, blk_ZQ
     &              , blk_ZT, blk_ZW, CC, Cd, Cd10, cff, cff1, cff2
     &              , Ch, Ch10, charn, Clam, Cp, Ct, Ct10, Cwet, delQ
     &              , delQc, delT, delTc, delW, diffh, diffw, e_sat
     &              , EminusP, emmiss, eps, Fc, g, Hcool, Hl, Hlb, Hlv
     &              , Hlw, Hs, Hsb, Hscale, Hscale2, Hsr, L, L10, lambd
     &              , LHeat, lhflx, PairM, pi, Q, Qair, Qbouy, Qcool
     &              , Qpsi, Qsea, Qstar, r3, RH, rhoAir, rhoref, rhoSea
     &              , rhow, Ri, Ribcu, Rr, scff, SHeat, shflx, SRad
     &              , StefBo, stflx, TairC, TairK, Taur, Taux, Tauy
     &              , Tcff, Tpsi, TseaC, TseaK, Tstar, twopi_inv, u10
     &              , upvel, Uwind, vap_p, VisAir, vonKar, Vwind
     &              , wet_bulb, Wgus, Wmag, Wpsi, Wspeed, Wstar
     &              , Zetu, Zo10, ZoL, ZoQ, ZoT, ZoT10, ZoW
        real(kind=rk), external :: bulk_psiu, bulk_psit
        integer Iter

        integer, parameter :: IterMax = 3
        logical, parameter :: COOLSKIN = .false.

        blk_Beta =   1.2
        blk_cpa  =1004.67
        blk_cpw  =   4.082793e3
        blk_dter =    .3
        blk_Rgas = 287.1
        blk_tcw  =    .6
        blk_visw =   1.e-6
        blk_Zabl = 600.
        Cp       =3983.212682927
        emmiss   =    .97
        g        =   9.81
        rhow     =1000.
        rhoref   =1025.
        stefbo   =   5.67e-8
        vonKar   =    .4
        eps = 1.e-20
        pi  = 4.*atan(1.)
        r3  = 1./3.
!
!  Compute Atmosphere-ocean fluxes using a bulk flux parameterization.
!
        Hscale    = rhoref*Cp
        Hscale2   = 1./(rhoref*Cp)
        twopi_inv = .5/pi
!
!  Input bulk parameterization fields.
!
        Wmag  = sqrt( uw*uw + vw*vw )
        PairM = Pair
        TairC = Tair
        TairK = TairC + 273.16
        TseaC = Tsea
        TseaK = TseaC + 273.16
        rhoSea= Rsea !*rhoref + 1000.
        RH    = Hair/100. ! Convert to fraction
        SRad  = srflx*Hscale
        Tcff  = 2.1e-5*( TseaC+3.2 )**0.79
        Scff  =  .026
!
!  Initialize.
!
        delTc = 0.
        delQc = 0.
!        LHeat = lhflx*Hscale
!        SHeat = shflx*Hscale
        Taur  = 0.
        Taux  = 0.
        Tauy  = 0.

!-----------------------------------------------------------------------
!  Compute net longwave radiation (W/m2), LRad.
!-----------------------------------------------------------------------
!
!  Use Berliand (1952) formula to calculate net longwave radiation.
!  The equation for saturation vapor pressure is from Gill (Atmosphere-
!  Ocean Dynamics, pp 606). Here the coefficient in the cloud term
!  is assumed constant, but it is a function of latitude varying from
!  1.0 at poles to 0.5 at the Equator).
!
        cff   = ( .7859 + .03477*TairC )/(1. + .00412*TairC )
        e_sat = 10.0**cff   ! saturation vapor pressure (hPa or mbar)
        vap_p = e_sat * RH       ! water vapor pressure (hPa or mbar)

        cff2  = TairK*TairK*TairK
        cff1  = cff2*TairK
!        LRad  =-emmiss*StefBo*
!     &              (cff1*( .39 - .05*sqrt(vap_p))*
!     &                    (1.   - .6823*cloud*cloud)+
!     &               cff2*4.*(TseaK-TairK))
!-----------------------------------------------------------------------
!  Compute specific humidities (kg/kg).
!
!    note that Qair(i) is the saturation specific humidity at Tair
!                 Q(i) is the actual specific humidity
!              Qsea(i) is the saturation specific humidity at Tsea
!
!          Saturation vapor pressure in mb is first computed and then
!          converted to specific humidity in kg/kg
!
!          The saturation vapor pressure is computed from Teten formula
!          using the approach of Buck (1981):
!
!          Esat(mb) = (1.0007_r8+3.46E-6_r8*PairM(mb))*6.1121_r8*
!                  EXP(17.502_r8*TairC(C)/(240.97_r8+TairC(C)))
!
!          The ambient vapor is found from the definition of the
!          Relative humidity:
!
!          RH = W/Ws*100 ~ E/Esat*100   E = RH/100*Esat if RH is in %
!                                       E = RH*Esat     if RH fractional
!
!          The specific humidity is then found using the relationship:
!
!          Q = 0.622 E/(P + (0.622-1)e)
!
!          Q(kg/kg) = 0.62197_r8*(E(mb)/(PairM(mb)-0.378_r8*E(mb)))
!
!-----------------------------------------------------------------------
!
!  Compute air saturation vapor pressure (mb), using Teten formula.
!
        cff = (1.0007 + 3.46e-6*PairM)*6.1121*
     &        exp( 17.502*TairC / (240.97+TairC) )
!
!  Compute specific humidity at Saturation, Qair (kg/kg).
!
        Qair = .62197*( cff / (PairM-.378*cff) )
!
!  Compute specific humidity, Q (kg/kg).
!
        cff = cff*RH                            !Vapor pres (mb)
        Q   = .62197*( cff / (PairM-.378*cff) ) !Spec hum (kg/kg)
!
!  Compute water saturation vapor pressure (mb), using Teten formula.
!
        cff = ( 1.0007 + 3.46e-6*PairM )*6.1121*
     &        exp( 17.502*TseaC / (240.97+TseaC) )
!
!  Vapor Pressure reduced for salinity (Kraus & Businger, 1994, pp 42).
!
        cff = cff*0.98
!
!  Compute Qsea (kg/kg) from vapor pressure.
!
        Qsea = 0.62197*( cff / (PairM-0.378*cff) )

!-----------------------------------------------------------------------
!  Compute Monin-Obukhov similarity parameters for wind (Wstar),
!  heat (Tstar), and moisture (Qstar), Liu et al. (1979).
!-----------------------------------------------------------------------
!
!  Moist air density (kg/m3).
!
        rhoAir = PairM*100.0 / ( blk_Rgas*TairK*(1.+.61*Q) )
!
!  Kinematic viscosity of dry air (m2/s), Andreas (1989).
!
        VisAir = 1.326e-5*
     &              (1.+TairC*(6.542e-3+TairC*
     &               (8.301e-6 - 4.84e-9*TairC)))
!
!  Compute latent heat of vaporization (J/kg) at sea surface, Hlv.
!
        Hlv = ( 2.501 - .00237*TseaC )*1.0e+6
!
!  Assume that wind is measured relative to sea surface and include
!  gustiness.
!
        Wgus = .5
        delW = sqrt( Wmag*Wmag + Wgus*Wgus )
        delQ = Qsea-Q
        delT = TseaC-TairC
!
!  Neutral coefficients.
!
        ZoW = 0.0001
        u10 = delW * log( 10./ZoW )/log( blk_ZW/ZoW )
        Wstar = .035*u10
        Zo10  = .011*Wstar*Wstar/g + .11*VisAir/Wstar
        Cd10  = ( vonKar / log( 10./Zo10 ) )**2
        Ch10  = .00115
        Ct10  = Ch10/sqrt(Cd10)
        ZoT10 = 10. / exp( vonKar/Ct10 )
        Cd    = ( vonKar / log( blk_ZW/Zo10 ) )**2
!
!  Compute Richardson number.
!
        Ct = vonKar / log( blk_ZT/ZoT10 )  ! T transfer coefficient
        CC = vonKar*Ct/Cd
        delTc = 0.
        delTc = blk_dter
        Ribcu = -blk_ZW/(blk_Zabl*0.004*blk_beta**3)
        Ri = -g*blk_ZW*((delT-delTc)+
     &                          .61*TairK*delQ)/
     &          (TairK*delW*delW)
        if ( Ri < 0. ) then
          Zetu = CC*Ri / ( 1. + Ri/Ribcu )   ! Unstable
        else
          Zetu = CC*Ri / ( 1. + 3.*Ri/CC )   ! Stable
        end if
        L10 = blk_ZW/Zetu
!
!  First guesses for Monon-Obukhov similarity scales.
!
        Wstar = delW*vonKar / ( log(blk_ZW/Zo10) -
     &                          bulk_psiu( blk_ZW/L10, pi ) )
        Tstar = -(delT-delTc)*vonKar/
     &           ( log( blk_ZT/ZoT10 )-
     &             bulk_psit( blk_ZT/L10, pi ) )
        Qstar = -(delQ-delQc)*vonKar/
     &           ( log(blk_ZQ/ZoT10 )-
     &             bulk_psit( blk_ZQ/L10, pi ) )
!
!  Modify Charnock for high wind speeds. The 0.125 factor below is for
!  1.0/(18.0-10.0).
!
        if ( delW > 18. ) then
          charn = .018
        elseif ( (10. < delW ).and.(delW >= 18. ) ) then
          charn = .011 + .125*(.018-.011)*(delW-10.)
        else
          charn = .011
        end if
!
!  Iterate until convergence. It usually converges within 3 iterations.
!
        do Iter = 1,IterMax

          ZoW = charn*Wstar*Wstar/g + .11*VisAir/(Wstar+eps)
          Rr  = ZoW*Wstar/VisAir
!
!  Compute Monin-Obukhov stability parameter, Z/L.
!
          ZoQ = min( 1.15e-4, 5.5e-5/Rr**.6 )
          ZoT = ZoQ
          ZoL = vonKar*g*blk_ZW*
     &             (Tstar*(1.+.61*Q)+
     &                        .61*TairK*Qstar)/
     &             (TairK*Wstar*Wstar*
     &             (1.+.61*Q)+eps)
          L = blk_ZW/(ZoL+eps)
!
!  Evaluate stability functions at Z/L.
!
          Wpsi = bulk_psiu( ZoL, pi )
          Tpsi = bulk_psit( blk_ZT/L, pi )
          Qpsi = bulk_psit( blk_ZQ/L, pi )
          if ( COOLSKIN ) then
            Cwet = .622*Hlv*Qsea/
     &              ( blk_Rgas*TseaK*TseaK )
            delQc = Cwet*delTc
          end if
!
!  Compute wind scaling parameters, Wstar.
!
          Wstar = max( eps, delW*vonKar/
     &                ( log( blk_ZW/ZoW ) - Wpsi ) )
          Tstar = -(delT-delTc)*vonKar/
     &                ( log( blk_ZT/ZoT ) - Tpsi )
          Qstar = -(delQ-delQc)*vonKar/
     &                ( log( blk_ZQ/ZoQ ) - Qpsi )
!
!  Compute gustiness in wind speed.
!
          Bf = -g/TairK*
     &         Wstar*( Tstar + .61*TairK*Qstar )
          if ( Bf > 0. ) then
            Wgus = blk_beta*( Bf*blk_Zabl )**r3
          else
            Wgus = .2
          end if
          delW = sqrt( Wmag*Wmag + Wgus*Wgus )

          if ( COOLSKIN ) then
!
!-----------------------------------------------------------------------
!  Cool Skin correction.
!-----------------------------------------------------------------------
!
!  Cool skin correction constants. Clam: part of Saunders constant
!  lambda; Cwet: slope of saturation vapor.
!
            Clam = 16.*g*blk_Cpw*( rhoSea*blk_visw )**3/
     &           ( blk_tcw*blk_tcw*rhoAir*rhoAir )
!
!  Set initial guesses for cool-skin layer thickness (Hcool).
!
            Hcool = 0.001
!
!  Backgound sensible and latent heat.
!
            Hsb = -rhoAir*blk_Cpa*Wstar*Tstar
            Hlb = -rhoAir*Hlv*Wstar*Qstar
!
!  Mean absoption in cool-skin layer.
!
            Fc = .065 + 11.*Hcool-
     &         (1.-exp(-Hcool*1250.))*6.6e-5/Hcool
!
!  Total cooling at the interface.
!
            Qcool = LRad+Hsb+Hlb-SRad*Fc
            Qbouy = Tcff*Qcool+Scff*Hlb*blk_Cpw/Hlv
!
!  Compute temperature and moisture change.
!
            if ( ( Qcool > 0. ).and.( Qbouy > 0. ) ) then
              lambd = 6./(1.+
     &               (Clam*Qbouy/(Wstar+eps)**4)**.75)**r3
              Hcool = lambd*blk_visw/(sqrt(rhoAir/rhoSea)*
     &                              Wstar+eps)
              delTc = Qcool*Hcool/blk_tcw
            else
              delTc = 0.
            end if
            delQc = Cwet*delTc
          end if

        end do
!-----------------------------------------------------------------------
!  Compute Atmosphere/Ocean fluxes.
!-----------------------------------------------------------------------
!
!
!  Compute transfer coefficients for momentum (Cd).
!
        Wspeed = sqrt( Wmag*Wmag + Wgus*Wgus )
        Cd = Wstar*Wstar / ( Wspeed*Wspeed+eps)
        Ch = Wstar*Tstar / (-Wspeed*delT+.0098*blk_ZT)
!
!  Compute turbulent sensible heat flux (W/m2), Hs.
!
        Hs = -blk_Cpa*rhoAir*Wstar*Tstar
!
!  Compute sensible heat flux (W/m2) due to rainfall (kg/m2/s), Hsr.
!
        diffw = 2.11e-5*(TairK/273.16)**1.94
        diffh =  .02411*(1.+TairC*
     &                    (3.309e-3 - 1.44e-6*TairC))/
     &          (rhoAir*blk_Cpa)
        cff = Qair*Hlv / (blk_Rgas*TairK*TairK)
        wet_bulb = 1. / (1.+.622*(cff*Hlv*diffw)/
     &                           (blk_Cpa*diffh))
        Hsr = rain*wet_bulb*blk_Cpw*
     &        ((TseaC-TairC)+(Qsea-Q)*Hlv/blk_Cpa)
        SHeat = Hs+Hsr
!
!  Compute turbulent latent heat flux (W/m2), Hl.
!
        Hl = -Hlv*rhoAir*Wstar*Qstar
!
!  Compute Webb correction (Webb effect) to latent heat flux, Hlw.
!
        upvel = -1.61*Wstar*Qstar-
     &          (1.+1.61*Q)*Wstar*Tstar/TairK
        Hlw = rhoAir*Hlv*upvel*Q
        LHeat = Hl+Hlw
!
!  Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
!
        Taur = .85*rain*Wmag
!
!  Compute wind stress components (N/m2), Tau.
!
        cff  = rhoAir*Cd*Wspeed
        Taux = (cff*Uwind + Taur*sign(1._rk,Uwind))
        Tauy = (cff*Vwind + Taur*sign(1._rk,Vwind))

!=======================================================================
!  Compute surface net heat flux and surface wind stress.
!=======================================================================
!
!  Compute kinematic, surface, net heat flux (degC m/s).  Notice that
!  the signs of latent and sensible fluxes are reversed because fluxes
!  calculated from the bulk formulations above are positive out of the
!  ocean.
!
!  For EMINUSP option,  EVAP = LHeat (W/m2) / Hlv (J/kg) = kg/m2/s
!                       PREC = rain = kg/m2/s
!
!  To convert these rates to m/s divide by freshwater density, rhow.
!
!  Note that when the air is undersaturated in water vapor (Q < Qsea)
!  the model will evaporate and LHeat > 0:
!
!                   LHeat positive out of the ocean
!                    evap positive out of the ocean
!
!  Note that if evaporating, the salt flux is positive
!        and if     raining, the salt flux is negative
!
!  Note that fresh water flux is positive out of the ocean and the
!  salt flux (stflx(isalt)) is positive into the ocean. It is converted
!  to (psu m/s) for stflx(isalt) in "set_vbc.F" or "ice_mk.h". The E-P
!  value is saved in variable EminusP for I/O purposes.
!
        cff = 1./rhow
        lrflx = LRad*Hscale2
        lhflx = LHeat !*Hscale2 ! do not invert latent...
        shflx = SHeat !*Hscale2 ! ...and sensible heat fluxes and convert them to [var*m/s] as well
        stflx =(srflx+lrflx+lhflx+shflx)
        evap  = LHeat/Hlv
        stflx = cff*(evap-rain)
        EminusP = stflx

      end ! subroutine
!----------------------------------------------------------------------
      subroutine coare35vn( u, zu, t, zt, rh, zq, P, ts, Rs, Rl, lat
     &                     ,zi, rain, cp, sigH, hsb, hlb, Evap )
!----------------------------------------------------------------------
! Input:
!
!     u = relative wind speed (m/s) at height zu(m)
!     t = bulk air temperature (degC) at height zt(m)
!    rh = relative humidity (%) at height zq(m)
!     P = surface air pressure (mb) (default = 1015)
!    ts = water temperature (degC) see jcool below
!    Rs = downward shortwave radiation (W/m^2) (default = 150)
!    Rl = downward longwave radiation (W/m^2) (default = 370)
!   lat = latitude (default = +45 N)
!    zi = PBL height (m) (default = 600m)
!  rain = rain rate (mm/hr)
!    cp = phase speed of dominant waves (m/s)
!  sigH =  significant wave height (m)
!
!
! Reference:
!
!  Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
!  Bulk parameterization of air sea fluxes: updates and verification for the
!  COARE algorithm, J. Climate, 16, 571-590.
!
        implicit none

        include 'realkind'

        integer*1, parameter :: jcool = 0

        real(kind=rk), external :: grv, psit_26, psiu_26
     &                           , psiu_40, qsat26sea, RHcalc

        real(kind=rk), intent(in) :: u,zu, t,zt, rh,zq, P, ts
     &                              ,Rs,Rl, lat, zi, rain, cp,sigH
        real(kind=rk)
     &                A, a1, a2, ad, Al, alfac, alq
     &              , B, Bd, be, Bf, Beta, bigc
     &              , CC, Cd, Cd10, cdhf, Cdn_10, Ce, Cen_10, Ch, Ch10
     &              , charn, charnC, charnS, charnW, Chn_10, cpa, cpv
     &              , cpw, cqhf, Ct, Ct10, cthf
     &              , dels, dq, dqer, dqs_dt, dt, dter, dtmp, dwat
     &              , Evap, fdg, gf, grav
     &              , hbb, hlb, hlwebb, hsb, hsbb
     &              , L, L10, lapse, Le
     &              , psi, psi10, psi10T, psirf, psirfQ, psirfT, psiT
     &              , Pv, Q, Q10, Q10N, Q10N2, qcol, QN, QN2, qout
     &              , Qrf, QrfN, QrfN2, Qs, qsr
     &              , RF, Rgas, RH10, rhoa, RHrf, rhodry, rhow
     &              , Ribcu, Ribu, Rnl, Rns, rr, S, S10, SSQ, SST
     &              , T10, T10N, T10N2, ta, tau, tcw, tdk, tkt, TN, TN2
     &              , Trf, TrfN, TrfN2, tsr, tssr, tvsr
     &              , u10, u10N, U10N2, ug, umax, UN, UN2, Urf, UrfN
     &              , UrfN2, usr, ut, visa, visw, von, wbar, wetc
     &              , xlamx, zet, zetu, zo, zo10, zoq, zoS, zot, zot10
     &              , zref, zrf_q, zrf_t, zrf_u
        integer*1 nits, i

        logical waveage, seastate

        waveage = .false.
        seastate = .false.

! convert rh to specific humidity (the functions below return g/kg)
        Qs = qsat26sea(ts,P)/1000.          ! surface water specific humidity [kg/kg]
        call qsat26air(t,P,rh, Q, Pv)       ! specific humidity of air [kg/kg]
        Q = Q/1000.

!-----------  set constants ----------------------------------------------
        zref = 10.
        Beta =  1.2                         ! Given as 1.25 in Fairall et al.(1996)
        von  =   .4                         ! Von Karman's "number"
        fdg  = 1.                           ! Turbulent Prandtl number
        tdk  = 273.16
        grav = grv(lat)

!-----------  air constants ----------------------------------------------
        Rgas   = 287.1                      ! Gas const. dry air [J/kg/K]
        Le     = ( 2.501 - .00237*ts )*1.e6 ! Latent Heat of vaporization [J/kg]
        cpa    = 1004.67                    ! Specific heat of dry air [J/kg/K] (Businger 1982)
        cpv    = cpa*( 1. + 0.84*Q )        ! Moist air - currently not used (Businger 1982)
        rhoa   =  P    *100./( Rgas*(t+tdk) * (1.+0.61*Q) ) ! Moist air density [kg/m^3]
        rhodry = (P-Pv)*100./( Rgas*(t+tdk)              )
        visa   = 1.326e-5*( 1.+t*(6.542e-3+t*(8.301e-6-4.84e-9*t)) ) !  Kinematic viscosity of dry air [m2/s], Andreas (1989).

!-----------  cool skin constants  ---------------------------------------
        Al   = 2.1e-5*( ts+3.2 )**0.79      ! Water thermal expansion coef.
        be   =  .026                        ! Salinity expansion coef.
        cpw  = 4000.                        ! Specific heat of water [J/kg/K]
        rhow = 1025.                        ! Seawater density
        visw = 1.e-6                        ! Kinematic viscosity of water [m^2/s]
        tcw  =  .6                          ! Thermal conductivity of water [W/m/K]
        bigc = 16.*grav*cpw*(rhow*visw)**3 / (tcw**2 * rhoa**2)
        wetc =  .622*Le*Qs / ( Rgas*(ts+tdk)**2 ) !correction for dq;slope of sat. vap.

!-----------  net radiation fluxes ---------------------------------------
        Rns =  .945 * Rs                    ! albedo correction
        Rnl =  .97*( 5.67e-8*( ts-.3*jcool+tdk )**4 - Rl) ! initial value

!----------------  begin bulk loop --------------------------------------------

!-----------  first guess ------------------------------------------------
!        du = u-us                          ! u is already relative to surface speed
        dt = ts-t-.0098*zt
        dq = Qs-Q
        ta = t+tdk
        ug = 0.5
        dter  = 0.3
        ut    = sqrt( u**2 + ug**2 )
        u10   = ut * log(10./1.e-4) / log(zu/1.e-4)
        usr   = 0.035 * u10
        zo10  = 0.011*usr**2/grav + 0.11*visa/usr
        Cd10  = ( von / log(10./zo10) )**2
        Ch10  = 0.00115
        Ct10  = Ch10 / sqrt(Cd10)
        zot10 = 10. / exp(von/Ct10)
        Cd    = ( von / log(zu/zo10) )**2
        Ct    =   von / log(zt/zot10)
        CC    =   von * Ct/Cd
        Ribcu = -zu/(zi*.004*Beta**3)
        Ribu  = -grav*zu/ta*( (dt-dter*jcool)+.61*ta*dq )/ut**2
        if ( Ribu < 0. ) then
          zetu = CC*Ribu/( 1. + Ribu/Ribcu )
        else
          zetu = CC*Ribu*( 1. + 3.*Ribu/CC )
        end if
        L10 = zu/zetu                         ! Monin-Obukhov length
        gf  = ut/u
        usr = ut*von / ( log(zu/zo10) - psiu_40(zu/L10) )
        tsr = -(dt-     dter*jcool)*von*fdg
     &                             /(log(zt/zot10)-psit_26(zt/L10))
        qsr = -(dq-wetc*dter*jcool)*von*fdg
     &                             /(log(zq/zot10)-psit_26(zq/L10))
        tkt = 0.001

!----------------------------------------------------------
!  The following gives the new formulation for the
!  Charnock variable
!----------------------------------------------------------

        charnC = 0.011
        umax   = 19.
        a1     =   .0017
        a2     =-  .0050

        charnC = a1*u10 + a2
        if ( u10 > umax ) charnC = a1*umax+a2
        if ( charnC < 0.011 ) charnC = 0.011

        A = 0.114   ! wave-age dependent coefficients
        B = 0.622

        Ad= 0.091   ! Sea-state/wave-age dependent coefficients
        Bd= 2.0

        charnW = A*(usr/cp)**B
        zoS    = sigH*Ad*(usr/cp)**Bd
        charnS = zoS*grav/usr/usr

        charn = 0.011 !*ones(N,1);
        if ( ut > 18. ) then
          charn = .018
        elseif ( ut > 10. ) then
          charn = 0.011 + (ut-10.)/(18.-10.)*(0.018-0.011)
        end if
!        charnC = charn ! Fix?

        nits = 10   ! number of iterations

        if ( zetu > 50. ) nits = 1
!--------------  bulk loop --------------------------------------------------

        do i = 1, nits

!          zet = von*grav*zu/ta*(tsr+.61*ta*qsr/(1+.61*Q))/(usr**2)
          zet = von*grav*zu*(tsr*(1.+.61*Q))
          if (waveage) then
            if (seastate) then
              charn = charnS
            else
              charn = charnW
            end if
          else
            charn = charnC
          end if
          L  = zu/zet
          zo = charn*usr**2/grav + .11*visa/usr  ! surface roughness
!          if (zo<1.d-10) zo = 1.d-10
          rr = zo*usr/visa
          zoq= min(1.6e-4, 5.8e-5/rr**.72)       ! These thermal roughness lengths give Stanton and
          zot= zoq                               ! Dalton numbers that closely approximate COARE 3.0
          cdhf = von    /(log(zu/zo) - psiu_26(zu/L))
          cqhf = von*fdg/(log(zq/zoq)- psit_26(zq/L))
          cthf = von*fdg/(log(zt/zot)- psit_26(zt/L))
          if (isnan(cdhf)) then
            print *, "zo:   ", charn, usr, grav, visa
            print *, "cdhf: ", von, zu, zo, L
            print *, "logs: ", log(zu/zo), psiu_26(zu/L), ut
          end if
          usr  = ut*cdhf
          qsr  =-(dq-wetc*dter*jcool)*cqhf
          tsr  =-(dt-dter*jcool)*cthf
          tvsr = tsr+0.61*ta*qsr
          tssr = tsr+0.51*ta*qsr
          Bf   =-grav/ta*usr*tvsr
          ug   = 0.2
          if ( Bf > 0. ) ug = Beta*(Bf*zi)**.333
          ut = sqrt(u**2 + ug**2)
          gf = ut/u
          hsb=-rhoa*cpa*usr*tsr
          hlb=-rhoa*Le*usr*qsr
          qout = Rnl+hsb+hlb
          dels = Rns*(0.065+11.*tkt-6.6e-5/tkt*(1.-exp(-tkt/8.0e-4))) ! 1./0.0008 = 1250.
          qcol = qout-dels
          alq  = Al*qcol + be*hlb*cpw/Le
          if ( alq > 0. ) then
            xlamx = 6./(1.+(bigc*alq/usr**4)**0.75)**0.333
            tkt   = xlamx*visw/(sqrt(rhoa/rhow)*usr)
          else
            xlamx= 6.
            tkt = min(0.01, xlamx*visw/(sqrt(rhoa/rhow)*usr))
          end if
          dter = qcol*tkt/tcw
          dqer = wetc*dter
          Rnl  = 0.97*(5.67e-8*(ts-dter*jcool+tdk)**4-Rl) ! update dter
          u10N = usr/von/gf*log(10/zo)
          charnC = a1*u10N+a2
          if ( u10N > umax ) charnC = a1*umax+a2
          charnW = A*(usr/cp)**B
          zoS = sigH*Ad*(usr/cp)**Bd - .11*visa/usr
          charnS = zoS*grav/usr/usr
          if ( charnC < 0.011 ) charnC = 0.011

        end do

!----------------  compute fluxes  --------------------------------------------
        tau = rhoa*usr*usr/gf       ! wind stress [N/m^2]
        hsb =-rhoa*cpa*usr*tsr      ! sensible heat flux [W/m^2]
        hlb =-rhoa*Le*usr*qsr       ! latent heat flux [W/m^2]
        hbb =-rhoa*cpa*usr*tvsr     ! buoyancy flux
        hsbb=-rhoa*cpa*usr*tssr     ! sonic heat flux
        wbar=1.61*hlb/Le/(1.+1.61*Q)/rhoa + hsb/rhoa/cpa/ta
        hlwebb = hlb + rhoa*wbar*Q*Le
        Evap   = hlwebb/Le
!        hlb = hlb + hlwebb ! ?????? Webb correction to latent heat flux already in ef via zoq/rr function so return hlwebb

!-----  compute transfer coeffs relative to ut @ meas. ht  --------------------
        Cd = tau/rhoa/ut/max(.1,u)
        Ch =-usr*tsr/ut/(dt-dter*jcool)
        Ce =-usr*qsr/ut/(dq-dqer*jcool)

!---  compute 10-m neutral coeff relative to ut (output if needed) ------------
        Cdn_10 = 1000.*von**2     / log(10./zo)**2
        Chn_10 = 1000.*von**2*fdg / log(10./zo)/log(10./zot)
        Cen_10 = 1000.*von**2*fdg / log(10./zo)/log(10./zoq)

!---  compute 10-m neutral coeff relative to ut (output if needed) ------------
!  Find the stability functions
!---------------------------------
        zrf_u = 10.             ! User defined reference heights
        zrf_t = 10.
        zrf_q = 10.
        psi   = psiu_26(zu/L)
        psi10 = psiu_26(10./L)
        psirf = psiu_26(zrf_u/L)
        psiT  = psit_26(zt/L)
        psi10T= psit_26(10./L)
        psirfT= psit_26(zrf_t/L)
        psirfQ= psit_26(zrf_q/L)
        gf    = ut/u

!---------------------------------------------------------
!  Determine the wind speeds relative to ocean surface
!  Note that usr is the friction velocity that includes
!  gustiness usr = sqrt(Cd) S, which is equation (18) in
!  Fairall et al. (1996)
!---------------------------------------------------------
        S = ut
        S10 = S + usr/von*(log(10./zu)-psi10+psi)
        U10 = S10/gf
! or U10 = U + usr/von/gf*(log(10./zu)-psi10+psi)
        Urf  = U   + usr/von/gf*(log(zrf_u/zu)-psirf+psi)
        UN   = U   + psi*usr/von/gf
        U10N = U10 + psi10*usr/von/gf
        UrfN = Urf + psirf*usr/von/gf

        UN2 = usr/von/gf*log(zu/zo)
        U10N2 = usr/von/gf*log(10./zo)
        UrfN2 = usr/von/gf*log(zrf_u/zo)

!******** rain heat flux (save to use if desired) *****************************
        dwat = 2.11e-5*((t+tdk)/tdk)**1.94  ! water vapour diffusivity
        dtmp =(1. + 3.309e-3*t - 1.44e-6*t*t)*0.02411/(rhoa*cpa) ! heat diffusivity
        dqs_dt = Q*Le/(Rgas*(t+tdk)**2.)  ! Clausius-Clapeyron
        alfac= 1./(1.+0.622*(dqs_dt*Le*dwat)/(cpa*dtmp))  ! wet bulb factor
        RF = rain*alfac*cpw*
     &      ((ts-t-dter*jcool)+(Qs-Q-dqer*jcool)*Le/cpa)/3600. ! [W/m2]??? 1/3600 here is to convert rain in mm/h to kg/m2/s.
        hsb = hsb+RF

        lapse = grav/cpa
        SST   = ts-dter*jcool

        T10 = T + tsr/von*(log(10./zt)-psi10T+psiT) + lapse*(zt-10.)
        Trf = T + tsr/von*(log(zrf_t/zt)-psirfT+psiT) + lapse*(zt-zrf_t)
        TN  = T + psiT*tsr/von
        T10N = T10 + psi10T*tsr/von
        TrfN = Trf + psirfT*tsr/von

        TN2   = SST + tsr/von*log(zt/zot)-lapse*zt
        T10N2 = SST + tsr/von*log(10./zot)-lapse*10
        TrfN2 = SST + tsr/von*log(zrf_t/zot)-lapse*zrf_t

        dqer = wetc*dter*jcool
        SSQ  = Qs-dqer
        SSQ  = SSQ*1000.
        Q   = Q*1000.
        qsr = qsr*1000.
        Q10 = Q + qsr/von*(log(10./zq)-psi10T+psiT)
        Qrf = Q + qsr/von*(log(zrf_q/zq)-psirfQ+psiT)
        QN  = Q + psiT*qsr/von/sqrt(gf)
        Q10N = Q10 + psi10T*qsr/von
        QrfN = Qrf + psirfQ*qsr/von

        QN2 = SSQ + qsr/von*log(zq/zoq)
        Q10N2 = SSQ + qsr/von*log(10./zoq)
        QrfN2 = SSQ + qsr/von*log(zrf_q/zoq)
        RHrf  = RHcalc(Trf,P,Qrf/1000.)
        RH10  = RHcalc(T10,P,Q10/1000.)

      end !subroutine coare25vn

!----------------------------------------------------------------------
      function qsat26sea(T,P)
!----------------------------------------------------------------------
! computes surface saturation specific humidity [g/kg]
! given T [degC] and P [mb]
!----------------------------------------------------------------------
        implicit none

        include 'realkind'

        real(kind=rk), external   :: bucksat
        real(kind=rk), intent(in) :: T, P
        real(kind=rk) ex, es, qsat26sea

        ex = bucksat(T,P)
        es = .98 * ex ! reduction at sea surface
        qsat26sea = 622.*es/( P - 0.378*es )

      end !function qsat26sea
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      function bucksat(T,P)
!----------------------------------------------------------------------
! computes saturation vapor pressure [mb]
! given T [degC] and P [mb]
!----------------------------------------------------------------------
        implicit none

        include 'realkind'

        real(kind=rk), intent(in) :: T, P
        real(kind=rk) bucksat
        bucksat = 6.1121*exp(17.502*T/(T+240.97))*(1.0007+3.46e-6*P)

      end
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      subroutine qsat26air(T,P,rh,q,em)
!----------------------------------------------------------------------
! computes surface saturation specific humidity [g/kg]
! given T [degC] and P [mb]
!----------------------------------------------------------------------
        implicit none

        include 'realkind'

        real(kind=rk), external   :: bucksat
        real(kind=rk), intent(in) :: T, P, rh
        real(kind=rk), intent(out):: q, em
        real(kind=rk) es

        es = bucksat(T,P)
        em = .01 * rh * es
        q  = 622.*es/( P - 0.378*es )

      end !function qsat26sea
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      function grv(lat)
!----------------------------------------------------------------------
! computes g [m/sec^2] given lat in deg
!----------------------------------------------------------------------
        implicit none

        include 'realkind'

        real(kind=rk), intent(in) :: lat
        real(kind=rk)             :: grv

        real(kind=rk), parameter  ::
     &                  gamma = 9.7803267715
     &                , c1    =  .0052790414
     &                , c2    =  .0000232718
     &                , c3    =  .0000001262
     &                , c4    =  .0000000007
        real(kind=rk) phi, pi, x

        pi = 4.*atan(1.)

        phi = lat*pi/180.
        x   = sin(phi)**2
        grv = gamma*( 1+x*( c1+x*( c2+x*( c3+x*c4 ) ) ) )

      end ! function grv
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      function psiu_40(zet)
!----------------------------------------------------------------------
! computes velocity structure function
!----------------------------------------------------------------------
        implicit none

        include 'realkind'

        real(kind=rk), intent(in) :: zet
        real(kind=rk)             :: psiu_40

        real(kind=rk), parameter ::
     &                  a = 1.
     &                 ,b = 3./4.
     &                 ,c = 5.
     &                 ,d =  .35
        real(kind=rk) dzet, f, psik, psic, x

        dzet = min(50., 0.35*zet)             ! stable
        psiu_40 = -(a*zet+b*(zet-c/d)*exp(-dzet)+b*c/d)
        if ( zet < 0. ) then                  ! unstable
          x = ( 1.-18.*zet)**.25
          psik = 2. *log((1.+x     )/2.) + log((1.+x*x)/2.)
     &                             - 2.*atan(x)+2.*atan(1.)
          x = ( 1.-10.*zet)**0.3333
          psic = 1.5*log((1.+x+x**2)/3.) - sqrt(3.)*
     &                                    atan((1.+2.*x)/sqrt(3.))
     &                                   + 4.*atan(1.)/sqrt(3.)
          f = zet**2 / ( 1+zet**2 )
          psiu_40 = (1-f)*psik + f*psic
        end if

      end !function psiu_40
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      function psit_26(zet)
!----------------------------------------------------------------------
! computes temperature structure function
!----------------------------------------------------------------------
        implicit none

        include 'realkind'

        real(kind=rk), intent(in) :: zet
        real(kind=rk)             :: psit_26

        real(kind=rk) dzet, f, psik, psic, x

        dzet = min(50., 0.35*zet)               ! stable
        if ( zet < 0. ) then                  ! unstable
          x = (1.-15.  *zet)**0.5
          psik = 2. *log((1+x     )/2.)
          x = (1.-34.15*zet)**0.3333
          psic = 1.5*log((1+x+x**2)/3.) - sqrt(3.)*
     &                                   atan((1.+2.*x)/sqrt(3.))
     &                                  + 4.*atan(1.)/sqrt(3.)
          f = zet**2/(1.+zet**2)
          psit_26 = (1-f)*psik + f*psic
        else
          psit_26 = -( (1+0.6667*zet)**1.5
     &          + 0.6667*(zet-14.28)*exp(-dzet) + 8.525 )
        end if

      end
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      function psiu_26(zet)
!----------------------------------------------------------------------
! computes velocity structure function
!----------------------------------------------------------------------
        implicit none

        include 'realkind'

        real(kind=rk), intent(in) :: zet
        real(kind=rk)             :: psiu_26

        real(kind=rk), parameter ::
     &                  a =  .7
     &                 ,b = 3./4.
     &                 ,c = 5.
     &                 ,d =  .35
        real(kind=rk) dzet, f, psik, psic, x

        dzet = min(50., 0.35*zet)               ! stable
        if ( zet < 0. ) then                   ! unstable
          x = (1.-15.  *zet)**0.25
          psik = 2. *log((1.+x     )/2.) + log((1.+x*x)/2.)
     &                             - 2.*atan(x)+2.*atan(1.)
          x = (1.-10.15*zet)**0.3333
          psic = 1.5*log((1.+x+x**2)/3.) - sqrt(3.)*
     &                                    atan((1.+2.*x)/sqrt(3.))
     &                                   + 4.*atan(1.)/sqrt(3.)
          f = zet**2 / ( 1+zet**2 )
          psiu_26 = (1.-f)*psik + f*psic
        else
          psiu_26 =-(a*zet+b*(zet-c/d)*exp(-dzet)+b*c/d)
        end if

      end !function psiu_26
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      function RHcalc(T,P,Q)
!----------------------------------------------------------------------
! computes relative humidity given T,P, & Q
!----------------------------------------------------------------------
        implicit none

        include 'realkind'

        real(kind=rk), intent(in) :: T, P, Q
        real(kind=rk)             :: RHcalc

        real(kind=rk) es, em

        es = 6.1121*exp(17.502*T/(T+240.97))*(1.0007+3.46e-6*P)
        em = Q*P/(0.378*Q+0.622)
        RHcalc = 100.*em/es

      end

      function bulk_psiu(ZoL, pi)
!
!=======================================================================
!                                                                      !
!  This function evaluates the stability function for  wind speed      !
!  by matching Kansas  and free convection forms.  The convective      !
!  form follows Fairall et al. (1996) with profile constants from      !
!  Grachev et al. (2000) BLM.  The  stable  form is from Beljaars      !
!  and Holtslag (1991).                                                !
!                                                                      !
!=======================================================================
!
        implicit none

        include 'realkind'

        real(rk) bulk_psiu
!
!  Imported variable declarations.
!
        real(rk), intent(in) :: ZoL, pi
!
!  Local variable declarations.
!
        real(rk), parameter :: r3 = 1./3.

        real(rk) :: Fw, cff, psic, psik, x, y
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Unstable conditions.
!
        if ( ZoL < 0. ) then
          x = (1.-15.*ZoL)**.25
          psik = 2.*log(.5*(1.+x))+
     &              log(.5*(1.+x*x))-
     &           2.*atan(x)+.5*pi
!
!  For very unstable conditions, use free-convection (Fairall).
!
          cff = sqrt(3.)
          y = (1.-10.15*ZoL)**r3
          psic = 1.5*log(r3*(1.+y+y*y))-
     &       cff*atan((1.+2.*y)/cff)+pi/cff
!
!  Match Kansas and free-convection forms with weighting Fw.
!
          cff = ZoL*ZoL
          Fw  = cff/(1.+cff)
          bulk_psiu=(1.-Fw)*psik+Fw*psic
!
!  Stable conditions.
!
        else
          cff = min( 50., .35*ZoL )
          bulk_psiu = -((1.+ZoL)+.6667*(ZoL-14.28)/
     &            exp(cff)+8.525)
        end if
        return
      end !function bulk_psiu

      function bulk_psit (ZoL, pi)
!
!=======================================================================
!                                                                      !
!  This function evaluates the  stability function  for moisture and   !
!  heat by matching Kansas and free convection forms. The convective   !
!  form follows Fairall et al. (1996) with  profile  constants  from   !
!  Grachev et al. (2000) BLM.  The stable form is from  Beljaars and   !
!  and Holtslag (1991).                                                !
!
!=======================================================================
!                                                                      !
!
        implicit none

        include 'realkind'

        real(rk) bulk_psit
!
!  Imported variable declarations.
!
        real(rk), intent(in) :: ZoL, pi
!
!  Local variable declarations.
!
        real(rk), parameter :: r3 = 1./3.

        real(rk) :: Fw, cff, psic, psik, x, y
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Unstable conditions.
!
        if ( ZoL < 0. ) then
          x = (1.-15.*ZoL)**.5
          psik = 2.*log(.5*(1.+x))
!
!  For very unstable conditions, use free-convection (Fairall).
!
          cff = sqrt(3.)
          y = (1.-34.15*ZoL)**r3
          psic = 1.5*log(r3*(1.+y+y*y))-
     &       cff*atan((1.+2.*y)/cff)+pi/cff
!
!  Match Kansas and free-convection forms with weighting Fw.
!
          cff = ZoL*ZoL
          Fw  = cff/(1.+cff)
          bulk_psit = (1.-Fw)*psik+Fw*psic
!
!  Stable conditions.
!
        else
          cff = min( 50., .35*ZoL )
          bulk_psit = -((1.+2.*ZoL)**1.5+
     &            0.6667*(ZoL-14.28)/EXP(cff)+8.525)
        end if

        return
      end !function bulk_psit