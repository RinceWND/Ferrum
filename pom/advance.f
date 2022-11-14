! advance.f

! advance POM

!______________________________________________________________________
!
      subroutine advance
!----------------------------------------------------------------------
!  Advances model a step
!----------------------------------------------------------------------
! called by: pom               [pom.f]
!
! calls    : check_nan         [advance.f]
!            check_nan_2d      [advance.f]
!            check_velocity    [advance.f]
!            ice_advance       [seaice.f]
!            lateral_viscosity [advance.f]
!            mode_interaction  [advance.f]
!            mode_external     [advance.f]
!            mode_internal     [advance.f]
!            print_section     [advance.f]
!            store_mean        [advance.f]
!            store_surf_mean   [advance.f]
!            update_bc         [advance.f]
!            update_time       [globals.f90]
!______________________________________________________________________
!
        use config   , only: spinup, use_ice
        use grid     , only: update_vcoord
        use model_run, only: iext, isplit
     &                     , update_time
        use seaice

! advance POM 1 step in time
        implicit none


! get time
        call update_time

! set time dependent boundary conditions
        if ( .not.spinup ) call update_bc

! set lateral viscosity
        call lateral_viscosity

! form vertical averages of 3-D fields for use in external (2-D) mode
        call mode_interaction

! external (2-D) mode calculation
        do iext = 1, isplit
!          write(6, "('==',i3.0)") iext
          call mode_external
          call check_nan_2d  !fhx:tide:debug
!          if (use_ice) call ice_advance
        end do

! refresh vertical coordinates
        call update_vcoord

! internal (3-D) mode calculation
        call mode_internal

! print section
        call print_section

! check nan
        call check_nan

! store mean 2010/4/26
        call store_mean

! store SURF mean
        call store_surf_mean

! write output
!      call write_output( dtime )

! write restart
!      if(mod(iint,irestart).eq.0) call write_restart_pnetcdf

! check CFL condition
        call check_velocity


      end ! subroutine advance
!
!______________________________________________________________________
!
      subroutine get_time
!----------------------------------------------------------------------
!  Returns the model time
!----------------------------------------------------------------------
! called by: [NO CALLS]
!______________________________________________________________________
!
      use model_run, only: dti, iint, ramp, time, time0

      implicit none


      time=dti*float(iint)/86400.+time0
      ramp=1.
!      if(lramp) then
!        ramp=time/period
!        if(ramp.gt.1.e0) ramp=1.e0
!      else
!        ramp=1.e0
!      endif


      end ! subroutine get_time
!
!______________________________________________________________________
!
      subroutine update_bc
!----------------------------------------------------------------------
!  Sets time-dependent boundary conditions
!----------------------------------------------------------------------
! called by: advance [advance.f]
!
! calls    : step    [air]
!            step    [bry]
!            step    [clim]
!            step    [seaice]
!______________________________________________________________________
!
      use air      , only: air_step => step
      use bry      , only: bry_step => step
      use clim     , only: clm_step => step
      use seaice   , only: ice_step => step
      use river    , only: riv_step => step
      use model_run, only: dtime

      implicit none


      call clm_step( dtime )
      call ice_step( dtime )
      call air_step( dtime )
      call bry_step( dtime )
      call riv_step( dtime )


      end ! subroutine update_bc
!
!______________________________________________________________________
!
      subroutine lateral_viscosity
!----------------------------------------------------------------------
!  Sets the lateral viscosity
!----------------------------------------------------------------------
! called by: advance        [advance.f]
!
! calls    : advct          [solver.f]
!            exchange3d_mpi [parallel_mpi.f]
!            pgscheme       [advance.f]
!______________________________________________________________________
!
      use config     , only: aam_init, horcon, mode, n1d, npg, rk
     &                     , pressure_gradient
      use bry        , only: aamfrz, USE_SPONGE
      use glob_const , only: MODE_BAROTROPIC
      use glob_domain, only: im, imm1, jm, jmm1, kmm1 ,my_task
      use grid       , only: dx, dy
      use glob_ocean , only: a, aam, aamfac, c, ee, u, v   ,drhox,dt

      implicit none

      integer i,j,k


! if mode=2 then initial values of aam2d are used. If one wishes
! to use Smagorinsky lateral viscosity and diffusion for an
! external (2-D) mode calculation, then appropiate code can be
! adapted from that below and installed just before the end of the
! "if(mode.eq.2)" loop in subroutine advave

! calculate Smagorinsky lateral viscosity:
! ( hor visc = horcon*dx*dy*sqrt((du/dx)**2+(dv/dy)**2
!                                +.5*(du/dy+dv/dx)**2) )
      if ( mode /= MODE_BAROTROPIC ) then

        call advct(a,c,ee)

        call pressure_gradient

!lyo:scs1d:
        if ( n1d /= 0 ) then
          aam(:,:,:) = aam_init
        else
          do k = 1, kmm1
            do j = 2, jmm1
              do i = 2, imm1
                aam(i,j,k) = horcon*dx(i,j)*dy(i,j)*aamfac(i,j)       !fhx:incmix
     &                      *sqrt( ( (u(i+1,j,k)-u(i,j,k))/dx(i,j) )**2
     &                           + ( (v(i,j+1,k)-v(i,j,k))/dy(i,j) )**2
     &                           +.5_rk*( .25_rk*( u(i  ,j+1,k)
     &                                           + u(i+1,j+1,k)
     &                                           - u(i  ,j-1,k)
     &                                           - u(i+1,j-1,k) )
     &                                          /dy(i,j)
     &                                  + .25_rk*( v(i+1,j  ,k)
     &                                           + v(i+1,j+1,k)
     &                                           - v(i-1,j  ,k)
     &                                           - v(i-1,j+1,k) )
     &                                          /dx(i,j) )**2 )
                if ( USE_SPONGE ) then
                  aam(i,j,k) = aam(i,j,k)*( 1._rk + aamfrz(i,j) )
                end if
              end do
            end do
          end do

!          if (my_task==0) then
!            i = 106
!            j = 4
!            print *, ">>> EXT_LAT_VISC_AAM: ", my_task
!            print *, "aam:   ", aam(i,j,1)
!    !        stop -1
!          end if

        end if !lyo:scs1d:
!
! create sponge zones
!        do k=1,kbm1
!          do j=2,jmm1
!            do i=2,imm1
!              aam(i,j,k)=aam(i,j,k)+1000*exp(-(j_global(j)-2)*1.5)
!     $                    +1000*exp((j_global(j)-jm_global+1)*1.5)
!            end do
!          end do
!        end do

        call exchange3d_mpi(aam(:,:,1:kmm1),im,jm,kmm1)
!        if (my_task==0) then
!          i = 106
!          j = 4
!          print *, ">>> EXT_LAT_VISC_AAM_POST_XCHNG: ", my_task
!          print *, "aam:   ", aam(i,j,1)
!  !        stop -1
!        end if

      end if


      end ! subroutine lateral_viscosity
!
!______________________________________________________________________
!
      subroutine mode_interaction
!----------------------------------------------------------------------
!  Forms vertical averages of 3-D fields for use in external (2-D) mode
!----------------------------------------------------------------------
! called by: advance [advance.f]
!
! calls    : advave  [solver.f]
!______________________________________________________________________
!
      use config     , only: mode
      use glob_const , only: MODE_BAROTROPIC
      use glob_domain, only: im, jm, kmm1 ,my_task
      use grid       , only: dz, fsm, h
      use glob_ocean , only: aam  , aam2d, advua, advva, advx , advy
     &                     , adx2d, ady2d, d    , drhox, drhoy, drx2d
     &                     , dry2d, egf  , el   , et   , tps  , ua
     &                     , utf  , va   , vtf
      use model_run  , only: isp2i, ispi

      implicit none

      integer i,j,k


      if ( mode /= MODE_BAROTROPIC ) then

        adx2d = 0.
        ady2d = 0.
        drx2d = 0.
        dry2d = 0.
        aam2d = 0.

!        if (my_task==0) then
!          i = 106
!          j = 4
!          print *, ">>> MODE_INT_NOT_BAR: ", my_task
!          print *, "advct:   ", advx(i,j,1), advy(i,j,1)
!          print *, "drho:    ", drhox(i,j,1), drhoy(i,j,1)
!          print *, "h:       ", h(i,j), et(i,j)
!          print *, "aam(k=1):", aam(i,j,1)
!  !        stop -1
!        end if

        do k = 1, kmm1

          adx2d = adx2d +  advx(:,:,k)
          ady2d = ady2d +  advy(:,:,k)
          drx2d = drx2d + drhox(:,:,k)
          dry2d = dry2d + drhoy(:,:,k)
          aam2d = aam2d +   aam(:,:,k)*dz(:,:,k)*fsm(:,:,k)

        end do
        aam2d = aam2d/( h + et )

        call advave(tps)

!        if (my_task==0) then
!          i = 106
!          j = 4
!          print *, ">>> MODE_INT_NOT_BAR_POST_ADVAVE: ", my_task
!          print *, "advave:  ", advua(i,j), advva(i,j)
!          print *, "adv2d:   ", adx2d(i,j), ady2d(i,j)
!          print *, "aam2d:   ", aam2d(i,j)
!  !        stop -1
!        end if

        adx2d = adx2d - advua
        ady2d = ady2d - advva

      end if

      egf = el*ispi

      do j = 1, jm
        do i = 2, im
          utf(i,j) = ua(i,j)*( d(i,j) + d(i-1,j) )*isp2i
        end do
      end do
      do j = 2, jm
        do i = 1, im
          vtf(i,j) = va(i,j)*( d(i,j) + d(i,j-1) )*isp2i
        end do
      end do


      end ! subroutine mode_interaction
!
!______________________________________________________________________
!
      subroutine mode_external
!----------------------------------------------------------------------
!  Calculates the external (2-D) mode
!----------------------------------------------------------------------
! called by: advance        [advance.f]
!
! calls    : advave         [solver.f]
!            exchange2d_mpi [solver.f]
!            bc_vel_ext     [bry.f90]
!            bc_zeta        [bry.f90]
!______________________________________________________________________
!
      use module_time
      use air        , only: e_atmos, vfluxf, wusurf, wvsurf
      use bry        , only: apply_tide, bc_vel_ext, bc_zeta ! TODO: Move apply_tide to tide
      use config     , only: alpha, cbcmax, cbcmin, hc, ispadv, smoth
     &                     , use_tide, wadsmoth, z0b, zsh,  mode
      use glob_const , only: grav, Kappa
      use glob_domain, only: im, imm1, jm, jmm1, kmm1   , my_task
      use grid       , only: art, aru, arv, cor, dx, dy, dzb, fsm, h, kb
     &                     , zz
      use glob_ocean
      use tide       , only: tide_ua, tide_va, tide_advance => step
      use model_run  , only: dte, dte2, iext, isp2i, ispi, isplit,dtime
     &                     , iint,iext

      implicit none

      integer i,j


      do j = 2, jm
        do i = 2, im
          fluxua(i,j) = .25_rk*(  d(i,j) +  d(i-1,j) )
     &                        *( dy(i,j) + dy(i-1,j) )*ua(i,j)
          fluxva(i,j) = .25_rk*(  d(i,j) +  d(i,j-1) )
     &                        *( dx(i,j) + dx(i,j-1) )*va(i,j)
        end do
      end do

!      if (my_task==0) then
!        i = 106
!        j = 3
!        print *, ">>> EXT_FLX: ", iext, my_task
!        print *, "flux:    ", fluxua(i,j), fluxva(i,j)
!        print *, "column:  ", d(i,j), d(i-1,j), d(i,j-1)
!        print *, "current: ", ua(i,j), va(i,j)
!!        stop -1
!      end if

! NOTE addition of surface freshwater flux, w(i,j,1)=vflux, compared
! with pom98.f. See also modifications to subroutine vertvl
      do j = 2, jmm1
        do i = 2, imm1
          elf(i,j) = elb(i,j)
     &             - dte2*( fluxua(i+1,j  ) - fluxua(i,j)
     &                    + fluxva(i  ,j+1) - fluxva(i,j) )/art(i,j)
     &             - vfluxf(i,j) ! vfluxf should already be normalised ny area
                                 ! FIXME: vfluxf should also be multiplied by dte2?
!          elf(i,j) = max( elf(i,j), hc*0.9_rk-h(i,j) )
        end do
      end do

!      if (my_task==0) then
!        i = 106
!        j = 3
!        print *, ">>> EXT_EL_PRE_BRY: ", iext, my_task
!        print *, "elev:    ", elf(i,j), elb(i,j)
!        print *, "fluxua:  ", fluxua(i,j), fluxua(i+1,j)
!        print *, "fluxva:  ", fluxva(i,j), fluxva(i,j+1)
!        print *, "vflux:   ", vfluxf(i,j)
!        print *, "h:       ", h(i,j), hc
!!        stop -1
!      end if

      call bc_zeta ! bcond(1)

      call exchange2d_mpi(elf,im,jm)

!      if (my_task==0) then
!        i = 106
!        j = 3
!        print *, ">>> EXT_EL_POST_BRY_XCNG: ", iext, my_task
!        print *, "elev:    ", elf(i,j), elb(i,j)
!        print *, "fluxua:  ", fluxua(i,j), fluxua(i+1,j)
!        print *, "fluxva:  ", fluxva(i,j), fluxva(i,j+1)
!        print *, "vflux:   ", vfluxf(i,j)
!        print *, "h:       ", h(i,j), hc
!!        stop -1
!      end if

! Update WaD mask !TODO: Manage subsequent layers for GCS
!      df = h + elf*fsm(:,:,1)
!      wdm = int(fsm(:,:,1),1)
!      where ( df <= hc )
!        wdm = 0
!      end where

!      if ( wadsmoth > 0._rk ) then
!        do j = 2,jmm1
!        do i = 2,imm1
!          if ( wdm(i,j) > ( wdm(i-1,j) + wdm(i+1,j)                     &
!     &                    + wdm(i,j-1) + wdm(i,j+1) )                   &
!     &         .and. ( df(i,j) <= hc*(1._rk+wadsmoth) ) ) then
!            wdm(i,j) = 0
!          end if
!        end do
!        end do
!        call exchange2d_byte_mpi(wdm,im,jm)
!      end if

      if ( mod(iext,ispadv) == 0 ) call advave(tps)

      do j = 2, jmm1
        do i = 2, im
          uaf(i,j) = adx2d(i,j) + advua(i,j)
     &             - aru(i,j)*.25_rk
     &              *( cor(i  ,j)*d(i  ,j)*(va(i  ,j+1)+va(i  ,j))
     &               + cor(i-1,j)*d(i-1,j)*(va(i-1,j+1)+va(i-1,j)) )
     &             + .25_rk*grav*(dy(i,j)+dy(i-1,j))
     &                          *(d (i,j)+d (i-1,j))
     &              *( (1._rk - 2._rk*alpha)*(el (i,j)-el (i-1,j))
     &               +                alpha *(elb(i,j)-elb(i-1,j)
     &                                       +elf(i,j)-elf(i-1,j))
     &               + (e_atmos(i,j)-e_atmos(i-1,j)) )
     &             + drx2d(i,j) + aru(i,j)*(wusurf(i,j)-wubot(i,j))
        end do
      end do

!      if (my_task==0) then
!        i = 106
!        j = 3
!        print *, ">>> EXT_UAF_INTERMEDIATE: ", iext, my_task
!        print *, "curr:    ", uaf(i,j)
!        print *, "curr_va: ", va(i,j), va(i,j+1), va(i-1,j+1), va(i-1,j)
!        print *, "adx2d:   ", adx2d(i,j)
!        print *, "advua:   ", advua(i,j)
!        print *, "elb:     ", d(i,j), d(i-1,j)
!        print *, "elb:     ", elb(i,j), elb(i-1,j)
!        print *, "el:      ", el(i,j), el(i-1,j)
!        print *, "elf:     ", elf(i,j), elf(i-1,j)
!        print *, "drx2d:   ", drx2d(i,j)
!        print *, "atmos:   ", e_atmos(i,j), e_atmos(i-1,j)
!        print *, "stress:  ", wusurf(i,j), wubot(i,j)
!!        stop -1
!      end if

      do j = 2, jmm1
        do i = 2, im
          uaf(i,j) = ( ( h(i,j) + elb(i,j) + h(i-1,j) + elb(i-1,j) )
     &                *aru(i,j)*uab(i,j)
     &               - 4._rk*dte*uaf(i,j) )
     &              /( ( h(i,j) + elf(i,j) + h(i-1,j) + elf(i-1,j) )
     &                *aru(i,j) )
        end do
      end do

!      if (my_task==0) then
!        i = 106
!        j = 3
!        print *, ">>> EXT_UAF: ", iext, my_task
!        print *, "curr:    ", uaf(i,j), uab(i,j)
!        print *, "elb:     ", elb(i,j), elb(i-1,j)
!        print *, "elf:     ", elf(i,j), elf(i-1,j)
!!        stop -1
!      end if

      do j = 2, jm
        do i = 2, imm1
          vaf(i,j) = ady2d(i,j) + advva(i,j)
     &             + arv(i,j)*.25_rk
     &              *( cor(i,j  )*d(i,j  )*(ua(i+1,j  )+ua(i,j  ))
     &               + cor(i,j-1)*d(i,j-1)*(ua(i+1,j-1)+ua(i,j-1)) )
     &             + .25_rk*grav*(dx(i,j)+dx(i,j-1))
     &                          *(d (i,j)+d (i,j-1))
     &              *( (1._rk - 2._rk*alpha)*(el (i,j)-el (i,j-1))
     &               +                alpha *(elb(i,j)-elb(i,j-1)
     &                                       +elf(i,j)-elf(i,j-1))
     &               + (e_atmos(i,j)-e_atmos(i,j-1)) )
     &             + dry2d(i,j) + arv(i,j)*(wvsurf(i,j)-wvbot(i,j))
        end do
      end do

      do j = 2, jm
        do i = 2, imm1
          vaf(i,j) = ( ( h(i,j) + elb(i,j) + h(i,j-1) + elb(i,j-1) )
     &                *vab(i,j)*arv(i,j)
     &               - 4._rk*dte*vaf(i,j) )
     &              /( ( h(i,j) + elf(i,j) + h(i,j-1) + elf(i,j-1) )
     &                *arv(i,j) )
        end do
      end do

!      if (my_task==0) then
!        i = 106
!        j = 3
!        print *, ">>> EXT_VAF: ", iext, my_task
!        print *, "curr:    ", vaf(i,j), vab(i,j)
!        print *, "elb:     ", elb(i,j), elb(i,j-1)
!        print *, "elf:     ", elf(i,j), elf(i,j-1)
!!        stop -1
!      end if

      if ( use_tide ) call tide_advance( dtime ) ! update tide boundaries before applying boundary conditions

      call bc_vel_ext ! bcond(2)

      call exchange2d_mpi(uaf,im,jm)
      call exchange2d_mpi(vaf,im,jm)

!      if (my_task==0) then
!        i = 106
!        j = 3
!        print *, ">>> EXT_UAF_POST_BRY_XCHNG: ", iext, my_task
!        print *, "curr:    ", uaf(i,j), vaf(i,j)
!!        stop -1
!      end if

      if     ( iext == (isplit-2) ) then

        etf = .25_rk*smoth*elf

      elseif ( iext == (isplit-1) ) then

        etf = etf + .5_rk*(1._rk-.5_rk*smoth)*elf

      elseif ( iext ==  isplit    ) then

        etf = ( etf + .5_rk*elf )*fsm(:,:,1)

      end if

! Restrict fluxes on dry-to-wet directions
!      do j = 1, jm
!      do i = 2, im
!        if ( .5_rk*( df(i,j) + df(i-1,j) ) <= hc ) then
!          uaf(i,j) = 0._rk
!        else
!          if ( wdm(i-1,j)==0._rk .and. uaf(i,j)>0._rk ) uaf(i,j) = 0._rk
!          if ( wdm(i  ,j)==0._rk .and. uaf(i,j)<0._rk ) uaf(i,j) = 0._rk
!        end if
!      end do
!      end do
!      call exchange2d_mpi(uaf,im,jm)
!
!      do j = 2, jm
!      do i = 1, im
!        if ( .5_rk*( df(i,j) + df(i,j-1) ) <= hc ) then
!          vaf(i,j) = 0._rk
!        else
!          if ( wdm(i,j-1)==0._rk .and. vaf(i,j)>0._rk ) vaf(i,j) = 0._rk
!          if ( wdm(i,j  )==0._rk .and. vaf(i,j)<0._rk ) vaf(i,j) = 0._rk
!        end if
!      end do
!      end do
!      call exchange2d_mpi(vaf,im,jm)

!      if (my_task==0) then
!        i = 106
!        j = 3
!        print *, ">>> EXT_UAF_POST_WAD: ", iext, my_task
!        print *, "curr:    ", uaf(i,j), vaf(i,j)
!        print *, "mask:    ", wdm(i,j), wdm(i-1,j), wdm(i,j-1)
!        print *, "df:      ", df(i,j), df(i-1,j), df(i,j-1)
!!        stop -1
!      end if

! apply filter to remove time split
      ua = ua + .5_rk*smoth*( uab + uaf - 2._rk*ua )
      va = va + .5_rk*smoth*( vab + vaf - 2._rk*va )
      el = el + .5_rk*smoth*( elb + elf - 2._rk*el )

      elb = el
      el  = elf
      d   = h + el
      uab = ua
      ua  = uaf
      vab = va
      va  = vaf

!      if (my_task==0) then
!        i = 106
!        j = 3
!        print *, ">>> EXT_UAF_ADVANCE: ", iext, my_task
!        print *, "curr:    ", uaf(i,j), vaf(i,j)
!        print *, "elev:    ", elf(i,j)
!!        stop -1
!      end if

! DELETE: Debug
!      do j = 1, jm
!      do i = 1, im
!        if ( d(i,j) < 0. ) then
!          print *, "NEGATIVE WCD @ ", i, j
!          print *, wdm(i,j)
!          d(i,j) = 0.01_rk
!          wdm(i,j) = 0
!          print *, uaf(i,j), uaf(i+1,j)
!          print *, vaf(i,j), vaf(i+1,j)
!          !stop -1
!        end if
!      end do
!      end do

! update bottom friction
      call bottom_friction

      if ( iext /= isplit ) then

        egf = egf + el*ispi

        do j = 1, jm
          do i = 2, im
            utf(i,j) = utf(i,j) + ua(i,j)*( d(i,j) + d(i-1,j) )*isp2i
          end do
        end do
        do j = 2, jm
          do i = 1, im
            vtf(i,j) = vtf(i,j) + va(i,j)*( d(i,j) + d(i,j-1) )*isp2i
          end do
        end do

        call exchange2d_mpi(utf,im,jm) ! FIXME: Unnecessary?
        call exchange2d_mpi(vtf,im,jm)

      end if


      end ! subroutine mode_external
!______________________________________________________________________
!
      subroutine mode_internal
!----------------------------------------------------------------------
!  Calculates the internal (3-D) mode
!----------------------------------------------------------------------
! called by: advance        [advance.f]
!
! calls    : advq           [solver.f]
!            advt1          [solver.f]
!            advt2          [solver.f]
!            advu           [solver.f]
!            advv           [solver.f]
!            bc_turb        [bry.f90]
!            bc_vel_int     [bry.f90]
!            bc_vel_vert    [bry.f90]
!            dens           [solver.f]
!            exchange3d_mpi [parallel_mpi.f]
!            profq          [solver.f]
!            proft          [solver.f]
!            profu          [solver.f]
!            profv          [solver.f]
!______________________________________________________________________
!
      use air        , only: vfluxb, vfluxf, wssurf, wtsurf
      use bry        , only: bc_ts, bc_turb, bc_vel_int, bc_vel_vert
      use clim       , only: tclim, sclim, relax_to_clim
      use glob_const , only: MODE_BAROTROPIC, MODE_DIAGNOSTIC, rk, small
      use config     , only: hc, hhi, mode , nadv, nbcs, nbct
     &                     , do_restart, smoth, s_hi, s_lo, t_hi, t_lo
      use glob_domain
      use grid       , only: fsm, dum, dvm, dz, dzb, dzf, h
      use glob_ocean
      use model_run

      implicit none

      real(rk), dimension(im,jm) :: col
      real(rk)                      col_h
      integer                       i,j,k


      if ( mode /= MODE_BAROTROPIC ) then
      if ( ( iint>2 .and. .not.do_restart )
     &               .or.      do_restart ) then

! adjust u(z) and v(z) such that depth average of (u,v) = (ua,va)
        tps = 0._rk
        col = 0._rk

        do k = 1, kmm1
          do j = 1, jm
            do i = 2, im
              col_h = .5_rk*( dz(i,j,k) + dz(i-1,j,k) )*dum(i,j,k)
              col(i,j) = col(i,j) + col_h
              tps(i,j) = tps(i,j) + u(i,j,k)*col_h
            end do
          end do
        end do

        do k = 1, kmm1
          do j = 1, jm
            do i = 2, im
              if ( dum(i,j,1) > 0._rk ) then
                u(i,j,k) = u(i,j,k)
     &                   + ( .5_rk*( utb(i,j) + utf(i,j) ) - tps(i,j) )
     &                    /col(i,j)
              end if
              u(i,j,k) = u(i,j,k)*dum(i,j,k)
            end do
          end do
        end do

        tps = 0._rk
        col = 0._rk

        do k = 1, kmm1
          do j = 2, jm
            do i = 1, im
              col_h = .5_rk*( dz(i,j,k) + dz(i,j-1,k) )*dvm(i,j,k)
              col(i,j) = col(i,j) + col_h
              tps(i,j) = tps(i,j) + v(i,j,k)*col_h
            end do
          end do
        end do

        do k = 1, kmm1
          do j = 2, jm
            do i = 1, im
              if ( dvm(i,j,1) > 0._rk ) then
                v(i,j,k) = v(i,j,k)
     &                   + ( .5_rk*( vtb(i,j) + vtf(i,j) ) - tps(i,j) )
     &                    /col(i,j)
              end if
              v(i,j,k) = v(i,j,k)*dvm(i,j,k)
            end do
          end do
        end do

!eda: do drifter assimilation
!eda: this was the place where Lin et al. assimilate drifters
!eda: /home/xil/OIpsLag/gomc27_test87d_hcast2me_pseudo1.f
!     IF(MOD(IINT,16).EQ.0) THEN   ! 16 = 3.*3600./dti, every 3 hours
!--      IF(iint.eq.1 .or. MOD(IINT,16).EQ.0) THEN
!--      IF(MOD(IINT,IASSIM).EQ.0) THEN
!        call assimdrf_OIpsLag(time, itime1, mins, sec,
!    1     IM, JM, KB, u, v, Nx, Ny, beta, alon, alat, zz, D,
!    2     igs, ige, jgs, jge, ndrfmax,
!    3     ub, vb, dz, DrfDir)
!     endif

! calculate w from u, v, dt (h+et), etf and etb
        call vertvl(a,c)

        call bc_vel_vert ! bcond(5)

        call exchange3d_mpi(w,im,jm,km)

! set uf and vf to zero
        uf = 0._rk
        vf = 0._rk

! calculate q2f and q2lf using uf, vf, a and c as temporary variables
        call advq(q2b,q2,uf,a,c)
        call advq(q2lb,q2l,vf,a,c)
        call profq(a,c,tps,dtef)

! an attempt to prevent underflow (DEBUG)
!        where ( q2l  < .5_rk*small ) q2l  = .5_rk*small
!        where ( q2lb < .5_rk*small ) q2lb = .5_rk*small

        call bc_turb ! bcond(6)

        call exchange3d_mpi(uf(:,:,2:kmm1),im,jm,kmm2)
        call exchange3d_mpi(vf(:,:,2:kmm1),im,jm,kmm2)

        q2  = q2  + .5_rk*smoth*( q2b  + uf - 2._rk*q2  )
        q2l = q2l + .5_rk*smoth*( q2lb + vf - 2._rk*q2l )
        q2b = q2
        q2  = uf
        q2lb= q2l
        q2l = vf


! calculate tf and sf using uf, vf, a and c as temporary variables
!        if( mode /= 4 .and. ( iint > 2 .or. do_restart ) ) then
        if ( mode /= MODE_DIAGNOSTIC ) then

          if     ( nadv == 1 ) then

            call advt1(tb,t,tclim,uf,a,c,'T')
            call advt1(sb,s,sclim,vf,a,c,'S')

          elseif ( nadv == 2 ) then

            call advt2(tb,tclim,uf,a,c,'T')
            call advt2(sb,sclim,vf,a,c,'S')

          else

            error_status = 1
            print *, '(/''Error: invalid value for nadv'')'

          end if

          call proft(uf,wtsurf,tsurf,nbct,tps)
          call proft(vf,wssurf,ssurf,nbcs,tps)

          call relax_to_clim( uf, vf )

          call bc_ts ! bcond(4)

          if ( t_lo > -999. ) then
            where ( uf < t_lo ) uf = t_lo
          end if
          if ( t_hi <  999. ) then
            where ( uf > t_hi ) uf = t_hi
          end if

          if ( s_lo > -999. ) then
            where ( vf < s_lo ) vf = s_lo
          end if
          if ( s_hi <  999. ) then
            where ( vf > s_hi ) vf = s_hi
          end if

          call exchange3d_mpi(uf(:,:,1:kmm1),im,jm,kmm1)
          call exchange3d_mpi(vf(:,:,1:kmm1),im,jm,kmm1)

! Relax T/S for dry cells
          do k = 1, kmm1
            where ( wdm == 0 )
              uf(:,:,k) = wet_relx1*tb(:,:,k) + wet_relx2*tsurf(:,:)
     &                                         *fsm(:,:,1)
              vf(:,:,k) = wet_relx1*sb(:,:,k) + wet_relx2*ssurf(:,:)
     &                                         *fsm(:,:,1)
            end where
          end do

          call exchange3d_mpi(uf(:,:,1:kmm1),im,jm,kmm1)
          call exchange3d_mpi(vf(:,:,1:kmm1),im,jm,kmm1)

          ! TODO: Make sure layer thickness cannot go infinitesimal
          where ( dz > 0. )
            t = t + .5_rk*smoth*( dzb*tb + dzf*uf -2._rk*dz*t )/dz
            s = s + .5_rk*smoth*( dzb*sb + dzf*vf -2._rk*dz*s )/dz
          end where

          tb = t
          t  = uf
          sb = s
          s  = vf

          call dens(s,t,rho)

        end if

! calculate uf and vf
        call advu
        call advv
        call profu
        call profv

        call bc_vel_int ! bcond(3)

! Mask currents
!        do k = 1, kmm1
!          where ( wdm(2:im,:)*wdm(1:imm1,:) == 0 )
!     &          uf(2:im,:,k) = 0._rk
!          where ( wdm(:,2:jm)*wdm(:,1:jmm1) == 0 )
!     &          vf(:,2:jm,k) = 0._rk
!        end do

        call exchange3d_mpi(uf(:,:,1:kmm1),im,jm,kmm1)
        call exchange3d_mpi(vf(:,:,1:kmm1),im,jm,kmm1)

        tps = 0._rk

        do k = 1, kmm1
          do j = 1, jm
            do i = 2, im
              tps(i,j) = tps(i,j)
     &                 + ( uf(i,j,k) + ub(i,j,k) - 2._rk*u(i,j,k) )
     &                  *( dz(i,j,k) + dz(i-1,j,k) )*dum(i,j,k)
            end do
          end do
        end do

        tps(2:im,1:jm) = tps(2:im,1:jm)
     &                  /( h(2:im  ,1:jm) + et(2:im  ,1:jm)
     &                   + h(1:imm1,1:jm) + et(1:imm1,1:jm) )

        do k = 1, kmm1
          do j = 1, jm
            do i = 2, im
              u(i,j,k) = u(i,j,k)
     &                 + .5_rk*smoth*( uf(i,j,k) + ub(i,j,k)
     &                               - 2._rk*u(i,j,k) - tps(i,j) )
            end do
          end do
        end do

        tps = 0._rk

        do k = 1, kmm1
          do j = 2, jm
            do i = 1, im
              tps(i,j) = tps(i,j)
     &                 + ( vf(i,j,k) + vb(i,j,k) - 2._rk*v(i,j,k) )
     &                  *( dz(i,j,k) + dz(i,j-1,k) )*dvm(i,j,k)
            end do
          end do
        end do

        tps(1:im,2:jm) = tps(1:im,2:jm)
     &                  /( h(1:im,2:jm  ) + et(1:im,2:jm  )
     &                   + h(1:im,1:jmm1) + et(1:im,1:jmm1) )

        do k = 1, kmm1
          do j = 1, jm
            do i = 1, im
              v(i,j,k) = v(i,j,k)
     &                 + .5_rk*smoth*( vf(i,j,k) + vb(i,j,k)
     &                               - 2._rk*v(i,j,k) - tps(i,j) )
            end do
          end do
        end do

        ub = u
        u  = uf
        vb = v
        v  = vf

        call geopotential_vertical_velocity

      end if
      end if

      egb = egf
      etb = et
      et  = etf
      dt  = h + et
      utb = utf
      vtb = vtf
      vfluxb = vfluxf

      dzb = dz
      dz  = dzf


      end ! subroutine mode_internal
!
!______________________________________________________________________
!
      subroutine geopotential_vertical_velocity
!----------------------------------------------------------------------
!  Calculates real (geopotential) vertical velocity as `wr`
!----------------------------------------------------------------------
! called by: mode_internal  [advance.f]
!
! calls    : exchange3d_mpi [parallel_mpi.f]
!______________________________________________________________________
!
        use glob_const , only: rk
        use glob_domain, only: im, imm1, jm, jmm1, km, kmm1
     &                       , n_east, n_north, n_south, n_west
        use grid       , only: dx, dy, fsm, zz
        use glob_ocean , only: dt, et, etb, etf, tps, u, v, w, wr
        use model_run  , only: dti2

        implicit none

        integer  i, j, k
        real(rk) dxr, dxl, dyt, dyb


        wr = 0.

        do k = 1, kmm1

          tps = zz(:,:,k)*dt + et
 
          do j = 2, jmm1
            do i = 2, imm1
              dxr = 2._rk/(dx(i+1,j)+dx(i  ,j))
              dxl = 2._rk/(dx(i  ,j)+dx(i-1,j))
              dyt = 2._rk/(dy(i,j+1)+dy(i,j  ))
              dyb = 2._rk/(dy(i,j  )+dy(i,j-1))
              wr(i,j,k) = .5_rk*(w(i,j,k)+w(i,j,k+1))
     $                  + .5_rk*
     $                   ( u(i+1,j,k)*(tps(i+1,j)-tps(i  ,j))*dxr
     $                   + u(i  ,j,k)*(tps(i  ,j)-tps(i-1,j))*dxl
     $                   + v(i,j+1,k)*(tps(i,j+1)-tps(i,j  ))*dyt
     $                   + v(i,j  ,k)*(tps(i,j  )-tps(i,j-1))*dyb )
     $                  + (1._rk+zz(i,j,k))*(etf(i,j)-etb(i,j))/dti2
            end do
          end do

        end do

        call exchange3d_mpi(wr(:,:,1:kmm1),im,jm,kmm1)

        do k = 1, km
          do i = 1, im
            if(n_south.eq.-1) wr(i,1,k)=wr(i,2,k)
            if(n_north.eq.-1) wr(i,jm,k)=wr(i,jmm1,k)
          end do
        end do
        do k = 1, km
          do j = 1, jm
            if(n_west.eq.-1) wr(1,j,k)=wr(2,j,k)
            if(n_east.eq.-1) wr(im,j,k)=wr(imm1,j,k)
          end do
        end do

        wr = fsm*wr


      end ! subroutine geopotential_vertical_velocity
!
!______________________________________________________________________
!
      subroutine print_section
!----------------------------------------------------------------------
!  Prints output
!----------------------------------------------------------------------
! called by: advance      [advance.f]
!
! calls    : bcast0i_mpi  [parallel_mpi.f]
!            finalize_mpi [parallel_mpi.f]
!            sum0d_mpi    [parallel_mpi.f]
!            sum0i_mpi    [parallel_mpi.f]
!______________________________________________________________________
!
      use config     , only: mode, sbias, tbias
      use glob_const , only: rk, MODE_DIAGNOSTIC
      use glob_domain
      use grid       , only: art, dz, fsm
      use glob_ocean , only: et, etb, dt, sb, tb, u, ub, v, vb, w, wdm
      use glob_out   , only: iprint
      use model_run

      implicit none

      real(rk), dimension(im,jm) :: d_vol
      real(rk) area_tot, vol_tot, kin_ave
      real(rk) elev_ave, temp_ave, salt_ave
      real(rk) e_resid, u_resid
      integer i,j,k


      if ( mod(iint,iprint) == 0 ) then

! print time
        if ( is_master ) print '(/
     $    ''=========================================================''
     $    /''time ='',f9.4,'', iint ='',i8,'', iext ='',i8,
     $    '', iprint ='',i8)', time,iint,iext,iprint

! check for errors
        call   sum0i_mpi(error_status,master_task)
        call bcast0i_mpi(error_status,master_task)
        if ( error_status /= 0 ) then
          if ( is_master ) print *, 'POM terminated with error'
          call finalize_mpi
          stop
        end if

! local averages
        vol_tot  = 0._rk
        area_tot = 0._rk
        kin_ave  = 0._rk
        temp_ave = 0._rk
        salt_ave = 0._rk
        elev_ave = 0._rk
! diagnostic residuals
        e_resid = 0.
        u_resid = 0.

        do k=1,kmm1
          d_vol    = art*dt*dz(:,:,k)*wdm
          vol_tot  = vol_tot + sum(d_vol)
          kin_ave  = kin_ave + sum( (u(:,:,k)*u(:,:,k)
     &                              +v(:,:,k)*v(:,:,k))*d_vol )
          temp_ave = temp_ave + sum(tb(:,:,k)*d_vol)
          salt_ave = salt_ave + sum(sb(:,:,k)*d_vol)
        end do

        area_tot = sum( art*wdm )
        elev_ave = sum( et*art )


        call sum0d_mpi( temp_ave, master_task )
        call sum0d_mpi( salt_ave, master_task )
        call sum0d_mpi( elev_ave, master_task )
        call sum0d_mpi(  vol_tot, master_task )
        call sum0d_mpi( area_tot, master_task )
        call sum0d_mpi(  kin_ave, master_task )

        if ( mode == MODE_DIAGNOSTIC ) then
          u_resid = max( maxval(abs(u-ub)), maxval(abs(v-vb)) )
          e_resid = maxval(abs(et-etb))
          call max0d_mpi( u_resid, master_task )
          call max0d_mpi( e_resid, master_task )
        end if

! print averages
        if ( is_master ) then

          temp_ave = temp_ave / vol_tot
          salt_ave = salt_ave / vol_tot
          elev_ave = elev_ave / area_tot
          kin_ave  = .5_rk*kin_ave / vol_tot

          if ( mode /= MODE_DIAGNOSTIC ) then
            print '(3(2(a,e15.8),a,/))'
     &        , "AVG: et = ", elev_ave        , " m  , KE = "
     &        , kin_ave         , " m^2/s^2;"
     &        , "   : tb = ", temp_ave + tbias, " deg, sb = "
     &        , salt_ave + sbias, " psu    ;"
     &        , "TOT: vol = ", vol_tot, " m^3, area = "
     &        , area_tot, " m^2. "
          else
            print '(a,e15.8,2(a,f11.8),a)'
     &        , "residual: elev = ", e_resid, " m, velocity = "
     &        , u_resid, " m/s"
          end if

          call findpsi

        end if

      end if


      end ! subroutine print_section
!
!______________________________________________________________________
!
      subroutine check_velocity
!----------------------------------------------------------------------
!  Checks if velocity condition is violated
!----------------------------------------------------------------------
! called by: advance [advance.f]
!______________________________________________________________________
!
      use config     , only: vmaxl
      use glob_const , only: rk
      use glob_domain
      use model_run  , only: iext, iint, time
      use glob_ocean , only: vaf
      use glob_out   , only: iprint

      implicit none

      integer  i,j,imax,jmax
      real(rk) vamax


      vamax = 0.

      do j=1,jm
        do i=1,im
          if ( abs(vaf(i,j)) /= vamax ) then
            vamax = abs(vaf(i,j))
            imax = i
            jmax = j
          end if
        end do
      end do

      if ( vamax > vmaxl ) then
        if ( error_status == 0 ) print '(/
     $    ''Error: velocity condition violated @ processor '',i3,/
     $    ''time ='',f9.4,
     $    '', iint ='',i8,'', iext ='',i8,'', iprint ='',i8,/
     $    ''vamax ='',e12.3,''   imax,jmax ='',2i5)'
     $    , my_task,time,iint,iext,iprint,vamax,imax,jmax
        error_status = 1
      end if


      end ! subroutine check_velocity
!
!______________________________________________________________________
!
      subroutine store_mean
!----------------------------------------------------------------------
!  Stores averages for further output (if enabled)
!----------------------------------------------------------------------
! called by: advance [advance.f]
!______________________________________________________________________
!
      use air        , only: wssurf, wtsurf, wusurf, wvsurf, swrad
      use glob_domain, only: km
      use glob_ocean , only: aam, cbc  , elb  , kh
     &                     , kmt, rho  , s    , t
     &                     , u  , uab  , v    , vab
     &                     , w  , wubot, wvbot
      use glob_out

      implicit none


      uab_mean    = uab_mean    + uab
      vab_mean    = vab_mean    + vab
      elb_mean    = elb_mean    + elb
      wusurf_mean = wusurf_mean + wusurf
      wvsurf_mean = wvsurf_mean + wvsurf
      wtsurf_mean = wtsurf_mean + wtsurf
      wssurf_mean = wssurf_mean + wssurf
      swrad_mean  = swrad_mean  + swrad
      u(:,:,km)   = wubot(:,:)        !fhx:20110318:store wvbot
      v(:,:,km)   = wvbot(:,:)        !fhx:20110318:store wvbot
      u_mean      = u_mean      + u
      v_mean      = v_mean      + v
      w_mean      = w_mean      + w
      t_mean      = t_mean      + t
      s_mean      = s_mean      + s
      rho_mean    = rho_mean    + rho
      kh_mean     = kh_mean     + kh
      aam(:,:,km) = cbc(:,:)          !lyo:20110315:botwavedrag:store cbc
      aam_mean    = aam_mean    + aam

      num = num + 1


      end ! subroutine store_mean
!
!______________________________________________________________________
!
      subroutine store_surf_mean
!----------------------------------------------------------------------
!  Stores averages for surface output
!----------------------------------------------------------------------
! called by: advance [advance.f]
!______________________________________________________________________
!
      use air       , only: uwsrf, vwsrf
      use glob_ocean, only: elb, u, v
      use glob_out

      implicit none


      usrf_mean  = usrf_mean  + u(:,:,1)
      vsrf_mean  = vsrf_mean  + v(:,:,1)
      elsrf_mean = elsrf_mean + elb
      uwsrf_mean = uwsrf_mean + uwsrf
      vwsrf_mean = vwsrf_mean + vwsrf

      nums = nums + 1


      end ! subroutine store_surf_mean
!______________________________________________________________________
!
      subroutine check_nan
!----------------------------------------------------------------------
!  Checks if NaNs present
!----------------------------------------------------------------------
! called by: advance    [advance.f]
!
! calls    : detect_nan [advance.f]
!______________________________________________________________________
!
      use glob_ocean, only: s, t, u, v

      implicit none


      call detect_nan( u, "u" )
      call detect_nan( v, "v" )
      call detect_nan( t, "t" )
      call detect_nan( s, "s" )


      end ! subroutine check_nan
!______________________________________________________________________
!
      subroutine detect_nan( var, varname )
!----------------------------------------------------------------------
!  Checks an array for NaNs
!----------------------------------------------------------------------
! called by: check_nan [advance.f]
!______________________________________________________________________
!
      use ieee_arithmetic, only: ieee_is_nan
      use glob_const , only: rk
      use glob_domain, only: i_global, im, j_global, jm, km
      use grid       , only: h

      implicit none

      real(rk)        , intent(in) :: var(im,jm,km)
      character(len=*), intent(in) :: varname

      integer i, j, k, num_nan


      num_nan = 0

      do k=1,km
        do j=1,jm
          do i=1,im
            if ( var(i,j,k) == var(i,j,k)+1 ) then
              print '(2a,3i4,2f12.4)'
     &            , "detect nan : ",varname
     &            , i_global(i), j_global(j), k
     &            , var(i,j,k), h(i,j)
              if ( k == 1 ) num_nan = num_nan + 1
            end if
          end do
        end do
      end do

      if ( num_nan /= 0 ) then
        print '(2a,2(a,i6))'
     &      , " detect_nan : ", varname
     &      , "j_global(1) = ", j_global(1)
     &      , ",   num_nan = ", num_nan
!         call finalize_mpi
        stop
      end if


      end ! subroutine detect_nan
!
!______________________________________________________________________
!fhx:tide:debug
      subroutine check_nan_2d
!----------------------------------------------------------------------
!  Checks for NaNs present in 2D
!----------------------------------------------------------------------
! called by: advance       [advance.f]
!
! calls    : detect_nan_2d [advance.f]
!______________________________________________________________________
!
      use glob_ocean, only: elf, uaf, vaf

      implicit none


      call detect_nan_2d( uaf, "uaf" )
      call detect_nan_2d( vaf, "vaf" )
      call detect_nan_2d( elf, "elf" )


      end ! subroutine check_nan_2d
!
!______________________________________________________________________
!fhx:tide;debug
      subroutine detect_nan_2d( var, varname )
!----------------------------------------------------------------------
!  Checks a 2D array for NaNs
!----------------------------------------------------------------------
! called by: check_nan_2d [advance.f]
!______________________________________________________________________
!
      use air
      use glob_const , only: rk
      use glob_domain, only: i_global, im, j_global, jm
      use grid       , only: fsm, h
      use model_run  , only: time
      use glob_ocean

      implicit none

      integer i, j, num_nan
      real(kind=rk), intent(in) :: var(im,jm)
      character(len=*),intent(in)  :: varname
!      logical isnanf


      num_nan = 0

      do j=1,jm
        do i=1,im
          if ( var(i,j) == var(i,j)+1 .or. var(i,j) /= var(i,j) ) then
            print '(2a,2i4,3f12.4)',
     $            "detect nan : ",varname,
     $            i_global(i),j_global(j),
     $            var(i,j),h(i,j),time
                  num_nan = num_nan + 1
          end if
        end do
      end do

      if ( num_nan /= 0 ) then

        print'(2a,2(a,i6))',
     $        " detect_nan : ", varname,
     $        "j_global(1) = ", j_global(1),
     $        ",   num_nan = ", num_nan
!         call finalize_mpi
         stop

      end if


      end ! subroutine detect_nan_2d
!
!______________________________________________________________________
!
      subroutine pgscheme(npg)
!----------------------------------------------------------------------
!  Redirects to a proper PGF scheme
!----------------------------------------------------------------------
! called by: lateral_viscosity [advance.f]
!            update_initial    [initialize.f]
!
! calls    : baropg            [solver.f]
!            baropg_lin        [solver.f]
!            baropg_mcc        [solver.f]
!            baropg_shch       [solver.f]
!            baropg_song_std   [solver.f]
!______________________________________________________________________
!
        implicit none

        integer, intent(in) :: npg


        select case (npg)
          case (1)
            call baropg
          case (2)
            call baropg_mcc
          case (3)
            call baropg_lin
          case (4)
            call baropg_song_std
          case (5)
            call baropg_shch
          case default
            call baropg_mcc
        end select


      end ! subroutine pgscheme
!
!______________________________________________________________________
!
      subroutine findpsi
!----------------------------------------------------------------------
!  Calculates the stream function, first assuming zero on the southern
! boundary and then, using the values on the western boundary,
! the stream function is calculated again. If the elevation field
! is near steady state, the two calculations should agree;
! otherwise not.
!  TODO: Modify for MPI
!----------------------------------------------------------------------
! called by: TODO: fill in
!
! calls    : TODO: fiil in
!______________________________________________________________________
!
      use glob_const , only: rk
      use glob_domain, only: im, imm1, jm, jmm1
      use grid       , only: dx, dy
      use glob_ocean , only: d, uab, vab

      implicit none

      real(rk), dimension(im,jm) :: psi_u, psi_v
      integer i, j


      psi_u = 0.
      psi_v = 0.

! Sweep northward:
      do j = 2, jmm1
        do i = 2, im
          psi_u(i,j+1) = psi_u(i,j)
     &                 + .25_rk*uab(i,j)*(  d(i,j) +  d(i-1,j) )
     &                                  *( dy(i,j) + dy(i-1,j) )
        end do
      end do

! Sweep eastward:
      do j = 2, jm
        do i = 2, imm1
          psi_v(i+1,j) = psi_v(i,j)
     &                 - .25_rk*vab(i,j)*(  d(i,j) +  d(i,j-1) )
     &                                  *( dx(i,j) + dx(i,j-1) )
        end do
      end do

! Calculate the difference
      psi_u = abs( psi_u - psi_v )

      print *, "Streamfunction stability index: ", maxval(psi_u)


      end subroutine findpsi
