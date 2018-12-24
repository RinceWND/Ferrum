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
        use config   , only: calc_ice
        use model_run, only: iext, isplit
     &                     , update_time
        use seaice

! advance POM 1 step in time
        implicit none


! get time
        call update_time

! set time dependent boundary conditions
        call update_bc

! set lateral viscosity
        call lateral_viscosity

! form vertical averages of 3-D fields for use in external (2-D) mode
        call mode_interaction

! external (2-D) mode calculation
        do iext=1,isplit
          call mode_external
          call check_nan_2d  !fhx:tide:debug
          if (calc_ice) call ice_advance
        end do

! internal (3-D) mode calculation
        call mode_internal

! print section
        call print_section

! check nan
        call check_nan

! store mean 2010/4/26
        call store_mean

! store SURF mean
        call store_surf_mean  !fhx:20110131:

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
!______________________________________________________________________
!
      use air      , only: air_step => step
      use bry      , only: bry_step => step
      use model_run, only: dtime

      implicit none


      call air_step( dtime )
      call bry_step( dtime )


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
      use config     , only: aam_init, horcon, mode, n1d, npg
      use glob_domain, only: im, imm1, jm, jmm1, kbm1
      use glob_grid  , only: dx, dy
      use glob_misc  , only: aamfac
      use glob_ocean , only: a, aam, c, ee, u, v

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
      if ( mode /= 2 ) then

        call advct(a,c,ee)

        call pgscheme(npg)

!lyo:scs1d:
        if ( n1d /= 0 ) then
          aam(:,:,:) = aam_init
        else
          do k=1,kbm1
            do j=2,jmm1
              do i=2,imm1
                aam(i,j,k) = horcon*dx(i,j)*dy(i,j)*aamfac(i,j)       !fhx:incmix
     $                    *sqrt( ((u(i+1,j,k)-u(i,j,k))/dx(i,j))**2
     $                          +((v(i,j+1,k)-v(i,j,k))/dy(i,j))**2
     $                    +.5*(.25*(u(i,j+1,k)+u(i+1,j+1,k)
     $                                 -u(i,j-1,k)-u(i+1,j-1,k))
     $                    /dy(i,j)
     $                    +.25*(v(i+1,j,k)+v(i+1,j+1,k)
     $                           -v(i-1,j,k)-v(i-1,j+1,k))
     $                    /dx(i,j)) **2)
              end do
            end do
          end do

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

        call exchange3d_mpi(aam(:,:,1:kbm1),im,jm,kbm1)

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
      use glob_domain, only: im, jm, kbm1
      use glob_grid  , only: dz
      use glob_ocean
      use model_run  , only: isp2i, ispi

      implicit none

      integer i,j,k


      if ( mode /= 2 ) then

        adx2d = 0.
        ady2d = 0.
        drx2d = 0.
        dry2d = 0.
        aam2d = 0.

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              adx2d(i,j) = adx2d(i,j) +  advx(i,j,k)*dz(k)
              ady2d(i,j) = ady2d(i,j) +  advy(i,j,k)*dz(k)
              drx2d(i,j) = drx2d(i,j) + drhox(i,j,k)*dz(k)
              dry2d(i,j) = dry2d(i,j) + drhoy(i,j,k)*dz(k)
              aam2d(i,j) = aam2d(i,j) +   aam(i,j,k)*dz(k)
            end do
          end do
        end do

        call advave(tps)

        adx2d = adx2d - advua
        ady2d = ady2d - advva

      end if

      egf = el*ispi

      do j=1,jm
        do i=2,im
          utf(i,j)=ua(i,j)*(d(i,j)+d(i-1,j))*isp2i
        end do
      end do
      do j=2,jm
        do i=1,im
          vtf(i,j)=va(i,j)*(d(i,j)+d(i,j-1))*isp2i
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
      use air        , only: e_atmos, vfluxf, wusurf, wvsurf
      use bry        , only: bc_vel_ext, bc_zeta
      use config     , only: alpha, ispadv, smoth
      use glob_const , only: grav
      use glob_domain, only: im, imm1, jm, jmm1!, my_task
      use glob_grid  , only: art, aru, arv, cor, dx, dy, fsm, h
      use glob_ocean
      use model_run  , only: dte, dte2, iext, isp2i, ispi, isplit!,iint

      implicit none

      integer i,j


      do j=2,jm
        do i=2,im
          fluxua(i,j)=.25*(d(i,j)+d(i-1,j))
     $                 *(dy(i,j)+dy(i-1,j))*ua(i,j)
          fluxva(i,j)=.25*(d(i,j)+d(i,j-1))
     $                 *(dx(i,j)+dx(i,j-1))*va(i,j)
        end do
      end do

! NOTE addition of surface freshwater flux, w(i,j,1)=vflux, compared
! with pom98.f. See also modifications to subroutine vertvl
      do j=2,jmm1
        do i=2,imm1
          elf(i,j)=elb(i,j)
     $              +dte2*(-(fluxua(i+1,j)-fluxua(i,j)
     $                      +fluxva(i,j+1)-fluxva(i,j))/art(i,j)
     $                      -vfluxf(i,j))
        end do
      end do

      call bc_zeta ! bcond(1)

      call exchange2d_mpi(elf,im,jm)

      if ( mod(iext,ispadv) == 0 ) call advave(tps)

      do j=2,jmm1
        do i=2,im
          uaf(i,j)=adx2d(i,j)+advua(i,j)
     $              -aru(i,j)*.25
     $                *(cor(i,j)*d(i,j)*(va(i,j+1)+va(i,j))
     $                 +cor(i-1,j)*d(i-1,j)*(va(i-1,j+1)+va(i-1,j)))
     $              +.25*grav*(dy(i,j)+dy(i-1,j))
     $                *(d(i,j)+d(i-1,j))
     $                *((1.-2.*alpha)
     $                   *(el(i,j)-el(i-1,j))
     $                  +alpha*(elb(i,j)-elb(i-1,j)
     $                         +elf(i,j)-elf(i-1,j))
     $                  +(e_atmos(i,j)-e_atmos(i-1,j)))
     $              +drx2d(i,j)+aru(i,j)*(wusurf(i,j)-wubot(i,j))
        end do
      end do

      do j=2,jmm1
        do i=2,im
          uaf(i,j)=((h(i,j)+elb(i,j)+h(i-1,j)+elb(i-1,j))
     $                *aru(i,j)*uab(i,j)
     $              -4.*dte*uaf(i,j))
     $             /((h(i,j)+elf(i,j)+h(i-1,j)+elf(i-1,j))
     $                 *aru(i,j))
        end do
      end do

      do j=2,jm
        do i=2,imm1
          vaf(i,j)=ady2d(i,j)+advva(i,j)
     $              +arv(i,j)*.25
     $                *(cor(i,j)*d(i,j)*(ua(i+1,j)+ua(i,j))
     $                 +cor(i,j-1)*d(i,j-1)*(ua(i+1,j-1)+ua(i,j-1)))
     $              +.25*grav*(dx(i,j)+dx(i,j-1))
     $                *(d(i,j)+d(i,j-1))
     $                *((1.-2.*alpha)*(el(i,j)-el(i,j-1))
     $                  +alpha*(elb(i,j)-elb(i,j-1)
     $                         +elf(i,j)-elf(i,j-1))
     $                  +(e_atmos(i,j)-e_atmos(i,j-1)))
     $              +dry2d(i,j)+arv(i,j)*(wvsurf(i,j)-wvbot(i,j))
        end do
      end do

      do j=2,jm
        do i=2,imm1
          vaf(i,j)=((h(i,j)+elb(i,j)+h(i,j-1)+elb(i,j-1))
     $                *vab(i,j)*arv(i,j)
     $              -4.*dte*vaf(i,j))
     $             /((h(i,j)+elf(i,j)+h(i,j-1)+elf(i,j-1))
     $                 *arv(i,j))
        end do
      end do

      call bc_vel_ext ! bcond(2)

      call exchange2d_mpi(uaf,im,jm)
      call exchange2d_mpi(vaf,im,jm)

      if     ( iext == (isplit-2) ) then

        etf = .25*smoth*elf

      elseif ( iext == (isplit-1) ) then

        etf = etf + .5*(1.-.5*smoth)*elf

      elseif ( iext ==  isplit    ) then

        etf = ( etf + .5*elf )*fsm

      end if

! apply filter to remove time split
      ua = ua + .5*smoth*( uab + uaf - 2.*ua )
      va = va + .5*smoth*( vab + vaf - 2.*va )
      el = el + .5*smoth*( elb + elf - 2.*el )

      elb = el
      el  = elf
      d   = h + el
      uab = ua
      ua  = uaf
      vab = va
      va  = vaf

      if ( iext /= isplit ) then

        egf = egf + el*ispi

        do j=1,jm
          do i=2,im
            utf(i,j)=utf(i,j)+ua(i,j)*(d(i,j)+d(i-1,j))*isp2i
          end do
        end do
        do j=2,jm
          do i=1,im
            vtf(i,j)=vtf(i,j)+va(i,j)*(d(i,j)+d(i,j-1))*isp2i
          end do
        end do

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
      use glob_const , only: rk, small
      use config     , only: mode, nadv, nbcs, nbct, do_restart
     &                     , smoth, t_lo, t_hi
      use glob_domain
      use glob_grid  , only: dx, dy, dz, fsm, h, zz
      use glob_ocean
      use model_run

      implicit none

      integer i,j,k
      real(rk) dxr,dxl,dyt,dyb


      if ( (iint/=1 .or. time0/=0.) .and. mode/=2 ) then

! adjust u(z) and v(z) such that depth average of (u,v) = (ua,va)
        tps = 0.

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              tps(i,j)=tps(i,j)+u(i,j,k)*dz(k)
            end do
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=2,im
              u(i,j,k)=(u(i,j,k)-tps(i,j))+
     $                 (utb(i,j)+utf(i,j))/(dt(i,j)+dt(i-1,j))
            end do
          end do
        end do

        do j=1,jm
          do i=1,im
            tps(i,j)=0.
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              tps(i,j)=tps(i,j)+v(i,j,k)*dz(k)
            end do
          end do
        end do

        do k=1,kbm1
          do j=2,jm
            do i=1,im
              v(i,j,k)=(v(i,j,k)-tps(i,j))+
     $                 (vtb(i,j)+vtf(i,j))/(dt(i,j)+dt(i,j-1))
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

        call exchange3d_mpi(w,im,jm,kb)

! set uf and vf to zero
        uf = 0.
        vf = 0.

! calculate q2f and q2lf using uf, vf, a and c as temporary variables
        call advq(q2b,q2,uf,a,c)
        call advq(q2lb,q2l,vf,a,c)
        call profq(a,c,tps,dtef)

! an attempt to prevent underflow (DEBUG)
        where(q2l.lt..5*small) q2l = .5*small
        where(q2lb.lt..5*small) q2lb = .5*small

        call bc_turb ! bcond(6)

        call exchange3d_mpi(uf(:,:,2:kbm1),im,jm,kbm2)
        call exchange3d_mpi(vf(:,:,2:kbm1),im,jm,kbm2)

        q2  = q2  + .5*smoth*( q2b  -2.*q2  + uf )
        q2l = q2l + .5*smoth*( q2lb -2.*q2l + vf )
        q2b = q2
        q2  = uf
        q2lb= q2l
        q2l = vf


! calculate tf and sf using uf, vf, a and c as temporary variables
!        if( mode /= 4 .and. ( iint > 2 .or. do_restart ) ) then
        if( mode /= 4 ) then

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


          call bc_ts ! bcond(4)
          if ( t_lo > -999. ) then
            where ( uf < t_lo ) uf = t_lo
          end if
          if ( t_hi <  999. ) then
            where ( uf > t_hi ) uf = t_hi
          end if


          call exchange3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
          call exchange3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

          t = t + .5*smoth*( tb + uf -2.*t )
          s = s + .5*smoth*( sb + vf -2.*s )
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

        call exchange3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
        call exchange3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

        tps = 0.

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              tps(i,j)=tps(i,j)
     $                  +(uf(i,j,k)+ub(i,j,k)-2.*u(i,j,k))*dz(k)
            end do
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              u(i,j,k)=u(i,j,k)
     $                  +.5*smoth*(uf(i,j,k)+ub(i,j,k)
     $                               -2.*u(i,j,k)-tps(i,j))
            end do
          end do
        end do

        tps = 0.

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              tps(i,j)=tps(i,j)
     $                  +(vf(i,j,k)+vb(i,j,k)-2.*v(i,j,k))*dz(k)
            end do
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              v(i,j,k)=v(i,j,k)
     $                  +.5*smoth*(vf(i,j,k)+vb(i,j,k)
     $                               -2.*v(i,j,k)-tps(i,j))
            end do
          end do
        end do

        ub = u
        u  = uf
        vb = v
        v  = vf

      end if

      egb = egf
      etb = et
      et  = etf
      dt  = h + et
      utb = utf
      vtb = vtf
      vfluxb = vfluxf

! calculate real w as wr
      wr = 0.

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tps(i,j) = zz(k)*dt(i,j) + et(i,j)
          end do
        end do
        do j=2,jmm1
          do i=2,imm1
            dxr=2.0/(dx(i+1,j)+dx(i,j))
            dxl=2.0/(dx(i,j)+dx(i-1,j))
            dyt=2.0/(dy(i,j+1)+dy(i,j))
            dyb=2.0/(dy(i,j)+dy(i,j-1))
            wr(i,j,k)=0.5*(w(i,j,k)+w(i,j,k+1))+0.5*
     $                (u(i+1,j,k)*(tps(i+1,j)-tps(i,j))*dxr+
     $                 u(i,j,k)*(tps(i,j)-tps(i-1,j))*dxl+
     $                 v(i,j+1,k)*(tps(i,j+1)-tps(i,j))*dyt+
     $                 v(i,j,k)*(tps(i,j)-tps(i,j-1))*dyb)
     $                +(1.0+zz(k))*(etf(i,j)-etb(i,j))/dti2
          end do
        end do
      end do

      call exchange3d_mpi(wr(:,:,1:kbm1),im,jm,kbm1)

      do k=1,kb
        do i=1,im
          if(n_south.eq.-1) wr(i,1,k)=wr(i,2,k)
          if(n_north.eq.-1) wr(i,jm,k)=wr(i,jmm1,k)
        end do
      end do
      do k=1,kb
        do j=1,jm
          if(n_west.eq.-1) wr(1,j,k)=wr(2,j,k)
          if(n_east.eq.-1) wr(im,j,k)=wr(imm1,j,k)
        end do
      end do

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            wr(i,j,k)=fsm(i,j)*wr(i,j,k)
          end do
        end do
      end do


      end ! subroutine mode_internal
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
!            finalize_mpi [mpi]
!            sum0d_mpi    [parallel_mpi.f]
!            sum0i_mpi    [parallel_mpi.f]
!______________________________________________________________________
!
      use config     , only: sbias, tbias
      use glob_const , only: rk
      use glob_domain
      use glob_grid  , only: art, dz, fsm
      use glob_ocean , only: et, dt, sb, tb
      use glob_out   , only: iprint
      use model_run

      implicit none

      real(rk) area_tot, vol_tot, d_vol
      real(rk) elev_ave, temp_ave, salt_ave
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
        vol_tot  = 0.
        area_tot = 0.
        temp_ave = 0.
        salt_ave = 0.
        elev_ave = 0.

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              d_vol    = art(i,j) * dt(i,j) * dz(k) * fsm(i,j)
              vol_tot  = vol_tot + d_vol
              temp_ave = temp_ave + tb(i,j,k)*d_vol
              salt_ave = salt_ave + sb(i,j,k)*d_vol
            end do
          end do
        end do

        do j=1,jm
          do i=1,im
            area_tot = area_tot + art(i,j)
            elev_ave = elev_ave + et(i,j) * art(i,j)
          end do
        end do


        call sum0d_mpi( temp_ave, master_task )
        call sum0d_mpi( salt_ave, master_task )
        call sum0d_mpi( elev_ave, master_task )
        call sum0d_mpi(  vol_tot, master_task )
        call sum0d_mpi( area_tot, master_task )

! print averages
        if ( is_master ) then

          temp_ave = temp_ave / vol_tot
          salt_ave = salt_ave / vol_tot
          elev_ave = elev_ave / area_tot
          print '(a,e15.8,2(a,f11.8),a)'
     &        , "mean ; et = ", elev_ave, " m, tb = "
     &        , temp_ave + tbias, " deg, sb = "
     &        , salt_ave + sbias, " psu"

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
      use air        , only: wssurf, wtsurf, wusurf, wvsurf
      use glob_domain, only: kb
      use glob_ocean , only: aam, cbc, elb, kh, rho, s, t, u, uab
     &                     , v, vab, w, wubot, wvbot
      use glob_out

      implicit none


      uab_mean    = uab_mean    + uab
      vab_mean    = vab_mean    + vab
      elb_mean    = elb_mean    + elb
      wusurf_mean = wusurf_mean + wusurf
      wvsurf_mean = wvsurf_mean + wvsurf
      wtsurf_mean = wtsurf_mean + wtsurf
      wssurf_mean = wssurf_mean + wssurf
      u(:,:,kb)   = wubot(:,:)        !fhx:20110318:store wvbot
      v(:,:,kb)   = wvbot(:,:)        !fhx:20110318:store wvbot
      u_mean      = u_mean      + u
      v_mean      = v_mean      + v
      w_mean      = w_mean      + w
      t_mean      = t_mean      + t
      s_mean      = s_mean      + s
      rho_mean    = rho_mean    + rho
      kh_mean     = kh_mean     + kh
      aam(:,:,kb) = cbc(:,:)          !lyo:20110315:botwavedrag:store cbc
      km_mean     = km_mean     + aam !lyo:20110202:save aam inst. of km

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
!
!_______________________________________________________________________
!      subroutine write_output( d_in )
!
!      use module_time
!
!      implicit none
!      include 'pom.h'
!
!      type(date), intent(in) :: d_in
!
!      integer i,j,k
!      real(kind=rk) u_tmp, v_tmp
!
!      if(netcdf_file.ne.'nonetcdf' .and. mod(iint,iprint).eq.0) then
!
!
!         uab_mean    = uab_mean    / real ( num )
!         vab_mean    = vab_mean    / real ( num )
!         elb_mean    = elb_mean    / real ( num )
!         wusurf_mean = wusurf_mean / real ( num )
!         wvsurf_mean = wvsurf_mean / real ( num )
!         wtsurf_mean = wtsurf_mean / real ( num )
!         wssurf_mean = wssurf_mean / real ( num )
!         u_mean      = u_mean      / real ( num )
!         v_mean      = v_mean      / real ( num )
!         w_mean      = w_mean      / real ( num )
!         t_mean      = t_mean      / real ( num )
!         s_mean      = s_mean      / real ( num )
!         rho_mean    = rho_mean    / real ( num )
!         kh_mean     = kh_mean     / real ( num )
!         km_mean     = km_mean     / real ( num )
!
!
!
!!         if ( my_task == 41 )
!!     $        print*, im/2,jm/2,rot(im/2,jm/2),
!!     $        uab_mean(im/2,jm/2),vab_mean(im/2,jm/2)
!!         do j = 1, jm
!!            do i = 1, im
!!               u_tmp = uab_mean(i,j)
!!               v_tmp = vab_mean(i,j)
!!               uab_mean(i,j)
!!     $              = u_tmp * cos( rot(i,j) * deg2rad )
!!     $              - v_tmp * sin( rot(i,j) * deg2rad )
!!               vab_mean(i,j)
!!     $              = u_tmp * sin( rot(i,j) * deg2rad )
!!     $              + v_tmp * cos( rot(i,j) * deg2rad )
!!            enddo
!!         enddo
!!         if ( my_task == 41 )
!!     $        print*, im/2,jm/2,
!!     $        cos(rot(im/2,jm/2)*deg2rad),
!!     $        uab_mean(im/2,jm/2),vab_mean(im/2,jm/2)
!
!
!!         do j = 1, jm
!!            do i = 1, im
!!               u_tmp = wusurf_mean(i,j)
!!               v_tmp = wvsurf_mean(i,j)
!!               wusurf_mean(i,j)
!!     $              = u_tmp * cos( rot(i,j) * deg2rad )
!!     $              - v_tmp * sin( rot(i,j) * deg2rad )
!!               wvsurf_mean(i,j)
!!     $              = u_tmp * sin( rot(i,j) * deg2rad )
!!     $              + v_tmp * cos( rot(i,j) * deg2rad )
!!            enddo
!!         enddo
!!         do k=1,kbm1
!!            do j = 1, jm
!!               do i = 1, im
!!                  u_tmp = u_mean(i,j,k)
!!                  v_tmp = v_mean(i,j,k)
!!                  u_mean(i,j,k)
!!     $                 = u_tmp * cos( rot(i,j) * deg2rad )
!!     $                 - v_tmp * sin( rot(i,j) * deg2rad )
!!                  v_mean(i,j,k)
!!     $                 = u_tmp * sin( rot(i,j) * deg2rad )
!!     $                 + v_tmp * cos( rot(i,j) * deg2rad )
!!               enddo
!!            enddo
!!         enddo
!
!
!         write( filename, '("out/",2a,".nc")' )
!     $        trim( netcdf_file ), date2str( d_in )
!
!         call write_output_pnetcdf( filename )
!
!         uab_mean    = 0.0
!         vab_mean    = 0.0
!         elb_mean    = 0.0
!         wusurf_mean = 0.0
!         wvsurf_mean = 0.0
!         wtsurf_mean = 0.0
!         wssurf_mean = 0.0
!         u_mean      = 0.0
!         v_mean      = 0.0
!         w_mean      = 0.0
!         t_mean      = 0.0
!         s_mean      = 0.0
!         rho_mean    = 0.0
!         kh_mean     = 0.0
!         km_mean     = 0.0
!
!         num = 0
!
!      endif
!
!      return
!      end
!
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
      use glob_const , only: rk
      use glob_domain, only: i_global, im, j_global, jm, kb
      use glob_grid  , only: h

      implicit none

      real(rk)        , intent(in) :: var(im,jm,kb)
      character(len=*), intent(in) :: varname

      integer i, j, k, num_nan


      num_nan = 0

      do k=1,kb
        do j=1,jm
          do i=1,im
            if ( isnan(var(i,j,k)) ) then
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
      use glob_grid  , only: h
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
          if ( isnan(var(i,j)) ) then
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
      endif


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
