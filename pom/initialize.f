! initialize.f

!  Initializes POM: defines constants, reads initial values, grid and
! initial conditions.

!______________________________________________________________________
!
      subroutine initialize
!----------------------------------------------------------------------
!  Initializes POM.
!______________________________________________________________________
!
      use config
      use glob_const , only: initialize_constants, rk
      use glob_domain
      use grid       , only: east_e, north_e, h
      use glob_misc  , only: aamfac
      use glob_ocean , only: d, dt, el, et
      use glob_out   , only: iout, iprint
      use model_run  , only: dtime, initialize_time, time, time0

      implicit none

      integer, parameter :: TEST_CASE = 0 !1


! Read domain distribution (TODO: get rid of hardcoded parameters)
      call read_domain_dist( 'domain.nml' )

! Read configuration files
      call read_config( 'config.nml' )
!______________________________________________________________________
!  `read_config` initializes all parameters to their default values
! prior to reading them.

! initialize the MPI execution environment and create communicator for
!internal POM communications
      call initialize_mpi

! distribute the model domain across processors
      call distribute_mpi

! distribute the coarse wind/assim data domain across processors
      call distribute_mpi_coarse   !fhx:20101206

! decide on interpolation feasibility (sets calc_interp flag)
      call check_interpolation

! Initialize various constants
      call initialize_constants

! Allocate and initialize arrays
      call initialize_arrays

! Initialize general IO and config
! [TODO] This routine should initialize output, not io
!      call initialize_io( 'config.nml' )

! Initialize time
      if ( do_restart ) then
        call initialize_time( init_file )
      else
        call initialize_time( "" )
      end if
!______________________________________________________________________
!  Uses `restart_file` string to offset time.
!  If `do_restart` is .false. this string is forced to be empty, thus
! the offseting is ignored.

! Print configuration
      call print_config

! Initialize grid
      call initialize_grid( TEST_CASE )

! The model is initialized. Below is some updates. TODO: Separate?
      call msg_print("MODEL CORE INITIALIZED", 1, "")

! Initialize modules
      call initialize_modules

! read previous (initial) record [TODO: NOTE! This step should be done before update_initial for clim to get rmean. But what happens to boundary?]
      call modules_initial_step( dtime )

! read initial and lateral boundary conditions
      call initial_conditions

! read laterel boundary conditions
      call lateral_boundary_conditions !ayumi:20100407

! read M2 tidal amplitude & phase
      if ( calc_tide ) call read_tide      !fhx:tide

! update initial conditions and set the remaining conditions
      call update_initial

! calculate the bottom friction coefficient
      call bottom_friction

!lyo:pac10:Keep this to increase aam near (lono,lato)
!     xs=1.5; ys=1.5; fak=0.5; !lyo:pac10:Specified in pom.h
!     lono=-74.59; lato=36.70; !Correspond to lyo:coarse io=234; jo=185
      if ( lono/=999. .and. lato/=999. ) then
!set lono.OR.lato=999. to skip
        call incmix(aamfac,im,jm,1,lono,lato,east_e,north_e,xs,ys,fak)
      end if

! read restart data from a previous run
      if ( do_restart ) then
        call read_restart_pnetcdf
! update elevation
        d  = h+el
        dt = h+et
! update time
        time = time0
        if ( append_output ) iout = int(time/iprint)
      end if

! write grid and initial conditions
      if (output_flag == 1) then
        if ( netcdf_file /= 'nonetcdf' )
     $      call write_output_pnetcdf(
     $      "out/"//trim(netcdf_file)//".nc")
      end if

      if ( SURF_flag == 2 ) then  !fhx:20110131: flag=2 for initial SURF output
        if ( netcdf_file /= 'nonetcdf' )
     $      call write_SURF_pnetcdf(
     $      "out/SRF."//trim(netcdf_file)//".nc")
      end if


! check for errors
      call   sum0i_mpi(error_status,master_task)
      call bcast0i_mpi(error_status,master_task)
      if ( error_status /= 0 ) then
        if ( is_master ) print '(/a/a)'
     &                       , "Something went wrong during init"
     &                       , "POM terminated with error"
        call finalize_mpi
        stop
      end if

      call msg_print("MODEL STATE UPDATED", 1, "")


      end ! subroutine initialize
!
!______________________________________________________________________
!
      subroutine initialize_arrays
!----------------------------------------------------------------------
!  Allocate arrays and initialize them for safety.
!______________________________________________________________________
!
      use config     , only: ntide
      use glob_const , only: rk
      use glob_domain, only: im, im_coarse, jm, jm_coarse, kb, kbm1
      use grid
      use glob_misc
      use glob_ocean
      use glob_out

      implicit none

      integer k

      iout = 0

      iouts = 0 !fhx:20110131:


      allocate(
     &  dz(kb), dzz(kb), z(kb), zz(kb)
     & )

      allocate(
     &  aam2d(im,jm)   ,
     $  aamfac(im,jm)  ,    !fhx:incmix
     $  advua(im,jm)   ,
     $  advva(im,jm)   ,
     $  adx2d(im,jm)   ,
     $  ady2d(im,jm)   ,
     $  art(im,jm)     ,
     $  aru(im,jm)     ,
     $  arv(im,jm)     ,
     $  cbc(im,jm)     ,
     $  icb(im,jm)     ,    !:rwnd
     $  ice(im,jm)     ,    !:rwnd
     $  icf(im,jm)     ,    !:rwnd
     $  cor(im,jm)     ,
     $  d(im,jm)       ,
     $  drx2d(im,jm)   ,
     $  dry2d(im,jm)   ,
     $  dt(im,jm)      ,
     $  dum(im,jm)     ,
     $  dvm(im,jm)     ,
     $  dx(im,jm)      ,
     $  dy(im,jm)      ,
     $  east_c(im,jm)  ,
     $  east_e(im,jm)  ,
     $  east_u(im,jm)  ,
     $  east_v(im,jm)  ,
     $  egb(im,jm)     ,
     $  egf(im,jm)     ,
     $  el(im,jm)      ,
     $  elb(im,jm)     ,
     $  elf(im,jm)     ,
     $  et(im,jm)      ,
     $  etb(im,jm)     ,
     $  etf(im,jm)     ,
     $  fluxua(im,jm)  ,
     $  fluxva(im,jm)  ,
     $  fsm(im,jm)     ,
     $  h(im,jm)       ,
     $  hi(im,jm)      ,    !:rwnd
     $  hz(im,jm,kb)   ,
     $  north_c(im,jm) ,
     $  north_e(im,jm) ,
     $  north_u(im,jm) ,
     $  north_v(im,jm) ,
     $  psi(im,jm)     ,
     $  rot(im,jm)     ,
     $  ssurf(im,jm)   ,
     $  tauiwu(im,jm)  ,    !:rwnd
     $  tauiwv(im,jm)  ,    !:rwnd
     $  tps(im,jm)     ,
     $  tsurf(im,jm)   ,
     $  ua(im,jm)      ,
     $  uab(im,jm)     ,
     $  uaf(im,jm)     ,
     $  uib(im,jm)     ,    !:rwnd
     $  ui(im,jm)      ,    !:rwnd
     $  uif(im,jm)     ,    !:rwnd
     $  utb(im,jm)     ,
     $  utf(im,jm)     ,
     $  va(im,jm)      ,
     $  vab(im,jm)     ,
     $  vaf(im,jm)     ,
     $  vib(im,jm)     ,    !:rwnd
     $  vi(im,jm)      ,    !:rwnd
     $  vif(im,jm)     ,    !:rwnd
     $  vtb(im,jm)     ,
     $  vtf(im,jm)     ,
     $  wubot(im,jm)   ,
     $  wvbot(im,jm)   ,
     $  alon_coarse(im_coarse,jm_coarse) ,
     $  alat_coarse(im_coarse,jm_coarse) ,
     $  mask_coarse(im_coarse,jm_coarse)
     & )

      allocate(
     &  aam(im,jm,kb)  ,
     $  advx(im,jm,kb) ,
     $  advy(im,jm,kb) ,
     $  a(im,jm,kb)    ,
     $  c(im,jm,kb)    ,
     $  drhox(im,jm,kb),
     $  drhoy(im,jm,kb),
     $  dtef(im,jm,kb) ,
     $  ee(im,jm,kb)   ,
     $  gg(im,jm,kb)   ,
     $  kh(im,jm,kb)   ,
     $  km(im,jm,kb)   ,
     $  kq(im,jm,kb)   ,
     $  l(im,jm,kb)    ,
     $  q2b(im,jm,kb)  ,
     $  q2(im,jm,kb)   ,
     $  q2lb(im,jm,kb) ,
     $  q2l(im,jm,kb)  ,
     $  rho(im,jm,kb)  ,
     $  sb(im,jm,kb)   ,
     $  s(im,jm,kb)    ,
     $  tb(im,jm,kb)   ,
     $  t(im,jm,kb)    ,
     $  ub(im,jm,kb)   ,
     $  uf(im,jm,kb)   ,
     $  u(im,jm,kb)    ,
     $  vb(im,jm,kb)   ,
     $  vf(im,jm,kb)   ,
     $  v(im,jm,kb)    ,
     $  w(im,jm,kb)    ,
     $  wr(im,jm,kb)   ,
     $  zflux(im,jm,kb)
     & )

      allocate(
     &  uab_mean(im,jm)    ,
     $  vab_mean(im,jm)    ,
     $  elb_mean(im,jm)    ,
     $  wusurf_mean(im,jm) ,
     $  wvsurf_mean(im,jm) ,
     $  wtsurf_mean(im,jm) ,
     $  wssurf_mean(im,jm)
     & )

      allocate(
     $  u_mean(im,jm,kb)   ,
     $  v_mean(im,jm,kb)   ,
     $  w_mean(im,jm,kb)   ,
     $  t_mean(im,jm,kb)   ,
     $  s_mean(im,jm,kb)   ,
     $  rho_mean(im,jm,kb) ,
     $  kh_mean(im,jm,kb)  ,
     $  km_mean(im,jm,kb)
     & )

! TODO: move to seaice and tide modules
!      allocate(
!     $  cibe(jm)       ,
!     $  cibn(im)       ,
!     $  cibs(im)       ,
!     $  cibw(jm)       ,
!     $  ampe(jm,ntide) ,
!     $  phae(jm,ntide) ,
!     $  amue(jm,ntide) ,
!     $  phue(jm,ntide)
!     & )


      allocate(
     $  usrf_mean(im,jm)    ,
     $  vsrf_mean(im,jm)    ,
     $  elsrf_mean(im,jm)   ,
     $  uwsrf_mean(im,jm)   ,
     $  vwsrf_mean(im,jm)   ,
     & )

! 2-D and 3-D arrays
      elb    = 0.
      elf    = 0.
      cor    = 0.
      tauiwu = 0.
      tauiwv = 0.
      drx2d  = 0.
      dry2d  = 0.
      wubot  = 0.     !lyo:20110315:botwavedrag:
      wvbot  = 0.     !lyo:20110315:botwavedrag:
      aamfac = 1.     !fhx:incmix
      hi     =  .35

      ub(:,:,1:kbm1) = 0.
      vb(:,:,1:kbm1) = 0.

      do k = 1, kb
        hz(:,:,kb) = -h*zz(k)
      end do


      end ! subroutine initialize_arrays
!
!______________________________________________________________________
!
      subroutine initialize_grid( TEST_CASE )
!----------------------------------------------------------------------
!  Reads of generates grid.
!______________________________________________________________________
!
      use glob_domain, only: im, imm1, jm, jmm1
      use grid

      implicit none

      integer, intent(in) :: TEST_CASE


      call initialize_mod( 'config.nml' )

      if ( TEST_CASE == 0 ) then
! read in grid data
        call read_grid
      else
        call make_grid(TEST_CASE)
      end if

!lyo:ecs:
!     Close northern boundary:
!!        if(n_north.eq.-1) then
!          do j=jm-5,jm
!!           do j=jm-0,jm
!!           fsm(:,j)=0.0
!!           h(:,j)=1.0
!!          enddo
!!        endif
!

! derive u- and v-mask from rho-mask
      dum(2:im,:) = fsm(2:im,:)*fsm(1:imm1,:)
      dum(1   ,:) = dum(2   ,:)
      dvm(:,2:jm) = fsm(:,2:jm)*fsm(:,1:jmm1)
      dvm(:,1   ) = dvm(:,2   )

      call exchange2d_mpi(  h,im,jm)
      call exchange2d_mpi(fsm,im,jm)
      call exchange2d_mpi(dum,im,jm)
      call exchange2d_mpi(dvm,im,jm)

      call check_cfl_min


      end ! subroutine initialize_grid
!
!______________________________________________________________________
!
      subroutine update_grid
!----------------------------------------------------------------------
!  Set up vertical and horizontal grid, topography, areas and masks.
!______________________________________________________________________
!
      use grid
      use glob_ocean , only: d, dt, el, et

      implicit none


! read grid
      call read_grid
!______________________________________________________________________
!      The followings are read in read_grid_pnetcdf:
!     z,zz,dx,dy
!     east_u,east_v,east_e,east_c
!     north_u,north_v,north_e,north_c
!     rot,h,fsm



! initialize water column depths
      d  = h+el
      dt = h+et


      end
!
!______________________________________________________________________
!
      subroutine initial_conditions
!----------------------------------------------------------------------
!  Set up initial and lateral boundary conditions.
!______________________________________________________________________
!
      use bry
      use glob_bry
      use config     , only: do_restart, init_file, init_record
      use glob_const , only: rk
      use glob_domain
      use grid       , only: dz
      use glob_ocean , only: el, elb, rho, s, sb, t, tb
     &                     , uab, ub, vab, vb

      implicit none

      integer k


! Read initial fields
      if ( .not.do_restart ) then
        call read_initial( trim(init_file), init_record
     &                   , elb, sb, tb, ub, vb )

! Derive barotropic velocities
         uab = 0.
         vab = 0.
         do k = 1, kbm1
           uab = uab + ub(:,:,k)*dz(k)
           vab = vab + vb(:,:,k)*dz(k)
         end do
!      do k = 1,kbm1
!        sb(:,:,k) = sb(:,:,k) * fsm
!        tb(:,:,k) = tb(:,:,k) * fsm
!      end do

! Calculate initial density.
        call dens(sb,tb,rho)

!      uab = .6*dum
!      if (n_south==-1) then
!        dvm(:,1)=0.
!        dum(:,1)=0.
!      end if

!      if(n_east.eq.-1) then
!        do j=1,jm
!          do i=1,im
!           aam2d(i,j)=100.e0*(1.+10.e0*exp(-(im-i)/10.e0))
!          end do
!        end do
!      end if

        t(:,:,1:kbm1) = tb(:,:,1:kbm1)
        s(:,:,1:kbm1) = sb(:,:,1:kbm1)
        el = elb

      end if
! set thermodynamic boundary conditions (for the seamount problem, and
! other possible applications, lateral thermodynamic boundary conditions
! are set equal to the initial conditions and are held constant
! thereafter - users may create variable boundary conditions)
!lyo:pac10:Comment out - replace by tobe etc below:
!     do k=1,kbm1
!       do j=1,jm
!         tbe(j,k) = tb(im,j,k) * fsm(im,j) !lyo:20110202:use initial
!         sbe(j,k) = sb(im,j,k) * fsm(im,j) ! t,sb instead of t,sclim
!         tbw(j,k) = tb( 1,j,k) * fsm( 1,j) !lyo:20110202:add west
!         sbw(j,k) = sb( 1,j,k) * fsm( 1,j)
!       end do
!       do i=1,im
!         tbn(i,k) = tb(i,jm,k) * fsm(i,jm) !lyo:20110202:add north
!         sbn(i,k) = sb(i,jm,k) * fsm(i,jm)
!         tbs(i,k) = tb(i, 1,k) * fsm(i, 1) !lyo:20110202:add south
!         sbs(i,k) = sb(i, 1,k) * fsm(i, 1)
!       end do
!     enddo

!lyo:pac10:beg:From /lustre/ltfs/scratch/Yu-Lin.Chang/sbPOM/stcc_exp/dx10_wind_5pe-5_kb350/pom/initialize.f; NOT from Alu's stcc code which does not have below
!lyo:pac10:end:

!      call read_ice_pnetcdf( "ice.19790102.nc", icb )

! Set initial conditions for modules.
      call initial_conditions_boundary


      end
!
!______________________________________________________________________
!lyo:pac10:beg:Here thro *end: replaced subr.lateral_boundary_conditions
      subroutine lateral_boundary_conditions
!----------------------------------------------------------------------
!  Read lateral boundary conditions.
!  Transport at eastern boundary for PROFS in GOM.
! ayumi 2010/4/7
!______________________________________________________________________
!
      use bry        , only: periodic_bc
      use glob_const , only: rk
      use glob_domain, only: i_global, im, im_global
     &                     , j_global, jm, jm_global, master_task
      use grid       , only: cor

      implicit none

      logical :: here, judge_inout !lyo:scs1d:
      integer :: i,j,ic,jc         !lyo:scs1d:
      real(kind=rk)    :: corcon            !lyo:scs1d:
!      character(len=120) in_file        !eda:uvforce

!     call read_uabe_pnetcdf(uabe)
!eda:uvforce
!      if (.not. calc_uvforce) then
!        write(in_file,'(a)') "bc.nc"
!        call read_bc_pnetcdf(uabe, uabw, vabs, vabn, in_file, 1)
!      endif

! Periodic in "x" and/or "y"?  !lyo:20110224:alu:stcc:
!     iperx.ne.0 if x-periodic; ipery.ne.0 if y-periodic               !
!     iperx(y) < 0 if south/north (west/east) walls are free-slip      !
!     iperx=-1; ipery= 0 !--> x-periodic & free-slip S/N boundaries
!     iperx= 0; ipery= 0 !--> x-periodic & free-slip S/N boundaries
!lyo:scs1d:moved to namelist: pom.nml
!lyo:scs1d:cannot be beta-plane if double-periodic (note cor(*,*)
!     was previously defined in "call read_grid")
      if ( periodic_bc%x .and. periodic_bc%y ) then
      ic = (im_global+1)/2; jc = (jm_global+1)/2
         here = judge_inout( ic, jc,
     $                       i_global(1), i_global(im),
     $                       j_global(1), j_global(jm) )
         if ( here ) then
            corcon=cor(ic-i_global(1)+1,jc-j_global(1)+1)
         endif
         call bcast(corcon,1,master_task)
         do j=1,jm; do i=1,im !f-plane:
            cor(i,j)=corcon
            enddo; enddo
         endif !if (iperx.eq.1 .and. ipery.eq.1) then


      end
!
!______________________________________________________________________
!
      subroutine bfrz(mw,me,ms,mn,nw,ne,ns,nn,im,jm,nu,frz)
!----------------------------------------------------------------------
!lyo:20110224:alu:stcc:
!lyo:Modified from subroutine assimfrs in
!     /wrk/newshoni/hunglu/model/sbPOM/stcc_ideal/
!     stcc_alu_30TSrelx_60aam_timescle1d/pom/pom.f
!     which was modified from:
!     /archive/lyo/arctic/beringsea/codes_and_runscripts/
!     bering07_hass50m_faccof1p2_dtfac.f; see also
!     /wrk/aden/lyo/crwu/assimssh/how_to_put_assimssh_in_pom.txt
!lyo:Modified from /home/lyo/gom/hindcast/fgrid/gomn117.f subr.NFRSASIM
!----------------------------------------------------------------------
!     calculate boundary flow-relaxation array "frz"
!----------------------------------------------------------------------
!     nu          = unit# (fort.nu) for ASCII printout
!     mw,me,ms,mn = #buffer grid-pnts near west, east, south & north
!                   boundaries where assimilation frz=1.0; recommended
!                   buffer 100~200km;  m?=0 for no buffer;
!                   As a precaution, program stops if 0<m?<=3 (or <0)
!     nw,ne,ns,nn = n_west,n_east,n_south,n_north
!                   is =-1 if "bfrz" is being called by processor that
!                   shares the west, east, south or north boundary
!                   Just set all to -1 for a single-processor run
!     frz         = 1 at boundaries and tapered (tanh) =0  interior
!
!                ... l.oey --- lyo@princeton.edu (Jan/2008)
!----------------------------------------------------------------------

      use glob_const , only: rk

      implicit none

      integer, intent(in) :: mw,me,ms,mn,im,jm,nu
      integer, intent(in) :: nw,ne,ns,nn
      real(kind=rk), dimension(im,jm), intent(out) :: frz
      integer :: mmid,i,j,ii,jj
      integer, parameter :: my_task=0
      real(kind=rk) :: c,tanhm


      frz(:,:)=0.0       !Initialize interior

!     West:
      if(nw.eq.-1) then
       if (mw.gt.3) then  !west buffer: needs at least 4 pts
         frz(1,:)=1.0; mmid=mw/2; c=5./float(mw-mmid)
         do i=2,mw
           ii=i
           tanhm=0.5*(1.-tanh(float(i-mmid)*c))
           frz(ii,:)=tanhm
         enddo
       elseif(mw.eq.0.or.mw.eq.1) then
!      do nothing, i.e. frz remains = 0 for mw=0 or 1
       else !mw=2 or 3 or <0
        write(nu,'(''Stopped in subr.bfrz, mw ='',i4)') mw
        write( *,'(''Stopped in bfrz. proc# ,mw ='',2i4)') my_task,mw
        stop
       endif
      endif
!
!     East:
      if(ne.eq.-1) then
       if (me.gt.3) then  !east buffer:
         frz(im,:)=1.0; mmid=me/2; c=5./float(me-mmid)
         do i=2,me
           ii=im-i+1
           tanhm=0.5*(1.-tanh(float(i-mmid)*c))
           frz(ii,:)=tanhm
         enddo
       elseif(me.eq.0.or.me.eq.1) then
!      do nothing, i.e. frz remains = 0 for me=0 or 1
       else !me=2 or 3 or <0
        write(nu,'(''Stopped in subr.bfrz, me ='',i4)') me
        write( *,'(''Stopped in bfrz. proc# ,me ='',2i4)') my_task,me
        stop
       endif
      endif
!
!     South:
      if(ns.eq.-1) then
       if (ms.gt.3) then  !south buffer:
         frz(:,1)=1.0; mmid=ms/2; c=5./float(ms-mmid)
         do j=2,ms
           jj=j
           tanhm=0.5*(1.-tanh(float(j-mmid)*c))
           do i=1,im !lyo:debug:lyo:20110224:fxu:stcc:delete "if (nw.eq.-1.."
           frz(i,jj)=max(tanhm,frz(i,jj)) !takes care of SW & SE corners
           enddo
         enddo
       elseif(ms.eq.0.or.ms.eq.1) then
!      do nothing, i.e. frz remains = 0 for ms=0 or 1
       else
        write(nu,'(''Stopped in subr.bfrz, ms ='',i4)') ms
        write( *,'(''Stopped in bfrz. proc# ,ms ='',2i4)') my_task,ms
        stop
       endif
      endif
!
!     North:
      if(nn.eq.-1) then
       if (mn.gt.3) then  !north buffer:
         frz(:,jm)=1.0; mmid=mn/2; c=5./float(mn-mmid)
         do j=2,mn
           jj=jm-j+1
           tanhm=0.5*(1.-tanh(float(j-mmid)*c))
           do i=1,im !lyo:debug:lyo:20110224:fxu:stcc:delete "if (nw.eq.-1.."
           frz(i,jj)=max(tanhm,frz(i,jj)) !takes care of NW & NE corners
           enddo
         enddo
       elseif(mn.eq.0.or.mn.eq.1) then
!      do nothing, i.e. frz remains = 0 for mn=0 or 1
       else
        write(nu,'(''Stopped in subr.bfrz, mn ='',i4)') mn
        write( *,'(''Stopped in bfrz. proc# ,mn ='',2i4)') my_task,mn
        stop
       endif
      endif


      end
!lyo:pac10:end:
!
!______________________________________________________________________
!
      subroutine read_tide
!----------------------------------------------------------------------
!fhx:tide:read tidal amplitude & phase at the eastern boundary for PROFS
!______________________________________________________________________
!
        use glob_const , only: rk


!      call read_tide_east_pnetcdf(ampe,phae)
        call read_tide_east_pnetcdf(ampe,phae,amue,phue)

!      if(my_task.eq.0)print*,ampe(10,1),phae(10,1),amue(10,1),phue(10,1)
!      if(my_task.eq.0)print*,ampe(10,2),phae(10,2),amue(10,2),phue(10,2)


      end
!fhx:tide:read_tide end
!
!______________________________________________________________________
!
      subroutine update_initial
!----------------------------------------------------------------------
!  Updates the initial conditions and sets the remaining
! initial conditions.
!______________________________________________________________________
!
      use air        , only: vfluxf
      use config     , only: aam_init, npg
      use glob_const , only: rk, SMALL
      use glob_domain, only: im, is_master, jm, kb, kbm1
      use grid       , only: dz, h
      use glob_ocean , only: aam, d, drhox, drhoy, drx2d, dry2d, dt
     &                     , el, elb, et, etb, etf
     &                     , kh, km, kq, l, q2, q2b, q2l, q2lb
     &                     , rho, s, sb, t, tb, u, ua, ub, uab
     &                     , v, va, vb, vab, w

      implicit none

      integer i,j,k


      ua = uab
      va = vab
      el = elb
      et = etb
      etf= et
      d  = h+el
      dt = h+et
      w(:,:,1) = vfluxf

      do k=1,kb
        l(:,:,k) = .1*dt
      end do

      q2b = SMALL
      q2lb= l*q2b
      kh  = l*sqrt(q2b)
      km  = kh
      kq  = kh
      aam = aam_init

      do k=1,kbm1
        do i=1,im
          do j=1,jm
            q2(i,j,k)=q2b(i,j,k)
            q2l(i,j,k)=q2lb(i,j,k)
            t(i,j,k)=tb(i,j,k)
            s(i,j,k)=sb(i,j,k)
            u(i,j,k)=ub(i,j,k)
            v(i,j,k)=vb(i,j,k)
          end do
        end do
      end do

      if ( is_master ) then
        if ( npg>5 .and. npg<0 ) then
          print '(/"[!] Invalid value for Pressure Gradient Scheme")'
          print '("   Defaulting to McCalpin 4th order (npg=2)")'
        end if
      end if

      call dens(s,t,rho)

      call pgscheme(npg)

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
            dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
          end do
        end do
      end do


      end
!______________________________________________________________________
!
      subroutine initialize_modules
!----------------------------------------------------------------------
!  Initialize all neccessary modules.
!______________________________________________________________________
!
      use air   , only: initialize_air => initialize_mod
      use bry   , only: initialize_bry => initialize_mod
      use clim  , only: initialize_clm => initialize_mod
      use seaice, only: initialize_ice => initialize_mod
      use tide  , only: initialize_tide=> initialize_mod

      implicit none


      call initialize_clm( 'config.nml' )
      call initialize_bry( 'config.nml' )
      call initialize_air( 'config.nml' )
      call initialize_ice( 'config.nml' )
      call initialize_tide('config.nml' )


      end ! subroutine initialize_modules
!
!______________________________________________________________________
!
      subroutine bottom_friction
!----------------------------------------------------------------------
!  Calculates the bottom friction coefficient.
!______________________________________________________________________
!
      use config     , only: cbcmax, cbcmin, z0b
      use glob_const , only: Kappa, rk
      use glob_domain, only: im, jm, kbm1
      use grid       , only: h, zz
      use glob_ocean , only: cbc

      implicit none

      integer i,j


! calculate bottom friction
      do j=1,jm
        do i=1,im
!lyo:correct:cbc(i,j)=(Kappa/log((1.+zz(kbm1))*h(i,j)/z0b))**2 !lyo:bug:
          cbc(i,j)=(Kappa/log(1.+(1.0+zz(kbm1))*h(i,j)/z0b))**2
          cbc(i,j)=max(cbcmin,cbc(i,j))
! if the following is invoked, then it is probable that the wrong
! choice of z0b or vertical spacing has been made:
          cbc(i,j)=min(cbcmax,cbc(i,j))
        end do
      end do


      end ! subroutine bottom_friction
!
!______________________________________________________________________
!
      subroutine ztosig(zs,tb,zz,h,t,im,jm,ks,kb,
     $                                   n_west,n_east,n_south,n_north)
!----------------------------------------------------------------------
!  Interpolates vertically from z- to sigma-coordinates.
!______________________________________________________________________
!
      use glob_const , only: rk

      implicit none

      integer , intent(in ) :: im, jm, ks, kb
      real(rk), intent(in ) :: zs(ks), tb(im,jm,ks), zz(kb)
     &                       , h(im,jm)
      real(rk), intent(out) :: t(im,jm,kb)
      integer , intent(in ) :: n_west,n_east,n_south,n_north

      integer  i,j,k
      real(rk) tmax, tin(ks), tout(kb), zzh(kb)


      do j=2,jm-1
        do i=2,im-1
          if ( h(i,j) > 1. ) then
! special interp on z-lev for cases of no data because h smoothing
            do k=1,ks
              tin(k) = tb(i,j,k)
              if ( zs(k)<=h(i,j) .and. tin(k)<.01 ) then
                tmax = max(tb(i-1,j,k),tb(i+1,j,k),
     &                     tb(i,j-1,k),tb(i,j+1,k))
                tin(k) = tmax
              end if
              if ( tin(k)<.01 .and. k/=1 ) tin(k) = tin(k-1)
            end do

            do k=1,kb
              zzh(k) = -zz(k)*h(i,j)
            end do

! vertical spline interp
            call splinc(zs,tin,ks,2.e30_rk,2.e30_rk,zzh,tout,kb)

            do k=1,kb
              t(i,j,k) = tout(k)
            end do

          end if
        end do
      end do
      call exchange3d_mpi(t,im,jm,kb)

! boundaries
      do k=1,kb
        do j=1,jm
          if ( n_west == -1 ) t( 1,j,k) = t(   2,j,k)
          if ( n_east == -1 ) t(im,j,k) = t(im-1,j,k)
        end do
        do i=1,im
          if ( n_south == -1 ) t(i, 1,k) = t(i,   2,k)
          if ( n_north == -1 ) t(i,jm,k) = t(i,jm-1,k)
        end do
      end do


      end
!
!______________________________________________________________________
!
      subroutine splinc(x,y,n,yp1,ypn,xnew,ynew,m)
!----------------------------------------------------------------------
!  Interpolate using splines.
!______________________________________________________________________
!
      use glob_const , only: rk

      implicit none

      integer, parameter :: nmax = 300

      integer               , intent(in ) :: n, m
      real(rk)              , intent(in ) :: yp1, ypn
      real(rk), dimension(n), intent(in ) :: x,y
      real(rk), dimension(m), intent(out) :: xnew,ynew

      integer                      i, k
      real(rk)                     p, qn, sig, un
      real(rk), dimension(nmax) :: y2,u


      if ( yp1 > .99e30 ) then
        y2(1) = 0.
        u(1)  = 0.
      else
        y2(1) = -.5
        u(1)  = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      end if

      do i=2,n-1
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*y2(i-1)+2.
        y2(i) = (sig-1.)/p
        u(i) = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     $        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do

      if ( ypn > .99e30 ) then
        qn = 0.
        un = 0.
      else
        qn =  .5
        un = (3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      end if

      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
        y2(k) = y2(k)*y2(k+1)+u(k)
      end do

      do i=1,m
        call splint(x,y,y2,n,xnew(i),ynew(i))
      end do


      end
!
!______________________________________________________________________
!
      subroutine splint(xa,ya,y2a,n,x,y)
!----------------------------------------------------------------------
!  Spline interpolation? [ DESCRIPTION NOT PROVIDED ]
!______________________________________________________________________
!
      use glob_const , only: rk
      use glob_domain, only: error_status

      implicit none

      integer               , intent(in ) :: n
      real(rk), dimension(n), intent(in ) :: xa, ya, y2a
      real(rk)              , intent(out) :: x, y

      integer k, khi, klo
      real(rk) a, b, h


      klo = 1
      khi = n
      do while ( khi-klo > 1 )
        k = (khi+klo)/2
        if ( xa(k) > x ) then
          khi = k
        else
          klo = k
        end if
      end do
      h = xa(khi) - xa(klo)
      if ( h == 0. )  then
        error_status = 1
        print '(/a)', 'Error: bad xa input in splint'
      end if
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo)+b*ya(khi)+
     &       ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.


      end ! subroutine splint
!
!______________________________________________________________________
!fhx:incmix:add subroutine incmix
      subroutine incmix(aam,im,jm,kb,lono,lato,x,y,xs,ys,fac)
!----------------------------------------------------------------------
!  Increase 'aam' by '(1+fac)*' at (lono,lato), then taper off to '1'
! in gaussian-manner for (x-xo, y-yo) > xs and ys.
!
!     Inputs: aam,x,y,xs,ys & fac
!     Output: aam is modified
!______________________________________________________________________
!
      use glob_const , only: rk

      implicit none

      integer                      , intent(in   ) :: im, jm, kb
      real(rk)                     , intent(in   ) :: lono, lato
     &                                              , xs, ys, fac
      real(rk), dimension(im,jm,kb), intent(inout) :: aam
      real(rk), dimension(im,jm   ), intent(in   ) :: x, y

      integer  i,j,k
      real(rk) factor,expon


!      print*,'incmix:', lono,lato
      do k=1,kb
        do j=1,jm
          do i=1,im
            factor=0.0
!      expon=((x(i,j)-x(io,jo))/xs)**2+((y(i,j)-y(io,jo))/ys)**2
            expon = ((x(i,j)-lono)/xs)**2+((y(i,j)-lato)/ys)**2
            if ( expon <= 10. ) then
              factor = fac*exp(-expon)
            end if
            aam(i,j,k) = aam(i,j,k)*(1.+factor)
          end do
        end do
      end do


      end
!lyo:pac10:exp013:
!______________________________________________________________________
!
      pure logical function judge_inout( i_in, j_in,
     &                                   imin_in, imax_in,
     &                                   jmin_in, jmax_in )
!----------------------------------------------------------------------
!  If the processor has ( i_in, j_in ) in its local domain.
!______________________________________________________________________
!
      implicit none

      integer, intent(in) ::    i_in,    j_in
     &                     , imin_in, imax_in
     &                     , jmin_in, jmax_in


      if (       ( i_in >= imin_in .and. i_in <= imax_in )
     &     .and. ( j_in >= jmin_in .and. j_in <= jmax_in ) ) then

         judge_inout = .true.

      else

         judge_inout = .false.

      endif


      end ! function judge_inout
!
!______________________________________________________________________
!
      subroutine make_grid( TEST_CASE )
!----------------------------------------------------------------------
!  Generate vertical and horizontal grid, topography, areas and masks.
!______________________________________________________________________
!
      use glob_const , only: DEG2RAD, rk
      use glob_domain
      use grid
      use glob_ocean , only: d, dt, el, et

      implicit none

      integer, intent(in) :: TEST_CASE

      integer i,j,k

      if ( is_master ) then
        call msg_print("IDEALIZED CASE: "//char(40+TEST_CASE),1,"")
      end if 

! generate grid
      do k = 1,kb
        z(k)  = -real(k-1)/real(kb-1)
      end do
      do k = 1,kb-1
        zz(k) = .5*(z(k+1)+z(k))
      end do
      zz(kb) = 2.*zz(kb-1)-zz(kb-2)

      do k=1,kb-1
        dz(k) = z(k)- z(k+1)
        dzz(k)=zz(k)-zz(k+1)
      end do
      dz(kb) = dz(kb-1)
      dzz(kb)=dzz(kb-1)


      do j = 1,jm
        do i = 1,im
          east_e(i,j) = 100.
     &        + 180.*real(i_global(i)-1)/real(im_global-1)
          north_e(i,j)= -80.
     &        + 160.*real(j_global(j)-1)/real(jm_global-1)
        end do
      end do

      do j = 1,jm
        do i = 2,im
          east_u(i,j) = .5*( east_e(i,j) + east_e(i-1,j) )
          north_u(i,j)= .5*( north_e(i,j)+ north_e(i-1,j))
        end do
      end do
      call exchange2d_mpi( east_u, im,jm)
      call exchange2d_mpi(north_u, im,jm)
      if (n_west==-1) then
        east_u(1,:) = 2.*east_u(2,:) - east_u(3,:)
        north_u(1,:)= 2.*north_u(2,:)- north_u(3,:)
      end if
      do j = 2,jm
        do i = 1,im
          east_v(i,j) = .5*( east_e(i,j) + east_e(i,j-1) )
          north_v(i,j)= .5*( north_e(i,j)+ north_e(i,j-1))
        end do
      end do
      call exchange2d_mpi( east_v, im,jm)
      call exchange2d_mpi(north_v, im,jm)
      if (n_south==-1) then
      east_v(:,1) = 2.*east_v(:,2) - east_v(:,3)
      north_v(:,1)= 2.*north_v(:,2)- north_v(:,3)
      end if

      do j = 1,jm
        do i = 2,im
          dx(i,j) = 111800.*cos(north_e(i,j)*DEG2RAD)
     &                     *abs(east_u(i,j)-east_u(i-1,j))
        end do
      end do
      do j = 2,jm
        do i = 1,im
          dy(i,j) = 111800.*cos(north_e(i,j)*DEG2RAD)
     &                     *abs(north_v(i,j)-north_v(i,j-1))
        end do
      end do
      call exchange2d_mpi(dx, im,jm)
      call exchange2d_mpi(dy, im,jm)
      if (n_west==-1) then
        dx(1,:) = dx(2,:)
        dy(1,:) = dy(2,:)
      end if
      if (n_south==-1) then
        dx(:,1) = dx(:,2)
        dy(:,1) = dy(:,2)
      end if

      rot = 0.
      do j = 1,jm
        do i = 1,im
          h(i,j)=100.-1000.*(tanh(5.*(real(i_global(i))
     &                               -real(im_global)/6.)
     &                              /(real(im_global)/6.
     &                               -real(im_global)))
     &                      +tanh(5.*(real(j_global(j))
     &                               -real(jm_global)/6.)
     &                              /(real(jm_global)/6.
     &                               -real(jm_global))))
          !h(i,j)=h(i,j)*(1.-real(i_global(i))/real(im_global))
          h(i,j)=h(i,j)+1000.*(tanh(5.*(real(im_global
     &                                      -i_global(i))
     &                               +real(im_global)/6.)
     &                              /(real(im_global)/6.
     &                               +real(im_global)))
     &                      +tanh(5.*(real(jm_global-
     &                                     j_global(j))
     &                               +real(jm_global)/6.)
     &                              /(real(jm_global)/6.
     &                               +real(jm_global))))
        end do
      end do

      fsm = 1.
      if (n_west==-1) fsm( 1,:) = 0.
      if (n_east==-1) fsm(im,:) = 0.
      if (n_south==-1) fsm(:, 1) = 0.
      if (n_north==-1) fsm(:,jm) = 0.

      dvm = fsm
      dum = fsm
      do j=1,jm-1
        do i=1,im
          if (fsm(i,j)==0..and.fsm(i,j+1)/=0.) dvm(i,j+1) = 0.
        end do
      end do
      do j=1,jm
        do i=1,im-1
          if (fsm(i,j)==0..and.fsm(i+1,j)/=0.) dum(i+1,j) = 0.
        end do
      end do
      call exchange2d_mpi(dum,im,jm)
      call exchange2d_mpi(dvm,im,jm)
!     The followings are read in read_grid_pnetcdf:
!     z,zz,dx,dy
!     east_u,east_v,east_e,east_c
!     north_u,north_v,north_e,north_c
!     rot,h,fsm,dum,dvm


      call print_grid_vert

! set up Coriolis parameter
        do j=1,jm
          do i=1,im
            cor(i,j)=2.*7.29e-5*sin(north_e(i,j)*DEG2RAD)
          end do
        end do

! inertial period for temporal filter
!      period=(2.e0*pi)/abs(cor(im/2,jm/2))/86400.e0

! calculate areas of "t" and "s" cells
      art = dx*dy

! calculate areas of "u" and "v" cells
      do j=2,jm
        do i=2,im
          aru(i,j)=.25*(dx(i,j)+dx(i-1,j))*(dy(i,j)+dy(i-1,j))
          arv(i,j)=.25*(dx(i,j)+dx(i,j-1))*(dy(i,j)+dy(i,j-1))
        end do
      end do
      call exchange2d_mpi(aru,im,jm)
      call exchange2d_mpi(arv,im,jm)

      if (n_west.eq.-1) then
        aru(1,:)=aru(2,:)
        arv(1,:)=arv(2,:)
      end if

      if (n_south.eq.-1) then
        aru(:,1)=aru(:,2)
        arv(:,1)=arv(:,2)
      end if

      d  = h + el
      dt = h + et


      end
!______________________________________________________________________
!
      subroutine print_grid_vert
!----------------------------------------------------------------------
!  Prints vertical sigma distribution.
!______________________________________________________________________

      use glob_domain, only: is_master, kb
      use grid       , only: z, zz, dz, dzz

      implicit none

      integer k

! print vertical grid information
      if ( is_master ) then
        print '(/a)', "=========== VERTICAL DISCRETIZATION ==========="
        print '(/3x,a,9x,a,8x,a,8x,a,8x,a)', 'k','z','zz','dz','dzz'
        do k=1,kb
          print '(1x,i5,4f10.3)', k,z(k),zz(k),dz(k),dzz(k)
        end do
      end if

      end ! subroutine

!______________________________________________________________________
! fhx:interp_flag
! fhx:fgrid
      subroutine check_interpolation
!----------------------------------------------------------------------
!  Decides on wether to use interpolation or not.
!______________________________________________________________________

      use config     , only: calc_interp
      use glob_domain, only: im_global, im_global_coarse, is_master
     &                     , jm_global, jm_global_coarse
     &                     , x_division, y_division

      implicit none

      if (  ((im_global_coarse-2)*x_division == (im_global-2))
     $ .and.((jm_global_coarse-2)*y_division == (jm_global-2)) ) then

        calc_interp = .TRUE.

      else if ( (im_global_coarse == im_global)
     $     .and.(jm_global_coarse == jm_global) ) then

        calc_interp = .FALSE.

      else

        if ( is_master ) print '(a//a)'
     &    , 'Incompatible number of *_global and *_global_coarse'
     &    , 'POM terminated with error'
        stop

      end if


      end ! subroutine

!______________________________________________________________________
!
      subroutine check_cfl_min
!----------------------------------------------------------------------
!  Estimates minimum barotropic timestep (sec)
!______________________________________________________________________
!
        use glob_const , only: grav, rk, SMALL
        use glob_domain, only: im, jm, master_task, my_task, POM_COMM
        use grid       , only: dx, dy, fsm, h
        use model_run  , only: dte
        use mpi        , only: MPI_DOUBLE, MPI_MIN, MPI_REAL
     &                       , mpi_reduce

        implicit none

        real(rk), dimension(im,jm) :: cfl
        real(rk)                      cflmin, tmp
        integer                       ierr, MPI_RK
        character(len=128)            msg, desc


        if ( rk == 8 ) then
          MPI_RK = MPI_DOUBLE
        else
          MPI_RK = MPI_REAL
        end if

        cfl = 0.
        cflmin = 1.d10

        cfl = fsm*.5 / sqrt(1./dx**2+1./dy**2) / sqrt(grav*(h+SMALL))

        cflmin = minval(cfl, cfl>0.)

        call mpi_reduce( cflmin, tmp, 1, MPI_RK, mpi_min
     &                 , master_task, POM_COMM, ierr     )

        if ( my_task == master_task ) then
          if ( cflmin < dte ) then
            write(msg, '(a,f5.2,a)')
     &                 "Timestep (", dte, ") is too large!"
            write(desc, '(a,f5.2,a)')
     &              "You are strongly advised to make dte smaller than "
     &                                                       ,cflmin,"."
            call msg_print(msg, 2, desc)
          end if
        end if


      end ! subroutine check_cfl_min
!______________________________________________________________________
!
      subroutine modules_initial_step( dtime )

        use air        , only: air_init => init
        use bry        , only: bry_init => init
        use clim       , only: clm_init => init
        use seaice     , only: ice_init => init
        use tide       , only: tide_init=> init
        use module_time

        implicit none

        type(date), intent(in) :: dtime


        call clm_init( dtime )
        call air_init( dtime )
        call bry_init( dtime )
        call ice_init( dtime )
        call tide_init(dtime )


      end ! subroutine
!!
!!______________________________________________________________________
!!
!      subroutine write_output_init_pnetcdf( out_file )
!!----------------------------------------------------------------------
!!  Write initial state output file.
!!______________________________________________________________________
!!
!      use air        , only: wssurf, wtsurf, wusurf, wvsurf
!      use bry
!      use config     , only: mode, title, use_air
!      use glob_domain
!      use grid
!      use glob_misc
!      use glob_ocean
!      use model_run  , only: time, time_start
!      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
!      use pnetcdf    , only: nf90mpi_close     , nf90mpi_create     &
!                           , nf90mpi_def_dim   , nf90mpi_enddef     &
!                           , nf90mpi_get_att   , nf90mpi_inq_varid  &
!                           , nf90mpi_open      , nf90mpi_put_att    &
!                           , nfmpi_put_vara_all                     &
!                           , nfmpi_put_vara_real_all                &
!                           , NF_64BIT_OFFSET, NF_CLOBBER            &
!                           , NF_FLOAT       , NF_GLOBAL             &
!                           , NF_NOERR       , NF_NOWRITE            &
!                           , NF_UNLIMITED   , NF_WRITE
!
!      implicit none
!
!      character(len=*), intent(in) :: out_file  ! Output filename
!
!      integer time_dimid, x_dimid, y_dimid, z_dimid
!      integer  aamfrz_varid,    dum_varid,    dvm_varid,     dx_varid &
!            ,      dy_varid, east_c_varid, east_e_varid, east_u_varid &
!            ,  east_v_varid,    ele_varid,    eln_varid,    els_varid &
!            ,     elw_varid,    frz_varid,    fsm_varid,      h_varid &
!            ,      kh_varid,     km_varid,north_c_varid,north_e_varid &
!            , north_u_varid,north_v_varid,    rho_varid,    rot_varid &
!            ,       s_varid,      t_varid,   time_varid, wssurf_varid &
!            ,  wtsurf_varid, wusurf_varid, wvsurf_varid,      z_varid &
!            ,      zz_varid
!
!      integer, dimension(4)      :: vdims
!      integer(MPI_OFFSET_KIND)                       &
!             , dimension(4)      :: start(4),edge(4)
!
!      integer ncid
!
!      real(kind=4), dimension(im      ) :: out1x
!      real(kind=4), dimension(   jm   ) :: out1y
!      real(kind=4), dimension(      kb) :: out1z
!      real(kind=4), dimension(im,jm   ) :: out2
!      real(kind=4), dimension(im,jm,kb) :: out3
!
!
!      out1x = 0.
!      out1y = 0.
!      out1z = 0.
!      out2 = 0.
!      out3 = 0.
!
!      call msg_print("", 6, "Writing file `"//trim(out_file)//"`")
!
!      call check( nf90mpi_create                  &
!                    ( POM_COMM, trim(out_file)    &
!                    , NF_CLOBBER+NF_64BIT_OFFSET  &
!                    , MPI_INFO_NULL, ncid )       &
!                , 'nf_create: '//out_file )
!
!! define global attributes
!      call check( nf90mpi_put_att                            &
!                  ( ncid, NF_GLOBAL, 'title', trim(title) )  &
!                , 'nf_put_att: title @ '//out_file )
!
!      call check( nf90mpi_put_att                                  &
!                  ( ncid, NF_GLOBAL, 'description', 'Init file' )  &
!                , 'nf_put_att: description @ '//out_file )
!
!! define dimensions
!      call check( nf90mpi_def_dim                                   &
!                  ( ncid, 'time', int(NF_UNLIMITED,8), time_dimid ) &
!                , 'nf_def_dim: time @ '//out_file )
!      call check( nf90mpi_def_dim                                   &
!                  ( ncid,    'z', int(          kb,8),    z_dimid ) &
!                , 'nf_def_dim: z @ '//out_file )
!      call check( nf90mpi_def_dim                                   &
!                  ( ncid,    'y', int(   jm_global,8),    y_dimid ) &
!                , 'nf_def_dim: y @ '//out_file )
!      call check( nf90mpi_def_dim                                   &
!                  ( ncid,    'x', int(   im_global,8),    x_dimid ) &
!                , 'nf_def_dim: x @ '//out_file )
!
!! define variables and their attributes
!      vdims(1) = time_dimid
!      time_varid = def_var_pnetcdf(ncid,'time',1,vdims,NF_FLOAT  &
!                              ,'time','days since '//time_start  &
!                              ,-1,0.,' ',.false.)
!
!! Skip baroclinic output if model run is barotropic
!      if ( mode /= 2 ) then
!
!        vdims(1) = z_dimid
!        z_varid = def_var_pnetcdf(ncid,'sigma',1,vdims,NF_FLOAT  &
!                            ,'sigma of cell face','sigma_level'  &
!                            ,-1,0.,' ',.false.)
!        call check( nf90mpi_put_att( ncid, z_varid               &
!                                   , 'standard_name'             &
!                                   , 'ocean_sigma_coordinate' )  &
!                  , 'nf_put_att: stdname @ z @ '//out_file )
!        call check( nf90mpi_put_att( ncid, z_varid                  &
!                                   , 'formula_terms'                &
!                                   , 'sigma: z eta: elb depth: h' ) &
!                  , 'nf_put_att: formterms @ z @ '//out_file )
!
!        zz_varid = def_var_pnetcdf(ncid,'zz',1,vdims,NF_FLOAT  &
!                        ,'sigma of cell centre','sigma_level'  &
!                        ,-1,0.,' ',.false.)
!        call check( nf90mpi_put_att( ncid, zz_varid              &
!                                   , 'standard_name'             &
!                                   , 'ocean_sigma_coordinate' )  &
!                  , 'nf_put_att: stdname @ zz @ '//out_file )
!        call check( nf90mpi_put_att( ncid, zz_varid                  &
!                                   , 'formula_terms'                 &
!                                   , 'sigma: zz eta: elb depth: h' ) &
!                  , 'nf_put_att: formterms @ zz @ '//out_file )
!
!      end if
!
!
!      vdims(1) = x_dimid
!      els_varid = def_var_pnetcdf(ncid,'els',1,vdims,NF_FLOAT  &
!                                ,'South elevation BC','meter'  &
!                                ,-1,0.,'east_e',.true.)
!      eln_varid = def_var_pnetcdf(ncid,'eln',1,vdims,NF_FLOAT  &
!                                ,'North elevation BC','meter'  &
!                                ,-1,0.,'east_e',.true.)
!      vdims(1) = y_dimid
!      ele_varid = def_var_pnetcdf(ncid,'ele',1,vdims,NF_FLOAT  &
!                                 ,'East elevation BC','meter'  &
!                                 ,-1,0.,'north_e',.true.)
!      elw_varid = def_var_pnetcdf(ncid,'elw',1,vdims,NF_FLOAT  &
!                                 ,'West elevation BC','meter'  &
!                                 ,-1,0.,'north_e',.true.)
!
!
!      vdims(1) = x_dimid
!      vdims(2) = y_dimid
!      dx_varid = def_var_pnetcdf(ncid,'dx',2,vdims,NF_FLOAT  &
!                             ,'grid increment in x','meter'  &
!                             ,-1,0.,'east_e north_e',.true.)
!      dy_varid = def_var_pnetcdf(ncid,'dy',2,vdims,NF_FLOAT  &
!                             ,'grid increment in y','meter'  &
!                             ,-1,0.,'east_e north_e',.true.)
!! TODO: Make units of easting and northing flexible (degree/meter)
!      east_u_varid = def_var_pnetcdf(ncid,'east_u',2,vdims,NF_FLOAT  &
!                                    ,'easting of u-points','degree'  &
!                                    ,-1,0.,'east_u north_u',.true.)
!      east_v_varid = def_var_pnetcdf(ncid,'east_v',2,vdims,NF_FLOAT  &
!                                    ,'easting of v-points','degree'  &
!                                    ,-1,0.,'east_v north_v',.true.)
!      east_e_varid = def_var_pnetcdf(ncid,'east_e',2,vdims,NF_FLOAT  &
!                            ,'easting of elevation points','degree'  &
!                            ,-1,0.,'east_e north_e',.true.)
!      east_c_varid = def_var_pnetcdf(ncid,'east_c',2,vdims,NF_FLOAT  &
!                                ,'easting of cell corners','degree'  &
!                                ,-1,0.,'east_c north_c',.true.)
!      north_u_varid= def_var_pnetcdf(ncid,'north_u',2,vdims,NF_FLOAT &
!                                    ,'northing of u-points','degree' &
!                                    ,-1,0.,'east_u north_u',.true.)
!      north_v_varid= def_var_pnetcdf(ncid,'north_v',2,vdims,NF_FLOAT &
!                                    ,'northing of v-points','degree' &
!                                    ,-1,0.,'east_v north_v',.true.)
!      north_e_varid= def_var_pnetcdf(ncid,'north_e',2,vdims,NF_FLOAT &
!                            ,'northing of elevation points','degree' &
!                            ,-1,0.,'east_e north_e',.true.)
!      north_c_varid= def_var_pnetcdf(ncid,'north_c',2,vdims,NF_FLOAT &
!                                ,'northing of cell corners','degree' &
!                                ,-1,0.,'east_c north_c',.true.)
!      rot_varid = def_var_pnetcdf(ncid,'rot',2,vdims,NF_FLOAT  &
!               ,'Rotation angle of x-axis wrt. east','degree'  &
!               ,-1,0.,'east_e north_e',.true.)
!      h_varid = def_var_pnetcdf(ncid,'h',2,vdims,NF_FLOAT  &
!                       ,'undisturbed water depth','metre'  &
!                       ,-1,0.,'east_e north_e',.true.)
!      fsm_varid = def_var_pnetcdf(ncid,'fsm',2,vdims,NF_FLOAT  &
!                         ,'free surface mask','dimensionless'  &
!                         ,-1,0.,'east_e north_e',.true.)
!      dum_varid = def_var_pnetcdf(ncid,'dum',2,vdims,NF_FLOAT  &
!                           ,'u-velocity mask','dimensionless'  &
!                           ,-1,0.,'east_u north_u',.true.)
!      dvm_varid = def_var_pnetcdf(ncid,'dvm',2,vdims,NF_FLOAT  &
!                           ,'v-velocity mask','dimensionless'  &
!                           ,-1,0.,'east_v north_v',.true.)
!      frz_varid = def_var_pnetcdf(ncid,'frz',2,vdims,NF_FLOAT  &
!                                ,'frz coeffi','dimensionless'  &
!                                ,-1,0.,'east_v north_v',.true.)
!      aamfrz_varid = def_var_pnetcdf(ncid,'aamfrz',2,vdims,NF_FLOAT  &
!                             ,'aamfrz coefficients','dimensionless'  &
!                             ,-1,0.,'east_v north_v',.true.)
!
!! If atmospheric forcing is enabled...
!      if ( use_air ) then
!        vdims(1) = x_dimid
!        vdims(2) = y_dimid
!        vdims(3) = time_dimid
!        wusurf_varid=def_var_pnetcdf(ncid,'wusurf',3,vdims,NF_FLOAT  &
!                                 ,'x-momentum flux','metre^2/sec^2'  &
!                                 ,-1,0.,'east_u north_u',.true.)
!        wvsurf_varid=def_var_pnetcdf(ncid,'wvsurf',3,vdims,NF_FLOAT  &
!                                 ,'y-momentum flux','metre^2/sec^2'  &
!                                 ,-1,0.,'east_v north_v',.true.)
!        wtsurf_varid=def_var_pnetcdf(ncid,'wtsurf',3,vdims,NF_FLOAT  &
!                                     ,'temperature flux','deg m/s'   &
!                                     ,-1,0.,'east_e north_e',.true.)
!        wssurf_varid=def_var_pnetcdf(ncid,'wssurf',3,vdims,NF_FLOAT  &
!                                     ,'salinity flux','psu m/s'      &
!                                     ,-1,0.,'east_e north_e',.true.)
!      end if
!
!! Skip baroclinic output...
!      if ( mode /= 2 ) then
!
!        vdims(1) = x_dimid
!        vdims(2) = y_dimid
!        vdims(3) = z_dimid
!        vdims(4) = time_dimid
!        t_varid = def_var_pnetcdf(ncid,'t',4,vdims,NF_FLOAT  &
!                          ,'potential temperature','degC'    &
!                          ,-1,0.,'east_e north_e zz',.true.)
!        s_varid = def_var_pnetcdf(ncid,'s',4,vdims,NF_FLOAT  &
!                          ,'salinity x rho / rhoref','PSS'   &
!                          ,-1,0.,'east_e north_e zz',.true.)
!        rho_varid = def_var_pnetcdf(ncid,'rho',4,vdims,NF_FLOAT  &
!                       ,'(density-1000)/rhoref','dimensionless'  &
!                       ,-1,0.,'east_e north_e zz',.true.)
!        kh_varid = def_var_pnetcdf(ncid,'kh',4,vdims,NF_FLOAT  &
!                         ,'vertical diffusivity','metre2/sec'  &
!                         ,-1,0.,'east_e north_e zz',.true.)
!        km_varid = def_var_pnetcdf(ncid,'km',4,vdims,NF_FLOAT  &
!                           ,'vertical viscosity','metre2/sec'  &
!                           ,-1,0.,'east_e north_e zz',.true.)
!
!      end if
!
!! end definitions
!      call check( nf90mpi_enddef(ncid)                     &
!                , 'nf_enddef: init file @ '//out_file )
!
!
!! write data
!      start(1) = 1
!      edge(1) = 1
!      out1z = real(time,4)
!      call check( nfmpi_put_vara_real_all                  &
!                  ( ncid,time_varid,start,edge,out1z )     &
!                , 'nf_put_vara_real:time @ '//out_file )
!
!! Skip baroclinic output...
!      if ( mode /= 2 ) then
!        start(1) = 1
!        edge(1) = kb
!        out1z = real(z,4)
!        call check( nfmpi_put_vara_real_all               &
!                    ( ncid,z_varid,start,edge,out1z )      &
!                  , 'nf_put_var_real: z @ '//out_file )
!        out1z = real(zz,4)
!        call check( nfmpi_put_vara_real_all                &
!                    ( ncid,zz_varid,start,edge,out1z )      &
!                  , 'nf_put_var_real: zz @ '//out_file )
!      end if
!
!! East
!      if ( n_east == -1 ) then
!        start(1) = j_global(1)
!        edge(1) = jm
!        out1y = real(el_bry%est(1,:),4)
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,ele_varid,start,edge,out1y )  &
!                  , 'nf_put_var_real: ele @ '//out_file )
!      else
!        start(1) = 1
!        edge(1) = 0
!        out1y = 0.
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,ele_varid,start,edge,out1y )  &
!                  , 'nf_put_var_real: ele (dummy) @ '//out_file )
!      end if
!! West
!      if ( n_west == -1 ) then
!        start(1) = j_global(1)
!        edge(1) = jm
!        out1y = real(el_bry%wst(1,:),4)
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,elw_varid,start,edge,out1y )  &
!                  , 'nf_put_var_real: elw @ '//out_file )
!      else
!        start(1) = 1
!        edge(1) = 0
!        out1y = 0.
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,elw_varid,start,edge,out1y )  &
!                  , 'nf_put_var_real: elw (dummy) @ '//out_file )
!      end if
!! South
!      if ( n_south == -1 ) then
!        start(1) = i_global(1)
!        edge(1) = im
!        out1x = real(el_bry%sth(:,1),4)
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,els_varid,start,edge,out1x )  &
!                  , 'nf_put_var_real: els @ '//out_file )
!      else
!        start(1) = 1
!        edge(1) = 0
!        out1x = 0.
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,els_varid,start,edge,out1x )  &
!                  , 'nf_put_var_real: els (dummy) @ '//out_file )
!      end if
!! North
!      if ( n_north == -1 ) then
!        start(1) = i_global(1)
!        edge(1) = im
!        out1x = real(el_bry%nth(:,1),4)
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,eln_varid,start,edge,out1x )  &
!                  , 'nf_put_var_real: eln @ '//out_file )
!      else
!        start(1) = 1
!        edge(1) = 0
!        out1x = 0.
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,eln_varid,start,edge,out1x )  &
!                  , 'nf_put_var_real: eln @ (dummy)'//out_file )
!      end if
!
!
!      start(1) = i_global(1)
!      start(2) = j_global(1)
!      edge(1) = im
!      edge(2) = jm
!      out2 = real(dx,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,dx_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: dx @ '//out_file )
!      out2 = real(dy,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,dy_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: dy @ '//out_file )
!      out2 = real(east_u,4)
!      call check( nfmpi_put_vara_real_all                    &
!                  ( ncid,east_u_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: east_u @ '//out_file )
!      out2 = real(east_v,4)
!      call check( nfmpi_put_vara_real_all                    &
!                  ( ncid,east_v_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: east_v @ '//out_file )
!      out2 = real(east_e,4)
!      call check( nfmpi_put_vara_real_all                    &
!                  ( ncid,east_e_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: east_e @ '//out_file )
!      out2 = real(east_c,4)
!      call check( nfmpi_put_vara_real_all                    &
!                  ( ncid,east_c_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: east_c @ '//out_file )
!      out2 = real(north_u,4)
!      call check( nfmpi_put_vara_real_all                     &
!                  ( ncid,north_u_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: north_u @ '//out_file )
!      out2 = real(north_v,4)
!      call check( nfmpi_put_vara_real_all                     &
!                  ( ncid,north_v_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: north_v @ '//out_file )
!      out2 = real(north_e,4)
!      call check( nfmpi_put_vara_real_all                     &
!                  ( ncid,north_e_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: north_e @ '//out_file )
!      out2 = real(north_c,4)
!      call check( nfmpi_put_vara_real_all                     &
!                  ( ncid,north_c_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: north_c @ '//out_file )
!      out2 = real(rot,4)
!      call check( nfmpi_put_vara_real_all                 &
!                  ( ncid,rot_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: rot @ '//out_file )
!      out2 = real(h,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,h_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: h @ '//out_file )
!      out2 = real(fsm,4)
!      call check( nfmpi_put_vara_real_all                 &
!                  ( ncid,fsm_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: fsm @ '//out_file )
!      out2 = real(dum,4)
!      call check( nfmpi_put_vara_real_all                 &
!                  ( ncid,dum_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: dum @ '//out_file )
!      out2 = real(dvm,4)
!      call check( nfmpi_put_vara_real_all                 &
!                  ( ncid,dvm_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: dvm @ '//out_file )
!
!      if ( USE_SPONGE ) then
!        out2 = real(frz,4)
!        call check( nfmpi_put_vara_real_all                 &
!                    ( ncid,frz_varid,start,edge,out2 )      &
!                  , 'nf_put_var_real: frz @ '//out_file )
!        out2 = real(aamfrz,4)
!        call check( nfmpi_put_vara_real_all                    &
!                    ( ncid,aamfrz_varid,start,edge,out2 )      &
!                  , 'nf_put_var_real: aamfrz @ '//out_file )
!      end if
!
!! Write only if atmospheric forcing is enabled.
!      if ( use_air ) then
!
!        start(1) = i_global(1)
!        start(2) = j_global(1)
!        start(3) = 1
!        edge(1) = im
!        edge(2) = jm
!        edge(3) = 1
!        out2 = real(wusurf,4)
!        call check( nfmpi_put_vara_real_all                     &
!                    ( ncid,wusurf_varid,start,edge,out2 )       &
!                  , 'nf_put_vara_real: wusurf @ '//out_file )
!        out2 = real(wvsurf,4)
!        call check( nfmpi_put_vara_real_all                     &
!                    ( ncid,wvsurf_varid,start,edge,out2 )       &
!                  , 'nf_put_vara_real: wvsurf @ '//out_file )
!        out2 = real(wtsurf,4)
!        call check( nfmpi_put_vara_real_all                     &
!                    ( ncid,wtsurf_varid,start,edge,out2 )       &
!                  , 'nf_put_vara_real: wtsurf @ '//out_file )
!        out2 = real(wssurf,4)
!        call check( nfmpi_put_vara_real_all                     &
!                    ( ncid,wssurf_varid,start,edge,out2 )       &
!                  , 'nf_put_vara_real: wssurf @ '//out_file )
!
!      end if
!
!! Skip baroclinic output...
!      if ( mode /= 2 ) then
!
!        start(1) = i_global(1)
!        start(2) = j_global(1)
!        start(3) = 1
!        start(4) = 1
!        edge(1) = im
!        edge(2) = jm
!        edge(3) = kb
!        edge(4) = 1
!        out3 = real(tb,4)
!        call check( nfmpi_put_vara_real_all               &
!                    ( ncid,t_varid,start,edge,out3 )      &
!                  , 'nf_put_vara_real: t @ '//out_file )
!        out3 = real(sb,4)
!        call check( nfmpi_put_vara_real_all               &
!                    ( ncid,s_varid,start,edge,out3 )      &
!                  , 'nf_put_vara_real: s @ '//out_file )
!        out3 = real(rho,4)
!        call check( nfmpi_put_vara_real_all                 &
!                    ( ncid,rho_varid,start,edge,out3 )      &
!                  , 'nf_put_vara_real: rho @ '//out_file )
!        out3 = real(kh,4)
!        call check( nfmpi_put_vara_real_all                &
!                    ( ncid,kh_varid,start,edge,out3 )      &
!                  , 'nf_put_vara_real: kh @ '//out_file )
!        out3 = real(km,4)
!        call check( nfmpi_put_vara_real_all                &
!                    ( ncid,km_varid,start,edge,out3 )      &
!                  , 'nf_put_vara_real: km @ '//out_file )
!
!      end if
!
!! close file
!      call check( nf90mpi_close(ncid)                   &
!                , 'nf_close: output:  @ '//out_file )
!
!
!    end ! subroutine write_output_init_pnetcdf
!!
!!______________________________________________________________________
!!
!    subroutine write_debug_pnetcdf( out_file )
!!----------------------------------------------------------------------
!!  Write initial state output file.
!!______________________________________________________________________
!!
!!      use air        , only: wssurf, wtsurf, wusurf, wvsurf
!      use bry        , only: el_bry
!      use config     , only: mode, title, use_air
!      use glob_domain
!      use grid
!      use glob_misc
!      use glob_ocean
!      use tide
!      use model_run  , only: time, time_start
!      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
!      use pnetcdf    , only: nf90mpi_close     , nf90mpi_create     &
!                           , nf90mpi_def_dim   , nf90mpi_enddef     &
!                           , nf90mpi_get_att   , nf90mpi_inq_varid  &
!                           , nf90mpi_open      , nf90mpi_put_att    &
!                           , nfmpi_put_vara_all                     &
!                           , nfmpi_put_vara_real_all                &
!                           , NF_64BIT_OFFSET, NF_CLOBBER            &
!                           , NF_FLOAT       , NF_GLOBAL             &
!                           , NF_NOERR       , NF_NOWRITE            &
!                           , NF_UNLIMITED   , NF_WRITE
!
!      implicit none
!
!      character(len=*), intent(in) :: out_file  ! Output filename
!
!      integer time_dimid, x_dimid, y_dimid, z_dimid, cons_dimid
!      integer  aam2d_varid,    aam_varid,  advua_varid, advva_varid  &
!            ,   advx_varid,   advy_varid,  adx2d_varid, ady2d_varid  &
!            ,  drhox_varid,  drhoy_varid,  drx2d_varid, dry2d_varid  &
!            ,e_atmos_varid,    egb_varid,    egf_varid,    el_varid  &
!            ,    elb_varid,    ele_varid,    eln_varid,   elf_varid  &
!            ,    els_varid,    elw_varid,     et_varid,   etb_varid  &
!            ,    etf_varid, fluxua_varid, fluxva_varid,    kh_varid  &
!            ,     km_varid,     kq_varid,      l_varid,    q2_varid  &
!            ,    q2b_varid,    q2l_varid,   q2lb_varid,   rho_varid  &
!            ,      s_varid,     sb_varid,      t_varid,  time_varid  &
!            ,     tb_varid,      u_varid,     ua_varid,   uab_varid  &
!            ,    uaf_varid,     ub_varid,     uf_varid,     v_varid  &
!            ,     va_varid,     vb_varid,    vab_varid,   vaf_varid  &
!            ,     vf_varid,      w_varid,     wr_varid,     z_varid  &
!            ,     zz_varid, el_amp_varid
!
!      integer, dimension(4)      :: vdims
!      integer(MPI_OFFSET_KIND)                       &
!             , dimension(4)      :: start(4),edge(4)
!
!      integer ncid
!
!      real(kind=4), dimension(im      ) :: out1x
!      real(kind=4), dimension(   jm   ) :: out1y
!      real(kind=4), dimension(      kb) :: out1z
!      real(kind=4), dimension(im,jm   ) :: out2
!      real(kind=4), dimension(im,jm,kb) :: out3
!
!
!      out1x = 0.
!      out1y = 0.
!      out1z = 0.
!      out2 = 0.
!      out3 = 0.
!
!      call msg_print("", 6, "Writing debug `"//trim(out_file)//"`")
!
!      call check( nf90mpi_create                  &
!                    ( POM_COMM, trim(out_file)    &
!                    , NF_CLOBBER+NF_64BIT_OFFSET  &
!                    , MPI_INFO_NULL, ncid )       &
!                , 'nf_create: '//out_file )
!
!! define global attributes
!      call check( nf90mpi_put_att                            &
!                  ( ncid, NF_GLOBAL, 'title', trim(title) )  &
!                , 'nf_put_att: title @ '//out_file )
!
!      call check( nf90mpi_put_att                                   &
!                  ( ncid, NF_GLOBAL, 'description', 'Debug file' )  &
!                , 'nf_put_att: description @ '//out_file )
!
!! define dimensions
!      call check( nf90mpi_def_dim                                   &
!                  ( ncid, 'time', int(NF_UNLIMITED,8), time_dimid ) &
!                , 'nf_def_dim: time @ '//out_file )
!      call check( nf90mpi_def_dim                                   &
!                  ( ncid,  'nct', int(       ncons,8), cons_dimid ) &
!                , 'nf_def_dim: constituents @ '//out_file )
!      call check( nf90mpi_def_dim                                   &
!                  ( ncid,    'z', int(          kb,8),    z_dimid ) &
!                , 'nf_def_dim: z @ '//out_file )
!      call check( nf90mpi_def_dim                                   &
!                  ( ncid,    'y', int(   jm_global,8),    y_dimid ) &
!                , 'nf_def_dim: y @ '//out_file )
!      call check( nf90mpi_def_dim                                   &
!                  ( ncid,    'x', int(   im_global,8),    x_dimid ) &
!                , 'nf_def_dim: x @ '//out_file )
!
!! define variables and their attributes
!      vdims(1) = time_dimid
!      time_varid = def_var_pnetcdf(ncid,'time',1,vdims,NF_FLOAT  &
!                              ,'time','days since '//time_start  &
!                              ,-1,0.,' ',.false.)
!
!      vdims(1) = z_dimid
!      z_varid = def_var_pnetcdf(ncid,'sigma',1,vdims,NF_FLOAT  &
!                          ,'sigma of cell face','sigma_level'  &
!                          ,-1,0.,' ',.false.)
!      call check( nf90mpi_put_att( ncid, z_varid               &
!                                 , 'standard_name'             &
!                                 , 'ocean_sigma_coordinate' )  &
!                , 'nf_put_att: stdname @ z @ '//out_file )
!      call check( nf90mpi_put_att( ncid, z_varid                  &
!                                 , 'formula_terms'                &
!                                 , 'sigma: z eta: elb depth: h' ) &
!                , 'nf_put_att: formterms @ z @ '//out_file )
!
!      zz_varid = def_var_pnetcdf(ncid,'zz',1,vdims,NF_FLOAT  &
!                      ,'sigma of cell centre','sigma_level'  &
!                      ,-1,0.,' ',.false.)
!      call check( nf90mpi_put_att( ncid, zz_varid              &
!                                 , 'standard_name'             &
!                                 , 'ocean_sigma_coordinate' )  &
!                , 'nf_put_att: stdname @ zz @ '//out_file )
!      call check( nf90mpi_put_att( ncid, zz_varid                  &
!                                 , 'formula_terms'                 &
!                                 , 'sigma: zz eta: elb depth: h' ) &
!                , 'nf_put_att: formterms @ zz @ '//out_file )
!
!
!      vdims(1) = x_dimid
!      els_varid = def_var_pnetcdf(ncid,'els',1,vdims,NF_FLOAT  &
!                                ,'South elevation BC','meter'  &
!                                ,-1,0.,'east_e',.true.)
!      eln_varid = def_var_pnetcdf(ncid,'eln',1,vdims,NF_FLOAT  &
!                                ,'North elevation BC','meter'  &
!                                ,-1,0.,'east_e',.true.)
!      vdims(1) = y_dimid
!      ele_varid = def_var_pnetcdf(ncid,'ele',1,vdims,NF_FLOAT  &
!                                 ,'East elevation BC','meter'  &
!                                 ,-1,0.,'north_e',.true.)
!      elw_varid = def_var_pnetcdf(ncid,'elw',1,vdims,NF_FLOAT  &
!                                 ,'West elevation BC','meter'  &
!                                 ,-1,0.,'north_e',.true.)
!
!
!      vdims(1) = x_dimid
!      vdims(2) = y_dimid
!      aam2d_varid = def_var_pnetcdf(ncid,'aam2d',2,vdims,NF_FLOAT  &
!                             ,'vertical average of aam',''  &
!                             ,-1,0.,'east_e north_e',.true.)
!      advua_varid = def_var_pnetcdf(ncid,'advua',2,vdims,NF_FLOAT      &
!                   ,'sum of the 2nd, 3rd and 4th terms in eq (18)',''  &
!                   ,-1,0.,'east_e north_e',.true.)
!      advva_varid = def_var_pnetcdf(ncid,'advva',2,vdims,NF_FLOAT  &
!                   ,'sum of the 2nd, 3rd and 4th terms in eq (19)',''  &
!                   ,-1,0.,'east_u north_u',.true.)
!      adx2d_varid = def_var_pnetcdf(ncid,'adx2d',2,vdims,NF_FLOAT   &
!                                   ,'vertical integral of advx',''  &
!                                   ,-1,0.,'east_v north_v',.true.)
!      ady2d_varid = def_var_pnetcdf(ncid,'ady2d',2,vdims,NF_FLOAT   &
!                                   ,'vertical integral of advy',''  &
!                                   ,-1,0.,'east_e north_e',.true.)
!      drx2d_varid = def_var_pnetcdf(ncid,'drx2d',2,vdims,NF_FLOAT    &
!                                   ,'vertical integral of drhox',''  &
!                                   ,-1,0.,'east_c north_c',.true.)
!      dry2d_varid = def_var_pnetcdf(ncid,'dry2d',2,vdims,NF_FLOAT   &
!                                   ,'vertical integral of drhoy','' &
!                                   ,-1,0.,'east_u north_u',.true.)
!      e_atmos_varid = def_var_pnetcdf(ncid,'e_atmos',2,vdims,NF_FLOAT &
!                                     ,'atmospheric pressure',''       &
!                                     ,-1,0.,'east_v north_v',.true.)
!      egb_varid = def_var_pnetcdf(ncid,'egb',2,vdims,NF_FLOAT          &
!         ,'surface elevation use for pressure gradient at time n-1','' &
!         ,-1,0.,'east_e north_e',.true.)
!      egf_varid = def_var_pnetcdf(ncid,'egf',2,vdims,NF_FLOAT &
!         ,'surface elevation use for pressure gradient at time n+1','' &
!         ,-1,0.,'east_c north_c',.true.)
!      el_varid = def_var_pnetcdf(ncid,'el',2,vdims,NF_FLOAT  &
!         ,'surface elevation used in the external mode at time n',''  &
!         ,-1,0.,'east_e north_e',.true.)
!      elb_varid = def_var_pnetcdf(ncid,'elb',2,vdims,NF_FLOAT  &
!         ,'surface elevation used in the external mode at time n-1','' &
!         ,-1,0.,'east_e north_e',.true.)
!      elf_varid = def_var_pnetcdf(ncid,'elf',2,vdims,NF_FLOAT  &
!         ,'surface elevation used in the external mode at time n+1','' &
!         ,-1,0.,'east_e north_e',.true.)
!      et_varid = def_var_pnetcdf(ncid,'et',2,vdims,NF_FLOAT  &
!         ,'surface elevation used in the internal mode at time n',''  &
!         ,-1,0.,'east_u north_u',.true.)
!      etb_varid = def_var_pnetcdf(ncid,'etb',2,vdims,NF_FLOAT  &
!         ,'surface elevation used in the internal mode at time n-1','' &
!         ,-1,0.,'east_v north_v',.true.)
!      etf_varid = def_var_pnetcdf(ncid,'etf',2,vdims,NF_FLOAT  &
!         ,'surface elevation used in the internal mode at time n+1','' &
!         ,-1,0.,'east_v north_v',.true.)
!      fluxua_varid = def_var_pnetcdf(ncid,'fluxua',2,vdims,NF_FLOAT  &
!                             ,'fluxua','dimensionless'  &
!                             ,-1,0.,'east_v north_v',.true.)
!      fluxva_varid = def_var_pnetcdf(ncid,'fluxva',2,vdims,NF_FLOAT  &
!                             ,'fluxva','dimensionless'  &
!                             ,-1,0.,'east_v north_v',.true.)
!      ua_varid = def_var_pnetcdf(ncid,'ua',2,vdims,NF_FLOAT         &
!                                ,'vertical mean of u at time n',''  &
!                                ,-1,0.,'east_v north_v',.true.)
!      uab_varid = def_var_pnetcdf(ncid,'uab',2,vdims,NF_FLOAT          &
!                                 ,'vertical mean of u at time n-1',''  &
!                                 ,-1,0.,'east_v north_v',.true.)
!      uaf_varid = def_var_pnetcdf(ncid,'uaf',2,vdims,NF_FLOAT          &
!                                 ,'vertical mean of u at time n+1',''  &
!                                 ,-1,0.,'east_v north_v',.true.)
!      va_varid = def_var_pnetcdf(ncid,'va',2,vdims,NF_FLOAT         &
!                                ,'vertical mean of v at time n',''  &
!                                ,-1,0.,'east_v north_v',.true.)
!      vab_varid = def_var_pnetcdf(ncid,'vab',2,vdims,NF_FLOAT          &
!                                 ,'vertical mean of v at time n-1',''  &
!                                 ,-1,0.,'east_v north_v',.true.)
!      vaf_varid = def_var_pnetcdf(ncid,'vaf',2,vdims,NF_FLOAT          &
!                                 ,'vertical mean of v at time n+1',''  &
!                                 ,-1,0.,'east_v north_v',.true.)
!
!      vdims(1) = x_dimid
!      vdims(2) = y_dimid
!      vdims(3) = cons_dimid
!      el_amp_varid = def_var_pnetcdf(ncid,'el_amp',3,vdims,NF_FLOAT    &
!                                 ,'elevation amplitude',''  &
!                                 ,-1,0.,'east_e north_e nct',.true.)
!
!      vdims(1) = x_dimid
!      vdims(2) = y_dimid
!      vdims(3) = z_dimid
!      vdims(4) = time_dimid
!      aam_varid = def_var_pnetcdf(ncid,'aam',4,vdims,NF_FLOAT   &
!                          ,'horizontal kinematic viscosity',''  &
!                          ,-1,0.,'east_e north_e zz',.true.)
!      advx_varid = def_var_pnetcdf(ncid,'advx',4,vdims,NF_FLOAT     &
!                  ,'x-horizontal advection and diffusion terms',''  &
!                  ,-1,0.,'east_e north_e zz',.true.)
!      advy_varid = def_var_pnetcdf(ncid,'advy',4,vdims,NF_FLOAT     &
!                  ,'y-horizontal advection and diffusion terms',''  &
!                  ,-1,0.,'east_e north_e zz',.true.)
!      drhox_varid = def_var_pnetcdf(ncid,'drhox',4,vdims,NF_FLOAT     &
!               ,'x-component of the internal baroclinic pressure',''  &
!               ,-1,0.,'east_e north_e zz',.true.)
!      drhoy_varid = def_var_pnetcdf(ncid,'drhoy',4,vdims,NF_FLOAT     &
!               ,'y-component of the internal baroclinic pressure',''  &
!               ,-1,0.,'east_e north_e zz',.true.)
!      kq_varid = def_var_pnetcdf(ncid,'kq',4,vdims,NF_FLOAT     &
!                ,'vertical turbulent something',''  &
!                ,-1,0.,'east_e north_e zz',.true.)
!      l_varid = def_var_pnetcdf(ncid,'l',4,vdims,NF_FLOAT     &
!                ,'turbulence length scale',''  &
!                ,-1,0.,'east_e north_e zz',.true.)
!      q2b_varid = def_var_pnetcdf(ncid,'q2b',4,vdims,NF_FLOAT     &
!                ,'twice the turbulent kinetic energy at time n-1',''  &
!                ,-1,0.,'east_e north_e zz',.true.)
!      q2_varid = def_var_pnetcdf(ncid,'q2',4,vdims,NF_FLOAT     &
!                ,'twice the turbulent kinetic energy at time n',''  &
!                ,-1,0.,'east_e north_e zz',.true.)
!      q2lb_varid = def_var_pnetcdf(ncid,'q2lb',4,vdims,NF_FLOAT     &
!                ,'q2 x l at time n-1',''  &
!                ,-1,0.,'east_e north_e zz',.true.)
!      q2l_varid = def_var_pnetcdf(ncid,'q2l',4,vdims,NF_FLOAT     &
!                ,'q2 x l at time n',''  &
!                ,-1,0.,'east_e north_e zz',.true.)
!      t_varid = def_var_pnetcdf(ncid,'t',4,vdims,NF_FLOAT  &
!                        ,'potential temperature','degC'    &
!                        ,-1,0.,'east_e north_e zz',.true.)
!      tb_varid = def_var_pnetcdf(ncid,'tb',4,vdims,NF_FLOAT    &
!                        ,'potential temperature at n-1','degC' &
!                        ,-1,0.,'east_e north_e zz',.true.)
!      s_varid = def_var_pnetcdf(ncid,'s',4,vdims,NF_FLOAT  &
!                        ,'salinity x rho / rhoref','PSS'   &
!                        ,-1,0.,'east_e north_e zz',.true.)
!      sb_varid = def_var_pnetcdf(ncid,'sb',4,vdims,NF_FLOAT      &
!                        ,'salinity x rho / rhoref at n-1','PSS'  &
!                        ,-1,0.,'east_e north_e zz',.true.)
!      rho_varid = def_var_pnetcdf(ncid,'rho',4,vdims,NF_FLOAT  &
!                     ,'(density-1000)/rhoref','dimensionless'  &
!                     ,-1,0.,'east_e north_e zz',.true.)
!      kh_varid = def_var_pnetcdf(ncid,'kh',4,vdims,NF_FLOAT  &
!                       ,'vertical diffusivity','metre2/sec'  &
!                       ,-1,0.,'east_e north_e zz',.true.)
!      km_varid = def_var_pnetcdf(ncid,'km',4,vdims,NF_FLOAT  &
!                         ,'vertical viscosity','metre2/sec'  &
!                         ,-1,0.,'east_e north_e zz',.true.)
!      u_varid = def_var_pnetcdf(ncid,'u',4,vdims,NF_FLOAT         &
!                        ,'horizontal velocity in x at time n',''  &
!                        ,-1,0.,'east_e north_e zz',.true.)
!      ub_varid = def_var_pnetcdf(ncid,'ub',4,vdims,NF_FLOAT         &
!                        ,'horizontal velocity in x at time n-1',''  &
!                        ,-1,0.,'east_e north_e zz',.true.)
!      uf_varid = def_var_pnetcdf(ncid,'uf',4,vdims,NF_FLOAT         &
!                        ,'horizontal velocity in x at time n-1',''  &
!                        ,-1,0.,'east_e north_e zz',.true.)
!      v_varid = def_var_pnetcdf(ncid,'v',4,vdims,NF_FLOAT         &
!                        ,'horizontal velocity in x at time n',''  &
!                        ,-1,0.,'east_e north_e zz',.true.)
!      vb_varid = def_var_pnetcdf(ncid,'vb',4,vdims,NF_FLOAT         &
!                        ,'horizontal velocity in x at time n-1',''  &
!                        ,-1,0.,'east_e north_e zz',.true.)
!      vf_varid = def_var_pnetcdf(ncid,'vf',4,vdims,NF_FLOAT         &
!                        ,'horizontal velocity in x at time n-1',''  &
!                        ,-1,0.,'east_e north_e zz',.true.)
!      w_varid = def_var_pnetcdf(ncid,'w',4,vdims,NF_FLOAT         &
!                        ,'sigma coordinate vertical velocity',''  &
!                        ,-1,0.,'east_e north_e zz',.true.)
!      wr_varid = def_var_pnetcdf(ncid,'wr',4,vdims,NF_FLOAT          &
!                        ,'real (z coordinate) vertical velocity',''  &
!                        ,-1,0.,'east_e north_e zz',.true.)
!
!! end definitions
!      call check( nf90mpi_enddef(ncid)                     &
!                , 'nf_enddef: debug file @ '//out_file )
!
!
!! write data
!      start(1) = 1
!      edge(1) = 1
!      out1z = real(time,4)
!      call check( nfmpi_put_vara_real_all                  &
!                  ( ncid,time_varid,start,edge,out1z )     &
!                , 'nf_put_vara_real:time @ '//out_file )
!
!      start(1) = 1
!      edge(1) = kb
!      out1z = real(z,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,z_varid,start,edge,out1z )      &
!                , 'nf_put_var_real: z @ '//out_file )
!      out1z = real(zz,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,zz_varid,start,edge,out1z )      &
!                , 'nf_put_var_real: zz @ '//out_file )
!
!! East
!      if ( n_east == -1 ) then
!        start(1) = j_global(1)
!        edge(1) = jm
!        out1y = real(el_bry%est(1,:),4)
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,ele_varid,start,edge,out1y )  &
!                  , 'nf_put_var_real: ele @ '//out_file )
!      else
!        start(1) = 1
!        edge(1) = 0
!        out1y = 0.
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,ele_varid,start,edge,out1y )  &
!                  , 'nf_put_var_real: ele (dummy) @ '//out_file )
!      end if
!! West
!      if ( n_west == -1 ) then
!        start(1) = j_global(1)
!        edge(1) = jm
!        out1y = real(el_bry%wst(1,:),4)
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,elw_varid,start,edge,out1y )  &
!                  , 'nf_put_var_real: elw @ '//out_file )
!      else
!        start(1) = 1
!        edge(1) = 0
!        out1y = 0.
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,elw_varid,start,edge,out1y )  &
!                  , 'nf_put_var_real: elw (dummy) @ '//out_file )
!      end if
!! South
!      if ( n_south == -1 ) then
!        start(1) = i_global(1)
!        edge(1) = im
!        out1x = real(el_bry%sth(:,1),4)
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,els_varid,start,edge,out1x )  &
!                  , 'nf_put_var_real: els @ '//out_file )
!      else
!        start(1) = 1
!        edge(1) = 0
!        out1x = 0.
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,els_varid,start,edge,out1x )  &
!                  , 'nf_put_var_real: els (dummy) @ '//out_file )
!      end if
!! North
!      if ( n_north == -1 ) then
!        start(1) = i_global(1)
!        edge(1) = im
!        out1x = real(el_bry%nth(:,1),4)
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,eln_varid,start,edge,out1x )  &
!                  , 'nf_put_var_real: eln @ '//out_file )
!      else
!        start(1) = 1
!        edge(1) = 0
!        out1x = 0.
!        call check( nfmpi_put_vara_real_all              &
!                    ( ncid,eln_varid,start,edge,out1x )  &
!                  , 'nf_put_var_real: eln @ (dummy)'//out_file )
!      end if
!
!      call exchange2d_mpi(aam2d,im,jm)
!      call exchange2d_mpi(advua,im,jm)
!      call exchange2d_mpi(advva,im,jm)
!      call exchange2d_mpi(adx2d,im,jm)
!      call exchange2d_mpi(ady2d,im,jm)
!      call exchange2d_mpi(drx2d,im,jm)
!      call exchange2d_mpi(dry2d,im,jm)
!      call exchange2d_mpi(egb,im,jm)
!      call exchange2d_mpi(egf,im,jm)
!      call exchange2d_mpi(el,im,jm)
!      call exchange2d_mpi(elb,im,jm)
!      call exchange2d_mpi(elf,im,jm)
!      call exchange2d_mpi(et,im,jm)
!      call exchange2d_mpi(etb,im,jm)
!      call exchange2d_mpi(etb,im,jm)
!      call exchange2d_mpi(fluxua,im,jm)
!      call exchange2d_mpi(fluxva,im,jm)
!      call exchange2d_mpi(ua,im,jm)
!      call exchange2d_mpi(uab,im,jm)
!      call exchange2d_mpi(uaf,im,jm)
!      call exchange2d_mpi(va,im,jm)
!      call exchange2d_mpi(vab,im,jm)
!      call exchange2d_mpi(vaf,im,jm)
!
!      call exchange3d_mpi(aam,im,jm,kb)
!      call exchange3d_mpi(advx,im,jm,kb)
!      call exchange3d_mpi(advy,im,jm,kb)
!      call exchange3d_mpi(drhox,im,jm,kb)
!      call exchange3d_mpi(drhoy,im,jm,kb)
!      call exchange3d_mpi(kq,im,jm,kb)
!      call exchange3d_mpi(l,im,jm,kb)
!      call exchange3d_mpi(q2b,im,jm,kb)
!      call exchange3d_mpi(q2,im,jm,kb)
!      call exchange3d_mpi(q2l,im,jm,kb)
!      call exchange3d_mpi(q2lb,im,jm,kb)
!      call exchange3d_mpi(t,im,jm,kb)
!      call exchange3d_mpi(tb,im,jm,kb)
!      call exchange3d_mpi(s,im,jm,kb)
!      call exchange3d_mpi(sb,im,jm,kb)
!      call exchange3d_mpi(rho,im,jm,kb)
!      call exchange3d_mpi(kh,im,jm,kb)
!      call exchange3d_mpi(km,im,jm,kb)
!      call exchange3d_mpi(u,im,jm,kb)
!      call exchange3d_mpi(ub,im,jm,kb)
!      call exchange3d_mpi(uf,im,jm,kb)
!      call exchange3d_mpi(v,im,jm,kb)
!      call exchange3d_mpi(vb,im,jm,kb)
!      call exchange3d_mpi(vf,im,jm,kb)
!      call exchange3d_mpi(w,im,jm,kb)
!      call exchange3d_mpi(wr,im,jm,kb)
!
!      call exchange3d_mpi(el_amp,im,jm,ncons)
!
!      start(1) = i_global(1)
!      start(2) = j_global(1)
!      edge(1) = im
!      edge(2) = jm
!      out2 = real(aam2d,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,aam2d_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: aam2d @ '//out_file )
!      out2 = real(advua,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,advua_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: advua @ '//out_file )
!      out2 = real(advva,4)
!      call check( nfmpi_put_vara_real_all                    &
!                  ( ncid,advva_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: advva @ '//out_file )
!      out2 = real(adx2d,4)
!      call check( nfmpi_put_vara_real_all                    &
!                  ( ncid,adx2d_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: adx2d @ '//out_file )
!      out2 = real(ady2d,4)
!      call check( nfmpi_put_vara_real_all                    &
!                  ( ncid,ady2d_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: ady2d @ '//out_file )
!      out2 = real(drx2d,4)
!      call check( nfmpi_put_vara_real_all                    &
!                  ( ncid,drx2d_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: drx2d @ '//out_file )
!      out2 = real(dry2d,4)
!      call check( nfmpi_put_vara_real_all                     &
!                  ( ncid,dry2d_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: dry2d @ '//out_file )
!      out2 = real(egb,4)
!      call check( nfmpi_put_vara_real_all                     &
!                  ( ncid,egb_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: egb @ '//out_file )
!      out2 = real(egf,4)
!      call check( nfmpi_put_vara_real_all                     &
!                  ( ncid,egf_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: egf @ '//out_file )
!      out2 = real(el,4)
!      call check( nfmpi_put_vara_real_all                     &
!                  ( ncid,el_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: el @ '//out_file )
!      out2 = real(elb,4)
!      call check( nfmpi_put_vara_real_all                 &
!                  ( ncid,elb_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: elb @ '//out_file )
!      out2 = real(elf,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,elf_varid,start,edge,out2 )    &
!                , 'nf_put_var_real: elf @ '//out_file )
!      out2 = real(et,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,et_varid,start,edge,out2 )     &
!                , 'nf_put_var_real: et @ '//out_file )
!      out2 = real(etb,4)
!      call check( nfmpi_put_vara_real_all                 &
!                  ( ncid,etb_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: etb @ '//out_file )
!      out2 = real(etf,4)
!      call check( nfmpi_put_vara_real_all                 &
!                  ( ncid,etf_varid,start,edge,out2 )      &
!                , 'nf_put_var_real: etf @ '//out_file )
!      out2 = real(fluxua,4)
!      call check( nfmpi_put_vara_real_all                  &
!                  ( ncid,fluxua_varid,start,edge,out2 )    &
!                , 'nf_put_var_real: fluxua @ '//out_file )
!      out2 = real(fluxva,4)
!      call check( nfmpi_put_vara_real_all                  &
!                  ( ncid,fluxva_varid,start,edge,out2 )    &
!                , 'nf_put_var_real: fluxva @ '//out_file )
!      out2 = real(ua,4)
!      call check( nfmpi_put_vara_real_all              &
!                  ( ncid,ua_varid,start,edge,out2 )    &
!                , 'nf_put_var_real: ua @ '//out_file )
!      out2 = real(uab,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,uab_varid,start,edge,out2 )    &
!                , 'nf_put_var_real: uab @ '//out_file )
!      out2 = real(uaf,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,uaf_varid,start,edge,out2 )    &
!                , 'nf_put_var_real: uaf @ '//out_file )
!      out2 = real(va,4)
!      call check( nfmpi_put_vara_real_all              &
!                  ( ncid,va_varid,start,edge,out2 )    &
!                , 'nf_put_var_real: va @ '//out_file )
!      out2 = real(vab,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,vab_varid,start,edge,out2 )    &
!                , 'nf_put_var_real: vab @ '//out_file )
!      out2 = real(vaf,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,vaf_varid,start,edge,out2 )    &
!                , 'nf_put_var_real: vaf @ '//out_file )
!
!      start(1) = i_global(1)
!      start(2) = j_global(1)
!      start(3) = 1
!      start(4) = 1
!      edge(1) = im
!      edge(2) = jm
!      edge(3) = ncons
!      edge(4) = 1
!      out3 = real(el_amp,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,el_amp_varid,start,edge,out3 )    &
!                , 'nf_put_var_real: el_amp @ '//out_file )
!
!      start(1) = i_global(1)
!      start(2) = j_global(1)
!      start(3) = 1
!      start(4) = 1
!      edge(1) = im
!      edge(2) = jm
!      edge(3) = kb
!      edge(4) = 1
!      out3 = real(aam,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,aam_varid,start,edge,out3 )     &
!                , 'nf_put_vara_real: aam @ '//out_file )
!      out3 = real(advx,4)
!      call check( nfmpi_put_vara_real_all                 &
!                  ( ncid,advx_varid,start,edge,out3 )     &
!                , 'nf_put_vara_real: advx @ '//out_file )
!      out3 = real(advy,4)
!      call check( nfmpi_put_vara_real_all                 &
!                  ( ncid,advy_varid,start,edge,out3 )     &
!                , 'nf_put_vara_real: advy @ '//out_file )
!      out3 = real(drhox,4)
!      call check( nfmpi_put_vara_real_all                  &
!                  ( ncid,drhox_varid,start,edge,out3 )     &
!                , 'nf_put_vara_real: drhox @ '//out_file )
!      out3 = real(drhoy,4)
!      call check( nfmpi_put_vara_real_all                  &
!                  ( ncid,drhoy_varid,start,edge,out3 )     &
!                , 'nf_put_vara_real: drhoy @ '//out_file )
!      out3 = real(kq,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,kq_varid,start,edge,out3 )     &
!                , 'nf_put_vara_real: kq @ '//out_file )
!      out3 = real(l,4)
!      call check( nfmpi_put_vara_real_all              &
!                  ( ncid,l_varid,start,edge,out3 )     &
!                , 'nf_put_vara_real: l @ '//out_file )
!      out3 = real(q2b,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,q2b_varid,start,edge,out3 )     &
!                , 'nf_put_vara_real: q2b @ '//out_file )
!      out3 = real(q2,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,q2_varid,start,edge,out3 )     &
!                , 'nf_put_vara_real: q2 @ '//out_file )
!      out3 = real(q2lb,4)
!      call check( nfmpi_put_vara_real_all                 &
!                  ( ncid,q2lb_varid,start,edge,out3 )     &
!                , 'nf_put_vara_real: q2lb @ '//out_file )
!      out3 = real(q2l,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,q2l_varid,start,edge,out3 )     &
!                , 'nf_put_vara_real: q2l @ '//out_file )
!      out3 = real(t,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,t_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: t @ '//out_file )
!      out3 = real(tb,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,tb_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: tb @ '//out_file )
!      out3 = real(s,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,s_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: s @ '//out_file )
!      out3 = real(sb,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,sb_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: sb @ '//out_file )
!      out3 = real(rho,4)
!      call check( nfmpi_put_vara_real_all                 &
!                  ( ncid,rho_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: rho @ '//out_file )
!      out3 = real(kh,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,kh_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: kh @ '//out_file )
!      out3 = real(km,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,km_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: km @ '//out_file )
!      out3 = real(ub,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,ub_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: ub @ '//out_file )
!      out3 = real(u,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,u_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: u @ '//out_file )
!      out3 = real(uf,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,uf_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: uf @ '//out_file )
!      out3 = real(vb,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,vb_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: vb @ '//out_file )
!      out3 = real(v,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,v_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: v @ '//out_file )
!      out3 = real(vf,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,vf_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: vf @ '//out_file )
!      out3 = real(w,4)
!      call check( nfmpi_put_vara_real_all               &
!                  ( ncid,w_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: w @ '//out_file )
!      out3 = real(wr,4)
!      call check( nfmpi_put_vara_real_all                &
!                  ( ncid,wr_varid,start,edge,out3 )      &
!                , 'nf_put_vara_real: wr @ '//out_file )
!
!! close file
!      call check( nf90mpi_close(ncid)                   &
!                , 'nf_close: output:  @ '//out_file )
!
!
!    end ! subroutine write_debug_pnetcdf
!!
!!______________________________________________________________________
!!
!    subroutine read_initial_ts_pnetcdf( temp, salt, ssh, n )
!!----------------------------------------------------------------------
!!  Reads initial temperature, salinity and elevation.
!! TODO: Add u and v. Add check for single record to read if not climatology.
!!______________________________________________________________________
!!
!      use glob_domain
!!      use grid       , only: fsm
!      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
!      use pnetcdf
!
!      implicit none
!
!      integer, external :: get_var_real_2d &
!                         , get_var_real_3d
!
!      integer , intent(in)  :: n
!      real(rk), intent(out) :: temp(im,jm,kb),salt(im,jm,kb),ssh(im,jm)
!
!      integer                  ncid, status
!      integer                  el_varid, sb_varid, tb_varid
!      integer(MPI_OFFSET_KIND) start(4),edge(4)
!
!
!      start = 1
!      edge  = 1
!
!!  Open netcdf file.
!      if ( is_master )                                              &
!           print '(''reading file `'',a,''`''/)', trim(ic_path)
!      call check( nf90mpi_open(                        &
!                    POM_COMM, ic_path, NF_NOWRITE  &
!                  , MPI_INFO_NULL, ncid )              &
!                , "nfmpi_open: "//ic_path )
!
!!  Get variables. [ TODO: parameterize varnames ]
!      call check( nf90mpi_inq_varid( ncid, 'tclim', tb_varid )  &
!                , "nfmpi_inq_varid: t_init @ "//ic_path )
!      call check( nf90mpi_inq_varid( ncid, 'sclim', sb_varid )  &
!                , "nfmpi_inq_varid: s_init @ "//ic_path )
!
!      status = nf90mpi_inq_varid( ncid, 'eclim', el_varid )
!      if ( status == NF_NOERR ) then
!! Get elevation.
!        start(1) = i_global(1)
!        start(2) = j_global(1)
!        start(3) = n
!        edge(1) = im
!        edge(2) = jm
!        edge(3) = 1
!        call check( get_var_real_2d( ncid, el_varid, start,edge, ssh ) &
!                  , "nfmpi_get_vara_real_all"//ic_path )
!      else
!        if ( status == -49 ) then
!          if ( is_master ) then
!            call msg_print("", 2                                       &
!                        , "Missing elevation data - falling back to 0.")
!          end if
!          ssh = 0.
!        else
!          call check(status, 'nfmpi_inq_varid: el_init @ '//ic_path)
!        end if
!      end if
!
!      start(1) = i_global(1)
!      start(2) = j_global(1)
!      start(3) = 1
!      start(4) = n
!      edge(1) = im
!      edge(2) = jm
!      edge(3) = kb
!      edge(4) = 1
!      call check( get_var_real_3d(ncid,tb_varid,start,edge,temp)  &
!                , "nf_get_var : temp @ "//ic_path )
!      call check( get_var_real_3d(ncid,sb_varid,start,edge,salt)  &
!                , "nf_get_var : salt @ "//ic_path )
!
!!! TODO: move out
!!      do k=1,kb-1
!!        where (fsm==0.)
!!          temp(:,:,k) = 0.
!!          salt(:,:,k) = 0.
!!        end where
!!      end do
!!      temp(:,:,kb) = temp(:,:,kb-1)
!!      salt(:,:,kb) = salt(:,:,kb-1)
!!      elb = elb*fsm
!
!!  Close file.
!      call check( nf90mpi_close(ncid), "nf_close @ "//ic_path)
!
!
!    end ! subroutine read_initial_ts_pnetcdf
!!
!!______________________________________________________________________
!!
!    subroutine read_ts_z_pnetcdf( temp, salt, ks, n, path )
!!----------------------------------------------------------------------
!!  Reads initial temperature, salinity and elevation.
!! TODO: Add u and v. Add check for single record to read if not climatology.
!!______________________________________________________________________
!!
!      use glob_domain
!      use grid       , only: zz, h
!      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
!      use pnetcdf
!
!      implicit none
!
!      integer, external :: get_var_real_2d &
!                         , get_var_real_3d
!
!      integer                      , intent(in   ) :: ks, n
!      character(len=*)             , intent(in   ) :: path
!      real(rk), dimension(im,jm,kb), intent(  out) :: temp,salt
!
!      integer                  ncid
!      integer                  z_varid, sb_varid, tb_varid
!      integer(MPI_OFFSET_KIND) start(4),edge(4)
!      real(rk)                 tz(im,jm,ks),sz(im,jm,ks),z_z(ks)
!
!
!      start = 1
!      edge  = 1
!
!!  Open netcdf file.
!      if ( is_master )                                              &
!           print '(''reading file `'',a,''`''/)', trim(ic_path)
!      call check( nf90mpi_open(                        &
!                    POM_COMM, path, NF_NOWRITE  &
!                  , MPI_INFO_NULL, ncid )              &
!                , "nfmpi_open: "//path )
!
!!  Get variables. [ TODO: parameterize varnames ]
!      call check( nf90mpi_inq_varid( ncid, 't', tb_varid )  &
!                , "nfmpi_inq_varid: t_init @ "//ic_path )
!      call check( nf90mpi_inq_varid( ncid, 's', sb_varid )  &
!                , "nfmpi_inq_varid: s_init @ "//ic_path )
!      call check( nf90mpi_inq_varid( ncid, 'z', z_varid )  &
!                , "nfmpi_inq_varid: z_init @ "//ic_path )
!
!      start(1) = 1
!      edge(1) = ks
!      call check( nf90mpi_get_var_all(ncid,z_varid,z_z,start,edge)  &
!                , "nf_get_var : z @ "//ic_path )
!
!      start(1) = i_global(1)
!      start(2) = j_global(1)
!      start(3) = 1
!      start(4) = n
!      edge(1) = im
!      edge(2) = jm
!      edge(3) = ks
!      edge(4) = 1
!      call check( get_var_real_3d(ncid,tb_varid,start,edge,tz)  &
!                , "nf_get_var : temp @ "//ic_path )
!      call check( get_var_real_3d(ncid,sb_varid,start,edge,sz)  &
!                , "nf_get_var : salt @ "//ic_path )
!
!!  Close file.
!      call check( nf90mpi_close(ncid), "nf_close @ "//ic_path)
!
!      call ztosig(z_z,tz,zz,h,temp,im,jm,ks,kb,n_west,n_east,n_south,n_north)
!      call ztosig(z_z,sz,zz,h,salt,im,jm,ks,kb,n_west,n_east,n_south,n_north)
!
!
!    end ! subroutine read_ts_z_pnetcdf
!!
!!______________________________________________________________________
!!
!    subroutine out_debug( out_file )
!!----------------------------------------------------------------------
!!  Wrapper for debug output procedure.
!!______________________________________________________________________
!!
!      implicit none
!
!      character(len=*), intent(in) :: out_file
!
!      character(len=4), parameter :: pfx = "dbg."
!
!
!      call write_debug_pnetcdf(      trim(out_path)    &
!                              //          pfx          &
!                              //     trim(out_file)    &
!                              //"."//trim(FORMAT_EXT) )
!
!
!    end ! subroutine out_debug
!______________________________________________________________________
!
      subroutine read_initial( init_file, rec, elev, salt, temp, u, v )
!----------------------------------------------------------------------
!  Reads initial fields.
!______________________________________________________________________
!
      use config     , only: el_name, s_name, t_name, u_name, v_name
      use glob_const , only: PATH_LEN, rk
      use glob_domain, only: i_global, im, is_master, j_global, jm, kb
      use io         , only: check, is_error, file_close, file_open
     &                     , var_read
      use mpi        , only: MPI_OFFSET_KIND

      implicit none

      character(*)                     , intent(in   ) :: init_file
      integer                          , intent(in   ) :: rec
      real(rk)    , dimension(im,jm   ), intent(inout) :: elev
      real(rk)    , dimension(im,jm,kb), intent(inout) :: salt, temp
     &                                                  , u   , v

      integer(MPI_OFFSET_KIND), dimension(4) :: start, stride
      integer                                   ncid, status


      if ( is_master ) then
        print '(/"Reading initial file ",a)', init_file
      end if

      ncid = file_open( init_file )

      start  = [ i_global(1), j_global(1),  1, rec ]
      stride = [ im         , jm         ,  1,   1 ]

      status = var_read( ncid, el_name, elev, start, stride )
      if ( is_error(status) ) then
        elev = 0.
        if ( is_master ) then
          print *, "[!] Elevation read fail... Fallback. "
!          print *, " @ initialize.f:read_initial" ! For debug purposes
        end if
      end if

      stride = [ im         , jm         , kb,   1 ]

      if ( s_name == '' ) then
        salt = 33.
      else
       call check( var_read( ncid, s_name, salt, start, stride )
     &           , " @ initialize.f:read_initial:var_read:"//s_name )
      end if

      status = var_read( ncid, t_name, temp, start, stride )
      if ( is_error(status) ) then
        temp = 10.
        if ( is_master ) then
          print *, "[!] Temperature read fail... Fallback. "
!          print *, " @ initialize.f:read_initial" ! For debug purposes
        end if
      end if

      status = var_read( ncid, u_name, u, start, stride )
      if ( is_error(status) ) then
        u = 0.
        if ( is_master ) then
          print *, "[!] Meridional velocity read fail... Fallback. "
!          print *, " @ initialize.f:read_initial" ! For debug purposes
        end if
      end if
      u = u/100.

      status = var_read( ncid, v_name, v, start, stride )
      if ( is_error(status) ) then
        v = 0.
        if ( is_master ) then
          print *, "[!] Zonal velocity read fail... Fallback. "
!          print *, " @ initialize.f:read_initial" ! For debug purposes
        end if
      end if
      v = v/100.

! close file
      ncid = file_close( ncid )


      end ! subroutine read_initial
!
!______________________________________________________________________
!
      subroutine out_init( out_file )
!----------------------------------------------------------------------
!  Wrapper for output procedure.
!______________________________________________________________________
!
      use air        , only: air_out_init_define => out_init_define
      use bry        , only: bry_out_init_define => out_init_define
      use config     , only: title
      use glob_const , only: PATH_LEN, rk
      use glob_domain, only: i_global, im, im_global, is_master
     &                     , j_global, jm, jm_global, kb
      use grid       , only: grid_out_init_define => out_init_define
     &                     , grid_out_write       => out_write
      use io
      use model_run  , only: time, time_start
      use mpi        , only: MPI_OFFSET_KIND
      use glob_ocean , only: kh, km, rho, sb, tb
      use pnetcdf    , only: NF90_FLOAT,NF_GLOBAL,NF90_DOUBLE,NF_DOUBLE
     &                     ,nf90mpi_inq_varid,nf90mpi_put_var_all
     &                     ,nfmpi_put_vara_double_all

      implicit none

      character(len=*), intent(in) :: out_file

      integer    ncid, time_dimid, varid
     &      , x_dimid,    y_dimid, z_dimid
      integer vdims(4)
      integer(MPI_OFFSET_KIND), dimension(4) :: start, stride

      real(rk), dimension(im,jm)    :: out1


      out1 = 0.

      if ( is_master ) then
        write(*,'(/''writing initial file '',a)') out_file
        print *, "NF_DOUBLE  : ", NF_DOUBLE
        print *, "NF90_DOUBLE: ", NF90_DOUBLE
      end if

      ncid = file_create( out_file )
! define global attributes
      call att_write( ncid, -1, 'title', trim(title) )
      call att_write( ncid, -1, 'description', "Init file" )

! define dimensions
      time_dimid = dim_define( ncid, 'time',         -1_8  )
      z_dimid    = dim_define( ncid, 'z', int(kb       ,8) )
      y_dimid    = dim_define( ncid, 'y', int(jm_global,8) )
      x_dimid    = dim_define( ncid, 'x', int(im_global,8) )

! define variables and their attributes
      vdims(1) = time_dimid
      varid = var_define( ncid, 'time', 1, vdims, NF90_FLOAT
     &                  , 'time', "days since "//time_start
     &                  , -1, 0., ' ', .false. )

! call modules
      call grid_out_init_define( ncid, x_dimid, y_dimid, z_dimid )
      call  bry_out_init_define( ncid, x_dimid, y_dimid )
      call  air_out_init_define( ncid, x_dimid, y_dimid, time_dimid )

      vdims = [ x_dimid, y_dimid, z_dimid, time_dimid ]
      varid = var_define( ncid, 't', 3, vdims, NF90_FLOAT
     &                  , "potential temperature", "K"
     &                  , 0, 100., "east_e north_e zz", .true. )
      varid = var_define( ncid, 's', 3, vdims, NF90_FLOAT
     &                  , "salinity x rho / rhoref", "pss"
     &                  , 1, 0., "east_e north_e zz", .true. )
      varid = var_define( ncid, 'rho', 3, vdims, NF90_FLOAT
     &                  , "(density-1000)/rhoref", "dimensionless"
     &                  , 0, 0., "east_e north_e zz", .true. )
      varid = var_define( ncid, 'kh', 3, vdims, NF90_FLOAT
     &                  , "vertical diffusivity", "m^2 s^-1"
     &                  , -1, 0., "east_e north_e zz", .true. )
      varid = var_define( ncid, 'km', 3, vdims, NF90_FLOAT
     &                  , "vertical viscosity", "m^2 s^-1"
     &                  , -1, 0., "east_e north_e zz", .true. )

! end definitions
      call file_end_definition( ncid )


! write data
      start(1)  = 1
      stride(1) = 1
      out1(1,1) = time
      call var_write( ncid, 'time', out1(1:1,1), start, stride )

!      call bry_out_write( ncid )
      call grid_out_write( ncid )

      start  = [ i_global(1), j_global(1),  1, 1 ]
      stride = [ im         , jm         , kb, 1 ]

      call var_write( ncid, 's'  ,  sb, start, stride )
      call var_write( ncid, 't'  ,  tb, start, stride )
      call var_write( ncid, 'rho', rho, start, stride )
      call var_write( ncid, 'kh' ,  kh, start, stride )
      call var_write( ncid, 'km' ,  km, start, stride )

! close file
      ncid = file_close( ncid )


      end ! subroutine out_init
      
!      elemental real(rk) function convert_units( var, from, to )
!
!      use config, only: rk
!
!      implicit none
!
!      real(rk)    , intent(in) :: var
!      character(*), intent(in) :: from, to
!
!      integer                    , parameter :: N = 7
!      character(10), dimension(N), parameter :: unit = [
!     &             "mm ", "cm ", "dm ", "m  ", "km ", "Pa ", "hPa" ]
!      real(rk), dimension(N,N), parameter :: conv =
!     &          [ [ 1.   , 1.e-1, 1.e-2, 1.e-3, 1.e-6, 0., 0. ]
!     &          , [ 1.e+1, 1.   , 1.e-1, 1.e-2, 1.e-5, 0., 0. ]
!     &          , [ 1.e+2, 1.e+1, 1.   , 1.e-1, 1.e-4, 0., 0. ]
!     &          , [ 1.e+3, 1.e+2, 1.e-1, 1.   , 1.e-3, 0., 0. ]
!     &          , [ 1.e+6, 1.e+5, 1.e+4, 1.e+3, 1.   , 0., 0. ]
!     &          , [ 0.   , 0.   , 0.   , 0.   , 0.   , 1.   , 1.e-2 ]
!     &          , [ 0.   , 0.   , 0.   , 0.   , 0.   , 1.e+2, 1.    ]
!     &          ]
!      integer f_ind, t_ind
!
!      f_ind = 0
!      t_ind = 0
!      f_ind = findloc( unit, from, 1 )
!      t_ind = findloc( unit, to  , 1 )
!      
!      if ( f_ind /= 0 .and. t_ind /= 0 ) then
!        convert_units == var*conv(f_ind,t_ind)
!      end if
!          
!      end ! function convert_units