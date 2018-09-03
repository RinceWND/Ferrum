! initialize.f

! initialize POM: define constant, read initial values, grid and initial
! conditions

!______________________________________________________________________
!
      subroutine initialize
!----------------------------------------------------------------------
!  Initializes POM.
!______________________________________________________________________

      use bry        , only: aamfrz, frz
      use glob_config
      use glob_const , only: SEC2DAY
      use glob_domain
      use glob_grid  , only: dum, dvm, east_e, north_e, fsm, h
      use glob_misc  , only: aamfac
      use glob_ocean , only: d, dt, el, et
      use glob_time  , only: dti, time, time0

      implicit none

      integer i,j
      integer, parameter :: TEST_CASE = 1
      real(kind=rk) rdisp !lyo:20110224:alu:stcc:dispensible real variable!lyo:pac10:
!     integer io_global,jo_global,io,jo !fhx:incmix!lyo:pac10:delete
!     real(kind=rk) xs,ys,fak!,lono,lato!lyo:pac10:moved to pom.h: !fhx:incmix


! Read domain distribution (TODO: get rid of hardcoded parameter)
      call read_domain_dist( 'domain.nml' )

! initialize the MPI execution environment and create communicator for
!internal POM communications
      call initialize_mpi

! distribute the model domain across processors
      call distribute_mpi

! distribute the coarse wind/assim data domain across processors
      call distribute_mpi_coarse   !fhx:20101206

! fhx:interp_flag
! fhx:fgrid
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

! allocate and initialize arrays (TODO: allocate them before or after initializing MPI?)
      call initialize_arrays

! read input values and define constants
      call read_input

! read in grid data
      if ( TEST_CASE == 0 ) then
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
      do j=1,jm
        do i=2,im
          dum(i,j) = fsm(i,j)*fsm(i-1,j)
        end do
      end do
      dum(1,:) = dum(2,:)
      do j=2,jm
        do i=1,im
          dvm(i,j) = fsm(i,j)*fsm(i,j-1)
        end do
      end do
      dvm(:,1) = dvm(:,2)

      call exchange2d_mpi(  h,im,jm)
      call exchange2d_mpi(fsm,im,jm)
      call exchange2d_mpi(dum,im,jm)
      call exchange2d_mpi(dvm,im,jm)

! read initial and lateral boundary conditions
      call initial_conditions

! read laterel boundary conditions
      call lateral_boundary_conditions !ayumi:20100407

!lyo:pac10:beg:
! set lateral boundary FRZ & SPONGE layer parameters !lyo:20110224:alu:stcc:
!      call bfrz(   nfw,   nfe,    nfs,    nfn,
!     $          n_west,n_east,n_south,n_north,
!     $              im,    jm,      6,    frz)       !FRZ:
!     For FRZ in bcond: T[n+1]={ (1-frz)*T[n] + frz*Tobs },
!     the inverse-relax-time = frz/dti,
!     so that since "bfrz" gives frz=1 @boundaries,
!     frz=1 (=dti/dti) -> inverse-relax-time = 1/dti @boundaries, and
!     frz=dti/86400.   -> inverse-relax-time = 1/86400 @boundaries etc.
      rdisp=dti*SEC2DAY; frz(:,:)=frz(:,:)*rdisp !lyo:stcc:mar_:dec_:
!     rdisp=1.;         frz(:,:)=frz(:,:)*rdisp !lyo:stcc:
!     For aam is increased by (1.+aamfrz(i,j)) in subr lateral_viscosity
      call bfrz(     0,     0,      0,      0,
     $          n_west,n_east,n_south,n_north,
     $              im,    jm,      6, aamfrz)       !SPONGE:
      rdisp=0.; aamfrz(:,:)=aamfrz(:,:)*rdisp
!lyo:pac10:end:

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
      if ( nread_rst /= 0 ) then
        call read_restart_pnetcdf
! update elevation
        d  = h+el
        dt = h+et
! update time
        time=time0
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

      call msg_print("MODEL CORE INITIALIZED", 1, "")


      return

      end

!______________________________________________________________________
!
      subroutine read_input
!----------------------------------------------------------------------
!  Reads input values and defines constants.
!______________________________________________________________________

      use glob_config
      use glob_const
      use glob_domain, only: is_master
      use glob_misc  , only: iprints
      use glob_time
      use module_time

      implicit none

      include 'io.h'
      include 'bulk.h'

      namelist/pom_nml/ title,netcdf_file,mode,nadv,nitera,sw,npg,dte, !fhx:Toni:npg
     $                  isplit,time_start,nread_rst,read_rst_file
     $                 ,write_rst,write_rst_file,days,prtd1,prtd2
     $                 ,iperx,ipery,n1d !lyo:scs1d:add iper* & n1d in pom.nml
                                        !n1d .ne. 0 for 1d simulation
     $                 ,windf           !lyo: windf = ccmp, ecmw, or gfsw etc
     $                 ,nbct,nbcs,tprni,umol,z0b,ntp

      namelist/switch_nml/
     $     calc_bulk, calc_wind, calc_tsforce, calc_river, calc_assim,
     $     calc_assimdrf,         !eda
     $     calc_tsurf_mc, calc_tide,!fhx:mcsst; fhx:tide
     $     calc_uvforce,         !eda:uvforce
     $     calc_ice,
     $     output_flag, SURF_flag, monthly_flag !fhx:20110131:
      type(date) ::  dtime, dtime0

      namelist/sensitivity_nml/ sf_bf, sf_hf, sf_wi
      namelist/misc_nml/ spinup, t_lo, t_hi
      namelist/input_files_nml/
     &  grid_path, bc_path, clim_path, mean_path
     &, bulk_path, flux_path
      namelist/input_vars_nml/
     &  dlwr_name, lht_name, lwr_name, prat_name, pres_name
     &, rhum_name, sht_name, sst_name, swr_name, tair_name
     &, tcld_name, umf_name, uwnd_name, vmf_name, vwnd_name
     &, wgst_name
      namelist/bulk_nml/ calc_swr, use_coare, lwrad_formula

! Input of filenames and constants

! Logical for inertial ramp (.true. if inertial ramp to be applied
! to wind stress and baroclinic forcing, otherwise .false.)
!      lramp=.false.

! Reference density (recommended values: 1025 for seawater,
! 1000 for freswater; S.I. units):
      rhoref=1025.

! Temperature bias (deg. C)
      tbias=0.

! Salinity bias
      sbias=0.

! Bottom roughness (metres)
      z0b=.01

! Minimum bottom friction coeff.
      cbcmin=.0025

! Maximum bottom friction coeff.
      cbcmax=5. !1.e0   !lyo:20110315:botwavedrag:

! Smagorinsky diffusivity coeff.
      horcon=0.2  !lyo:pac10:exp004: !=0.1e0

! Inverse horizontal turbulent Prandtl number (ah/am; dimensionless):
! NOTE that tprni=0.e0 yields zero horizontal diffusivity!
      tprni=.2

! Background viscosity used in subroutines profq, proft, profu and
! profv (S.I. units):
      umol=2.e-5

! Maximum magnitude of vaf (used in check that essentially tests
! for CFL violation):
      vmaxl=5. !lyo:debug:100.e0

! Maximum allowable value of:
!   <difference of depths>/<sum of depths>
! for two adjacent cells (dimensionless). This is used in subroutine
! slpmax. If >= 1, then slpmax is not applied:
      slmax=2.

! Water type, used in subroutine proft.
!    ntp    Jerlov water type
!     1            i
!     2            ia
!     3            ib
!     4            ii
!     5            iii
      ntp=2

! Surface temperature boundary condition, used in subroutine proft:
!    nbct   prescribed    prescribed   short wave
!           temperature      flux      penetration
!     1        no           yes           no
!     2        no           yes           yes
!     3        yes          no            no
!     4        yes          no            yes
      nbct=1

! Surface salinity boundary condition, used in subroutine proft:
!    nbcs   prescribed    prescribed
!            salinity      flux
!     1        no           yes
!     3        yes          no
! NOTE that only 1 and 3 are allowed for salinity.
      nbcs=1

! Step interval during which external (2-D) mode advective terms are
! not updated (dimensionless):
      ispadv=5

! Constant in temporal filter used to prevent solution splitting
! (dimensionless):
      smoth=0.10

! Weight used for surface slope term in external (2-D) dynamic
! equation (a value of alpha = 0.e0 is perfectly acceptable, but the
! value, alpha=.225e0 permits a longer time step):
      alpha=0.225

! Initial value of aam:
      aam_init=500.

! Bulk default parameters
!
! Use air-sea bulk calculations
      calc_bulk = .true.

! Calculate shortwave radiation according to Reed
      calc_swr  = .true.

! Do not use COARE algorithm for heat fluxes estimation by default
      use_coare = .false.


! Other parameters and flags
      sf_bf = 1.
      sf_hf = 1.
      sf_wi = 1.
      t_lo = -999.
      t_hi =  999.
      iperx = 0
      ipery = 0

! read input namelist
      open(73,file='pom.nml',status='old')
      read(73,nml=pom_nml)
      close(73)

! read main switches
      open(73,file='switch.nml',status='old')
      read(73,nml=switch_nml)
      read(73,nml=misc_nml)
      read(73,nml=sensitivity_nml)
      close(73)

! read main switches
      open(73,file='io.nml',status='old')
      read(73,nml=input_files_nml)
      read(73,nml=input_vars_nml)
      close(73)

! read bulk configuration
      open(73,file='bulk.nml',status='old')
      read(73,nml=bulk_nml)
      close(73)

! End of input of constants

! Initialize various constants
      call initialize_constants

      dti=dte*float(isplit)
      dte2=dte*2
      dti2=dti*2

      iend=max0(nint(days*24.*3600./dti),2)
      iprint=max(nint(prtd1*24.*3600./dti),1)
      irestart=nint(write_rst*24.*3600./dti)
      iprints=nint(prtd2*24.*3600./dti) !fhx:20110131:add 3hrly output

      ispi=1./float(isplit)
      isp2i=1./(2.*float(isplit))

! Initialise time
! Do not offset time if not restarting
      if ( nread_rst /= 0 ) then

        dtime0 = str2date( time_start(1:19) )
        dtime = str2date( read_rst_file(1:13)//":"//
     &  read_rst_file(15:16)//":"//read_rst_file(18:19) )

        time0 = real( (dtime-dtime0)/86400, rk )

      else

        time0 = 0.

      end if

!      if (is_master) then
!         print*, 'dtime',dtime
!         print*, 'dtime0',dtime0
!         print*, time0
!      endif

!      time0=0.e0
      time=time0 !lyo:20110224:alu:stcc:!lyo:pac10:

! print initial summary
      if ( is_master ) then
        print '(/'' title      = '',a40)', title
        print '(/'' mode       = '',i10)', mode
        print '('' nadv       = '',i10)', nadv
        print '('' nitera     = '',i10)', nitera
        print '('' sw         = '',f10.4)', sw
        print '('' npg        = '',i10)', npg    !fhx:Toni:npg
        print '('' nread_rst  = '',i10)', nread_rst
        print '('' write_rst  = '',f10.4)', write_rst
        print '('' irestart   = '',i10)', irestart
        print '('' dte        = '',f10.2)', dte
        print '('' dti        = '',f10.1)', dti
        print '('' isplit     = '',i10)', isplit
        print '('' time_start = '',a26)', time_start
        print '('' days       = '',f10.4)', days
        print '('' iend       = '',i10)', iend
        print '('' prtd1      = '',f10.4)', prtd1
        print '('' iprint     = '',i10)', iprint
        print '('' prtd2      = '',f10.4)', prtd2
        print '('' rhoref     = '',f10.3)', rhoref
        print '('' tbias      = '',f10.3)', tbias
        print '('' sbias      = '',f10.3)', sbias
        print '('' grav       = '',f10.4)', grav
        print '('' Kappa      = '',f10.4)', Kappa
        print '('' z0b        = '',f10.6)', z0b
        print '('' cbcmin     = '',f10.6)', cbcmin
        print '('' cbcmax     = '',f10.6)', cbcmax
        print '('' horcon     = '',f10.3)', horcon
        print '('' tprni      = '',f10.4)', tprni
        print '('' umol       = '',f10.4)', umol
        print '('' vmaxl      = '',f10.4)', vmaxl
        print '('' slmax      = '',f10.4)', slmax
        print '('' lono,lato  = '',2f10.4)', lono,lato
        print '('' xs,ys,fak  = '',3f10.4)', xs,ys,fak
        print '('' iperx,ipery= '',2i10)', iperx,ipery
        print '('' ntp        = '',i10)', ntp
        print '('' nbct       = '',i10)', nbct
        print '('' nbcs       = '',i10)', nbcs
        print '('' ispadv     = '',i10)', ispadv
        print '('' smoth      = '',f10.4)', smoth
        print '('' alpha      = '',f10.4)', alpha
!        print '('' lramp      = '',l10)', lramp
        print '('' calc_wind      = '',l2)', calc_wind
        print '('' calc_tsforce   = '',l2)', calc_tsforce
        print '('' calc_river     = '',l2)', calc_river
        print '('' calc_assim     = '',l2)', calc_assim
        print '('' calc_assimdrf  = '',l2)', calc_assimdrf !eda
        print '('' calc_uvforce   = '',l2)', calc_uvforce !eda:uvforce
        print '('' calc_tsurf_mc  = '',l2)', calc_tsurf_mc !fhx:mcsst
        print '('' calc_tide      = '',l2)', calc_tide     !fhx:tide
        print '('' calc_interp    = '',l2)', calc_interp   !fhx:interp_flag
        print '('' output_flag    = '',i2)', output_flag  !fhx:20110131:
        print '('' SURF_flag      = '',i2)', SURF_flag    !fhx:20110131:
        print '(/'' Sensitivity:'')'
        print '(''   heat flux     = '',f10.4)', sf_hf
        print '(''   wind speed    = '',f10.4)', sf_wi
        print '(''   lat.velocities= '',f10.4)', sf_bf
        if ( calc_bulk ) then
          print '(/'' Air-Sea Bulk: [ ENABLED ]'')'
          print '(''   calc swrad    = '',l2)', calc_swr
          print '(''   use COARE     = '',l2)', use_coare
          print '(''   longwave form.= '',i2)', lwrad_formula
        else
          print '(/'' Air-Sea Bulk: [ DISABLED ]'')'
        end if
      end if

      return

      end

!______________________________________________________________________
!
      subroutine initialize_arrays
!----------------------------------------------------------------------
!  Allocate arrays and initialize them for safety.
!______________________________________________________________________

!     ayumi 2010/5/14
!      use river, only : totq
      use bry        , only: initialize_boundary
      use glob_atmos
      use glob_bry
      use glob_config, only: ntide, rk
      use glob_domain, only: im, im_coarse, jm, jm_coarse, kb, kbm1
      use glob_grid
      use glob_misc
      use glob_ocean
      use glob_out

      implicit none

      integer i,j,k

!     ayumi 2010/5/14
!      totq = 0.e0

!     ayumi 2010/6/13
      iout = 0

      iouts = 0 !fhx:20110131:


      call initialize_boundary

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
     $  e_atmos(im,jm) ,
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
     $  north_c(im,jm) ,
     $  north_e(im,jm) ,
     $  north_u(im,jm) ,
     $  north_v(im,jm) ,
     $  psi(im,jm)     ,
     $  rot(im,jm)     ,
     $  ssurf(im,jm)   ,
     $  swrad(im,jm)   ,
     $  vfluxb(im,jm)  ,
     $  tauiwu(im,jm)  ,    !:rwnd
     $  tauiwv(im,jm)  ,    !:rwnd
     $  tps(im,jm)     ,
     $  tsurf(im,jm)   ,
     $  ua(im,jm)      ,
     $  vfluxf(im,jm)  ,
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
     $  wssurf(im,jm)  ,
     $  wtsurf(im,jm)  ,
     $  wubot(im,jm)   ,
     $  wusurf(im,jm)  ,
     $  wvbot(im,jm)   ,
     $  wvsurf(im,jm)  ,
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
     $  rmean(im,jm,kb),
     $  sb(im,jm,kb)   ,
     $  sclim(im,jm,kb),
     $  s(im,jm,kb)    ,
     $  tb(im,jm,kb)   ,
     $  tclim(im,jm,kb),
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

      allocate(
     $  cibe(jm)       ,
     $  cibn(im)       ,
     $  cibs(im)       ,
     $  cibw(jm)       ,
     $  ampe(jm,ntide) ,
     $  phae(jm,ntide) ,
     $  amue(jm,ntide) ,
     $  phue(jm,ntide)
     & )


      allocate(
     $  usrf_mean(im,jm)    ,
     $  vsrf_mean(im,jm)    ,
     $  elsrf_mean(im,jm)   ,
     $  uwsrf_mean(im,jm)   ,
     $  vwsrf_mean(im,jm)   ,
     $  uwsrf(im,jm)        ,   ! wind data for SURF mean output
     $  vwsrf(im,jm)
     & )

! 2-D and 3-D arrays
      do j=1,jm
        do i=1,im
          cor(i,j)    =0.
          e_atmos(i,j)=0.
          vfluxb(i,j) =0.
          vfluxf(i,j) =0.
          wusurf(i,j) =0.
          wvsurf(i,j) =0.
          tauiwu(i,j) =0.
          tauiwv(i,j) =0.
          wtsurf(i,j) =0.
          wssurf(i,j) =0.
          swrad(i,j)  =0.
          drx2d(i,j)  =0.
          dry2d(i,j)  =0.
          uwsrf(i,j)  =0.     !fhx:20110131:wind u for SURF output
          vwsrf(i,j)  =0.     !fhx:20110131:wind v for SURF output
          wubot(i,j)  =0.     !lyo:20110315:botwavedrag:
          wvbot(i,j)  =0.     !lyo:20110315:botwavedrag:
          aamfac(i,j) =1.     !fhx:incmix
          hi(i,j)=.35
        end do
      end do

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            ub(i,j,k)=0.
            vb(i,j,k)=0.
          end do
        end do
      end do



      return
      end

!______________________________________________________________________
!
      subroutine read_grid
!----------------------------------------------------------------------
!  Set up vertical and horizontal grid, topography, areas and masks.
!______________________________________________________________________

      use glob_config, only: n1d, rk
      use glob_const , only: DEG2RAD, Ohm, SMALL
      use glob_domain, only: im, is_master, jm, kb, n_south, n_west
      use glob_grid
      use glob_ocean , only: d, dt, el, et

      implicit none

      integer i,j,k

! read grid
      call read_grid_pnetcdf
!______________________________________________________________________
!      The followings are read in read_grid_pnetcdf:
!     z,zz,dx,dy
!     east_u,east_v,east_e,east_c
!     north_u,north_v,north_e,north_c
!     rot,h,fsm

! define "u" and "v" masks from "rho" mask
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

! modify SMALL to zero distances to 1 meter
      where( dx < SMALL ) dx = 1.
      where( dy < SMALL ) dy = 1.

!fhx:debug:close north boundary for NWATL
!       if(n_north.eq.-1) then
!          fsm(:,jm)=0.0
!          h(:,jm)=1.0
!       endif
!
!lyo:scs1d:beg:Artificially increase dx & dy for 1d-simulations
      if (n1d.ne.0) then
        do j=1,jm; do i=1,im
          dx(i,j)=1.e9; dy(i,j)=1.e9
        enddo; enddo
      endif
!lyo:scs1d:end:
!
! derive vertical grid variables
!      z(kb) = -1. ! :DIRTY HACK!!! (for stupid grid)
!      zz(kb-1) = .5*(z(kb-1)+z(kb)) ! :
!      zz(kb) = 2.*zz(kb-1)-zz(kb-2) ! :
!!      if (my_task==1) fsm(:,jm-8:jm) = 0.
!      h = 1500.
      do k=1,kb-1
        dz(k) = z(k)- z(k+1)
        dzz(k)=zz(k)-zz(k+1)
      end do
      dz(kb) = dz(kb-1) !=0. !lyo:20110202 =0 is dangerous but checked
      dzz(kb)=dzz(kb-1) !=0. !thro' code - and found to be ok.

! print vertical grid information
      if ( is_master ) then
        print '(/2x,a,7x,a,9x,a,9x,a,9x,a)', 'k','z','zz','dz','dzz'
        do k=1,kb
          print '(1x,i5,4f10.3)', k,z(k),zz(k),dz(k),dzz(k)
        end do
      end if

! set up Coriolis parameter
      cor = 2.*Ohm*sin(north_e*DEG2RAD)

! inertial period for temporal filter
!      period=(2.e0*pi)/abs(cor(im/2,jm/2))/86400.e0

! calculate areas of "rho" cells
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

      if (n_west == -1) then
        do j=1,jm
          aru(1,j)=aru(2,j)
          arv(1,j)=arv(2,j)
        end do
      end if

      if (n_south == -1) then
        do i=1,im
          aru(i,1)=aru(i,2)
          arv(i,1)=arv(i,2)
        end do
      end if

! initialize water column depths
      d  = h+el
      dt = h+et

      return

      end

!______________________________________________________________________
!
      subroutine initial_conditions
!----------------------------------------------------------------------
!  Set up initial and lateral boundary conditions.
!______________________________________________________________________

      use bry
      use glob_bry
      use glob_config, only: read_rst_file, rk
      use glob_domain
      use glob_grid  , only: fsm
      use glob_ocean , only: elb, rmean, rho, s, sb, sclim
     &                     , t, tb, tclim

      implicit none

      integer i,j,k, n
      integer :: ii, jj                !lyo:pac10:
      character(len=120) netcdf_ic_file     !lyo:20110202:
      logical :: fexist                !lyo:20110202:

      read(read_rst_file, '(5x,i2)') n

! read initial temperature and salinity from ic file
!      call read_initial_ts_pnetcdf(kb,tb,sb)
      write(netcdf_ic_file,'(a)') "./in/tsclim/ts_clim.nc"
      inquire(file=trim(netcdf_ic_file),exist=fexist)

      if (fexist) then
!        call read_clim_ts_pnetcdf(tb,sb,n)
        call read_initial_ts_pnetcdf(tb,sb,elb,n)

! read annual-mean, xy-ave t,sclim if avail !lyo:20110202:
        write(netcdf_ic_file,'(a)') "./in/tsclim/ts_mean.nc"
        inquire(file=trim(netcdf_ic_file),exist=fexist)
        if(fexist) then  !annual-mean xy-ave t,sclim
          call read_mean_ts_pnetcdf(tclim,sclim,n)
        else             !annual-mean *clim* - dangerous for rmean
!      call read_clim_ts_pnetcdf_obs(tclim,sclim,rmean,n)
          call read_clim_ts_pnetcdf(tclim,sclim,n)
        endif

      else
        call msg_print("Failed reading clim...", 2, "")
        tb = 15.
        sb = 33.
        tclim = 15.
        sclim = 33.
      end if
! calc. initial density
      call dens(sb,tb,rho)

! calc. rmean; should use xy-ave t,sclim; ok if tsforce is called later
      call dens(sclim,tclim,rmean)

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            t(i,j,k)=tb(i,j,k)
            s(i,j,k)=sb(i,j,k)
          end do
        end do
      end do

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
      do k=1,kb
         do j=1,jm
            do i=1,nfw
               tobw(i,j,k) = tclim( i, j, k) * fsm( i, j)
               sobw(i,j,k) = sclim( i, j, k) * fsm( i, j)
               enddo
            do i=1,nfe
               ii=im-i+1
               tobe(i,j,k) = tclim(ii, j, k) * fsm(ii, j)
               sobe(i,j,k) = sclim(ii, j, k) * fsm(ii, j)
               enddo
               tbw(j,k) = tobw(1,j,k); sbw(j,k) = sobw(1,j,k)
               tbe(j,k) = tobe(1,j,k); sbe(j,k) = sobe(1,j,k)
          enddo
         do i=1,im
            do j=1,nfs
               tobs(i,j,k) = tclim( i, j, k) * fsm( i, j)
               sobs(i,j,k) = sclim( i, j, k) * fsm( i, j)
               enddo
            do j=1,nfn
               jj=jm-j+1
               tobn(i,j,k) = tclim( i,jj, k) * fsm( i,jj)
               sobn(i,j,k) = sclim( i,jj, k) * fsm( i,jj)
               enddo
               tbs(i,k) = tobs(i,1,k); sbs(i,k) = sobs(i,1,k)
               tbn(i,k) = tobn(i,1,k); sbn(i,k) = sobn(i,1,k)
          enddo
      enddo
!lyo:pac10:end:

!      call read_ice_pnetcdf( "ice.19790102.nc", icb )

      return
      end

!______________________________________________________________________
!lyo:pac10:beg:Here thro *end: replaced subr.lateral_boundary_conditions
      subroutine lateral_boundary_conditions
!----------------------------------------------------------------------
!  Read lateral boundary conditions.
!  Transport at eastern boundary for PROFS in GOM.
! ayumi 2010/4/7
!______________________________________________________________________

      use bry        , only: RFE, RFN, RFS, RFW
      use glob_config, only: iperx, ipery, rk
      use glob_domain, only: i_global, im, im_global
     &                     , j_global, jm, jm_global, master_task
      use glob_grid  , only: cor

      implicit none

      logical :: here, judge_inout !lyo:scs1d:
      integer :: i,j,ic,jc         !lyo:scs1d:
      real(kind=rk)    :: corcon            !lyo:scs1d:
!      character(len=120) in_file        !eda:uvforce

      namelist/bry_nml/ rfn, rfe, rfs, rfw

!      if ( n_north == -1 ) eln = elb( :,jm)
!      if ( n_east  == -1 ) ele = elb(im,: )
!      if ( n_south == -1 ) els = elb( :, 1)
!      if ( n_west  == -1 ) elw = elb( 1,: )
!     call read_uabe_pnetcdf(uabe)
!eda:uvforce
!      if (.not. calc_uvforce) then
!        write(in_file,'(a)') "bc.nc"
!        call read_bc_pnetcdf(uabe, uabw, vabs, vabn, in_file, 1)
!      endif

!     Radiation factors for use in subroutine bcond !alu:20101216
!      rfe=0.; rfw=0.; rfn=0.; rfs=0. !=1 Flather; =0 clamped
      rfe=1.; rfw=1.; rfn=1.; rfs=1. !=1 Flather; =0 clamped
      open(73,file='bry.nml',status='old')
      read(73,nml=bry_nml)
      close(73)

! Periodic in "x" and/or "y"?  !lyo:20110224:alu:stcc:
!     iperx.ne.0 if x-periodic; ipery.ne.0 if y-periodic               !
!     iperx(y) < 0 if south/north (west/east) walls are free-slip      !
!     iperx=-1; ipery= 0 !--> x-periodic & free-slip S/N boundaries
!     iperx= 0; ipery= 0 !--> x-periodic & free-slip S/N boundaries
!lyo:scs1d:moved to namelist: pom.nml
!lyo:scs1d:cannot be beta-plane if double-periodic (note cor(*,*)
!     was previously defined in "call read_grid")
      if (iperx.eq.1 .and. ipery.eq.1) then
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

      return
      end

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

      use glob_config, only: rk

      implicit none

      integer, intent(in) :: mw,me,ms,mn,im,jm,nu
      integer, intent(in) :: nw,ne,ns,nn
      real(kind=rk), dimension(im,jm), intent(out) :: frz
      integer :: mmid,i,j,ii,jj
      integer, parameter :: my_task=0
      real(kind=rk) :: c,tanhm
!
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
!
      return
      end
!lyo:pac10:end:

!______________________________________________________________________
!
      subroutine read_tide
!----------------------------------------------------------------------
!fhx:tide:read tidal amplitude & phase at the eastern boundary for PROFS
!______________________________________________________________________

        use glob_config, only: rk

!      call read_tide_east_pnetcdf(ampe,phae)
        call read_tide_east_pnetcdf(ampe,phae,amue,phue)

!      if(my_task.eq.0)print*,ampe(10,1),phae(10,1),amue(10,1),phue(10,1)
!      if(my_task.eq.0)print*,ampe(10,2),phae(10,2),amue(10,2),phue(10,2)
        return

      end
!fhx:tide:read_tide end

!______________________________________________________________________
!
      subroutine update_initial
!----------------------------------------------------------------------
!  Updates the initial conditions and sets the remaining
! initial conditions.
!______________________________________________________________________

      use glob_atmos , only: vfluxf
      use glob_config, only: aam_init, npg, rk
      use glob_const , only: SMALL
      use glob_domain, only: im, is_master, jm, kb, kbm1
      use glob_grid  , only: dz, h
      use glob_ocean , only: aam, d, drhox, drhoy, drx2d, dry2d, dt
     &                     , el, elb, et, etb, etf
     &                     , kh, km, kq, l, q2, q2b, q2l, q2lb
     &                     , s, sb, t, tb, u, ua, ub, uab
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

      call pgscheme(npg)

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
            dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
          end do
        end do
      end do


      return

      end

!______________________________________________________________________
!
      subroutine bottom_friction
!----------------------------------------------------------------------
!  Calculates the bottom friction coefficient.
!______________________________________________________________________

      use glob_config, only: cbcmax, cbcmin, rk, z0b
      use glob_const , only: Kappa
      use glob_domain, only: im, jm, kbm1
      use glob_grid  , only: h, zz
      use glob_ocean , only: cbc

      implicit none

      integer i,j
! calculate bottom friction
      do i=1,im
        do j=1,jm
!lyo:correct:cbc(i,j)=(Kappa/log((1.+zz(kbm1))*h(i,j)/z0b))**2 !lyo:bug:
          cbc(i,j)=(Kappa/log(1.+(1.0+zz(kbm1))*h(i,j)/z0b))**2
          cbc(i,j)=max(cbcmin,cbc(i,j))
! if the following is invoked, then it is probable that the wrong
! choice of z0b or vertical spacing has been made:
          cbc(i,j)=min(cbcmax,cbc(i,j))
        end do
      end do
      return
      end

!______________________________________________________________________
!
      subroutine ztosig(zs,tb,zz,h,t,im,jm,ks,kb,
     $                                   n_west,n_east,n_south,n_north)
!----------------------------------------------------------------------
!  Interpolate vertically.
!______________________________________________________________________

      use glob_config, only: rk

      implicit none

      integer , intent(in ) :: im, jm, ks, kb
      real(rk), intent(in ) :: zs(ks), tb(im,jm,ks), zz(kb)
     &                       , h(im,jm)
      real(rk), intent(out) :: t(im,jm,kb)
      integer , intent(in ) :: n_west,n_east,n_south,n_north

      real(rk) tmax, tin(ks), tout(kb), zzh(kb)
      integer  i,j,k

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

      return

      end

!______________________________________________________________________
!
      subroutine splinc(x,y,n,yp1,ypn,xnew,ynew,m)
!----------------------------------------------------------------------
!  Interpolate using splines.
!______________________________________________________________________

      use glob_config, only: rk

      implicit none

      integer, parameter :: nmax = 300

      integer               , intent(in ) :: n, m
      real(rk)              , intent(in ) :: yp1, ypn
      real(rk), dimension(n), intent(in ) :: x,y
      real(rk), dimension(m), intent(out) :: xnew,ynew

      integer  i, k
      real(rk), dimension(nmax) :: y2,u
      real(rk) p, qn, sig, un

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

      return

      end

!______________________________________________________________________
!
      subroutine splint(xa,ya,y2a,n,x,y)
!----------------------------------------------------------------------
!  Spline interpolation? [ DESCRIPTION NOT PROVIDED ]
!______________________________________________________________________

      use glob_config, only: rk
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
      return

      end

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

      use glob_config, only: rk

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

      return

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

      end function judge_inout

!______________________________________________________________________
!
      subroutine make_grid( TEST_CASE )
!----------------------------------------------------------------------
!  Generate vertical and horizontal grid, topography, areas and masks.
!______________________________________________________________________

      use glob_config, only: rk
      use glob_const , only: DEG2RAD
      use glob_domain
      use glob_grid
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
      do j = 1,jm
        do i = 2,im
          east_v(i,j) = .5*( east_e(i,j) + east_e(i-1,j) )
          north_v(i,j)= .5*( north_e(i,j)+ north_e(i-1,j))
        end do
      end do
      call exchange2d_mpi( east_v, im,jm)
      call exchange2d_mpi(north_v, im,jm)
      if (n_west==-1) then
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


! print vertical grid information
      if ( is_master ) then
        print '(/2x,a,7x,a,9x,a,9x,a,9x,a)', 'k','z','zz','dz','dzz'
        do k=1,kb
          print '(1x,i5,4f10.3)', k,z(k),zz(k),dz(k),dzz(k)
        end do
      end if

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

      return

      end
