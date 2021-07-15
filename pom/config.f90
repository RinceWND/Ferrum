module config

  use glob_const, only: PATH_LEN, RK, VAR_LEN

  implicit none

  public

!----------------------------------------------------------------------
! Input variables
!----------------------------------------------------------------------
  logical            &
    do_restart       & ! flag to start from restart file
  , append_output      ! write to an existing output file

  integer            &
    mode             & ! calculation mode
  , ntp              & ! water type
  , ispadv           & ! step interval for updating external advective terms
  , nadv             & ! advection scheme
  , nbct             & ! surface temperature boundary condition
  , nbcs             & ! surface salinity boundary condition
  , nitera           & ! maximum number of iterations for Smolarkiewicz scheme
  , npg              & ! pressure gradient scheme
  , iperx            & !periodic boundary condition in x direction !alu:stcc
  , ipery            & !periodic boundary condition in y direction !lyo:scs1d:
  , n1d                !n1d .ne. 0 for 1d-simulation !lyo:scs1d:

  real(kind=rk)      &
    alpha            & ! weight for surface slope term in external eq
  , sbias            & ! salinity bias
  , slmax            & 
  , tbias            & ! temperature bias
  , tprni            & ! inverse horizontal turbulent Prandtl number
  , umol             & ! background viscosity
  , vmaxl            & ! max vaf used to test for model blow-up
  , aam_init         & ! initial value of aam
  , cbcmax           & ! maximum bottom friction coefficient
  , cbcmin           & ! minimum bottom friction coefficient
  , horcon           & ! smagorinsky diffusivity coefficient
  , smoth            & ! constant to prevent solution splitting
  , sw               & ! smoothing parameter for Smolarkiewicz scheme
  , z0b              & ! bottom roughness
  , lono,lato        & ! lon,lat where aam*(1+fak) larger influence (xs,ys)
  , xs,ys,fak          ! set lono or lato=999.     to skip !lyo:pac10:
! , period           & ! inertial period

  character(VAR_LEN)  &
    title

  character(PATH_LEN) &
    netcdf_file       & ! output netcdf filename
  , initial_file        ! restart filename to read from

  procedure(PROC_INTERFACE), pointer :: pressure_gradient, advection

  abstract interface
    subroutine PROC_INTERFACE
    end subroutine PROC_INTERFACE
  end interface

!  parameter( lono=999.0,lato=999.0, xs=1.5,ys=1.5, fak=0.5)


!----------------------------------------------------------------------
! Module switches
!----------------------------------------------------------------------
  logical            &
    calc_bulk        & ! use bulk formulations for surface fluxes
  , calc_wind        & ! use surface wind stress BC
  , calc_tsforce     & ! use surface heat flux BC
  , calc_river       & ! use river discharge       [ NOT TESTED ]
  , calc_assim       & ! use SST assimilation      [ NOT TESTED ]
  , calc_assimdrf    & ! use drifter assimilation  [ NOT TESTED ]
  , calc_interp      & ! use interpolation (wind?) [ NOT TESTED ]
  , calc_tsurf_mc    & ! use MCSST                 [ NOT TESTED ]
  , calc_tide        & ! use tides                 [ NOT TESTED ]
  , calc_uvforce     & ! use lateral BC for velocities
  , calc_ice         & ! use simple ice drag model [ TODO: TEST ]
  , spinup           & ! spinup flag - freeze initial forcing
  , use_air          & ! apply surface forcing
  , use_bry          & ! apply active lateral forcing
  , use_ice          & ! use simple ice drag model [ TODO: TEST ]
  , use_river        & ! use river discharge
  , use_tide           ! use tides


!----------------------------------------------------------------------
! Strings
!----------------------------------------------------------------------
  character(4) &
    windf      ! TODO: remove this path

  logical          &
    output_means     ! write variables averaged over time strting from previous output

  integer          &
    output_flag    &
  , SURF_flag      & ! TODO: move to potential output module?.. (glob_out rework)
  , monthly_flag

  real(rk)         &
    sf_bf          & ! boundary flux (barotropic velocities) factor
  , sf_hf          & ! heat flux factor
  , sf_wi            ! wind velocity factor

  real(rk)         &
    s_hi, s_lo     & ! thresholds for salinity and temperature
  , t_hi, t_lo

!----------------------------------------------------------------------
! Tide parameters
!----------------------------------------------------------------------
  integer            &
    ntide              ! number of tidal components

  parameter( ntide=1 )


  contains

!______________________________________________________________________
!
    subroutine parameters_init
!----------------------------------------------------------------------
!  Initializes parameters with their default values.
!______________________________________________________________________
!
      use model_run, only: days, dte, isplit, time_start

      implicit none

! Input of filenames and constants

! Title of model run:
      title = "Numerical experiment"

! Output filename: [ TODO: handle it in io module? ]
      netcdf_file = "results"

! Model calculation mode:
!   0 : No calculation [ NOT IMPLEMENTED ]
!   2 : Barotropic 2-D mode
!   3 : Baroclinic 3-D mode (default)
!   4 : Diagnostic mode with T and S held fixed
      mode = 3

! Advection scheme: [ TODO: Move to solver module (TODO: create solver module) ]
!   1 : Simple upstream scheme (???)
!   2 : Smolarkiewicz scheme
      nadv = 2

! Number of maximum iterations for Smolarkiewicz scheme.
! Usual "safe" number is 4.
!   In the current implementation the scheme convergence is checked
!  against |eps|, so it may take fewer steps than the set maximum:
      nitera = 10

! Smoothing factor for Smolarkiewicz scheme [ TODO: move to a separate module? ]
      sw = .75

! Pressure gradient scheme:
!   1 : Simple POM scheme
!   2 : 4th order McCalpin (default)
!   3 : Lin         [ TESTS NEEDED ]
!   4 : Song        [ TESTS NEEDED ]
!   5 : Shchepetkin [ TESTS NEEDED ]
      npg = 2

! Logical for inertial ramp (.true. if inertial ramp to be applied
! to wind stress and baroclinic forcing, otherwise .false.)
!      lramp=.false.

! Temperature bias (deg. C)
      tbias = 0.

! Salinity bias
      sbias = 0.

! Bottom roughness (metres)
      z0b = .01

! Minimum bottom friction coeff.
      cbcmin = .0025

! Maximum bottom friction coeff.
      cbcmax = 5. !1.e0   !lyo:20110315:botwavedrag:

! Smagorinsky diffusivity coeff.
      horcon = .2  !lyo:pac10:exp004: !=0.1e0

! Inverse horizontal turbulent Prandtl number (ah/am; dimensionless):
! NOTE that tprni=0.e0 yields zero horizontal diffusivity!
      tprni = .2

! Background viscosity used in subroutines profq, proft, profu and
! profv (S.I. units):
      umol = 2.e-5

! Maximum magnitude of vaf (used in check that essentially tests
! for CFL violation):
      vmaxl = 5. !lyo:debug:100.e0

! Maximum allowable value of:
!   <difference of depths>/<sum of depths>
! for two adjacent cells (dimensionless). This is used in subroutine
! slpmax. If >= 1, then slpmax is not applied:
      slmax = 2.

! Water type, used in subroutine proft.
!    ntp    Jerlov water type
!     1            i
!     2            ia
!     3            ib
!     4            ii
!     5            iii
      ntp = 2

! Surface temperature boundary condition, used in subroutine proft:
!    nbct   prescribed    prescribed   short wave
!           temperature      flux      penetration
!     1        no           yes           no
!     2        no           yes           yes
!     3        yes          no            no
!     4        yes          no            yes
      nbct = 1

! Surface salinity boundary condition, used in subroutine proft:
!    nbcs   prescribed    prescribed
!            salinity      flux
!     1        no           yes
!     3        yes          no
! NOTE that only 1 and 3 are allowed for salinity.
      nbcs = 1

! Step interval during which external (2-D) mode advective terms are
! not updated (dimensionless):
      ispadv = 5

! Default mode splitting step value:
      isplit = 30

! Constant in temporal filter used to prevent solution splitting
! (dimensionless):
      smoth = .10

! Starting datetime string:
      time_start = "0001-01-01 00:00:00"

! Run duration (days):
      days = 0.

! External time step (sec):
      dte = 60.

! Weight used for surface slope term in external (2-D) dynamic
! equation (a value of alpha = 0.e0 is perfectly acceptable, but the
! value, alpha=.225e0 permits a longer time step):
      alpha = .225

! Initial value of aam:
      aam_init = 500.

! Thresholds for temperature and salinity
      t_hi =  999.
      t_lo =   -2. ! allow water to overcool a bit
      s_hi =  999.
      s_lo =    0.

! Mixing point modifiers
      fak  =    .5
      lato = 999.
      lono = 999.
      xs   =   1.5
      ys   =   1.5

! Module switches
      USE_AIR   = .false.
      USE_BRY   = .false.
      USE_ICE   = .false.
      USE_RIVER = .false.
      USE_TIDE  = .false.

! Restart flags
      do_restart    = .false.
      append_output = .false.

      output_flag  = 1
      output_means = .false.

    end ! subroutine
!
!______________________________________________________________________
!
    subroutine get_pgproc(id, proc)
!----------------------------------------------------------------------
!  Selects specified routine for pressure gradient.
!______________________________________________________________________
!
      implicit none

      procedure(), pointer, intent(out) :: proc
      integer             , intent(in ) :: id

      external :: baropg, baropg_mcc, baropg_lin, baropg_song_std, baropg_shch


      select case (id)
        case (1)
          proc => baropg
        case (2)
          proc => baropg_mcc
        case (3)
          proc => baropg_lin
        case (4)
          proc => baropg_song_std
        case (5)
          proc => baropg_shch
        case default
          proc => baropg_mcc
      end select


    end subroutine
!
!______________________________________________________________________
!
    subroutine read_config( config )
!----------------------------------------------------------------------
!  Reads input values and defines constants.
!______________________________________________________________________
!
      use glob_const , only: rk
      use glob_out   , only: prtd1, prtd2, swtch , write_rst
      use model_run  , only: days , dte  , isplit, initialize_model_run  &
                           , time_start

      implicit none


      character(*), intent(in) :: config

!      include 'io.h'
!      include 'bulk.h'

      namelist/run_nml/                          &
        days      , dte  , do_restart  , isplit  &
      , time_start, title, initial_file, spinup

      namelist/setup_nml/                     &
        aam_init, alpha , cbcmax, cbcmin      &
      , fak     , horcon, ispadv, lato        &
      , lono    , mode  , nadv  , nbcs        &
      , nbct    , nitera, npg   , ntp         &
      , s_hi    , s_lo  , sbias , smoth       &
      , sw      , t_hi  , t_lo  , tbias       &
      , tprni   , umol  , vmaxl , xs          &
      , ys      , z0b

      namelist/output_nml/                                     &
        append_output, monthly_flag, netcdf_file, output_flag  &
      , output_means , prtd1       , prtd2      , swtch        &
      , write_rst    , SURF_flag

      namelist/modules_nml/                             &
        USE_AIR, USE_BRY, USE_ICE, USE_RIVER, USE_TIDE

      namelist/sensitivity_nml/ sf_bf, sf_hf, sf_wi


!  Initialize parameters first.
      call parameters_init

!  Read input namelist.
      open ( 73, file = config, status = 'old' )
      read ( 73, nml =     run_nml )
      rewind( 73 )
      read ( 73, nml =  output_nml )
      rewind( 73 )
      read ( 73, nml =   setup_nml )
      rewind( 73 )
      read ( 73, nml = modules_nml )
      close( 73 )

!  Initialise timestep-related variables
      call initialize_model_run

! Set pressure gradient and advection schemes procedures
      call get_pgproc (npg ,pressure_gradient)

! End of input of constants
      call print_config


    end ! subroutine read_config
!
!______________________________________________________________________
!
    subroutine print_config
!----------------------------------------------------------------------
!  Prints model configuration.
!______________________________________________________________________
!
      use glob_const
      use glob_domain, only: is_master
      use glob_out   , only: iprint, irestart, prtd1, prtd2, write_rst
      use model_run

      implicit none


      if ( .not.is_master ) return

! Print initial (core) summary
      call msg_print("CONFIGURATION READ", 1, "")
      print '(/a)', "---- EXPERIMENT ---------------------------------"
      print '( "  ",a40)', title
      print '(a/)', "-------------------------------------------------"
      print '(/a)', "---- CONFIGURATION ------------------------------"
      print '(/" MODE                : ",i10)' , mode
      print '(/" Experiment duration : ",f10.4)', days
      print '(" Start time          : ",a26)'  , time_start
      print '(/" External timestep   : ",f10.2)', dte
      print '(" Internal timestep   : ",f10.1)', dti
      print '(" INT/EXT ratio       : ",i10)'  , isplit
      print '(" Total model steps   : ",i10)'  , iend
      print '(/" Advection scheme    : ",i10)' , nadv
      print '(" Max iterations      : ",i10)'  , nitera
      print '(" Smolar. smoothing   : ",f10.4)', sw
      print '(/" PGScheme            : ",i10)'  , npg
      print '(" Output interval     : ",f10.4)', prtd1
      print '(" Output steps        : ",i10)'  , iprint
      print '(" AUX. output steps   : ",f10.4)', prtd2
      print '(/" IS RESTART?         : ",l2)'   , do_restart
      print '(" Restart interval    : ",f10.4)', write_rst
      print '(" Restart steps       : ",i10)'  , irestart
      print '(/" Temperature bias    : ",f10.3)', tbias
      print '(" Salinity bias       : ",f10.3)', sbias
      print '(/" Bottom roughness    : ",f10.6)', z0b
      print '(" Min bottom friction : ",f10.6)', cbcmin
      print '(" Max bottom friction : ",f10.6)', cbcmax
      print '(/" Smag. diffusivity   : ",f10.3)', horcon
      print '(" Turb.Pr.inversed    : ",f10.4)', tprni
      print '(" Background visc.    : ",f10.4)', umol
      print '(/" Velocity limit      : ",f10.4)', vmaxl
      print '(" Slope factor[UNUSED]: ",f10.4)', slmax
!      print '('' lono,lato  = '',2f10.4)', lono,lato
!      print '('' xs,ys,fak  = '',3f10.4)', xs,ys,fak
!      print '('' iperx,ipery= '',2i10)', iperx,ipery
      print '(/" Water type          : ",i10)'  , ntp
      print '(" Surface temp BC     : ",i10)'  , nbct
      print '(" Surface salt BC     : ",i10)'  , nbcs
      print '(/" Int.adv. skip steps : ",i10)'  , ispadv
      print '(" Temporal filter     : ",f10.4)', smoth
      print '(" Surf. slope weight  : ",f10.4)', alpha
!      print '('' lramp      = '',l10)', lramp
      print '(/" calc_wind           : ",l2)', calc_wind
      print '(" calc_tsforce        : ",l2)', calc_tsforce
      print '(" calc_river          : ",l2)', calc_river
      print '(" calc_assim          : ",l2)', calc_assim
      print '(" calc_assimdrf       : ",l2)', calc_assimdrf !eda
      print '(" calc_uvforce        : ",l2)', calc_uvforce !eda:uvforce
      print '(" calc_tsurf_mc       : ",l2)', calc_tsurf_mc !fhx:mcsst
      print '(" calc_tide           : ",l2)', calc_tide     !fhx:tide
      print '(" calc_interp         : ",l2)', calc_interp   !fhx:interp_flag
      print '(" output_flag         : ",i2)', output_flag  !fhx:20110131:
      print '(" SURF_flag           : ",i2)', SURF_flag    !fhx:20110131:
! TODO: move below output to their respected modules.
!        print '(/'' Sensitivity:'')'
!        print '(''   heat flux     = '',f10.4)', sf_hf
!        print '(''   wind speed    = '',f10.4)', sf_wi
!        print '(''   lat.velocities= '',f10.4)', sf_bf
!        if ( calc_bulk ) then
!          print '(/'' Air-Sea Bulk: [ ENABLED ]'')'
!          print '(''   calc swrad    = '',l2)', calc_swr
!          print '(''   use COARE     = '',l2)', use_coare
!          print '(''   longwave form.= '',i2)', lwrad_formula
!        else
!          print '(/'' Air-Sea Bulk: [ DISABLED ]'')'
!        end if
!      end if


    end ! subroutine print_config


end module config
