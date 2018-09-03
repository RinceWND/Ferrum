module glob_config

  implicit none

  public

! Get precision parameter
  include 'realkind'

!----------------------------------------------------------------------
! Input variables
!----------------------------------------------------------------------
  integer            &
    mode             & ! calculation mode
  , ntp              & ! water type
  , ispadv           & ! step interval for updating external advective terms
  , isplit           & ! dti/dte
  , iprint           & ! interval in iint at which variables are printed
  , nadv             & ! advection scheme
  , nbct             & ! surface temperature boundary condition
  , nbcs             & ! surface salinity boundary condition
  , nitera           & ! maximum number of iterations for Smolarkiewicz scheme
  , npg              & ! pressure gradient scheme
  , nread_rst        & ! index to start from restart file
  , iperx            & !periodic boundary condition in x direction !alu:stcc
  , ipery            & !periodic boundary condition in y direction !lyo:scs1d:
  , n1d                !n1d .ne. 0 for 1d-simulation !lyo:scs1d:

  real(kind=rk)      &
    alpha            & ! weight for surface slope term in external eq
  , rhoref           & ! reference density
  , sbias            & ! salinity bias
  , slmax            & 
  , tbias            & ! temperature bias
  , tprni            & ! inverse horizontal turbulent Prandtl number
  , umol             & ! background viscosity
  , vmaxl            & ! max vaf used to test for model blow-up
  , write_rst        & 
  , aam_init         & ! initial value of aam
  , cbcmax           & ! maximum bottom friction coefficient
  , cbcmin           & ! minimum bottom friction coefficient
  , days             & ! run duration in days
  , horcon           & ! smagorinsky diffusivity coefficient
  , prtd1            & ! initial print interval (days)
  , prtd2            & ! final print interval (days)
  , smoth            & ! constant to prevent solution splitting
  , sw               & ! smoothing parameter for Smolarkiewicz scheme
  , z0b              & ! bottom roughness
  , lono,lato        & ! lon,lat where aam*(1+fak) larger influence (xs,ys)
  , xs,ys,fak          ! set lono or lato=999.     to skip !lyo:pac10:
! , period           & ! inertial period

  character(len=26)  &
    time_start         ! date and time of start of initial run of model

  character(len=40)  &
    source           & ! TODO: Remove; unused var
  , title

  character(len=120) &
    netcdf_file      & ! output netcdf filename
  , read_rst_file    & ! restart filename
  , write_rst_file     ! output restart filename (TODO: Remove; unused var)


  parameter( lono=999.0,lato=999.0, xs=1.5,ys=1.5, fak=0.5)


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
  , spinup             ! spinup flag - freeze initial forcing


!----------------------------------------------------------------------
! Strings
!----------------------------------------------------------------------
  character(len=4) &
    windf      ! TODO: remove this path

  integer          &
    output_flag    &
  , SURF_flag      &
  , monthly_flag

  real(kind=rk)    &
    sf_bf          & ! boundary flux (barotropic velocities) factor
  , sf_hf          & ! heat flux factor
  , sf_wi            ! wind velocity factor

  real(kind=rk)    &
    t_lo, t_hi

!----------------------------------------------------------------------
! Tide parameters
!----------------------------------------------------------------------
  integer            &
    ntide              ! number of tidal components

  parameter( ntide=1 )

end module glob_config


module glob_domain

  implicit none

  public

!----------------------------------------------------------------------
! Grid size parameters
!----------------------------------------------------------------------
  integer            &
    im_global        & ! number of global grid points in x
  , jm_global        & ! number of global grid points in y
  , kb               & ! number of grid points in z
  , im_local         & ! number of local grid points in x
  , jm_local         & ! number of local grid points in y
  , im_global_coarse & ! number of global grid points in x for coarse grids
  , jm_global_coarse & ! number of global grid points in y for coarse grids
  , im_local_coarse  &
  , jm_local_coarse  &
  , x_division       & ! number of divisions from coarse to fine grids in x
  , y_division       & ! number of divisions from coarse to fine grids in y
  , n_proc             ! number of processors

!______________________________________________________________________
!   Correct values for im_local and jm_local are found using:
!    n_proc=(im_global-2)/(im_local-2)*(jm_global-2)/(jm_local-2)
!   Values higher than necessary will not cause the code to fail, but
! will allocate more memory than is necessary. Values that are too
! low will cause the code to exit.
!   x_divison and y_division can change according to the nest
! requirement; must >=2 for otherwise "if" statements
! in initialize.f won't work

!----------------------------------------------------------------------
! Efective grid size
!----------------------------------------------------------------------
  integer            &
    im               & ! number of grid points used in local x domains
  , imm1             & ! im-1
  , imm2             & ! im-2
  , jm               & ! number of grid points used in local y domains
  , jmm1             & ! jm-1
  , jmm2             & ! jm-2
  , kbm1             & ! kb-1
  , kbm2             & ! kb-2
  , im_coarse        & ! number of coarse grid points used in local x domains
  , jm_coarse          ! number of coarse grid points used in local y domains

! _____________________________________________________________________
!   Note that im and jm may be different between local domains
! im and jm are equal to or lower than im_local and jm_local,
! depending on the use of correct values for im_local and jm_local


!----------------------------------------------------------------------
! Parallel variables
!----------------------------------------------------------------------
  integer           &
    error_status    &
  , my_task         & ! actual parallel processor ID
  , master_task     & ! master processor ID
  , n_west          & ! western parallel processor ID
  , n_east          & ! eastern parallel processor ID
  , n_south         & ! southern parallel processor ID
  , n_north         & ! northern parallel processor ID
  , pom_comm        & ! POM model MPI group communicator
  , pom_comm_coarse   ! satellite data MPI group communicator

  logical is_master   ! flag for master task

  integer               &
      , allocatable     &
      , dimension(:) :: & 
    i_global        & ! global i index for each point in local domain
  , i_global_coarse & ! global i index for each point in local domain
  , j_global        & ! global j index for each point in local domain
  , j_global_coarse   ! global j index for each point in local domain


  contains

!----------------------------------------------------------------------
! Read domain distribution
!----------------------------------------------------------------------
    subroutine read_domain_dist( dist )

      implicit none

      character(len=*), intent(in) :: dist

      namelist/domain_nml/ &
        im_global          &
      , jm_global          &
      , im_local           &
      , jm_local           &
      , im_global_coarse   &
      , jm_global_coarse   &
      , im_local_coarse    &
      , jm_local_coarse    &
      , x_division         &
      , y_division         &
      , kb, n_proc

      x_division = 2
      y_division = 1

      im_global_coarse = 0
      jm_global_coarse = 0
      im_local_coarse = 0
      jm_local_coarse = 0

! read distribution namelist
      open ( 10, file = dist      , status = 'old' )
      read ( 10,  nml = domain_nml )
      close( 10 )

      if ( im_global_coarse == 0 ) im_global_coarse = im_global
      if ( jm_global_coarse == 0 ) jm_global_coarse = jm_global
      if ( im_local_coarse == 0 ) im_local_coarse = im_local
      if ( jm_local_coarse == 0 ) jm_local_coarse = jm_local

! allocate vriables for domain distribution
      allocate(                          &
       i_global(im_local)                &
     , j_global(jm_local)                &
     , i_global_coarse(im_local_coarse)  &
     , j_global_coarse(jm_local_coarse)  &
      )

    end subroutine ! read_domain

end module glob_domain


module glob_time

  use glob_config, only: rk

  implicit none

  public

!----------------------------------------------------------------------
! Internal time-related variables
!----------------------------------------------------------------------
  integer            &
    iint             & ! internal mode step counter
  , iend             & ! total internal mode time steps
  , iext             & ! external mode step counter
  , irestart           ! restart file writing time step (in iint)

  real(kind=rk)      &
    dte              & ! external (2-D) time step (s)
  , dte2             & ! 2*dte
  , dti              & ! internal (3-D) time step (s)
  , dti2             & ! 2*dti
  , ispi             & ! dte/dti
  , isp2i            & ! dte/(2*dti)
  , ramp             & ! inertial ramp
  , time             & ! model time (days)
  , time0              ! initial time (days)

end module glob_time


module glob_const

  use glob_config, only: rk

  implicit none

  public


  real(kind=rk)      &
    deg2rad          & ! pi/180.
  , grav             & ! gravitational acceleration
  , kappa            & ! von Karman's constant
  , ohm              & ! Earth's rotation angular frequency 
  , pi               & ! pi
  , sec2day          & ! 1./86400.
  , small              ! small value

  contains

    subroutine initialize_constants

      pi      = atan(1.)*4.  ! PI
      deg2rad = pi/180.      ! degrees to radians conversion factor
      grav    = 9.806        ! gravity constant (m/s^2)
      small   = 1.e-10       ! small value
      kappa   = 0.4          ! VonKarman's constant
      ohm     = 7.29e-5      ! angular frequency of Earth
      sec2day = 1./86400.    ! seconds to days conversion factor

    end subroutine

end module glob_const


module glob_grid

  use glob_config, only: rk

  implicit none

  public

!----------------------------------------------------------------------
! Vertical discretization
!----------------------------------------------------------------------
  real(kind=rk)              &
       , allocatable         &
       , dimension(:)     :: &
    dz               & ! z(k)-z(k+1)
  , dzz              & ! zz(k)-zz(k+1)
  , z                & ! sigma coordinate from z=0 (surface) to z=-1 (bottom)
  , zz                 ! sigma coordinate, intermediate between z

!----------------------------------------------------------------------
! Horizontal discretization
!----------------------------------------------------------------------
  real(kind=rk)              &
       , allocatable         &
       , dimension(:,:)   :: &
    art              & ! cell area centered on T grid points
  , aru              & ! cell area centered on U grid points
  , arv              & ! cell area centered on V grid points
  , cor              & ! coriolis parameter
  , dum              & ! mask for u velocity
  , dvm              & ! mask for v velocity
  , dx               & ! grid spacing in x
  , dy               & ! grid spacing in y
  , east_c           & ! horizontal coordinate of cell corner points in x
  , east_e           & ! horizontal coordinate of elevation points in x
  , east_u           & ! horizontal coordinate of U points in x
  , east_v           & ! horizontal coordinate of V points in x
  , fsm              & ! mask for scalar variables
  , h                & ! bottom depth
  , north_c          & ! horizontal coordinate of cell corner points in y
  , north_e          & ! horizontal coordinate of elevation points in y
  , north_u          & ! horizontal coordinate of U points in y
  , north_v          & ! horizontal coordinate of V points in y
  , rot                ! rotation angle

end module glob_grid


module glob_ocean

  use glob_config, only: rk

  implicit none

  public

!----------------------------------------------------------------------
! Ocean state 2D arrays
!----------------------------------------------------------------------
  real(kind=rk)              &
       , allocatable         &
       , dimension(:,:)   :: &
    aam2d            & ! vertical average of aam
  , advua            & ! sum of the 2nd, 3rd and 4th terms in eq (18)
  , advva            & ! sum of the 2nd, 3rd and 4th terms in eq (19)
  , adx2d            & ! vertical integral of advx
  , ady2d            & ! vertical integral of advy
  , cbc              & ! bottom friction coefficient
  , d                & ! h+el
  , drx2d            & ! vertical integral of drhox
  , dry2d            & ! vertical integral of drhoy
  , dt               & ! h+et
  , egb              & ! surface elevation use for pressure gradient at time n-1
  , egf              & ! surface elevation use for pressure gradient at time n+1
  , el               & ! surface elevation used in the external mode at time n
  , elb              & ! surface elevation used in the external mode at time n-1
  , elf              & ! surface elevation used in the external mode at time n+1
  , et               & ! surface elevation used in the internal mode at time n
  , etb              & ! surface elevation used in the internal mode at time n-1
  , etf              & ! surface elevation used in the internal mode at time n+1
  , fluxua           & ! [int.use] horizontal x-flux
  , fluxva           & ! [int.use] horizontal y-flux
  , psi              & 
  , ssurf            & 
  , tps              & 
  , tsurf            & 
  , ua               & ! vertical mean of u at time n
  , uab              & ! vertical mean of u at time n-1
  , uaf              & ! vertical mean of u at time n+1
  , utb              & ! ua time averaged over the interval dti at time n-1
  , utf              & ! ua time averaged over the interval dti at time n+1
  , va               & ! vertical mean of v at time n
  , vab              & ! vertical mean of v at time n-1
  , vaf              & ! vertical mean of v at time n+1
  , vtb              & ! va time averaged over the interval dti at time n-1
  , vtf              & ! va time averaged over the interval dti at time n+1
  , wubot            & ! x-momentum flux at the bottom
  , wvbot              ! y-momentum flux at the bottom

!----------------------------------------------------------------------
! Ocean state 3D arrays
!----------------------------------------------------------------------
  real(kind=rk)              &
       , allocatable         &
       , dimension(:,:,:) :: &
    aam              & ! horizontal kinematic viscosity
  , advx             & ! x-horizontal advection and diffusion terms
  , advy             & ! y-horizontal advection and diffusion terms
  , a                & 
  , c                & 
  , drhox            & ! x-component of the internal baroclinic pressure
  , drhoy            & ! y-component of the internal baroclinic pressure
  , dtef             & 
  , ee               & 
  , gg               & 
  , kh               & ! vertical diffusivity
  , km               & ! vertical kinematic viscosity
  , kq               & 
  , l                & ! turbulence length scale
  , q2b              & ! twice the turbulent kinetic energy at time n-1
  , q2               & ! twice the turbulent kinetic energy at time n
  , q2lb             & ! q2 x l at time n-1
  , q2l              & ! q2 x l at time n
  , rho              & ! density
  , rmean            & ! horizontally averaged density
  , sb               & ! salinity at time n-1
  , sclim            & ! horizontally averaged salinity
  , s                & ! salinity at time n
  , tb               & ! temperature at time n-1
  , tclim            & ! horizontally averaged temperature
  , t                & ! temperature at time n
  , ub               & ! horizontal velocity in x at time n-1
  , uf               & ! horizontal velocity in x at time n+1
  , u                & ! horizontal velocity in x at time n
  , vb               & ! horizontal velocity in y at time n-1
  , vf               & ! horizontal velocity in y at time n+1
  , v                & ! horizontal velocity in y at time n
  , w                & ! sigma coordinate vertical velocity
  , wr               & ! real (z coordinate) vertical velocity
  , zflux


end module glob_ocean



module glob_atmos

  use glob_config, only: rk

  implicit none

  public


!----------------------------------------------------------------------
! Air-Sea interaction related arrays
!----------------------------------------------------------------------
  real(kind=rk)              &
       , allocatable         &
       , dimension(:,:)   :: &
    e_atmos          & ! atmospheric pressure
  , swrad            & ! short wave radiation incident on the ocean surface
  , uwsrf            & ! wind speed in x-direction
  , vfluxb           & ! volume flux through water column surface at time n-1
  , vfluxf           & ! volume flux through water column surface at time n+1
  , vwsrf            & ! wind speed in y-direction
  , wssurf           & ! <ws(0)> salinity flux at the surface
  , wtsurf           & ! <wt(0)> temperature flux at the surface
  , wusurf           & ! <wu(0)> momentum flux at the surface
  , wvsurf             ! <wv(0)> momentum flux at the surface


end module glob_atmos


module glob_misc

  use glob_config, only: rk

  implicit none

  public

!----------------------------------------------------------------------
! Scalars
!----------------------------------------------------------------------
  integer          &
    nums           &
  , iprints        &
  , iouts

!----------------------------------------------------------------------
! 2-D arrays
!----------------------------------------------------------------------
  real(kind=rk)              &
       , allocatable         &
       , dimension(:,:)   :: &
    aamfac           & ! aam factor for subroutine incmix
  , icb              & !:sea ice concentration at time n-1
  , ice              & !:sea ice concentration
  , icf              & !:sea ice concentration at time n+1
  , hi               & ! sea ice thickness
  , tauiwu           & ! momentum flux through the ice-water interface
  , tauiwv           & ! momentum flux through the ice-water interface
  , uib              & ! sea ice u-velocity at time n-1
  , ui               & ! sea ice u-velocity
  , uif              & ! sea ice u-velocity at time n+1
  , vib              & ! sea ice v-velocity at time n-1
  , vi               & ! sea ice v-velocity
  , vif              & ! sea ice v-velocity at time n+1
  , alon_coarse      & ! elevation points in x for wind and satellite data
  , alat_coarse      & ! elevation points in y for wind and satellite data
  , mask_coarse        ! mask for scalar variables for wind and satellite data


end module glob_misc


module glob_out

  use glob_config, only: rk

  implicit none

  public

!----------------------------------------------------------------------
! Averages for 2D variables' output
!----------------------------------------------------------------------
  real(kind=rk)              &
       , allocatable         &
       , dimension(:,:)   :: &
! Standard output
    elb_mean         &
  , uab_mean         &
  , vab_mean         &
  , wssurf_mean      &
  , wtsurf_mean      &
  , wusurf_mean      &
  , wvsurf_mean      &
! Surface output
  , elsrf_mean       &
  , usrf_mean        &
  , uwsrf_mean       & ! wind x-velocity output field?
  , vsrf_mean        &
  , vwsrf_mean         ! wind y-velocity output field?


!----------------------------------------------------------------------
! Averages for 3D variables' output
!----------------------------------------------------------------------
  real(kind=rk)              &
       , allocatable         &
       , dimension(:,:,:) :: &
    kh_mean          &
  , km_mean          &
  , rho_mean         &
  , s_mean           &
  , t_mean           &
  , u_mean           &
  , v_mean           &
  , w_mean

  integer &
    num, iout

end module glob_out


module glob_bry

  use glob_config, only: rk

  implicit none

  public

!----------------------------------------------------------------------
! 1D boundary arrays
!----------------------------------------------------------------------
  real(kind=rk)              &
       , allocatable         &
       , dimension(:)     :: &
    cibe      & ! sea ice concentration at the eastern open boundary
  , cibn      & ! sea ice concentration at the northern open boundary
  , cibs      & ! sea ice concentration at the southern open boundary
  , cibw        ! sea ice concentration at the western open boundary

!----------------------------------------------------------------------
! 2D boundary arrays
!----------------------------------------------------------------------
  real(kind=rk)              &
       , allocatable         &
       , dimension(:,:)   :: &
    ampe      & ! M2/K1 eta amplitude at the eastern open boundary
  , amue      & ! M2/K1 UA amplitude at the eastern open boundary
  , phae      & ! M2/K1 eta phase at the eastern open boundary
  , phue        ! M2/K1 UA phase at the eastern open boundary


end module glob_bry
