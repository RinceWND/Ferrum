module glob_const

  implicit none

  public

! Get precision parameter
  include 'realkind'   ! integer, parameter :: rk = 8
                       ! or 4... 16 is not supported by NetCDF.

  integer, parameter :: PATH_LEN = 120  &
                      ,  VAR_LEN =  40

  real(rk)           &
    c2k              & ! Celcius to Kelvin
  , cpw              & ! seawater specific heat
  , deg2rad          & ! pi/180.
  , grav             & ! gravitational acceleration
  , kappa            & ! von Karman's constant
  , ohm              & ! Earth's rotation angular frequency 
  , pi               & ! pi
  , rhoref           & ! reference density
  , rho_cpw          & ! rho*Cpw
  , sec2day          & ! 1./86400.
  , small              ! small value

  contains

    subroutine initialize_constants !( constants_override_nml ) : TODO

      c2k     = 273.16        _rk  ! Celcius to Kelvin offset
      cpw     = 3986.         _rk  ! Specific heat of water (J/kg)
      pi      = atan(1._rk)*4._rk  ! PI
      deg2rad = pi/180._rk         ! degrees to radians conversion factor
      grav    = 9.806         _rk  ! gravity constant (m/s^2)
      small   = 1.e-10        _rk  ! small value
      kappa   = 0.4           _rk  ! VonKarman's constant
      ohm     = 7.29e-5       _rk  ! angular frequency of Earth
      rhoref  = 1025.         _rk  ! recommended values: 1025 for seawater,
                                   !                     1000 for freshwater
      rho_cpw = rhoref*cpw         ! Seawater density times the specific heat of seawater
      sec2day = 1._rk/86400._rk    ! seconds to days conversion factor

    end subroutine

end module glob_const


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
! Effective grid size
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
    i_global            & ! global i index for each point in local domain
  , i_global_coarse     & ! global i index for each point in local domain
  , j_global            & ! global j index for each point in local domain
  , j_global_coarse       ! global j index for each point in local domain


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
      if ( im_local_coarse  == 0 ) im_local_coarse  = im_local
      if ( jm_local_coarse  == 0 ) jm_local_coarse  = jm_local

! allocate variables for domain distribution
      allocate(                          &
       i_global(im_local)                &
     , j_global(jm_local)                &
     , i_global_coarse(im_local_coarse)  &
     , j_global_coarse(jm_local_coarse)  &
      )

    end subroutine ! read_domain

end module glob_domain


module glob_out

  use glob_const, only: rk

  implicit none

  public

!----------------------------------------------------------------------
! Scalars
!----------------------------------------------------------------------
  integer          &
    nums           &
  , iprint         & ! interval [iint] at which variables are printed
  , iprints        & ! interval [iint] for writing SURF?
  , irestart       & ! restart file writing time step [in iint]
  , iouts

  real(rk)         &
    prtd1          & ! output interval (days)
  , prtd2          & ! output interval for SURF.* file (days)
  , write_rst        ! restart output interval [days]

!----------------------------------------------------------------------
! Averages for 2D variables' output
!----------------------------------------------------------------------
  real(rk)                   &
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
  real(rk)                   &
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

!TODO: Move to a separate file
module model_run

  use glob_const , only: rk
  use module_time

  implicit none

  public

!----------------------------------------------------------------------
! Internal time-related variables
!----------------------------------------------------------------------
  integer            &
    iint             & ! internal mode step counter
  , iend             & ! total internal mode time steps
  , iext             & ! external mode step counter
  , isplit           & ! dti/dte
  , mid_in_month     & ! middle of a current month (s)
  , mid_in_nbr       & ! middle of a neighbour month (s)
  , sec_of_month     & ! current second of a month (s)
  , sec_of_year        ! seconds since the beginning of year

  real(rk)           &
    days             & ! run duration in days
  , dte              & ! external (2-D) time step (s)
  , dte2             & ! 2*dte
  , dti              & ! internal (3-D) time step (s)
  , dti2             & ! 2*dti
  , ispi             & ! dte/dti
  , isp2i            & ! dte/(2*dti)
  , ramp             & ! inertial ramp
  , time             & ! model time (days)
  , time0              ! initial time (days)

  character(26)      &
    time_start       & ! date and time of start of initial run of model
  , datetime           ! current model datetime string

  type(date)         &
    dtime              ! current model datetime variable

  integer(2)                &
       , dimension(0:12) :: &  ! number of days in month (0 is December)
    days_in_month


  contains

    subroutine initialize_time( time_string )

      use glob_out, only: iprint, iprints, irestart   &
                        , prtd1 , prtd2  , write_rst

      implicit none

      character(len=*), intent(in) :: time_string

      type(date) dtime_offset

      days_in_month = int( (/31,31,28,31,30,31,30,31,31,30,31,30,31/), 2 )

      ramp = 1.

! Determine internal timestep
      dti  = dte*real(isplit)
      dte2 = dte*2._rk
      dti2 = dti*2._rk

! Define number of steps and output intervals
      iend    = max( nint(       days*86400._rk/dti), 2 )
      iprint  = max( nint(      prtd1*86400._rk/dti), 1 )
      irestart= max( nint(  write_rst*86400._rk/dti), 1 )
      iprints = max( nint(      prtd2*86400._rk/dti), 1 )

      ispi = 1._rk/       real(isplit)
      isp2i= 1._rk/(2._rk*real(isplit))

! Initialize time
! Do not offset time if not restarting
      if ( time_string /= '' ) then

        dtime_offset = str2date(  time_start(1:19) )
        dtime        = str2date( time_string(1:19) )

        time0 = real( (dtime-dtime_offset)/86400 ) ! Why should this be integer division?

      else

        dtime = str2date( time_start )

        time0 = 0._rk

      end if

      time = time0

      call update_time


    end subroutine

    subroutine update_time

      implicit none

      time = dti*real(iint)/86400._rk + time0

      if( is_leap( dtime%year ) ) then
        days_in_month(2) = 29
      else
        days_in_month(2) = 28
      end if

      mid_in_month = days_in_month(dtime%month  )*43200 !*24*3600/2

      sec_of_month = (dtime%day-1)*24*3600      &
                    + dtime%hour     *3600      &
                    + dtime%min*60 + dtime%sec

      if ( sec_of_month <= mid_in_month ) then
        mid_in_nbr = days_in_month(     dtime%month-1      )*43200
      else
        mid_in_nbr = days_in_month(mod( dtime%month+1, 12 ))*43200
      end if

      sec_of_year = sec_of_month                                   &
                  + ( sum(days_in_month(1:dtime%month-1)) )*86400
!______________________________________________________________________
!  var(1:0) subset is expected to return 0.
! gfortran does this. Not sure about ifort or others.

!      print *, sec_of_year, ":", sec_of_month, "|", mid_in_month
!      sec_of_year = seconds_of_year(dtime) ! TODO: is this faster?

!      ramp = 1.

!      if(lramp) then
!        ramp=time/period
!        if(ramp.gt.1.e0) ramp=1.e0
!      else
!        ramp=1.e0
!      endif

    end ! subroutine


end module model_run


module glob_grid

  use glob_const, only: rk

  implicit none

  public

!----------------------------------------------------------------------
! Vertical discretisation
!----------------------------------------------------------------------
  real(rk)                   &
       , allocatable         &
       , dimension(:)     :: &
    dz               & ! z(k)-z(k+1)
  , dzz              & ! zz(k)-zz(k+1)
  , z                & ! sigma coordinate from z=0 (surface) to z=-1 (bottom)
  , zz                 ! sigma coordinate, intermediate between z

!----------------------------------------------------------------------
! Horizontal discretisation
!----------------------------------------------------------------------
  real(rk)                   &
       , allocatable         &
       , dimension(:,:)   :: &
    art              & ! cell area centered on T grid points
  , aru              & ! cell area centered on U grid points
  , arv              & ! cell area centered on V grid points
  , cor              & ! Coriolis parameter
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

  use glob_const, only: rk

  implicit none

  public

!----------------------------------------------------------------------
! Ocean state 2D arrays
!----------------------------------------------------------------------
  real(rk)                   &
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
  , egb              & ! surface elevation used for pressure gradient at time n-1
  , egf              & ! surface elevation used for pressure gradient at time n+1
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
  real(rk)                   &
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
  , hz               & ! model cell depth [m]
  , kh               & ! vertical diffusivity
  , km               & ! vertical kinematic viscosity
  , kq               & 
  , l                & ! turbulence length scale
  , q2b              & ! twice the turbulent kinetic energy at time n-1
  , q2               & ! twice the turbulent kinetic energy at time n
  , q2lb             & ! q2 x l at time n-1
  , q2l              & ! q2 x l at time n
  , rho              & ! density
  , sb               & ! salinity at time n-1
  , s                & ! salinity at time n
  , tb               & ! temperature at time n-1
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


module glob_misc

  use glob_const, only: rk

  implicit none

  public

!----------------------------------------------------------------------
! 2-D arrays
!----------------------------------------------------------------------
  real(rk)                   &
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


module glob_bry

  use glob_const, only: rk

  implicit none

  public

!----------------------------------------------------------------------
! 1D boundary arrays
!----------------------------------------------------------------------
  real(rk)                   &
       , allocatable         &
       , dimension(:)     :: &
    cibe      & ! sea ice concentration at the eastern open boundary
  , cibn      & ! sea ice concentration at the northern open boundary
  , cibs      & ! sea ice concentration at the southern open boundary
  , cibw        ! sea ice concentration at the western open boundary

!----------------------------------------------------------------------
! 2D boundary arrays
!----------------------------------------------------------------------
  real(rk)                   &
       , allocatable         &
       , dimension(:,:)   :: &
    ampe      & ! M2/K1 eta amplitude at the eastern open boundary
  , amue      & ! M2/K1 UA amplitude at the eastern open boundary
  , phae      & ! M2/K1 eta phase at the eastern open boundary
  , phue        ! M2/K1 UA phase at the eastern open boundary


end module glob_bry
