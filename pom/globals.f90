module glob_const

  implicit none

  public

! Get precision parameter
  integer, parameter :: rk = 8  ! or 4... 16 is not supported by NetCDF.

  integer, parameter :: PATH_LEN = 120  &
                      ,  VAR_LEN =  40

  integer(1), parameter :: MODE_BAROTROPIC = 2  & ! Barotropic 2-D mode
                         , MODE_BAROCLINIC = 3  & ! Baroclinic 3-D mode
                         , MODE_DIAGNOSTIC = 4    ! Diagnostic 3-D mode with time-invariant temperature and salinity

  integer            &
    MPI_RK             ! real precision for mpi routines

  real(rk)           &
    c2k              & ! Celcius to Kelvin
  , cpw              & ! seawater specific heat
  , deg2rad          & ! pi/180.
  , grav             & ! gravitational acceleration
  , kappa            & ! von Karman's constant
  , ohm              & ! Earth's rotation angular frequency 
  , pi               & ! pi
  , rad2deg          & ! 180./pi
  , rhoref           & ! reference density
  , rho_cpw          & ! rho*Cpw
  , sec2day          & ! 1./86400.
  , small              ! small value

  contains

    subroutine initialize_constants !( constants_override_nml ) : TODO

      c2k     = 273.15        _rk  ! Celcius to Kelvin offset
      cpw     = 3986.         _rk  ! Specific heat of water (J/kg)
      pi      = atan(1._rk)*4._rk  ! PI
      rad2deg = 180._rk/pi         ! radians to degrees conversion factor
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

      call print_config


    end subroutine ! read_domain


    subroutine print_config

      implicit none

      character(48) str


      write( str, '("Size: [",i5," x",i5," x",i4," ]")' ) im_global, jm_global, kb
      call msg_print( "DOMAIN READ", 1, str )


    end subroutine print_config


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
  , iswtch         & ! interval [iint] to switch prtd1 to prtd2 for main output
  , irestart       & ! restart file writing time step [in iint]
  , iouts

  real(rk)         &
    prtd1          & ! output interval (days)
  , prtd2          & ! output interval for SURF.* file (days)
  , swtch          & ! time step interval(days)  to switch from prtd1 to prtd2
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
  , swrad_mean       &
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
    aam_mean          &
  , kh_mean          &
  , km_mean          &
  , rho_mean         &
  , s_mean           &
  , t_mean           &
  , u_mean           &
  , v_mean           &
  , w_mean

  integer &
    num, out_record

  contains

    subroutine allocate_arrays

      use glob_domain, only: im, jm, kb

      implicit none


      num        = 0
      out_record = 1

      swtch      = huge(swtch)

      allocate(             &
        elb_mean   (im,jm)  &
      , swrad_mean (im,jm)  &
      , uab_mean   (im,jm)  &
      , vab_mean   (im,jm)  &
      , wssurf_mean(im,jm)  &
      , wtsurf_mean(im,jm)  &
      , wusurf_mean(im,jm)  &
      , wvsurf_mean(im,jm)  &
      )

         elb_mean = 0.
       swrad_mean = 0.
         uab_mean = 0.
         vab_mean = 0.
      wssurf_mean = 0.
      wtsurf_mean = 0.
      wusurf_mean = 0.
      wvsurf_mean = 0.

      allocate(             &
        aam_mean(im,jm,kb)  &
      , kh_mean (im,jm,kb)  &
      , km_mean (im,jm,kb)  &
      , rho_mean(im,jm,kb)  &
      , s_mean  (im,jm,kb)  &
      , t_mean  (im,jm,kb)  &
      , u_mean  (im,jm,kb)  &
      , v_mean  (im,jm,kb)  &
      , w_mean  (im,jm,kb)  &
      )

      aam_mean = 0.
       kh_mean = 0.
       km_mean = 0.
      rho_mean = 0.
        s_mean = 0.
        t_mean = 0.
        u_mean = 0.
        v_mean = 0.
        w_mean = 0.

      allocate(            &
        elsrf_mean(im,jm)  &
      , usrf_mean (im,jm)  &
      , uwsrf_mean(im,jm)  &
      , vsrf_mean (im,jm)  &
      , vwsrf_mean(im,jm)  &
      )

      elsrf_mean = 0.
       usrf_mean = 0.
      uwsrf_mean = 0.
       vsrf_mean = 0.
      vwsrf_mean = 0.


    end ! subroutine allocate_arrays

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

    subroutine initialize_model_run

      use glob_out, only: iprint, iprints, irestart, iswtch     &
                        , prtd1 , prtd2  , swtch   , write_rst

      implicit none

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
      iswtch  =      nint(      swtch*86400._rk/dti)

      ispi = 1._rk/       real(isplit)
      isp2i= 1._rk/(2._rk*real(isplit))


    end subroutine

    subroutine update_time

      use glob_out, only: iprint, iswtch, prtd2

      implicit none

      time = dti*real(iint)/86400._rk + time0
      if ( iint >= iswtch ) iprint = nint( prtd2*86400._rk/dti )

      if( is_leap( dtime%year ) ) then
        days_in_month(2) = 29
      else
        days_in_month(2) = 28
      end if

      mid_in_month = days_in_month( dtime%month )*43200 !*24*3600/2

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
  , aamfac           & ! horizontal viscosity factor
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

  contains
!______________________________________________________________________
!
    subroutine allocate_arrays
!----------------------------------------------------------------------
!  Allocates module arrays and initializes (some of) them
!______________________________________________________________________
!
      use glob_domain, only: im, im_coarse, jm, jm_coarse, kb, kbm1

      implicit none


      allocate(        &
        aam2d (im,jm)  &
      , aamfac(im,jm)  &
      , advua (im,jm)  &
      , advva (im,jm)  &
      , adx2d (im,jm)  &
      , ady2d (im,jm)  &
      , cbc   (im,jm)  &
      , d     (im,jm)  &
      , drx2d (im,jm)  &
      , dry2d (im,jm)  &
      , dt    (im,jm)  &
      , egb   (im,jm)  &
      , egf   (im,jm)  &
      , el    (im,jm)  &
      , elb   (im,jm)  &
      , elf   (im,jm)  &
      , et    (im,jm)  &
      , etb   (im,jm)  &
      , etf   (im,jm)  &
      , fluxua(im,jm)  &
      , fluxva(im,jm)  &
      , psi   (im,jm)  &
      , ssurf (im,jm)  &
      , tps   (im,jm)  &
      , tsurf (im,jm)  &
      , ua    (im,jm)  &
      , uab   (im,jm)  &
      , uaf   (im,jm)  &
      , utb   (im,jm)  &
      , utf   (im,jm)  &
      , va    (im,jm)  &
      , vab   (im,jm)  &
      , vaf   (im,jm)  &
      , vtb   (im,jm)  &
      , vtf   (im,jm)  &
      , wubot (im,jm)  &
      , wvbot (im,jm)  &
      ! , alon_coarse(im_coarse,jm_coarse) &
      ! , alat_coarse(im_coarse,jm_coarse) &
      ! , mask_coarse(im_coarse,jm_coarse) &
      )

      allocate(          &
        aam  (im,jm,kb)  &
      , advx (im,jm,kb)  &
      , advy (im,jm,kb)  &
      , a    (im,jm,kb)  &
      , c    (im,jm,kb)  &
      , drhox(im,jm,kb)  &
      , drhoy(im,jm,kb)  &
      , dtef (im,jm,kb)  &
      , ee   (im,jm,kb)  &
      , gg   (im,jm,kb)  &
      , hz   (im,jm,kb)  &
      , kh   (im,jm,kb)  &
      , km   (im,jm,kb)  &
      , kq   (im,jm,kb)  &
      , l    (im,jm,kb)  &
      , q2b  (im,jm,kb)  &
      , q2   (im,jm,kb)  &
      , q2lb (im,jm,kb)  &
      , q2l  (im,jm,kb)  &
      , rho  (im,jm,kb)  &
      , sb   (im,jm,kb)  &
      , s    (im,jm,kb)  &
      , tb   (im,jm,kb)  &
      , t    (im,jm,kb)  &
      , ub   (im,jm,kb)  &
      , uf   (im,jm,kb)  &
      , u    (im,jm,kb)  &
      , vb   (im,jm,kb)  &
      , vf   (im,jm,kb)  &
      , v    (im,jm,kb)  &
      , w    (im,jm,kb)  &
      , wr   (im,jm,kb)  &
      , zflux(im,jm,kb)  &
      )

! Initialize arrays
      aamfac   = 1.
      drhox    = 0.
      drhoy    = 0.
      u        = 0.
      uab      = 0.
      uaf      = 0.
      ub       = 0.
      uf       = 0.
      v        = 0.
      vab      = 0.
      vaf      = 0.
      vb       = 0.
      vf       = 0.
      w        = 0.


      end ! subroutine initialize_arrays

end module glob_ocean
