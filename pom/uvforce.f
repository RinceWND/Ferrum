      module uvforce

      implicit none

      private

      public :: uvforce_init, uvforce_main 

      include 'pom.h'

!     days in month
      integer :: mday(0:12) = (/31, 31, 28, 31, 30, 31, 30,               
     $                          31, 31, 30, 31, 30, 31/)
                             

!      real(kind=rk), dimension( im_local ) ::
!     $     vabs_a, vabn_a, vabs_b, vabn_b
!
!      real(kind=rk), dimension( jm_local ) ::
!     $     uabe_a, uabw_a, uabe_b, uabw_b

      real(kind=rk), dimension( im_local, kb ) ::
     $     vbs_a, vbn_a, vbs_b, vbn_b

      real(kind=rk), dimension( jm_local, kb ) ::
     $     ube_a, ubw_a, ube_b, ubw_b

    
      integer :: mon_a, mon_b, sec_in_month, mid_in_month
      integer :: k
      character(len=7) :: infile_a, infile_b
      real(kind=rk) :: aa


      contains

!==============================================================
! Initialization variables for U & V boundary condition 
!--------------------------------------------------------------
      subroutine uvforce_init( d_in )

      use module_time

      implicit none

!     intent(in)
      type(date), intent(in) :: d_in 


!     check leap year
      if( ( mod( d_in%year, 4 ) .eq. 0 
     $     .and. mod( d_in%year, 100 ) .ne. 0 ) 
     $     .or. mod( d_in%year, 400 ) .eq. 0 ) then
         mday(2) = 29
      else
         mday(2) = 28
      endif
                      
!     current time [sec] from the beginning of the month.
      
      sec_in_month = d_in%day * 24 * 3600  
     $             + d_in%hour * 3600 + d_in%min * 60 + d_in%sec


!     mid-point [sec] in the month.
      
      mid_in_month = int( real( mday( d_in%month ) )/ 2.0  
     $             * 24. * 3600. )


!     decide between which two months.
 
!     former half in the month.    
      if ( sec_in_month .le. mid_in_month ) then
         mon_a = d_in%month - 1
!     latter half in the month
      else
         mon_a = d_in%month 
      endif
      
      mon_b = mon_a + 1

      if ( mon_a ==  0 ) mon_a = 12
      if ( mon_b == 13 ) mon_b =  1


!     data open.
         write( infile_a, '( "bc",i2.2,".nc" )' ) mon_a
         call read_bc_pnetcdf
     $        ( ube_a, ubw_a, vbs_a, vbn_a, "bc.nc", mon_a )

!         write( infile_b, '( "bc",i2.2,".nc" )' ) mon_b
         call read_bc_pnetcdf
     $        ( ube_b, ubw_b, vbs_b, vbn_b, "bc.nc", mon_b )
     
!         call read_bc_pnetcdf_obs
!     $        ( ube_a, ubw_a, vbs_a, vbn_a, "bry_old.nc", mon_a )
!
!         call read_bc_pnetcdf_obs
!     $        ( ube_b, ubw_b, vbs_b, vbn_b, "bry_old.nc", mon_b )
!     TODO: Treat water transport increase properly in case of baroclinic velocities input
!           The below won't work for vertically non-homogeneous velocity field
         ube_a = sf_bf*ube_a
         ubw_a = sf_bf*ubw_a
         vbs_a = sf_bf*vbs_a
         vbn_a = sf_bf*vbn_a

      if ( my_task == master_task ) 
     $        write(*,'(/a/)') "---------- uvforce_init."


      return

      end subroutine uvforce_init
!--------------------------------------------------------------

!==============================================================
! Read & time-interpolation of U & V boundary condition
!--------------------------------------------------------------
      subroutine uvforce_main( d_in )

      use module_time

      implicit none
      
      ! intent(in)
      type(date), intent(in) :: d_in
      


!     check leap year. 
      if( ( mod( d_in%year, 4 ) .eq. 0 
     $     .and. mod( d_in%year, 100 ) .ne. 0 ) 
     $     .or. mod( d_in%year, 400 ) .eq. 0 ) then
         mday(2) = 29
      else
         mday(2) = 28
      endif

                      
!     current time [sec] from the beginning of the month.
      
      sec_in_month = d_in%day * 24 * 3600  
     $             + d_in%hour * 3600 + d_in%min * 60 + d_in%sec


!     mid-point [sec] in the month.
      
      mid_in_month = int( real( mday( d_in%month ) )/ 2. * 24.*3600. )


!     decide between which two months.
 
!     former half in the month.    
      if ( sec_in_month .le. mid_in_month ) then
         mon_a = d_in%month - 1
         aa = real( sec_in_month ) 
     $      + real( mday( mon_a ) )/ 2. * 24. * 3600.

!     latter half in the month
      else
         mon_a = d_in%month 
         aa = real( sec_in_month ) 
     $      - real( mday( mon_a ) )/ 2. * 24. * 3600.
      endif
      
      mon_b = mon_a + 1

      if ( mon_a ==  0 ) mon_a = 12
      if ( mon_b == 13 ) mon_b =  1

      aa = aa / 
     $     ( real( mday( mon_a ) + mday( mon_b ) ) / 2. * 24. * 3600. )            



!     data open.

      if ( sec_in_month > mid_in_month  .and.
     $     sec_in_month - mid_in_month <= int( dti ) ) then

!         uabe_a = uabe_b
!         uabw_a = uabw_b
!         vabs_a = vabs_b
!         vabn_a = vabn_b
         
         ube_a = ube_b
         ubw_a = ubw_b
         vbs_a = vbs_b
         vbn_a = vbn_b


         write( infile_b, '( "bc",i2.2,".nc" )' ) mon_b
         call read_bc_pnetcdf
     $        ( ube_b, ubw_b, vbs_b, vbn_b, "bc.nc", mon_b )
     
         ube_b = sf_bf*ube_b
         ubw_b = sf_bf*ubw_b
         vbs_b = sf_bf*vbs_b
         vbn_b = sf_bf*vbn_b

!         call read_bc_pnetcdf_obs
!     $        ( ube_b, ubw_b, vbs_b, vbn_b, "bry_old.nc", mon_b )


      endif
     

!     time interpolation.

      ube = ( 1.0 - aa ) * ube_a + aa * ube_b
      ubw = ( 1.0 - aa ) * ubw_a + aa * ubw_b
      vbs = ( 1.0 - aa ) * vbs_a + aa * vbs_b
      vbn = ( 1.0 - aa ) * vbn_a + aa * vbn_b
      
!     integrate vertically
      uabe = 0.
      uabw = 0.
      vabs = 0.
      vabn = 0.
      do k=1,kb
        uabe = uabe + ube(:,k)*dz(k)
        uabw = uabw + ubw(:,k)*dz(k)
        vabs = vabs + vbs(:,k)*dz(k)
        vabn = vabn + vbn(:,k)*dz(k)
      end do
      ! TODO: Check barotropic velocities!!!

      return

      end subroutine uvforce_main
!--------------------------------------------------------------

      end module uvforce
