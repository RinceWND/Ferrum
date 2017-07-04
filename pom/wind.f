      module wind

      implicit none

      private

      public :: wind_init, wind_main

      include 'pom.h'

      real(kind=rk), dimension( im_local, jm_local ) :: 
     $  uwnd_a, vwnd_a, uwnd_b, vwnd_b, uwnd_fine, vwnd_fine

!lyo:pac10:more efficient:Comment out the followings
!     real(kind=rk), dimension( im_local, jm_local ) :: !fhx:20110131:
!    $ uwind_surf, vwind_surf

      real(kind=rk), dimension( im_local_coarse, jm_local_coarse ) :: 
     $  uwnd_a_coarse, vwnd_a_coarse, uwnd_b_coarse, vwnd_b_coarse,
     $  uwnd_coarse, vwnd_coarse

      real(kind=rk), dimension( im_local_coarse, jm_local_coarse, 4 ) :: 
     $  uwnd_buf_coarse, vwnd_buf_coarse

      real(kind=rk), dimension( im_local, jm_local, 4 ) :: 
     $  uwnd_buf, vwnd_buf 

      integer :: nsec(4)=(/ 0, 6*3600, 12*3600, 18*3600 /)
      integer :: sec_in_day, i, j, k
      
      real(kind=rk), parameter :: rhoa = 1.22, rhow = 1025.0
      real(kind=rk) :: aa, uwnd, vwnd, cda, uvabs!, rdisp !lyo:pac10:add rdisp

      character*16 :: infile


      contains

!==================================================================
! Initialize variables for wind 
!------------------------------------------------------------------
      subroutine wind_init( d_in )
      
      use module_time
      use interp
      
      implicit none

      ! intent(in)
      type(date), intent(in) :: d_in

      ! local     
      type(date) :: d_tmp, d_tmp2

      logical :: lexist

      integer n

      ! initialize
      n = 1

      wusurf = 0.0
      wvsurf = 0.0
      

      ! wind stresses in initial state

      if ( calc_wind ) then

         sec_in_day = d_in%hour*3600 + d_in%min*60 + d_in%sec

!!*      write( infile, '( //trim(windf)//"_",i4.4,2i2.2,".nc" )' )
!!*  $        d_in%year, d_in%month, d_in%day     
!         write( infile, '( a4,"_",i4.4,2i2.2,".nc" )' )
!     $        windf, d_in%year, d_in%month, d_in%day
         write( infile, '( a3,".",i4.4,".nc" )' )
     $        "mfl", d_in%year
!!       write( infile, '( "gfsw_",i4.4,2i2.2,".nc" )' )   ! fhx:read gfsw wind
!!    $        d_in%year, d_in%month, d_in%day     

! fhx: check wind data exists,is ~exist, set wind to be 0.10/26/2010  
         inquire(file='in/'//trim(windf)//'/'//trim(infile),
     $           exist=lexist)   
!!       inquire(file='in/gfsw/'//trim(infile),exist=lexist)   

         if(lexist) then
         d_tmp = str2date("1979-01-01 00:00:00")
         d_tmp%year = d_in%year
         n = int(dif_date(d_in, d_tmp)/86400.)*4+1

         call read_wind_pnetcdfc
     $             ( uwnd_buf_coarse, vwnd_buf_coarse, trim(infile), n )

         else
             if ( my_task == master_task ) then
                  write(*,'(/2a)') 
     $              "missing wind data at wind_init : "
     $              , trim(infile)
             endif
             uwnd_buf_coarse = 0.
             vwnd_buf_coarse = 0.         
            
         endif

         write(*,*) minval(uwnd_buf_coarse),maxval(uwnd_buf_coarse)

! interpolation
      if(calc_interp) then ! fhx:interp_flag:add flag for interp fgrid.  

       do k=1,4  
        uwnd_coarse = uwnd_buf_coarse( :, :, k )
        vwnd_coarse = vwnd_buf_coarse( :, :, k )

        call interp_mask_2d(uwnd_coarse,1,east_e,north_e,fsm,uwnd_fine)
        call interp_mask_2d(vwnd_coarse,2,east_e,north_e,fsm,vwnd_fine)
!fhx: after interpolation, i=1,3 and j=1,3 seem to have problems
      if (n_west.eq.-1) then
!         print*, uwnd_fine(1,1),uwnd_fine(1,2),uwnd_fine(2,1)
!         print*, vwnd_fine(1,1),vwnd_fine(1,2),vwnd_fine(2,1)
!         print*, ''
       do i=1,2
         uwnd_fine(i,:)=uwnd_fine(3,:)      
         vwnd_fine(i,:)=vwnd_fine(3,:)
       enddo

      endif

      if(n_south.eq.-1) then
       do j=1,2
         uwnd_fine(:,j)=uwnd_fine(:,3)       
         vwnd_fine(:,j)=vwnd_fine(:,3)        
       enddo 
      endif
         uwnd_buf(:,:,k)=uwnd_fine
         vwnd_buf(:,:,k)=vwnd_fine
 
       enddo ! k=1,4 
              
       else

          uwnd_buf = reshape(uwnd_buf_coarse, shape(uwnd_buf))    
          vwnd_buf = reshape(vwnd_buf_coarse, shape(vwnd_buf))  
  
       endif !  if(calc_interp) then fhx:interp_flag
                
         do i=1,3

            if (  sec_in_day >= nsec(i) .and.
     $            sec_in_day <  nsec(i+1)  ) then               

            uwnd_a = uwnd_buf( :, :, i )
            vwnd_a = vwnd_buf( :, :, i )
            uwnd_b = uwnd_buf( :, :, i+1 )
            vwnd_b = vwnd_buf( :, :, i+1 )
           
               exit
               
            endif   

         enddo
         
         if (  sec_in_day >= nsec(4) ) then            

            uwnd_a = uwnd_buf( :, :, 4 )
            vwnd_a = vwnd_buf( :, :, 4 )
            
            d_tmp2 = d_in + 24 * 3600
            n = int(dif_date(d_tmp2, d_tmp)/86400.)*4+1

            if (inc_leap(d_in%year).eq.1) then
              n = modulo(int(dif_date(d_tmp2, d_tmp)/86400.)*4, 1463)
            else
              n = modulo(int(dif_date(d_tmp2, d_tmp)/86400.)*4, 1459)
            end if

!!*      write( infile, '( //trim(windf)//"_",i4.4,2i2.2,".nc" )' )
!!*  $        d_in%year, d_in%month, d_in%day     
!         write( infile, '( a4,"_",i4.4,2i2.2,".nc" )' )
!debug     $        windf, d_in%year, d_in%month, d_in%day     
!     $        windf, d_tmp%year, d_tmp%month, d_tmp%day
         write( infile, '( a3,".",i4.4,".nc" )' )
     $        "mfl", d_tmp2%year
!!       write( infile, '( "gfsw_",i4.4,2i2.2,".nc" )' )   ! fhx:read gfsw wind
!!   $        d_in%year, d_in%month, d_in%day     
!            write( infile, '( "gfs_",i4.4,2i2.2,".nc" )' ) ! fhx:read gfs wind
!     $           d_tmp%year, d_tmp%month, d_tmp%day 

! fhx: check wind data exists,is ~exist, set wind to be 0.10/26/2010             
         inquire(file='in/'//trim(windf)//'/'//trim(infile),
     $           exist=lexist)   
!!       inquire(file='in/gfsw/'//trim(infile),exist=lexist)   

         if(lexist) then

            call read_wind_pnetcdfc
     $         ( uwnd_buf_coarse, vwnd_buf_coarse, trim(infile), n )

         else
              if ( my_task == master_task ) then
                 write(*,'(/2a)') 
     $              "missing wind data at wind_init : "
     $              , trim(infile)
               endif
             uwnd_buf_coarse = 0.
             vwnd_buf_coarse = 0.
             
         endif

! interpolation
      if(calc_interp) then ! fhx:interp_flag:add flag for interp fgrid.  

       do k=1,4  
        uwnd_coarse = uwnd_buf_coarse( :, :, k )
        vwnd_coarse = vwnd_buf_coarse( :, :, k )
        call interp_mask_2d(uwnd_coarse,1,east_e,north_e,fsm,uwnd_fine)
        call interp_mask_2d(vwnd_coarse,2,east_e,north_e,fsm,vwnd_fine)
!fhx: after interpolation, i=1,2 and j=1,2 seem to have problems
      if (n_west.eq.-1) then
       do i=1,2
         uwnd_fine(i,:)=uwnd_fine(3,:)      
         vwnd_fine(i,:)=vwnd_fine(3,:)
       enddo
      endif

      if(n_south.eq.-1) then
       do j=1,2
         uwnd_fine(:,j)=uwnd_fine(:,3)       
         vwnd_fine(:,j)=vwnd_fine(:,3)        
       enddo 
      endif
         uwnd_buf(:,:,k)=uwnd_fine
         vwnd_buf(:,:,k)=vwnd_fine
 
       enddo ! k=1,4 

       else

          uwnd_buf = reshape(uwnd_buf_coarse, shape(uwnd_buf))    
          vwnd_buf = reshape(vwnd_buf_coarse, shape(vwnd_buf))  
  
       endif !  if(calc_interp) then !fhx:interp_flag

            uwnd_b = uwnd_buf( :, :, 1 )
            vwnd_b = vwnd_buf( :, :, 1 )

         endif

      endif   !if ( calc_wind ) then


      if ( my_task == master_task ) then
         write(*,'(/a/)') "----------- wind_init."
      endif

      return

      end subroutine wind_init
!-----------------------------------------------------------------

!=================================================================
! Read & time-interpolate wind data. 
! Calclate wind stress wusurf and wvsurf.
!-----------------------------------------------------------------
      subroutine wind_main( d_in )

      use module_time
      use interp

      implicit none
      

      ! intent(in)
      type(date), intent(in) :: d_in
!      logical, intent(in) :: initial

      ! local     
      type(date) :: d_tmp
      integer n

      logical :: lexist     

      sec_in_day = d_in%hour*3600 + d_in%min*60 + d_in%sec

      d_tmp = str2date("1979-01-01 00:00:00")
      d_tmp%year = d_in%year
      n = int(dif_date(d_in, d_tmp)/86400.)

      do i=1,3
         
         if (  sec_in_day >= nsec(i) .and.
     $         sec_in_day - nsec(i) < int( dti )  ) then
            
            uwnd_a = uwnd_buf( :, :, i )
            vwnd_a = vwnd_buf( :, :, i )
            uwnd_b = uwnd_buf( :, :, i+1 )
            vwnd_b = vwnd_buf( :, :, i+1 )

            exit
            
         endif
         
      enddo


         
      if (  sec_in_day >= nsec(4) .and.
     $      sec_in_day - nsec(4) < int ( dti )  ) then 

            uwnd_a = uwnd_buf( :, :, 4 )
            vwnd_a = vwnd_buf( :, :, 4 )
         
         d_tmp = d_in + 24 * 3600

!!*      write( infile, '( //trim(windf)//"_",i4.4,2i2.2,".nc" )' )
!!*  $        d_in%year, d_in%month, d_in%day     
!         write( infile, '( a4,"_",i4.4,2i2.2,".nc" )' )
!debug     $        windf, d_in%year, d_in%month, d_in%day    
!     $        windf, d_tmp%year, d_tmp%month, d_tmp%day
         write( infile, '( a3,".",i4.4,".nc" )' )
     $        "mfl", d_tmp%year
!!       write( infile, '( "gfsw_",i4.4,2i2.2,".nc" )' )   ! fhx:read gfsw wind
!!   $        d_in%year, d_in%month, d_in%day     
!         write( infile, '( "gfs_",i4.4,2i2.2,".nc" )' )  ! fhx:read gfs wind
!     $        d_tmp%year, d_tmp%month, d_tmp%day

! fhx: check wind data exists,is ~exist, set wind to be 0.10/26/2010         
         inquire(file='in/'//trim(windf)//'/'//trim(infile),
     $           exist=lexist)   
!!       inquire(file='in/gfsw/'//trim(infile),exist=lexist)   

         if(lexist) then

         call read_wind_pnetcdfc
     $        ( uwnd_buf_coarse, vwnd_buf_coarse, trim(infile), n*4+1 )
         
         else
            if ( my_task == master_task ) then
                 write(*,'(/2a)') 
     $              "missing wind data at wind_main : "
     $              , trim(infile)
             endif
             uwnd_buf_coarse = 0.
             vwnd_buf_coarse = 0.
             
         endif

! interpolation
      if(calc_interp) then ! fhx:interp_flag:add flag for interp fgrid.  

       do k=1,4  
        uwnd_coarse = uwnd_buf_coarse( :, :, k )
        vwnd_coarse = vwnd_buf_coarse( :, :, k )
        call interp_mask_2d(uwnd_coarse,1,east_e,north_e,fsm,uwnd_fine)
        call interp_mask_2d(vwnd_coarse,2,east_e,north_e,fsm,vwnd_fine)
!fhx: after interpolation, i=1,2 and j=1,2 seem to have problems
      if (n_west.eq.-1) then
       do i=1,2
         uwnd_fine(i,:)=uwnd_fine(3,:)      
         vwnd_fine(i,:)=vwnd_fine(3,:)
       enddo
      endif

      if(n_south.eq.-1) then
       do j=1,2
         uwnd_fine(:,j)=uwnd_fine(:,3)       
         vwnd_fine(:,j)=vwnd_fine(:,3)        
       enddo 
      endif
         uwnd_buf(:,:,k)=uwnd_fine
         vwnd_buf(:,:,k)=vwnd_fine
 
       enddo ! k=1,4 

       else

          uwnd_buf = reshape(uwnd_buf_coarse, shape(uwnd_buf))    
          vwnd_buf = reshape(vwnd_buf_coarse, shape(vwnd_buf))          
  
       endif !  if(calc_interp) then fhx:interp_flag  
         
           uwnd_b = uwnd_buf( :, :, 1 )
           vwnd_b = vwnd_buf( :, :, 1 )

      endif
 

      aa = real( mod( sec_in_day, 6 * 3600 ) ) / ( 6. * 3600. )

!lyo:pac10:more efficient:Comment out the followings
! prepare wind data for SURF output !fhx:20110131:
!       uwind_surf = ( 1.0 - aa ) * uwnd_a + aa * uwnd_b
!       vwind_surf = ( 1.0 - aa ) * vwnd_a + aa * vwnd_b
!       uwsrf = uwind_surf
!       vwsrf = vwind_surf

      do j=1,jm
         do i=1,im

            uwnd = ( 1.0 - aa ) * uwnd_a( i, j ) + aa * uwnd_b( i, j )
            vwnd = ( 1.0 - aa ) * vwnd_a( i, j ) + aa * vwnd_b( i, j )

!lyo:pac10:exp016:reduce wind to zero west of 129
!           rdisp=0.5*(1.+tanh((east_e(i,j)-129.0)*0.5))
!           uwnd=uwnd*rdisp; vwnd=vwnd*rdisp
!lyo:pac10:more efficient:
            uwsrf(i,j)=uwnd; vwsrf(i,j)=vwnd;

            uvabs = sqrt( uwnd**2 + vwnd**2 )


!     test ayumi 2010/7/20
!
!            hmin = 15. ; wndmax = 15.
!            if ( h(i,j) <= hmin .and. uvabs > wndmax ) then
!               uwnd = uwnd * wndmax / uvabs
!               vwnd = vwnd * wndmax / uvabs
!               uvabs = wndmax
!            endif


            if (     uvabs <= 11.0 ) then
               cda = 0.0012
            elseif ( uvabs <= 19.0 ) then
               cda = 0.00049 + 0.000065*uvabs
            elseif ( uvabs <= 100.0 ) then
               cda = 0.001364 + 0.0000234*uvabs - 2.31579e-7*uvabs**2
            else
               cda = 0.00138821   !modify cda = 0.000138821 ---> 0.00138821
            endif

            wusurf( i, j ) = - rhoa / rhow * cda * uvabs * uwnd !- uwnd / rhow
            wvsurf( i, j ) = - rhoa / rhow * cda * uvabs * vwnd !- vwnd / rhow

         enddo
      enddo

!     wind wave-induced enhanced bottom drag !lyo:20110315:botwavedrag:4lines
!       if ( calc_botwavedrag ) call botwavedrag(...) !future use calc_..
        call botwavedrag (im,jm,fsm,wusurf,wvsurf,
     $  0.0314159,          !kp=2.*pi/200.=0.0314159
     $  wubot,wvbot,h,zz(kbm1),z0b,cbcmin,cbcmax,cbc)
      

      
      return

      end subroutine wind_main
!-----------------------------------------------------------------

      end module wind
