      module assim

      use config     , only: calc_interp
      use glob_const , only: rk
      use glob_domain, only: im, jm
      use glob_grid  , only: east_e, north_e
      use module_time

      implicit none

      private

      public :: assim_init, assim_main, assim_store_ssha

      integer :: iassim

      real(kind=rk), allocatable, dimension(:,:) ::
     $     ssha_a, ssha_b, ssha
     $   , frs, frs_coarse
     $   , ssha_a_coarse, ssha_b_coarse!, ssha_coarse

      real(kind=rk), allocatable, dimension(:,:,:) ::
     $     tav, fac, cof
     $   , tav_coarse, fac_coarse, cof_coarse
     

      real(kind=rk), parameter :: errval = -999.

      real(kind=rk) :: nassim

      type(date) :: d1, d2, dstart, d_buf

      character(len=19) ::
     $     dstart_str = "1992-10-22_00:00:00"

      character(len=30) :: infile

      logical, allocatable, dimension(:,:) ::
     $     mask_msla

      integer ::  i, j
      integer :: count_days            !fhx:debug


      namelist/assim_nml/ nassim

      contains

!==================================================================
! Initialize 
!------------------------------------------------------------------
      subroutine assim_init( d_in )

      use bry        , only: frz, USE_SPONGE
      use config     , only: calc_assim
      use glob_domain, only: im, im_coarse, is_master
     &                     , jm, jm_coarse, kb
     &                     , n_south, n_west
      use model_run  , only: dti
      use interp

      implicit none

      ! intent(in)
      type(date), intent(in) :: d_in
    
      ! local
      integer :: status, access
!      logical :: lexist

      allocate(
     &  ssha_a(im,jm), ssha_b(im,jm), ssha(im,jm), frs(im,jm)
     &, ssha_a_coarse(im_coarse,jm_coarse)
     &, ssha_b_coarse(im_coarse,jm_coarse)
     &, frs_coarse(im_coarse,jm_coarse)
     &, tav(im,jm,kb), fac(im,jm,kb), cof(im,jm,kb)
     &, tav_coarse(im_coarse,jm_coarse,kb)
     &, fac_coarse(im_coarse,jm_coarse,kb)
     &, cof_coarse(im_coarse,jm_coarse,kb)
     &, mask_msla(im,jm)
     & )

      if ( calc_assim ) then

         dstart = str2date( dstart_str )
                     
!     read namelist
         open(73,file='switch.nml',status='old')
         read(73, nml=assim_nml)
         close(73)

         iassim = max( int( nassim * 24 * 3600 / dti ), 1 )
      
         if ( is_master ) 
     $        print '(/a,2(i4,a))', "assimilation interval : ",
     $        int(nassim), "days, every  ", iassim, " time steps."

!     read information for data assimilation
      
         call read_assiminfo_pnetcdfc
     $    ( frs_coarse, tav_coarse, fac_coarse, cof_coarse )

! fhx: interpolation to fine grid. 12/13/2010
      if(calc_interp) then ! fhx:interp_flag:add flag for interp fgrid.

        call interp_mask_2d(frs_coarse,0,east_e,north_e,frs)
        call interp_mask_3d(tav_coarse,0,east_e,north_e,tav)
        call interp_mask_3d(fac_coarse,0,east_e,north_e,fac)
        call interp_mask_3d(cof_coarse,0,east_e,north_e,cof)

      if (n_west.eq.-1) then

!        print*,frs(1,1),frs(2,2),frs(3,3),frs(4,4)
!        print*,fac(1,1,1),fac(2,2,2),fac(3,3,3),fac(4,4,4)
!        print*,''
       do i=1,3
         frs(i,:)=frs(4,:)
         tav(i,:,:)=tav(4,:,:)         
         fac(i,:,:)=fac(4,:,:)
         cof(i,:,:)=cof(4,:,:)     
       enddo
      endif

      if(n_south.eq.-1) then
       do j=1,3
         frs(:,j)=frs(:,4)
         tav(:,j,:)=tav(:,4,:)         
         fac(:,j,:)=fac(:,4,:)
         cof(:,j,:)=cof(:,4,:)     
       enddo 
      endif

      else

      frs = reshape(frs_coarse,shape(frs))
      tav = reshape(tav_coarse,shape(tav))
      fac = reshape(fac_coarse,shape(fac))
      cof = reshape(cof_coarse,shape(cof))

      endif  ! if(calc_interp) then !fhx:interp_flag

! fhx: test assim with missing msla*.nc. 10/26/2010
 
        
         if ( dstart > d_in ) then

            if ( is_master ) then
               print '(/2a)'
     &             , "Before assimilation period. d_in : "
     &             ,  date2str(d_in)
               print '(2a)', "no assimilation until "
     &             , dstart_str
            end if
            
            d1 = dstart
            d2 = dstart

         else

 
      !     decide the backward date.

      ! set initial d_buf 10-20days before d_in.
 
            d_buf = dstart + 
     $    max(( int(( d_in - dstart ) / ( 10 * 86400 ) ) - 1 ) 
     $    * ( 10 * 86400 ), 0 ) 
             
      
            do while ( d_in >= d_buf )
               write( infile, 
     $              '( "./in/assim/msla_",i4.4,2i2.2,".nc" )' ) 
     $              d_buf%year, d_buf%month, d_buf%day
               status = access( trim(infile), 'r' )         
               if ( status == 0 ) then
                  d1 = d_buf
               else
                 if ( is_master ) then
                   print '(/2a)'
     $                 , "missing satellite data at assim_init : "
     $                 , trim(infile)
                 endif

!               go to 200   !fhx:debug
               endif       ! fhx: return if missing msla_*.nc. 10/26/2010  
               d_buf = d_buf + 24 * 3600
            end do

            
      !     decide the forward date.

            d_buf = d1 + 24 * 3600
            status = -1
            count_days = 0      !fhx:debug
            do while ( status /=0 )
               write( infile, 
     $              '( "./in/assim/msla_",i4.4,2i2.2,".nc" )' ) 
     $              d_buf%year, d_buf%month, d_buf%day
               status = access( trim(infile), 'r' )         
               if ( status == 0 ) then 
                   d2 = d_buf
               else
!                 if ( my_task == master_task ) then
!                   write(*,'(/2a)') 
!     $              "missing satellite data at assim_init : "
!     $              , trim(infile)
!                 endif
               count_days = count_days+1     !fhx:debug
               if(count_days>12)  go to 200  !fhx:debug
               endif       ! fhx: return if missing msla_*.nc. 10/26/2010 
               d_buf = d_buf + 24 * 3600
            end do

         endif
        
         
!     read data
     
         write( infile, '( "msla_",i4.4,2i2.2,".nc" )' ) 
     $        d1%year, d1%month, d1%day

         call read_msla_pnetcdfc( ssha_a_coarse, trim(infile) )    

         write( infile, '( "msla_",i4.4,2i2.2,".nc" )' ) 
     $        d2%year, d2%month, d2%day

         call read_msla_pnetcdfc( ssha_b_coarse, trim(infile) )

! interpolate after read in ssha. fhx:12/13/2010  
      if(calc_interp) then ! fhx:interp_flag:add flag for interp fgrid. 

        call interp_mask_2d(ssha_a_coarse,0,east_e,north_e,ssha_a)
        call interp_mask_2d(ssha_b_coarse,0,east_e,north_e,ssha_b)
      if (n_west.eq.-1) then
       do i=1,3
         ssha_a(i,:)=ssha_a(4,:)
         ssha_b(i,:)=ssha_b(4,:)
       enddo
      endif

      if(n_south.eq.-1) then
       do j=1,3
         ssha_a(:,j)=ssha_a(:,4)
         ssha_b(:,j)=ssha_b(:,4)
       enddo 
      endif

      else

      ssha_a = reshape(ssha_a_coarse,shape(ssha_a))
      ssha_b = reshape(ssha_b_coarse,shape(ssha_b))

      endif  ! if(calc_interp) then !fhx:interp_flag

      
         if ( is_master ) then
            write(*,'(/4a)')  "assim backward & forward : "
     $           , date2str(d1)," ", date2str(d2)
            write(*,'(/a/)')  "---------- assim_init. "
         endif
         
      else
         

         if ( is_master ) then
            write(*,'(/a)')  "no assimilation. "
            write(*,'(/a/)')  "---------- assim_init. "
         endif

      endif



!      if ( n_east == -1 ) 
!     $     print*, frs(im/2,jm/2), tav(im/2,jm/2,1), 
!     $     fac(im/2,jm/2,1), cof(im/2,jm/2,1)

!200   return
!yoyo:200 return
200   continue

      if ( USE_SPONGE ) then
        frs = frz !lyo:pac10:exp032:quick fix for frs
                          !instead of reading from assiminfo.nc
                          !which may have frs=1 everywhere
      else
        frs = 0.
      end if

      return

      end subroutine assim_init
!==================================================================


!==================================================================
! Main
!------------------------------------------------------------------
      subroutine assim_main( d_in )

!      use glob_config, only: 
      use glob_domain, only: i_global, im, j_global, jm, kb, kbm1
     &                     , my_task, n_proc, n_south, n_west
      use glob_grid  , only: h, zz
      use glob_ocean , only: t, tb, tclim
      use model_run  , only: iint
      use module_time
      use interp

      implicit none

      ! intent(in)
      type(date), intent(in) :: d_in

      ! local
      real(kind=rk) :: aa
      integer :: i, j ,k 
      integer :: sec_d1din, sec_d1d2
      integer :: status, access      

      real(kind=rk), parameter :: dtfac = 1.0 / 40.0    
      real(kind=rk) :: cfrs, corml, corhh, corbd
      real(kind=rk) :: c2, pt, cfg2, zlev, tobs
      
!      logical :: lexist

!     no assimilation before DSTART or mod( iint, iassim ) /= 0

      if ( dstart > d_in .or. mod(iint-1, iassim) /= 0 ) return
!      if ( dstart > d_in .or. mod(iint, iassim) /= 0 ) return


!     do assimilation

      if ( d_in >= d2 ) then   

         
      !  renew the backward and forward data

      ! set initial d_buf 10-20days before d_in.
 
         d_buf = dstart + 
     $    max(( int(( d_in - dstart ) / ( 10 * 86400 ) ) - 1 ) 
     $    * ( 10 * 86400 ), 0) 
         
            do while ( d_in >= d_buf )
               write( infile, 
     $              '( "./in/assim/msla_",i4.4,2i2.2,".nc" )' ) 
     $              d_buf%year, d_buf%month, d_buf%day
               status = access( trim(infile), 'r' )   
               if ( status == 0 ) then 
                 d1 = d_buf
               else
!                 if ( my_task == master_task ) then
!                   write(*,'(/2a)') 
!     $              "missing satellite data at assim_main : "
!     $              , trim(infile)
!                 endif
!                d1 = d_buf 
               endif   ! fhx: return if missing msla_*.nc. 10/26/2010 
                d_buf = d_buf + 24 * 3600
            end do
            
      !     decide the forward date.

            d_buf = d1 + 24 * 3600
            status = -1
            count_days = 0          !fhx:debug
            do while ( status /=0 )
               write( infile, 
     $              '( "./in/assim/msla_",i4.4,2i2.2,".nc" )' ) 
     $              d_buf%year, d_buf%month, d_buf%day
               status = access( trim(infile), 'r' )     
               
               if ( status == 0 ) then
                 d2 = d_buf
               else

               count_days = count_days+1     !fhx:debug
!                if ( my_task == master_task ) then
!                  print*,'count_days:',count_days
!                  write(*,'(/2a)') 
!     $              "missing satellite data at assim_main : "
!     $              , trim(infile)
!                endif

!                d2 = d_buf 

                if(count_days>12) go to 300  !fhx:debug

               endif   ! fhx: return if missing msla_*.nc. 10/26/2010 

               d_buf = d_buf + 24 * 3600
            end do

!     read data
   
         write( infile, '( "msla_",i4.4,2i2.2,".nc" )' ) 
     $        d1%year, d1%month, d1%day

         call read_msla_pnetcdfc( ssha_a_coarse, trim(infile) )

         write( infile, '( "msla_",i4.4,2i2.2,".nc" )' ) 
     $        d2%year, d2%month, d2%day

         call read_msla_pnetcdfc( ssha_b_coarse, trim(infile) )

! interpolate after read in ssha. fhx:12/13/2010
       if(calc_interp) then !fhx:interp_flag:add flag for interp fgrid. 

        call interp_mask_2d(ssha_a_coarse,0,east_e,north_e,ssha_a)
        call interp_mask_2d(ssha_b_coarse,0,east_e,north_e,ssha_b)
       if (n_west.eq.-1) then
       do i=1,3
         ssha_a(i,:)=ssha_a(4,:)
         ssha_b(i,:)=ssha_b(4,:)
       enddo
       endif

       if(n_south.eq.-1) then
       do j=1,3
         ssha_a(:,j)=ssha_a(:,4)
         ssha_b(:,j)=ssha_b(:,4)
       enddo 
       endif
      else

      ssha_a = reshape(ssha_a_coarse,shape(ssha_a))
      ssha_b = reshape(ssha_b_coarse,shape(ssha_b))

      endif   ! if(calc_interp) then !fhx:interp_flag
        
      endif   ! if ( d_in >= d2 ) then   

      
         
!     seconds from the previous data date.

      sec_d1din = d_in - d1 


!     time interval between the data.    
    
      sec_d1d2 = d2 -d1

      aa = real( sec_d1din ) / real( sec_d1d2 )


!     error mask
     
      where ( ssha_a == errval .or. ssha_b == errval ) 
         mask_msla = .false.
         ssha = 0.0
      elsewhere
         mask_msla = .true. 
         ssha = ( 1.0 - aa ) * ssha_a + aa * ssha_b
      endwhere
      

!     unit chage [cm] --> [m]
      
      ssha = ssha * 0.01

     
!      if ( my_task == 20 ) 
!     $     write(*,'(a,i3,5f11.4)') date2str(d_in),my_task, aa,(1.-aa)
!     $     ,ssha_a(im_local/2,jm_local/2),ssha_b(im_local/2,jm_local/2)
!     $     ,ssha(im_local/2,jm_local/2)*100.
      

!     use monthly climatology instead of the model mean
      tav(:,:,1:kbm1) = tclim(:,:,1:kbm1)

!     before
      if ( my_task == n_proc/2-1 )  
     $     print'(/2(a,f12.5),a,l2)',
     $     "before : tb = ", tb(im/2,jm/2,1), 
     $     ", t = ", t(im/2,jm/2,1),
     $     ", mask_msla = ",mask_msla(im/2,jm/2)



!     calculation of coefficients & analysya t, tb

      do j = 1,jm
         do i = 1,im
            if ( mask_msla(i,j) ) then
               do k = 1, kbm1

                  ! determines coefficients
                  corml = 1.0; corhh = 1.0; corbd = 1.0

                  ! taper corml to zero 1000-1500m, 
                  ! then no assimilation below 1500m.
                  zlev = -1.0 * zz(k) * h(i,j)
                  if ( zlev > 1000.0 ) 
     $                 corml = 1.0 - ( zlev - 1000.0 ) / 500.0
                  if ( zlev >= 1500.0 ) corml = 0.0


                  ! no assimilation on shelves, taper to one 200-500m
                  if ( h(i,j) <= 200.0 ) corhh = 0.0
                  if ( h(i,j) > 200.0 .and. h(i,j) < 500.0 )
     $                 corhh = 1.0 - ( 500.0 - h(i,j) ) / 300.0

                  ! assume some interpolation erroe exists so max = 0.95
                  corhh = min( corhh, 0.95_rk )

                  ! cfrs suppresses assimilation near boundaries
                  ! frs = 1 near eastern boundary
                  cfrs = 1.0 - frs(i,j)


                  ! calculate the observed value
                  tobs = tav(i,j,k) + fac(i,j,k) * ssha(i,j)

                  ! calculate the analysis value 
                  c2 = cof(i,j,k)**2
                  cfg2 = corbd * corml * corhh * 2.0 * dtfac * c2
                  pt = cfg2 / ( 1.0 + cfg2 - c2 )
                  t(i,j,k) = t(i,j,k) 
     $                 + cfrs * pt * ( tobs - t(i,j,k) ) 
                  tb(i,j,k) = tb(i,j,k) 
     $                 + cfrs * pt * ( tobs - tb(i,j,k) )
               enddo
            endif
         enddo
      enddo



      if ( my_task == n_proc/2-1 ) then 
         write(*,'(/2(a,f12.5),a,2i4)') 
     $        "assim-ssha done: tb@kb/4 = ", tb(im/2,jm/2,kb/4), 
     $        ", t@kb/4 = ", t(im/2,jm/2,kb/4),
     $        ", (i,j) = ", i_global(im/2),j_global(jm/2)
         
         
         write(*,'(/a/)')  "---------- assim_main. "

      endif

300   return

      end subroutine assim_main
!==================================================================

!==================================================================
! store ssha at t(k=kb) which is not used.
! note that ssha is at the last assimilation time.
!------------------------------------------------------------------
      subroutine assim_store_ssha( var, varname, d_in )

      use glob_domain, only: im, is_master, jm
      use module_time

      real(kind=rk), dimension( im, jm )
     $     , intent(inout) :: var
      character(len=*) , intent(in) :: varname
      type(date), intent(in) :: d_in

      if ( dstart > d_in  ) then

         return

      else

         var = ssha

         if ( is_master ) then
            write(*,'(/2a)')  
     $           " ssha used in last assmilation is stored to "
     $           , varname
            write(*,'(/a/)')  "---------- assim_store_ssha. "
      
         endif
      
      return

      endif

      end subroutine assim_store_ssha              
!==================================================================

!==================================================================
! calculate weighting coefficient frs
!------------------------------------------------------------------
!      subroutine calc_frs( frs_out )
!      
!      real(kind=rk), dimension( im_local, jm_local ), intent(out) :: frs_out
!
!      frs_out = 0.0
!
! ....      
!      end function calc_frs
!
!==================================================================

      end module assim
