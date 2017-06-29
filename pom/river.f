      module river
      
      implicit none

      private
      
      public :: river_init, river_main, totq

      include 'pom.h'

      integer, parameter :: nr_gom = 33 , nr_mab = 23
      
      real(kind=rk), dimension( nr_gom ) :: riv_gom_b, riv_gom_f
      real(kind=rk), dimension( nr_mab ) :: riv_mab_b, riv_mab_f

      real(kind=rk) :: totq
      real(kind=rk) :: ffac

      character(len=19) :: 
     $     dstart_river_str   = "1992-01-01_00:00:00",
!     $     dend_river_gom_str = "2000-09-30_00:00:00" ,  !fhx:2011/01/18
     $     dend_river_gom_str = "2010-10-31_00:00:00" ,   !fhx:2011/01/18
     $     dend_river_mab_str = "2008-12-31_00:00:00" 

      character(len=50) ::
     $     rivers_gom_file = './in/rivers_gom.bin',
     $     rivers_mab_file = './in/rivers_mab.bin'

      integer :: i_gom( nr_gom ) =
     $(/ 117, 107, 100,  99,  90,  82,  81,  80,  75,  74, 
     $    71,  70,  72,  65,  63,  64,  62,  61,  34,  32, 
     $    29,  28,  30,  22,  22,  18,  13,  10,  10,   6,
     $     6,  62,  46 /)
      integer :: i_gomf(nr_gom)

      integer :: j_gom( nr_gom ) =
     $(/ 131, 145, 144, 145, 167, 171, 171, 171, 177, 177, 
     $   182, 182, 179, 186, 183, 184, 182, 182, 198, 201,
     $   201, 201, 202, 194, 195, 183, 175, 166, 165, 151,
     $   109, 151, 179 /)
      integer :: j_gomf(nr_gom)

      integer :: i_mab( nr_mab ) = 
     $(/ 166, 174, 176, 191, 211, 216, 219, 224, 227, 232, 
     $   233, 246, 238, 249, 268, 287, 291, 301, 317, 323,
     $   339, 343, 357 /)
      integer :: i_mabf( nr_mab )
 
      integer :: j_mab( nr_mab ) = 
     $(/ 166, 167, 171, 174, 175, 175, 223, 221, 237, 245,
     $   239, 246, 231, 248, 241, 232, 230, 237, 244, 243,
     $   239, 236, 251 /)
      integer :: j_mabf( nr_mab )
!fhx:river_gen:beg
      integer, parameter :: numc = 7  !# of coefficients 
      real(kind=rk) :: T1, T2, T3               !periods of harmonic analysis
      real(kind=rk) :: river_time              !Julian days to calculte river discharge
!fhx:river_gen:end
! Belows are river points in the previous serial code. 
! In the present points above, 1pts added in i&j direction
! according to 1pt increase of grid number in i&j direction.
! Also, j_gom is added by 14. 
! This was done inside the previous serial code.
! 
!      integer :: i_gom( nr_gom ) =
!     $(/ 116, 106,  99,  98,  89,  81,  80,  79,  74,  73, 
!     $    70,  69,  71,  64,  62,  63,  61,  60,  45,  33, 
!     $    31,  28,  27,  29,  21,  21,  17,  12,   9,   9,
!     $     5,   5,  61,  45 /)

!      integer :: j_gom( nr_gom ) =
!     $(/ 116, 130, 129, 130, 152, 156, 156, 156, 162, 162, 
!     $   167, 167, 164, 171, 168, 169, 167, 167, 164, 183, 
!     $   186, 186, 186, 187, 179, 180, 168, 160, 151, 150, 
!     $   136, 94, 136, 164 /)

!      integer :: i_mab( nr_mab ) = 
!     $(/ 165, 173, 175, 190, 210, 215, 218, 223, 226, 231, 
!     $   232, 245, 237, 248, 267, 286, 290, 300, 316, 322,
!     $   338, 342, 356 /)

!      integer :: j_mab( nr_mab ) = 
!     $(/ 165, 166, 170, 173, 174, 174, 222, 220, 236, 244,
!     $   238, 245, 230, 247, 240, 231, 229, 236, 243, 242,
!     $   238, 235, 250 /)


      contains
      
!=================================================================
! Initialize variables for river
!-----------------------------------------------------------------
      subroutine river_init( d_in )

      use module_time

      implicit none

!     intent (in)
      type( date ), intent( in ) :: d_in

!     initialize

      totq = 0.0

!     river fluxes in an initial state
      if ( calc_river ) call river_main( d_in, .true. )

      if ( my_task == master_task ) 
     $     write(*,'(/a/)') "---------- river_init."


      return

      end subroutine river_init
!-----------------------------------------------------------------

!=================================================================
! Main routine for river
! Calculate total volume flux by river: totq
!-----------------------------------------------------------------
      subroutine river_main( d_in, initial )
     
      use module_time

      implicit none

!     intent (in)

      type( date ), intent( in ) :: d_in
      logical, intent( in ) :: initial

      real(kind=rk), dimension(im,jm) :: vol_gom, vol_mab
      real(kind=rk) :: totq_gom, totq_mab

!     d_in     : present date
!     initial  : .true. for initializarion 
!     vol_gom  : grided volume flux for GOM [m3/s]
!     vol_mab  : grided volume flux for MAB [m3/s]
!     totq_gom : total volume flux for GOM [m3/s]
!     totq_mab : total volume flux for MAB [m3/s]

      vol_gom = 0.0; vol_mab = 0.0

!     River flux [m3/s] Gulf of Mexico.

      call river_gom( d_in, initial, vol_gom, totq_gom  )


!     River flux [m3/s] Mid Atlantic Bight.

      call river_mab( d_in, initial, vol_mab, totq_mab )


!     Inflow == downward ( negative ) velocity.
!     unit change [m3/s] -> [m/s]

      vfluxf = -1.0 * ( vol_gom + vol_mab ) / art
      

!     total discharge [m3/s]

      totq = totq_gom + totq_mab


!      print*, totq_gom, totq_mab, date2str(d_in), my_task


      return

      end subroutine river_main
!-------------------------------------------------------------

!=============================================================
! Read & time-interpolate river data for GOM
!-------------------------------------------------------------
      subroutine river_gom( d_in, initial, vol_out, totq_out )

      use module_time

      implicit none

      include 'pom.h'

!     intent(in)
      logical, intent( in ) :: initial
      type( date ), intent( in ) :: d_in

!     intent(out)
      real(kind=rk), dimension(im,jm), intent( out ) :: vol_out
      real(kind=rk) :: totq_out

!     local
      integer :: i, sec_in_day, elapsed_days, n_b, n_f
      real(kind=rk) :: wtime, aa, riv( nr_gom ) 
!      logical :: here, judge_inout
      logical :: here
      type( date ) :: dstart_river, dend_river
            

      dstart_river = str2date( dstart_river_str )
      dend_river   = str2date( dend_river_gom_str )


      if ( d_in >= dstart_river  
     $     .and. dend_river > d_in  ) then

         sec_in_day = d_in%hour * 3600 + d_in%min * 60 + d_in%sec

      if ( initial 
     $     .or. sec_in_day < int(dti) ) then 

         elapsed_days = ( d_in - dstart_river ) / 86400 
         
         n_b = elapsed_days + 1
         n_f = n_b + 1
         
         open(10, file=trim(rivers_gom_file),form='unformatted',
     $        access='direct', recl = (nr_gom+1)*4, status ='old')
!         if ( my_task == 0 )  print'(3a)', 
!     $        'reading file ',trim(rivers_gom_file),' .'
         read(10, rec=n_b ) wtime, riv_gom_b
         read(10, rec=n_f ) wtime, riv_gom_f
         close(10)
         
      endif
      
!     time interpolation
      aa = real( sec_in_day ) / 86400.
      riv = ( 1.0 - aa ) * riv_gom_b + aa * riv_gom_f 

!fhx:river_gen:beg:2011/10/19
      else if (d_in >= dend_river) then
            if(1.eq.1) then
              sec_in_day = d_in%hour * 3600 + d_in%min * 60 + d_in%sec
              elapsed_days = ( d_in - dstart_river ) / 86400 
              river_time = real(elapsed_days)+real(sec_in_day)/86400.-4.!fhx:adjust time
              call river_gen_gom(river_time,riv)
            else
              open(10, file=trim(rivers_gom_file),form='unformatted',
!fhx: put a steady river inflow. 2010/12/29
     $        access='direct', recl = (nr_gom+1)*4, status ='old')
!             read(10, rec=3196) wtime, riv   ! last record of old available river      
              read(10, rec=6879) wtime, riv   ! last record of new available river.fhx, 2011/1/18.
              close(10)
            endif !if(1.eq.1) then 
      endif
!fhx:river_gen:end
!     total discharge [m3/s]
      totq_out = sum( riv )

!             if ( my_task == master_task ) then
!                  print*,'steady rivers in GOM :'
!                  print*, totq_out,riv(1),riv(10)
!             endif

      if(calc_interp) then !fhx:interp_flag
     
!fhx:locate i_gom&j_gom to i_gomf&j_gomf in fine-grid locations.2010/12/07
! factor to evenly distribute the coarse-grid river dischange
      ffac=1./(x_division*y_division)  
    
      do i = 1, nr_gom
!         i_gomf(i) = (i_gom(i)-2+1)*x_division
!         j_gomf(i) = (j_gom(i)-2+1)*y_division
!
!fhx:fgrid:+2 to adjust river after adding dummy points at south and west boundaries
         i_gomf(i) = (i_gom(i)-2+1)*x_division + 2 
         j_gomf(i) = (j_gom(i)-2+1)*y_division + 2 !fhx:fgrid
!(i,j)    
         here = judge_inout( i_gomf(i), j_gomf(i), 
     $                       i_global(1), i_global(im),
     $                       j_global(1), j_global(jm) )
         
         if ( here ) then
            vol_out( i_gomf(i)-i_global(1) + 1, 
     $           j_gomf(i)-j_global(1) + 1 ) = riv(i)*ffac 
            
!            print*, i, i_gomf(i), j_gomf(i), riv(i)       

         endif
!(i-1,j-1)
         here = judge_inout( i_gomf(i)-1, j_gomf(i)-1, 
     $                       i_global(1), i_global(im),
     $                       j_global(1), j_global(jm) )
         
         if ( here ) then
            vol_out( i_gomf(i)-i_global(1), 
     $           j_gomf(i)-j_global(1)) = riv(i)*ffac 
            
!            print*, i, i_gomf(i)-1, j_gomf(i)-1, riv(i)       

         endif
!(i-1,j)
         here = judge_inout( i_gomf(i)-1, j_gomf(i), 
     $                       i_global(1), i_global(im),
     $                       j_global(1), j_global(jm) )
         
         if ( here ) then
            vol_out( i_gomf(i)-i_global(1), 
     $           j_gomf(i)-j_global(1) + 1 ) = riv(i)*ffac 
            
!            print*, i, i_gomf(i)-1, j_gomf(i), riv(i)       

         endif
!(i,j-1)
         here = judge_inout( i_gomf(i), j_gomf(i)-1, 
     $                       i_global(1), i_global(im),
     $                       j_global(1), j_global(jm) )
         
         if ( here ) then
            vol_out( i_gomf(i)-i_global(1) + 1, 
     $           j_gomf(i)-j_global(1) ) = riv(i)*ffac 
            
!            print*, i, i_gomf(i), j_gomf(i), riv(i)       

         endif

      enddo

      else

            do i = 1, nr_gom
    
         here = judge_inout( i_gom(i), j_gom(i), 
     $                       i_global(1), i_global(im),
     $                       j_global(1), j_global(jm) )
         
         if ( here ) then
            vol_out( i_gom(i)-i_global(1) + 1, 
     $           j_gom(i)-j_global(1) + 1 ) = riv(i) 
            
!            print*, i, i_gom(i), j_gom(i), riv(i)       

         endif

      enddo

      call exchange2d_mpi(vol_out,im,jm)  !fhx:fgrid     

      endif ! if(calc_interp) then !fhx:interp_flag
      

!      endif
      
      return

      end subroutine river_gom
!-------------------------------------------------------------

!=============================================================
! Read & time-interpolate river data for MAB
!-------------------------------------------------------------
      subroutine river_mab( d_in, initial, vol_out, totq_out )

      use module_time

      implicit none

      include 'pom.h'

!     intent(in)
      logical, intent( in ) :: initial
      type( date ), intent( in ) :: d_in

!     intent(out)
      real(kind=rk), dimension(im,jm), intent( out ) :: vol_out
      real(kind=rk) :: totq_out

!     local
      integer :: i, sec_in_day, elapsed_days, n_b, n_f
      real(kind=rk) :: wtime, aa, riv( nr_mab )
!      logical :: here, judge_inout
      logical :: here
      type( date ) :: dstart_river, dend_river
     
      dstart_river = str2date( dstart_river_str )
      dend_river   = str2date( dend_river_mab_str )

      if ( d_in >= dstart_river  
     $     .and. dend_river > d_in  ) then

         sec_in_day = d_in%hour * 3600 + d_in%min * 60 + d_in%sec

      if ( initial
     $     .or. sec_in_day < int(dti) ) then 

         elapsed_days = ( d_in - dstart_river ) / 86400 
         
         n_b = elapsed_days + 1
         n_f = n_b + 1
        
         open(10, file=trim(rivers_mab_file),form='unformatted',
     $        access='direct', recl = (nr_mab+1)*4, status ='old')
!         if ( my_task == 0 )  print'(3a)', 
!     $        'reading file ',trim(rivers_mab_file),' .'
         read(10, rec=n_b ) wtime, riv_mab_b
         read(10, rec=n_f ) wtime, riv_mab_f
         close(10)
         
      endif
      
!     time interpolation
      aa = real( sec_in_day ) / 86400.
      riv = ( 1.0 - aa ) * riv_mab_b + aa * riv_mab_f 
      
!fhx:river_gen:beg:2011/10/19
      else if (d_in >= dend_river) then
            if(1.eq.1) then
              sec_in_day = d_in%hour * 3600 + d_in%min * 60 + d_in%sec
              elapsed_days = ( d_in - dstart_river ) / 86400 
              river_time = real(elapsed_days)+real(sec_in_day)/86400.-4.!fhx:adjust time
              call river_gen_mab(river_time,riv)
            endif
      endif
!fhx:river_gen:end
!     total discharge [m3/s]
      totq_out = sum( riv )

      if(calc_interp) then !fhx:interp_flag

!fhx:locate i_mab&j_mab to i_mabf&j_mabf in fine-grid locations.2010/12/07
! factor to evenly distribute the coarse-grid river dischange
      ffac=1./(x_division*y_division)  

      do i = 1, nr_mab
!         i_mabf(i) = (i_mab(i)-2+1)*x_division
!         j_mabf(i) = (j_mab(i)-2+1)*y_division
!
!fhx:fgrid:+2 to adjust river after adding dummy points at south and west boundaries
         i_mabf(i) = (i_mab(i)-2+1)*x_division + 2 
         j_mabf(i) = (j_mab(i)-2+1)*y_division + 2 !fhx:fgrid
!(i,j)
         here = judge_inout( i_mabf(i), j_mabf(i), 
     $                       i_global(1), i_global(im),
     $                       j_global(1), j_global(jm) )
         
         if ( here )  then

            vol_out( i_mabf(i)-i_global(1) + 1, 
     $           j_mabf(i)-j_global(1) + 1 ) = riv(i) 
!            print*, i, i_mabf(i), j_mabf(i), riv(i)
            
         endif
!(i-1,j-1)
         here = judge_inout( i_mabf(i)-1, j_mabf(i)-1, 
     $                       i_global(1), i_global(im),
     $                       j_global(1), j_global(jm) )
         
         if ( here )  then

            vol_out( i_mabf(i)-i_global(1), 
     $           j_mabf(i)-j_global(1) ) = riv(i) 
!            print*, i, i_mabf(i)-1, j_mabf(i)-1, riv(i)
            
         endif
!(i-1,j)
         here = judge_inout( i_mabf(i)-1, j_mabf(i), 
     $                       i_global(1), i_global(im),
     $                       j_global(1), j_global(jm) )
         
         if ( here )  then

            vol_out( i_mabf(i)-i_global(1), 
     $           j_mabf(i)-j_global(1) + 1 ) = riv(i) 
!            print*, i, i_mabf(i)-1, j_mabf(i), riv(i)
            
         endif
!(i,j-1)
         here = judge_inout( i_mabf(i), j_mabf(i)-1, 
     $                       i_global(1), i_global(im),
     $                       j_global(1), j_global(jm) )
         
         if ( here )  then

            vol_out( i_mabf(i)-i_global(1) + 1, 
     $           j_mabf(i)-j_global(1) ) = riv(i) 
!            print*, i, i_mabf(i), j_mabf(i)-1, riv(i)
            
         endif


      enddo

      else
      do i = 1, nr_mab

         here = judge_inout( i_mab(i), j_mab(i), 
     $                       i_global(1), i_global(im),
     $                       j_global(1), j_global(jm) )
         
         if ( here )  then

            vol_out( i_mab(i)-i_global(1) + 1, 
     $           j_mab(i)-j_global(1) + 1 ) = riv(i) 
!            print*, i, i_mab(i), j_mab(i), riv(i)
            
         endif


      enddo

      call exchange2d_mpi(vol_out,im,jm)  !fhx:fgrid

      endif  !if(calc_interp) then !fhx:interp_flag

      
!      endif  !fhx:river_gen
      
      return

      end subroutine river_mab
!-------------------------------------------------------------


!=============================================================
!calculate river dischage according to Julian days,
!    riv =a(1)+a(2)*sin(2*pi*t/T1)+a(3)*cos(2*pi*t/T1)+
!              a(4)*sin(2*pi*t/T2)+a(5)*cos(2*pi*t/T2)+
!              a(6)*sin(2*pi*t/T3)+a(7)*cos(2*pi*t/T3);
!-------------------------------------------------------------
       subroutine river_gen_gom(tt,riv_out)
       integer :: i
       real(kind=rk) tt
       real(kind=rk) riv_out(nr_gom)
       character(len=50) :: coef_file = './in/GOM_rivers_coef.bin'
       real(kind=rk), dimension(numc,nr_gom) :: a 

          T1=365.;T2=182.5;T3=91.25 !periods of harmonic analysis

          open(20, file=coef_file, form='unformatted',
     $     access='direct',recl = (numc*nr_gom)*4 , status = 'old')     
      if ( my_task == master_task )  print'(3a)',
     $      'reading file ',trim(coef_file),' .'
          read(20,rec=1) a
!          if ( my_task == master_task ) 
!     $     print*,'GOM rivers coef:', a(:,15)
          close(20)
 
      do i = 1,nr_gom
       riv_out(i) =a(1,i)+a(2,i)*sin(2*pi*tt/T1)+a(3,i)*cos(2*pi*tt/T1) 
     $                   +a(4,i)*sin(2*pi*tt/T2)+a(5,i)*cos(2*pi*tt/T2)
     $                   +a(6,i)*sin(2*pi*tt/T3)+a(7,i)*cos(2*pi*tt/T3)
      enddo

!          if ( my_task == master_task )  print*,tt,riv_out

       return
       end subroutine river_gen_gom

!-------------------------------------------------------------
       subroutine river_gen_mab(tt,riv_out)
       integer :: i
       real(kind=rk) tt
       real(kind=rk) riv_out(nr_mab)
       character(len=50) :: coef_file = './in/MAB_rivers_coef.bin'
       real(kind=rk), dimension(numc,nr_mab) :: a 

          T1=365.;T2=182.5;T3=91.25 !periods of harmonic analysis

          open(20, file=coef_file, form='unformatted',
     $     access='direct',recl = (numc*nr_mab)*4 , status = 'old')     
      if ( my_task == master_task )  print'(3a)',
     $      'reading file ',trim(coef_file),' .'
          read(20,rec=1) a
!          if ( my_task == master_task ) 
!     $     print*,'MAB rivers coef:', a(:,15)
          close(20)
 
      do i = 1,nr_mab
       riv_out(i) =a(1,i)+a(2,i)*sin(2*pi*tt/T1)+a(3,i)*cos(2*pi*tt/T1) 
     $                   +a(4,i)*sin(2*pi*tt/T2)+a(5,i)*cos(2*pi*tt/T2)
     $                   +a(6,i)*sin(2*pi*tt/T3)+a(7,i)*cos(2*pi*tt/T3)
      enddo

!          if ( my_task == master_task )  print*,tt,riv_out

       return
       end subroutine river_gen_mab

!=============================================================
! If the processor has ( i_in, j_in ) in its local domain.
!-------------------------------------------------------------
      logical function judge_inout( i_in, j_in, 
     $                              imin_in, imax_in, 
     $                              jmin_in, jmax_in ) 
      integer :: i_in, j_in, imin_in, imax_in, jmin_in, jmax_in
      
      if ( ( i_in >= imin_in .and. i_in <= imax_in )
     $     .and.( j_in >= jmin_in .and. j_in <= jmax_in ) ) then

         judge_inout = .true.

      else

         judge_inout = .false.

      endif

      end function judge_inout
!-------------------------------------------------------------

      end module river

