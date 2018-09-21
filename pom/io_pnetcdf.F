! Updated (in/grid/grid.nc & >2GB) io_pnetcdf.F combining FH+Alu
!     This version is for FH's case with AluUpperSST & AluUpperSSS
!     commented out
!
! io_pnetcdf.F !lyo:pac10:Incorporating Alu's version for pac10
!              !          in ncu's /home/lyo/gaea/alu/pac10

! parallel NetCDF I/O

!_______________________________________________________________________
      subroutine def_var_pnetcdf(ncid,name,nvdims,vdims,varid,vartype
     &                          ,long_name,units,nofill,fillval
     &                          ,coords,lcoords)
! defines variable in netcdf file
!        use mpi,     only: mpi_offset_kind, mpi_info_null
        use pnetcdf, only: nf90mpi_def_var, nf90mpi_put_att
     &                   , nf90mpi_def_var_fill
     &                   , NF_NOERR

        integer         , intent(in) :: vdims(nvdims)
        integer         , intent(in) :: ncid,nvdims,vartype,nofill
        logical         , intent(in) :: lcoords
        character(len=*), intent(in) :: name,long_name,units,coords
        real            , intent(in) :: fillval
        integer         , intent(out):: varid
        integer status

! define variable
        status = nf90mpi_def_var( ncid, name, vartype, vdims, varid )
        call handle_error_pnetcdf( 'nf_def_var', status )
! define fill value
        if ( nofill > -1 ) then
          status = nf90mpi_def_var_fill( ncid, varid, nofill, fillval )
          call handle_error_pnetcdf( 'nf_def_var_fill', status )
        end if
! define attributes
        status = nf90mpi_put_att( ncid, varid, 'long_name'
     &                           ,trim(long_name) )
        call handle_error_pnetcdf( 'nf_put_att : long_name', status )

        status = nf90mpi_put_att( ncid, varid, 'units'
     &                           ,trim(units) )
        call handle_error_pnetcdf( 'nf_put_att : units', status )
! add coordinates attribute, if necessary
        if ( lcoords ) then
          status = nf90mpi_put_att( ncid, varid, 'coordinates'
     &                             ,trim(coords) )
          call handle_error_pnetcdf( 'nf_put_att : coordinates'
     &                             , status )
        end if

        return

      end

!_______________________________________________________________________
      subroutine handle_error_pnetcdf(routine,status)

        use pnetcdf, only: nf90mpi_strerror, NF_NOERR

        implicit none
        include 'pom.h'

        integer, intent(in) :: status
        character(len=*), intent(in) :: routine

        if ( status /= NF_NOERR ) then
          error_status = 1
          if ( my_task == master_task ) then
             print '(/a,a)', 'Error: NetCDF routine ', routine
             print *, nf90mpi_strerror(status)
             stop
          endif
        end if

        return

      end

!_______________________________________________________________________
!      subroutine write_output_pnetcdf
! ayumi 2010/5/10
      subroutine write_output_pnetcdf0( netcdf_out_file )
! write output data
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      character(len=*), intent(in) :: netcdf_out_file
      character(len=120) str_tmp!, netcdf_out_file
      integer time_dimid,x_dimid,y_dimid,z_dimid
      integer z_varid,zz_varid,dx_varid,dy_varid,east_c_varid,
     $        east_e_varid,east_u_varid,east_v_varid,north_c_varid,
     $        north_e_varid,north_u_varid,north_v_varid,rot_varid,
     $        h_varid,fsm_varid,dum_varid,dvm_varid,
     $        wusurf_varid,wvsurf_varid,wtsurf_varid,wssurf_varid
      integer time_varid,uab_varid,vab_varid,elb_varid,u_varid,v_varid,
     $        w_varid,t_varid,s_varid,rho_varid,km_varid,kh_varid
      integer ncid,status
      integer vdims(4)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(4),edge(4)


! create netcdf file
!      nprint=(iint+time0*86400/dti)/iprint
!      write(netcdf_out_file,'(''out/'',a,''.'',i4.4,''.nc'')')
!     $                                          trim(netcdf_file),nprint
      if(my_task.eq.0) write(*,'(/''writing file '',a)')
     $     trim(netcdf_out_file)
      status=nfmpi_create(pom_comm,trim(netcdf_out_file),
     $     nf_clobber+nf_64bit_offset,mpi_info_null,ncid) !lyo:pac10:
      call handle_error_pnetcdf( 'nf_create: '//netcdf_out_file
     &                         , status )

! define global attributes
      length=len(trim(title))
      status=nfmpi_put_att_text(ncid,nf_global,'title',length,
     $                          trim(title))
      call handle_error_pnetcdf( 'nf_put_att_text', status )

      str_tmp='output file'
      length=len(trim(str_tmp))
      status=nfmpi_put_att_text(ncid,nf_global,'description',length,
     $                          trim(str_tmp))
      call handle_error_pnetcdf( 'nf_put_att_text', status )

! define dimensions
      length=1
      status=nfmpi_def_dim(ncid,'time',length,time_dimid)
      call handle_error_pnetcdf( 'nf_def_dim: time', status )
      length=kb
      status=nfmpi_def_dim(ncid,'z',length,z_dimid)
      call handle_error_pnetcdf( 'nf_def_dim: z', status )
      length=jm_global
      status=nfmpi_def_dim(ncid,'y',length,y_dimid)
      call handle_error_pnetcdf( 'nf_def_dim: y', status )
      length=im_global
      status=nfmpi_def_dim(ncid,'x',length,x_dimid)
      call handle_error_pnetcdf( 'nf_def_dim: x', status )

! define variables and their attributes
      vdims(1)=time_dimid
      str_tmp='days since '//time_start
      call def_var_pnetcdf(ncid,'time',1,vdims,time_varid,nf_float,
     $                     'time',str_tmp,
     $                     -1,0.,' ',.false.)

      vdims(1)=z_dimid
      call def_var_pnetcdf(ncid,'z',1,vdims,z_varid,nf_float,
     $                     'sigma of cell face','sigma_level',
     $                     -1,0.,' ',.false.)
      length=22
      status=nfmpi_put_att_text(ncid,z_varid,'standard_name',
     $                          length,'ocean_sigma_coordinate')
      length=26
      status=nfmpi_put_att_text(ncid,z_varid,'formula_terms',
     $                          length,'sigma: z eta: elb depth: h')
      call handle_error_pnetcdf( 'nf_put_att_text: z', status )

      call def_var_pnetcdf(ncid,'zz',1,vdims,zz_varid,nf_float,
     $                     'sigma of cell centre','sigma_level',
     $                     -1,0.,' ',.false.)
      length=22
      status=nfmpi_put_att_text(ncid,zz_varid,'standard_name',
     $                          length,'ocean_sigma_coordinate')
      length=27
      status=nfmpi_put_att_text(ncid,zz_varid,'formula_terms',
     $                          length,'sigma: zz eta: elb depth: h')
      call handle_error_pnetcdf( 'nf_put_att_text: zz', status )

      vdims(1)=x_dimid
      vdims(2)=y_dimid
      call def_var_pnetcdf(ncid,'dx',2,vdims,dx_varid,nf_float,
     $                     'grid increment in x','metre',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'dy',2,vdims,dy_varid,nf_float,
     $                     'grid increment in y','metre',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'east_u',2,vdims,east_u_varid,nf_float,
     $                     'easting of u-points','degree',
     $                     -1,0.,'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'east_v',2,vdims,east_v_varid,nf_float,
     $                     'easting of v-points','degree',
     $                     -1,0.,'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'east_e',2,vdims,east_e_varid,nf_float,
     $                     'easting of elevation points','degree',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'east_c',2,vdims,east_c_varid,nf_float,
     $                     'easting of cell corners','degree',
     $                     -1,0.,'east_c north_c',.true.)
      call def_var_pnetcdf(ncid,'north_u',2,vdims,north_u_varid,nf_float
     $                    ,'northing of u-points','degree',
     $                     -1,0.,'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'north_v',2,vdims,north_v_varid,nf_float
     $                    ,'northing of v-points','degree',
     $                     -1,0.,'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'north_e',2,vdims,north_e_varid,nf_float
     $                    ,'northing of elevation points','degree',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'north_c',2,vdims,north_c_varid,nf_float
     $                    ,'northing of cell corners','degree',
     $                     -1,0.,'east_c north_c',.true.)
      call def_var_pnetcdf(ncid,'rot',2,vdims,rot_varid,nf_float,
     $                    'Rotation angle of x-axis wrt. east','degree',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'h',2,vdims,h_varid,nf_float,
     $                     'undisturbed water depth','metre',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'fsm',2,vdims,fsm_varid,nf_float,
     $                     'free surface mask','dimensionless',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'dum',2,vdims,dum_varid,nf_float,
     $                     'u-velocity mask','dimensionless',
     $                     -1,0.,'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'dvm',2,vdims,dvm_varid,nf_float,
     $                     'v-velocity mask','dimensionless',
     $                     -1,0.,'east_v north_v',.true.)

      vdims(1)=x_dimid
      vdims(2)=y_dimid
      vdims(3)=time_dimid
      call def_var_pnetcdf(ncid,'uab',3,vdims,uab_varid,nf_float,
     $                     'depth-averaged u','metre/sec',
     $                     -1,0.,'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'vab',3,vdims,vab_varid,nf_float,
     $                     'depth-averaged v','metre/sec',
     $                     -1,0.,'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'elb',3,vdims,elb_varid,nf_float,
     $                     'surface elevation','metre',
     $                     0,-999.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'wusurf',3,vdims,wusurf_varid,nf_float,
     $                     'x-momentum flux','metre^2/sec^2',
     $                     -1,0.,'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'wvsurf',3,vdims,wvsurf_varid,nf_float,
     $                     'y-momentum flux','metre^2/sec^2',
     $                     -1,0.,'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'wtsurf',3,vdims,wtsurf_varid,nf_float,
     $                     'temperature flux','deg m/s',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'wssurf',3,vdims,wssurf_varid,nf_float,
     $                     'salinity flux','psu m/s',
     $                     -1,0.,'east_e north_e',.true.)
      vdims(1)=x_dimid
      vdims(2)=y_dimid
      vdims(3)=z_dimid
      vdims(4)=time_dimid
      call def_var_pnetcdf(ncid,'u',4,vdims,u_varid,nf_float,
     $                     'x-velocity','metre/sec',
     $                     -1,0.,'east_u north_u zz',.true.)
      call def_var_pnetcdf(ncid,'v',4,vdims,v_varid,nf_float,
     $                     'y-velocity','metre/sec',
     $                     -1,0.,'east_v north_v zz',.true.)
      call def_var_pnetcdf(ncid,'w',4,vdims,w_varid,nf_float,
     $                     'z-velocity','metre/sec',
     $                     -1,0.,'east_e north_e z',.true.)
      call def_var_pnetcdf(ncid,'t',4,vdims,t_varid,nf_float,
     $                     'potential temperature','K',
     $                     0,-999.,'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'s',4,vdims,s_varid,nf_float,
     $                     'salinity x rho / rhoref','PSS',
     $                     0,0.,'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'rho',4,vdims,rho_varid,nf_float,
     $                     '(density-1000)/rhoref','dimensionless',
     $                     -1,0.,'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'kh',4,vdims,kh_varid,nf_float,
     $                     'vertical diffusivity','metre2/sec',
     $                     -1,0.,'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'km',4,vdims,km_varid,nf_float,
     $                     'vertical viscosity','metre2/sec',
     $                     -1,0.,'east_e north_e zz',.true.)

! end definitions
      status=nfmpi_enddef(ncid)
      call handle_error_pnetcdf( 'nf_enddef: output', status )

! write data
      start(1)=1
      edge(1)=1
      status=nfmpi_put_vara_real_all(ncid,time_varid,start,edge
     $                                                  ,real(time,4))
      call handle_error_pnetcdf( 'nf_put_vara_real:time', status )

      start(1)=1
      edge(1)=kb
      status=nfmpi_put_vara_real_all(ncid,z_varid,start,edge
     $                                                     ,real(z,4))
      call handle_error_pnetcdf( 'nf_put_var_real: z', status )
      status=nfmpi_put_vara_real_all(ncid,zz_varid,start,edge
     $                                                    ,real(zz,4))
      call handle_error_pnetcdf( 'nf_put_var_real: zz', status )

      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      status=nfmpi_put_vara_real_all(ncid,dx_varid,start,edge
     $                                                    ,real(dx,4))
      call handle_error_pnetcdf( 'nf_put_var_real: dx', status )
      status=nfmpi_put_vara_real_all(ncid,dy_varid,start,edge
     $                                                    ,real(dy,4))
      call handle_error_pnetcdf( 'nf_put_var_real: dy', status )
      status=nfmpi_put_vara_real_all(ncid,east_u_varid,start,edge,
     $                                                 real(east_u,4))
      call handle_error_pnetcdf( 'nf_put_var_real: east_u', status )
      status=nfmpi_put_vara_real_all(ncid,east_v_varid,start,edge,
     $                                                 real(east_v,4))
      call handle_error_pnetcdf( 'nf_put_var_real: east_v', status )
      status=nfmpi_put_vara_real_all(ncid,east_e_varid,start,edge,
     $                                                 real(east_e,4))
      call handle_error_pnetcdf( 'nf_put_var_real: east_e', status )
      status=nfmpi_put_vara_real_all(ncid,east_c_varid,start,edge,
     $                                                 real(east_c,4))
      call handle_error_pnetcdf( 'nf_put_var_real: east_c', status )
      status=nfmpi_put_vara_real_all(ncid,north_u_varid,start,edge,
     $                                                real(north_u,4))
      call handle_error_pnetcdf( 'nf_put_var_real: north_u', status )
      status=nfmpi_put_vara_real_all(ncid,north_v_varid,start,edge,
     $                                                real(north_v,4))
      call handle_error_pnetcdf( 'nf_put_var_real: north_v', status )
      status=nfmpi_put_vara_real_all(ncid,north_e_varid,start,edge,
     $                                                real(north_e,4))
      call handle_error_pnetcdf( 'nf_put_var_real: north_e', status )
      status=nfmpi_put_vara_real_all(ncid,north_c_varid,start,edge,
     $                                                real(north_c,4))
      call handle_error_pnetcdf( 'nf_put_var_real: north_c', status )
      status=nfmpi_put_vara_real_all(ncid,rot_varid,start,edge
     $                                                   ,real(rot,4))
      call handle_error_pnetcdf( 'nf_put_var_real: rot', status )
      status=nfmpi_put_vara_real_all(ncid,h_varid,start,edge
     $                                                      ,real(h,4))
      call handle_error_pnetcdf( 'nf_put_var_real: h', status )
      status=nfmpi_put_vara_real_all(ncid,fsm_varid,start,edge
     $                                                    ,real(fsm,4))
      call handle_error_pnetcdf( 'nf_put_var_real: fsm', status )
      status=nfmpi_put_vara_real_all(ncid,dum_varid,start,edge
     $                                                    ,real(dum,4))
      call handle_error_pnetcdf( 'nf_put_var_real: dum', status )
      status=nfmpi_put_vara_real_all(ncid,dvm_varid,start,edge
     $                                                    ,real(dvm,4))
      call handle_error_pnetcdf( 'nf_put_var_real: dvm', status )

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=1
      status=nfmpi_put_vara_real_all(ncid,uab_varid,start,edge
     $                                               ,real(uab_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: uab', status )
      status=nfmpi_put_vara_real_all(ncid,vab_varid,start,edge
     $                                               ,real(vab_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: vab', status )
      status=nfmpi_put_vara_real_all(ncid,elb_varid,start,edge
     $                                               ,real(elb_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: elb', status )
      status=nfmpi_put_vara_real_all(
     $     ncid,wusurf_varid,start,edge,real(wusurf_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: wusurf', status )
      status=nfmpi_put_vara_real_all(
     $     ncid,wvsurf_varid,start,edge,real(wvsurf_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: wvsurf', status )
      status=nfmpi_put_vara_real_all(
     $     ncid,wtsurf_varid,start,edge,real(wtsurf_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: wtsurf', status )
      status=nfmpi_put_vara_real_all(
     $     ncid,wssurf_varid,start,edge,real(wssurf_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: wssurf', status )

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      start(4)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=kb
      edge(4)=1
      status=nfmpi_put_vara_real_all(ncid,u_varid,start,edge
     $                                                ,real(u_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: u', status )
      status=nfmpi_put_vara_real_all(ncid,v_varid,start,edge
     $                                                ,real(v_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: v', status )
      status=nfmpi_put_vara_real_all(ncid,w_varid,start,edge
     $                                                ,real(w_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: w', status )
      status=nfmpi_put_vara_real_all(ncid,t_varid,start,edge
     $                                                ,real(t_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: t', status )
      status=nfmpi_put_vara_real_all(ncid,s_varid,start,edge
     $                                                ,real(s_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: s', status )
      status=nfmpi_put_vara_real_all(ncid,rho_varid,start,edge
     $                                              ,real(rho_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: rho', status )
      status=nfmpi_put_vara_real_all(ncid,kh_varid,start,edge
     $                                               ,real(kh_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: kh', status )
      status=nfmpi_put_vara_real_all(ncid,km_varid,start,edge
     $                                               ,real(km_mean,4))
      call handle_error_pnetcdf( 'nf_put_vara_real: km', status )

! close file
      status = nfmpi_close(ncid)
      call handle_error_pnetcdf( 'nf_close: output', status )

      return
      end

!_______________________________________________________________________
!      subroutine write_restart_pnetcdf
      subroutine write_restart_pnetcdf( netcdf_out_file )
! write data for a seamless restart
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      character(len=*), intent(in) :: netcdf_out_file
      character(len=120) str_tmp
      integer nprint
      integer ncid,status
      integer time_dimid,x_dimid,y_dimid,z_dimid
      integer time_varid,wubot_varid,wvbot_varid,aam2d_varid,ua_varid,
     $        wusurf_varid, wvsurf_varid,
     $        uab_varid,va_varid,vab_varid,el_varid,elb_varid,et_varid,
     $        etb_varid,egb_varid,utb_varid,vtb_varid,u_varid,ub_varid,
     $        w_varid,v_varid,vb_varid,t_varid,tb_varid,s_varid,
     $        sb_varid,rho_varid,adx2d_varid,ady2d_varid,advua_varid,
     $        advva_varid,km_varid,kh_varid,kq_varid,l_varid,q2_varid,
     $        q2b_varid,aam_varid,q2l_varid,q2lb_varid
      integer vdims(3)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(3),edge(3)
      logical spinup
      namelist/misc_nml/ spinup, t_lo, t_hi

      open(73, file='switch.nml',status='old')
      read(73, nml=misc_nml)
      close(73)

! create netcdf restart file
      nprint=int((iint+time0*86400./dti)/irestart)
!      write(netcdf_out_file,'(''out/'',a,''.'',i4.4,''.nc'')')
!     $     trim(write_rst_file),nprint
      if(my_task.eq.0) write(*,'(/''writing file '',a)')
     $ trim(netcdf_out_file)
!      status=nfmpi_create(pom_comm,netcdf_out_file,nf_64bit_offset,
      status=nfmpi_create(pom_comm,netcdf_out_file,
     $     nf_clobber+nf_64bit_offset,mpi_info_null,ncid) !lyo:pac10:
      call handle_error_pnetcdf( 'nf_create', status )

! define global attributes
      length=len(trim(title))
      status=nfmpi_put_att_text(ncid,nf_global,'title',length,
     $                          trim(title))
      call handle_error_pnetcdf( 'nf_put_att_text', status )

      str_tmp='restart file'
      length=len(trim(str_tmp))
      status=nfmpi_put_att_text(ncid,nf_global,'description',length,
     $                          trim(str_tmp))
      call handle_error_pnetcdf( 'nf_put_att_text', status )

! define dimensions
      length=1
      status=nfmpi_def_dim(ncid,'time',length,time_dimid)
      call handle_error_pnetcdf( 'nf_def_dim', status )

      length=kb
      status=nfmpi_def_dim(ncid,'z',length,z_dimid)
      call handle_error_pnetcdf( 'nf_def_dim', status )

      length=jm_global
      status=nfmpi_def_dim(ncid,'y',length,y_dimid)
      call handle_error_pnetcdf( 'nf_def_dim', status )

      length=im_global
      status=nfmpi_def_dim(ncid,'x',length,x_dimid)
      call handle_error_pnetcdf( 'nf_def_dim', status )

! define variables and their attributes:
      vdims(1)=time_dimid
      str_tmp='days since '//time_start
      call def_var_pnetcdf(ncid,'time',1,vdims,time_varid,nf_double,
     $                     'time',str_tmp,-1,0.,' ',.false.)

      vdims(1)=x_dimid
      vdims(2)=y_dimid
      call def_var_pnetcdf(ncid,'wusurf',2,vdims,wusurf_varid,nf_double,
     $                     'x-momentum flux at the surface',
     $                     'metre^2/sec^2',-1,0.,
     $                     'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'wvsurf',2,vdims,wvsurf_varid,nf_double,
     $                     'y-momentum flux at the surface',
     $                     'metre^2/sec^2',-1,0.,
     $                     'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'wubot',2,vdims,wubot_varid,nf_double,
     $                     'x-momentum flux at the bottom',
     $                     'metre^2/sec^2',-1,0.,
     $                     'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'wvbot',2,vdims,wvbot_varid,nf_double,
     $                     'y-momentum flux at the bottom',
     $                     'metre^2/sec^2',-1,0.,
     $                     'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'aam2d',2,vdims,aam2d_varid,nf_double,
     $                     'vertical average of aam',
     $                     'metre^2/sec',-1,0.,
     $                     'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'ua',2,vdims,ua_varid,nf_double,
     $                     'vertical mean of u',
     $                     'metre/sec',-1,0.,
     $                     'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'uab',2,vdims,uab_varid,nf_double,
     $                     'vertical mean of u at time -dt',
     $                     'metre/sec',-1,0.,
     $                     'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'va',2,vdims,va_varid,nf_double,
     $                     'vertical mean of v',
     $                     'metre/sec',-1,0.,
     $                     'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'vab',2,vdims,vab_varid,nf_double,
     $                     'vertical mean of v at time -dt',
     $                     'metre/sec',-1,0.,
     $                     'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'el',2,vdims,el_varid,nf_double,
     $                     'surface elevation in external mode',
     $                     'metre',-1,0.,
     $                     'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'elb',2,vdims,elb_varid,nf_double,
     $                     'surface elevation in external mode at -dt',
     $                     'metre',-1,0.,
     $                     'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'et',2,vdims,et_varid,nf_double,
     $                     'surface elevation in internal mode',
     $                     'metre',-1,0.,
     $                     'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'etb',2,vdims,etb_varid,nf_double,
     $                     'surface elevation in internal mode at -dt',
     $                     'metre',-1,0.,
     $                     'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'egb',2,vdims,egb_varid,nf_double,
     $                     'surface elevation for pres. grad. at -dt',
     $                     'metre',-1,0.,
     $                     'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'utb',2,vdims,utb_varid,nf_double,
     $                     'ua time averaged over dti',
     $                     'metre/sec',-1,0.,
     $                     'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'vtb',2,vdims,vtb_varid,nf_double,
     $                     'va time averaged over dti',
     $                     'metre/sec',-1,0.,
     $                     'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'adx2d',2,vdims,adx2d_varid,nf_double,
     $                     'vertical integral of advx',
     $                     '-',-1,0.,
     $                     'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'ady2d',2,vdims,ady2d_varid,nf_double,
     $                     'vertical integral of advy',
     $                     '-',-1,0.,
     $                     'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'advua',2,vdims,advua_varid,nf_double,
     $                     'sum of 2nd, 3rd and 4th terms in eq (18)',
     $                     '-',-1,0.,
     $                     'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'advva',2,vdims,advva_varid,nf_double,
     $                     'sum of 2nd, 3rd and 4th terms in eq (19)',
     $                     '-',-1,0.,
     $                     'east_v north_v',.true.)

      vdims(1)=x_dimid
      vdims(2)=y_dimid
      vdims(3)=z_dimid
      call def_var_pnetcdf(ncid,'u',3,vdims,u_varid,nf_double,
     $                     'x-velocity',
     $                     'metre/sec',-1,0.,
     $                     'east_u north_u zz',.true.)
      call def_var_pnetcdf(ncid,'ub',3,vdims,ub_varid,nf_double,
     $                     'x-velocity at time -dt',
     $                     'metre/sec',-1,0.,
     $                     'east_u north_u zz',.true.)
      call def_var_pnetcdf(ncid,'v',3,vdims,v_varid,nf_double,
     $                     'y-velocity',
     $                     'metre/sec',-1,0.,
     $                     'east_v north_v zz',.true.)
      call def_var_pnetcdf(ncid,'vb',3,vdims,vb_varid,nf_double,
     $                     'y-velocity at time -dt',
     $                     'metre/sec',-1,0.,
     $                     'east_v north_v zz',.true.)
      call def_var_pnetcdf(ncid,'w',3,vdims,w_varid,nf_double,
     $                     'sigma-velocity',
     $                     'metre/sec',-1,0.,
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'t',3,vdims,t_varid,nf_double,
     $                     'potential temperature',
     $                     'K',-1,0.,
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'tb',3,vdims,tb_varid,nf_double,
     $                     'potential temperature at time -dt',
     $                     'K',-1,0.,
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'s',3,vdims,s_varid,nf_double,
     $                     'salinity x rho / rhoref',
     $                     'PSS',-1,0.,
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'sb',3,vdims,sb_varid,nf_double,
     $                     'salinity x rho / rhoref at time -dt',
     $                     'PSS',-1,0.,
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'rho',3,vdims,rho_varid,nf_double,
     $                     '(density-1000)/rhoref',
     $                     'dimensionless',-1,0.,
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'km',3,vdims,km_varid,nf_double,
     $                     'vertical kinematic viscosity',
     $                     'metre^2/sec',-1,0.,
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'kh',3,vdims,kh_varid,nf_double,
     $                     'vertical diffusivity',
     $                     'metre^2/sec',-1,0.,
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'kq',3,vdims,kq_varid,nf_double,
     $                     'kq',
     $                     'metre^2/sec',-1,0.,
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'l',3,vdims,l_varid,nf_double,
     $                     'turbulence length scale',
     $                     '-',-1,0.,
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'q2',3,vdims,q2_varid,nf_double,
     $                     'twice the turbulent kinetic energy',
     $                     'metre^2/sec^2',-1,0.,
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'q2b',3,vdims,q2b_varid,nf_double,
     $                     'twice the turbulent kinetic energy at -dt',
     $                     'metre^2/sec^2',-1,0.,
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'aam',3,vdims,aam_varid,nf_double,
     $                     'horizontal kinematic viscosity',
     $                     'metre^2/sec',-1,0.,
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'q2l',3,vdims,q2l_varid,nf_double,
     $                     'q2 x l',
     $                     'metre^3/sec^2',-1,0.,
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'q2lb',3,vdims,q2lb_varid,nf_double,
     $                     'q2 x l at time -dt',
     $                     'metre^3/sec^2',-1,0.,
     $                     'east_e north_e zz',.true.)

! end definitions
      status=nfmpi_enddef(ncid)
      call handle_error_pnetcdf( 'nf_enddef', status )

! write data
      start(1)=1
      edge(1)=1
      if (spinup) then
        status=nfmpi_put_vara_double_all(ncid,time_varid,start,edge
     &                                                   ,real(time0,8))
      else
        status=nfmpi_put_vara_double_all(ncid,time_varid,start,edge
     &                                                    ,real(time,8))
      end if
      call handle_error_pnetcdf( 'nf_put_vara_double', status )

      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      status=nfmpi_put_vara_double_all(
     $     ncid,wusurf_varid,start,edge,real(wusurf,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(
     $     ncid,wvsurf_varid,start,edge,real(wvsurf,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,wubot_varid,start,edge
     $ ,real(wubot,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,wvbot_varid,start,edge
     $ ,real(wvbot,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,aam2d_varid,start,edge
     $ ,real(aam2d,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,ua_varid,start,edge
     &                                                      ,real(ua,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,uab_varid,start,edge
     &                                                     ,real(uab,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,va_varid,start,edge
     &                                                      ,real(va,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,vab_varid,start,edge
     &                                                     ,real(vab,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,el_varid,start,edge
     &                                                      ,real(el,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,elb_varid,start,edge
     &                                                     ,real(elb,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,et_varid,start,edge
     &                                                      ,real(et,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,etb_varid,start,edge
     &                                                     ,real(etb,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,egb_varid,start,edge
     &                                                     ,real(egb,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,utb_varid,start,edge
     &                                                     ,real(utb,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,vtb_varid,start,edge
     &                                                     ,real(vtb,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,adx2d_varid,start,edge
     $ ,real(adx2d,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,ady2d_varid,start,edge
     $ ,real(ady2d,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,advua_varid,start,edge
     $ ,real(advua,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,advva_varid,start,edge
     $ ,real(advva,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=kb
      status=nfmpi_put_vara_double_all(ncid,u_varid,start,edge
     &                                                       ,real(u,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,ub_varid,start,edge
     &                                                      ,real(ub,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,v_varid,start,edge
     &                                                       ,real(v,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,vb_varid,start,edge
     &                                                      ,real(vb,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,w_varid,start,edge
     &                                                       ,real(w,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,t_varid,start,edge
     &                                                       ,real(t,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,tb_varid,start,edge
     &                                                      ,real(tb,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,s_varid,start,edge
     &                                                       ,real(s,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,sb_varid,start,edge
     &                                                      ,real(sb,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,rho_varid,start,edge
     &                                                     ,real(rho,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,km_varid,start,edge
     &                                                      ,real(km,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,kh_varid,start,edge
     &                                                      ,real(kh,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,kq_varid,start,edge
     &                                                      ,real(kq,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,l_varid,start,edge
     &                                                       ,real(l,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,q2_varid,start,edge
     &                                                      ,real(q2,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,q2b_varid,start,edge
     &                                                     ,real(q2b,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,aam_varid,start,edge
     &                                                     ,real(aam,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,q2l_varid,start,edge
     &                                                     ,real(q2l,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )
      status=nfmpi_put_vara_double_all(ncid,q2lb_varid,start,edge
     $ ,real(q2lb,8))
      call handle_error_pnetcdf( 'nf_put_vara_real', status )

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf( 'nf_close', status )

      return
      end

!_______________________________________________________________________
      subroutine read_grid_pnetcdf
! read grid data

      use mpi
      use pnetcdf

      implicit none

      include 'pom.h'

      integer get_var_real_1d
     &       ,get_var_real_2d
     &       ,get_var_real_3d

      character(len=120) netcdf_grid_file
      integer z_varid,zz_varid,dx_varid,dy_varid,east_c_varid,
     $        east_e_varid,east_u_varid,east_v_varid,north_c_varid,
     $        north_e_varid,north_u_varid,north_v_varid,rot_varid,
     $        h_varid,fsm_varid!,dum_varid,dvm_varid
      integer ncid,status,dimn
      integer i,j
      integer(mpi_offset_kind) start(4),edge(4)
!      real(kind=rk) minh

      start = 1
      edge  = 1

! open netcdf file
      write(netcdf_grid_file,'(a)') "./in/grid/grid.nc"

      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_grid_file)
      status=nfmpi_open(pom_comm,netcdf_grid_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_grid_file
     $                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'z',z_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: z',  status )
      status=nfmpi_inq_varid(ncid,'zz',zz_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: zz', status )
      status=nfmpi_inq_varid(ncid,'dx',dx_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: dx', status )
      status=nfmpi_inq_varid(ncid,'dy',dy_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: dy', status )
      status=nfmpi_inq_varid(ncid,'lon_u',east_u_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: east_u', status )
      status=nfmpi_inq_varid(ncid,'lon_v',east_v_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: east_v'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'lon_e',east_e_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: east_e'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'lon_c',east_c_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: east_c'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'lat_u',north_u_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: north_u'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'lat_v',north_v_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: north_v'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'lat_e',north_e_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: north_e'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'lat_c',north_c_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: north_c'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'rot',rot_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: rot'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'h',h_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: h'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'fsm',fsm_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: fsm'
     $                          ,status)
!      status=nfmpi_inq_varid(ncid,'dum',dum_varid)
!      call handle_error_pnetcdf('nfmpi_inq_varid: dum'
!     $                          ,status)
!      status=nfmpi_inq_varid(ncid,'dvm',dvm_varid)
!      call handle_error_pnetcdf('nfmpi_inq_varid: dvm'
!     $                          ,status)
! get data
      status = nfmpi_inq_varndims(ncid, z_varid, dimn)
      call handle_error_pnetcdf('inq_vardimn: z' ,status)

      if (dimn==1) then
        start(1)=1
        edge(1)=kb
        status = get_var_real_1d(ncid, z_varid,start,edge, z)
        call handle_error_pnetcdf('get_var_real: z' ,status)
        status = get_var_real_1d(ncid,zz_varid,start,edge,zz)
        call handle_error_pnetcdf('get_var_real: zz',status)
      else if (dimn==3) then
! If z and zz are 3-dimensional get just the (1,1) cell distribution
        start(1)=1
        start(2)=1
        start(3)=1
        edge(1)=1
        edge(2)=1
        edge(3)=kb
        status = get_var_real_3d(ncid, z_varid,start,edge, z)
        call handle_error_pnetcdf('get_var_real: z' ,status)
        status = get_var_real_3d(ncid,zz_varid,start,edge,zz)
        call handle_error_pnetcdf('get_var_real: zz',status)
      end if

      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      status = get_var_real_2d(ncid,dx_varid,start,edge,dx)
      call handle_error_pnetcdf('get_var_real: dx',status)
      status = get_var_real_2d(ncid,dy_varid,start,edge,dy)
      call handle_error_pnetcdf('get_var_real: dy',status)
      status = get_var_real_2d(ncid,east_u_varid,start,edge,east_u)
      call handle_error_pnetcdf('get_var_real: east_u',status)
      status = get_var_real_2d(ncid,east_v_varid,start,edge,east_v)
      call handle_error_pnetcdf('get_var_real: east_v',status)
      status = get_var_real_2d(ncid,east_e_varid,start,edge,east_e)
      call handle_error_pnetcdf('get_var_real: east_e',status)
      status = get_var_real_2d(ncid,east_c_varid,start,edge,east_c)
      call handle_error_pnetcdf('get_var_real: east_c',status)
      status = get_var_real_2d(ncid,north_u_varid,start,edge,north_u)
      call handle_error_pnetcdf('get_var_real: north_u',status)
      status = get_var_real_2d(ncid,north_v_varid,start,edge,north_v)
      call handle_error_pnetcdf('get_var_real: north_v',status)
      status = get_var_real_2d(ncid,north_e_varid,start,edge,north_e)
      call handle_error_pnetcdf('get_var_real: north_e',status)
      status = get_var_real_2d(ncid,north_c_varid,start,edge,north_c)
      call handle_error_pnetcdf('get_var_real: north_c',status)
      status = get_var_real_2d(ncid,rot_varid,start,edge,rot)
      call handle_error_pnetcdf('get_var_real: rot',status)
      status = get_var_real_2d(ncid,h_varid,start,edge,h)
      call handle_error_pnetcdf('get_var_real: h',status)
      status = get_var_real_2d(ncid,fsm_varid,start,edge,fsm)
      call handle_error_pnetcdf('get_var_real: fsm',status)

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
      call exchange2d_mpi(dum,im_local,jm_local)
      call exchange2d_mpi(dvm,im_local,jm_local)

! close file:
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close: grid',status)

      return
      end
!_______________________________________________________________________
      subroutine read_grid_pnetcdf_obs
! read grid data
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'

      integer get_var_real_1d
     &       ,get_var_real_2d

      character(len=120) netcdf_grid_file
      integer z_varid,zz_varid,dx_varid,dy_varid,east_e_varid,
     $        north_e_varid,rot_varid,h_varid,fsm_varid
      integer ncid,status
      integer i,j
      integer(mpi_offset_kind) start(4),edge(4)

      start = 1
      edge  = 1

! open netcdf file
      write(netcdf_grid_file,'(a)') "./in/grid/grid_old.nc"

      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_grid_file)
      status=nfmpi_open(pom_comm,netcdf_grid_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_grid_file
     &                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'z',z_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: z'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'zz',zz_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: zz'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'dx',dx_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: dx'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'dy',dy_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: dy'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'longitude',east_e_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: east_e'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'latitude',north_e_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: north_e'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'rot',rot_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: rot'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'h',h_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: h'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'fsm',fsm_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: fsm'
     $                          ,status)
! get data
      start(1)=1
      edge(1)=kb
      status = get_var_real_1d(ncid, z_varid,start,edge, z)
      call handle_error_pnetcdf('get_var_real: z' ,status)
      status = get_var_real_1d(ncid,zz_varid,start,edge,zz)
      call handle_error_pnetcdf('get_var_real: zz',status)

      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      status = get_var_real_2d(ncid,dx_varid,start,edge,dx)
      call handle_error_pnetcdf('get_var_real: dx',status)
      status = get_var_real_2d(ncid,dy_varid,start,edge,dy)
      call handle_error_pnetcdf('get_var_real: dy',status)
      status = get_var_real_2d(ncid,east_e_varid,start,edge,east_e)
      call handle_error_pnetcdf('get_var_real: east_e',status)
      status = get_var_real_2d(ncid,north_e_varid,start,edge,north_e)
      call handle_error_pnetcdf('get_var_real: north_e',status)
      status = get_var_real_2d(ncid,rot_varid,start,edge,rot)
      call handle_error_pnetcdf('get_var_real: rot',status)
      status = get_var_real_2d(ncid,h_varid,start,edge,h)
      call handle_error_pnetcdf('get_var_real: h',status)
      status = get_var_real_2d(ncid,fsm_varid,start,edge,fsm)
      call handle_error_pnetcdf('get_var_real: fsm',status)
!
!     u-points:
!
      do j=1,jm-1
        do i=1,im
           east_u(i,j) = ( east_e(i,j)+ east_e(i,j+1))/2.
          north_u(i,j) = (north_e(i,j)+north_e(i,j+1))/2.
        end do
      end do
!
!     Extrapolate ends:
!
      do i=1,im
         east_u(i,jm) = ( east_e(i,jm)*3.- east_e(i,jm-1))/2.
        north_u(i,jm) = (north_e(i,jm)*3.-north_e(i,jm-1))/2.
      end do
!
!     v-points:
!
      do j=1,jm
        do i=1,im-1
           east_v(i,j) = ( east_e(i,j)+ east_e(i+1,j))/2.
          north_v(i,j) = (north_e(i,j)+north_e(i+1,j))/2.
        end do
      end do
!
!     Extrapolate ends:
!
      do j=1,jm
         east_v(im,j) = ( east_e(im,j)*3.- east_e(im-1,j))/2.
        north_v(im,j) = (north_e(im,j)*3.-north_e(im-1,j))/2.
      end do
!
!     c-points:
!
      do j=1,jm
        do i=1,im-1
           east_c(i,j) = ( east_v(i,j)+ east_v(i+1,j))/2.
          north_c(i,j) = (north_u(i,j)+north_u(i+1,j))/2.
        end do
      end do
!
!     Extrapolate ends:
!
      do j=1,jm
         east_c(im,j) = ( east_v(im,j)*3.- east_v(im-1,j))/2.
        north_c(im,j) = (north_u(im,j)*3.-north_u(im-1,j))/2.
      end do

      write(*,*) my_task,"::::",minval(h),maxval(h)
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
      call exchange2d_mpi(dum,im_local,jm_local)
      call exchange2d_mpi(dvm,im_local,jm_local)

! close file:
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close: grid',status)

      return
      end
!_______________________________________________________________________
      subroutine read_grid_pnetcdf0  ! fhx:called by subroutine "interp_init"
! read original grid data for wind and satellite data
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      character(len=120) netcdf_grid_file
      integer alon_coarse_varid,alat_coarse_varid,mask_coarse_varid
      integer ncid,status
      integer(mpi_offset_kind) start(2),edge(2)

! open netcdf file, /wrk/newshoni/fayumi/work/mpi-profs/netcdf_files/grid/grid_i402j252_hmin12m.nc
      write(netcdf_grid_file,'(a)') "./in/grid_i402j252_hmin12m.nc"

      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_grid_file)
      status=nfmpi_open(pom_comm_coarse,netcdf_grid_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_grid_file
     &                         , status )

! get variables

      status=nfmpi_inq_varid(ncid,'east_e',alon_coarse_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: east_e',  status )

      status=nfmpi_inq_varid(ncid,'north_e',alat_coarse_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: north_e', status )

      status=nfmpi_inq_varid(ncid,'fsm',mask_coarse_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: fsm',     status )

! get data


      start(1)=i_global_coarse(1)
      start(2)=j_global_coarse(1)
      edge(1)=im_coarse
      edge(2)=jm_coarse

      status=nfmpi_get_vara_real_all(ncid,alon_coarse_varid,start,edge,
     $                               alon_coarse)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all', status )

      status=nfmpi_get_vara_real_all(ncid,alat_coarse_varid,start,edge,
     $                               alat_coarse)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all', status )

      status=nfmpi_get_vara_real_all(ncid,mask_coarse_varid,start,edge,
     $                               mask_coarse)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all', status )


! close file:
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf( 'nf_close: grid', status )

      return
      end

!_______________________________________________________________________
      subroutine read_restart_pnetcdf
! read data for a seamless restart
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer get_var_real_2d, get_var_real_3d
      character(len=120) netcdf_in_file
      integer ncid,status
      integer wubot_varid,wvbot_varid,aam2d_varid,ua_varid,
     $        uab_varid,va_varid,vab_varid,el_varid,elb_varid,et_varid,
     $        etb_varid,egb_varid,utb_varid,vtb_varid,u_varid,ub_varid,
     $        w_varid,v_varid,vb_varid,t_varid,tb_varid,s_varid,
     $        sb_varid,rho_varid,adx2d_varid,ady2d_varid,advua_varid,
     $        advva_varid,km_varid,kh_varid,kq_varid,l_varid,q2_varid,
     $        q2b_varid,aam_varid,q2l_varid,q2lb_varid
! ayumi 2010/4/8
      integer wusurf_varid, wvsurf_varid
      integer(mpi_offset_kind) start(4),edge(4)

      start = 1
      edge  = 1

! open netcdf restart file
      write(netcdf_in_file,'(''in/'',a)') trim(read_rst_file)
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $ trim(netcdf_in_file)
      status=nfmpi_open(pom_comm,netcdf_in_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//trim(read_rst_file)
     &                         , status )

! get variables
! ayumi 2010/4/8 --------------------------------------------
!      status=nfmpi_inq_varid(ncid,'time',time_varid)
!      call handle_error_pnetcdf('nfmpi_inq_varid: time'
!     $                          ,status)
      status=nfmpi_inq_varid(ncid,'wusurf',wusurf_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: wusurf', status )
      status=nfmpi_inq_varid(ncid,'wvsurf',wvsurf_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: wvsurf', status )
!------------------------------------------------------------
      status=nfmpi_inq_varid(ncid,'wubot',wubot_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: wubot', status )
      status=nfmpi_inq_varid(ncid,'wvbot',wvbot_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: wvbot', status )
      status=nfmpi_inq_varid(ncid,'aam2d',aam2d_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: aam2d', status )
      status=nfmpi_inq_varid(ncid,'ua',ua_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: ua',  status )
      status=nfmpi_inq_varid(ncid,'uab',uab_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: uab', status )
      status=nfmpi_inq_varid(ncid,'va',va_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: va',  status )
      status=nfmpi_inq_varid(ncid,'vab',vab_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: vab', status )
      status=nfmpi_inq_varid(ncid,'el',el_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: el',  status )
      status=nfmpi_inq_varid(ncid,'elb',elb_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: elb', status )
      status=nfmpi_inq_varid(ncid,'et',et_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: et',  status )
      status=nfmpi_inq_varid(ncid,'etb',etb_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: etb', status )
      status=nfmpi_inq_varid(ncid,'egb',egb_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: egb', status )
      status=nfmpi_inq_varid(ncid,'utb',utb_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: utb', status )
      status=nfmpi_inq_varid(ncid,'vtb',vtb_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: vtb', status )
      status=nfmpi_inq_varid(ncid,'u',u_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: u',  status )
      status=nfmpi_inq_varid(ncid,'ub',ub_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: ub', status )
      status=nfmpi_inq_varid(ncid,'v',v_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: v',  status )
      status=nfmpi_inq_varid(ncid,'vb',vb_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: vb', status )
      status=nfmpi_inq_varid(ncid,'w',w_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: w',  status )
      status=nfmpi_inq_varid(ncid,'t',t_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: t',  status )
      status=nfmpi_inq_varid(ncid,'tb',tb_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: tb', status )
      status=nfmpi_inq_varid(ncid,'s',s_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: s',  status )
      status=nfmpi_inq_varid(ncid,'sb',sb_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: sb', status )
      status=nfmpi_inq_varid(ncid,'rho',rho_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: rho', status )
      status=nfmpi_inq_varid(ncid,'adx2d',adx2d_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: adx2d', status )
      status=nfmpi_inq_varid(ncid,'ady2d',ady2d_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: ady2d', status )
      status=nfmpi_inq_varid(ncid,'advua',advua_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: advua', status )
      status=nfmpi_inq_varid(ncid,'advva',advva_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: advva', status )
      status=nfmpi_inq_varid(ncid,'km',km_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: km', status )
      status=nfmpi_inq_varid(ncid,'kh',kh_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: kh', status )
      status=nfmpi_inq_varid(ncid,'kq',kq_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: kq', status )
      status=nfmpi_inq_varid(ncid,'l',l_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: l',  status )
      status=nfmpi_inq_varid(ncid,'q2',q2_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: q2',  status )
      status=nfmpi_inq_varid(ncid,'q2b',q2b_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: q2b', status )
      status=nfmpi_inq_varid(ncid,'aam',aam_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: aam', status )
      status=nfmpi_inq_varid(ncid,'q2l',q2l_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: q2l',  status )
      status=nfmpi_inq_varid(ncid,'q2lb',q2lb_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: q2lb', status )

! get data
! ayumi 2010/4/8
!      start(1)=1
!      edge(1)=1
!      status=nfmpi_get_vara_real_all(ncid,time_varid,start,edge,time0)
!      call handle_error_pnetcdf('nfmpi_get_vara_real_all',
!     $                          status,nf_noerr)

      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
! ayumi 2010/4/8 --------------------------------------------------------
      status = get_var_real_2d(ncid,wusurf_varid,start,edge,wusurf)
      call handle_error_pnetcdf('get_var_real: wusurf',status)
      status = get_var_real_2d(ncid,wvsurf_varid,start,edge,wvsurf)
      call handle_error_pnetcdf('get_var_real: wvsurf',status)
!------------------------------------------------------------------------
      status = get_var_real_2d(ncid,wubot_varid,start,edge,wubot)
      call handle_error_pnetcdf('get_var_real: wubot',status)
      status = get_var_real_2d(ncid,wvbot_varid,start,edge,wvbot)
      call handle_error_pnetcdf('get_var_real: wvbot',status)
      status = get_var_real_2d(ncid,aam2d_varid,start,edge,aam2d)
      call handle_error_pnetcdf('get_var_real: aam2d',status)
      status = get_var_real_2d(ncid,ua_varid,start,edge,ua)
      call handle_error_pnetcdf('get_var_real: ua',status)
      status = get_var_real_2d(ncid,uab_varid,start,edge,uab)
      call handle_error_pnetcdf('get_var_real: uab',status)
      status = get_var_real_2d(ncid,va_varid,start,edge,va)
      call handle_error_pnetcdf('get_var_real: va',status)
      status = get_var_real_2d(ncid,vab_varid,start,edge,vab)
      call handle_error_pnetcdf('get_var_real: vab',status)
      status = get_var_real_2d(ncid,el_varid,start,edge,el)
      call handle_error_pnetcdf('get_var_real: el',status)
      status = get_var_real_2d(ncid,elb_varid,start,edge,elb)
      call handle_error_pnetcdf('get_var_real: elb',status)
      status = get_var_real_2d(ncid,et_varid,start,edge,et)
      call handle_error_pnetcdf('get_var_real: et',status)
      status = get_var_real_2d(ncid,etb_varid,start,edge,etb)
      call handle_error_pnetcdf('get_var_real: etb',status)
      status = get_var_real_2d(ncid,egb_varid,start,edge,egb)
      call handle_error_pnetcdf('get_var_real: egb',status)
      status = get_var_real_2d(ncid,utb_varid,start,edge,utb)
      call handle_error_pnetcdf('get_var_real: utb',status)
      status = get_var_real_2d(ncid,vtb_varid,start,edge,vtb)
      call handle_error_pnetcdf('get_var_real: vtb',status)
      status = get_var_real_2d(ncid,adx2d_varid,start,edge,adx2d)
      call handle_error_pnetcdf('get_var_real: adx2d',status)
      status = get_var_real_2d(ncid,ady2d_varid,start,edge,ady2d)
      call handle_error_pnetcdf('get_var_real: ady2d',status)
      status = get_var_real_2d(ncid,advua_varid,start,edge,advua)
      call handle_error_pnetcdf('get_var_real: advua',status)
      status = get_var_real_2d(ncid,advva_varid,start,edge,advva)
      call handle_error_pnetcdf('get_var_real: advva',status)

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=kb
      status = get_var_real_3d(ncid,u_varid,start,edge,u)
      call handle_error_pnetcdf('get_var_real: u',status)
      status = get_var_real_3d(ncid,ub_varid,start,edge,ub)
      call handle_error_pnetcdf('get_var_real: ub',status)
      status = get_var_real_3d(ncid,v_varid,start,edge,v)
      call handle_error_pnetcdf('get_var_real: v',status)
      status = get_var_real_3d(ncid,vb_varid,start,edge,vb)
      call handle_error_pnetcdf('get_var_real: vb',status)
      status = get_var_real_3d(ncid,w_varid,start,edge,w)
      call handle_error_pnetcdf('get_var_real: w',status)
      status = get_var_real_3d(ncid,t_varid,start,edge,t)
      call handle_error_pnetcdf('get_var_real: t',status)
      status = get_var_real_3d(ncid,tb_varid,start,edge,tb)
      call handle_error_pnetcdf('get_var_real: tb',status)
      status = get_var_real_3d(ncid,s_varid,start,edge,s)
      call handle_error_pnetcdf('get_var_real: s',status)
      status = get_var_real_3d(ncid,sb_varid,start,edge,sb)
      call handle_error_pnetcdf('get_var_real: sb',status)
      status = get_var_real_3d(ncid,rho_varid,start,edge,rho)
      call handle_error_pnetcdf('get_var_real: rho',status)
      status = get_var_real_3d(ncid,km_varid,start,edge,km)
      call handle_error_pnetcdf('get_var_real: km',status)
      status = get_var_real_3d(ncid,kh_varid,start,edge,kh)
      call handle_error_pnetcdf('get_var_real: kh',status)
      status = get_var_real_3d(ncid,kq_varid,start,edge,kq)
      call handle_error_pnetcdf('get_var_real: kq',status)
      status = get_var_real_3d(ncid,l_varid,start,edge,l)
      call handle_error_pnetcdf('get_var_real: l',status)
      status = get_var_real_3d(ncid,q2_varid,start,edge,q2)
      call handle_error_pnetcdf('get_var_real: q2',status)
      status = get_var_real_3d(ncid,q2b_varid,start,edge,q2b)
      call handle_error_pnetcdf('get_var_real: q2b',status)
      status = get_var_real_3d(ncid,aam_varid,start,edge,aam)
      call handle_error_pnetcdf('get_var_real: aam',status)
      status = get_var_real_3d(ncid,q2l_varid,start,edge,q2l)
      call handle_error_pnetcdf('get_var_real: q2l',status)
      status = get_var_real_3d(ncid,q2lb_varid,start,edge,q2lb)
      call handle_error_pnetcdf('get_var_real: q2bl',status)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

! initialize elevation
      d  = h+el
      dt = h+et

! initialize time
      time=time0

      return
      end

!_______________________________________________________________________
      subroutine read_initial_ts_pnetcdf(temp,salt,ssh,n)
! read initial temperature and salinity
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer get_var_real_2d, get_var_real_3d
      integer k, n
      real(kind=rk) temp(im,jm,kb),salt(im,jm,kb),ssh(im,jm)
      character(len=120) netcdf_ic_file
      integer ncid,status
      integer tb_varid,sb_varid,el_varid
      integer(mpi_offset_kind) start(4),edge(4)

      start = 1
      edge  = 1

! open netcdf ic file
      write(netcdf_ic_file,'(a)') "./in/tsclim/ts_clim.nc"

      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_ic_file)
      status=nfmpi_open(pom_comm,netcdf_ic_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_ic_file
     $                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'tclim',tb_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: t_init'
     $                         , status )
      status=nfmpi_inq_varid(ncid,'sclim',sb_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: s_init'
     $                         , status )
      status=nfmpi_inq_varid(ncid,'eclim',el_varid)

      if ( status == nf_noerr ) then
! get data
        start(1)=i_global(1)
        start(2)=j_global(1)
        start(3)=n
        edge(1)=im
        edge(2)=jm
        edge(3)=1
        status = get_var_real_2d(ncid,el_varid,start,edge,ssh)
        call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     &                           , status )
      else
        if ( status == -49 ) then
          if ( my_task == master_task ) then
            print *, "missing elevation data - falling back to 0"
          end if
          ssh = 0.
        else
          call handle_error_pnetcdf( 'nfmpi_inq_varid: el_init'
     &                             , status )
        end if
      end if

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      start(4)=n
      edge(1)=im
      edge(2)=jm
      edge(3)=kb
      edge(4)=1
      status = get_var_real_3d(ncid,tb_varid,start,edge,temp)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     &                         , status )
      status = get_var_real_3d(ncid,sb_varid,start,edge,salt)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     &                         , status )
      do k=1,kb-1
        where (fsm==0.)
          temp(:,:,k) = 0.
          salt(:,:,k) = 0.
        end where
      end do
      temp(:,:,kb) = temp(:,:,kb-1)
      salt(:,:,kb) = salt(:,:,kb-1)
      elb = elb*fsm

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end
!_______________________________________________________________________
      subroutine read_clim_ts_pnetcdf(temp,salt,n)
! read annual-mean temperature and salinity
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer get_var_real_3d
      integer k,n
!      real(kind=rk) temp(im,jm,k),salt(im,jm,k)
      real(kind=rk) temp(im,jm,kb),salt(im,jm,kb)
      character(len=120) netcdf_ic_file
      integer ncid,status
      integer tclim_varid,sclim_varid
      integer(mpi_offset_kind) start(4),edge(4)

! open netcdf ic file
      write(netcdf_ic_file,'(a)') "./in/tsclim/ts_clim.nc"

      if(my_task.eq.master_task) write(*,'(/''reading file '',a)')
     $  trim(netcdf_ic_file)
      status=nfmpi_open(pom_comm,netcdf_ic_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_ic_file
     &                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'tclim',tclim_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: tclim'
     &                         , status )
      status=nfmpi_inq_varid(ncid,'sclim',sclim_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: sclim'
     &                         , status )

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      start(4)=n
      edge(1)=im
      edge(2)=jm
      edge(3)=kb
      edge(4)=1
      status = get_var_real_3d(ncid,tclim_varid,start,edge,temp)
      call handle_error_pnetcdf('get_var_real: tclim',status)
      status = get_var_real_3d(ncid,sclim_varid,start,edge,salt)
      call handle_error_pnetcdf('get_var_real: sclim',status)

      do k=1,kbm1
        where (fsm==0.)
          temp(:,:,k) = 0.
          salt(:,:,k) = 0.
        end where
      end do
      temp(:,:,kb) = temp(:,:,kbm1)
      salt(:,:,kb) = salt(:,:,kbm1)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end
!_______________________________________________________________________
      subroutine read_clim_ts_pnetcdf_obs(temp,salt,rm,n)
! read annual-mean temperature and salinity
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer get_var_real_3d
      integer k,n
!      real(kind=rk) temp(im,jm,k),salt(im,jm,k)
      real(kind=rk) temp(im,jm,kb),salt(im,jm,kb),rm(im,jm,kb)
      character(len=120) netcdf_ic_file
      integer ncid,status
      integer tclim_varid,sclim_varid,rm_varid
      integer(mpi_offset_kind) start(4),edge(4)

! open netcdf ic file
      write(netcdf_ic_file,'(a)') "./in/tsclim/ts_clim_old.nc"

      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_ic_file)
      status=nfmpi_open(pom_comm,netcdf_ic_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_ic_file
     &                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'tclim',tclim_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: tclim'
     &                         , status )
      status=nfmpi_inq_varid(ncid,'sclim',sclim_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: sclim'
     &                         , status )
      status=nfmpi_inq_varid(ncid,'Rmean',rm_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: rmean'
     &                         , status )

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      start(4)=n
      edge(1)=im
      edge(2)=jm
      edge(3)=kb
      edge(4)=1
      status = get_var_real_3d(ncid,tclim_varid,start,edge,temp)
      call handle_error_pnetcdf('get_var_real: tclim',status)
      status = get_var_real_3d(ncid,sclim_varid,start,edge,salt)
      call handle_error_pnetcdf('get_var_real: sclim',status)
      status = get_var_real_3d(ncid,rm_varid,start,edge,rm)
      call handle_error_pnetcdf('get_var_real: rmean',status)

      do k=1,kb
        where (fsm==0.)
          temp(:,:,k) = 0.
          salt(:,:,k) = 0.
            rm(:,:,k) = 0.
        end where
      end do

        rm(:,:,kb) =   rm(:,:,kbm1)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end

!_______________________________________________________________________
      subroutine read_mean_ts_pnetcdf(temp,salt)  !lyo:20110202
! read xy-averaged, annual-mean temp & salt

      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer get_var_real_3d
!      real(kind=rk) temp(im,jm,k),salt(im,jm,k)
      real(kind=rk) temp(im,jm,kb),salt(im,jm,kb)
      character(len=120) netcdf_ic_file
      integer ncid,status
      integer tmean_varid,smean_varid
      integer(mpi_offset_kind) start(4),edge(4)

      start = 1
      edge  = 1
!
! open netcdf ic file !lyo:try first *mean* - the correct one
      write(netcdf_ic_file,'(a)') "./in/tsclim/ts_mean.nc"

      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_ic_file)
      status=nfmpi_open(pom_comm,netcdf_ic_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_ic_file
     &                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'tmean',tmean_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: t'
     &                         , status )
      status=nfmpi_inq_varid(ncid,'smean',smean_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: s'
     &                         , status )

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=kb
      status = get_var_real_3d(ncid,tmean_varid,start,edge,temp)
      call handle_error_pnetcdf('get_var_real: tmean',status)
      status = get_var_real_3d(ncid,smean_varid,start,edge,salt)
      call handle_error_pnetcdf('get_var_real: smean',status)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end

!_______________________________________________________________________
      subroutine read_mean_ts_z_pnetcdf(temp,salt,ks,n)
! read xy-averaged, annual-mean temp & salt

      use mpi
      use pnetcdf

      implicit none

      include 'pom.h'
      integer get_var_real_1d, get_var_real_2d, ks
      real(kind=rk) t1(kb),s1(kb),zs(ks)
!      real(kind=rk) tz(im,jm,ks),sz(im,jm,ks)
      real(kind=rk) temp(im,jm,kb),salt(im,jm,kb)
      character(len=120) netcdf_ic_file
      integer ncid,status, i,j, n
      integer tmean_varid,smean_varid,z_varid
      integer(mpi_offset_kind) start(4),edge(4)

      start = 1
      edge  = 1
!
! open netcdf ic file !lyo:try first *mean* - the correct one
      write(netcdf_ic_file,'(a)') "./in/tsclim/ts_mean_z.nc"

      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_ic_file)
      status=nfmpi_open(pom_comm,netcdf_ic_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_ic_file
     &                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'z',z_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: z'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'tmean_s',tmean_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: t'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'smean_s',smean_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: s'
     $                          ,status)

! get data
      start(1)=1
      edge(1)=ks
      status = get_var_real_1d(ncid,z_varid,start,edge,zs)
      call handle_error_pnetcdf('get_var_real: zs',status)
      start(1)=1
      start(2)=n
      edge(1)=kb
      edge(2)=1
      print *, my_task, "::", n, ks
      status = get_var_real_2d(ncid,tmean_varid,start,edge,t1)
      call handle_error_pnetcdf('get_var_real: tmean',status)
      status = get_var_real_2d(ncid,smean_varid,start,edge,s1)
      call handle_error_pnetcdf('get_var_real: smean',status)
      print *, my_task, "::", t1
      do j=1,jm
        do i=1,im
!          tz(i,j,:) = t1
!          sz(i,j,:) = s1
          temp(i,j,:) = t1
          salt(i,j,:) = s1
        end do
      end do

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

!      call ztosig(zs,tz,zz,h,temp,im,jm,ks,kb,
!     $                                  n_west,n_east,n_south,n_north)
!      call ztosig(zs,sz,zz,h,salt,im,jm,ks,kb,
!     $                                  n_west,n_east,n_south,n_north)

      return
      end

!_______________________________________________________________________
      subroutine read_boundary_conditions_pnetcdf !lyo:pac10:not used!
! read boundary conditions
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      character(len=120) netcdf_bc_file
      integer ncid,status
      integer tw_varid,sw_varid,uw_varid,ew_varid
      integer te_varid,se_varid,ue_varid,ee_varid
      integer ts_varid,ss_varid,vs_varid,es_varid
      integer tn_varid,sn_varid,vn_varid,en_varid
      integer(mpi_offset_kind) start(2),edge(2)

! open netcdf ic file
      write(netcdf_bc_file,'(''in/'',a,''.bc.nc'')')
     $  trim(netcdf_file)

      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_bc_file)
      status=nfmpi_open(pom_comm,netcdf_bc_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_bc_file
     &                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'tw',tw_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: tw',status)
      status=nfmpi_inq_varid(ncid,'sw',sw_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: sw',status)
      status=nfmpi_inq_varid(ncid,'uw',uw_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: uw',status)
      status=nfmpi_inq_varid(ncid,'ew',ew_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: ew',status)
      status=nfmpi_inq_varid(ncid,'te',te_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: te',status)
      status=nfmpi_inq_varid(ncid,'se',se_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: se',status)
      status=nfmpi_inq_varid(ncid,'ue',ue_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: ue',status)
      status=nfmpi_inq_varid(ncid,'ee',ee_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: ee',status)
      status=nfmpi_inq_varid(ncid,'ts',ts_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: ts',status)
      status=nfmpi_inq_varid(ncid,'ss',ss_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: ss',status)
      status=nfmpi_inq_varid(ncid,'vs',vs_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: vs',status)
      status=nfmpi_inq_varid(ncid,'es',es_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: es',status)
      status=nfmpi_inq_varid(ncid,'tn',tn_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: tn',status)
      status=nfmpi_inq_varid(ncid,'sn',sn_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: sn',status)
      status=nfmpi_inq_varid(ncid,'vn',vn_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: vn',status)
      status=nfmpi_inq_varid(ncid,'en',en_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: en',status)

! get data
      start(1)=j_global(1)
      edge(1)=jm
      status=nfmpi_get_vara_real_all(ncid,uw_varid,start,edge,uabw)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: uabw'
     $                         , status )
      status=nfmpi_get_vara_real_all(ncid,ew_varid,start,edge,elw)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: elw'
     $                         , status )
      status=nfmpi_get_vara_real_all(ncid,ue_varid,start,edge,uabe)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: uabe'
     $                         , status )
      status=nfmpi_get_vara_real_all(ncid,ee_varid,start,edge,ele)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: ele'
     $                         , status )

      start(1)=i_global(1)
      edge(1)=im
      status=nfmpi_get_vara_real_all(ncid,vs_varid,start,edge,vabs)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: vabs'
     $                         , status )
      status=nfmpi_get_vara_real_all(ncid,es_varid,start,edge,els)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: els'
     $                         , status )
      status=nfmpi_get_vara_real_all(ncid,vn_varid,start,edge,vabn)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: vabn'
     $                         , status )
      status=nfmpi_get_vara_real_all(ncid,en_varid,start,edge,eln)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: eln'
     $                         , status )

      start(1)=j_global(1)
      start(2)=1
      edge(1)=jm
      edge(2)=kb
      status=nfmpi_get_vara_real_all(ncid,tw_varid,start,edge,tbw)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: tbw'
     $                         , status )
      status=nfmpi_get_vara_real_all(ncid,sw_varid,start,edge,sbw)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: sbw'
     $                         , status )
      status=nfmpi_get_vara_real_all(ncid,te_varid,start,edge,tbe)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: tbe'
     $                         , status )
      status=nfmpi_get_vara_real_all(ncid,se_varid,start,edge,sbe)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: sbe'
     $                         , status )

      start(1)=i_global(1)
      start(2)=1
      edge(1)=im
      edge(2)=kb
      status=nfmpi_get_vara_real_all(ncid,ts_varid,start,edge,tbs)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: tbs'
     $                         , status )
      status=nfmpi_get_vara_real_all(ncid,ss_varid,start,edge,sbs)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: sbs'
     $                         , status )
      status=nfmpi_get_vara_real_all(ncid,tn_varid,start,edge,tbn)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: tbn'
     $                         , status )
      status=nfmpi_get_vara_real_all(ncid,sn_varid,start,edge,sbn)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all: sbn'
     $                         , status )

! close file
      status = nfmpi_close(ncid)
      call handle_error_pnetcdf( 'nf_close', status )

      return
      end
!_______________________________________________________________________
      subroutine read_wind_pnetcdf0(n,wu,wv) !lyo:pac10:not used!
! read wind stress
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer n
      real(kind=rk) wu(im,jm),wv(im,jm)
      character(len=120) netcdf_wind_file
      integer ncid,status
      integer wu_varid,wv_varid
      integer(mpi_offset_kind) start(2),edge(2)

! open netcdf wind file
      write(netcdf_wind_file,'(''in/'',a,''.wind.'',i4.4,''.nc'')')
     $                                               trim(netcdf_file),n
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $                                            trim(netcdf_wind_file)
      status=nfmpi_open(pom_comm,netcdf_wind_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_wind_file
     $                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'wu',wu_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: wu'
     $                         , status )
      status=nfmpi_inq_varid(ncid,'wv',wv_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: wv'
     $                         , status )

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      status=nfmpi_get_vara_real_all(ncid,wu_varid,start,edge,wu)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     $                         , status )
      status=nfmpi_get_vara_real_all(ncid,wv_varid,start,edge,wv)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     $                         , status )

! close file
      status = nfmpi_close(ncid)
      call handle_error_pnetcdf( 'nf_close', status )

      return
      end
!_______________________________________________________________________




! followings are for MAB model
! ayumi since 2010/4/?
!_______________________________________________________________________
      subroutine read_tsclim_monthly_pnetcdf(tm,sm,tc,sc,ht,in_file,n)
! read initial temperature and salinity
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer get_var_real_3d, get_var_real_2d
      real(kind=rk) tm(im,jm,kb),sm(im,jm,kb),tc(im,jm,kb),sc(im,jm,kb)
      real(kind=rk) ht(im,jm)
      character(len=120) netcdf_ic_file
      character(len=*) in_file
      integer ncid,status,n,k
      integer tm_varid,sm_varid,tc_varid,sc_varid,ht_varid
      integer(mpi_offset_kind) start(4),edge(4)

! open netcdf ic file
      write(netcdf_ic_file,'(''in/tsclim/'',a)')
     $     trim(in_file)
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_ic_file)
      status=nfmpi_open(pom_comm,netcdf_ic_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      if (status/=NF_NOERR) then
          if (my_task==0) write(*,*) 'Failed reading clim! Fallback...'
          tm = 15.
          sm = 33.
          tc = 15.
          sc = 33.
          return
      end if
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_ic_file
     &                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'tmean',tm_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: tm'
     $                         , status )
      status=nfmpi_inq_varid(ncid,'smean',sm_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: sm'
     $                         , status )
      status=nfmpi_inq_varid(ncid,'tclim',tc_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: tc'
     $                         , status )
      status=nfmpi_inq_varid(ncid,'sclim',sc_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: sc'
     $                         , status )
      status=nfmpi_inq_varid(ncid,'eclim',ht_varid)
      if ( status /= nf_noerr ) then
        if ( my_task == master_task ) then
          print *, "Failed reading elevation"
        end if
      end if

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      start(4)=n
      edge(1)=im
      edge(2)=jm
      edge(3)=kb
      edge(4)=1
      status = get_var_real_3d(ncid,tm_varid,start,edge,tm)
      call handle_error_pnetcdf('get_var_real: tm',status)
      status = get_var_real_3d(ncid,sm_varid,start,edge,sm)
      call handle_error_pnetcdf('get_var_real: sm',status)
      status = get_var_real_3d(ncid,tc_varid,start,edge,tc)
      call handle_error_pnetcdf('get_var_real: tc',status)
      status = get_var_real_3d(ncid,sc_varid,start,edge,sc)
      call handle_error_pnetcdf('get_var_real: sc',status)
      start(3)=n
      edge(3)=1
      status = get_var_real_2d(ncid,ht_varid,start,edge,ht)
      if ( status /= nf_noerr ) then
        if ( my_task == master_task ) then
          print *, "Falling back to zero"
        end if
        el = 0.
      end if

      do k=1,kbm1
        where (fsm(:,:)==0.)
          tm(:,:,k) = 0.
          sm(:,:,k) = 0.
          tc(:,:,k) = 0.
          sc(:,:,k) = 0.
          el        = 0.
        end where
      end do

!      tm(:,:,kb) = tm(:,:,kb-1)
!      sm(:,:,kb) = sm(:,:,kb-1)
!      tc(:,:,kb) = tc(:,:,kb-1)
!      sc(:,:,kb) = sc(:,:,kb-1)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end

!_______________________________________________________________________
      subroutine read_tsclim_monthly_pnetcdf_obs(rm,tc,sc,in_file,n)
! read initial temperature and salinity
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer get_var_real_3d
      real(kind=rk) rm(im,jm,kb),tc(im,jm,kb),sc(im,jm,kb)
      character(len=120) netcdf_ic_file
      character(len=*) in_file
      integer ncid,status,n,k
      integer rm_varid,tc_varid,sc_varid
      integer(mpi_offset_kind) start(4),edge(4)

! open netcdf ic file
      write(netcdf_ic_file,'(''in/tsclim/'',a)')
     $     trim(in_file)
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_ic_file)
      status=nfmpi_open(pom_comm,netcdf_ic_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_ic_file
     $                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'Rmean',rm_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: rm', status )
      status=nfmpi_inq_varid(ncid,'tclim',tc_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: tc', status )
      status=nfmpi_inq_varid(ncid,'sclim',sc_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: sc', status )

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      start(4)=n
      edge(1)=im
      edge(2)=jm
      edge(3)=kb
      edge(4)=1
      status = get_var_real_3d(ncid,rm_varid,start,edge,rm)
      call handle_error_pnetcdf('get_var_real: rm',status)
      status = get_var_real_3d(ncid,tc_varid,start,edge,tc)
      call handle_error_pnetcdf('get_var_real: tc',status)
      status = get_var_real_3d(ncid,sc_varid,start,edge,sc)
      call handle_error_pnetcdf('get_var_real: sc',status)

      do k=1,kbm1
        where (fsm(:,:)==0.)
          rm(:,:,k) = 0.
          tc(:,:,k) = 0.
          sc(:,:,k) = 0.
        end where
      end do

      rm(:,:,kb) = rm(:,:,kb-1)
      tc(:,:,kb) = tc(:,:,kb-1)
      sc(:,:,kb) = sc(:,:,kb-1)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end

!_______________________________________________________________________
      subroutine read_wind_monthly_pnetcdf(uw,vw,in_file,n)
! read initial temperature and salinity
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer get_var_real_2d
      real(kind=rk) uw(im,jm),vw(im,jm)
      character(len=120) netcdf_ic_file
      character(len=*) in_file
      integer ncid,status,n
      integer uw_varid,vw_varid
      integer(mpi_offset_kind) start(4),edge(4)

! open netcdf ic file
      write(netcdf_ic_file,'(''in/tsclim/'',a)')
     $     trim(in_file)
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_ic_file)
      status=nfmpi_open(pom_comm,netcdf_ic_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_ic_file
     $                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'wU',uw_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: uw'
     &                         , status )
      status=nfmpi_inq_varid(ncid,'wV',vw_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: vw'
     &                         , status )

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=n
      edge(1)=im
      edge(2)=jm
      edge(3)=1
      status = get_var_real_2d(ncid,uw_varid,start,edge,uw)
      call handle_error_pnetcdf('get_var_real: uw',status)
      status = get_var_real_2d(ncid,vw_varid,start,edge,vw)
      call handle_error_pnetcdf('get_var_real: vw',status)

      where (fsm==0.)
        uw = 0.
        vw = 0.
      end where

! close file
      status = nfmpi_close(ncid)
      call handle_error_pnetcdf( 'nf_close', status )

      return
      end

!_______________________________________________________________________
      subroutine read_ssts_monthly_pnetcdf(sst,sss,in_file)
! read wind stress
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      real(kind=rk) sst(im,jm),sss(im,jm)
      character(len=9) in_file
      character(len=120) netcdf_ssts_file
      integer ncid,status
      integer sst_varid,sss_varid
      integer(mpi_offset_kind) start(2),edge(2)

! open netcdf wind file
      write(netcdf_ssts_file,'(''in/ssts/'',a)') trim(in_file)
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $                                            trim(netcdf_ssts_file)
      status=nfmpi_open(pom_comm,netcdf_ssts_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_ssts_file
     &                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'sst',sst_varid)
!     status=nfmpi_inq_varid(ncid,'SST',sst_varid) !lyo:pac10:AluUpperSST
      call handle_error_pnetcdf( 'nfmpi_inq_varid: sst', status )
      status=nfmpi_inq_varid(ncid,'sss',sss_varid)
!     status=nfmpi_inq_varid(ncid,'SSS',sss_varid) !lyo:pac10:AluUpperSSS
      call handle_error_pnetcdf( 'nfmpi_inq_varid: sss', status )

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      status=nfmpi_get_vara_real_all(ncid,sst_varid,start,edge,sst)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     $                         , status )
      status=nfmpi_get_vara_real_all(ncid,sss_varid,start,edge,sss)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     $                         , status )

! close file
      status = nfmpi_close(ncid)
      call handle_error_pnetcdf( 'nf_close', status )

      return
      end
!_______________________________________________________________________

!      subroutine read_uabe_pnetcdf(uabe_out) !lyo:pac10:use stcc version
      subroutine read_bc_pnetcdf(ube_out, ubw_out,
     &                           vbs_out, vbn_out, in_file, n) !eda: read all boundaries

      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      character(len=64) v_uabe,v_uabw,v_vabs,v_vabn
      namelist/ncinf_nml/ v_uabe,v_uabw,v_vabs,v_vabn
      integer get_var_real_3d
      character(len=120) netcdf_bc_file
      integer ncid,status
      integer uabe_varid, uabw_varid, vabs_varid, vabn_varid
      integer(mpi_offset_kind) start(4),edge(4)
      real(kind=rk) ube_out(jm,kb), ubw_out(jm,kb)
     $   , vbs_out(im,kb), vbn_out(im,kb)
      logical :: fexist !lyo:20110224:alu:stcc:
      character(len=*), intent(in) :: in_file
      integer, intent(in) :: n

      open(73,file='ncinf.nml',status='old')
      read(73,nml=ncinf_nml)
      close(73)

      start = 1
      edge  = 1
!
! eda:uvforce
      write(netcdf_bc_file,'(''in/bc/'',a)')
     $     trim(in_file)

!      write(netcdf_bc_file,'(a)') "./in/bc.nc"

      inquire(file=trim(netcdf_bc_file),exist=fexist)
!
      if(fexist) then
!
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_bc_file)
      status=nfmpi_open(pom_comm,netcdf_bc_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_bc_file
     $                         , status )
! get variables
      status=nfmpi_inq_varid(ncid,trim(v_uabe),uabe_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: uabe',status)
      status=nfmpi_inq_varid(ncid,trim(v_uabw),uabw_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: uabw',status)
      status=nfmpi_inq_varid(ncid,trim(v_vabs),vabs_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: vabs',status)
      status=nfmpi_inq_varid(ncid,trim(v_vabn),vabn_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: vabn',status)
! get data
      start(1)=j_global(1)
      start(2)=1
      start(3)=n
      edge(1)=jm
      edge(2)=kb
      edge(3)=1
      status = get_var_real_3d(ncid,uabe_varid,start,edge,ube_out)
      call handle_error_pnetcdf('get_var_real_3d: uabe',status)
      status = get_var_real_3d(ncid,uabw_varid,start,edge,ubw_out)
      call handle_error_pnetcdf('get_var_real_3d: uabw',status)

      start(1)=i_global(1)
      start(2)=1
      start(3)=n
      edge(1)=im
      edge(2)=kb
      edge(3)=1
      status = get_var_real_3d(ncid,vabs_varid,start,edge,vbs_out)
      call handle_error_pnetcdf('get_var_real_3d: vabs',status)
      status = get_var_real_3d(ncid,vabn_varid,start,edge,vbn_out)
      call handle_error_pnetcdf('get_var_real_3d: vabn',status)

! close file
      status = nfmpi_close(ncid)
      call handle_error_pnetcdf( 'nf_close', status )
!
      else !lyo:no, does not exist, give warning...

      if(my_task.eq.master_task) then
!----------------------------------------------------------------------!
         write(6,'(/''****Warning from subr. read_bc_pnetcdf****'')')
         write(6,'(/''bc file '',a)')trim(netcdf_bc_file)
         write(6,'(/''does not exist, and bc set to zero'')')
         ube_out(:,:)=0.0 ; ubw_out(:,:)=0.0
         vbs_out(:,:)=0.0 ; vbn_out(:,:)=0.0
      end if

      endif

      return
      end
!_______________________________________________________________________
      subroutine read_bc_transport_pnetcdf(trans,y,x, in_file, n)

      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer get_var_real_1d
      character(len=120) netcdf_bc_file
      integer ncid,status
      integer trans_varid, y_varid, x_varid, t_dimid
      integer(mpi_offset_kind) start(4),edge(4)
      logical :: fexist
      integer, intent(in) :: n
      integer n_bry
      real(kind=rk), intent(out) :: trans(10), y(2,10), x(2,10)
      character(len=*), intent(in) :: in_file

      start = 1
      edge  = 1
!
! eda:uvforce
      write(netcdf_bc_file,'(''in/bc/'',a)')
     $     trim(in_file)

!      write(netcdf_bc_file,'(a)') "./in/bc.nc"
      inquire(file=trim(netcdf_bc_file),exist=fexist)
!
      if(fexist) then
!
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_bc_file)
      status=nfmpi_open(pom_comm,netcdf_bc_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_bc_file
     $                         , status )

      status = nfmpi_inq_dimid(ncid, 'num', t_dimid)
      call handle_error_pnetcdf('nfmpi_inq_dimid: num',status)
      status = nfmpi_inq_dim(ncid, t_dimid, netcdf_bc_file, n_bry)
      call handle_error_pnetcdf('nfmpi_inq_dim:num',status)
! get variables
      status=nfmpi_inq_varid(ncid,'transport',trans_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid:trns',status)
      status=nfmpi_inq_varid(ncid,'y_bounds',y_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid:ybnd',status)
      status=nfmpi_inq_varid(ncid,'x_bounds',x_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid:xbnd',status)
! get data
      start(1)=1
      start(2)=n
      edge(1)=n_bry
      edge(2)=1
      status = get_var_real_1d(ncid,trans_varid,start,edge
     &                       ,trans(1:n_bry))
      call handle_error_pnetcdf('get_var_real_1:trans',status)
      start(1)=1
      start(2)=1
      edge(1)=n_bry
      edge(2)=1
      status = get_var_real_1d(ncid,y_varid,start,edge,y(1,1:n_bry))
      call handle_error_pnetcdf('get_var_real_1d: y.1',status)
      status = get_var_real_1d(ncid,x_varid,start,edge,x(1,1:n_bry))
      call handle_error_pnetcdf('get_var_real_1d: x.1',status)
      start(1)=1
      start(2)=2
      edge(1)=n_bry
      edge(2)=1
      status = get_var_real_1d(ncid,y_varid,start,edge,y(2,1:n_bry))
      call handle_error_pnetcdf('get_var_real_1d: y.2',status)
      status = get_var_real_1d(ncid,x_varid,start,edge,x(2,1:n_bry))
      call handle_error_pnetcdf('get_var_real_1d: x.2',status)

! convert Sv to m^3/s
      trans = trans * 1.d6

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)
!
      else !lyo:no, does not exist, give warning...

      if(my_task.eq.master_task) then
!----------------------------------------------------------------------!
         write(6,'(/''****Warning from subr. read_bc_transport****'')')
         write(6,'(/''bc file '',a)')trim(netcdf_bc_file)
         write(6,'(/''does not exist, and bc set to zero'')')
         trans = 0.
         y = 1.
         x = 1.
      end if

      endif

      return
      end
!_______________________________________________________________________
      subroutine read_tide_east_pnetcdf
     $           (ampe_out,phae_out,amue_out,phue_out)
!fhx:tide:read tide at the eastern boundary for PROFS
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      character(len=120) netcdf_bc_file
      integer ncid,status
      integer ampe_varid,phae_varid,amue_varid,phue_varid
      integer(mpi_offset_kind) start(2),edge(2)
      real(kind=rk) ampe_out(jm,ntide),phae_out(jm,ntide)
      real(kind=rk) amue_out(jm,ntide),phue_out(jm,ntide)

! open netcdf ic file
      write(netcdf_bc_file,'(a)') "./in/tide.nc"

      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_bc_file)
      status=nfmpi_open(pom_comm,netcdf_bc_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_bc_file
     $                         , status )


! get variables
      status=nfmpi_inq_varid(ncid,'ampe',ampe_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: ampe',status)

      status=nfmpi_inq_varid(ncid,'phae',phae_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: phae',status)

      status=nfmpi_inq_varid(ncid,'amue',amue_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: amue',status)

      status=nfmpi_inq_varid(ncid,'phue',phue_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: phue',status)


! get data
      start(1)=j_global(1)
      start(2)=1
      edge(1)=jm
      edge(2)=ntide
      status=nfmpi_get_vara_real_all(
     $ ncid,ampe_varid,start,edge,ampe_out)
      call handle_error_pnetcdf(
     $ 'nfmpi_get_vara_real_all: ampe',status)

      status=nfmpi_get_vara_real_all(
     $ ncid,phae_varid,start,edge,phae_out)
      call handle_error_pnetcdf(
     $ 'nfmpi_get_vara_real_all: phae',status)

      status=nfmpi_get_vara_real_all(
     $ ncid,amue_varid,start,edge,amue_out)
      call handle_error_pnetcdf(
     $ 'nfmpi_get_vara_real_all: amue',status)

      status=nfmpi_get_vara_real_all(
     $ ncid,phue_varid,start,edge,phue_out)
      call handle_error_pnetcdf(
     $ 'nfmpi_get_vara_real_all: phue',status)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end
!_______________________________________________________________________

!_______________________________________________________________________
      subroutine read_mflx_pnetcdf( uflx, vflx, filename, n ) !lyo:pac10:use stcc version !lyo:pac10:not used!
! read wind stress
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer get_var_real_3d
      integer, parameter :: wind_npd = 4
      character(len=*), intent(in) :: filename
      real(kind=rk),dimension(im,jm,wind_npd),intent(out) ::
     &                                            uflx, vflx
      character :: netcdf_mflx_file*24
      integer ncid,status,n
      integer uflx_varid,vflx_varid
      integer(mpi_offset_kind) start(4),edge(4)

      start = 1
      edge  = 1

! open netcdf wind file
      netcdf_mflx_file="in/"//windf//"/"//trim(filename)
      if(my_task.eq.0)
     &       write(*,'(/''reading wind stress '',a,'' @ '',i4)')
     &                                  trim(netcdf_mflx_file), n
      status=nfmpi_open(pom_comm,netcdf_mflx_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_mflx_file
     $                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'uflx',uflx_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: uflx', status )
      status=nfmpi_inq_varid(ncid,'vflx',vflx_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: vflx', status )

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=n
      edge(1)=im
      edge(2)=jm
      edge(3)=wind_npd

      status = get_var_real_3d(ncid,uflx_varid,start,edge,uflx)
      call handle_error_pnetcdf('get_var_real: uflx',status)
      status = get_var_real_3d(ncid,vflx_varid,start,edge,vflx)
      call handle_error_pnetcdf('get_var_real: vflx',status)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end
!_______________________________________________________________________
!_______________________________________________________________________
      subroutine read_wind_pnetcdfc( uwnd_out, vwnd_out, filename, n )
! read wind stress
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      include 'io.h'

      integer get_var_real_3d!,k
      integer, parameter :: wind_npd = 4
      real(kind=rk) uwnd_out(im_coarse,jm_coarse,wind_npd)
      real(kind=rk) vwnd_out(im_coarse,jm_coarse,wind_npd)
!      character :: filename*15, netcdf_wind_file*23
      character(len=*)   filename
      character(len=256) netcdf_read_file
      integer ncid,status,n
      integer uwnd_varid,vwnd_varid
      integer(mpi_offset_kind) start(4),edge(4)
      logical lexist

      start = 1
      edge  = 1

! open netcdf wind file
      netcdf_read_file = trim(bulk_path)//trim(filename)

      inquire( file=netcdf_read_file, exist=lexist )

      if ( .not.lexist ) then

        if ( my_task == master_task ) then
          print *, "[!] Missing wind data file at `"
     &             , trim(netcdf_read_file), "`"
        end if
        return

      end if

      if ( my_task == master_task ) then
        write( *,'(/''reading wind '',a,'' @ '',i4)' )
     &                                    trim(netcdf_read_file), n
      end if
      status = nfmpi_open( pom_comm_coarse, netcdf_read_file
     &                   , nf_nowrite, mpi_info_null, ncid )
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_read_file
     &                         , status )

! get variables
      status = nfmpi_inq_varid( ncid, trim(uwnd_name), uwnd_varid )
      call handle_error_pnetcdf( 'nfmpi_inq_varid: uwnd'
     &                         , status )
      status = nfmpi_inq_varid( ncid, trim(vwnd_name), vwnd_varid )
      call handle_error_pnetcdf( 'nfmpi_inq_varid: vwnd'
     &                         , status)

! get data
      start(1)=i_global_coarse(1)
      start(2)=j_global_coarse(1)
      start(3)=n
      edge(1)=im_coarse
      edge(2)=jm_coarse
      edge(3)=wind_npd

      status = get_var_real_3d(ncid,uwnd_varid,start,edge,uwnd_out)
      call handle_error_pnetcdf( 'get_var_real: '//trim(uwnd_name)
     &                         , status )

      status = get_var_real_3d(ncid,vwnd_varid,start,edge,vwnd_out)
      call handle_error_pnetcdf( 'get_var_real: '//trim(vwnd_name)
     &                         , status )

!      do k=1,wind_npd
!        uwnd_out(2:imm1,2:jmm1,k) = uwnd_out(2:imm1,2:jmm1,k)*.25
!     &                        *(fsm(2:imm1,3:jm)+fsm(3:im,2:jmm1)
!     &                         +fsm(2:imm1,1:jmm2)+fsm(1:imm2,2:jmm1))
!        vwnd_out(2:imm1,2:jmm1,k) = vwnd_out(2:imm1,2:jmm1,k)*.25
!     &                        *(fsm(2:imm1,3:jm)+fsm(3:im,2:jmm1)
!     &                         +fsm(2:imm1,1:jmm2)+fsm(1:imm2,2:jmm1))
!        uwnd_out(1,2:jmm1,k) = uwnd_out(1,2:jmm1,k)/3.
!     &                        *(fsm(1,1:jmm2)+fsm(2,2:jmm1)+fsm(1,3:jm))
!        vwnd_out(1,2:jmm1,k) = vwnd_out(1,2:jmm1,k)/3.
!     &                        *(fsm(1,1:jmm2)+fsm(2,2:jmm1)+fsm(1,3:jm))
!        uwnd_out(im,2:jmm1,k) = uwnd_out(im,2:jmm1,k)/3.
!     &                   *(fsm(im,1:jmm2)+fsm(imm1,2:jmm1)+fsm(im,3:jm))
!        vwnd_out(im,2:jmm1,k) = vwnd_out(im,2:jmm1,k)/3.
!     &                   *(fsm(im,1:jmm2)+fsm(imm1,2:jmm1)+fsm(im,3:jm))
!        uwnd_out(2:imm1,1,k) = uwnd_out(2:imm1,1,k)/3.
!     &                        *(fsm(1:imm2,1)+fsm(2:imm1,2)+fsm(3:im,1))
!        vwnd_out(2:imm1,1,k) = vwnd_out(2:imm1,1,k)/3.
!     &                        *(fsm(1:imm2,1)+fsm(2:imm1,2)+fsm(3:im,1))
!        uwnd_out(2:imm1,jm,k) = uwnd_out(2:imm1,jm,k)/3.
!     &                   *(fsm(1:imm2,jm)+fsm(2:imm1,jmm1)+fsm(3:im,jm))
!        vwnd_out(2:imm1,jm,k) = vwnd_out(2:imm1,jm,k)/3.
!     &                   *(fsm(1:imm2,jm)+fsm(2:imm1,jmm1)+fsm(3:im,jm))
!        uwnd_out(1,1,k) = uwnd_out(1,1,k)*.5*(fsm(1,2)+fsm(2,1))
!        vwnd_out(1,1,k) = vwnd_out(1,1,k)*.5*(fsm(1,2)+fsm(2,1))
!        uwnd_out(im,1,k) = uwnd_out(im,1,k)*.5*(fsm(im,2)+fsm(imm1,1))
!        vwnd_out(im,1,k) = vwnd_out(im,1,k)*.5*(fsm(im,2)+fsm(imm1,1))
!        uwnd_out(1,jm,k) = uwnd_out(1,jm,k)*.5*(fsm(1,jmm1)+fsm(2,jm))
!        vwnd_out(1,jm,k) = vwnd_out(1,jm,k)*.5*(fsm(1,jmm1)+fsm(2,jm))
!        uwnd_out(im,jm,k) = uwnd_out(im,jm,k)*.5
!     &                     *(fsm(im,jmm1)+fsm(imm1,jm))
!        vwnd_out(im,jm,k) = vwnd_out(im,jm,k)*.5
!     &                     *(fsm(im,jmm1)+fsm(imm1,jm))
!      end do

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end

!_______________________________________________________________________
      subroutine read_heat_pnetcdf( sh, lh, sr, lr, dlr, filename, n )
! read heat fluxes
        use mpi,     only: mpi_offset_kind, mpi_info_null
        use pnetcdf, only: nf90mpi_open, nf90mpi_inq_varid
     &                   , nf90mpi_close
     &                   , NF_NOERR, NF_NOWRITE

        implicit none

        include 'pom.h'
        include 'io.h'

        integer get_var_real_2d
        real(rk), dimension(im,jm), intent(out) :: sh, lh, sr, lr, dlr
        character(len=*),           intent(in)  :: filename
        integer,                    intent(in)  :: n

        character(len=64) netcdf_heat_file
        integer ncid,status
        integer sh_varid, lh_varid, sr_varid, lr_varid, dlr_varid
        integer(mpi_offset_kind) start(4),edge(4)

        start = 1
        edge  = 1

! open netcdf file
        netcdf_heat_file = trim(flux_path)//trim(filename)//".nc"

        if ( my_task == master_task )
     &       write(*,'(/''reading heat '',a,'' @ '',i4)')
     &                                    trim(netcdf_heat_file),n
        status = nf90mpi_open( pom_comm, netcdf_heat_file, NF_NOWRITE
     &                       , mpi_info_null, ncid )
        call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_heat_file
     &                           , status )

! get variables
        status = nf90mpi_inq_varid( ncid, sht_name, sh_varid )
        call handle_error_pnetcdf(
     &         'nfmpi_inq_varid: '//trim(sht_name),status)
        status = nf90mpi_inq_varid( ncid, lht_name, lh_varid )
        call handle_error_pnetcdf(
     &         'nfmpi_inq_varid: '//trim(lht_name),status)
        status = nf90mpi_inq_varid( ncid, swr_name, sr_varid )
        call handle_error_pnetcdf(
     &         'nfmpi_inq_varid: '//trim(swr_name),status)
        status = nf90mpi_inq_varid( ncid, lwr_name, lr_varid )
        call handle_error_pnetcdf(
     &         'nfmpi_inq_varid: '//trim(lwr_name),status)
        status = nf90mpi_inq_varid( ncid, dlwr_name, dlr_varid )
        call handle_error_pnetcdf(
     &        'nfmpi_inq_varid: '//trim(dlwr_name),status)
! get data
        start(1) = i_global(1)
        start(2) = j_global(1)
        start(3) = n
        edge(1) = im
        edge(2) = jm
        edge(3) = 1

        status = get_var_real_2d(ncid,sh_varid,start,edge,sh)
        call handle_error_pnetcdf(
     &               'get_var_real: '//trim(sht_name),status)

        status = get_var_real_2d(ncid,lh_varid,start,edge,lh)
        call handle_error_pnetcdf(
     &               'get_var_real: '//trim(lht_name),status)

        status = get_var_real_2d(ncid,sr_varid,start,edge,sr)
        call handle_error_pnetcdf(
     &               'get_var_real: '//trim(swr_name),status)

        status = get_var_real_2d(ncid,lr_varid,start,edge,lr)
        call handle_error_pnetcdf(
     &               'get_var_real: '//trim(lwr_name),status)

        status = get_var_real_2d(ncid,dlr_varid,start,edge,dlr)
        call handle_error_pnetcdf(
     &              'get_var_real: '//trim(dlwr_name),status)

! close file
        status = nf90mpi_close(ncid)
        call handle_error_pnetcdf( 'nf_close', status )

        return

      end
!_______________________________________________________________________
      subroutine read_surf_corr_pnetcdf( dtemp, dsalt, filename, n )
! read heat fluxes
        use mpi,     only: mpi_offset_kind, mpi_info_null
        use pnetcdf, only: nf90mpi_open, nf90mpi_inq_varid
     &                   , nf90mpi_close
     &                   , NF_NOERR, NF_NOWRITE

        implicit none

        include 'pom.h'
        include 'io.h'

        integer get_var_real_2d
        real(rk), dimension(im,jm), intent(out) :: dtemp, dsalt
        character(len=*),           intent(in)  :: filename
        integer,                    intent(in)  :: n

        character(len=64) netcdf_tmp_file
        integer ncid,status
        integer dtemp_varid, dsalt_varid
        integer(mpi_offset_kind) start(4),edge(4)

        start = 1
        edge  = 1

! open netcdf file
        netcdf_tmp_file = 'in/'//trim(filename)

        if ( my_task == master_task )
     &       write(*,'(/''reading heat '',a,'' @ '',i4)')
     &                                    trim(netcdf_tmp_file),n
        status = nf90mpi_open( pom_comm, netcdf_tmp_file, NF_NOWRITE
     &                       , mpi_info_null, ncid )
        call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_tmp_file
     &                           , status )

! get variables
        status = nf90mpi_inq_varid(ncid,"dt",dtemp_varid)
        call handle_error_pnetcdf( 'nfmpi_inq_varid: dtemp', status )
        status = nf90mpi_inq_varid(ncid,"ds",dsalt_varid)
        call handle_error_pnetcdf( 'nfmpi_inq_varid: dsalt', status )
! get data
        start(1) = i_global(1)
        start(2) = j_global(1)
        start(3) = n
        edge(1) = im
        edge(2) = jm
        edge(3) = 1

        status = get_var_real_2d(ncid,dtemp_varid,start,edge,dtemp)
        call handle_error_pnetcdf(
     &               'get_var_real: dtemp',status)

        status = get_var_real_2d(ncid,dsalt_varid,start,edge,dsalt)
        call handle_error_pnetcdf(
     &               'get_var_real: dsalt',status)

! close file
        status = nf90mpi_close(ncid)
        call handle_error_pnetcdf( 'nf_close', status )

        return

      end
!_______________________________________________________________________
      subroutine read_heat_aux_pnetcdf( shm, prt, cld, filename, n )
! read wind stress
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer get_var_real_2d
      integer n
      real(kind=rk), dimension(im,jm) :: shm, prt, cld
      character(len=14) filename
      character(len=24) netcdf_heat_aux_file
      integer ncid,status
      integer shm_varid, prt_varid, cld_varid
      integer(mpi_offset_kind) start(4),edge(4)

      start = 1
      edge  = 1


! open netcdf wind file
      netcdf_heat_aux_file="in/heat/"//trim(filename)
      if(my_task.eq.0) write(*,'(/''reading aux heat '',a,'' @ '',i4)')
     $                                    trim(netcdf_heat_aux_file),n
      status=nfmpi_open(pom_comm,netcdf_heat_aux_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_heat_aux_file
     &                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'rhum',shm_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: rhum',status)
      status=nfmpi_inq_varid(ncid,'rain',prt_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: prate',status)
      status=nfmpi_inq_varid(ncid,'cloud',cld_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: cloud',status)

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=n
      edge(1)=im
      edge(2)=jm
      edge(3)=1

      status = get_var_real_2d(ncid,shm_varid,start,edge,shm)
      call handle_error_pnetcdf('get_var_real: rhum',status)
      shm = shm*fsm

      status = get_var_real_2d(ncid,prt_varid,start,edge,prt)
      call handle_error_pnetcdf('get_var_real: prate',status)
      prt = prt*fsm

      status = get_var_real_2d(ncid,cld_varid,start,edge,cld)
      call handle_error_pnetcdf('get_var_real: cloud',status)
      cld = cld*fsm

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end
!_______________________________________________________________________
      subroutine read_ncep_bulk_pnetcdf( pres, tair, rhum, rain, cloud
     $                                  ,uwnd, vwnd, tskn, filename, n )
! read bulk data
        use mpi,     only: mpi_offset_kind, mpi_info_null
        use pnetcdf, only: nf90mpi_open, nf90mpi_inq_varid
     &                   , nf90mpi_get_att, nf90mpi_close
     &                   , NF_NOERR, NF_NOWRITE

      implicit none

      include 'pom.h'
      include 'io.h'

      integer get_var_real_2d
      real(rk), dimension(im,jm), intent(out) :: pres, tair, rhum, rain
     &                                         , cloud, uwnd, vwnd, tskn
      character(len=*),           intent(in)  :: filename
      integer,                    intent(in)  :: n

      character(len=256) netcdf_read_file, att_val
      integer ncid,status
      integer pres_varid, tair_varid, rhum_varid, rain_varid
     $       ,cloud_varid, uwnd_varid, vwnd_varid, tskn_varid
      integer(mpi_offset_kind) start(4),edge(4)
      logical lexist


      start = 1
      edge  = 1

! open netcdf wind file
      netcdf_read_file = trim(bulk_path)//trim(filename)
      inquire( file=netcdf_read_file, exist=lexist )

      if ( .not.lexist ) then

        if ( my_task == master_task ) then
          print *, "[!] Missing surface bulk data file at `"
     &             , trim(netcdf_read_file), "`"
        end if
        return

      end if

      if ( my_task == master_task ) then
        write(*,'(/''reading bulk data '',a,'' @ '',i4)')
     $                                    trim(netcdf_read_file),n
      end if

      status = nf90mpi_open( pom_comm, netcdf_read_file, nf_nowrite,
     &                       mpi_info_null, ncid )
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_read_file
     &                         , status )

! get variables
      status = nf90mpi_inq_varid( ncid, trim(pres_name), pres_varid )
      call handle_error_pnetcdf('nfmpi_inq_varid: '//trim(pres_name)
     $                          ,status)
      status = nf90mpi_inq_varid( ncid, trim(tair_name), tair_varid )
      call handle_error_pnetcdf('nfmpi_inq_varid: '//trim(tair_name)
     $                          ,status)
      status = nf90mpi_inq_varid( ncid, trim(rhum_name), rhum_varid )
      call handle_error_pnetcdf('nfmpi_inq_varid: '//trim(rhum_name)
     $                          ,status)
      status = nf90mpi_inq_varid( ncid, trim(prat_name), rain_varid )
      call handle_error_pnetcdf('nfmpi_inq_varid: '//trim(prat_name)
     $                          ,status)
      status = nf90mpi_inq_varid( ncid, trim(tcld_name), cloud_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: '//trim(tcld_name)
     $                          ,status)
      status = nf90mpi_inq_varid( ncid, trim(uwnd_name), uwnd_varid )
      call handle_error_pnetcdf('nfmpi_inq_varid: '//trim(uwnd_name)
     $                          ,status)
      status = nf90mpi_inq_varid( ncid, trim(vwnd_name), vwnd_varid )
      call handle_error_pnetcdf('nfmpi_inq_varid: '//trim(vwnd_name)
     $                          ,status)
      status = nf90mpi_inq_varid( ncid, trim(sst_name), tskn_varid )
      call handle_error_pnetcdf('nfmpi_inq_varid: '//trim(sst_name)
     $                          ,status)

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=n
      edge(1)=im
      edge(2)=jm
      edge(3)=1

      status = get_var_real_2d( ncid, pres_varid, start,edge, pres )
      call handle_error_pnetcdf('get_var_real: '//trim(pres_name)
     &                          ,status)
      status = nf90mpi_get_att( ncid, pres_varid, 'units', att_val )
      if ( status == nf_noerr ) then
        if ( trim(att_val) == 'Pa' ) pres = pres/100.
      end if

      status = get_var_real_2d( ncid, tair_varid, start,edge, tair )
      call handle_error_pnetcdf('get_var_real: '//trim(tair_name)
     &                          ,status)
      status = nf90mpi_get_att( ncid, tair_varid, 'units', att_val )
      if ( status == nf_noerr ) then
        if ( index(trim(att_val), 'K') /= 0 ) tair = tair - 273.15
      end if

      status = get_var_real_2d( ncid, rhum_varid, start,edge, rhum )
      call handle_error_pnetcdf('get_var_real: '//trim(rhum_name)
     &                          ,status)

      status = get_var_real_2d( ncid, rain_varid, start,edge, rain )
      call handle_error_pnetcdf('get_var_real: '//trim(prat_name)
     &                          ,status)

      status = get_var_real_2d( ncid, cloud_varid, start,edge, cloud )
      call handle_error_pnetcdf('get_var_real: '//trim(tcld_name)
     &                          ,status)
      where( cloud < 0. ) cloud = 0.

      status = get_var_real_2d( ncid, uwnd_varid, start,edge, uwnd )
      call handle_error_pnetcdf('get_var_real: '//trim(uwnd_name)
     &                          ,status)

      status = get_var_real_2d( ncid, vwnd_varid, start,edge, vwnd )
      call handle_error_pnetcdf('get_var_real: '//trim(vwnd_name)
     &                          ,status)

      status = get_var_real_2d( ncid, tskn_varid, start,edge, tskn )
      call handle_error_pnetcdf('get_var_real: '//trim(pres_name)
     &                          ,status)
      status = nf90mpi_get_att( ncid, tskn_varid, 'units', att_val )
      if ( status == nf_noerr ) then
        if ( index(trim(att_val), 'K') /= 0 ) tskn = tskn - 273.15
      end if

! close file
      status = nf90mpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end
!_______________________________________________________________________
      subroutine read_ice_pnetcdf( ci, filename )
! read wind stress
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer get_var_real_2d
      real(kind=rk), dimension(im,jm) :: ci
      character(len=15) filename
      character(len=24) netcdf_sice_file
      integer ncid,status
      integer ci_varid
      integer(mpi_offset_kind) start(4),edge(4)

      start = 1
      edge  = 1


! open netcdf wind file
      netcdf_sice_file="in/sice/"//trim(filename)
      if(my_task.eq.0) write(*,'(/''reading ice '',a)')
     $                                            trim(netcdf_sice_file)
      status=nfmpi_open(pom_comm_coarse,netcdf_sice_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_sice_file
     $                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'ci',ci_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: sice'
     $                         , status )
! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=1

      status = get_var_real_2d(ncid,ci_varid,start,edge,ci)
      call handle_error_pnetcdf('get_var_real: sice',status)
!      where(ci<=minval(ci,ci>0.)) ci = 0.
      ci = fsm*ci/100.
!      write(*,*) my_task,"::",minval(ci, ci>0.),maxval(ci)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end

!_______________________________________________________________________
      subroutine read_msla_pnetcdf( ssha_out, filename )
! read mean sea level anomaly data
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      real(kind=rk)ssha_out(im,jm)
      real(kind=rk)errval
      character :: filename*16, netcdf_msla_file*25
      integer ncid,status
      integer ssha_varid
      integer(mpi_offset_kind) start(3),edge(3)


! open netcdf msla file
      write(netcdf_msla_file,'(''in/assim/'',a )') trim(filename)
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $                                    trim(netcdf_msla_file)
      status=nfmpi_open(pom_comm,netcdf_msla_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_msla_file
     $                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'ssha',ssha_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: ssha', status )
      status=nfmpi_get_att_real(ncid, ssha_varid,
     $     'missing_value', errval )

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=1

      status=nfmpi_get_vara_real_all(
     $     ncid,ssha_varid,start,edge,ssha_out)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     &                         , status )


! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end
!_______________________________________________________________________
!_______________________________________________________________________
      subroutine read_msla_pnetcdfc( ssha_out, filename )
! read mean sea level anomaly data
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      real(kind=rk)ssha_out(im_coarse,jm_coarse)
      real(kind=rk)errval
      character :: filename*16, netcdf_msla_file*25
      integer ncid,status
      integer ssha_varid
      integer(mpi_offset_kind) start(3),edge(3)


! open netcdf msla file
      write(netcdf_msla_file,'(''in/assim/'',a )') trim(filename)
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $                                    trim(netcdf_msla_file)
      status=nfmpi_open(pom_comm_coarse,netcdf_msla_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_msla_file
     &                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'ssha',ssha_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: ssha'
     $                         , status )
      status=nfmpi_get_att_real(ncid, ssha_varid,
     $     'missing_value', errval )
!      if(my_task.eq.0) write(*,'(a,f12.3)')
!     $                                  " missing_value = ", errval


! get data
      start(1)=i_global_coarse(1)
      start(2)=j_global_coarse(1)
      start(3)=1
      edge(1)=im_coarse
      edge(2)=jm_coarse
      edge(3)=1

      status=nfmpi_get_vara_real_all(
     $     ncid,ssha_varid,start,edge,ssha_out)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     $                         , status )


! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end
!_______________________________________________________________________

!_______________________________________________________________________
      subroutine read_assiminfo_pnetcdf
     $     ( frs_out, tav_out, fac_out, cof_out )
! read information for data assimilation
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      real(kind=rk), dimension(im_local,jm_local) :: frs_out
      real(kind=rk), dimension(im_local,jm_local,kb)
     $     :: tav_out, fac_out, cof_out
      character :: netcdf_assiminfo_file*21
      integer ncid,status
      integer frs_varid, tav_varid, fac_varid, cof_varid
      integer(mpi_offset_kind) start(3),edge(3)


! open netcdf assiminfo file
      write(netcdf_assiminfo_file,'(''in/assim/'',a )') "assiminfo.nc"
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $                                    trim(netcdf_assiminfo_file)
      status=nfmpi_open(pom_comm,netcdf_assiminfo_file,nf_nowrite,
     $                  mpi_info_null,ncid)

! get variables
      status=nfmpi_inq_varid(ncid,'frs',frs_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: frs', status )
      status=nfmpi_inq_varid(ncid,'tav',tav_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: tav', status )
      status=nfmpi_inq_varid(ncid,'fac',fac_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: fac', status )
      status=nfmpi_inq_varid(ncid,'cof',cof_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: cof', status )



! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=1

      status=nfmpi_get_vara_real_all(
     $     ncid,frs_varid,start,edge,frs_out)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     $                         , status )

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=kb

      status=nfmpi_get_vara_real_all(
     $     ncid,tav_varid,start,edge,tav_out)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     $                         , status )

      status=nfmpi_get_vara_real_all(
     $     ncid,fac_varid,start,edge,fac_out)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     $                         , status )

      status=nfmpi_get_vara_real_all(
     $     ncid,cof_varid,start,edge,cof_out)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     $                         , status )


! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end
!_______________________________________________________________________
!_______________________________________________________________________
      subroutine read_assiminfo_pnetcdfc
     $     ( frs_out, tav_out, fac_out, cof_out )
! read information for data assimilation
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      real(kind=rk),dimension(im_local_coarse,jm_local_coarse)::frs_out
      real(kind=rk),dimension(im_local_coarse,jm_local_coarse,kb)
     $     :: tav_out, fac_out, cof_out
      character :: netcdf_assiminfo_file*21
      integer ncid,status
      integer frs_varid, tav_varid, fac_varid, cof_varid
      integer(mpi_offset_kind) start(3),edge(3)


! open netcdf assiminfo file
      write(netcdf_assiminfo_file,'(''in/assim/'',a )') "assiminfo.nc"
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $                                    trim(netcdf_assiminfo_file)
      status=nfmpi_open(pom_comm_coarse,netcdf_assiminfo_file,
     $                 nf_nowrite,mpi_info_null,ncid)

! get variables
      status=nfmpi_inq_varid(ncid,'frs',frs_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: frs'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'tav',tav_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: tav'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'fac',fac_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: fac'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'cof',cof_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: cof'
     $                          ,status)



! get data
      start(1)=i_global_coarse(1)
      start(2)=j_global_coarse(1)
      start(3)=1
      edge(1)=im_coarse
      edge(2)=jm_coarse
      edge(3)=1

      status=nfmpi_get_vara_real_all(
     $     ncid,frs_varid,start,edge,frs_out)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     &                         , status )

      start(1)=i_global_coarse(1)
      start(2)=j_global_coarse(1)
      start(3)=1
      edge(1)=im_coarse
      edge(2)=jm_coarse
      edge(3)=kb

      status=nfmpi_get_vara_real_all(
     $     ncid,tav_varid,start,edge,tav_out)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     $                         , status )

      status=nfmpi_get_vara_real_all(
     $     ncid,fac_varid,start,edge,fac_out)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     $                         , status )

      status=nfmpi_get_vara_real_all(
     $     ncid,cof_varid,start,edge,cof_out)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     $                         , status )


! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end
!_______________________________________________________________________
!_______________________________________________________________________
!fhx:mcsst:start
      subroutine read_mcsst_pnetcdf( mcsst_out, filename )
! read mean MCSST data
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      real(kind=rk) mcsst_out(im_local,jm_local)
      real(kind=rk) errval
      character :: filename*17, netcdf_mcsst_file*26
      integer ncid,status
      integer mcsst_varid
      integer(mpi_offset_kind) start(3),edge(3)


! open netcdf mcsst file
      write(netcdf_mcsst_file,'(''in/mcsst/'',a )') trim(filename)
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $                                    trim(netcdf_mcsst_file)
      status=nfmpi_open(pom_comm,netcdf_mcsst_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nfmpi_open: '//netcdf_mcsst_file
     &                         , status )

! get variables
      status=nfmpi_inq_varid(ncid,'mcsst',mcsst_varid)
      call handle_error_pnetcdf( 'nfmpi_inq_varid: mcsst'
     &                         , status )
      status=nfmpi_get_att_real(ncid, mcsst_varid,
     $     'missing_value', errval )
!      if(my_task.eq.0) write(*,'(a,f12.3)')
!     $                                  " missing_value = ", errval


! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=1

      status=nfmpi_get_vara_real_all(
     $     ncid,mcsst_varid,start,edge,mcsst_out)
      call handle_error_pnetcdf( 'nfmpi_get_vara_real_all'
     $                         , status )


! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status)

      return
      end
!fhx:mcsst:end
!_______________________________________________________________________
!_______________________________________________________________________

      subroutine write_output_pnetcdf( netcdf_out_file ) !lyo:pac10:use stcc version
! write output data
        use mpi,     only: MPI_OFFSET_KIND, MPI_INFO_NULL
        use pnetcdf, only: nf90mpi_close    , nf90mpi_create
     &                   , nf90mpi_def_dim  , nf90mpi_enddef
     &                   , nf90mpi_inq_varid, nf90mpi_get_att
     &                   , nf90mpi_open     , nf90mpi_put_att
     &                   , nfmpi_put_vara_all,nfmpi_put_vara_real_all
     &                   , NF_64BIT_OFFSET, NF_CLOBBER
     &                   , NF_FLOAT       , NF_GLOBAL
     &                   , NF_NOERR       , NF_NOWRITE
     &                   , NF_UNLIMITED   , NF_WRITE

        implicit none

        include 'pom.h'

        character(len=*), intent(in) :: netcdf_out_file

        character(len=120) str_tmp!, netcdf_out_file
        integer time_dimid, x_dimid, y_dimid, z_dimid
        integer       z_varid,     zz_varid,     dx_varid,     dy_varid
     &         , east_c_varid, east_e_varid, east_u_varid, east_v_varid
     &         ,north_c_varid,north_e_varid,north_u_varid,north_v_varid
     &         ,    rot_varid,      h_varid,    fsm_varid,    dum_varid
     &         ,    dvm_varid,  uhtfl_varid, wusurf_varid, wvsurf_varid
     &         , wtsurf_varid, wssurf_varid,   icec_varid,     ui_varid
     &         , vi_varid    ,   time_varid,    uab_varid,    vab_varid
     &         ,    elb_varid,      u_varid,      v_varid,      w_varid
     &         ,      t_varid,      s_varid,    rho_varid,     km_varid
     &         ,     kh_varid
        integer ncid,status
        integer vdims(4)
        integer(MPI_OFFSET_KIND) start(4),edge(4)

        real(kind=4), dimension(      kb) :: out1
        real(kind=4), dimension(im,jm   ) :: out2
        real(kind=4), dimension(im,jm,kb) :: out3


        out1 = 0.
        out2 = 0.
        out3 = 0.

        iout = iout + 1

        if ( iout == 1  ) then

!     create netcdf file
          if ( my_task == master_task )
     &         write(*,'(/''writing file '',a)') trim(netcdf_out_file)
          status = nf90mpi_create( pom_comm, trim(netcdf_out_file)
     &           , NF_CLOBBER+NF_64BIT_OFFSET, MPI_INFO_NULL, ncid )
          call handle_error_pnetcdf( 'nf_create: '//netcdf_out_file
     &                             , status )

! define global attributes
          status = nf90mpi_put_att( ncid, NF_GLOBAL
     &                            , 'title', trim(title) )
          call handle_error_pnetcdf('nf_put_att: title',status)

          status = nf90mpi_put_att( ncid, NF_GLOBAL
     &                            , 'description', 'Output file' )
          call handle_error_pnetcdf('nf_put_att: description',status)

! define dimensions
          status = nf90mpi_def_dim(ncid,'time'
     &                            ,int(NF_UNLIMITED,8),time_dimid)
          call handle_error_pnetcdf('nf_def_dim: time',status)
          status = nf90mpi_def_dim(ncid,   'z'
     &                            ,int(          kb,8),   z_dimid)
          call handle_error_pnetcdf('nf_def_dim: z'   ,status)
          status = nf90mpi_def_dim(ncid,   'y'
     &                            ,int(   jm_global,8),   y_dimid)
          call handle_error_pnetcdf('nf_def_dim: y'   ,status)
          status = nf90mpi_def_dim(ncid,   'x'
     &                            ,int(   im_global,8),   x_dimid)
          call handle_error_pnetcdf('nf_def_dim: x'   ,status)

! define variables and their attributes
          vdims(1) = time_dimid
          str_tmp  = 'days since '//time_start
          call def_var_pnetcdf(ncid,'time',1,vdims,time_varid,NF_FLOAT
     &                        ,'time',str_tmp
     &                        ,-1,0.,' ',.false.)

          vdims(1) = z_dimid
          call def_var_pnetcdf(ncid,   'z',1,vdims,   z_varid,NF_FLOAT
     &                        ,'sigma of cell face','sigma_level'
     &                        ,-1,0.,' ',.false.)
          status = nf90mpi_put_att( ncid, z_varid
     &                            , 'standard_name'
     &                            , 'ocean_sigma_coordinate' )
          call handle_error_pnetcdf('nf_put_att: z@std_name',status)
          status = nf90mpi_put_att( ncid, z_varid
     &                            , 'formula_terms'
     &                            , 'sigma: z eta: elb depth: h' )
          call handle_error_pnetcdf('nf_put_att: z@frm_terms',status)

          call def_var_pnetcdf(ncid,  'zz',1,vdims,  zz_varid,NF_FLOAT
     &                        ,'sigma of cell centre','sigma_level'
     &                        ,-1,0.,' ',.false.)
          status = nf90mpi_put_att( ncid, zz_varid
     &                            , 'standard_name'
     &                            , 'ocean_sigma_coordinate' )
          call handle_error_pnetcdf('nf_put_att: zz@std_name',status)
          status = nf90mpi_put_att( ncid, zz_varid
     &                            , 'formula_terms'
     &                            , 'sigma: zz eta: elb depth: h' )
          call handle_error_pnetcdf('nf_put_att: zz@frm_terms',status)

          vdims(1) = x_dimid
          vdims(2) = y_dimid
          call def_var_pnetcdf(ncid,'dx',2,vdims,dx_varid,NF_FLOAT
     &                        ,'grid increment in x','metre'
     &                        ,-1,0.,'east_e north_e',.true.)
          call def_var_pnetcdf(ncid,'dy',2,vdims,dy_varid,NF_FLOAT
     &                        ,'grid increment in y','metre'
     &                        ,-1,0.,'east_e north_e',.true.)
          call def_var_pnetcdf(ncid,'east_u',2,vdims
     &                        ,east_u_varid,NF_FLOAT
     &                        ,'easting of u-points','degree'
     &                        ,-1,0.,'east_u north_u',.true.)
          call def_var_pnetcdf(ncid,'east_v',2,vdims
     &                        ,east_v_varid,NF_FLOAT
     &                        ,'easting of v-points','degree'
     &                        ,-1,0.,'east_v north_v',.true.)
          call def_var_pnetcdf(ncid,'east_e',2,vdims
     &                        ,east_e_varid,NF_FLOAT
     &                        ,'easting of elevation points','degree'
     &                        ,-1,0.,'east_e north_e',.true.)
          call def_var_pnetcdf(ncid,'east_c',2,vdims
     &                        ,east_c_varid,NF_FLOAT
     &                        ,'easting of cell corners','degree'
     &                        ,-1,0.,'east_c north_c',.true.)
          call def_var_pnetcdf(ncid,'north_u',2,vdims
     &                        ,north_u_varid,NF_FLOAT
     &                        ,'northing of u-points','degree'
     &                        ,-1,0.,'east_u north_u',.true.)
          call def_var_pnetcdf(ncid,'north_v',2,vdims
     &                        ,north_v_varid,NF_FLOAT
     &                        ,'northing of v-points','degree'
     &                        ,-1,0.,'east_v north_v',.true.)
          call def_var_pnetcdf(ncid,'north_e',2,vdims
     &                        ,north_e_varid,NF_FLOAT
     &                        ,'northing of elevation points','degree'
     &                        ,-1,0.,'east_e north_e',.true.)
          call def_var_pnetcdf(ncid,'north_c',2,vdims
     &                        ,north_c_varid,NF_FLOAT
     &                        ,'northing of cell corners','degree'
     &                        ,-1,0.,'east_c north_c',.true.)
          call def_var_pnetcdf(ncid,'rot',2,vdims
     &                        ,rot_varid,NF_FLOAT
     &                        ,'Rotation angle of x-axis wrt. east'
     &                        ,'degree'
     &                        ,-1,0.,'east_e north_e',.true.)
          call def_var_pnetcdf(ncid,'h',2,vdims
     &                        ,h_varid,NF_FLOAT
     &                        ,'undisturbed water depth','metre'
     &                        ,-1,0.,'east_e north_e',.true.)
          call def_var_pnetcdf(ncid,'fsm',2,vdims
     &                        ,fsm_varid,NF_FLOAT
     &                        ,'free surface mask','dimensionless'
     &                        ,-1,0.,'east_e north_e',.true.)
          call def_var_pnetcdf(ncid,'dum',2,vdims
     &                        ,dum_varid,NF_FLOAT
     &                        ,'u-velocity mask','dimensionless'
     &                        ,-1,0.,'east_u north_u',.true.)
          call def_var_pnetcdf(ncid,'dvm',2,vdims
     &                        ,dvm_varid,NF_FLOAT
     &                        ,'v-velocity mask','dimensionless'
     &                        ,-1,0.,'east_v north_v',.true.)

          vdims(1) = x_dimid
          vdims(2) = y_dimid
          vdims(3) = time_dimid
          call def_var_pnetcdf(ncid,'uab',3,vdims
     &                        ,uab_varid,NF_FLOAT
     &                        ,'depth-averaged u','metre/sec'
     &                        ,-1,0.,'east_u north_u',.true.)
          call def_var_pnetcdf(ncid,'vab',3,vdims
     &                        ,vab_varid,NF_FLOAT
     &                        ,'depth-averaged v','metre/sec'
     &                        ,-1,0.,'east_v north_v',.true.)
          call def_var_pnetcdf(ncid,'elb',3,vdims
     &                        ,elb_varid,NF_FLOAT
     &                        ,'surface elevation','metre'
     &                        ,0,-999.,'east_e north_e',.true.)
          call def_var_pnetcdf(ncid,'wusurf',3,vdims
     &                        ,wusurf_varid,NF_FLOAT
     &                        ,'x-momentum flux','metre^2/sec^2'
     &                        ,-1,0.,'east_u north_u',.true.)
          call def_var_pnetcdf(ncid,'wvsurf',3,vdims
     &                        ,wvsurf_varid,NF_FLOAT
     &                        ,'y-momentum flux','metre^2/sec^2'
     &                        ,-1,0.,'east_v north_v',.true.)
          call def_var_pnetcdf(ncid,'wtsurf',3,vdims
     &                        ,wtsurf_varid,NF_FLOAT
     &                        ,'temperature flux','deg m/s'
     &                        ,-1,0.,'east_e north_e',.true.)
          call def_var_pnetcdf(ncid,'swrad',3,vdims
     &                        ,uhtfl_varid,NF_FLOAT
     &                        ,'upward net heat flux','K m/s'
     &                        ,-1,0.,'east_e north_e',.true.)
          call def_var_pnetcdf(ncid,'wssurf',3,vdims
     &                        ,wssurf_varid,NF_FLOAT
     &                        ,'salinity flux','psu m/s'
     &                        ,-1,0.,'east_e north_e',.true.)
          if ( calc_ice ) then
            call def_var_pnetcdf(ncid,'icec',3,vdims
     &                          ,icec_varid,NF_FLOAT
     &                          ,'sea ice concentration'
     &                          ,'fraction'
     &                          ,0,0.,'east_e north_e',.true.)
            call def_var_pnetcdf(ncid,'ui',3,vdims
     &                          ,ui_varid,NF_FLOAT
     &                          ,'sea ice x-velocity','m/s'
     &                          ,-1,0.,'east_u north_u',.true.)
            call def_var_pnetcdf(ncid,'vi',3,vdims
     &                         ,vi_varid,NF_FLOAT
     &                         ,'sea ice y-velocity','m/s'
     &                         ,-1,0.,'east_v north_v',.true.)
          end if

          if ( mode/= 2 ) then
            vdims(1) = x_dimid
            vdims(2) = y_dimid
            vdims(3) = z_dimid
            vdims(4) = time_dimid
            call def_var_pnetcdf(ncid,'u',4,vdims
     &                          ,u_varid,NF_FLOAT
     &                          ,'x-velocity','metre/sec'
     &                          ,-1,0.,'east_u north_u zz',.true.)
            call def_var_pnetcdf(ncid,'v',4,vdims
     &                          ,v_varid,NF_FLOAT
     &                          ,'y-velocity','metre/sec'
     &                          ,-1,0.,'east_v north_v zz',.true.)
            call def_var_pnetcdf(ncid,'w',4,vdims
     &                          ,w_varid,NF_FLOAT
     &                          ,'z-velocity','metre/sec'
     &                          ,-1,0.,'east_e north_e z',.true.)
            call def_var_pnetcdf(ncid,'t',4,vdims
     &                          ,t_varid,NF_FLOAT
     &                          ,'potential temperature','degC'
     &                          ,0,-999.,'east_e north_e zz',.true.)
            call def_var_pnetcdf(ncid,'s',4,vdims
     &                          ,s_varid,NF_FLOAT
     &                          ,'salinity x rho / rhoref','PSS'
     &                          ,0,0.,'east_e north_e zz',.true.)
            call def_var_pnetcdf(ncid,'rho',4,vdims
     &                          ,rho_varid,NF_FLOAT
     &                          ,'(density-1000)/rhoref'
     &                          ,'dimensionless'
     &                          ,-1,0.,'east_e north_e zz',.true.)
            call def_var_pnetcdf(ncid,'kh',4,vdims
     &                          ,kh_varid,NF_FLOAT
     &                          ,'vertical diffusivity','metre2/sec'
     &                          ,-1,0.,'east_e north_e zz',.true.)
            call def_var_pnetcdf(ncid,'km',4,vdims
     &                          ,km_varid,NF_FLOAT
     &                          ,'vertical viscosity','metre2/sec'
     &                          ,-1,0.,'east_e north_e zz',.true.)
          end if

! end definitions
          status = nf90mpi_enddef(ncid)
          call handle_error_pnetcdf('nf_enddef: output',status)

! write data
!      start(1)=1
!      edge(1)=1
!      status=nfmpi_put_vara_real_all(ncid,time_varid,start,edge,time)
!      call handle_error_pnetcdf('nf_put_vara_real:time',status)

          start(1) = 1
          edge(1) = kb

          out1 = real(z,4)
          status=nfmpi_put_vara_real_all(ncid, z_varid,start,edge,out1)
          call handle_error_pnetcdf('nf_put_vara: z',status)

          out1 = real(zz,4)
          status=nfmpi_put_vara_real_all(ncid,zz_varid,start,edge,out1)
          call handle_error_pnetcdf('nf_put_vara: zz',status)


          start(1) = i_global(1)
          start(2) = j_global(1)
          edge(1) = im
          edge(2) = jm

          out2 = real(dx,4)
          status=nfmpi_put_vara_real_all(ncid,dx_varid,start,edge,out2)
          call handle_error_pnetcdf('nf_put_var_real: dx',status)
          out2 = real(dy,4)
          status=nfmpi_put_vara_real_all(ncid,dy_varid,start,edge,out2)
          call handle_error_pnetcdf('nf_put_var_real: dy',status)
          out2 = real(east_u,4)
          status=nfmpi_put_vara_real_all(ncid,east_u_varid,start,edge
     &                                                           ,out2)
          call handle_error_pnetcdf('nf_put_var_real: east_u',status)
          out2 = real(east_v,4)
          status=nfmpi_put_vara_real_all(ncid,east_v_varid,start,edge
     &                                                           ,out2)
          call handle_error_pnetcdf('nf_put_var_real: east_v',status)
          out2 = real(east_e,4)
          status=nfmpi_put_vara_real_all(ncid,east_e_varid,start,edge
     &                                                           ,out2)
          call handle_error_pnetcdf('nf_put_var_real: east_e',status)
          out2 = real(east_c,4)
          status=nfmpi_put_vara_real_all(ncid,east_c_varid,start,edge
     &                                                           ,out2)
          call handle_error_pnetcdf('nf_put_var_real: east_c',status)
          out2 = real(north_u,4)
          status=nfmpi_put_vara_real_all(ncid,north_u_varid,start,edge
     &                                                           ,out2)
          call handle_error_pnetcdf('nf_put_var_real: north_u',status)
          out2 = real(north_v,4)
          status=nfmpi_put_vara_real_all(ncid,north_v_varid,start,edge
     &                                                           ,out2)
          call handle_error_pnetcdf('nf_put_var_real: north_v',status)
          out2 = real(north_e,4)
          status=nfmpi_put_vara_real_all(ncid,north_e_varid,start,edge
     &                                                           ,out2)
          call handle_error_pnetcdf('nf_put_var_real: north_e',status)
          out2 = real(north_c,4)
          status=nfmpi_put_vara_real_all(ncid,north_c_varid,start,edge
     &                                                           ,out2)
          call handle_error_pnetcdf('nf_put_var_real: north_c',status)
          out2 = real(rot,4)
          status=nfmpi_put_vara_real_all(ncid,rot_varid,start,edge
     &                                                           ,out2)
          call handle_error_pnetcdf('nf_put_var_real: rot',status)
          out2 = real(h,4)
          status=nfmpi_put_vara_real_all(ncid,h_varid,start,edge,out2)
          call handle_error_pnetcdf('nf_put_var_real: h',status)
          out2 = real(fsm,4)
          status=nfmpi_put_vara_real_all(ncid,fsm_varid,start,edge
     &                                                           ,out2)
          call handle_error_pnetcdf('nf_put_var_real: fsm',status)
          out2 = real(dum,4)
          status=nfmpi_put_vara_real_all(ncid,dum_varid,start,edge
     &                                                           ,out2)
          call handle_error_pnetcdf('nf_put_var_real: dum',status)
          out2 = real(dvm,4)
          status=nfmpi_put_vara_real_all(ncid,dvm_varid,start,edge
     &                                                           ,out2)
          call handle_error_pnetcdf('nf_put_var_real: dvm',status)


        else


!     open netcdf file
          if ( my_task == master_task )
     &         write(*,'(/''writing file '',a)') trim(netcdf_out_file)
          status = nf90mpi_open( pom_comm, trim(netcdf_out_file)
     &                         , NF_WRITE, MPI_INFO_NULL, ncid )
          call handle_error_pnetcdf( 'nf_open: '//netcdf_out_file
     &                             , status )


!!     inquire dimensions
!
!      status=nfmpi_inq_dimid(ncid,'time',time_dimid)
!      call handle_error_pnetcdf('nfmpi_inq_dim: time',status)
!      status=nfmpi_inq_dim(ncid,'z',z_dimid)
!      call handle_error_pnetcdf('nfmpi_inq_dim: z',status)
!      status=nfmpi_inq_dim(ncid,'y',y_dimid)
!      call handle_error_pnetcdf('nfmpi_inq_dim: y',status)
!      status=nfmpi_inq_dim(ncid,'x',x_dimid)
!     call handle_error_pnetcdf('nfmpi_inq_dim: x',status)



!     inqure variables
          status = nf90mpi_inq_varid(ncid,'time',time_varid)
          call handle_error_pnetcdf( 'nfmpi_inq_varid: time'
     &                             , status )
          status = nf90mpi_inq_varid(ncid,'uab',uab_varid)
          call handle_error_pnetcdf( 'nfmpi_inq_varid: uab'
     &                             , status )
          status = nf90mpi_inq_varid(ncid,'vab',vab_varid)
          call handle_error_pnetcdf('nfmpi_inq_varid: vab'
     &                             , status )
          status = nf90mpi_inq_varid(ncid,'elb',elb_varid)
          call handle_error_pnetcdf( 'nfmpi_inq_varid: elb'
     &                             , status )
          status = nf90mpi_inq_varid(ncid,'wusurf',wusurf_varid)
          call handle_error_pnetcdf( 'nfmpi_inq_varid: wusurf'
     &                             , status )
          status = nf90mpi_inq_varid(ncid,'wvsurf',wvsurf_varid)
          call handle_error_pnetcdf('nfmpi_inq_varid: wvsurf'
     &                             , status )
          status = nf90mpi_inq_varid(ncid,'wtsurf',wtsurf_varid)
          call handle_error_pnetcdf('nfmpi_inq_varid: wtsurf'
     &                             , status )
          status = nf90mpi_inq_varid(ncid,'swrad',uhtfl_varid)
          call handle_error_pnetcdf('nfmpi_inq_varid: swrad'
     &                             , status )
          status = nf90mpi_inq_varid(ncid,'wssurf',wssurf_varid)
          call handle_error_pnetcdf('nfmpi_inq_varid: wssurf'
     &                             , status )
          if ( calc_ice ) then
            status = nf90mpi_inq_varid(ncid,'icec',icec_varid)
            call handle_error_pnetcdf('nfmpi_inq_varid: icec'
     &                               , status )
            status = nf90mpi_inq_varid(ncid,'ui',ui_varid)
            call handle_error_pnetcdf('nfmpi_inq_varid: ui'
     &                               , status )
            status = nf90mpi_inq_varid(ncid,'vi',vi_varid)
            call handle_error_pnetcdf('nfmpi_inq_varid: vi'
     &                               , status )
          end if
          
          if ( mode /= 2 ) then
            status = nf90mpi_inq_varid(ncid,'u',u_varid)
            call handle_error_pnetcdf('nfmpi_inq_varid: u'
     &                               , status )
            status = nf90mpi_inq_varid(ncid,'v',v_varid)
            call handle_error_pnetcdf('nfmpi_inq_varid: v'
     &                               , status )
            status = nf90mpi_inq_varid(ncid,'w',w_varid)
            call handle_error_pnetcdf('nfmpi_inq_varid: w'
     &                               , status )
            status = nf90mpi_inq_varid(ncid,'t',t_varid)
            call handle_error_pnetcdf('nfmpi_inq_varid: t'
     &                               , status )
            status = nf90mpi_inq_varid(ncid,'s',s_varid)
            call handle_error_pnetcdf('nfmpi_inq_varid: s'
     &                               , status )
            status = nf90mpi_inq_varid(ncid,'rho',rho_varid)
            call handle_error_pnetcdf('nfmpi_inq_varid: rho'
     &                               , status )
            status = nf90mpi_inq_varid(ncid,'kh',kh_varid)
            call handle_error_pnetcdf('nfmpi_inq_varid: kh'
     &                               , status )
            status = nf90mpi_inq_varid(ncid,'km',km_varid)
            call handle_error_pnetcdf('nfmpi_inq_varid: km'
     &                               , status )
          end if

        end if !if ( iout == 1  ) then...

!     write data

        if ( iout == 1  ) then !lyo:20110224:stcc:use orig t, not t_mean

          start(1) = iout
          edge(1) = 1
          out1 = real(time,4)
          status = nfmpi_put_vara_real_all(ncid,time_varid,start,edge
     &                                                            ,out1)
          call handle_error_pnetcdf('nf_put_vara_real:time',status)

          start(1) = i_global(1)
          start(2) = j_global(1)
          start(3) = iout
          edge(1) = im
          edge(2) = jm
          edge(3) = 1
          out2 = real(uab,4)
          status=nfmpi_put_vara_real_all(ncid,uab_varid,start,edge,out2)
          call handle_error_pnetcdf('nf_put_vara_real: uab',status)
          out2 = real(vab,4)
          status=nfmpi_put_vara_real_all(ncid,vab_varid,start,edge,out2)
          call handle_error_pnetcdf('nf_put_vara_real: vab',status)
          out2 = real(elb,4)
          status=nfmpi_put_vara_real_all(ncid,elb_varid,start,edge,out2)
          call handle_error_pnetcdf('nf_put_vara_real: elb',status)
          out2 = real(wusurf,4)
          status=nfmpi_put_vara_real_all(ncid,wusurf_varid,start,edge
     &                                                            ,out2)
          call handle_error_pnetcdf('nf_put_vara_real: wusurf',status)
          out2 = real(wvsurf,4)
          status=nfmpi_put_vara_real_all(ncid,wvsurf_varid,start,edge
     &                                                            ,out2)
          call handle_error_pnetcdf('nf_put_vara_real: wvsurf',status)
          out2 = real(wtsurf,4)
          status=nfmpi_put_vara_real_all(ncid,wtsurf_varid,start,edge
     &                                                            ,out2)
          call handle_error_pnetcdf('nf_put_vara_real: wtsurf',status)
          out2 = real(swrad,4)
          status=nfmpi_put_vara_real_all(ncid,uhtfl_varid,start,edge
     &                                                            ,out2)
          call handle_error_pnetcdf('nf_put_vara_real: swrad',status)
          out2 = real(wssurf,4)
          status=nfmpi_put_vara_real_all(ncid,wssurf_varid,start,edge
     &                                                            ,out2)
          call handle_error_pnetcdf('nf_put_vara_real: wssurf',status)
          if ( calc_ice ) then
            out2 = real(ice,4)
            status=nfmpi_put_vara_real_all(ncid,icec_varid,start,edge
     &                                                            ,out2)
            call handle_error_pnetcdf('nf_put_vara_real: icec',status)
            out2 = real(ui,4)
            status=nfmpi_put_vara_real_all(ncid,ui_varid,start,edge
     &                                                            ,out2)
            call handle_error_pnetcdf('nf_put_vara_real: ui',status)
            out2 = real(vi,4)
            status=nfmpi_put_vara_real_all(ncid,vi_varid,start,edge
     &                                                            ,out2)
            call handle_error_pnetcdf('nf_put_vara_real: vi',status)
          end if

          start(1) = i_global(1)
          start(2) = j_global(1)
          start(3) = 1
          start(4) = iout
          edge(1) = im
          edge(2) = jm
          edge(3) = kb
          edge(4) = 1

          out3 = real(u,4)
          status=nfmpi_put_vara_real_all(ncid,u_varid,start,edge,out3)
          call handle_error_pnetcdf('nf_put_vara_real: u',status)
          out3 = real(v,4)
          status=nfmpi_put_vara_real_all(ncid,v_varid,start,edge,out3)
          call handle_error_pnetcdf('nf_put_vara_real: v',status)
          out3 = real(w,4)
          status=nfmpi_put_vara_real_all(ncid,w_varid,start,edge,out3)
          call handle_error_pnetcdf('nf_put_vara_real: w',status)
          out3 = real(t,4)
          status=nfmpi_put_vara_real_all(ncid,t_varid,start,edge,out3)
          call handle_error_pnetcdf('nf_put_vara_real: t',status)
          out3 = real(s,4)
          status=nfmpi_put_vara_real_all(ncid,s_varid,start,edge,out3)
          call handle_error_pnetcdf('nf_put_vara_real: s',status)
          out3 = real(rho,4)
          status=nfmpi_put_vara_real_all(ncid,rho_varid,start,edge,out3)
          call handle_error_pnetcdf('nf_put_vara_real: rho',status)
          out3 = real(kh,4)
          status=nfmpi_put_vara_real_all(ncid,kh_varid,start,edge,out3)
          call handle_error_pnetcdf('nf_put_vara_real: kh',status)
          out3 = real(km,4)
          status=nfmpi_put_vara_real_all(ncid,km_varid,start,edge,out3)
          call handle_error_pnetcdf('nf_put_vara_real: km',status)

        else !lyo:20110224:stcc:save *_mean

          start(1) = iout
          edge(1) = 1

          out1 = real(time,4)
          status=nfmpi_put_vara_real_all(ncid,time_varid,start,edge
     &                                                            ,out1)
          call handle_error_pnetcdf('nf_put_vara_real:time',status)

          start(1) = i_global(1)
          start(2) = j_global(1)
          start(3) = iout
          edge(1) = im
          edge(2) = jm
          edge(3) = 1

          out2 = real(uab_mean,4)
          status=nfmpi_put_vara_real_all(ncid,uab_varid,start,edge,out2)
          call handle_error_pnetcdf('nf_put_vara_real: uab',status)
          out2 = real(vab_mean,4)
          status=nfmpi_put_vara_real_all(ncid,vab_varid,start,edge,out2)
          call handle_error_pnetcdf('nf_put_vara_real: vab',status)
          out2 = real(elb_mean,4)
          status=nfmpi_put_vara_real_all(ncid,elb_varid,start,edge,out2)
          call handle_error_pnetcdf('nf_put_vara_real: elb',status)
          out2 = real(wusurf_mean,4)
          status=nfmpi_put_vara_real_all(ncid,wusurf_varid,start,edge
     &                                                            ,out2)
          call handle_error_pnetcdf('nf_put_vara_real: wusurf',status)
          out2 = real(wvsurf_mean,4)
          status=nfmpi_put_vara_real_all(ncid,wvsurf_varid,start,edge
     &                                                            ,out2)
          call handle_error_pnetcdf('nf_put_vara_real: wvsurf',status)
          out2 = real(wtsurf_mean,4)
          status=nfmpi_put_vara_real_all(ncid,wtsurf_varid,start,edge
     &                                                            ,out2)
          call handle_error_pnetcdf('nf_put_vara_real: wtsurf',status)
          out2 = real(swrad,4)
          status=nfmpi_put_vara_real_all(ncid,uhtfl_varid,start,edge
     &                                                            ,out2)
          call handle_error_pnetcdf('nf_put_vara_real: swrad',status)
          out2 = real(wssurf_mean,4)
          status=nfmpi_put_vara_real_all(ncid,wssurf_varid,start,edge
     &                                                            ,out2)
          call handle_error_pnetcdf('nf_put_vara_real: wssurf',status)
          if ( calc_ice ) then
            out2 = real(ice,4)
            status=nfmpi_put_vara_real_all(ncid,icec_varid,start,edge
     &                                                            ,out2)
            call handle_error_pnetcdf('nf_put_vara_real: icec',status)
            out2 = real(ui,4)
            status=nfmpi_put_vara_real_all(ncid,ui_varid,start,edge
     &                                                            ,out2)
            call handle_error_pnetcdf('nf_put_vara_real: ui',status)
            out2 = real(vi,4)
            status=nfmpi_put_vara_real_all(ncid,vi_varid,start,edge
     &                                                            ,out2)
            call handle_error_pnetcdf('nf_put_vara_real: vi',status)
          end if

          if ( mode /= 2 ) then
            start(1) = i_global(1)
            start(2) = j_global(1)
            start(3) = 1
            start(4) = iout
            edge(1) = im
            edge(2) = jm
            edge(3) = kb
            edge(4) = 1

            out3 = real(u_mean,4)
            status=nfmpi_put_vara_real_all(ncid,u_varid,start,edge,out3)
            call handle_error_pnetcdf('nf_put_vara_real: u',status)
            out3 = real(v_mean,4)
            status=nfmpi_put_vara_real_all(ncid,v_varid,start,edge,out3)
            call handle_error_pnetcdf('nf_put_vara_real: v',status)
            out3 = real(w_mean,4)
            status=nfmpi_put_vara_real_all(ncid,w_varid,start,edge,out3)
            call handle_error_pnetcdf('nf_put_vara_real: w',status)
            out3 = real(t_mean,4)
            status=nfmpi_put_vara_real_all(ncid,t_varid,start,edge,out3)
            call handle_error_pnetcdf('nf_put_vara_real: t',status)
            out3 = real(s_mean,4)
            status=nfmpi_put_vara_real_all(ncid,s_varid,start,edge,out3)
            call handle_error_pnetcdf('nf_put_vara_real: s',status)
            out3 = real(rho_mean,4)
            status=nfmpi_put_vara_real_all(ncid,rho_varid,start,edge
     &                                                            ,out3)
            call handle_error_pnetcdf('nf_put_vara_real: rho',status)
            out3 = real(kh_mean,4)
            status=nfmpi_put_vara_real_all(ncid,kh_varid,start,edge
     &                                                            ,out3)
            call handle_error_pnetcdf('nf_put_vara_real: kh',status)
            out3 = real(km_mean,4)
            status=nfmpi_put_vara_real_all(ncid,km_varid,start,edge
     &                                                            ,out3)
            call handle_error_pnetcdf('nf_put_vara_real: km',status)
          end if

      end if !if ( iout == 1  ) then... !lyo:20110224:stcc:

! close file
      status = nf90mpi_close(ncid)
      call handle_error_pnetcdf( 'nf_close: output', status )

      return
      end

!_______________________________________________________________________
      subroutine write_SURF_pnetcdf( netcdf_out_file )!fhx:20110131:new subr
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      character(len=*), intent(in) :: netcdf_out_file
      character(len=120) str_tmp!, netcdf_out_file
      integer time_dimid,x_dimid,y_dimid,z_dimid
      integer z_varid,zz_varid,
     $        east_e_varid,
     $        north_e_varid,
     $        h_varid,fsm_varid
      integer time_varid
      integer usrf_varid,vsrf_varid,elsrf_varid,uwsrf_varid,vwsrf_varid
      integer ncid,status
      integer vdims(4)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(4),edge(4)


      iouts = iouts + 1

      if ( iouts == 1  ) then

!     create netcdf file
         if(my_task.eq.0) write(*,'(/''writing file '',a)')
     $        trim(netcdf_out_file)
         status=nfmpi_create(pom_comm,trim(netcdf_out_file),
     $     nf_clobber+nf_64bit_offset,mpi_info_null,ncid) !lyo:pac10:
         call handle_error_pnetcdf( 'nf_create: '//netcdf_out_file
     &                            , status )

! define global attributes
      length=len(trim(title))
      status=nfmpi_put_att_text(ncid,nf_global,'title',length,
     $                          trim(title))
      call handle_error_pnetcdf('nf_put_att_text',status)

      str_tmp='output file'
      length=len(trim(str_tmp))
      status=nfmpi_put_att_text(ncid,nf_global,'description',length,
     $                          trim(str_tmp))
      call handle_error_pnetcdf('nf_put_att_text',status)

! define dimensions
      length=nf_unlimited
      status=nfmpi_def_dim(ncid,'time',length,time_dimid)
      call handle_error_pnetcdf('nf_def_dim: time',status)
      length=kb
      status=nfmpi_def_dim(ncid,'z',length,z_dimid)
      call handle_error_pnetcdf('nf_def_dim: z',status)
      length=jm_global
      status=nfmpi_def_dim(ncid,'y',length,y_dimid)
      call handle_error_pnetcdf('nf_def_dim: y',status)
      length=im_global
      status=nfmpi_def_dim(ncid,'x',length,x_dimid)
      call handle_error_pnetcdf('nf_def_dim: x',status)

! define variables and their attributes
      vdims(1)=time_dimid
      str_tmp='days since '//time_start
      call def_var_pnetcdf(ncid,'time',1,vdims,time_varid,nf_float,
     $                     'time',str_tmp,
     $                     -1,0.,' ',.false.)

      vdims(1)=z_dimid
      call def_var_pnetcdf(ncid,'z',1,vdims,z_varid,nf_float,
     $                     'sigma of cell face','sigma_level',
     $                     -1,0.,' ',.false.)
      length=22
      status=nfmpi_put_att_text(ncid,z_varid,'standard_name',
     $                          length,'ocean_sigma_coordinate')
      length=26
      status=nfmpi_put_att_text(ncid,z_varid,'formula_terms',
     $                          length,'sigma: z eta: elb depth: h')
      call handle_error_pnetcdf('nf_put_att_text: z',status)

      call def_var_pnetcdf(ncid,'zz',1,vdims,zz_varid,nf_float,
     $                     'sigma of cell centre','sigma_level',
     $                     -1,0.,' ',.false.)
      length=22
      status=nfmpi_put_att_text(ncid,zz_varid,'standard_name',
     $                          length,'ocean_sigma_coordinate')
      length=27
      status=nfmpi_put_att_text(ncid,zz_varid,'formula_terms',
     $                          length,'sigma: zz eta: elb depth: h')
      call handle_error_pnetcdf('nf_put_att_text: zz',status)

      vdims(1)=x_dimid
      vdims(2)=y_dimid

      call def_var_pnetcdf(ncid,'alon',2,vdims,east_e_varid,nf_float,
     $                     'easting of elevation points','degree',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'alat',2,vdims,north_e_varid,nf_float,
     $                     'northing of elevation points','degree',
     $                     -1,0.,'east_e north_e',.true.)

      call def_var_pnetcdf(ncid,'h',2,vdims,h_varid,nf_float,
     $                     'undisturbed water depth','metre',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'fsm',2,vdims,fsm_varid,nf_float,
     $                     'free surface mask','dimensionless',
     $                     -1,0.,'east_e north_e',.true.)


      vdims(1)=x_dimid
      vdims(2)=y_dimid
      vdims(3)=time_dimid
      call def_var_pnetcdf(ncid,'usrf',3,vdims,usrf_varid,nf_float,
     $                     'surface u','metre/sec',
     $                     -1,0.,'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'vsrf',3,vdims,vsrf_varid,nf_float,
     $                     'surface v','metre/sec',
     $                     -1,0.,'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'elsrf',3,vdims,elsrf_varid,nf_float,
     $                     'surface elevation','metre',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'uwsrf',3,vdims,uwsrf_varid,nf_float,
     $                     'wind u','metre/sec',
     $                     -1,0.,'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'vwsrf',3,vdims,vwsrf_varid,nf_float,
     $                     'wind v','metre/sec',
     $                     -1,0.,'east_v north_v',.true.)


! end definitions
      status=nfmpi_enddef(ncid)
      call handle_error_pnetcdf('nf_enddef: output',status)

! write data
!      start(1)=1
!      edge(1)=1
!      status=nfmpi_put_vara_real_all(ncid,time_varid,start,edge,time)
!      call handle_error_pnetcdf('nf_put_vara_real:time',status)

      start(1)=1
      edge(1)=kb
      status=nfmpi_put_vara_real_all(ncid,z_varid,start,edge,real(z))
      call handle_error_pnetcdf('nf_put_var_real: z',status)
      status=nfmpi_put_vara_real_all(ncid,zz_varid,start,edge,real(zz))
      call handle_error_pnetcdf('nf_put_var_real: zz',status)

      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm

      status=nfmpi_put_vara_real_all(ncid,east_e_varid,start,edge,
     $                               real(east_e))
      call handle_error_pnetcdf('nf_put_var_real: east_e',status)
      status=nfmpi_put_vara_real_all(ncid,north_e_varid,start,edge,
     $                               real(north_e))
      call handle_error_pnetcdf('nf_put_var_real: north_e',status)
      status=nfmpi_put_vara_real_all(ncid,h_varid,start,edge,real(h))
      call handle_error_pnetcdf('nf_put_var_real: h',status)
      status=nfmpi_put_vara_real_all(ncid,fsm_varid,start,edge
     $                                                       ,real(fsm))
      call handle_error_pnetcdf('nf_put_var_real: fsm',status)

      else

!     open netcdf file
         if(my_task.eq.0) write(*,'(/''writing file '',a)')
     $        trim(netcdf_out_file)
      status=nfmpi_open(pom_comm,trim(netcdf_out_file),nf_write,
     $     mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nf_open: '//netcdf_out_file
     &                         , status )

!!     inquire dimensions
!
!      status=nfmpi_inq_dimid(ncid,'time',time_dimid)
!      call handle_error_pnetcdf('nfmpi_inq_dim: time',status)
!      status=nfmpi_inq_dim(ncid,'z',z_dimid)
!      call handle_error_pnetcdf('nfmpi_inq_dim: z',status)
!      status=nfmpi_inq_dim(ncid,'y',y_dimid)
!      call handle_error_pnetcdf('nfmpi_inq_dim: y',status)
!      status=nfmpi_inq_dim(ncid,'x',x_dimid)
!     call handle_error_pnetcdf('nfmpi_inq_dim: x',status)



!     inqure variables
      status=nfmpi_inq_varid(ncid,'time',time_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: time'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'usrf',usrf_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: usrf'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'vsrf',vsrf_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: vsrf'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'elsrf',elsrf_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: elsrf'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'uwsrf',uwsrf_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: uwsrf'
     $                          ,status)
      status=nfmpi_inq_varid(ncid,'vwsrf',vwsrf_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: vwsrf'
     $                          ,status)

      endif

!     write data

      start(1)=iouts
      edge(1)=1
      status=nfmpi_put_vara_real_all(ncid,time_varid,start,edge
     $                                                      ,real(time))
      call handle_error_pnetcdf('nf_put_vara_real:time',status)

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=iouts
      edge(1)=im
      edge(2)=jm
      edge(3)=1
      status=nfmpi_put_vara_real_all(
     $      ncid,usrf_varid,start,edge,real(usrf_mean))
      call handle_error_pnetcdf(
     $     'nf_put_vara_real: usrf',status)
      status=nfmpi_put_vara_real_all(
     $     ncid,vsrf_varid,start,edge,real(vsrf_mean))
      call handle_error_pnetcdf(
     $     'nf_put_vara_real: vsrf',status)
      status=nfmpi_put_vara_real_all(
     $      ncid,elsrf_varid,start,edge,real(elsrf_mean))
      call handle_error_pnetcdf(
     $      'nf_put_vara_real: elsrf',status)
      status=nfmpi_put_vara_real_all(
     $     ncid,uwsrf_varid,start,edge,real(uwsrf_mean))
      call handle_error_pnetcdf(
     $     'nf_put_vara_real: uwsrf',status)
      status=nfmpi_put_vara_real_all(
     $     ncid,vwsrf_varid,start,edge,real(vwsrf_mean))
      call handle_error_pnetcdf(
     $     'nf_put_vara_real: vwsrf',status)


! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close: output',status)

      return
      end

!_______________________________________________________________________
!_______________________________________________________________________
      subroutine write_output_init_pnetcdf( netcdf_out_file ) !lyo:20110224:alu:stcc:!lyo:pac10:add this subroutine
!alu:20100614:  write output data for initial fields
!alu:20110105: change the name and clean up a little bit
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      character(len=*), intent(in) :: netcdf_out_file
      character(len=120) str_tmp !,netcdf_out_file !lyo:20110224:alu:stcc:
      integer time_dimid,x_dimid,y_dimid,z_dimid
      integer time_varid
     $       ,z_varid,zz_varid,dx_varid,dy_varid
     $       , east_c_varid, east_e_varid, east_u_varid, east_v_varid
     $       ,north_c_varid,north_e_varid,north_u_varid,north_v_varid
     $       ,h_varid,fsm_varid,dum_varid,dvm_varid,rot_varid
     $       ,wusurf_varid,wvsurf_varid,wtsurf_varid,wssurf_varid
     $       ,t_varid,s_varid,rho_varid,km_varid,kh_varid
     $       ,frz_varid
     $       ,aamfrz_varid
      integer ele_varid,elw_varid,els_varid,eln_varid
      integer ncid,status,flags
      integer vdims(4)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(4),edge(4)


!     netcdf_out_file='./out/out.init.nc' !lyo:20110224:alu:stcc:
      if(my_task.eq.0) write(*,'(/''writing file '',a)')
     $     trim(netcdf_out_file)

      flags= NF_CLOBBER + NF_64BIT_OFFSET  !for large file( > 2GB) !lyo:pac10:
!     flags= NF_CLOBBER                     !for small file( < 2GB)

      status=nfmpi_create(pom_comm,trim(netcdf_out_file),flags,
     $     mpi_info_null,ncid)
      call handle_error_pnetcdf( 'nf_create: '//netcdf_out_file
     &                         , status )

! define global attributes
      length=len(trim(title))
      status=nfmpi_put_att_text(ncid,nf_global,'title',length,
     $                          trim(title))
      call handle_error_pnetcdf('nf_put_att_text',status)

      str_tmp='output file'
      length=len(trim(str_tmp))
      status=nfmpi_put_att_text(ncid,nf_global,'description',length,
     $                          trim(str_tmp))
      call handle_error_pnetcdf('nf_put_att_text',status)

! define dimensions
      length=1
      status=nfmpi_def_dim(ncid,'time',length,time_dimid)
      call handle_error_pnetcdf('nf_def_dim: time',status)
      length=kb
      status=nfmpi_def_dim(ncid,'z',length,z_dimid)
      call handle_error_pnetcdf('nf_def_dim: z',status)
      length=jm_global
      status=nfmpi_def_dim(ncid,'y',length,y_dimid)
      call handle_error_pnetcdf('nf_def_dim: y',status)
      length=im_global
      status=nfmpi_def_dim(ncid,'x',length,x_dimid)
      call handle_error_pnetcdf('nf_def_dim: x',status)

! define variables and their attributes
      vdims(1)=time_dimid
      str_tmp='days since '//time_start
      call def_var_pnetcdf(ncid,'time',1,vdims,time_varid,nf_float,
     $                     'time',str_tmp,
     $                     -1,0.,' ',.false.)

      vdims(1)=z_dimid
      call def_var_pnetcdf(ncid,'sigma',1,vdims,z_varid,nf_float,
     $                     'sigma of cell face','sigma_level',
     $                     -1,0.,' ',.false.)
      str_tmp='ocean_sigma_coordinate'
      length=len(trim(str_tmp))
      status=nfmpi_put_att_text(ncid,z_varid,'standard_name',
     $                          length,str_tmp)
      str_tmp='sigma: z eta: elb depth: h'
      length=len(trim(str_tmp))
      status=nfmpi_put_att_text(ncid,z_varid,'formula_terms',
     $                          length,str_tmp)
      call handle_error_pnetcdf('nf_put_att_text: z',status)

      call def_var_pnetcdf(ncid,'zz',1,vdims,zz_varid,nf_float,
     $                     'sigma of cell centre','sigma_level',
     $                     -1,0.,' ',.false.)
      str_tmp='ocean_sigma_coordinate'
      length=len(trim(str_tmp))
      status=nfmpi_put_att_text(ncid,zz_varid,'standard_name',
     $                          length,str_tmp)
      str_tmp='sigma: zz eta: elb depth: h'
      length=len(trim(str_tmp))
      status=nfmpi_put_att_text(ncid,zz_varid,'formula_terms',
     $                          length,str_tmp)
      call handle_error_pnetcdf('nf_put_att_text: zz',status)

!alu:add calculated elebc output
      vdims(1)=x_dimid
      call def_var_pnetcdf(ncid,'els',1,vdims,els_varid,nf_float, !South
     $                      'South elebc','metre',
     $                      -1,0.,'east_e',.true.)
      call def_var_pnetcdf(ncid,'eln',1,vdims,eln_varid,nf_float, !North
     $                      'Nouth elebc','metre',
     $                      -1,0.,'east_e',.true.)
      vdims(1)=y_dimid
      call def_var_pnetcdf(ncid,'ele',1,vdims,ele_varid,nf_float, !East
     $                      'East elebc','metre',
     $                      -1,0.,'north_e',.true.)
      call def_var_pnetcdf(ncid,'elw',1,vdims,elw_varid,nf_float, !West
     $                      'West elebc','metre',
     $                      -1,0.,'north_e',.true.)


      vdims(1)=x_dimid
      vdims(2)=y_dimid
      call def_var_pnetcdf(ncid,'dx',2,vdims,dx_varid,nf_float,
     $                     'grid increment in x','metre',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'dy',2,vdims,dy_varid,nf_float,
     $                     'grid increment in y','metre',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'east_u',2,vdims,east_u_varid,nf_float,
     $                     'easting of u-points','degree',
     $                     -1,0.,'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'east_v',2,vdims,east_v_varid,nf_float,
     $                     'easting of v-points','degree',
     $                     -1,0.,'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'east_e',2,vdims,east_e_varid,nf_float,
     $                     'easting of elevation points','degree',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'east_c',2,vdims,east_c_varid,nf_float,
     $                     'easting of cell corners','degree',
     $                     -1,0.,'east_c north_c',.true.)
      call def_var_pnetcdf(ncid,'north_u',2,vdims,north_u_varid,nf_float
     $                    ,'northing of u-points','degree',
     $                     -1,0.,'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'north_v',2,vdims,north_v_varid,nf_float
     $                    ,'northing of v-points','degree',
     $                     -1,0.,'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'north_e',2,vdims,north_e_varid,nf_float
     $                    ,'northing of elevation points','degree',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'north_c',2,vdims,north_c_varid,nf_float
     $                    ,'northing of cell corners','degree',
     $                     -1,0.,'east_c north_c',.true.)
      call def_var_pnetcdf(ncid,'rot',2,vdims,rot_varid,nf_float,
     $                    'Rotation angle of x-axis wrt. east','degree',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'h',2,vdims,h_varid,nf_float,
     $                     'undisturbed water depth','metre',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'fsm',2,vdims,fsm_varid,nf_float,
     $                     'free surface mask','dimensionless',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'dum',2,vdims,dum_varid,nf_float,
     $                     'u-velocity mask','dimensionless',
     $                     -1,0.,'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'dvm',2,vdims,dvm_varid,nf_float,
     $                     'v-velocity mask','dimensionless',
     $                     -1,0.,'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'frz',2,vdims,frz_varid,nf_float,
     $                     'frz coeffi','dimensionless',
     $                     -1,0.,'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'aamfrz',2,vdims,aamfrz_varid,nf_float,
     $                     'aamfrz coeffi','dimensionless',
     $                     -1,0.,'east_v north_v',.true.)

      vdims(1)=x_dimid
      vdims(2)=y_dimid
      vdims(3)=time_dimid
      call def_var_pnetcdf(ncid,'wusurf',3,vdims,wusurf_varid,nf_float,
     $                     'x-momentum flux','metre^2/sec^2',
     $                     -1,0.,'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'wvsurf',3,vdims,wvsurf_varid,nf_float,
     $                     'y-momentum flux','metre^2/sec^2',
     $                     -1,0.,'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'wtsurf',3,vdims,wtsurf_varid,nf_float,
     $                     'temperature flux','deg m/s',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'wssurf',3,vdims,wssurf_varid,nf_float,
     $                     'salinity flux','psu m/s',
     $                     -1,0.,'east_e north_e',.true.)

      vdims(1)=x_dimid
      vdims(2)=y_dimid
      vdims(3)=z_dimid
      vdims(4)=time_dimid
      call def_var_pnetcdf(ncid,'t',4,vdims,t_varid,nf_float,
     $                     'potential temperature','K',
     $                     -1,0.,'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'s',4,vdims,s_varid,nf_float,
     $                     'salinity x rho / rhoref','PSS',
     $                     -1,0.,'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'rho',4,vdims,rho_varid,nf_float,
     $                     '(density-1000)/rhoref','dimensionless',
     $                     -1,0.,'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'kh',4,vdims,kh_varid,nf_float,
     $                     'vertical diffusivity','metre2/sec',
     $                     -1,0.,'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'km',4,vdims,km_varid,nf_float,
     $                     'vertical viscosity','metre2/sec',
     $                     -1,0.,'east_e north_e zz',.true.)

! end definitions
      status=nfmpi_enddef(ncid)
      call handle_error_pnetcdf('nf_enddef: output',status)


! write data
      start(1)=1
      edge(1)=1
      status=nfmpi_put_vara_real_all(ncid,time_varid,start,edge
     &                                                    ,real(time,4))
      call handle_error_pnetcdf('nf_put_vara_real:time',status)

      start(1)=1
      edge(1)=kb
      status=nfmpi_put_vara_real_all(ncid,z_varid,start,edge
     &                                                       ,real(z,4))
      call handle_error_pnetcdf('nf_put_var_real: z',status)
      status=nfmpi_put_vara_real_all(ncid,zz_varid,start,edge
     &                                                      ,real(zz,4))
      call handle_error_pnetcdf('nf_put_var_real: zz',status)

!alu:add elebc
!--East
       start(1)=j_global(1)
       edge(1)=jm
       status=nfmpi_put_vara_real_all(ncid,ele_varid,start,edge
     &                                                     ,real(ele,4))
       call handle_error_pnetcdf('nf_put_var_real: ele',status)
!--West
       start(1)=j_global(1)
       edge(1)=jm
       status=nfmpi_put_vara_real_all(ncid,elw_varid,start,edge
     &                                                     ,real(elw,4))
       call handle_error_pnetcdf('nf_put_var_real: elw',status)
!--South
       start(1)=i_global(1)
       edge(1)=im
       status=nfmpi_put_vara_real_all(ncid,els_varid,start,edge
     &                                                     ,real(els,4))
       call handle_error_pnetcdf('nf_put_var_real: els',status)
!--North
       start(1)=i_global(1)
       edge(1)=im
       status=nfmpi_put_vara_real_all(ncid,eln_varid,start,edge
     &                                                     ,real(eln,4))
       call handle_error_pnetcdf('nf_put_var_real: eln',status)



      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      status=nfmpi_put_vara_real_all(ncid,dx_varid,start,edge
     &                                                      ,real(dx,4))
      call handle_error_pnetcdf('nf_put_var_real: dx',status)
      status=nfmpi_put_vara_real_all(ncid,dy_varid,start,edge
     &                                                      ,real(dy,4))
      call handle_error_pnetcdf('nf_put_var_real: dy',status)
      status=nfmpi_put_vara_real_all(ncid,east_u_varid,start,edge
     &                                                  ,real(east_u,4))
      call handle_error_pnetcdf('nf_put_var_real: east_u',status)
      status=nfmpi_put_vara_real_all(ncid,east_v_varid,start,edge
     &                                                  ,real(east_v,4))
      call handle_error_pnetcdf('nf_put_var_real: east_v',status)
      status=nfmpi_put_vara_real_all(ncid,east_e_varid,start,edge
     &                                                  ,real(east_e,4))
      call handle_error_pnetcdf('nf_put_var_real: east_e',status)
      status=nfmpi_put_vara_real_all(ncid,east_c_varid,start,edge
     &                                                  ,real(east_c,4))
      call handle_error_pnetcdf('nf_put_var_real: east_c',status)
      status=nfmpi_put_vara_real_all(ncid,north_u_varid,start,edge
     &                                                 ,real(north_u,4))
      call handle_error_pnetcdf('nf_put_var_real: north_u',status)
      status=nfmpi_put_vara_real_all(ncid,north_v_varid,start,edge
     &                                                 ,real(north_v,4))
      call handle_error_pnetcdf('nf_put_var_real: north_v',status)
      status=nfmpi_put_vara_real_all(ncid,north_e_varid,start,edge
     &                                                 ,real(north_e,4))
      call handle_error_pnetcdf('nf_put_var_real: north_e',status)
      status=nfmpi_put_vara_real_all(ncid,north_c_varid,start,edge
     &                                                 ,real(north_c,4))
      call handle_error_pnetcdf('nf_put_var_real: north_c',status)
      status=nfmpi_put_vara_real_all(ncid,rot_varid,start,edge
     &                                                     ,real(rot,4))
      call handle_error_pnetcdf('nf_put_var_real: rot',status)
      status=nfmpi_put_vara_real_all(ncid,h_varid,start,edge,real(h,4))
      call handle_error_pnetcdf('nf_put_var_real: h',status)
      status=nfmpi_put_vara_real_all(ncid,fsm_varid,start,edge
     &                                                     ,real(fsm,4))
      call handle_error_pnetcdf('nf_put_var_real: fsm',status)
      status=nfmpi_put_vara_real_all(ncid,dum_varid,start,edge
     &                                                     ,real(dum,4))
      call handle_error_pnetcdf('nf_put_var_real: dum',status)
      status=nfmpi_put_vara_real_all(ncid,dvm_varid,start,edge
     &                                                     ,real(dvm,4))
      call handle_error_pnetcdf('nf_put_var_real: dvm',status)
      status=nfmpi_put_vara_real_all(ncid,frz_varid,start,edge
     &                                                     ,real(frz,4))
      call handle_error_pnetcdf('nf_put_var_real: frz',status)
      status=nfmpi_put_vara_real_all(ncid,aamfrz_varid,start,edge
     &                                                  ,real(aamfrz,4))
      call handle_error_pnetcdf('nf_put_var_real: aamfrz',status)

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=1
      status=nfmpi_put_vara_real_all(
     $     ncid,wusurf_varid,start,edge,real(wusurf,4))
      call handle_error_pnetcdf('nf_put_vara_real: wusurf',status)
      status=nfmpi_put_vara_real_all(
     $     ncid,wvsurf_varid,start,edge,real(wvsurf,4))
      call handle_error_pnetcdf('nf_put_vara_real: wvsurf',status)
      status=nfmpi_put_vara_real_all(
     $     ncid,wtsurf_varid,start,edge,real(wtsurf,4))
      call handle_error_pnetcdf('nf_put_vara_real: wtsurf',status)
      status=nfmpi_put_vara_real_all(
     $     ncid,wssurf_varid,start,edge,real(wssurf,4))
      call handle_error_pnetcdf('nf_put_vara_real: wssurf',status)


      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      start(4)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=kb
      edge(4)=1
      status=nfmpi_put_vara_real_all(ncid,t_varid,start,edge,real(tb,4))
      call handle_error_pnetcdf('nf_put_vara_real: t',status)
      status=nfmpi_put_vara_real_all(ncid,s_varid,start,edge,real(sb,4))
      call handle_error_pnetcdf('nf_put_vara_real: s',status)
      status=nfmpi_put_vara_real_all(ncid,rho_varid,start,edge
     &                                                     ,real(rho,4))
      call handle_error_pnetcdf('nf_put_vara_real: rho',status)
      status=nfmpi_put_vara_real_all(ncid,kh_varid,start,edge
     &                                                      ,real(kh,4))
      call handle_error_pnetcdf('nf_put_vara_real: kh',status)
      status=nfmpi_put_vara_real_all(ncid,km_varid,start,edge
     &                                                      ,real(km,4))
      call handle_error_pnetcdf('nf_put_vara_real: km',status)


! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close: output: '
     $                           ,status)

      return
      end
!_______________________________________________________________________
      subroutine write_debug_pnetcdf( netcdf_out_file, flu, flv )
        implicit none
# include "mpif.h"
# include "pnetcdf.inc"
        include 'pom.h'
        character(len=*), intent(in) :: netcdf_out_file
        character(len=120) str_tmp !,netcdf_out_file !lyo:20110224:alu:stcc:
        integer time_dimid,x_dimid,y_dimid,z_dimid
        integer time_varid,z_varid,zz_varid,icec_varid,ui_varid,vi_varid
        integer uabe_varid,uabw_varid,vabs_varid,vabn_varid,hi_varid
        integer tbe_varid,tbw_varid,tbs_varid,tbn_varid,uf_varid
        integer sbe_varid,sbw_varid,sbs_varid,sbn_varid,vf_varid
        integer km_varid,kh_varid
        integer ncid,status,flags
        integer vdims(4)
        integer(mpi_offset_kind) length
        integer(mpi_offset_kind) start(4),edge(4)
        real(kind=rk), dimension(im,jm), intent(in) :: flu, flv

        if(my_task.eq.0) write(*,'(/''writing file '',a)')
     $     trim(netcdf_out_file)

        flags= NF_CLOBBER + NF_64BIT_OFFSET

        status=nfmpi_create(pom_comm,trim(netcdf_out_file),flags,
     $     mpi_info_null,ncid)
        call handle_error_pnetcdf( 'nf_create: '//netcdf_out_file
     &                           , status )

! define global attributes
        length=len(trim(title))
        status=nfmpi_put_att_text(ncid,nf_global,'title',length,
     $                          trim(title))
        call handle_error_pnetcdf('nf_put_att_text',status)

        str_tmp='output file'
        length=len(trim(str_tmp))
        status=nfmpi_put_att_text(ncid,nf_global,'description',length,
     $                          trim(str_tmp))
        call handle_error_pnetcdf('nf_put_att_text',status)

! define dimensions
        length=1
        status=nfmpi_def_dim(ncid,'time',length,time_dimid)
        call handle_error_pnetcdf('nf_def_dim: time',status)
        length=kb
        status=nfmpi_def_dim(ncid,'z',length,z_dimid)
        call handle_error_pnetcdf('nf_def_dim: z',status)
        length=jm_global
        status=nfmpi_def_dim(ncid,'y',length,y_dimid)
        call handle_error_pnetcdf('nf_def_dim: y',status)
        length=im_global
        status=nfmpi_def_dim(ncid,'x',length,x_dimid)
        call handle_error_pnetcdf('nf_def_dim: x',status)

! define variables and their attributes
        vdims(1)=time_dimid
        str_tmp='days since '//time_start
        call def_var_pnetcdf(ncid,'time',1,vdims,time_varid,nf_float,
     $                     'time',str_tmp,
     $                     -1,0.,' ',.false.)

        vdims(1)=z_dimid

!alu:add calculated elebc output
        vdims(1)=x_dimid
        call def_var_pnetcdf(ncid,'vabs',1,vdims,vabs_varid,nf_float, !South
     $                      'South barotropic boundary velocity','m/s',
     $                      -1,0.,'east_e',.true.)
        call def_var_pnetcdf(ncid,'vabn',1,vdims,vabn_varid,nf_float, !North
     $                      'Nouth barotropic boundary velocity','m/s',
     $                      -1,0.,'east_e',.true.)
        vdims(1)=y_dimid
        call def_var_pnetcdf(ncid,'uabe',1,vdims,uabe_varid,nf_float, !East
     $                      'East barotropic boundary velocity','m/s',
     $                      -1,0.,'north_e',.true.)
        call def_var_pnetcdf(ncid,'uabw',1,vdims,uabw_varid,nf_float, !West
     $                      'West barotropic boundary velocity','m/s',
     $                      -1,0.,'north_e',.true.)


      vdims(1)=y_dimid
      vdims(2)=z_dimid
      call def_var_pnetcdf(ncid,'tbe',2,vdims,tbe_varid,nf_float,
     $                     'east boundary temp','degree',
     $                     -1,0.,'north_e sigma',.true.)
      call def_var_pnetcdf(ncid,'sbe',2,vdims,sbe_varid,nf_float,
     $                     'east boundary salt','pps',
     $                     -1,0.,'north_e sigma',.true.)
      call def_var_pnetcdf(ncid,'tbw',2,vdims,tbw_varid,nf_float,
     $                     'west boundary temp','degree',
     $                     -1,0.,'north_e sigma',.true.)
      call def_var_pnetcdf(ncid,'sbw',2,vdims,sbw_varid,nf_float,
     $                     'west boundary salt','pps',
     $                     -1,0.,'north_e sigma',.true.)
      vdims(1)=x_dimid
      call def_var_pnetcdf(ncid,'tbs',2,vdims,tbs_varid,nf_float,
     $                     'south boundary temp','degree',
     $                     -1,0.,'east_e sigma',.true.)
      call def_var_pnetcdf(ncid,'sbs',2,vdims,sbs_varid,nf_float,
     $                     'south boundary salt','pps',
     $                     -1,0.,'east_e sigma',.true.)
      call def_var_pnetcdf(ncid,'tbn',2,vdims,tbn_varid,nf_float
     $                    ,'north boundary temp','degree',
     $                     -1,0.,'east_e sigma',.true.)
      call def_var_pnetcdf(ncid,'sbn',2,vdims,sbn_varid,nf_float
     $                    ,'north boundary salt','pps',
     $                     -1,0.,'east_v north_v',.true.)
!      call def_var_pnetcdf(ncid,'north_e',2,vdims,north_e_varid,nf_float
!     $                    ,'northing of elevation points','degree',
!     $                     'east_e north_e',.true.)
!      call def_var_pnetcdf(ncid,'north_c',2,vdims,north_c_varid,nf_float
!     $                    ,'northing of cell corners','degree',
!     $                     'east_c north_c',.true.)
!      call def_var_pnetcdf(ncid,'rot',2,vdims,rot_varid,nf_float,
!     $                    'Rotation angle of x-axis wrt. east','degree',
!     $                     'east_e north_e',.true.)
!      call def_var_pnetcdf(ncid,'h',2,vdims,h_varid,nf_float,
!     $                     'undisturbed water depth','metre',
!     $                     'east_e north_e',.true.)
!      call def_var_pnetcdf(ncid,'fsm',2,vdims,fsm_varid,nf_float,
!     $                     'free surface mask','dimensionless',
!     $                     'east_e north_e',.true.)
!      call def_var_pnetcdf(ncid,'dum',2,vdims,dum_varid,nf_float,
!     $                     'u-velocity mask','dimensionless',
!     $                     'east_u north_u',.true.)
!      call def_var_pnetcdf(ncid,'dvm',2,vdims,dvm_varid,nf_float,
!     $                     'v-velocity mask','dimensionless',
!     $                     'east_v north_v',.true.)
!      call def_var_pnetcdf(ncid,'frz',2,vdims,frz_varid,nf_float,
!     $                     'frz coeffi','dimensionless',
!     $                     'east_v north_v',.true.)
!      call def_var_pnetcdf(ncid,'aamfrz',2,vdims,aamfrz_varid,nf_float,
!     $                     'aamfrz coeffi','dimensionless',
!     $                     'east_v north_v',.true.)
!
      vdims(1)=x_dimid
      vdims(2)=y_dimid
      vdims(3)=time_dimid
      call def_var_pnetcdf(ncid,'ice_c',3,vdims,icec_varid,nf_float,
     $                     'Sea ice concentration','',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'ice_u',3,vdims,ui_varid,nf_float,
     $                     'Sea ice u-velocity','metre/sec',
     $                     -1,0.,'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'ice_v',3,vdims,vi_varid,nf_float,
     $                     'Sea ice v-velocity','metre/sec',
     $                     -1,0.,'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'hi',3,vdims,hi_varid,nf_float,
     $                     'sea ice thickness','metre',
     $                     -1,0.,'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'flx_u',3,vdims,uf_varid,nf_float,
     $                     'Sea ice u-flux','conc/sec',
     $                     -1,0.,'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'flx_v',3,vdims,vf_varid,nf_float,
     $                     'Sea ice v-flux','conc/sec',
     $                     -1,0.,'east_v north_v',.true.)
!
      vdims(1)=x_dimid
      vdims(2)=y_dimid
      vdims(3)=z_dimid
      vdims(4)=time_dimid
!      call def_var_pnetcdf(ncid,'t',4,vdims,t_varid,nf_float,
!     $                     'potential temperature','K',
!     $                     'east_e north_e zz',.true.)
!      call def_var_pnetcdf(ncid,'s',4,vdims,s_varid,nf_float,
!     $                     'salinity x rho / rhoref','PSS',
!     $                     'east_e north_e zz',.true.)
!      call def_var_pnetcdf(ncid,'rho',4,vdims,rho_varid,nf_float,
!     $                     '(density-1000)/rhoref','dimensionless',
!     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'drhox',4,vdims,kh_varid,nf_float,
     $                     'vertical diffusivity','metre2/sec',
     $                     -1,0.,'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'drhoy',4,vdims,km_varid,nf_float,
     $                     'vertical viscosity','metre2/sec',
     $                     -1,0.,'east_e north_e zz',.true.)

! end definitions
      status=nfmpi_enddef(ncid)
      call handle_error_pnetcdf('nf_enddef: output',status)


! write data
      start(1)=1
      edge(1)=1
      status=nfmpi_put_vara_real_all(ncid,time_varid,start,edge
     &                                                    ,real(time,4))
      call handle_error_pnetcdf('nf_put_vara_real:time',status)

!      start(1)=1
!      edge(1)=kb
!      status=nfmpi_put_vara_real_all(ncid,z_varid,start,edge
!     &                                                       ,real(z,4))
!      call handle_error_pnetcdf('nf_put_var_real: z',status)
!      status=nfmpi_put_vara_real_all(ncid,zz_varid,start,edge
!     &                                                      ,real(zz,4))
!      call handle_error_pnetcdf('nf_put_var_real: zz',status)

!alu:add elebc
!--East
       start(1)=j_global(1)
       edge(1)=jm
       status=nfmpi_put_vara_real_all(ncid,uabe_varid,start,edge
     &                                                    ,real(uabe,4))
       call handle_error_pnetcdf('nf_put_var_real:uabe',status)
!--West
       start(1)=j_global(1)
       edge(1)=jm
       status=nfmpi_put_vara_real_all(ncid,uabw_varid,start,edge
     &                                                    ,real(uabw,4))
       call handle_error_pnetcdf('nf_put_var_real:uabw',status)
!--South
       start(1)=i_global(1)
       edge(1)=im
       status=nfmpi_put_vara_real_all(ncid,vabs_varid,start,edge
     &                                                    ,real(vabs,4))
       call handle_error_pnetcdf('nf_put_var_real:vabs',status)
!--North
       start(1)=i_global(1)
       edge(1)=im
       status=nfmpi_put_vara_real_all(ncid,vabn_varid,start,edge
     &                                                    ,real(vabn,4))
       call handle_error_pnetcdf('nf_put_var_real:vabn',status)



      start(1)=j_global(1)
      start(2)=1
      edge(1)=jm
      edge(2)=kb
      status=nfmpi_put_vara_real_all(ncid,tbe_varid,start,edge
     &                                                     ,real(tbe,4))
      call handle_error_pnetcdf('nf_put_var_real: tbe',status)
      status=nfmpi_put_vara_real_all(ncid,sbe_varid,start,edge
     &                                                     ,real(sbe,4))
      call handle_error_pnetcdf('nf_put_var_real: sbe',status)
      status=nfmpi_put_vara_real_all(ncid,tbw_varid,start,edge
     &                                                     ,real(tbw,4))
      call handle_error_pnetcdf('nf_put_var_real: tbw',status)
      status=nfmpi_put_vara_real_all(ncid,sbw_varid,start,edge
     &                                                     ,real(sbw,4))
      call handle_error_pnetcdf('nf_put_var_real: sbw',status)
      start(1)=i_global(1)
      edge(1)=im
      status=nfmpi_put_vara_real_all(ncid,tbs_varid,start,edge
     &                                                     ,real(tbs,4))
      call handle_error_pnetcdf('nf_put_var_real: tbs',status)
      status=nfmpi_put_vara_real_all(ncid,sbs_varid,start,edge
     &                                                     ,real(sbs,4))
      call handle_error_pnetcdf('nf_put_var_real: sbs',status)
      status=nfmpi_put_vara_real_all(ncid,tbn_varid,start,edge
     &                                                     ,real(tbn,4))
      call handle_error_pnetcdf('nf_put_var_real: tbn',status)
      status=nfmpi_put_vara_real_all(ncid,sbn_varid,start,edge
     &                                                     ,real(sbn,4))
      call handle_error_pnetcdf('nf_put_var_real: sbn',status)
!      status=nfmpi_put_vara_real_all(ncid,north_e_varid,start,edge
!     &                                                 ,real(north_e,4))
!      call handle_error_pnetcdf('nf_put_var_real: north_e',status
!     &                                                        ,nf_noerr)
!      status=nfmpi_put_vara_real_all(ncid,north_c_varid,start,edge
!     &                                                 ,real(north_c,4))
!      call handle_error_pnetcdf('nf_put_var_real: north_c',status
!     &                                                        ,nf_noerr)
!      status=nfmpi_put_vara_real_all(ncid,rot_varid,start,edge
!     &                                                     ,real(rot,4))
!      call handle_error_pnetcdf('nf_put_var_real: rot',status)
!      status=nfmpi_put_vara_real_all(ncid,h_varid,start,edge,real(h,4))
!      call handle_error_pnetcdf('nf_put_var_real: h',status)
!      status=nfmpi_put_vara_real_all(ncid,fsm_varid,start,edge
!     &                                                     ,real(fsm,4))
!      call handle_error_pnetcdf('nf_put_var_real: fsm',status)
!      status=nfmpi_put_vara_real_all(ncid,dum_varid,start,edge
!     &                                                     ,real(dum,4))
!      call handle_error_pnetcdf('nf_put_var_real: dum',status)
!      status=nfmpi_put_vara_real_all(ncid,dvm_varid,start,edge
!     &                                                     ,real(dvm,4))
!      call handle_error_pnetcdf('nf_put_var_real: dvm',status)
!      status=nfmpi_put_vara_real_all(ncid,frz_varid,start,edge
!     &                                                     ,real(frz,4))
!      call handle_error_pnetcdf('nf_put_var_real: frz',status)
!      status=nfmpi_put_vara_real_all(ncid,aamfrz_varid,start,edge
!     &                                                  ,real(aamfrz,4))
!      call handle_error_pnetcdf('nf_put_var_real: aamfrz',status
!     &                                                        ,nf_noerr)
!
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=1
      status=nfmpi_put_vara_real_all(
     $     ncid,icec_varid,start,edge,real(icb,4))
      call handle_error_pnetcdf(
     $     'nf_put_vara_real: icec',status)
      status=nfmpi_put_vara_real_all(
     $     ncid,ui_varid,start,edge,real(ui,4))
      call handle_error_pnetcdf(
     $     'nf_put_vara_real: iceu',status)
      status=nfmpi_put_vara_real_all(
     $     ncid,vi_varid,start,edge,real(vi,4))
      call handle_error_pnetcdf(
     $     'nf_put_vara_real: icev',status)
      status=nfmpi_put_vara_real_all(
     $     ncid,hi_varid,start,edge,real(hi,4))
      call handle_error_pnetcdf(
     $     'nf_put_vara_real: hi',status)
      status=nfmpi_put_vara_real_all(
     $     ncid,uf_varid,start,edge,real(flu,4))
      call handle_error_pnetcdf(
     $     'nf_put_vara_real: flxu',status)
      status=nfmpi_put_vara_real_all(
     $     ncid,vf_varid,start,edge,real(flv,4))
      call handle_error_pnetcdf(
     $     'nf_put_vara_real: flxv',status)
!
!
      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      start(4)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=kb
      edge(4)=1
!      status=nfmpi_put_vara_real_all(ncid,t_varid,start,edge,real(tb,4))
!      call handle_error_pnetcdf('nf_put_vara_real: t',status)
!      status=nfmpi_put_vara_real_all(ncid,s_varid,start,edge,real(sb,4))
!      call handle_error_pnetcdf('nf_put_vara_real: s',status)
!      status=nfmpi_put_vara_real_all(ncid,rho_varid,start,edge
!     &                                                     ,real(rho,4))
!      call handle_error_pnetcdf('nf_put_vara_real: rho',status)
      status=nfmpi_put_vara_real_all(ncid,kh_varid,start,edge
     &                                                   ,real(drhox,4))
      call handle_error_pnetcdf('nf_put_vara_real: kh',status)
      status=nfmpi_put_vara_real_all(ncid,km_varid,start,edge
     &                                                   ,real(drhoy,4))
      call handle_error_pnetcdf('nf_put_vara_real: km',status)


! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close: output: '
     $                           ,status)

      return
      end
!_______________________________________________________________________
!      subroutine write_sflx( netcdf_out_file, ta,sh,ra,ts,pr,cl )
!
!        use mpi
!        use pnetcdf
!
!        implicit none
!
!        include 'pom.h'
!        character(len=*), intent(in) :: netcdf_out_file
!        character(len=120) str_tmp
!        integer time_dimid,x_dimid,y_dimid
!        integer time_varid,tair_varid,rhum_varid,rain_varid
!        integer tskn_varid,pres_varid,tcld_varid
!        integer ncid,status,flags
!        integer vdims(4)
!        integer(mpi_offset_kind) length
!        integer(mpi_offset_kind) start(4),edge(4)
!        real(kind=rk),dimension(im,jm),intent(inout)::ta,sh,ra,ts,pr,cl
!        real(kind=rk),dimension(1) :: tarr
!
!!        if(my_task.eq.0) write(*,'(/''writing file '',a)')
!!     $     trim(netcdf_out_file)
!
!        flags= nf90_clobber + nf90_64bit_offset
!
!        status=nf90mpi_create(pom_comm,trim(netcdf_out_file),flags,
!     $     mpi_info_null,ncid)
!        call handle_error_pnetcdf('nf_create: '//netcdf_out_file,
!     $                          status,nf_noerr)
!
!! define global attributes
!        status=nf90mpi_put_att(ncid,nf_global,'title',trim(title))
!        call handle_error_pnetcdf('nf_put_att',status)
!
!        str_tmp='output file'
!        status=nf90mpi_put_att(ncid,nf_global,'description',
!     $                          trim(str_tmp))
!        call handle_error_pnetcdf('nf_put_att_text',status)
!
!! define dimensions
!        length=1
!        status=nf90mpi_def_dim(ncid,'time',length,time_dimid)
!        call handle_error_pnetcdf('nf_def_dim: time',status)
!        length=jm_global
!        status=nf90mpi_def_dim(ncid,'y',length,y_dimid)
!        call handle_error_pnetcdf('nf_def_dim: y',status)
!        length=im_global
!        status=nf90mpi_def_dim(ncid,'x',length,x_dimid)
!        call handle_error_pnetcdf('nf_def_dim: x',status)
!
!! define variables and their attributes
!        vdims(1)=time_dimid
!        str_tmp='days since '//time_start
!        call def_var_pnetcdf(ncid,'time',1,vdims,time_varid,nf90_double,
!     $                     'time',str_tmp,
!     $                     ' ',.false.)
!
!      vdims(1)=x_dimid
!      vdims(2)=y_dimid
!      vdims(3)=time_dimid
!      status=nf90mpi_def_var(ncid,'tair',nf90_double,vdims(1:3)
!     &                      ,tair_varid)
!      status=nf90mpi_def_var(ncid,'rhum',nf90_double,vdims(1:3)
!     &                      ,rhum_varid)
!      status=nf90mpi_def_var(ncid,'rain',nf90_double,vdims(1:3)
!     &                      ,rain_varid)
!      status=nf90mpi_def_var(ncid,'tskn',nf90_double,vdims(1:3)
!     &                      ,tskn_varid)
!      status=nf90mpi_def_var(ncid,'pres',nf90_double,vdims(1:3)
!     &                      ,pres_varid)
!      status=nf90mpi_def_var(ncid,'tcld',nf90_double,vdims(1:3)
!     &                      ,tcld_varid)
!! end definitions
!      status=nf90mpi_enddef(ncid)
!      call handle_error_pnetcdf('nf_enddef: output',status)
!
!
!! write data
!
!      start(1)=1
!      edge(1)=1
!      tarr(1) = time
!      status=nfmpi_put_vara_double_all(
!     $     ncid,time_varid,start,edge,tarr)
!      call handle_error_pnetcdf(
!     $     'nf_put_vara_real: time',status)
!
!      start(1)=i_global(1)
!      start(2)=j_global(1)
!      start(3)=1
!      edge(1)=im
!      edge(2)=jm
!      edge(3)=1
!      status=nfmpi_put_vara_double_all(
!     $     ncid,tair_varid,start,edge,ta)
!      call handle_error_pnetcdf(
!     $     'nf_put_vara_real: tair',status)
!      status=nfmpi_put_vara_double_all(
!     $     ncid,rhum_varid,start,edge,sh)
!      call handle_error_pnetcdf(
!     $     'nf_put_vara_real: rhum',status)
!      status=nfmpi_put_vara_double_all(
!     $     ncid,rain_varid,start,edge,ra)
!      call handle_error_pnetcdf(
!     $     'nf_put_vara_real: rain',status)
!      status=nfmpi_put_vara_double_all(
!     $     ncid,tskn_varid,start,edge,ts)
!      call handle_error_pnetcdf(
!     $     'nf_put_vara_real: tskn',status)
!      status=nfmpi_put_vara_double_all(
!     $     ncid,pres_varid,start,edge,pr)
!      call handle_error_pnetcdf(
!     $     'nf_put_vara_real: pres',status)
!      status=nfmpi_put_vara_double_all(
!     $     ncid,tcld_varid,start,edge,cl)
!      call handle_error_pnetcdf(
!     $     'nf_put_vara_real: tcld',status)
!!
!! close file
!      status=nfmpi_close(ncid)
!      call handle_error_pnetcdf('nf_close: output: '
!     $                           ,status)
!
!      return
!      end
!_______________________________________________________________________

      integer function get_var_real_1d(ncid, varid, start, edge, var)

        implicit none
# include "mpif.h"
# include "pnetcdf.inc"
        include 'realkind'


        integer                              ,intent(in)  :: ncid, varid
        integer(mpi_offset_kind),dimension(4),intent(in)  :: start, edge
        real(kind=rk),     dimension(edge(1)),intent(out) :: var

        integer       status, vtype
        integer(mpi_offset_kind) k,l,m,n
        real(kind=4)  f4, o4
        real(kind=8)  f8, o8
        real(kind=rk) fk, ok

        integer(kind=1), dimension(:), allocatable :: i1  ! byte
        integer(kind=2), dimension(:), allocatable :: i2  ! short
        integer(kind=4), dimension(:), allocatable :: i4  ! int
!        integer(kind=8), dimension(:), allocatable :: i8  ! long (not supported by netcdf?)
        real(kind=4)   , dimension(:), allocatable :: r4  ! float
        real(kind=8)   , dimension(:), allocatable :: r8  ! double
!        real(kind=16)  , dimension(:), allocatable :: rF  ! quadruple (not supported by netcdf?)

        fk = 1. ! scale factor
        ok = 0. ! add offset
        k = edge(1)
        l = edge(2)
        m = edge(3)
        n = edge(4)

        status = nfmpi_inq_vartype(ncid, varid, vtype)

! BYTE to REAL
        if (vtype==NF_BYTE) then
          allocate(i1(k))
          status=nfmpi_get_vara_int1_all(ncid,varid,start,edge,i1)
          var = 1.
          where (i1==0) var = 0.
          deallocate(i1)
! SHORT to REAL
        else if (vtype==NF_SHORT) then
! unpack
          status = nfmpi_inq_atttype(ncid, varid, "add_offset", vtype)
          if (vtype==NF_FLOAT) then
            status = nfmpi_get_att_real(ncid,varid,"add_offset",o4)
            ok = real(o4)
          else if (vtype==NF_DOUBLE) then
            status = nfmpi_get_att_double(ncid,varid,"add_offset",o8)
            ok = real(o8)
          else
            call handle_error_pnetcdf( 'Failed reading `add_offset`'
     $                               , status )
          end if
          status = nfmpi_inq_atttype(ncid, varid, "scale_factor", vtype)
          if (vtype==NF_FLOAT) then
            status = nfmpi_get_att_real(ncid,varid,"scale_factor",f4)
            fk = real(f4)
          else if (vtype==NF_DOUBLE) then
            status = nfmpi_get_att_double(ncid,varid,"scale_factor",f8)
            fk = real(f8)
          else
            call handle_error_pnetcdf( 'Failed reading `scale_factor`'
     $                               , status )
          end if
          allocate(i2(k))
          status=nfmpi_get_vara_int2_all(ncid,varid,start,edge,i2)
          var = real(i2)*fk+ok
          deallocate(i2)
        else if (vtype==NF_INT) then
          allocate(i4(k))
          status=nfmpi_get_vara_int_all(ncid,varid,start,edge,i4)
          var = real(i4)
          deallocate(i4)
        else if (vtype==NF_FLOAT) then
          allocate(r4(k))
          status=nfmpi_get_vara_real_all(ncid,varid,start,edge,r4)
          var = real(r4)
          deallocate(r4)
        else if (vtype==NF_DOUBLE) then
          allocate(r8(k))
          status=nfmpi_get_vara_double_all(ncid,varid,start,edge,r8)
          var = real(r8)
          deallocate(r8)
        end if
!      call handle_error_pnetcdf('nfmpi_get_vara_real_all',
!     $                          status,nf_noerr)
        get_var_real_1d = status

      end function



      integer function get_var_real_2d(ncid, varid, start, edge, var)

        implicit none
# include "mpif.h"
# include "pnetcdf.inc"
        include 'realkind'

        integer                              ,intent(in)  :: ncid, varid
        integer(mpi_offset_kind),dimension(4),intent(in)  :: start, edge
        real(kind=rk),
     &             dimension(edge(1),edge(2)),intent(out) :: var

        integer       status, vtype
        integer(mpi_offset_kind) k,l,m,n
        real(kind=4)  f4, o4, m4
        real(kind=8)  f8, o8, m8
        real(kind=rk) fk, ok, mv

        integer(kind=1), dimension(:,:), allocatable :: i1  ! byte
        integer(kind=2), dimension(:,:), allocatable :: i2  ! short
        integer(kind=4), dimension(:,:), allocatable :: i4  ! int
!        integer(kind=8), dimension(:,:), allocatable :: i8  ! long (not supported by netcdf?)
        real(kind=4)   , dimension(:,:), allocatable :: r4  ! float
        real(kind=8)   , dimension(:,:), allocatable :: r8  ! double
!        real(kind=16)  , dimension(:,:), allocatable :: rF  ! quadruple (not supported by netcdf?)

        fk = 1. ! scale factor
        ok = 0. ! add offset
        mv = 9.96920996838687e+36 ! missing value
        k = edge(1)
        l = edge(2)
        m = edge(3)
        n = edge(4)

        status = nfmpi_inq_atttype(ncid, varid, "add_offset", vtype)
        if (status==NF_NOERR) then
          if (vtype==NF_FLOAT) then
            status = nfmpi_get_att_real(ncid,varid,"add_offset",o4)
            ok = real(o4)
          else if (vtype==NF_DOUBLE) then
            status = nfmpi_get_att_double(ncid,varid,"add_offset",o8)
            ok = real(o8)
          end if
        end if

        status = nfmpi_inq_atttype(ncid, varid, "scale_factor", vtype)
        if (status==NF_NOERR) then
          if (vtype==NF_FLOAT) then
            status = nfmpi_get_att_real(ncid,varid,"scale_factor",f4)
            fk = real(f4)
          else if (vtype==NF_DOUBLE) then
            status = nfmpi_get_att_double(ncid,varid,"scale_factor",f8)
            fk = real(f8)
          end if
        end if

        status = nfmpi_inq_vartype(ncid, varid, vtype)

! BYTE to REAL
        if (vtype==NF_BYTE) then
          allocate(i1(k,l))
          status=nfmpi_get_vara_int1_all(ncid,varid,start,edge,i1)
          var = 1.
          where (i1==0) var = 0.
          deallocate(i1)
! SHORT to REAL
        else if (vtype==NF_SHORT) then
          allocate(i2(k,l))
          status=nfmpi_get_vara_int2_all(ncid,varid,start,edge,i2)
          var = real(i2)
          deallocate(i2)
        else if (vtype==NF_INT) then
          allocate(i4(k,l))
          status=nfmpi_get_vara_int_all(ncid,varid,start,edge,i4)
          var = real(i4)
          deallocate(i4)
        else if (vtype==NF_FLOAT) then
          allocate(r4(k,l))
          status=nfmpi_get_vara_real_all(ncid,varid,start,edge,r4)
          var = real(r4)
          deallocate(r4)
        else if (vtype==NF_DOUBLE) then
          allocate(r8(k,l))
          status=nfmpi_get_vara_double_all(ncid,varid,start,edge,r8)
          var = real(r8)
          deallocate(r8)
        end if
        get_var_real_2d = status

        status = nfmpi_inq_atttype(ncid, varid, "_FillValue", vtype)
        if (status==NF_NOERR) then
          if (vtype==NF_FLOAT) then
            status = nfmpi_get_att_real(ncid,varid,"_FillValue",m4)
            mv = real(m4)
          else if (vtype==NF_DOUBLE) then
            status = nfmpi_get_att_double(ncid,varid,"_FillValue",m8)
            mv = real(m8)
          else
            call handle_error_pnetcdf( 'Failed reading `_FillValue`'
     &                               , status )
          end if

          where (var==mv) var = 0.
        end if
        status = nf_noerr

        var = var*fk+ok

      end function



      integer function get_var_real_3d(ncid, varid, start, edge, var)

        implicit none
# include "mpif.h"
# include "pnetcdf.inc"
        include 'realkind'

        integer                              ,intent(in)  :: ncid, varid
        integer(mpi_offset_kind),dimension(4),intent(in)  :: start, edge
        real(kind=rk),
     &     dimension(edge(1),edge(2),edge(3)),intent(out) :: var

        integer       status, vtype
        integer(mpi_offset_kind) k,l,m,n
        real(kind=4)  f4, o4, m4
        real(kind=8)  f8, o8, m8
        real(kind=rk) fk, ok, mv

        integer(kind=1), dimension(:,:,:), allocatable :: i1  ! byte
        integer(kind=2), dimension(:,:,:), allocatable :: i2  ! short
        integer(kind=4), dimension(:,:,:), allocatable :: i4  ! int
!        integer(kind=8), dimension(:,:,:), allocatable :: i8  ! long (not supported by netcdf?)
        real(kind=4)   , dimension(:,:,:), allocatable :: r4  ! float
        real(kind=8)   , dimension(:,:,:), allocatable :: r8  ! double
!        real(kind=16)  , dimension(:,:,:), allocatable :: rF  ! quadruple (not supported by netcdf?)

        fk = 1. ! scale factor
        ok = 0. ! add offset
        mv = 9.96920996838687e+36 ! missing value
        k = edge(1)
        l = edge(2)
        m = edge(3)
        n = edge(4)

        status = nfmpi_inq_atttype(ncid, varid, "add_offset", vtype)
        if (status==NF_NOERR) then
          if (vtype==NF_FLOAT) then
            status = nfmpi_get_att_real(ncid,varid,"add_offset",o4)
            ok = real(o4)
          else if (vtype==NF_DOUBLE) then
            status = nfmpi_get_att_double(ncid,varid,"add_offset",o8)
            ok = real(o8)
          end if
        end if

        status = nfmpi_inq_atttype(ncid, varid, "scale_factor", vtype)
        if (status==NF_NOERR) then
          if (vtype==NF_FLOAT) then
            status = nfmpi_get_att_real(ncid,varid,"scale_factor",f4)
            fk = real(f4)
          else if (vtype==NF_DOUBLE) then
            status = nfmpi_get_att_double(ncid,varid,"scale_factor",f8)
            fk = real(f8)
          end if
        end if

        status = nfmpi_inq_vartype(ncid, varid, vtype)
! BYTE to REAL
        if (vtype==NF_BYTE) then
          allocate(i1(k,l,m))
          status=nfmpi_get_vara_int1_all(ncid,varid,start,edge,i1)
          var = 1.
          where (i1==0) var = 0.
          deallocate(i1)
! SHORT to REAL
        else if (vtype==NF_SHORT) then
          allocate(i2(k,l,m))
          status=nfmpi_get_vara_int2_all(ncid,varid,start,edge,i2)
          var = real(i2)
          deallocate(i2)
        else if (vtype==NF_INT) then
          allocate(i4(k,l,m))
          status=nfmpi_get_vara_int_all(ncid,varid,start,edge,i4)
          var = real(i4)
          deallocate(i4)
        else if (vtype==NF_FLOAT) then
          allocate(r4(k,l,m))
          status=nfmpi_get_vara_real_all(ncid,varid,start,edge,r4)
          var = real(r4)
          deallocate(r4)
        else if (vtype==NF_DOUBLE) then
          if (rk/=8) then
            allocate(r8(k,l,m))
            status=nfmpi_get_vara_double_all(ncid,varid,start,edge,r8)
            var = real(r8)
            deallocate(r8)
          else
            status=nfmpi_get_vara_double_all(ncid,varid,start,edge,var)
          end if
        end if

        get_var_real_3d = status

        ! the below code is optional
        status = nfmpi_inq_atttype(ncid, varid, "_FillValue", vtype)
        if (status==NF_NOERR) then
          if (vtype==NF_FLOAT) then
            status = nfmpi_get_att_real(ncid,varid,"_FillValue",m4)
            mv = real(m4)
          else if (vtype==NF_DOUBLE) then
            status = nfmpi_get_att_double(ncid,varid,"_FillValue",m8)
            mv = real(m8)
          else
            call handle_error_pnetcdf( 'Failed reading `_FillValue`'
     &                               , status )
          end if

          where (var==mv) var = 0.
        end if

        if (fk/=1. .and. ok/=0.) then
          var = var*fk+ok
        end if

      end function
