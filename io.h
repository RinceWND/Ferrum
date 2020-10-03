      character(len=256)
     &  grid_path, bc_path, clim_path, mean_path
     &, bulk_path, flux_path

      common/io/
     &  grid_path, bc_path, clim_path, mean_path
     &, bulk_path, flux_path

! Input variables' names
      character(len=24)
     &  dlwr_name,  ! Downward longwave radiation
     &  lht_name,   ! Latent heat flux
     &  lwr_name,   ! Longwave net radiation
     &  prat_name,  ! Precipitation rate
     &  pres_name,  ! Atm. pressure
     &  rhum_name,  ! Relative humidity
     &  sht_name,   ! Sensible heat flux
     &  sst_name,   ! Sea surface temperature
     &  swr_name,   ! Shortwave net radiation
     &  tair_name,  ! Air temperature
     &  tcld_name,  ! Total cloud cover
     &  umf_name,   ! U-component of momentum flux
     &  uwnd_name,  ! U-component of wind
     &  vmf_name,   ! V-component of momentum flux
     &  vwnd_name,  ! V-component of wind
     &  wgst_name   ! Wind gustiness [NOT IMPLEMENTED]

      common/input_vars/
     &  dlwr_name, lht_name, lwr_name, prat_name, pres_name
     &, rhum_name, sht_name, sst_name, swr_name, tair_name
     &, tcld_name, umf_name, uwnd_name, vmf_name, vwnd_name
     &, wgst_name