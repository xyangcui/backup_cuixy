;-----------------------------------------------------------------
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_4D.ncl"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
begin
  
  nwa  = 100.
  nwb  = 20.

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/vs/"
  fq  = addfile(ddir+"ERA5_"+"q"+"level_daily_mean_1979-2022_r2.5.nc","r")
  ft  = addfile(ddir+"ERA5_"+"t"+"level_daily_mean_1979-2022_r2.5.nc","r")

  q = fq->q
  t = ft->t

  p = q&level*100

  theta_e = pot_temp_equiv(tofloat(p), t, q, 1, "q") 
  copy_VarCoords(q, theta_e)
  print("calc EPT done.")
  
  xa = theta_e

  do i = 0, dimsizes(q&level)-1
    xa(:,i,:,:) = iso_filter(theta_e(:,i,:,:),201,nwa,nwb,1)
  end do

  print("calc done. "+ systemfunc("date"))
  copy_VarMeta(theta_e,xa)

  create_nc_file(theta_e,q&time,q&level,q&lat,q&lon,ddir+\
                "ERA5_"+"ept"+"level_daily_mean_1979-2022_r2.5.nc","ept")

  create_nc_file(xa,q&time,q&level,q&lat,q&lon,ddir+\
                "ERA5_"+"ept"+"level_daily_anom_1979-2022_r2.5.nc","ept")

end