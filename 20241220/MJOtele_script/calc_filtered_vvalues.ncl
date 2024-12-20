;-----------------------------------------------------------------
;; this script is just used to filter values by removing seasonal
;; cycle. Lanczos filter.
;; time: 2024-05-13
;-----------------------------------------------------------------
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_4D.ncl"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
begin
  
  var  = var

  nwa  = 100.
  nwb  = 20.

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/vs/"
  fin  = addfile(ddir+"ERA5_"+var+"level_daily_mean_1979-2022_r2.5.nc","r")

  x  = fin->$var$(:,{1000:100},:,:)
  xa = x
  do i = 0, dimsizes(x&level)-1
    xa(:,i,:,:) = iso_filter(x(:,i,:,:),201,nwa,nwb,1)
  end do
  print("calc done. "+ systemfunc("date"))
  copy_VarMeta(x,xa)

  create_nc_file(xa,x&time,x&level,x&lat,x&lon,ddir+\
                "ERA5_"+var+"level_daily_anom_1979-2022_r2.5.nc",var)
  
end