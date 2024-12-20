;-----------------------------------------------------------------
;; this script is just used to filter values by removing seasonal
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
  fu  = addfile(ddir+"ERA5_"+"u"+"level_daily_mean_1979-2022_r2.5.nc","r")
  fv  = addfile(ddir+"ERA5_"+"v"+"level_daily_mean_1979-2022_r2.5.nc","r")
  fw  = addfile(ddir+"ERA5_"+"w"+"level_daily_mean_1979-2022_r2.5.nc","r")
  ft  = addfile(ddir+"ERA5_"+"t"+"level_daily_mean_1979-2022_r2.5.nc","r")
; statistic_stability & term2
  t = ft->t(:,{1000:100},:,:)
  w = fw->w(:,{1000:100},:,:)
  t1 = dim_avg_n(t, 0)
  t1 := dim_avg_n(t1, (/1,2/))
  p = t&level*100
  p@units = "Pa"
  ss = static_stability(p, t1, 0, 0)

  s = conform_dims(dimsizes(w), ss, 1)
  term2 = w*s
  delete([/w,ss,s,t1/])
  print("term2 done.")

;T tendency K/s
  term1 = center_finite_diff_n(t,24*60*60.,False,0,0)
  print("term1 done.")

; advection
  u = fu->u(:,{1000:100},:,:)
  v = fv->v(:,{1000:100},:,:)

  term3 = advect_variable_cfd(u, v, t, t&lat, t&lon, False, "VgradT", "K*m/s", 0)
  delete([/u,v/])
  print("term3 done.")

; sum.
  cp = 1004.64    ;J/(K-kg)
  q = cp*(term1 - term2 - term3) 

  q@units = "J/(kg-s)"
  q@long_name  = "precipitation_heating"
  q@short_name = "diabatic_heating"

  copy_VarCoords(t, q)
  delete(t)

  xa = q
  do i = 0, dimsizes(q&level)-1
    xa(:,i,:,:) = iso_filter(q(:,i,:,:),201,nwa,nwb,1)
  end do
  print("calc done. "+ systemfunc("date"))
  copy_VarMeta(q,xa)

  create_nc_file(q,q&time,q&level,q&lat,q&lon,ddir+\
                "ERA5_"+"h"+"level_daily_mean_1979-2022_r2.5.nc","h")

  create_nc_file(xa,q&time,q&level,q&lat,q&lon,ddir+\
                "ERA5_"+"h"+"level_daily_anom_1979-2022_r2.5.nc","h")
  
  exit
end