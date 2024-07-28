;-----------------------------------------------------------------
;; this script is just used to filter values by removing seasonal
;; cycle. Lanczos filter.
;; time: 2023-11-22
;-----------------------------------------------------------------
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"

;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
begin
  var  = "olr"
  nwa  = 100.
  nwb  = 20.

  ddir = "/home/sunming/data5/cuixy/global_jets/data"

  fin  = addfile(ddir+"/"+var+"_daily_mean_1979-2022.nc","r")

  xv   = fin->$var$(time|:,lat|:,lon|:)

  xa   = iso_filter(xv,201,nwa,nwb,0)

  copy_VarCoords(xv,xa)

  fout = addfile(ddir+"/"+var+"_daily_anom_1979-2022.nc","c")

  fout->$var$ = xa

end

