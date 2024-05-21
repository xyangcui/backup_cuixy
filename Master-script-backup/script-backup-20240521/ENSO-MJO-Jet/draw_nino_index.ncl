;-----------------------------------------------------------------
;; calculate iso of jet indices (10-90 days).
;; calculate std. of jet indices.
;; draw std. per year with nino index.
;; time: 2024-03-27 author: cuixy.
;-----------------------------------------------------------------

nwa  = 90.
nwb  = 10.
ymd_srt = 19790701
ymd_end = 20220630

;;---------------------------------------------
;; nino index.
;;---------------------------------------------
  ninodir = "/home/sunming/data5/cuixy/data/clim_indices/Nino_3.4_index.csv"
  pdir = "/home/sunming/data5/cuixy/ENSO-MJO-Jet/figure/bar_nino_index"

  delim  = ","
  ninof  = asciiread(ninodir,-1,"string")
  iyear  = str_get_field(ninof,1,delim)
  iJan   = str_get_field(ninof,2,delim)
  iFeb   = str_get_field(ninof,3,delim)
  iDec   = str_get_field(ninof,4,delim)

  iyear:= toint(iyear(1::))
  iJan := tofloat(iJan(1::))
  iFeb := tofloat(iFeb(1::))
  iDec := tofloat(iDec(1::))

  nino = new((/43/),float)
  nino = 0.
  do i = 0,42
    nino(i) = (iJan(i)+iFeb(i)+iDec(i))/3.
  end do

  std = dim_stddev(nino)

  res   = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False

  ;res@trYMinF   = -3.0   
  ;res@trYMaxF   =  3.0   

  res@tmXBMode    = "Explicit"
  ;res@tmXBValues  = levels
  res@tmXBLabels  = (/"79/80","89/90","99/00","09/10","19/20"/)

  res@tiYAxisString = ""     
  res@tiXAxisString = ""
  res@tiMainString  = "" 

  res@gsnYRefLine           = 0.                
  res@gsnXYBarChart         = True            
  res@gsnAboveYRefLineColor = "red"           
  res@gsnBelowYRefLineColor = "blue"          


  wks   = gsn_open_wks ("eps", pdir )  

  plot = new(1,graphic)

  plot = gsn_csm_xy(wks, ispan(1980, 2022, 1), nino, res)

  resp = True
  resp@gsLineColor  = "black"
  resp@gsLineThicknessF = 3.5

  dum1 = gsn_add_polyline(wks, plot, (/1970,2025/), (/0.7,0.7/), resp)
  dum2 = gsn_add_polyline(wks, plot, (/1970,2025/), (/-0.7,-0.7/), resp)

  draw(plot)

  ;i = 0
  ;do while (i.eq.0)
  ;  i = 0
  ;end do
