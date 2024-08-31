;--------------------------------------------------------
; draw uwind anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour. 
; waf vector
; This is for test sig. 95% by monte carlo method.
; time: 2023-12-18
; 2.5x2.5.
;--------------------------------------------------------

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  udir   = ddir+"/"+"ERA5_daily_u250_1979-2022_r2.5.nc"
  ucdir  = ddir+"/"+"ERA5_daily_u250_1979-2022_r2.5.nc"

  lats = 10
  latn = 60
  lonw = 0
  lone = 210
  ;lonw = 200
  ;lone = 350
;--------------------------------------------------------
; u250
;--------------------------------------------------------
  f  = addfile(ucdir,"r")

  time = f->time
  TIME = cd_calendar(time, 0)
  mon  = toint(TIME(:,1))
  nt   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  uc   = dim_avg_n_Wrap(f->u(nt,0,:,:),0)
;contour plot
  resc = True
  resc@gsnFrame     = False
  resc@gsnDraw      = False
  resc@gsnMaximize  = False
  resc@gsnAddCyclic = True
  resc@gsnContourNegLineDashPattern = 1
  resc@gsnLeftString      = ""
  resc@gsnRightString     = ""
  resc@cnMonoLineColor    = True
  resc@cnLinesOn          = True
  resc@cnLineLabelsOn     = False
  resc@cnLineColor        = "magenta";magenta
  resc@cnLineThicknessF   = 2.5
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels      = (/30.,45.,60.,75.,90./);(/30,35,40,45/)
  delete([/time,f,TIME,mon,nt/])
  print("draw clim.u complete. "+ systemfunc("date"))
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  fu = addfile(udir,"r")

  u  = fu->u(:,0,:,:)

  time = u&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  u := runave_n_Wrap(u, 5, 0, 0)
  print("complete data selling "+ systemfunc("date"))
;-------------------------------------------------------
;; selling MJO.
  delim  = ","
  rmmf   = asciiread(rmmdir,-1,"string")
  iyear  = str_get_field(rmmf,1,delim)
  imon   = str_get_field(rmmf,2,delim)
  iday   = str_get_field(rmmf,3,delim)
  phaset = str_get_field(rmmf,6,delim)
  ampt   = str_get_field(rmmf,7,delim)

  iyear  := toint(iyear(1::))
  imon   := toint(imon(1::))
  iday   := toint(iday(1::))
  phaset := toint(phaset(1::))
  ampt   := tofloat(ampt(1::))

  amp   = new((/dimsizes(time)/),float)
  phase = new((/dimsizes(time)/),integer)
  amp   = 0.
  phase = 0

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i).and.iday.eq.day(i))
    amp(i)   = ampt(nt)
    phase(i) = phaset(nt)
    delete(nt)
  end do

  nn = ispan(0, dimsizes(time)-1, 1)
  amp!0   = "timeind"
  phase!0 = "timeind"
  amp&timeind   = nn
  phase&timeind = nn
  delete([/nn,ampt,phaset,iyear,imon,iday/])
  print("complete RMM selling "+ systemfunc("date"))
;----------------------------------------------------------------------
;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
  printVarSummary(u)
  llu = new((/dimsizes(u&lat),dimsizes(u&lon)/),float)
  copy_VarCoords(u(0,:,:), llu)
  llu = 0.

  n1 = ind(year.eq.1993.and.mon.eq.1.and.day.eq.6)
  n2 = ind(year.eq.1993.and.mon.eq.1.and.day.eq.21) 

  llu = dim_avg_n(u(n1:n2,:,:),0) 

  ;do ntt = 1,4
  ;  nit = nind + ntt
  ;  nindd := array_append_record(nindd, nit, 0)
  ;  delete(nit)
  ;end do
  ;llu(0,:,:) = dim_avg_n(u(nindd,:,:), 0)
  ;delete([/nt,nind,ntt,nindd/])
  print("complete calc leadlag data. "+ systemfunc("date"))
;----------------------------------------------------------------------
;; resources
  cmap = read_colormap_file("MPL_RdBu")  ;MPL_RdBu
;shading u plot
  resu = True
  resu@gsnFrame           = False
  resu@gsnDraw            = False
  resu@gsnMaximize        = False
  resu@gsnAddCyclic       = False
  resu@gsnMajorLatSpacing = 15
  resu@gsnMajorLonSpacing = 30
  resu@gsnRightStringFontHeightF = 0.015 
  resu@gsnLeftStringFontHeightF  = 0.015

  resu@gsnRightStringOrthogonalPosF = 0.01
  resu@gsnLeftStringOrthogonalPosF  = 0.01

  resu@tmXTOn = False
  resu@tmYROn = False
  resu@tmXBLabelFont = 21
  resu@tmYLLabelFont = 21
  resu@tmXBLabelFontHeightF = 0.01
  resu@tmYLLabelFontHeightF = 0.01

  resu@mpMinLatF = lats
  resu@mpMaxLatF = latn
  resu@mpMinLonF = lonw
  resu@mpMaxLonF = lone
  resu@mpFillOn        = False
  resu@mpLandFillColor = "white"
  resu@mpGridAndLimbDrawOrder = "PostDraw"
  resu@mpGridLineColor        = "Gray60"

  resu@cnFillOn               = True
  resu@cnFillPalette          = "WhiteBlueGreenYellowRed"
  resu@cnLinesOn              = False
  resu@cnLineLabelsOn         = False
  resu@cnLevelSelectionMode   = "Manuallevels"
  resu@cnMinLevelValF         = 30;15;15
  resu@cnMaxLevelValF         = 90;45;90.
  resu@cnLevelSpacingF        = 15.;5;15.
  resu@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  resu@lbLabelBarOn           = True

  respoly = True
  respoly@gsLineColor      = "purple"
  respoly@gsLineThicknessF = 4.5

  print("complete res. "+ systemfunc("date"))
;----------------------------------------------------------------------
;; plot.
  wks = gsn_open_wks("x11", pdir+"/"+"Fig3_U__MEEAJS")

  plot = new(1,graphic)
  dum1 = plot
  Label  = (/"(a) 250hPa U@P2","(b) 250hPa U@P6"/)

  n = 0
  x  = llu
  resu@gsnLeftString  = "~F21~"+Label(n)
  plot = gsn_csm_contour_map(wks, x, resu)
  plotut  = gsn_csm_contour(wks, uc, resc)

  lonmin =  10
  lonmax =  70
  latmin = -100
  latmax =  100
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  respoly@gsLineDashPattern = 0
  dum2 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)

  lonmin =  110
  lonmax =  180
  latmin = -100
  latmax =  100
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  respoly@gsLineDashPattern = 1
  dum3 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)

  overlay(plot,plotut)
  overlay(plot,dum2)
  overlay(plot,dum3)
  draw(plot)

  i = 0
  do while(i.eq.0)
    i= 0
  end do