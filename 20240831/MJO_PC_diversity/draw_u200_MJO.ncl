
;---------------------------------------------------------
;function: 返回位相3对应的日期 格式：YYYYMMDD
;---------------------------------------------------------
function return_days(timt:integer,xdir)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,pt,at,ip
begin
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nf = asciiread(xdir,-1,"string")
  is = str_get_field(nf,1,delim)
  ie = str_get_field(nf,2,delim)
  is := toint(is(1::))
  ie := toint(ie(1::))
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ;event all days.
  nda := ndtooned(nda)
  nda := nda(ind(.not.ismissing(nda)))  ;index for z.
  return nda

end
;----------------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  udir = ddir+"/"+"ERA5_u200_daily_mean_1979-2022_r2.5.nc"
;;jet events
  ndir = ddir+"/"+"IOnP_date_omi.csv"
  fdir = ddir+"/"+"IOfP_date_omi.csv"
  sdir = ddir+"/"+"IOsP_date_omi.csv"
  idir = ddir+"/"+"IOiP_date_omi.csv"

  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"

  lats = 0
  latn = 90
;----------------------------------------------------------
; main code.
;----------------------------------------------------------
  f = addfile(udir,"r")
  u = f->u(:,{200},{lats:latn},:)
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  timi = cd_calendar(time,-2)
  u := runave_n_Wrap(u,5,0,0)
;;OMI
  delim = ","
  rmmf  = asciiread(omidir,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
  phaset = str_get_field(rmmf,4,delim)
  ampt = str_get_field(rmmf,7,delim)

  iyear := toint(iyear(1::))
  imon  := toint(imon(1::))
  iday  := toint(iday(1::))
  phaset := tointeger(phaset(1::))
  ampt := tofloat(ampt(1::))

  phase = new((/dimsizes(time)/),integer)
  amp   = new((/dimsizes(time)/),float)
  phase = 0
  amp   = 0.

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i)\
             .and.iday.eq.day(i))
    phase(i) = phaset(nt)
    amp(i) = ampt(nt)
    delete(nt)
  end do
  delete([/phaset,ampt,iyear,imon,iday/])
  print("complete OMI selling "+ systemfunc("date"))

  IOnP = return_days(timi,ndir) ;IOnP
  IOiP = return_days(timi,idir) ;IOiP
  IOfP = return_days(timi,fdir) ;IOfP
  IOsP = return_days(timi,sdir) ;IOsP

  un = dim_avg_n_Wrap(u(IOnP,:,:),0)
  ui = dim_avg_n_Wrap(u(IOiP,:,:),0)
  uf = dim_avg_n_Wrap(u(IOfP,:,:),0)
  us = dim_avg_n_Wrap(u(IOsP,:,:),0)

  uni = un - ui
  ufi = uf - ui
  usi = us - ui
  copy_VarCoords(ui, uni)
  copy_VarCoords(ui, ufi)
  copy_VarCoords(ui, usi)
; u250 resources.
;--plot
  cmap = read_colormap_file("MPL_RdBu") 

  resc = True
  resc@gsnFrame     = False
  resc@gsnDraw      = False
  resc@gsnMaximize  = False

  resc@gsnLeftString = ""
  resc@gsnRightString = ""

  resc@mpFillOn = True
  resc@mpLandFillColor = "white"
  resc@mpCenterLonF    = 180.
  resc@mpMinLatF       = 0
  resc@mpLandFillColor = "white"
  resc@mpGridAndLimbDrawOrder = "PostDraw"
  resc@mpGridLineColor        = "Gray60"

  resc@cnFillOn = False
  resc@cnLinesOn       = True
  resc@cnLineLabelsOn  = False
  resc@cnLineColor        = "magenta";magenta
  resc@cnLineThicknessF   = 2.5
  resc@cnInfoLabelOn   = False
  resc@cnInfoLabelPerimOn = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels  = (/20.,30.,40.,50.,60.,70./)

  res = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False
  res@gsnAddCyclic = True
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing  = 30

  res@gsnLeftStringOrthogonalPosF  = 0.012
  res@gsnRightStringOrthogonalPosF = 0.012
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnRightStringFontHeightF = 0.015
  res@tmXBLabelFontHeightF  = 0.012
  res@tmYLLabelFontHeightF  = 0.012

  res@mpFillOn = True
  res@mpLandFillColor = "white"
  res@mpCenterLonF    = 180.
  res@mpMinLatF       = 0

  res@cnFillOn       = True              
  res@cnFillPalette  = cmap(::-1,:)  
  res@cnLinesOn = False    
  res@cnLevelSelectionMode = "ManualLevels"     
  res@cnMinLevelValF       = -5         
  res@cnMaxLevelValF       =  5            
  res@cnLevelSpacingF      =  1
  res@cnLinesOn  = False
  res@cnInfoLabelPerimOn = False
  res@cnInfoLabelOns   = False

  res@lbLabelBarOn = False
  res@pmLabelBarOrthogonalPosF = 0.25
 
  res@tmXBLabelFont = 21
  res@tmYRLabelFont = 21

  res@vpWidthF = 0.5
  res@vpHeightF= 0.15

  wks = gsn_open_wks("x11","u200")

  label = (/"~F21~(a) IOiP","~F21~(b) IOnP","~F21~(c) IOfP","~F21~(d) IOsP"/)

  res@vpXF     = 0.2
  res@vpYF     = 0.7
  res@gsnLeftString  = label(1)
  plot2 = gsn_csm_contour_map(wks, uni, res)
  plot2c = gsn_csm_contour(wks, un, resc)
  overlay(plot2, plot2c)

  res@vpXF     = 0.2
  res@vpYF     = 0.5
  res@gsnLeftString  = label(2)
  plot3 = gsn_csm_contour_map(wks, ufi, res)
  plot3c = gsn_csm_contour(wks, uf, resc)
  overlay(plot3, plot3c)

  res@lbLabelBarOn = True
  res@vpXF     = 0.2
  res@vpYF     = 0.3
  res@gsnLeftString  = label(3)
  plot4 = gsn_csm_contour_map(wks, usi, res)
  plot4c = gsn_csm_contour(wks, us, resc)
  overlay(plot4, plot4c)

  resc@vpWidthF = 0.5
  resc@vpHeightF= 0.15
  resc@vpXF     = 0.2
  resc@vpYF     = 0.95
  resc@gsnLeftString  = label(0)
  plot1 = gsn_csm_contour_map(wks, ui, resc)

  draw(plot1)
  draw(plot2)
  draw(plot3)
  draw(plot4)

  i = 0
  do while (i.eq.0)
    i = 0
  end do