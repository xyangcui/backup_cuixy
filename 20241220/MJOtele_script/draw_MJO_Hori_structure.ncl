
;--------------------------------------------------------------------------
;function: 计算每种事件并返回矩阵的程序。
;--------------------------------------------------------------------------
function calc_pro(timt:integer,xdir,zt:float,pht:integer,ampt:float)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,igood,pt,at,ip,out,pct1
begin
  dims = dimsizes(zt)
  dim2 = dims(1)
  dim3 = dims(2)
  pct1 = new((/8,dim2,dim3/),float) 
  pct1 = 0.
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
  igood = ind(.not.ismissing(nda))
  nda := nda(igood)  ;index for z.
  pt = pht(nda)
  at = ampt(nda)
 ; pt  := pt(ind(at.gt.1.))
  ;nda := nda(ind(at.gt.1.))
  delete(igood)
  ;select amp >1. and phase.
  do m = 0, 7
    np = m+1
    ip = ind(pt.eq.m+1)
    ip@_FillValue = -999
    if (any(ip.eq.ip@_FillValue).eq."Missing")
      continue
    else
      pct1(m,:,:) = dim_avg_n(zt(nda(ip)+10,:,:),0)
      delete(ip)
    end if
  end do

  return pct1

end
;------------------------------------------------------------------

  lats = -15
  latn =  20

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  odir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
;;jet events
  ndir = ddir+"/"+"IOnP_date_omi.csv"
  fdir = ddir+"/"+"IOfP_date_omi.csv"
  sdir = ddir+"/"+"IOsP_date_omi.csv"
  idir = ddir+"/"+"IOiP_date_omi.csv"

  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"
;-----------------------------------------------------------------
; main code.
;-----------------------------------------------------------------
  f = addfile(odir,"r")
  olr = f->olr(:,{lats:latn},:)
  time = olr&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  timi = cd_calendar(time, -2) 
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

;;MJO events.
;;jet events (year,mon,day) (events,days,lat,lon)
  print("jet events begin. "+ systemfunc("date"))
  IOnP = calc_pro(timi,ndir,olr,phase,amp) ;IOnP
  IOiP = calc_pro(timi,idir,olr,phase,amp) ;IOiP
  IOfP = calc_pro(timi,fdir,olr,phase,amp) ;IOfP
  IOsP = calc_pro(timi,sdir,olr,phase,amp) ;IOsP
  print("jet events done. "+ systemfunc("date"))

  IOP = (IOnP+IOiP+IOfP+IOsP)/4.
  copy_VarCoords(olr(0,:,:),IOP(0,:,:))
  copy_VarCoords(olr(0,:,:),IOfP(0,:,:))
  copy_VarCoords(olr(0,:,:),IOsP(0,:,:))
  copy_VarCoords(olr(0,:,:),IOnP(0,:,:))
  copy_VarCoords(olr(0,:,:),IOiP(0,:,:))
;--------------------------------------------------------------------
  wks  = gsn_open_wks ("x11", "hov" )  

  cmap = read_colormap_file("BlueWhiteOrangeRed")         
;shading u plot
  res = True
  res@gsnFrame           = False
  res@gsnDraw            = False
  res@gsnMaximize        = False
  res@gsnAddCyclic       = True
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30
  ;res@gsnRightStringFontHeightF = 0.015 
  res@gsnLeftStringFontHeightF  = 0.015

  ;res@gsnRightStringOrthogonalPosF = 0.01
  res@gsnLeftStringOrthogonalPosF  = 0.01

  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBLabelFont = 21
  res@tmYLLabelFont = 21
  res@tmXBLabelFontHeightF = 0.01
  res@tmYLLabelFontHeightF = 0.01

  res@mpMinLatF = lats
  res@mpMaxLatF = latn
  ;res@mpMinLonF = lonw
  ;res@mpMaxLonF = lone
  res@mpCenterLonF = 120
  res@mpFillOn        = False
  res@mpLandFillColor = "white"
  res@mpGridAndLimbDrawOrder = "PostDraw"
  res@mpGridLineColor        = "Gray60"

  res@cnFillOn               = True
  res@cnFillPalette          = cmap
  res@cnLinesOn              = False
  res@cnLineLabelsOn         = False
  res@cnLevelSelectionMode   = "Manuallevels"
  res@cnMinLevelValF         = -20.
  res@cnMaxLevelValF         =  20.
  res@cnLevelSpacingF        =  4.
  res@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  res@lbLabelBarOn           = False

  label = (/"p1","p2","p3","p4","p5","p6","p7","p8"/)
  plot = new(8,graphic)

  pv = IOfP
  copy_VarCoords(IOfP,pv)
  do n = 0, 8-1
    res@gsnLeftString = "~F21~"+label(n)
    plot(n) = gsn_csm_contour_map(wks, pv(n,:,:), res) 
  end do

  resP = True
  resP@gsnMaximize  = True
  resP@gsnFrame     = False
  resP@gsnDraw      = True

  resP@pmLabelBarOrthogonalPosF = 0.0
  resP@pmLabelBarParallelPosF   = 0.0
  resP@pmLabelBarWidthF         = 0.6

  resP@gsnPanelLabelBar = True
  resP@lbTitleOn = True
  resP@lbTitleString = "[W m~S~-2~N~]"
  resP@lbTitlePosition  = "Right"
  resP@lbTitleDirection = "Across"
  resP@lbTitleFontHeightF = 0.015
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.013
  ;resP@lbTitleStringParallelPosF = -0.1

  gsn_panel(wks,plot,(/1,1,1,1,1,1,1,1/),resP)

  i = 0
  do while(i.eq.0)
    i = 0
  end do