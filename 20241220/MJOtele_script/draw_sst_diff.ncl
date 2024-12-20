;---------------------------------------------
;统计每个事件经过的月份，平均
;---------------------------------------------

;---------------------------------------------------------
;function: 返回第一天和最后一天的日期 格式化为：YYYYMM01
;---------------------------------------------------------
function return_days(timt:integer,xdir)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,pt,at,ip
begin
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nf = asciiread(xdir,-1,"string")
  is = str_get_field(nf,1,delim)
  ie = str_get_field(nf,2,delim)
  is := is(1::)
  ie := ie(1::)
  dtime_map = (/4,2,2/)
  is_dtime = str_split_by_length(is, dtime_map)
  ie_dtime = str_split_by_length(ie, dtime_map)
  is_dtime(:,2) = "01"
  ie_dtime(:,2) = "01"
  nda = new((/dimsizes(is),2/),integer)
  do i = 0,dimsizes(is)-1
    its = toint(str_concat(is_dtime(i,:)))
    ite = toint(str_concat(ie_dtime(i,:)))
    nda(i,0) = its
    nda(i,1) = ite
  end do

  return nda

end

  lats = -30
  latn =  30
  lonw =  40
  lone =  300

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  sstdir = ddir+"/"+"ERSSTv5.nc"
  udir = ddir+"/"+"ERA5_u200_daily_mean_1979-2022_r2.5.nc"
;;jet events
  fdir = ddir+"/"+"IOfP_date_omi.csv"
  sdir = ddir+"/"+"IOsP_date_omi.csv"

  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"
;---------------------------------------------
;main code.
;---------------------------------------------
; calc monthly anomaly.
  f = addfile(sstdir,"r")
  sst = f->sst(:,{lats:latn},{lonw:lone})
  clm = clmMonTLL(sst)
  sa  = calcMonAnomTLL(sst, clm)
;omi index
  fu = addfile(udir,"r")
  time = fu->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  timi = cd_calendar(time,-2)
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
  delete([/time,TIME,year,mon,day,phaset,ampt,iyear,imon,iday/])
  print("complete OMI selling "+ systemfunc("date"))

; select month (start,end) (events,2)
  IOfP = return_days(timi,fdir) ;IOfP
  IOsP = return_days(timi,sdir) ;IOsP
  delete(timi)

  time = sst&time
  timi = cd_calendar(time,-2)

  print("select date done. "+ systemfunc("date"))

;select month to composite.
  mP = ind(timi.ge.IOfP(0,0).and.timi.le.IOfP(0,1))
  do i = 1, dimsizes(IOfP(:,0))-1
    mt  = ind(timi.ge.IOfP(i,0).and.timi.le.IOfP(i,1))
    mP := array_append_record(mP, mt, 0)
    delete(mt)
  end do
  safP = dim_avg_n_Wrap(sa(mP,:,:), 0)
  delete(mP)

  mP = ind(timi.ge.IOsP(0,0).and.timi.le.IOsP(0,1))
  do i = 1, dimsizes(IOsP(:,0))-1
    mt  = ind(timi.ge.IOsP(i,0).and.timi.le.IOsP(i,1))
    mP := array_append_record(mP, mt, 0)
    delete(mt)
  end do
  sasP = dim_avg_n_Wrap(sa(mP,:,:), 0)
  delete(mP)

;;for statas.
  saf = new((/dimsizes(IOfP(:,0)),dimsizes(sa&lat),dimsizes(sa&lon)/),float)
  saf = 0.
  sas = saf
  do i = 0, dimsizes(IOsP(:,0))-1
    mt  = ind(timi.ge.IOsP(i,0).and.timi.le.IOsP(i,1))
    sas(i,:,:) = dim_avg_n_Wrap(sa(mt,:,:), 0)
    delete(mt)
  end do

  do i = 0, dimsizes(IOfP(:,0))-1
    mt  = ind(timi.ge.IOfP(i,0).and.timi.le.IOfP(i,1))
    saf(i,:,:) = dim_avg_n_Wrap(sa(mt,:,:), 0)
    delete(mt)
  end do

  copy_VarCoords(sa(0,:,:), saf(0,:,:))
  copy_VarCoords(sa(0,:,:), sas(0,:,:))

  fP_std=dim_variance_n_Wrap(saf, 0)
  fP_std=sqrt(fP_std/(dimsizes(IOfP(:,0)));;标准差
  fP_std=where(fP_std.eq.0,fP_std@_FillValue,fP_std)

  sP_std=dim_variance_n_Wrap(sas, 0)
  sP_std=sqrt(sP_std/dimsizes(IOsP(:,0)));;标准差
  sP_std=where(sP_std.eq.0,sP_std@_FillValue,sP_std)

  t_fP=safP/fP_std
  confi_fP=safP
  confi_fP=student_t(t_fP, dimsizes(IOfP(:,0))-1)

  t_sP=sasP/sP_std
  confi_sP=sasP
  confi_sP=student_t(t_sP, dimsizes(IOsP(:,0))-1)
;--------------------------------------------------------------------
  wks  = gsn_open_wks ("x11", "hov" )  
        
  res = True
  res@gsnFrame           = False
  res@gsnDraw            = False
  res@gsnMaximize        = False
  res@gsnAddCyclic       = False
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnLeftStringOrthogonalPosF  = 0.01

  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBLabelFont = 21
  res@tmYLLabelFont = 21
  res@tmXBLabelFontHeightF = 0.01
  res@tmYLLabelFontHeightF = 0.01

  res@mpMinLatF = lats
  res@mpMaxLatF = latn
  res@mpMinLonF = lonw
  res@mpMaxLonF = lone
  res@mpCenterLonF = 120
  res@mpFillOn        = False
  res@mpLandFillColor = "white"
  res@mpGridAndLimbDrawOrder = "PostDraw"
  res@mpGridLineColor        = "Gray60"

  res@cnFillOn               = True
  res@cnFillPalette          = "temp_diff_18lev"
  res@cnLinesOn              = False
  res@cnLineLabelsOn         = False
  res@cnLevelSelectionMode   = "Manuallevels"
  res@cnMinLevelValF         = -0.5
  res@cnMaxLevelValF         =  0.5
  res@cnLevelSpacingF        =  0.05
  res@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  res@lbLabelBarOn           = False

  opt                      =True
  opt@gsnShadeFillType     ="pattern"
  opt@gsnShadeLow          =17;打点
  opt@gsnAddCyclic         =True
  opt@cnFillDotSizeF       =0.1;改变点大小的默认值

  sres                     =True
  sres@gsnDraw             =False
  sres@gsnFrame            =False
  sres@cnLinesOn           =False
  sres@gsnLeftString       =""
  sres@gsnRightString      =""
  sres@cnLevelSelectionMode="ExplicitLevels"
  sres@cnLevels            =(/0.05,0.01/)
  sres@cnFillPalette       ="GMT_gray"
  sres@cnFillColors        =(/5,7,-1/)
  sres@cnLineLabelsOn      =False
  sres@cnInfoLabelOn       =False
  sres@lbLabelBarOn        =False
  sres@cnFillScaleF        = 1
  sres@cnFillDotSizeF      = 0.005


  label = (/"IOfP","IOsP"/)
  plot = new(2,graphic)
  base = plot

  n = 0
  res@gsnLeftString = "~F21~"+label(n)
  plot(n) = gsn_csm_contour_map(wks, safP, res) 
  n = 1
  res@gsnLeftString = "~F21~"+label(n)
  plot(n) = gsn_csm_contour_map(wks, sasP, res) 

  base(0)=gsn_csm_contour(wks,confi_fP, sres)
  base(0)=gsn_contour_shade(base(0),0.1,1,opt)
  overlay(plot(0),base(0))
  base(1)=gsn_csm_contour(wks,confi_sP, sres)
  base(1)=gsn_contour_shade(base(1),0.1,1,opt)
  overlay(plot(1),base(1))


  resP = True
  resP@gsnMaximize  = True
  resP@gsnFrame     = False
  resP@gsnDraw      = True

  resP@pmLabelBarOrthogonalPosF = 0.0
  resP@pmLabelBarParallelPosF   = 0.0
  resP@pmLabelBarWidthF         = 0.6

  resP@gsnPanelLabelBar = True
  resP@lbTitleOn = True
  resP@lbTitleString = ""
  resP@lbTitlePosition  = "Right"
  resP@lbTitleDirection = "Across"
  resP@lbTitleFontHeightF = 0.015
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.013
  ;resP@lbTitleStringParallelPosF = -0.1

  gsn_panel(wks,plot,(/1,1/),resP)

  i = 0
  do while(i.eq.0)
    i = 0
  end do