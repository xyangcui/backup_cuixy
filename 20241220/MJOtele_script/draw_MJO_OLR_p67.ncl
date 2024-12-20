
;--------------------------------------------------------------------------
;function: 计算每种事件并返回矩阵的程序。
;--------------------------------------------------------------------------
function calc_pro(timt:integer,xdir,zt:float,pht:integer,ampt:float)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,igood,pt,at,ip,out,pct1
begin
  dims = dimsizes(zt)
  dim2 = dims(1)
  dim3 = dims(2)
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nf = asciiread(xdir,-1,"string")
  is = str_get_field(nf,1,delim)
  ie = str_get_field(nf,2,delim)
  is := toint(is(1::))
  ie := toint(ie(1::))
  ind_time = ispan(0,dimsizes(timt)-1,1)
  ;event by event.
  ;pct1 = new((/dimsizes(is),dim2,dim3/),float) 
  ;pct1 = 0.
  ;do i = 0,dimsizes(is)-1
  ;  its = ind(timt.eq.is(i))
  ;  ite = ind(timt.eq.ie(i))
  ;  pt = pht(its:ite)
  ;  at = ampt(its:ite)
  ;  nda = ind_time(its:ite)
  ;  pt  := pt(ind(at.gt.1.))
  ;  nda := nda(ind(at.gt.1.))
  ;  ip = ind(pt.eq.6.or.pt.eq.7)
  ;  ip@_FillValue = -999
  ;  if (any(ip.eq.ip@_FillValue).eq."Missing")
  ;    delete([/ip,its,ite,at,pt,nda/])
  ;    continue
  ;  else
  ;    if (dimsizes(ip).eq.1)
  ;      pct1(i,:,:) = zt(nda(ip),:,:)
  ;    else
  ;      pct1(i,:,:) = dim_avg_n(zt(nda(ip),:,:),0)
  ;    end if
  ;  end if
  ;  delete([/ip,its,ite,at,pt,nda/])
  ;end do

  ;day by day.
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  nda := ndtooned(nda)
  igood = ind(.not.ismissing(nda))
  nda := nda(igood)  ;index for z.
  pt = pht(nda)
  at = ampt(nda)
  ;pt  := pt(ind(at.gt.1.))
  ;nda := nda(ind(at.gt.1.))
  delete(igood)
  ph = ind(pt.eq.6.or.pt.eq.7)
  nda := nda(ph)
  pct1 = zt(nda,:,:)

  return pct1

end
;------------------------------------------------------------------

  lats = -15
  latn =  20
  lonw =  60
  lone =  220

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  odir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
;;jet events
  fdir = ddir+"/"+"IOfP_date_omi.csv"
  sdir = ddir+"/"+"IOsP_date_omi.csv"

  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"
;-----------------------------------------------------------------
; main code.
;-----------------------------------------------------------------
  f = addfile(odir,"r")
  olr = f->olr(:,{lats:latn},{lonw:lone})
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
  IOsP = calc_pro(timi,sdir,olr,phase,amp) ;IOsP
  IOfP = calc_pro(timi,fdir,olr,phase,amp) ;IOfP
  print("jet events done. "+ systemfunc("date"))

  IOsP!0 = "time"
  IOsP!1 = "lat"
  IOsP!2 = "lon"

  IOfP!0 = "time"
  IOfP!1 = "lat"
  IOfP!2 = "lon"
;;value for plot.
  If = dim_avg_n(IOfP,0)
  Is = dim_avg_n(IOsP,0)
  Id = Is - If

  copy_VarCoords(olr(0,:,:),If)
  copy_VarCoords(olr(0,:,:),Is)
  copy_VarCoords(olr(0,:,:),Id)
;;for ttest.
  nf = dimsizes(IOsP(:,0,0))
  ns = dimsizes(IOfP(:,0,0))

;(t-test)
  doff = equiv_sample_size(IOfP(lat|:,lon|:,time|:), 0.1, 0)
  dff  =  toint(wgt_areaave(doff, 1., 1., 0))

  dofs = equiv_sample_size(IOsP(lat|:,lon|:,time|:), 0.1, 0)
  dfs  =  toint(wgt_areaave(dofs, 1., 1., 0))

  print(dff)
  print(dfs)
;;sst
  f_var=dim_variance_n(IOfP, 0)
  f_std=sqrt(f_var/nf)
  f_std=where(f_std.eq.0,f_std@_FillValue,f_std);

  t_f=If/f_std
  confi_f=If
  confi_f=student_t(t_f, dff)

  s_var=dim_variance_n(IOsP, 0)
  s_std=sqrt(s_var/ns)
  s_std=where(s_std.eq.0,s_std@_FillValue,s_std)

  t_s=Is/s_std
  confi_s=Is
  confi_s=student_t(t_s, dfs)

  copy_VarCoords(If, confi_f)
  copy_VarCoords(If, confi_s)
;;diff
  confi_d = ttest(If,f_var,dff, Is,s_var,dfs, True, False)
  alpha = 100.*(1. - confi_d) 

  Is = smth9_Wrap(Is, 0.5, 0.25, False)
  If = smth9_Wrap(If, 0.5, 0.25, False)
  Id = smth9_Wrap(Id, 0.5, 0.25, False)
;--------------------------------------------------------------------
  wks  = gsn_open_wks ("x11", "hov" )  

  cmap = read_colormap_file("BlueWhiteOrangeRed")         
;shading u plot
  res = True
  res@gsnFrame           = False
  res@gsnDraw            = False
  res@gsnMaximize        = False
  res@gsnAddCyclic       = False
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
  res@mpMinLonF = lonw
  res@mpMaxLonF = lone
  ;res@mpCenterLonF = 120
  res@mpFillOn        = False
  res@mpLandFillColor = "white"
  res@mpGridAndLimbDrawOrder = "PostDraw"
  res@mpGridLineColor        = "Gray60"

  res@cnFillOn               = True
  res@cnFillPalette          = cmap
  res@cnLinesOn              = False
  res@cnLineLabelsOn         = False
  res@cnLevelSelectionMode   = "Manuallevels"
  ;res@cnMinLevelValF         = -14.
  ;res@cnMaxLevelValF         =  14.
  ;res@cnLevelSpacingF        =  2.
  res@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  res@lbLabelBarOn           = False

  opt                      =True
  opt@gsnShadeFillType     ="pattern"
  opt@gsnShadeLow          =17;打点
  opt@gsnAddCyclic         =True
  opt@cnFillDotSizeF       =0.001

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
  sres@cnFillScaleF        = 1.3
  sres@cnFillDotSizeF      = 0.005

  label = (/"IOsP","IOfP","IOsP-IOfP"/)
  plot = new(3,graphic)
  base = plot

  n = 0
  res@gsnLeftString = "~F21~"+label(n)
  plot(n) = gsn_csm_contour_map(wks, Is, res) 
  base(n)=gsn_csm_contour(wks,confi_s, sres)
  base(n)=gsn_contour_shade(base(n),0.05,1,opt)
  overlay(plot(n),base(n))

  n = 1
  res@gsnLeftString = "~F21~"+label(n)
  plot(n) = gsn_csm_contour_map(wks, If, res) 
  base(n)=gsn_csm_contour(wks,confi_f, sres)
  base(n)=gsn_contour_shade(base(n),0.05,1,opt)
  overlay(plot(n),base(n))

  n = 2
  res@gsnLeftString = "~F21~"+label(n)
  plot(n) = gsn_csm_contour_map(wks, Id, res) 
  base(n)=gsn_csm_contour(wks,confi_d, sres)
  base(n)=gsn_contour_shade(base(n),0.05,1,opt)
  overlay(plot(n),base(n))

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

  gsn_panel(wks,(/plot(2),plot(1),plot(0)/),(/1,1,1/),resP)

  i = 0
  do while(i.eq.0)
    i = 0
  end do