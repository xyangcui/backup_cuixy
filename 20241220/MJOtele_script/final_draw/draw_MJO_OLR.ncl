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
  pct1 = new((/dimsizes(is),dim2,dim3/),float) 
  pct1 = 0.
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    pt = pht(its:ite)
    at = ampt(its:ite)
    nda = ind_time(its:ite)
    ;pt  := pt(ind(at.gt.1.))
    ;nda := nda(ind(at.gt.1.))
    ;ip = ind(pt.eq.3.or.pt.eq.2)
    ip = ind(pt.eq.5)
    ip@_FillValue = -999
    if (any(ip.eq.ip@_FillValue).eq."Missing")
      delete([/ip,its,ite,at,pt,nda/])
      continue
    else
      if (dimsizes(ip).eq.1)
        pct1(i,:,:) = zt(nda(ip)+5,:,:)
      else
        pct1(i,:,:) = dim_avg_n(zt(nda(ip)+5,:,:),0)
      end if
    end if
    delete([/ip,its,ite,at,pt,nda/])
  end do
  return pct1

end
;--------------------------------------------------------------------------
;function: 计算每种事件并返回矩阵的程序。
;--------------------------------------------------------------------------
function calc_pro7(timt:integer,xdir,zt:float,pht:integer,ampt:float)
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
  pct1 = new((/dimsizes(is),dim2,dim3/),float) 
  pct1 = 0.
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    pt = pht(its:ite)
    at = ampt(its:ite)
    nda = ind_time(its:ite)
    ;pt  := pt(ind(at.gt.1.))
    ;nda := nda(ind(at.gt.1.))
    ip = ind(pt.eq.6)
    ip@_FillValue = -999
    if (any(ip.eq.ip@_FillValue).eq."Missing")
      delete([/ip,its,ite,at,pt,nda/])
      continue
    else
      if (dimsizes(ip).eq.1)
        pct1(i,:,:) = zt(nda(ip)+10,:,:)
      else
        pct1(i,:,:) = dim_avg_n(zt(nda(ip)+0,:,:),0)
      end if
    end if
    delete([/ip,its,ite,at,pt,nda/])
  end do
  return pct1

end

;--------------------------------------------------------------------------
;function: 计算每种事件并返回矩阵的程序。
;--------------------------------------------------------------------------
function calc_olr(timt:integer,xdir,ztt:float,pht:integer,ampt:float)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,igood,pt,at,ip,out,pct1
begin
  lats = -15
  latn =  15
  lonw =  135
  lone =  165
  zt = wgt_areaave_Wrap(ztt(:,{lats:latn},{lonw:lone}), 1., 1., 0)
  ;zt = dim_standardize_n(ztd, 0, 0)
  dims = dimsizes(zt)

  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nf = asciiread(xdir,-1,"string")
  is = str_get_field(nf,1,delim)
  ie = str_get_field(nf,2,delim)
  is := toint(is(1::))
  ie := toint(ie(1::))
  ind_time = ispan(0,dimsizes(timt)-1,1)
  ;event by event.
  pct1 = new((/dimsizes(is),26/),float) 
  pct1 = 0.
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    pt  = pht(its:ite)
    nda = ind_time(its:ite)
    ip  = ind(pt.eq.5)
    ip@_FillValue = -999
    if (any(ip.eq.ip@_FillValue).eq."Missing")
      delete([/ip,its,ite,pt,nda/])
      continue
    else
      do j =0,25
        if (dimsizes(ip).eq.1)
          pct1(i,j) = zt(nda(ip)+j)
        else
          pct1(i,j) = dim_avg_n(zt(nda(ip)+j),0)
        end if
      end do

    end if
    delete([/ip,its,ite,pt,nda/])
  end do

  return pct1

end

;------------------------------------------------------------------
  lats = -25
  latn =  30
  lonw =  50
  lone =  210

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
  undef("time")
  time = olr&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  undef("timi")
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

  undef("phase")
  undef("amp")
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
  olr_iosp = calc_olr(timi,sdir,olr,phase,amp)
  olr_iofp = calc_olr(timi,fdir,olr,phase,amp)

  olr_iosp := dim_avg_n(olr_iosp,0)
  olr_iofp := dim_avg_n(olr_iofp,0)

  print(olr_iosp)
  print(olr_iofp)

  IOsP = calc_pro(timi,sdir,olr,phase,amp) ;IOsP
  IOfP = calc_pro(timi,fdir,olr,phase,amp) ;IOfP

  IOsP!0 = "time"
  IOsP!1 = "lat"
  IOsP!2 = "lon"

  IOfP!0 = "time"
  IOfP!1 = "lat"
  IOfP!2 = "lon"
;;value for plot.
  If = dim_avg_n(IOfP,0)
  Is = dim_avg_n(IOsP,0)

  copy_VarCoords(olr(0,:,:),If)
  copy_VarCoords(olr(0,:,:),Is)
;;for ttest.
  undef("nf")
  undef("ns")
  nf = dimsizes(IOsP(:,0,0))
  ns = dimsizes(IOfP(:,0,0))
;(t-test)
  doff = equiv_sample_size(IOfP(lat|:,lon|:,time|:), 0.1, 0)
  dff  = toint(wgt_areaave(doff, 1., 1., 0))

  dofs = equiv_sample_size(IOsP(lat|:,lon|:,time|:), 0.1, 0)
  dfs  = toint(wgt_areaave(dofs, 1., 1., 0))
;;sst
  f_var=dim_variance_n(IOfP, 0)
  f_std=sqrt(f_var/nf)
  f_std=where(f_std.eq.0,f_std@_FillValue,f_std);

  t_f=If/f_std
  confi_f=If
  confi_f=student_t(t_f, nf-1)

  s_var=dim_variance_n(IOsP, 0)
  s_std=sqrt(s_var/ns)
  s_std=where(s_std.eq.0,s_std@_FillValue,s_std)

  t_s=Is/s_std
  confi_s=Is
  confi_s=student_t(t_s, ns-1)

  copy_VarCoords(If, confi_f)
  copy_VarCoords(If, confi_s)

  Id = Is - If
  copy_VarCoords(Is,Id) 
  confi_d = ttest(Is, s_var, dfs, If, f_var, dff, True, False)
  copy_VarCoords(If, confi_d) 

  Is = smth9_Wrap(Is, 0.5, 0.25, False)
  If = smth9_Wrap(If, 0.5, 0.25, False)
;-------------------------------------------------------------------- 
  wks = gsn_open_wks("x11", "bar")
  cmap = read_colormap_file("BlueWhiteOrangeRed")         
;shading u plot
  res = True
  res@gsnFrame           = False
  res@gsnDraw            = False
  res@gsnMaximize        = False
  res@gsnAddCyclic       = False
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30
  res@gsnLeftStringFontHeightF  = 0.02
  res@gsnLeftStringOrthogonalPosF  = 0.01;-0.15

  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBLabelFont = 21
  res@tmYLLabelFont = 21
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015

  res@mpMinLatF = lats
  res@mpMaxLatF = latn
  res@mpMinLonF = lonw
  res@mpMaxLonF = lone
  res@mpCenterLonF = 140
  res@mpFillOn        = False
  res@mpLandFillColor = "white"
  res@mpGridAndLimbDrawOrder = "PostDraw"
  res@mpGridLineColor        = "Gray60"

  res@cnFillOn               = True
  res@cnFillPalette          = cmap
  res@cnLinesOn              = False
  res@cnLineLabelsOn         = False
  res@cnLevelSelectionMode   = "Manuallevels"
  res@cnMinLevelValF         = -20
  res@cnMaxLevelValF         =  20
  res@cnLevelSpacingF        =  4.
  res@cnLabelBarEndStyle     = "ExcludeOuterBoxes"

  res@lbTitleOn = True
  res@lbTitleString = "[W m~S~-2~N~]"
  res@lbTitlePosition  = "Right"
  res@lbTitleDirection = "Across"
  res@lbTitleFontHeightF = 0.015
  res@lbBoxLinesOn       = True
  res@lbLabelFontHeightF = 0.012
  res@pmLabelBarOrthogonalPosF = 0.15

  opt                      =True
  opt@gsnShadeFillType     ="pattern"
  opt@gsnShadeLow          =17;打点
  opt@gsnAddCyclic         =True
  opt@cnFillDotSizeF       =0.001
  opt@cnFillColor          ="gray64"

  sres                     =True
  sres@gsnDraw             =False
  sres@gsnFrame            =False
  sres@cnLinesOn           =False
  sres@gsnLeftString       =""
  sres@gsnRightString      =""
  sres@cnLevelSelectionMode="ExplicitLevels"
  sres@cnLevels            =(/0.15,0.1,0.05,0.01/)
  sres@cnFillPalette       ="GMT_gray"
  sres@cnFillColor         ="gray64"
  sres@cnLineLabelsOn      =False
  sres@cnInfoLabelOn       =False
  sres@lbLabelBarOn        =False
  sres@cnFillScaleF        = 1.5
  sres@cnFillDotSizeF      = 0.003

  n = 0
  res@vpXF = 0.15
  res@vpYF = 0.59
  res@lbLabelBarOn  = False
  res@gsnLeftString = "~F21~ (a) OLR (IOsP P6)"
  plots = gsn_csm_contour_map(wks, Is, res) 
  bases=gsn_csm_contour(wks,confi_s, sres)
  bases=gsn_contour_shade(bases,0.1,1,opt)
  overlay(plots,bases)

  n = 1
  res@vpXF = 0.15
  res@vpYF = 0.325
  res@lbLabelBarOn  = False
  res@gsnLeftString = "~F21~ (b) OLR (IOfP P6)"
  plotf = gsn_csm_contour_map(wks, If, res) 
  basef=gsn_csm_contour(wks,confi_f, sres)
  basef=gsn_contour_shade(basef,0.1,1,opt)
  overlay(plotf,basef)

  ;draw(plots)
  ;draw(plotf)

;P67.
  undef("IOsP")
  undef("IOfP")
  IOsP = calc_pro7(timi,sdir,olr,phase,amp) ;IOsP
  IOfP = calc_pro7(timi,fdir,olr,phase,amp) ;IOfP

  IOsP!0 = "time"
  IOsP!1 = "lat"
  IOsP!2 = "lon"

  IOfP!0 = "time"
  IOfP!1 = "lat"
  IOfP!2 = "lon"
;;value for plot.
  If = dim_avg_n(IOfP,0)
  Is = dim_avg_n(IOsP,0)

  copy_VarCoords(olr(0,:,:),If)
  copy_VarCoords(olr(0,:,:),Is)
;;for ttest.
  undef("nf")
  undef("ns")
  nf = dimsizes(IOsP(:,0,0))
  ns = dimsizes(IOfP(:,0,0))

;(t-test)
  doff = equiv_sample_size(IOfP(lat|:,lon|:,time|:), 0.1, 0)
  dff  =  toint(wgt_areaave(doff, 1., 1., 0))

  dofs = equiv_sample_size(IOsP(lat|:,lon|:,time|:), 0.1, 0)
  dfs  =  toint(wgt_areaave(dofs, 1., 1., 0))
;;sst
  f_var=dim_variance_n(IOfP, 0)
  f_std=sqrt(f_var/nf)
  f_std=where(f_std.eq.0,f_std@_FillValue,f_std);

  t_f=If/f_std
  confi_f=If
  confi_f=student_t(t_f, nf-1)

  s_var=dim_variance_n(IOsP, 0)
  s_std=sqrt(s_var/ns)
  s_std=where(s_std.eq.0,s_std@_FillValue,s_std)

  t_s=Is/s_std
  confi_s=Is
  confi_s=student_t(t_s, ns-1)

  copy_VarCoords(If, confi_f)
  copy_VarCoords(If, confi_s)

  Id = Is - If
  copy_VarCoords(Is,Id)
  confi_d = ttest(Is, s_var, dfs, If, f_var, dff, True, False)
  copy_VarCoords(If, confi_d)

  Is = smth9_Wrap(Is, 0.5, 0.25, False)
  If = smth9_Wrap(If, 0.5, 0.25, False)
;-------------------------------------------------------------------- 
  undef("cmap")
  cmap = read_colormap_file("BlueWhiteOrangeRed")         

  res@vpXF = 0.15+0.35
  res@vpYF = 0.59
  res@lbLabelBarOn  = False
  res@gsnLeftString = "~F21~ (c) OLR (IOsP P7)"
  plotst = gsn_csm_contour_map(wks, Is, res) 
  basest=gsn_csm_contour(wks,confi_s, sres)
  basest=gsn_contour_shade(basest,0.1,1,opt)
  overlay(plotst,basest)

  res@vpXF = 0.15+0.35
  res@vpYF = 0.325
  res@lbLabelBarOn  = False
  res@gsnLeftString = "~F21~ (d) OLR (IOfP P7)"
  plotft = gsn_csm_contour_map(wks, If, res) 
  baseft=gsn_csm_contour(wks,confi_f, sres)
  baseft=gsn_contour_shade(baseft,0.1,1,opt)
  overlay(plotft,baseft)

  resp = True
  resp@gsLineColor = "black"
  resp@gsLineThicknessF = 3.
  resp@gsLineDashPattern = 0

  polyx = (/130,170,170,130,130/)
  polyy = (/-15,-15,15,15,-15/)

  dum1 = gsn_add_polyline(wks, plots, polyx, polyy, resp)
  dum2 = gsn_add_polyline(wks, plotf, polyx, polyy, resp)
  dum3 = gsn_add_polyline(wks, plotst, polyx, polyy, resp)
  dum4 = gsn_add_polyline(wks, plotft, polyx, polyy, resp) 
  ;draw(plotst)
  ;draw(plotft)

  resP = True
  resP@gsnMaximize      = True
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelYWhiteSpacePercent = 2

  resP@lbLabelFont       = 21
  resP@lbBoxLinesOn      = True
  resP@gsnPanelRowSpec   = True

  resP@lbTitleOn = True
  resP@lbTitleString = "[W m~S~-2~N~]"
  resP@lbTitlePosition  = "Right"
  resP@lbTitleDirection = "Across"
  resP@lbTitleFontHeightF = 0.015
  resP@lbBoxLinesOn       = True
  resP@lbLabelFontHeightF = 0.012

  resP@pmLabelBarWidthF         =  0.6
  resP@pmLabelBarParallelPosF   =  0.005
  resP@pmLabelBarOrthogonalPosF =  0.

  gsn_panel(wks,(/plots,plotst,plotf,plotft/),(/2,2/),resP)

  i = 0 
  do while(i.eq.0)
    i = 0
  end do