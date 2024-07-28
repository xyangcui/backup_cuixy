load "/home/sunming/data5/cuixy/script/ncl_func/cal_WFsf.ncl"
;;----------------------------------------------------------
function siglevel(ovalue,mvalue)
begin

  rv = ovalue
  lati = rv&lat
  loni = rv&lon

  do i = 0, dimsizes(lati)-1
    do j = 0, dimsizes(loni)-1
      if(ovalue(i,j).lt.mvalue(0,i,j).or.ovalue(i,j).gt.mvalue(1,i,j)) then
        rv(i,j) = 2.
      else
        rv(i,j) = -2.
      end if
    end do
  end do

  return rv

end

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  phidir = ddir+"/"+"ERA5_daily_anom_phi250_1979-2022_2.5x2.5.nc"
  udir   = ddir+"/"+"ERA5_daily_anom_u250_1979-2022_2.5x2.5.nc"
  vdir   = ddir+"/"+"ERA5_daily_anom_v250_1979-2022_2.5x2.5.nc"
  tdir   = ddir+"/"+"ERA5_daily_anom_tvert_1979-2022_2.5x2.5.nc"

  ucdir  = ddir+"/"+"ERA5_daily_u250_1979-2022_r2.5.nc"

  sdir = "/home/sunming/data5/cuixy/global_jets/data/monte_carlo/Fig8_usft_sig.nc"

  lats = 0
  latn = 60
  lonw = 200
  lone = 350
;--------------------------------------------------------
; u250
;--------------------------------------------------------
  f    = addfile(ucdir,"r")

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
  resc@cnLevels      = (/30.,40.,46.,60.,70./)

  delete([/time,f,TIME,mon,nt/])

  print("draw clim.u complete. "+ systemfunc("date"))
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  fphi = addfile(phidir,"r")
  fu   = addfile(udir,"r")
  fv   = addfile(vdir,"r")
  ft = addfile(tdir,"r")

  t    = ft->t(:,{lats:latn},{lonw:lone})
  phi  = fphi->phi(:,{lats:latn},{lonw:lone})
  u    = fu->u(:,{lats:latn},{lonw:lone})
  v    = fv->v(:,{lats:latn},{lonw:lone})

  lat = phi&lat
  lon = phi&lon

  time = phi&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  print("complete data selling "+ systemfunc("date"))

  fsig = addfile(sdir,"r")
  mvalue = fsig->mvalue
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
;-------------------------------------------------------------

  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.2.and.day.eq.28)

  amp   := amp(n1:n2)
  phase := phase(n1:n2)
  mon   := mon(n1:n2)
  time  := time(n1:n2)

  nt     = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  amp   := amp(nt)
  phase := phase(nt)

  delete(nt)
  ind_DJF = amp&timeind   ; indice of DJF between 1979-2022 DJF.

;----------------------------------------------------------------------
;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
  llphi = new((/dimsizes(phi&lat),dimsizes(phi&lon)/),float)
  copy_VarCoords(phi(0,:,:), llphi(:,:))
  llphi = 0.

  llu = new((/dimsizes(u&lat),dimsizes(u&lon)/),float)
  copy_VarCoords(u(0,:,:), llu(:,:))
  llu = 0.

  llv = new((/dimsizes(v&lat),dimsizes(v&lon)/),float)
  copy_VarCoords(v(0,:,:), llv(:,:))
  llv = 0.

  llt = new((/dimsizes(t&lat),dimsizes(t&lon)/),float)
  copy_VarCoords(t(0,:,:), llt(:,:))
  llt = 0.

  nt = ind(phase.eq.8.and.amp.gt.1.)
  nind = ind_DJF(nt)
  ntt  = 3
  nindd = nind + ntt
  do ntt = 4,10
    nit = nind + ntt
    nindd := array_append_record(nindd, nit, 0)
    delete(nit)
  end do

  llphi = dim_avg_n(phi(nindd,:,:), 0)
  llu   = dim_avg_n(u(nindd,:,:), 0)
  llv   = dim_avg_n(v(nindd,:,:), 0)
  llt   = dim_avg_n(t(nindd,:,:), 0)
  delete([/nt,nind,ntt,nindd/])

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
  resu@mpFillOn        = False
  resu@mpLandFillColor = "white"
  resu@mpGridAndLimbDrawOrder = "PostDraw"
  resu@mpGridLineColor        = "Gray60"

  resu@cnFillOn               = True
  resu@cnFillPalette          = cmap(::-1,:);(66:146,:)
  resu@cnLinesOn              = False
  resu@cnLineLabelsOn         = False
  resu@cnLevelSelectionMode   = "Manuallevels"
  resu@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  resu@lbLabelBarOn           = False
;; vector res.
  resw = True
  resw@gsnFrame     = False
  resw@gsnDraw      = False
  resw@gsnMaximize  = False
  resw@gsnAddCyclic = False
  resw@vcGlyphStyle            = "LineArrow"
  resw@vcLineArrowThicknessF   = 2.5
  resw@vcRefMagnitudeF         = 5.
  resw@vcRefLengthF            = 0.08
  resw@vcLineArrowColor        = "mediumblue"
  resw@vcRefAnnoString2On      = False
  resw@vcRefAnnoString1        = "$VMG$" +" m s~S~-1~N~"
  resw@vcRefAnnoPerimOn        = False
  resw@vcRefAnnoFontHeightF    = 0.02
  resw@vcRefAnnoBackgroundColor = -1

  sres   = True 
  sres@cnLineDrawOrder    = "PostDraw"                    
  sres@gsnDraw            = False                          
  sres@gsnFrame           = False                        
  sres@cnLineLabelsOn     = False                              
  sres@cnLinesOn          = False                           
  sres@cnInfoLabelOn      = False
  sres@cnFillOn           = False                       
  sres@cnInfoLabelOn      = False
  sres@cnFillScaleF       = 0.5
  sres@cnFillDotSizeF     = 0.003
  sres@gsnLeftString      =" "
  sres@gsnRightString     =" "
  sres@cnMonoFillColor    = True
  sres@cnFillColor        = "gray"

  opt = True
  opt@gsnShadeFillType   = "pattern"
  opt@gsnShadeHigh       = 17
  opt@gsnShadeFillScaleF = 1.5

  respoly = True
  respoly@gsLineColor       = "purple"
  respoly@gsLineThicknessF  = 3.5
  respoly@gsLineDashPattern = 2.

  print("complete res. "+ systemfunc("date"))

  lonmin = -110
  lonmax = -40
  latmin = -100
  latmax =  100
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
;----------------------------------------------------------------------
  wks = gsn_open_wks("eps", pdir+"/"+"Fig7_U_SF_T_NAJS")

  plott   = new(1,graphic)
  plotphi = plott

  plotu   = plott
  plotu2  = plott
  plotu3  = plott

  plotv   = plott

  dum1    = plott
  dum2    = plott
  dum3    = plott

  Label  = (/"~F21~(a) 250hPa_U@P4","~F21~(a) 250hPa_U@P8"/)
  Label1 = (/"~F21~(a) 250hPa_SF_P4","~F21~(a) 250hPa_SF_UV_P8"/)
  Label2 = (/"~F21~(c) Temp_P4","~F21~(b) Temp_P8"/)

  xf = (/0.05,0.05/)
  n = 1
  x  = llphi
  y  = llu
  z  = llt
  l  = llv

  mvalueu = mvalue(0,1,:,:,:)
  mvaluev = mvalue(1,1,:,:,:)
  mvalues = mvalue(2,1,:,:,:)
  mvaluet = mvalue(3,1,:,:,:)

  xt = siglevel(x,mvalues)
  yt1 = siglevel(y,mvalueu)
  zt = siglevel(z,mvaluet)

  lt = mask(l,y.lt.mvalueu(0,:,:).or.y.gt.mvalueu(0,:,:).or.\
            l.lt.mvaluev(0,:,:).or.l.gt.mvaluev(1,:,:),True)
  yt = mask(y,y.lt.mvalueu(0,:,:).or.y.gt.mvalueu(0,:,:).or.\
            l.lt.mvaluev(0,:,:).or.l.gt.mvaluev(1,:,:),True)
      
  copy_VarCoords(y,lt)
  copy_VarCoords(y,yt)

  resu@lbLabelBarOn = True
  resu@pmLabelBarWidthF  = 0.4
  resu@pmLabelBarHeightF = 0.05
  resu@pmLabelBarParallelPosF   = 0.5
  resu@pmLabelBarOrthogonalPosF = 0.12
  resw@vcRefAnnoOn = True
  resw@vcRefAnnoOrthogonalPosF =  0.025
  resw@vcRefAnnoParallelPosF   =  0.95

  delete([/mvalueu,mvaluev,mvalues,mvaluet/])

  x = x*10^(-6)
  resu@vpWidthF  = 0.43
  resu@vpHeightF = 0.25
  resu@vpXF = xf(n)
  resu@vpYF = 0.98

  resu@mpMinLatF   = 20
  resu@mpMaxLatF   = 55
  resu@mpMinLonF   = -120
  resu@mpMaxLonF   = -30
  resu@gsnLeftString  = "~F21~"+Label(n)

  resu@cnMinLevelValF  = -4.
  resu@cnMaxLevelValF  =  6.
  resu@cnLevelSpacingF =  1.

  resu@lbTitleOn = True
  resu@lbTitleString = "[m s~S~-1~N~]"
  resu@lbTitlePosition  = "Right"
  resu@lbTitleDirection = "Across"
  resu@lbTitleFontHeightF = 0.015
  resu@lbBoxLinesOn       = True
  resu@gsnPanelRowSpec    = True
  resu@lbLabelFontHeightF = 0.013

  plotu = gsn_csm_contour_map(wks, y, resu)
  dum2  = gsn_csm_contour(wks,yt1,sres)
  dum2  = gsn_contour_shade(dum2,1.9,-999.,opt)
  overlay(plotu,dum2)
;; draw phi.
  resu@vpHeightF = 0.45
  resu@vpXF = xf(n)
  resu@vpYF = 0.9

  resu@mpMinLatF   =  0
  resu@mpMaxLatF   =  60
  resu@mpMinLonF   = -120
  resu@mpMaxLonF   = -30
  resu@gsnLeftString  = "~F21~"+Label1(n)

  resu@lbTitleString = "[10~S~6~N~ m~S~2~N~ s~S~-1~N~]"

  resu@cnMinLevelValF    = -4.
  resu@cnMaxLevelValF    =  6.
  resu@cnLevelSpacingF   =  1.

  plotu2 = gsn_csm_contour_map(wks, x, resu)
  plotv  = gsn_csm_vector(wks, yt(::2,::2), lt(::2,::2), resw)
  dum1   = gsn_csm_contour(wks,xt,sres)
  dum1  = gsn_contour_shade(dum1,1.9,-999.,opt)
  overlay(plotu2,dum1)
  overlay(plotu2,plotv)

; temp.
  resu@vpHeightF = 0.25
  resu@vpXF = xf(n)
  resu@vpYF = 0.4

  resu@mpMinLatF   =  10
  resu@mpMaxLatF   =  55
  resu@mpMinLonF   = -120
  resu@mpMaxLonF   = -30
  resu@gsnLeftString  = "~F21~"+Label2(n)

  resu@lbTitleString = "[K]"
  resu@lbLabelBarOn = True
  resu@cnMinLevelValF  = -0.6
  resu@cnMaxLevelValF  =  0.35
  resu@cnLevelSpacingF =  0.05

  resu@pmLabelBarWidthF  = 0.4
  resu@pmLabelBarHeightF = 0.05
  resu@pmLabelBarOrthogonalPosF = 0.12

  plotu3 = gsn_csm_contour_map(wks, z, resu)
  dum3  = gsn_csm_contour(wks,zt,sres)
  dum3   = gsn_contour_shade(dum3,1.9,-999.,opt)
  overlay(plotu3,dum3)

  plotut  = gsn_csm_contour(wks, uc, resc)
  plotut1 = gsn_csm_contour(wks, uc, resc)
  plotut2 = gsn_csm_contour(wks, uc, resc)
  overlay(plotu,plotut)
  overlay(plotu2,plotut1)
  overlay(plotu3,plotut2)

  dum4 = gsn_add_polyline(wks, plotu,  polyx, polyy, respoly)
  dum5 = gsn_add_polyline(wks, plotu2, polyx, polyy, respoly)
  dum6 = gsn_add_polyline(wks, plotu3, polyx, polyy, respoly)
  overlay(plotu ,dum4)
  overlay(plotu2,dum5)
  overlay(plotu3,dum6)

  draw(plotu2)
  draw(plotu3)

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

  phidir = ddir+"/"+"ERA5_daily_anom_phi250_1979-2022_2.5x2.5.nc"
  udir   = ddir+"/"+"ERA5_daily_anom_u250_1979-2022_2.5x2.5.nc"
  olrdir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
  ucdir  = ddir+"/"+"ERA5_daily_u250_1979-2022_r2.5.nc.nc"
  sigdir = ddir+"/"+"monte_carlo/Fig9_sf_sig.nc"

  lats = -20
  latn =  90
  lonw = 100
  lone = 360
;--------------------------------------------------------
;; for calc WAF.
  u1dir = ddir+"/"+"ERA5_daily_u250_1979-2022_r2.5.nc"
  undef("vdir")
  vdir = ddir+"/"+"ERA5_daily_v250_1979-2022_r2.5.nc"

  fu = addfile(u1dir,"r")
  fv = addfile(vdir,"r")

  u0 = fu->u(:,{250},{lats:latn},{lonw:lone})
  v0 = fv->v(:,{250},{lats:latn},{lonw:lone})

  undef("time")
  time = u0&time
  undef("lat")
  lat  = u0&lat
  undef("lon")
  lon  = u0&lon

  undef("TIME")
  TIME = cd_calendar(time,0)
  undef("year")
  year = toint(TIME(:,0))
  undef("mon")
  mon  = toint(TIME(:,1))
  undef("day")
  day  = toint(TIME(:,2))

  wgtl = filwgts_lanczos(201, 0, 1./100., -999., 1.)
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  fphi = addfile(phidir,"r")
  fu   = addfile(udir,"r")
  fv   = addfile(vdir,"r")
  folr = addfile(olrdir,"r")

  olr  = folr->olr(:,{lats:latn},{lonw:lone})
  undef("phi")
  phi  = fphi->phi(:,{lats:latn},{lonw:lone})
  undef("u")
  u    = fu->u(:,{lats:latn},{lonw:lone})
  print("complete data selling "+ systemfunc("date"))

;; sig
  fsig = addfile(sigdir,"r")
  undef("mvalue")
  mvalue = fsig->mvalue
;-------------------------------------------------------
;; selling MJO.

  delim = ","
  rmmf  = asciiread(rmmdir,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
  phaset = str_get_field(rmmf,6,delim)
  ampt   = str_get_field(rmmf,7,delim)

  iyear  := toint(iyear(1::))
  imon   := toint(imon(1::))
  iday   := toint(iday(1::))
  phaset := toint(phaset(1::))
  ampt   := tofloat(ampt(1::))

  undef("amp")
  amp   = new((/dimsizes(time)/),float)
  undef("phase")
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
;-------------------------------------------------------------

  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.2.and.day.eq.28)

  amp   := amp(n1:n2)
  phase := phase(n1:n2)
  mon   := mon(n1:n2)
  time  := time(n1:n2)
  u1 = u0(n1:n2,:,:)
  v1 = v0(n1:n2,:,:)

  undef("nt")
  nt     = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  um     = dim_avg_n_Wrap(u1(nt,:,:),0)
  vm     = dim_avg_n_Wrap(v1(nt,:,:),0)
  amp   := amp(nt)
  phase := phase(nt)

  delete(nt)
  ind_DJF = amp&timeind   ; indice of DJF between 1979-2022 DJF.
;----------------------------------------------------------------------
;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
  undef("llphi")
  llphi = new((/4,dimsizes(phi&lat),dimsizes(phi&lon)/),float)
  copy_VarCoords(phi(0,:,:), llphi(0,:,:))
  llphi = 0.

  undef("llu")
  llu = new((/4,dimsizes(u&lat),dimsizes(u&lon)/),float)
  copy_VarCoords(u(0,:,:), llu(0,:,:))
  llu = 0.

  llolr = new((/4,dimsizes(olr&lat),dimsizes(olr&lon)/),float)
  copy_VarCoords(olr(0,:,:), llolr(0,:,:))
  llolr = 0.

  lluw = new((/4,dimsizes(u&lat),dimsizes(u&lon)/),float)
  copy_VarCoords(u(0,:,:), lluw(0,:,:))
  lluw = 0.

  llvw = new((/4,dimsizes(u&lat),dimsizes(u&lon)/),float)
  copy_VarCoords(u(0,:,:), llvw(0,:,:))
  llvw = 0.

  prem = 25000./100000.
  nlati = ind(lat.ge.-20.and.lat.le.10)

  nt = ind(phase.eq.8.and.amp.gt.1.)
  nind = ind_DJF(nt)
  do j = 0, 3
    ntt   = -10 + j*5
    nindd = nind + ntt
    llphi(j,:,:) = dim_avg_n(phi(nindd,:,:), 0)
    llu(j,:,:)   = dim_avg_n(u(nindd,:,:), 0)
    llolr(j,:,:) = dim_avg_n(olr(nindd,:,:), 0)
    ;um = dim_avg_n(u0(nindd,:,:), 0)
    ;vm = dim_avg_n(v0(nindd,:,:), 0)
    WAFuv = cal_WFsf(llphi(j,:,:),um,vm,phi&lat,phi&lon)
    WAFuv(:,nlati,:) = 0.
    lluw(j,:,:) = WAFuv(0,:,:)*prem
    llvw(j,:,:) = WAFuv(1,:,:)*prem
    delete([/ntt,nindd,WAFuv/])
  end do
  delete([/nt,nind/])

  llvalue = sqrt(lluw*lluw+llvw*llvw)

  lluw = mask(lluw,llvalue.gt.0.1,True)
  llvw = mask(llvw,llvalue.gt.0.1,True)

  llphi = llphi*10^(-6)
  llolr = smth9_Wrap(llolr, 0.5, 0.25, False)

  print("complete calc leadlag data. "+ systemfunc("date"))
;----------------------------------------------------------------------
;; resources
  undef("cmap")
  cmap = read_colormap_file("BlueWhiteOrangeRed")
;shading u plot
  undef("rest")
  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = False
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 30
  rest@gsnRightStringFontHeightF = 0.025 
  rest@gsnLeftStringFontHeightF  = 0.025
  rest@gsnRightStringOrthogonalPosF = 0.01
  rest@gsnLeftStringOrthogonalPosF  = 0.01
  rest@gsnRightString = "~F21~250hPa_SF_WAF"


  rest@tmXTOn = False
  rest@tmYROn = False
  rest@tmXBLabelFont = 21
  rest@tmYLLabelFont = 21
  rest@tmXBLabelFontHeightF = 0.02
  rest@tmYLLabelFontHeightF = 0.02
  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
  rest@mpMinLatF   = lats
  rest@mpMaxLatF   = latn-5
  rest@mpMinLonF   = lonw+10
  rest@mpMaxLonF   = lone-10
  rest@mpLandFillColor = "white"

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap;(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -6.
  rest@cnMaxLevelValF         =  6.
  rest@cnLevelSpacingF        =  1.
  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = False

; contour olr plot.
  resolr = True
  resolr@gsnDraw      = False
  resolr@gsnMaximize  = False
  resolr@gsnFrame     = False
  resolr@gsnAddCyclic = False
  resolr@gsnContourNegLineDashPattern = 1
  resolr@cnMonoLineColor    = True
  resolr@cnLinesOn          = True
  resolr@cnLineLabelsOn     = False
  resolr@cnLineColor        = "red3"
  resolr@cnLineThicknessF   = 3.0
  resolr@cnInfoLabelPerimOn = False
  resolr@cnInfoLabelOn      = False
  resolr@cnLevelSelectionMode = "Explicitlevels"
  resolr@cnLevels      = (/-25,-20,-15,-10,-5,5,10,15,20,25/)

;; vector res.
  undef("resw")
  resw = True
  resw@gsnFrame     = False
  resw@gsnDraw      = False
  resw@gsnMaximize  = False
  resw@gsnAddCyclic = False
  resw@vcGlyphStyle            = "LineArrow"
  resw@vcLineArrowThicknessF   = 3.0
  resw@vcRefMagnitudeF         = 4
  resw@vcRefLengthF            = 0.05
  resw@vcLineArrowColor        = "black"
  resw@vcRefAnnoOn             = False
  resw@vcRefAnnoString2On      = False
  resw@vcRefAnnoString1On      = False
  resw@vcRefAnnoPerimSpaceF    = 0.5
  resw@vcRefAnnoPerimOn        = False
  resw@vcRefAnnoFontHeightF    = 0.015
  resw@vcRefAnnoJust           = "BottomRight"
  resw@vcRefAnnoBackgroundColor = -1

  resP = True
  resP@gsnMaximize      = True
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02
  resP@gsnPanelLeft     = 0.5
  resP@gsnPanelYWhiteSpacePercent = 3

  resP@pmLabelBarOrthogonalPosF =  0.012;-0.005
  resP@pmLabelBarWidthF         =  0.4
  resP@pmLabelBarParallelPosF   =  0.005;-0.1

  resP@lbLabelFont       = 21
  resP@lbBoxLinesOn      = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012

  undef("sres")
  sres   = True 
  sres@cnLineDrawOrder    = "PostDraw"                    
  sres@gsnDraw            = False                          
  sres@gsnFrame           = False                        
  sres@cnLineLabelsOn     = False                              
  sres@cnLinesOn          = False                           
  sres@cnInfoLabelOn      = False
  sres@cnFillOn           = False                       
  sres@cnInfoLabelOn      = False
  sres@cnFillScaleF       = 0.1
  sres@cnFillDotSizeF     = 0.0025
  sres@gsnLeftString      =" "
  sres@gsnRightString     =" "
  sres@cnMonoFillColor    = True
  sres@cnFillColor        = "gray"

  undef("opt")
  opt = True
  opt@gsnShadeFillType   = "pattern"
  opt@gsnShadeHigh       = 17
  opt@gsnShadeFillScaleF = 1.2

  undef("respoly")
  respoly = True
  respoly@gsLineThicknessF  = 3.5
  respoly@gsLineDashPattern = 2
  respoly@gsLineColor       = "purple"

  print("complete res. "+ systemfunc("date"))
;----------------------------------------------------------------------
;; plot.
  undef("plott")
  plott   = new(4,graphic)
  plotolr = plott
  plotvc  = plott
  undef("dum")
  dum     = plott
  undef("dum2")
  dum2    = plott
  undef("Label1")
  Label1 = (/"~F21~(c) day-10_P8","~F21~(d) day-5_P8",\
             "~F21~(e) day0_P8","~F21~(f) day5_P8"/)

  undef("x")
  undef("xt")
  do n = 0,3
    rest@gsnLeftString = "~F21~"+Label1(n)
    if (n.eq.3)then
      resw@vcRefAnnoOn             = True
      resw@vcRefAnnoString1On      = True
      resw@vcRefAnnoString1        = "$VMG$" +" ~F21~m~S~2~N~ s~S~-2~N~"
      resw@vcRefAnnoOrthogonalPosF = 0.1
      resw@vcRefAnnoParallelPosF   = -0.05
      resw@vcRefAnnoFontHeightF    = 0.02
    end if

    x  = llphi(n,:,:)*10^(6)
    xt = llphi(n,:,:)
    xt = 0.
    do i = 0, dimsizes(lat)-1
      do j = 0, dimsizes(lon)-1
        if(x(i,j).lt.mvalue(0,0,n,i,j).or.x(i,j).gt.mvalue(0,1,n,i,j)) then
          xt(i,j) = 2.
        else
          xt(i,j) = -2.
        end if
      end do
    end do

    plott(n)   = gsn_csm_contour_map(wks, llphi(n,:,:), rest)
    plotolr(n) = gsn_csm_contour(wks, llolr(n,:,:), resolr)
    plotvc(n)  = gsn_csm_vector(wks, lluw(n,::3,::3), llvw(n,::3,::3), resw)
    plotu  = gsn_csm_contour(wks, uc, resc)
    dum(n)  = gsn_csm_contour(wks,xt,sres)
    dum(n)  = gsn_contour_shade(dum(n),1.9,-999.,opt)
    dum2(n) = gsn_add_polyline(wks, plott(n), (/-110,-40,-40,-110,-110/), \
                              (/-90,-90,90,90,-90/), respoly)

    overlay(plott(n),dum2(n))
    overlay(plott(n),dum(n))
    overlay(plott(n),plotvc(n))
    overlay(plott(n),plotu)
    overlay(plott(n),plotolr(n))

    delete([/x,xt/])

  end do

  resP@lbTitleOn = True
  resP@lbTitleString = "[10~S~6~N~ m~S~2~N~ s~S~-1~N~]"
  resP@lbTitleFontHeightF = 0.013
  resP@lbTitlePosition = "Right"
  resP@lbTitleDirection = "Across"
  gsn_panel(wks,(/plott(0),plott(1),\
                  plott(2),plott(3)/),(/1,1,1,1/),resP)

  ;i = 1
  ;do while (i.eq.1)
  ;  i = 1
  ;end do