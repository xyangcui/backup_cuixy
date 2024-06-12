;--------------------------------------------------------
; draw uwind anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour. 
; waf vector
; This is for test sig. 95% by monte carlo method.
; time: 2023-12-18
; 2.5x2.5.
;--------------------------------------------------------

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

  olrdir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
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
  folr = addfile(olrdir,"r")
  ft = addfile(tdir,"r")

  t    = ft->t(:,{lats:latn},{lonw:lone})
  olr  = folr->olr(:,{lats:latn},{lonw:lone})
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
;; sig
  ;fsig = addfile(sigdir,"r")
  ;m_value = fsig->m_value

;----------------------------------------------------------------------
;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
  llphi = new((/2,dimsizes(phi&lat),dimsizes(phi&lon)/),float)
  copy_VarCoords(phi(0,:,:), llphi(0,:,:))
  llphi = 0.

  llu = new((/2,dimsizes(u&lat),dimsizes(u&lon)/),float)
  copy_VarCoords(u(0,:,:), llu(0,:,:))
  llu = 0.

  llv = new((/2,dimsizes(v&lat),dimsizes(v&lon)/),float)
  copy_VarCoords(v(0,:,:), llv(0,:,:))
  llv = 0.

  llt = new((/2,dimsizes(t&lat),dimsizes(t&lon)/),float)
  copy_VarCoords(t(0,:,:), llt(0,:,:))
  llt = 0.

  nt = ind(phase.eq.8.and.amp.gt.1.)
  nind = ind_DJF(nt)
  ntt  = 6
  nindd = nind + ntt
  do ntt = 7,15
    nit = nind + ntt
    nindd := array_append_record(nindd, nit, 0)
    delete(nit)
  end do

  llphi(1,:,:) = dim_avg_n(phi(nindd,:,:), 0)
  llu(1,:,:)   = dim_avg_n(u(nindd,:,:), 0)
  llv(1,:,:)   = dim_avg_n(v(nindd,:,:), 0)
  llt(1,:,:)   = dim_avg_n(t(nindd,:,:), 0)
  delete([/nt,nind,ntt,nindd/])

  nt = ind(phase.eq.4.and.amp.gt.1.)
  nind = ind_DJF(nt)
  ntt   = 1
  nindd = nind + ntt
  do ntt = 2,10
    nit = nind + ntt
    nindd := array_append_record(nindd, nit, 0)
    delete(nit)
  end do

  llphi(0,:,:) = dim_avg_n(phi(nindd,:,:), 0)
  llu(0,:,:)   = dim_avg_n(u(nindd,:,:), 0)
  llv(0,:,:)   = dim_avg_n(v(nindd,:,:), 0)
  llt(0,:,:)   = dim_avg_n(t(nindd,:,:), 0)
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
  resw@vcLineArrowColor        = "green3"
  resw@vcRefAnnoString2On      = False
  resw@vcRefAnnoString1        = "$VMG$" +" m s~S~-1~N~"
  resw@vcRefAnnoPerimOn        = False
  resw@vcRefAnnoFontHeightF    = 0.02
  resw@vcRefAnnoBackgroundColor = -1

  resP = True
  resP@gsnMaximize      = True
  resP@gsnFrame         = False
  resP@gsnDraw          = True

  resP@pmLabelBarOrthogonalPosF = 0.0
  resP@pmLabelBarParallelPosF   = 0.0
  resP@pmLabelBarWidthF         = 0.5

  resP@lbTitleFontHeightF= 0.012
  resP@lbBoxLinesOn      = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012

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
  respoly@gsLineColor      = "purple"
  respoly@gsLineThicknessF = 4.5

  print("complete res. "+ systemfunc("date"))

  lonmin = -110
  lonmax = -40
  latmin = -100
  latmax =  100
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
;----------------------------------------------------------------------
  wks = gsn_open_wks("x11", pdir+"/"+"Fig6_U_SF_T_NAJS")

  plott   = new(2,graphic)
  plotphi = plott

  plotu   = plott
  plotu2  = plott
  plotu3  = plott

  plotv   = plott

  dum1    = plott
  dum2    = plott
  dum3    = plott

  Label  = (/"~F21~(a) 250hPa_U@P4","~F21~(d) 250hPa_U@P8"/)
  Label1 = (/"~F21~(b) 250hPa_SF@P4","~F21~(e) 250hPa_SF@P8"/)
  Label2 = (/"~F21~(c) Temp@P4","~F21~(f) Temp@P8"/)

  xf = (/0.05,0.55/)

  do n = 0,1

    x  = llphi(n,:,:)
    y  = llu(n,:,:)
    z  = llt(n,:,:)
    l  = llv(n,:,:)

    if(n.eq.0)then
      mvalueu = mvalue(0,0,:,:,:)
      mvaluev = mvalue(1,0,:,:,:)
      mvalues = mvalue(2,0,:,:,:)
      mvaluet = mvalue(3,0,:,:,:)

      xt = siglevel(x,mvalues)
      yt1 = siglevel(y,mvalueu)
      zt = siglevel(z,mvaluet)

      lt = mask(l,y.lt.mvalueu(0,:,:).or.y.gt.mvalueu(0,:,:).or.\
                l.lt.mvaluev(0,:,:).or.l.gt.mvaluev(1,:,:),True)
      yt = mask(y,y.lt.mvalueu(0,:,:).or.y.gt.mvalueu(0,:,:).or.\
                l.lt.mvaluev(0,:,:).or.l.gt.mvaluev(1,:,:),True)

      copy_VarCoords(y,lt)
      copy_VarCoords(y,yt)

      resu@lbLabelBarOn = False
      resw@vcRefAnnoOn  = False

    else
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
      resu@pmLabelBarWidthF  = 0.5
      resu@pmLabelBarHeightF = 0.05
      resu@pmLabelBarParallelPosF   = -0.15
      resu@pmLabelBarOrthogonalPosF = 0.12
      resw@vcRefAnnoOn = True
      resw@vcRefAnnoOrthogonalPosF = -0.2
      resw@vcRefAnnoParallelPosF   =  0.6
    end if
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

    plotu(n) = gsn_csm_contour_map(wks, y, resu)
    dum2(n)  = gsn_csm_contour(wks,yt1,sres)
    dum2(n)  = gsn_contour_shade(dum2(n),1.9,-999.,opt)
    overlay(plotu(n),dum2(n))
;; draw phi.
    resu@vpHeightF = 0.45
    resu@vpXF = xf(n)
    resu@vpYF = 0.75

    resu@mpMinLatF   =  0
    resu@mpMaxLatF   =  60
    resu@mpMinLonF   = -120
    resu@mpMaxLonF   = -30
    resu@gsnLeftString  = "~F21~"+Label1(n)

    resu@lbTitleString = "[10~S~6~N~ m~S~2~N~ s~S~-1~N~]"

    resu@cnMinLevelValF    = -4.
    resu@cnMaxLevelValF    =  6.
    resu@cnLevelSpacingF   =  1.

    plotu2(n) = gsn_csm_contour_map(wks, x, resu)
    plotv(n)  = gsn_csm_vector(wks, yt(::2,::2), lt(::2,::2), resw)
    dum1(n)   = gsn_csm_contour(wks,xt,sres)
    dum1(n)   = gsn_contour_shade(dum1(n),1.9,-999.,opt)
    overlay(plotu2(n),dum1(n))
    overlay(plotu2(n),plotv(n))
; temp.
    resu@vpHeightF = 0.25
    resu@vpXF = xf(n)
    resu@vpYF = 0.3

    resu@mpMinLatF   =  10
    resu@mpMaxLatF   =  55
    resu@mpMinLonF   = -120
    resu@mpMaxLonF   = -30
    resu@gsnLeftString  = "~F21~"+Label2(n)

    resu@lbTitleString = "[K]"

    resu@cnMinLevelValF  = -0.7
    resu@cnMaxLevelValF  =  0.6
    resu@cnLevelSpacingF =  0.05

    plotu3(n) = gsn_csm_contour_map(wks, z, resu)
    dum3(n)   = gsn_csm_contour(wks,zt,sres)
    dum3(n)   = gsn_contour_shade(dum3(n),1.9,-999.,opt)
    overlay(plotu3(n),dum3(n))

    plotut  = gsn_csm_contour(wks, uc, resc)
    plotut1 = gsn_csm_contour(wks, uc, resc)
    plotut2 = gsn_csm_contour(wks, uc, resc)
    overlay(plotu(n),plotut)
    overlay(plotu2(n),plotut1)
    overlay(plotu3(n),plotut2)

    dum4 = gsn_add_polyline(wks, plotu(n),  polyx, polyy, respoly)
    dum5 = gsn_add_polyline(wks, plotu2(n), polyx, polyy, respoly)
    dum6 = gsn_add_polyline(wks, plotu3(n), polyx, polyy, respoly)
    overlay(plotu(n) ,dum4)
    overlay(plotu2(n),dum5)
    overlay(plotu3(n),dum6)

  end do

  respoly = True
  respoly@gsLineThicknessF = 3.5
  respoly@gsLineColor      = "gold"

  respoly@gsLineDashPattern = 0
  dumea = gsn_add_polyline(wks, plotu3(1), (/-105,-75,-75,-105,-105/), \
                          (/35,35,50,50,35/), respoly)

  respoly@gsLineDashPattern = 0
  dumea2 = gsn_add_polyline(wks,plotu3(0), (/-105,-75,-75,-105,-105/), \
                          (/35,35,50,50,35/), respoly)

  draw(plotu(0))
  draw(plotu(1))
  draw(plotu2(0))
  draw(plotu2(1))
  draw(plotu3(0))
  draw(plotu3(1))

  i = 0
  do while(i.eq.0)
    i= 0
  end do