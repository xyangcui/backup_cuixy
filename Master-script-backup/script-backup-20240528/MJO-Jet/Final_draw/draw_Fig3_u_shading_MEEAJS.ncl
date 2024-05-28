;--------------------------------------------------------
; draw uwind anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour. 
; This is for test sig. 0.05 by monte carlo method.
; time: 2024-05-15

;--------------------------------------------------------
;function: calculate sig values.
;--------------------------------------------------------
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
;-----------------------------------------------------------------------------
  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  udir  = ddir+"/"+"ERA5_daily_anom_u250_1979-2022_2.5x2.5.nc"
  ucdir = ddir+"/"+"ERA5_daily_u250_1979-2022_r2.5.nc"
  sdir  = ddir+"/"+"monte_carlo/Fig3_u_sig.nc"

  lats = 10
  latn = 50
  lonw = 0
  lone = 210
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
  resc@cnLevels      = (/30.,45.,60.,75./)

  delete([/time,f,TIME,mon,nt/])

  print("draw clim.u complete. "+ systemfunc("date"))
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------
  fu = addfile(udir,"r")

  u  = fu->u(:,{lats:latn},{lonw:lone})

  time = u&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  print("complete data selling "+ systemfunc("date"))

  fs = addfile(sdir,"r")
  mvalue1 = fs->mvalue;(:,:,{lats:latn},{lonw:lone})
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
  llu = new((/4,dimsizes(u&lat),dimsizes(u&lon)/),float)
  copy_VarCoords(u(0,:,:), llu(0,:,:))
  llu = 0.

  nt = ind(phase.eq.2.and.amp.gt.1.)
  nind = ind_DJF(nt)
  ntt  = 0
  nindd = nind + ntt
  do ntt = 1,4
    nit = nind + ntt
    nindd := array_append_record(nindd, nit, 0)
    delete(nit)
  end do
  llu(0,:,:) = dim_avg_n(u(nindd,:,:), 0)
  delete([/nt,nind,ntt,nindd/])

  nt = ind(phase.eq.6.and.amp.gt.1.)
  nind = ind_DJF(nt)
  ntt   = 0
  nindd = nind + ntt
  do ntt = 1,4
    nit = nind + ntt
    nindd := array_append_record(nindd, nit, 0)
    delete(nit)
  end do
  llu(1,:,:) = dim_avg_n(u(nindd,:,:), 0)
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

  resu@mpMinLatF = lats
  resu@mpMaxLatF = latn
  resu@mpMinLonF = lonw
  resu@mpMaxLonF = lone
  resu@mpFillOn        = False
  resu@mpLandFillColor = "white"
  resu@mpGridAndLimbDrawOrder = "PostDraw"
  resu@mpGridLineColor        = "Gray60"

  resu@cnFillOn               = True
  resu@cnFillPalette          = cmap(::-1,:);(66:146,:)
  resu@cnLinesOn              = False
  resu@cnLineLabelsOn         = False
  resu@cnLevelSelectionMode   = "Manuallevels"
  resu@cnMinLevelValF         = -5.
  resu@cnMaxLevelValF         =  5.
  resu@cnLevelSpacingF        =  1.
  resu@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  resu@lbLabelBarOn           = False

  resP = True
  resP@gsnMaximize  = True
  resP@gsnFrame     = False
  resP@gsnDraw      = True

  resP@pmLabelBarOrthogonalPosF = 0.0
  resP@pmLabelBarParallelPosF   = 0.0
  resP@pmLabelBarWidthF         = 0.6

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
;----------------------------------------------------------------------
;; plot.
  wks = gsn_open_wks("x11", pdir+"/"+"Fig3_U_MEEAJS")

  plot = new(2,graphic)
  dum1 = plot
  Label  = (/"(a) 250hPa_U@P2","(b) 250hPa_U@P6"/)

  do n = 0,1
    x  = llu(n,:,:)
    mvalue = mvalue1(:,n,:,:)
    xt = siglevel(x,mvalue)
    
    dum1(n) = gsn_csm_contour(wks,xt,sres)
    dum1(n) = gsn_contour_shade(dum1(n),1.9,-999.,opt)

    resu@gsnLeftString  = "~F21~"+Label(n)
    plot(n) = gsn_csm_contour_map(wks, x, resu)
    plotut  = gsn_csm_contour(wks, uc, resc)

    lonmin =  10
    lonmax =  70
    latmin = -100
    latmax =  100
    polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
    polyy = (/latmin,latmin,latmax,latmax,latmin/)
    respoly@gsLineDashPattern = 0
    dum2 = gsn_add_polyline(wks, plot(n), polyx, polyy, respoly)

    lonmin =  110
    lonmax =  180
    latmin = -100
    latmax =  100
    polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
    polyy = (/latmin,latmin,latmax,latmax,latmin/)
    respoly@gsLineDashPattern = 1
    dum3 = gsn_add_polyline(wks, plot(n), polyx, polyy, respoly)

    overlay(plot(n),plotut)
    overlay(plot(n),dum1(n))
    overlay(plot(n),dum2)
    overlay(plot(n),dum3)
  end do

  resP@gsnPanelLabelBar = True
  resP@lbTitleOn = True
  resP@lbTitleString = "[m/s]"
  resP@lbTitlePosition  = "Right"
  resP@lbTitleDirection = "Across"
  resP@lbTitleFontHeightF = 0.015
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.013
  resP@lbTitleStringParallelPosF = -0.1

  gsn_panel(wks,(/plot(0),plot(1)/),(/1,1/),resP);

  i = 0
  do while(i.eq.0)
    i= 0
  end do