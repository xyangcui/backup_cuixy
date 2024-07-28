;-------------------------------------------
;; draw MJO comp. u250 for three SWJ.
;; MEJS and EAJS (0-4 d) NAJS (10-14 d)
;; Information
;; time: 2024-05-13
;; version: v1
;-------------------------------------------

  lats = 10
  latn = 50
  lonw = 100
  lone = 200

  ddir   = "/home/sunming/data5/cuixy/ENSO-MJO-Jet/data"
  pdir   = "/home/sunming/data5/cuixy/ENSO-MJO-Jet/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
  ninodir = "/home/sunming/data5/cuixy/data/clim_indices/Nino_3.4_index.csv"

  udir   = ddir+"/"+"ERA5_u250_daily_anom_1979-2022_r1.nc"
  ucdir  = ddir+"/"+"u_DJF_mean_1979-2022_r2.5.nc"

;--------------------------------------------------------
;; Main code.
;--------------------------------------------------------

  f = addfile(udir,"r")

  u = f->u(:,{lats:latn},{lonw:lone})

  time = u&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  fuc = addfile(ucdir,"r")
  uc  = fuc->u(0,0,:,:)

;; MJO indices.
;-------------------------------------------------------
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

  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.2.and.day.eq.28)

  amp   := amp(n1:n2)
  phase := phase(n1:n2)
  mon   := mon(n1:n2)
  time  := time(n1:n2)
  year  := year(n1:n2)

  nt     = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  amp   := amp(nt)
  phase := phase(nt)
  year  := year(nt)
  mon   := mon(nt)

  delete(nt)
  ind_DJF = amp&timeind 

;; select for EI, LA or NE. Dec year + 1.
  delim  = ","
  ninof  = asciiread(ninodir,-1,"string")
  iyear  = str_get_field(ninof,1,delim)
  ievent = str_get_field(ninof,6,delim)

  iyear  := toint(iyear(1::))
  ievent := tofloat(ievent(1::))

;;EI
  nei = ind(ievent.eq.1)
  yei = iyear(nei)
;;LA
  nla = ind(ievent.eq.-1)
  yla = iyear(nla)
;;NE
  nne = ind(ievent.eq.0)
  yne = iyear(nne)

  nt = ind(mon.eq.12)
  year(nt) = year(nt)+1
  delete(nt)

;;for EI.
  ny = yla
  nc = ind(year.eq.ny(0))
  do i = 1, dimsizes(ny)-1
    ntemp = ind(year.eq.ny(i))
    nc := array_append_record(nc, ntemp, 0)
    delete(ntemp)
  end do
  indt   = ind_DJF(nc)
  phaset = phase(nc)
  ampt   = amp(nc)
;------------------------------------------------------
;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
  llu = new((/8,dimsizes(u&lat),dimsizes(u&lon)/),float)
  copy_VarCoords(u(0,:,:), llu(0,:,:))
  llu = 0.

  do i = 1,8

    nt = ind(phaset.eq.i.and.ampt.gt.1.)
    nind = indt(nt)
    ntt   = 0
    nindd = nind + ntt

    do ntt = 1,4
      nit = nind + ntt
      nindd := array_append_record(nindd, nit, 0)
      delete(nit)
    end do

    llu(i-1,:,:) = dim_avg_n(u(nindd,:,:), 0)
    delete([/nt,nind,ntt,nindd/])

  end do

  print("complete calc leadlag data. "+ systemfunc("date"))

;------------------------------------------------------------
;;plot
;; resources

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
  resc@cnLineColor        = "black";magenta
  resc@cnLineThicknessF   = 2.5
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels = (/30.,40.,46.,60.,70./)

  cmap = read_colormap_file("MPL_RdBu") 
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

  resu@mpMinLatF  = lats
  resu@mpMaxLatF  = latn
  resu@mpMinLonF  = lonw
  resu@mpMaxLonF  = lone
  resu@mpFillOn        = False
  resu@mpLandFillColor = "white"
  resu@mpGridAndLimbDrawOrder = "PostDraw"
  resu@mpGridLineColor        = "Gray60"

  resu@cnFillOn               = True
  resu@cnFillPalette          = cmap(::-1,:)
  resu@cnLinesOn              = False
  resu@cnLineLabelsOn         = False
  resu@cnLevelSelectionMode   = "Manuallevels"
  resu@cnMinLevelValF         = -6.
  resu@cnMaxLevelValF         =  6.
  resu@cnLevelSpacingF        =  1.
  resu@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  resu@lbLabelBarOn           = False

  resP = True
  resP@gsnMaximize  = True
  resP@gsnFrame     = False
  resP@gsnDraw      = True

  resP@pmLabelBarOrthogonalPosF = 0.0
  resP@pmLabelBarParallelPosF   = 0.0
  resP@pmLabelBarWidthF         = 0.5

  resP@lbTitleFontHeightF = 0.012
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012

  lonmin =  110
  lonmax =  180
  latmin = -100
  latmax =  100
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)

  respoly = True
  respoly@gsLineColor      = "purple"
  respoly@gsLineThicknessF = 4.5

  print("complete res. "+ systemfunc("date"))

  wks = gsn_open_wks("x11", pdir+"/"+"Fig_comp.U_MEJS")

  plotu = new(8,graphic)

  resu@gsnRightString = "250hPa u"

  label = (/"(a) phase1","(b) phase2","(c) phase3","(d) phase4",\
            "(e) phase5","(f) phase6","(g) phase7","(h) phase8"/)

  do i = 0,7
    resu@gsnLeftString = label(i)
    plotu(i) = gsn_csm_contour_map(wks, llu(i,:,:), resu)
    plotuc = gsn_csm_contour(wks,uc,resc)
    dum = gsn_add_polyline(wks, plotu(i), polyx, polyy, respoly)
    overlay(plotu(i),dum)
    overlay(plotu(i),plotuc)
  end do

  resP@gsnPanelLabelBar = True
  gsn_panel(wks,(/plotu(0),plotu(4),\
                  plotu(1),plotu(5),\
                  plotu(2),plotu(6),\
                  plotu(3),plotu(7)/),(/2,2,2,2/),resP)

  i = 0
  do while (i.eq.0)
    i = 0
  end do