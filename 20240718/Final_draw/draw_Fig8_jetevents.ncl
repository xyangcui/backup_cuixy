;----------------------------------------
;; select Jet events by PC1.
;; "candidate": consecutive 1 pentad 1.0.
;; Pre 12 days must grow and decay.
;; remove each year pre 12 days and late.
;----------------------------------------

;;EOF plot.
  lats =  10
  latn =  80
  lonw = -60
  lone =  30

  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  edir = ddir+"/EOF_ulevel.nc"

  fe = addfile(edir,"r")
  eof1 = fe->u(0,:,:,:)

  uv = dim_avg_n_Wrap(eof1,2)

  level = eof1&level
  lat = eof1&lat
  lon = eof1&lon

  cmap = read_colormap_file("MPL_RdBu")

  res              = True         
  res@gsnDraw      = False        
  res@gsnFrame     = False       
  res@gsnAddCyclic = False  
  res@gsnRightStringFontHeightF = 0.015
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnRightStringOrthogonalPosF = 0.01
  res@gsnLeftStringOrthogonalPosF  = 0.01 
  res@tmYLLabelFontHeightF = 0.012     
  res@tmXBLabelFontHeightF = 0.012 
  res@gsnMajorLatSpacing = 15
  ;res@gsnMajorLonSpacing= 

  res@mpFillOn     = False        
  res@mpMinLatF    = lats        
  res@mpMaxLatF    = latn
  res@mpMinLonF    = lonw
  res@mpMaxLonF    = lone
  res@cnFillOn      = True                
  res@cnFillPalette = cmap(::-1,:)
  res@cnLevelSelectionMode   = "Manuallevels"

  ;res@cnMinLevelValF         = -1.
  ;res@cnMaxLevelValF         =  1.
  ;res@cnLevelSpacingF        =  0.1 
  symMinMaxPlt(eof1, 11, False, res) 

  res@pmLabelBarOrthogonalPosF = 0.15

  resv   = True               
  resv@gsnDraw      = False        
  resv@gsnFrame     = False       
  resv@gsnAddCyclic = False
  resv@gsnRightStringFontHeightF = 0.015
  resv@gsnLeftStringFontHeightF  = 0.015
  resv@gsnRightStringOrthogonalPosF = 0.01
  resv@gsnLeftStringOrthogonalPosF  = 0.01 
  resv@tmYLLabelFontHeightF = 0.012     
  resv@tmXBLabelFontHeightF = 0.012 
  resv@gsnMajorLatSpacing = 15
  resv@cnLevelSelectionMode = "ManualLevels"      
  ;resv@cnLevelSpacingF      = 5.0                  
  ;resv@cnMinLevelValF       = -50.                 
  ;resv@cnMaxLevelValF       =  50.                 
  resv@cnLineLabelsOn       = False                 
  resv@cnFillOn             = True                 
  resv@cnFillPalette        = cmap(::-1,:)
  resv@tmYRMode             = "Automatic"  
  symMinMaxPlt(uv, 11, False, resv)  


  wks = gsn_open_wks("x11","eof") 


  res@vpXF = 0.15
  res@vpYF = 0.98
  res@vpWidthF  = 0.3
  res@vpHeightF = 0.3
  res@gsnLeftString = "~F21~(a) EOF1_250hPa"
  plot1 = gsn_csm_contour_map(wks,eof1({250},:,:),res)

  resv@vpXF = 0.55
  resv@vpYF = 0.95
  resv@vpWidthF  = 0.3
  resv@vpHeightF = 0.25
  resv@gsnLeftString = "~F21~(b) EOF1_levels"
  plot2 = gsn_csm_pres_hgt(wks, uv, resv )

  draw(plot1)
  draw(plot2)

  ;i = 0
  ;do while (i.eq.0)
  ;  i = 0
  ;end do

  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  indir = ddir+"/centered_date.csv"

  pdir = ddir+"/Ulevel_PC"
; time table.
  idir = ddir + "/ERA5_daily_mean_ulevel_1979-2022_r2.5.nc"
  f = addfile(idir,"r")

  u = f->u(:,{250},{lats:latn},{lonw:lone})
  printVarSummary(u)
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

;; selling events.
  delim  = ","
  rmmf   = asciiread(indir,-1,"string")
  iyear  = str_get_field(rmmf,1,delim)
  imon   = str_get_field(rmmf,2,delim)
  iday   = str_get_field(rmmf,3,delim)
  pct = str_get_field(rmmf,4,delim)

  iyear  := toint(iyear(1::))
  imon   := toint(imon(1::))
  iday   := toint(iday(1::))
  pct := toint(pct(1::))

  pev   = ind(pct.gt.0)
  pyear = iyear(pev)
  pmon  = imon(pev)
  pday  = iday(pev)

  pu = new((/dimsizes(pev),11,dimsizes(u&lat),dimsizes(u&lon)/),float)
  pu = 0.

  do i = 0, dimsizes(pev)-1
    nt = ind(year.eq.pyear(i).and.mon.eq.pmon(i).and.day.eq.pday(i))
    pu(i,:,:,:) = u(nt-5:nt+5,:,:)
    delete(nt)
  end do

  nev   = ind(pct.lt.0)
  nyear = iyear(nev)
  nmon  = imon(nev)
  nday  = iday(nev)

  nu = new((/dimsizes(nev),11,dimsizes(u&lat),dimsizes(u&lon)/),float)
  nu = 0.

  do i = 0, dimsizes(nev)-1
    nt = ind(year.eq.nyear(i).and.mon.eq.nmon(i).and.day.eq.nday(i))
    nu(i,:,:,:) = u(nt-5:nt+5,:,:)
    delete(nt)
  end do

  copy_VarCoords(u(0,:,:), pu(0,0,:,:))
  copy_VarCoords(u(0,:,:), nu(0,0,:,:))

  pu := dim_avg_n_Wrap(pu, 0)
  nu := dim_avg_n_Wrap(nu, 0)

;;plot
  plotp = new(3,graphic)                
  plotn = new(3,graphic)                

  undef("res")
  res              = True         
  res@gsnDraw      = False        
  res@gsnFrame     = False       
  res@gsnAddCyclic = False  
  res@gsnRightStringFontHeightF = 0.03
  res@gsnLeftStringFontHeightF  = 0.03
  res@gsnRightStringOrthogonalPosF = 0.01
  res@gsnLeftStringOrthogonalPosF  = 0.01   
  res@tmYLLabelFontHeightF = 0.025     
  res@tmXBLabelFontHeightF = 0.025
  res@gsnMajorLatSpacing = 15
  ;res@gsnMajorLonSpacing= 

  res@mpFillOn     = False        
  res@mpMinLatF    = lats        
  res@mpMaxLatF    = latn
  res@mpMinLonF    = lonw
  res@mpMaxLonF    = lone
                                       
  res = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False
  res@gsnAddCyclic = False
  res@gsnContourNegLineDashPattern = 1
  res@gsnLeftString      = ""
  res@gsnRightString     = ""
  res@cnMonoLineColor    = True
  res@cnLinesOn          = True
  res@cnLineLabelsOn     = False
  res@cnLineColor        = "black"
  res@cnLineThicknessF   = 2.0
  res@cnInfoLabelPerimOn = False
  res@cnInfoLabelOn      = False
  res@cnLevelSelectionMode = "Explicitlevels"
  res@cnLevels   = (/20.,25.,30.,35./)

  resP                    = True         
  resP@gsnMaximize        = True
  resP@gsnPanelYWhiteSpacePercent = 3     
  resP@gsnPanelTop = 0.63           

  label = (/"(c) pDay -5","(d) pDay 0","(e) pDay 5",\
            "(f) nDay -5","(g) nDay 0","(h) nDay 5"/)

  res@gsnLeftString = "~F21~"+label(0)
  plotp(0)= gsn_csm_contour_map(wks,pu(0,:,:),res)
  res@gsnLeftString = "~F21~"+label(1)
  plotp(1)= gsn_csm_contour_map(wks,pu(6,:,:),res)
  res@gsnLeftString = "~F21~"+label(2)
  plotp(2)= gsn_csm_contour_map(wks,pu(10,:,:),res)

  res@gsnLeftString = "~F21~"+label(3)
  plotn(0)= gsn_csm_contour_map(wks,nu(0,:,:),res)
  res@gsnLeftString = "~F21~"+label(4)
  plotn(1)= gsn_csm_contour_map(wks,nu(6,:,:),res)
  res@gsnLeftString = "~F21~"+label(5)
  plotn(2)= gsn_csm_contour_map(wks,nu(10,:,:),res)

  gsn_panel(wks,(/plotp(0),plotp(1),plotp(2),\
                  plotn(0),plotn(1),plotn(2)/),(/2,3/),resP)    

  i = 0
  do while (i.eq.0)
    i = 0
  end do