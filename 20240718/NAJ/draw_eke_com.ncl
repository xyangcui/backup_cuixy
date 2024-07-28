;----------------------------------------
;; select Jet events by PC1.
;; "candidate": consecutive 1 pentad 1.0.
;; Pre 12 days must grow and decay.
;; remove each year pre 12 days and late.
;----------------------------------------
  lats =  15
  latn =  80
  lonw = -120
  lone =  10

  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  indir = ddir+"/mhf_events.csv";"/centered_date.csv"

  pdir = ddir+"/Ulevel_PC"
; time table.
  idir = ddir + "/EKE250_daily_anom_1979-2022_r2.5.nc"
  f = addfile(idir,"r")

  u = f->EKE(:,{lats:latn},{lonw:lone})
  printVarSummary(u)
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.4.and.day.eq.1)

  mont = mon(n1:n2)
  ut = u(n1:n2,:,:)
  nt = ind(mont.eq.12.or.mont.eq.1.or.mont.eq.2)
  ut := ut(nt,{45:55},{-85:-55})

  ;rad  = 4.0*atan(1.0)/180.0
  ;clat = cos(ut&lat * rad)
  ;mind = wgt_areaave_Wrap(ut, clat, 1.0, 1)

  ;mind := dim_standardize_n_Wrap(mind, 0, 0)

  ;filename = ddir+"/eke_index.csv"
  ;system("rm -rf " + filename)
  ;write_table(filename,"w",[/mind/],"%5.2f")

  delete([/n1,n2,nt,ut/])
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
    pu(i,:,:,:) = u(nt-4:nt+6,:,:)
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
    nu(i,:,:,:) = u(nt-4:nt+6,:,:)
    delete(nt)
  end do

  copy_VarCoords(u(0,:,:), pu(0,0,:,:))
  copy_VarCoords(u(0,:,:), nu(0,0,:,:))

  pu := dim_avg_n_Wrap(pu, 0)
  nu := dim_avg_n_Wrap(nu, 0)

;;plot
  wks = gsn_open_wks("x11","eof")         
  plot = new(2,graphic)                
  plotc= plot         
  res              = True         
  res@gsnDraw      = False        
  res@gsnFrame     = False       
  res@gsnAddCyclic = False       
 
  res@mpFillOn     = False        
  res@mpMinLatF    = lats        
  res@mpMaxLatF    = latn
  res@mpMinLonF    = lonw
  res@mpMaxLonF    = lone
  res@cnFillOn      = True                
  res@cnFillPalette = "BlWhRe" 
  ;res@cnLevelSelectionMode   = "Manuallevels"
  ;res@cnMinLevelValF         = -1.
  ;res@cnMaxLevelValF         =  1.
  ;res@cnLevelSpacingF        =  0.1                                       
  symMinMaxPlt(pu, 11, False, res)   

  ;res = True
  ;res@gsnFrame     = False
  ;res@gsnDraw      = False
  ;res@gsnMaximize  = False
  ;res@gsnAddCyclic = False
  ;res@gsnContourNegLineDashPattern = 1
  ;res@gsnLeftString      = ""
  ;res@gsnRightString     = ""
  ;res@cnMonoLineColor    = True
  ;res@cnLinesOn          = True
  ;res@cnLineLabelsOn     = False
  ;res@cnLineColor        = "black"
  ;res@cnLineThicknessF   = 2.0
  ;res@cnInfoLabelPerimOn = False
  ;res@cnInfoLabelOn      = False
  ;res@cnLevelSelectionMode = "Explicitlevels"
  ;res@cnLevels   = (/20.,25.,30.,35./)

  resP                    = True         
  resP@gsnMaximize        = True         
  resP@gsnPanelLabelBar   = True        
;*******************************************
; first plot
;*******************************************
  ;nu = nu-pu
  ;copy_VarCoords(pu,nu)
  plot(0)=gsn_csm_contour_map(wks,dim_avg_n_Wrap(pu(6:8,:,:),0),res);
  plot(1)=gsn_csm_contour_map(wks,dim_avg_n_Wrap(nu(6:8,:,:),0),res) ;dim_avg_n_Wrap(nu(2:4,:,:),0)

  lonw = -70
  lone = -40
  lats =  35
  latn =  50
  polyx = (/lonw,lone,lone,lonw,lonw/)
  polyy = (/lats,lats,latn,latn,lats/)
  respoly = True
  respoly@gsLineDashPattern = 1
  dum3 = gsn_add_polyline(wks, plot(1), polyx, polyy, respoly)

  gsn_panel(wks,plot,(/2,1/),resP)  


  i = 0
  do while (i.eq.0)
    i = 0
  end do