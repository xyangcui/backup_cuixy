;-----------------------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  idir = ddir + "/ERA5_daily_anom_ulevel_1979-2022_r2.5.nc"

  lats =  10
  latn =  80
  lonw = -60
  lone =  30

  neof   = 3  ; number of EOFs
  optEOF = True       
  optEOF@jopt = 0 

  optETS = False
;------------------------------------------------------------
; Main code.
;------------------------------------------------------------
  f = addfile(idir,"r")
  u = f->u(:,:,{lats:latn},{lonw:lone})
  time = u&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.4.and.day.eq.1)

  mon:= mon(n1:n2)
  u  := u(n1:n2,:,:,:)
  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  u := u(nt,:,:,:)
  delete([/n1,n2,nt/])

;; weight for greater latitude range.
  rad    = 4.*atan(1.)/180.
  clat   = u&lat        
  clat   = sqrt( cos(rad*clat) )    

  wu   = u ; copy meta data
  wu   = u*tofloat( conform(u, clat, 2) )
  wu@long_name = "Wgt: "+wu@long_name

; Reorder (lat,lon,time) the *weighted* input data
; Access the area of interest via coordinate subscripting
  xw = wu(level|:,lat|:,lon|:,time|:)
  x  = wu(time|:,level|:,lat|:,lon|:)

  eof    = eofunc_Wrap(xw, neof, optEOF)      
  eof_ts = eofunc_ts_Wrap(xw, eof, optETS)

  printVarSummary( eof )         
  printVarSummary( eof_ts )

  ;eval = sqrt(eof@eval)
  ;eof  = eof*conform(eof, eval, 0)
  eof_ts = dim_standardize_n_Wrap(eof_ts, 0, 1)
  ;eof_ts = eof_ts/conform(eof_ts, eval, 0) ;standardize.
;; write PC1 in CSV.
  ;filename = ddir+"/Ulevel_PC.csv"
  ;system("rm -rf " + filename)
  ;write_table(filename,"w",[/eof_ts(0,:),eof_ts(1,:)/],"%5.2f,%5.2f")
;; REOF
  ;eof_rot = eofunc_varimax_Wrap( eof, 1 )
  ;printVarSummary( eof_rot )        
  ;print("eof_rot: min="+min(eof_rot)+"   max="+max(eof_rot) )
  ;eofunc_varimax_reorder( eof_rot )
  ;printVarSummary( eof_rot )  
  ;eof_ts = eofunc_ts_Wrap(xw, eof_rot, optETS) 

;---Open a new NetCDF file to write to
  fout_name = ddir + "/EOF_ulevel.nc"
  system("rm -f " + fout_name)
  fout = addfile(fout_name,"c")
  fout->u = eof

; PLOTS
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
  res@lbLabelBarOn  = False  
  res@cnLinesOn              = False
  res@cnLineLabelsOn         = False
  res@cnLevelSelectionMode   = "Manuallevels"
  res@cnMinLevelValF         = -1.
  res@cnMaxLevelValF         =  1.
  res@cnLevelSpacingF        =  0.1
  res@cnLabelBarEndStyle     = "ExcludeOuterBoxes"      
                                     
  symMinMaxPlt(eof, 16, False, res)   

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
  resc@cnLineColor        = "black"
  resc@cnLineThicknessF   = 2.0
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels   = (/15.,20.,25.,30.,35./)

  resP                    = True         
  resP@gsnMaximize        = True         
  resP@gsnPanelLabelBar   = True        
;*******************************************
; first plot
;*******************************************
  do n=0,1
    res@gsnLeftString  = "EOF "+(n+1)
    res@gsnRightString = sprintf("%5.1f", eof@pcvar(n)) +"%"
    plot(n)=gsn_csm_contour_map(wks,eof(n,{250},:,:),res)
  end do
  gsn_panel(wks,plot,(/2,1/),resP)    

  i = 0
  do while (i.eq.0)
    i = 0
  end do