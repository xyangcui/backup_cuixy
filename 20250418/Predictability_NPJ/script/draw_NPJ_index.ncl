

latn = 80
lats = 5
lonw = 100
lone = 240

ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"
;----------------------------------------------------------------
;main code
;----------------------------------------------------------------
  f = addfile(ddir+"/ERA5_daily_anom_u200_1982-2022.nc","r")
  fu= addfile(ddir+"/ERA5_daily_mean_u200_1979-2024.nc","r")

  TIME=cd_calendar(f->time,0)
  timi=cd_calendar(f->time,-2)
  mon = TIME(:,1)
  nt = ind((timi.ge.19821201.and.timi.le.20220331)\
           .and.(mon.eq.12.or.mon.eq.1.or.mon.eq.2.or.mon.eq.3) )
  print(dimsizes(nt))

  x = f->u(nt,{lats:latn},{lonw:lone})
 
  umt = fu->u(nt,{lats:latn},{lonw:lone})
  um := dim_avg_n_Wrap(umt, 0)

  delete([/TIME,timi,mon,nt/])
  timi=cd_calendar(x&time,-2)
  printVarSummary(timi)
  ;std.
  xstd = dim_stddev_n_Wrap(x, 0)

  ;NPJ index.
  nlatn = 40
  nlats = 30
  nlonw = 175
  nlone = 205

  xt = x(:,{nlatn:nlats},{nlonw:nlone})
  pi = atan(1.)*4
  clat = cos(pi/180.*xt&lat)
  NPJ  = wgt_areaave(xt, clat, 1.0, 0)
  NPJ_index = dim_standardize(NPJ, 0) ;standardized index.
  pos = ind(NPJ_index.gt.1.5)
  neg = ind(NPJ_index.lt.-1.5)
  xpos= dim_avg_n_Wrap(x(pos,:,:),0)
  xneg= dim_avg_n_Wrap(x(neg,:,:),0)
  ump = dim_avg_n_Wrap(umt(pos,:,:),0)
  umn = dim_avg_n_Wrap(umt(neg,:,:),0)

  ;system("rm -rf " + ddir+"/NPJ_index.csv")
  ;dq = str_get_dq()
  ;fields = (/"date","NPJ_index"/)
  ;header = [/str_join(fields,",")/]
  ;alist  = [/timi,NPJ_index/]  ; List of variables to write.
  ;format = "%d,%2f"            ; Format to use for each variable.

  ;write_table(ddir+"/NPJ_index.csv", "w", header, "%s")   ; Write header to CSV file.
  ;write_table(ddir+"/NPJ_index.csv", "a", alist, format)  ; Write mixed data to CSV file.

  wks  = gsn_open_wks("x11","regress")         ; send graphics to PNG file
  cmap = read_colormap_file("WhiteBlueGreenYellowRed")
  res                       = True     
  res@gsnMaximize           = False             
  res@gsnDraw               = False
  res@gsnAddCyclic          = False
  res@gsnFrame              = False
  res@gsnMajorLatSpacing = 10

  res@cnFillOn              = True             
  res@cnLinesOn             = False            
  res@cnLineLabelsOn        = False            

  res@mpFillOn              = False 
  res@mpMinLatF = 10;lats
  res@mpMaxLatF = 60;latn
  res@mpMinLonF = 120;lonw
  res@mpMaxLonF = 240;lone

  res@cnFillOn               = True
  res@cnFillPalette          = cmap;(::-1,:)
  res@cnLinesOn              = False
  res@cnLineLabelsOn         = False
  res@cnLevelSelectionMode   = "Manuallevels"
  res@cnMinLevelValF         = 6
  res@cnMaxLevelValF         = 20
  res@cnLevelSpacingF        = 1
  res@lbLabelBarOn           = True

  resc = True
  resc@gsnFrame     = False
  resc@gsnDraw      = False
  resc@gsnMaximize  = False
  resc@gsnAddCyclic = False
  resc@gsnContourNegLineDashPattern = 1
  resc@gsnLeftString      = ""
  resc@gsnRightString     = ""
  resc@cnMonoLineColor    = True
  resc@cnLinesOn          = True
  resc@cnLineLabelsOn     = False
  resc@cnLineColor        = "blue";magenta
  resc@cnLineThicknessF   = 3.0
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels  = (/30.,40.,50.,60./)  

  ;STD
  res@gsnLeftString ="(a) Standard deviation"
  res@gsnRightString= ""
  res@vpXF = 0.13
  res@vpYF = 0.9
  plot1 = gsn_csm_contour_map_ce(wks,xstd,res)
  plotu = gsn_csm_contour(wks,um,resc)
  overlay(plot1, plotu)
  polyx = (/175,175,205,205,175/)
  polyy = (/30, 40, 40, 30, 30/)
  resl = True
  resl@gsLineColor = "black"
  resl@gsLineThicknessF = 4.
  dum1 = gsn_add_polyline(wks, plot1, polyx, polyy, resl)

  ;pos events
  cmap := read_colormap_file("cmocean_balance")
  rest                       = True     
  rest@gsnMaximize           = False             
  rest@gsnDraw               = False
  rest@gsnAddCyclic          = False
  rest@gsnFrame              = False
  rest@gsnMajorLatSpacing = 10

  rest@cnFillOn              = True             
  rest@cnLinesOn             = False            
  rest@cnLineLabelsOn        = False       

  rest@mpFillOn  = False 
  rest@mpMinLatF = 5
  rest@mpMaxLatF = 70
  rest@mpMinLonF = 120
  rest@mpMaxLonF = 240

  rest@cnFillOn               = True
  rest@cnFillPalette         := cmap
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -32
  rest@cnMaxLevelValF         =  32
  rest@cnLevelSpacingF        =  4
  rest@lbLabelBarOn           =  False

  rest@vpXF = 0.08
  rest@vpYF = 0.4
  rest@vpWidthF = 0.4
  rest@vpHeightF= 0.25

  rest@gsnLeftString = "(b) Positive "+"("+dimsizes(pos)+" days)"
  rest@gsnRightString= ""
  plot2  = gsn_csm_contour_map_ce(wks,xpos,rest)
  plotu2 = gsn_csm_contour(wks,ump,resc)
  ;overlay(plot2, plotu2)
  ;neg events
  rest@lbLabelBarOn =  True
  rest@pmLabelBarWidthF = 0.6
  rest@pmLabelBarOrthogonalPosF =  0.15
  rest@pmLabelBarParallelPosF   = -0.12

  rest@vpXF = 0.55
  rest@vpYF = 0.4
  rest@vpWidthF = 0.4
  rest@vpHeightF= 0.25

  rest@gsnLeftString = "(c) Negative "+"("+dimsizes(neg)+" days)"
  plot3  = gsn_csm_contour_map_ce(wks,xneg,rest)
  plotu3 = gsn_csm_contour(wks,umn,resc)
  ;overlay(plot3, plotu3)

  draw(plot1)
  draw(plot2)
  draw(plot3)

  i=0
  do while(i.eq.0)
    i=0
  end do