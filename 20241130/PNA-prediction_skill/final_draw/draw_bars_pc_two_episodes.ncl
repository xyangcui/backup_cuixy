
undef("ErrorBar")
procedure ErrorBar(wks1,plotebar,xx,yy,yup,ydown,corr)
local wks1,xx,yy,yup,ydown,up1,down1,polyres,polyres2,\
dum1,plotebar
begin
  mstring =       "-"
  fontnum =       19
  xoffset =       0.0
  yoffset =       0.0
  aspect  =       1E-10
  size    =       1.2
  angle   =       0.0

  new_index = NhlNewMarker(wks1, mstring, fontnum, xoffset, yoffset,\
   aspect, size, angle)

  COLOR    = corr
  LThick   = 3
  MarkSiz  = 0.02
  Dpattern = 0
  polyres                       = True
  polyres@gsLineColor           = COLOR
  polyres@gsLineDashPattern     = Dpattern
  polyres@gsLineThicknessF      = LThick
  polyres@tfPolyDrawOrder       = "PostDraw"

  res_lines = True
  res_lines@gsLineDashPattern = Dpattern                   
  res_lines@gsLineThicknessF  = LThick                    
  res_lines@gsLineColor       = COLOR     

  dum1   = gsn_add_polyline(wks1,plotebar,(/ydown,yup/),(/xx,xx/),polyres)
  up1    = gsn_add_polyline(wks1,plotebar,(/yup,yup/),(/xx-0.1,xx+0.1/),res_lines)
  down1  = gsn_add_polyline(wks1,plotebar,(/ydown,ydown/),(/xx-0.1,xx+0.1/),res_lines)

  str1 = unique_string("dum")
  str2 = unique_string("up")
  str3 = unique_string("down")
  plotebar@$str1$ = dum1
  plotebar@$str2$ = up1
  plotebar@$str3$ = down1

end

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  idir = ddir+"/"+"pattern_consistency_omi_alldays.nc"
  adir = ddir+"/"+"pattern_consistency_omi_all_alldays.nc"

  fi = addfile(idir,"r")
  IO = fi->pc
  IO = IO * 100

  IOP = new((/2,8/),float) ;0 fp 1 sp

  IOP(0,:) = IO(2,:,0)
  IOP(1,:) = IO(3,:,0)

  res             = True  
  res@gsnDraw     = False
  res@gsnMaximize = False
  res@gsnFrame    = False
  res@gsnRightStringFontHeightF = 0.015
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnRightStringOrthogonalPosF = 0.0
  res@gsnLeftStringOrthogonalPosF  = 0.0

  res@tiYAxisString = "~F21~Frequency[%]"
  res@tiYAxisFontHeightF = 0.015
  res@tiXAxisFontHeightF = 0.015
  ;res@tiYAxisOffsetXF = 0.2

  res@tmYROn = False
  res@tmXTOn = False
  res@tmXBLabelFontHeightF = 0.015
  res@tmXBMode = "Explicit" 
  res@tmXBValues = (/0.25,0.6/)
  res@tmXBLabels = (/"P4","P8"/)
  res@tmYLLabelFontHeightF = 0.015

  res@vpYF = 0.9
  res@vpXF = 0.55   
  res@vpWidthF  = 0.4
  res@vpHeightF = 0.25

  res@trXMaxF = 0.82
  res@trXMinF = 0.02
  res@trYMinF = 10
  res@trYMaxF = 25
  res@gsnYRefLine  = 5

  res@gsnXYBarChart = True
  res@gsnXYBarChartBarWidth = 0.1

  wks  = gsn_open_wks ("x11", "bar" )

  odir = ddir+"/"+"boostrap_pc_omi_days.nc"
  fout = addfile(odir,"r")           
  pc   = fout->pc
  printVarSummary(pc)
  nrand = dimsizes(pc(0,:,0))
  pc_me = pc
  do i = 0,nrand-1
    pc_me(0,i,:) = pc_me(0,i,:)*100 - IOP(0,:)
    pc_me(1,i,:) = pc_me(1,i,:)*100 - IOP(1,:)
  end do
  pc_sorted = dim_pqsort_n(pc_me, 2, 1)

  ind_up = toint(nrand*0.95-1)
  ind_dn = toint(nrand*0.05-1)

  iofp_up = pc_me(0,ind_up,:)
  iofp_dn = pc_me(0,ind_dn,:)

  iosp_up = pc_me(1,ind_up,:)
  iosp_dn = pc_me(1,ind_dn,:)

  xf = fspan(.2,.55,2)
  xs = fspan(.3,.65,2)

  fcolor = "dodgerblue1"
  scolor = "firebrick1"
  res@tiXAxisString = ""
  res@gsnLeftString = "(b) Diff events"
  res@gsnXYBarChartColors = (/fcolor/)
  res@xyLineColor = fcolor  
  plot1 = gsn_csm_xy(wks,xf,(/IOP(0,3),IOP(0,7)/),res)   
  res@gsnXYBarChartColors = (/scolor/)  
  res@xyLineColor = scolor       
  plot2 = gsn_csm_xy(wks,xs,(/IOP(1,3),IOP(1,7)/),res)  
  ;do i = 0,7
  ;  ErrorBar(wks,plot1,xf(i),IOP(0,i),IOP(0,i)-iofp_dn(i),IOP(0,i)-iofp_up(i),"dodgerblue4")
  ;  ErrorBar(wks,plot2,xs(i),IOP(1,i),IOP(1,i)-iosp_dn(i),IOP(1,i)-iosp_up(i),"firebrick4")
  ;end do
  overlay(plot1,plot2)

  txres               = True 
  txres@txFontHeightF = 0.018
  txres@txFont = 22
;; EAJS
  txres@txFontColor = scolor
  text1 = gsn_add_text(wks, plot1, "IOsP", 0.15, 22, txres)
;; MEJS
  txres@txFontColor = fcolor
  text2 = gsn_add_text(wks, plot1, "IOfP", 0.15, 20, txres)
  draw(plot1)

  undef("IO")
  undef("IOP")
  lats = 20
  latn = 70

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  zdir = ddir+"/"+"ERA5_z500_daily_anom_1979-2022_r2.5.nc"
  fz = addfile(zdir,"r")
  z  = fz->z(0,10:18,0)

  fa = addfile(adir,"r")
  IO = fa->pc
  IO = IO * 100

  IOP = new((/9,31/),float)
  IOP(0:7,:) = IO(:,:)
  IOP(8,:)  = IO(0,:)  
  copy_VarCoords(z, IOP(:,0))
  IOP!1 = "time"
  IOP := IOP(time|:,lat|:)
  printVarSummary(IOP)
  lat = IOP&lat

  cmap  = read_colormap_file("MPL_hot")
  cmap2 = read_colormap_file("ViBlGrWhYeOrRe")

  undef("res")
  res             = True  
  res@gsnDraw     = False
  res@gsnMaximize = False
  res@gsnFrame    = False
  res@gsnRightStringFontHeightF = 0.015
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnRightStringOrthogonalPosF = 0.0
  res@gsnLeftStringOrthogonalPosF  = 0.0

  res@tmYLLabelFontHeightF = 0.012
  res@tmXBLabelFontHeightF = 0.012

  res@tmYROn = False
  res@tmXTOn = False

  res@tiXAxisString = "~F21~RMM phase"
  res@tiYAxisString = "~F21~Lag days"
  res@tiYAxisFontHeightF = 0.015
  res@tiXAxisFontHeightF = 0.015

  res@tmXBMode = "Explicit" 
  res@tmXBValues = z&lat
  res@tmXBLabels = (/"1","2","3","4","5","6","7","8","1" /)

  res@cnFillOn  = True              
  res@cnLinesOn = False    

  res@lbOrientation  = "Vertical" 
  res@lbLabelStride  = 2          
  res@cnFillOn       = True              
  res@cnFillPalette  = cmap(::-1,:)  
  res@cnLevelSelectionMode = "ManualLevels"     
  res@cnMinLevelValF       =  10;10         
  res@cnMaxLevelValF       =  15;25              
  res@cnLevelSpacingF      =  0.5;1

  res@vpWidthF  = 0.3
  res@vpHeightF = 0.5
  res@vpXF = 0.1
  res@vpYF = 0.9
  res@gsnLeftString = "~F21~(a) Pattern consistency(all events)"
  plot = gsn_csm_time_lat(wks, IOP, res) 
;Marker
  pmres                        =  True
  pmres@gsMarkerColor          = "white"       
  pmres@gsMarkerSizeF          =  0.01      
  pmres@gsLineThicknessF       =  6.
  pmres@gsMarkerIndex = 16
           
  pm1 = gsn_add_polymarker(wks, plot, lat(3), 0, pmres)
  pm2 = gsn_add_polymarker(wks, plot, lat(7), 0, pmres)

  draw(plot)

  opt         = True
  opt@bin_min = 10
  opt@bin_max = 40

  pdff = pdfx(pc(0,:,7)*100,25,opt)
  pdfs = pdfx(pc(1,:,7)*100,25,opt)
  pdf = new((/2,dimsizes(pdff)/),double)
  pdf(0,:) = pdff
  pdf(1,:) = pdfs
  pdf@bin_center = pdfs@bin_center

  resc  = True
  resc@gsnDraw = False
  resc@gsnFrame = False
  resc@gsnMaximize = False
  resc@gsnRightStringFontHeightF = 0.015
  resc@gsnLeftStringFontHeightF  = 0.015
  resc@gsnRightStringOrthogonalPosF = 0.0
  resc@gsnLeftStringOrthogonalPosF  = 0.0

  resc@tmYLLabelFontHeightF = 0.012
  resc@tmXBLabelFontHeightF = 0.012

  resc@tmYROn = False
  resc@tmXTOn = False
  resc@xyLineColors = (/fcolor,scolor/)
  resc@xyLineThicknesses = (/4,4/)
  resc@xyDashPatterns = (/0,0/)
  resc@gsnYRefLine = 0.

  resc@vpXF = 0.75+0.03
  resc@vpYF = 0.58
  resc@vpWidthF  = 0.2
  resc@vpHeightF = 0.2
  resc@gsnLeftString = "~F21~(d) P8"
  resc@tiYAxisString = ""
  resc@tiXAxisString = "~F21~Frequency[%]"
  resc@tiXAxisFontHeightF = 0.012
  resc@tiYAxisFontHeightF = 0.012
  plotc = gsn_csm_xy (wks, pdf@bin_center, pdf, resc)
  draw(plotc)

;p4
  opt         = True
  opt@bin_min = 10
  opt@bin_max = 40

  pdff = pdfx(pc(0,:,3)*100,25,opt)
  pdfs = pdfx(pc(1,:,3)*100,25,opt)
  pdf = new((/2,dimsizes(pdff)/),double)
  pdf(0,:) = pdff
  pdf(1,:) = pdfs
  pdf@bin_center = pdfs@bin_center

  resc@tiYAxisString = "~F21~PDF"
  resc@tiXAxisString = "~F21~Frequency[%]"
  resc@tiXAxisFontHeightF = 0.012
  resc@tiYAxisFontHeightF = 0.012
  resc@tiYAxisOffsetXF = 0.068
  resc@vpXF = 0.5+0.03   
  resc@gsnLeftString = "~F21~(c) P4"
  plotc = gsn_csm_xy (wks, pdf@bin_center, pdf, resc)
  draw(plotc)

  i = 0
  do while(i.eq.0)
    i = 0
  end do