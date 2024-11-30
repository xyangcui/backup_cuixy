  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  odir = ddir+"/"+"boostrap_pc_omi_days.nc"
  fout = addfile(odir,"r")           
  pc   = fout->pc
  printVarSummary(pc)

  fcolor = "dodgerblue1";chocolate1"
  scolor = "firebrick1";"deepskyblue1"

  resc  = True
  resc@gsnDraw = False
  resc@gsnFrame = False
  resc@gsnMaximize = False
  resc@gsnRightStringFontHeightF = 0.015
  resc@gsnLeftStringFontHeightF  = 0.015
  resc@gsnRightStringOrthogonalPosF = 0.0
  resc@gsnLeftStringOrthogonalPosF  = 0.0

  resc@tiYAxisString = "~F21~PDF"
  resc@tiXAxisString = "~F21~Frequency[%]"
  resc@tiXAxisFontHeightF = 0.015
  resc@tiYAxisFontHeightF = 0.015

  resc@tmYLLabelFontHeightF = 0.0135
  resc@tmXBLabelFontHeightF = 0.0135

  resc@tmYROn = False
  resc@tmXTOn = False
  resc@xyLineColors = (/fcolor,scolor/)
  resc@xyLineThicknesses = (/4,4/)
  resc@xyDashPatterns = (/0,0/)

  wks  = gsn_open_wks ("x11", "bar" )
  label = (/"~F21~(a) P1","~F21~(b) P2","~F21~(c) P3","",\
            "~F21~(d) P5","~F21~(e) P6","~F21~(f) P7",""/)
  plotc = new(6,graphic)
  resc@vpWidthF  = 0.35
  resc@vpHeightF = 0.216
;P1,2,3
  ;resc@vpYF = 0.45
  opt  = True
  opt@bin_min = 10
  opt@bin_max = 25
do i = 0,2
  pdff = pdfx(pc(0,:,i)*100,20,opt)
  pdfs = pdfx(pc(1,:,i)*100,20,opt)
  pdf = new((/2,dimsizes(pdff)/),double)
  pdf(0,:) = pdff
  pdf(1,:) = pdfs
  pdf@bin_center = pdfs@bin_center

  resc@gsnLeftString = label(i)
  plotc(i) = gsn_csm_xy(wks, pdf@bin_center, pdf, resc)
end do

do i = 3,5
  pdff = pdfx(pc(0,:,i+1)*100,20,opt)
  pdfs = pdfx(pc(1,:,i+1)*100,20,opt)
  pdf = new((/2,dimsizes(pdff)/),double)
  pdf(0,:) = pdff
  pdf(1,:) = pdfs
  pdf@bin_center = pdfs@bin_center

  resc@gsnLeftString = label(i+1)
  plotc(i) = gsn_csm_xy(wks, pdf@bin_center, pdf, resc)
end do

  resP = True
  resP@gsnPanelRowSpec = False

  gsn_panel(wks, plotc, (/2,3/), resP)

  i = 0
  do while(i.eq.0)
    i = 0
  end do