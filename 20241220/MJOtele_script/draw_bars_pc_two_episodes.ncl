

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  idir = ddir+"/"+"pattern_consistency_omi_alldays.nc"


  fi = addfile(idir,"r")
  IO = fi->pc
  IO = IO * 100

  IOP = new((/2,8,5/),float)

  IOP(0,:,:) = IO(2,:,1:5)
  IOP(1,:,:) = IO(3,:,1:5)

  IOP := dim_avg_n(IOP,2)

  wks  = gsn_open_wks ("x11", "bar" )  

  res             = True  
  res@gsnDraw     = False
  res@gsnMaximize = False
  res@gsnFrame    = False
  res@gsnRightStringFontHeightF = 0.015
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnRightStringOrthogonalPosF = 0.0
  res@gsnLeftStringOrthogonalPosF  = 0.0

  res@tiXAxisString = "~F21~MJO phase"
  res@tiYAxisString = "~F21~[%]"
  res@tiYAxisFontHeightF = 0.015
  res@tiXAxisFontHeightF = 0.015

  res@tmYROn = False
  res@tmXTOn = False
  res@tmYLLabelFontHeightF = 0.012
  res@tmXBLabelFontHeightF = 0.012
  res@tmXBMode = "Explicit" 
  res@tmXBValues = (/1,2,3,4,5,6,7,8/)
  res@tmXBLabels = (/"1","2","3","4","5","6","7","8"/)

  ;res@vpXF = 0.18
  ;res@vpYF = 0.95   
  ;res@vpWidthF  = 0.25
  ;res@vpHeightF = 0.35

  res@trXMaxF = 9
  res@trXMinF = 0
  res@trYMinF = 10
  res@trYMaxF = 23
  res@gsnYRefLine  = 10.


  res@gsnXYBarChart = True
  res@gsnXYBarChartBarWidth = 0.25

  res@gsnLeftString = "(a)"
  res@gsnXYBarChartColors = (/"red3"/)  
  plot1 = gsn_csm_xy(wks,fspan(.875,7.875,8),IOP(0,:),res)   
  res@gsnXYBarChartColors = (/"blue3"/)         
  plot2 = gsn_csm_xy(wks,fspan(1.125,8.125,8),IOP(1,:),res)  

  overlay(plot1,plot2)
  draw(plot1)

  i = 0
  do while(i.eq.0)
    i = 0
  end do