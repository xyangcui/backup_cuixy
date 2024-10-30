
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

  IOP = new((/2,8/),float)

  IOP(0,:) = IO(2,:,0)
  IOP(1,:) = IO(3,:,0)

  fa = addfile(adir,"r")
  ia = fa->pc(:,0)*100

  res             = True  
  res@gsnDraw     = False
  res@gsnMaximize = False
  res@gsnFrame    = False
  res@gsnRightStringFontHeightF = 0.015
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnRightStringOrthogonalPosF = 0.0
  res@gsnLeftStringOrthogonalPosF  = 0.0

  res@tiXAxisString = "~F21~Frequency[%]"
  res@tiXAxisFontHeightF = 0.015
  res@tiYAxisFontHeightF = 0.015
  res@tiXAxisOffsetXF = 0.2

  res@tmYROn = False
  res@tmXTOn = False
  res@tmYLLabelFontHeightF = 0.012
  res@tmYLMode = "Explicit" 
  res@tmYLValues = (/1,2,3,4,5,6,7,8/)
  res@tmYLLabels = (/"1","2","3","4","5","6","7","8"/)
  res@tmXBLabelFontHeightF = 0.012


  res@vpXF = 0.6
  res@vpYF = 0.8   
  res@vpWidthF  = 0.25
  res@vpHeightF = 0.4

  res@trYMaxF = 8.7
  res@trYMinF = 0.3
  res@trXMinF = 5
  res@trXMaxF = 30
  res@gsnXRefLine  = 5

  res@gsnXYBarChart = True
  res@gsnXYBarChartBarWidth = 0.35

  wks  = gsn_open_wks ("x11", "bar" )

  odir = ddir+"/"+"boostrap_pc_omi_days.nc"
  fout = addfile(odir,"r")           
  pc   = fout->pc
  sdf  = dim_stddev_n(pc(0,:,:), 0)*100
  sds  = dim_stddev_n(pc(1,:,:), 0)*100

  pc_sorted = dim_pqsort_n(pc, 2, 1)
  ind_up = toint(1000*0.975-1)
  ind_dn = toint(1000*0.025-1)

  iosp_up = pc(1,ind_up,:)*100 - IOP(1,:)
  iosp_dn = pc(1,ind_dn,:)*100 - IOP(1,:)

  iofp_up = pc(0,ind_up,:)*100 - IOP(0,:)
  iofp_dn = pc(0,ind_dn,:)*100 - IOP(0,:)

  print(iosp_dn)
  print(iofp_up)

  odir = ddir+"/"+"boostrap_pc_omi_all_alldays.nc"
  fout = addfile(odir,"r")           
  pca   = fout->pc
  sda  = dim_stddev_n(pca, 0)*100


  xf = fspan(.825,7.825,8)
  xs = fspan(1.175,8.175,8)

  fcolor = "lightskyblue";chocolate1"
  scolor = "lightpink";"deepskyblue1"
  res@tiXAxisString = ""
  res@gsnLeftString = "(b) Diff events"
  res@gsnXYBarChartColors = (/fcolor/)  
  plot1 = gsn_csm_xy(wks,IOP(0,:),xf,res)   
  res@gsnXYBarChartColors = (/scolor/)         
  plot2 = gsn_csm_xy(wks,IOP(1,:),xs,res)  
  do i = 0,7
    ErrorBar(wks,plot1,xf(i),IOP(0,i),IOP(0,i)+iofp_up(i),IOP(0,i)+iofp_dn(i),"lightskyblue3")
    ErrorBar(wks,plot2,xs(i),IOP(1,i),IOP(1,i)+iosp_up(i),IOP(1,i)+iosp_dn(i),"lightpink3")
  end do
  overlay(plot1,plot2)

  res@tmYLOn = True
  res@tiXAxisString = "~F21~Frequency[%]"
  res@tiYAxisString = "~F21~MJO phase"
  res@gsnXYBarChartBarWidth = 0.6
  res@trXMinF = 5
  res@trXMaxF = 20
  res@gsnXRefLine  = 5
  res@gsnLeftString = "(a) All events"
  res@gsnXYBarChartColors = (/"slategray4"/) 
  res@vpWidthF  = 0.4
  res@vpHeightF = 0.4  
  res@vpXF = res@vpXF - res@vpWidthF 
  plot = gsn_csm_xy(wks,ia,fspan(1,8,8),res) 

  do i = 0,7
    ErrorBar(wks,plot,i+1,ia(i),ia(i)+sda(i),ia(i)-sda(i),"black")
  end do

  txres               = True 
  txres@txFontHeightF = 0.018
  txres@txFont = 22
;; EAJS
  txres@txFontColor = scolor
  text1 = gsn_add_text(wks, plot1, "IOsP", 24, 7.1, txres)
;; MEJS
  txres@txFontColor = fcolor
  text2 = gsn_add_text(wks, plot1, "IOfP", 24, 6.5, txres)

  draw(plot)
  draw(plot1)

  i = 0
  do while(i.eq.0)
    i = 0
  end do