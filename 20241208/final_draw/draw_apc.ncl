;-------------------------------------------------------------------------
;function
;error bars.
;-------------------------------------------------------------------------
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

  dum1   = gsn_add_polyline(wks1,plotebar,(/xx,xx/),(/ydown,yup/),polyres)
  up1    = gsn_add_polyline(wks1,plotebar,(/xx-0.05,xx+0.05/),(/yup,yup/),res_lines)
  down1  = gsn_add_polyline(wks1,plotebar,(/xx-0.05,xx+0.05/),(/ydown,ydown/),res_lines)

  str1 = unique_string("dum")
  str2 = unique_string("up")
  str3 = unique_string("down")
  plotebar@$str1$ = dum1
  plotebar@$str2$ = up1
  plotebar@$str3$ = down1
end
;------------------------------------------------------------
;function
;取所有模式的APC值，连接到一起。
;------------------------------------------------------------
function re_apc_values(Sdir)
local fi,Sdir,xi,xi_dims,xout,fn,xn,xn_dims,xnout
begin
  fi = addfile(Sdir(0),"r")
  xi = fi->apc
  xi_dims = dimsizes(xi)
  xout = xi
  do i = 1,dimsizes(Sdir)-1
    fn = addfile(Sdir(i),"r")
    xn = fn->apc
    xn_dims = dimsizes(xn)
    xnout = xn
    xout := array_append_record(xout,xnout,0)
    delete([/fn,xn,xn_dims,xnout/])
  end do

  return xout
end
;------------------------------------------------------------------------
;function
;计算实际的APC。
;------------------------------------------------------------------------
;;calc phase values.
function calc_weekly_apc(iosp_pit)
local iosp_pit,iosp_pjt,iosp_pi,i,apc_diff,iop
begin
  iosp_pi = new((/dimsizes(iosp_pit(:,0)),4/),float)
  do i = 0, dimsizes(iosp_pit(:,0))-1
    iosp_pi(i,0) = dim_avg_n(iosp_pit(i,0+7*0:6+7*0),0)
    iosp_pi(i,1) = dim_avg_n(iosp_pit(i,0+7*1:6+7*1),0)
    iosp_pi(i,2) = dim_avg_n(iosp_pit(i,0+7*2:6+7*2),0)
    iosp_pi(i,3) = dim_avg_n(iosp_pit(i,0+7*3:6+7*3),0)
  end do

  iosp_pij = iosp_pi
  apc_diff = dim_avg_n(iosp_pij,0)

  return apc_diff

end

;--------------------------------------------------------------------------
ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/hindcast/"
;ecmwf.
model_name = "ecmwf"
sdir = new((/8,8/),string)
fdir = new((/8,8/),string)
adir = new((/8,8/),string)
do iyear = 2016,2023
  do i = 1,8
    sdir(iyear-2016,i-1) = ddir+model_name+"/"+iyear+"/apc/apc_IOsP"+"_p"+i+"_a1_ens.nc"
    fdir(iyear-2016,i-1) = ddir+model_name+"/"+iyear+"/apc/apc_IOfP"+"_p"+i+"_a1_ens.nc"
    adir(iyear-2016,i-1) = ddir+model_name+"/"+iyear+"/apc/apc_IO"+"_p"+i+"_a1_ens.nc"
  end do
end do
;----------------------------------------------------------------------------
;main code.
;----------------------------------------------------------------------------
  i = 1
  ecmwf_iosp_p1 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p1 = re_apc_values(fdir(:,i-1))
  ecmwf_io_p1 = re_apc_values(adir(:,i-1))

  i = 2
  ecmwf_iosp_p2 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p2 = re_apc_values(fdir(:,i-1))
  ecmwf_io_p2 = re_apc_values(adir(:,i-1))

  i = 3
  ecmwf_iosp_p3 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p3 = re_apc_values(fdir(:,i-1))
  ecmwf_io_p3 = re_apc_values(adir(:,i-1))

  i = 4
  ecmwf_iosp_p4 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p4 = re_apc_values(fdir(:,i-1))
  ecmwf_io_p4 = re_apc_values(adir(:,i-1))

  i = 5
  ecmwf_iosp_p5 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p5 = re_apc_values(fdir(:,i-1))
  ecmwf_io_p5 = re_apc_values(adir(:,i-1))

  i = 6
  ecmwf_iosp_p6 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p6 = re_apc_values(fdir(:,i-1))
  ecmwf_io_p6 = re_apc_values(adir(:,i-1))

  i = 7
  ecmwf_iosp_p7 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p7 = re_apc_values(fdir(:,i-1))
  ecmwf_io_p7 = re_apc_values(adir(:,i-1))

  i = 8
  ecmwf_iosp_p8 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p8 = re_apc_values(fdir(:,i-1))
  ecmwf_io_p8 = re_apc_values(adir(:,i-1))


;;iofp
  ecmwf_iofp_p  = array_append_record(ecmwf_iofp_p1,ecmwf_iofp_p2,0)
  ecmwf_iofp_p := array_append_record(ecmwf_iofp_p,ecmwf_iofp_p3,0)
  ecmwf_iofp_p := array_append_record(ecmwf_iofp_p,ecmwf_iofp_p4,0)
  ecmwf_iofp_p := array_append_record(ecmwf_iofp_p,ecmwf_iofp_p5,0)
  ecmwf_iofp_p := array_append_record(ecmwf_iofp_p,ecmwf_iofp_p6,0)
  ecmwf_iofp_p := array_append_record(ecmwf_iofp_p,ecmwf_iofp_p8,0)
;;iosp
  ecmwf_iosp_p  = array_append_record(ecmwf_iosp_p1,ecmwf_iosp_p2,0)
  ecmwf_iosp_p := array_append_record(ecmwf_iosp_p,ecmwf_iosp_p3,0)
  ecmwf_iosp_p := array_append_record(ecmwf_iosp_p,ecmwf_iosp_p4,0)
  ecmwf_iosp_p := array_append_record(ecmwf_iosp_p,ecmwf_iosp_p5,0)
  ecmwf_iosp_p := array_append_record(ecmwf_iosp_p,ecmwf_iosp_p6,0)
  ecmwf_iosp_p := array_append_record(ecmwf_iosp_p,ecmwf_iosp_p8,0)

;;io
  ecmwf_io_p  = array_append_record(ecmwf_io_p1,ecmwf_io_p2,0)
  ecmwf_io_p := array_append_record(ecmwf_io_p,ecmwf_io_p3,0)
  ecmwf_io_p := array_append_record(ecmwf_io_p,ecmwf_io_p4,0)
  ecmwf_io_p := array_append_record(ecmwf_io_p,ecmwf_io_p5,0)
  ecmwf_io_p := array_append_record(ecmwf_io_p,ecmwf_io_p6,0)
  ecmwf_io_p := array_append_record(ecmwf_io_p,ecmwf_io_p8,0)

;calc apc real values.

  ecmwf_apc = new((/3,8,4/),float)

;iosp
  i = 0
  ecmwf_apc(i,0,:) = calc_weekly_apc(ecmwf_iosp_p1)
  ecmwf_apc(i,1,:) = calc_weekly_apc(ecmwf_iosp_p2)
  ecmwf_apc(i,2,:) = calc_weekly_apc(ecmwf_iosp_p3)
  ecmwf_apc(i,3,:) = calc_weekly_apc(ecmwf_iosp_p4)
  ecmwf_apc(i,4,:) = calc_weekly_apc(ecmwf_iosp_p5)
  ecmwf_apc(i,5,:) = calc_weekly_apc(ecmwf_iosp_p6)
  ecmwf_apc(i,6,:) = calc_weekly_apc(ecmwf_iosp_p7)
  ecmwf_apc(i,7,:) = calc_weekly_apc(ecmwf_iosp_p8)

;iofp
  i = 2
  ecmwf_apc(i,0,:) = calc_weekly_apc(ecmwf_iofp_p1)
  ecmwf_apc(i,1,:) = calc_weekly_apc(ecmwf_iofp_p2)
  ecmwf_apc(i,2,:) = calc_weekly_apc(ecmwf_iofp_p3)
  ecmwf_apc(i,3,:) = calc_weekly_apc(ecmwf_iofp_p4)
  ecmwf_apc(i,4,:) = calc_weekly_apc(ecmwf_iofp_p5)
  ecmwf_apc(i,5,:) = calc_weekly_apc(ecmwf_iofp_p6)
  ecmwf_apc(i,6,:) = calc_weekly_apc(ecmwf_iofp_p7)
  ecmwf_apc(i,7,:) = calc_weekly_apc(ecmwf_iofp_p8)

;io
  i = 1
  ecmwf_apc(i,0,:) = calc_weekly_apc(ecmwf_io_p1)
  ecmwf_apc(i,1,:) = calc_weekly_apc(ecmwf_io_p2)
  ecmwf_apc(i,2,:) = calc_weekly_apc(ecmwf_io_p3)
  ecmwf_apc(i,3,:) = calc_weekly_apc(ecmwf_io_p4)
  ecmwf_apc(i,4,:) = calc_weekly_apc(ecmwf_io_p5)
  ecmwf_apc(i,5,:) = calc_weekly_apc(ecmwf_io_p6)
  ecmwf_apc(i,6,:) = calc_weekly_apc(ecmwf_io_p7)
  ecmwf_apc(i,7,:) = calc_weekly_apc(ecmwf_io_p8)

;;real apc.
  undef("res")
  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@gsnRightString = ""

  res@gsnRightStringFontHeightF = 0.015
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnRightStringOrthogonalPosF = 0.0
  res@gsnLeftStringOrthogonalPosF  = 0.0

  res@tmYROn = False
  res@tmXTOn = False
  res@tmYLLabelFontHeightF = 0.013
  res@tmXBMode = "Explicit" 
  res@tmXBValues = (/1,2,3,4,5,6,7,8/)
  res@tmXBLabels = (/"P1","P2","P3","P4","P5","P6","P7","P8"/)
  res@tmXBLabelFontHeightF = 0.015  

  res@tiYAxisString = "~F21~Z500 APC"
  res@tiYAxisFontHeightF = 0.015 

  res@gsnXYBarChart = True
  res@gsnXYBarChartBarWidth = 0.2
  res@gsnXYBarChartOutlineThicknessF = 0.01
         
  res@vpWidthF  = 0.4
  res@vpHeightF = 0.247 
  res@trXMinF = 0.3
  res@trXMaxF = 8.7


  wks = gsn_open_wks("x11", "bar")

  xf = ispan(1, 8, 1)
  xd = (/-0.2,0.,0.2/)
  colors = (/"red","blue","green"/)

  plot1 = new(3,graphic)
  plot2 = new(3,graphic)
  plot3 = new(3,graphic)

;week1
  res@tmYLOn = False
  res@tmXBOn = False
  res@gsnLeftString = ""
  res@tiYAxisString = ""

  res@trYMaxF = 0.9
  res@trYMinF = 0.4
  res@gsnYRefLine  = res@trYMinF
  res@tmYLFormat = "0@;*.1f"
  res@vpXF = 0.1
  res@vpYF = 0.97

  do i = 0,2
    if (i.eq.2) then
      res@tmYLOn = True
      res@tmXBOn = True
      res@gsnLeftString = "~F21~(a) week1"
      res@tiYAxisString = "~F21~Z500 APC"
    end if
    res@gsnXYBarChartColors = colors(i)
    res@xyLineColor = colors(i)
    plot1(i) = gsn_csm_xy(wks,xf+xd(i),ecmwf_apc(i,:,0),res)
  end do

;week2
  res@tmYLOn = False
  res@tmXBOn = False
  res@gsnLeftString = ""
  res@tiYAxisString = ""

  res@trYMaxF = 0.8
  res@trYMinF = 0.3
  res@tmYLFormat = "0@;*.1f"
  res@gsnYRefLine  = res@trYMinF
  res@vpXF = 0.1
  res@vpYF = 0.65
  do i = 0,2
    if (i.eq.2) then
      res@tmYLOn = True
      res@tmXBOn = True
      res@gsnLeftString = "~F21~(b) week2"
      res@tiYAxisString = "~F21~Z500 APC"
    end if
    res@gsnXYBarChartColors = colors(i)
    res@xyLineColor = colors(i)
    plot2(i) = gsn_csm_xy(wks,xf+xd(i),ecmwf_apc(i,:,1),res)
  end do

;week3
  res@tmYLOn = False
  res@tmXBOn = False
  res@gsnLeftString = ""
  res@tiYAxisString = ""

  res@trYMaxF = 0.6
  res@trYMinF = 0.1
  res@gsnYRefLine  = res@trYMinF
  res@vpXF = 0.1
  res@vpYF = 0.33

  do i = 0,2
    if (i.eq.2) then
      res@tmYLOn = True
      res@tmXBOn = True
      res@gsnLeftString = "~F21~(c) week3"
      res@tiYAxisString = "~F21~Z500 APC"
    end if
    res@gsnXYBarChartColors = colors(i)
    res@xyLineColor = colors(i)
    plot3(i) = gsn_csm_xy(wks,xf+xd(i),ecmwf_apc(i,:,2),res)
  end do

; Create some labelbars for a legend
  lbres                    = True          ; labelbar only resources
  lbres@vpWidthF           = 0.2           ; labelbar width
  lbres@vpHeightF          = 0.1           ; labelbar height
  lbres@lbBoxMajorExtentF  = 0.15          ; puts space between color boxes
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
  lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
  lbres@lbPerimOn          = False

  labels = (/"IOsP","All","IOfP"/)
  ;xpos   = (/0.22,0.37,0.52/)
  xpos   = (/0.05,0.2,0.35/)
  do i=0,2
    lbres@lbFillColors = colors(i)
    gsn_labelbar_ndc(wks,1,labels(i),xpos(i),0.07,lbres)
  end do 

;stats.
;week1
  fs = addfile(ddir+"boostrap_apc_w4.nc","r")
  sb = fs->apc ;(n,events,Phase)
  sbc = sb
  do i = 0,5000-1
    sbc(i,0,:) = sb(i,0,:) - ecmwf_apc(0,:,0)
    sbc(i,1,:) = sb(i,1,:) - ecmwf_apc(1,:,0)
    sbc(i,2,:) = sb(i,2,:) - ecmwf_apc(2,:,0)
  end do

  sb_sort = dim_pqsort_n(sbc,2,0)

  nup = toint(5000*0.95)-1
  ndn = toint(5000*0.05)-1

  sb_up = sbc(nup,:,:)
  sb_dn = sbc(ndn,:,:)
  print(sb_up)
  print(sb_dn)

;all events
  xf := (/1,2,3,4,5,6,7,8/)
  do i = 0,7
    ErrorBar(wks,plot1(1),xf(i),ecmwf_apc(1,i,0),\
             ecmwf_apc(1,i,0)-sb_dn(1,i),ecmwf_apc(1,i,0)+sb_dn(1,i),"black") ;ecmwf_apc_two_diff(i)-
  end do
;IOsP
  xf := (/1,2,3,4,5,6,7,8/) + xd(0)
  do i = 0,7
    ErrorBar(wks,plot1(0),xf(i),ecmwf_apc(0,i,0),\
             ecmwf_apc(0,i,0)-sb_dn(0,i),ecmwf_apc(0,i,0)+sb_dn(0,i),"black") ;ecmwf_apc_two_diff(i)-
  end do
;IOfP
  xf := (/1,2,3,4,5,6,7,8/) + xd(2)
  do i = 0,7
    ErrorBar(wks,plot1(2),xf(i),ecmwf_apc(2,i,0),\
             ecmwf_apc(2,i,0)-sb_dn(2,i),ecmwf_apc(2,i,0)+sb_dn(2,i),"black") ;ecmwf_apc_two_diff(i)-
  end do

  delete([/fs,sb,sbc,sb_sort,nup,ndn,sb_up,sb_dn,xf/])

;week2
  fs = addfile(ddir+"boostrap_apc_w2.nc","r")
  sb = fs->apc ;(n,events,Phase)
  sbc = sb
  do i = 0,5000-1
    sbc(i,0,:) = sb(i,0,:) - ecmwf_apc(0,:,1)
    sbc(i,1,:) = sb(i,1,:) - ecmwf_apc(1,:,1)
    sbc(i,2,:) = sb(i,2,:) - ecmwf_apc(2,:,1)
  end do

  sb_sort = dim_pqsort_n(sbc,2,0)

  nup = toint(5000*0.95)-1
  ndn = toint(5000*0.05)-1

  sb_up = sbc(nup,:,:)
  sb_dn = sbc(ndn,:,:)

;all events
  xf := (/1,2,3,4,5,6,7,8/)
  do i = 0,7
    ErrorBar(wks,plot2(1),xf(i),ecmwf_apc(1,i,1),\
             ecmwf_apc(1,i,1)-sb_dn(1,i),ecmwf_apc(1,i,1)+sb_dn(1,i),"black") ;ecmwf_apc_two_diff(i)-
  end do
;IOsP
  xf := (/1,2,3,4,5,6,7,8/) + xd(0)
  do i = 0,7
    ErrorBar(wks,plot2(0),xf(i),ecmwf_apc(0,i,1),\
             ecmwf_apc(0,i,1)-sb_dn(0,i),ecmwf_apc(0,i,1)+sb_dn(0,i),"black") ;ecmwf_apc_two_diff(i)-
  end do
;IOfP
  xf := (/1,2,3,4,5,6,7,8/) + xd(2)
  do i = 0,7
    ErrorBar(wks,plot2(2),xf(i),ecmwf_apc(2,i,1),\
             ecmwf_apc(2,i,1)-sb_dn(2,i),ecmwf_apc(2,i,1)+sb_dn(2,i),"black") ;ecmwf_apc_two_diff(i)-
  end do

  delete([/fs,sb,sbc,sb_sort,nup,ndn,sb_up,sb_dn,xf/])
;week3
  fs = addfile(ddir+"boostrap_apc_w3.nc","r")
  sb = fs->apc ;(n,events,Phase)
  sbc = sb
  do i = 0,5000-1
    sbc(i,0,:) = sb(i,0,:) - ecmwf_apc(0,:,2)
    sbc(i,1,:) = sb(i,1,:) - ecmwf_apc(1,:,2)
    sbc(i,2,:) = sb(i,2,:) - ecmwf_apc(2,:,2)
  end do

  sb_sort = dim_pqsort_n(sbc,2,0)

  nup = toint(5000*0.95)-1
  ndn = toint(5000*0.05)-1

  sb_up = sbc(nup,:,:)
  sb_dn = sbc(ndn,:,:)

;all events
  xf := (/1,2,3,4,5,6,7,8/)
  do i = 0,7
    ErrorBar(wks,plot3(1),xf(i),ecmwf_apc(1,i,2),\
             ecmwf_apc(1,i,2)-sb_dn(1,i),ecmwf_apc(1,i,2)+sb_dn(1,i),"black") ;ecmwf_apc_two_diff(i)-
  end do
;IOsP
  xf := (/1,2,3,4,5,6,7,8/) + xd(0)
  do i = 0,7
    ErrorBar(wks,plot3(0),xf(i),ecmwf_apc(0,i,2),\
             ecmwf_apc(0,i,2)-sb_dn(0,i),ecmwf_apc(0,i,2)+sb_dn(0,i),"black") ;ecmwf_apc_two_diff(i)-
  end do
;IOfP
  xf := (/1,2,3,4,5,6,7,8/) + xd(2)
  do i = 0,7
    ErrorBar(wks,plot3(2),xf(i),ecmwf_apc(2,i,2),\
             ecmwf_apc(2,i,2)-sb_dn(2,i),ecmwf_apc(2,i,2)+sb_dn(2,i),"black") ;ecmwf_apc_two_diff(i)-
  end do

  delete([/fs,sb,sbc,sb_sort,nup,ndn,sb_up,sb_dn,xf/])

  draw(plot1(0))
  draw(plot1(1))
  draw(plot1(2))

  draw(plot2(0))
  draw(plot2(1))
  draw(plot2(2))

  draw(plot3(0))
  draw(plot3(1))
  draw(plot3(2))

  i = 0
  do while(i.eq.0)
    i = 0
  end do
