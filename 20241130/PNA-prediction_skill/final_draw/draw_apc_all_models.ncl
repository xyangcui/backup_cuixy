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
  up1    = gsn_add_polyline(wks1,plotebar,(/xx-0.1,xx+0.1/),(/yup,yup/),res_lines)
  down1  = gsn_add_polyline(wks1,plotebar,(/xx-0.1,xx+0.1/),(/ydown,ydown/),res_lines)

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
function calc_weekly_apc(iosp_pit,iosp_pjt)
local iosp_pit,iosp_pjt,iosp_pi,i,apc_diff,iop
begin
  iosp_pi = new((/dimsizes(iosp_pit(:,0)),4/),float)
  do i = 0, dimsizes(iosp_pit(:,0))-1
    iosp_pi(i,0) = dim_avg_n(iosp_pit(i,0+7*0:6+7*0),0)
    iosp_pi(i,1) = dim_avg_n(iosp_pit(i,0+7*1:6+7*1),0)
    iosp_pi(i,2) = dim_avg_n(iosp_pit(i,0+7*2:6+7*2),0)
    iosp_pi(i,3) = dim_avg_n(iosp_pit(i,0+7*3:6+7*3),0)
  end do

  iosp_pj = new((/dimsizes(iosp_pjt(:,0)),4/),float)
  do i = 0, dimsizes(iosp_pjt(:,0))-1
    iosp_pj(i,0) = dim_avg_n(iosp_pjt(i,0+7*0:6+7*0),0)
    iosp_pj(i,1) = dim_avg_n(iosp_pjt(i,0+7*1:6+7*1),0)
    iosp_pj(i,2) = dim_avg_n(iosp_pjt(i,0+7*2:6+7*2),0)
    iosp_pj(i,3) = dim_avg_n(iosp_pjt(i,0+7*3:6+7*3),0)
  end do

  iosp_pij = array_append_record(iosp_pi,iosp_pj,0)
  apc_diff = dim_avg_n(iosp_pij,0)

  return apc_diff

end
;------------------------------------------------------------------------
;function
;计算每一周的APC的DIFF。
;------------------------------------------------------------------------
;;calc phase values.
function calc_weekly_apc_diff(iosp_pit,iosp_pjt,iofp_pit,iofp_pjt)
local iosp_pit,iosp_pjt,iofp_pit,iofp_pjt,iosp_pi,iofp_pi,i,apc_diff,iop
begin
  iosp_pi = new((/dimsizes(iosp_pit(:,0)),4/),float)
  iofp_pi = new((/dimsizes(iofp_pit(:,0)),4/),float)
  do i = 0, dimsizes(iosp_pit(:,0))-1
    iosp_pi(i,0) = dim_avg_n(iosp_pit(i,0+7*0:6+7*0),0)
    iosp_pi(i,1) = dim_avg_n(iosp_pit(i,0+7*1:6+7*1),0)
    iosp_pi(i,2) = dim_avg_n(iosp_pit(i,0+7*2:6+7*2),0)
    iosp_pi(i,3) = dim_avg_n(iosp_pit(i,0+7*3:6+7*3),0)
  end do

  do i = 0, dimsizes(iofp_pit(:,0))-1
    iofp_pi(i,0) = dim_avg_n(iofp_pit(i,0+7*0:6+7*0),0)
    iofp_pi(i,1) = dim_avg_n(iofp_pit(i,0+7*1:6+7*1),0)
    iofp_pi(i,2) = dim_avg_n(iofp_pit(i,0+7*2:6+7*2),0)
    iofp_pi(i,3) = dim_avg_n(iofp_pit(i,0+7*3:6+7*3),0)
  end do

  iosp_pj = new((/dimsizes(iosp_pjt(:,0)),4/),float)
  iofp_pj = new((/dimsizes(iofp_pjt(:,0)),4/),float)
  do i = 0, dimsizes(iosp_pjt(:,0))-1
    iosp_pj(i,0) = dim_avg_n(iosp_pjt(i,0+7*0:6+7*0),0)
    iosp_pj(i,1) = dim_avg_n(iosp_pjt(i,0+7*1:6+7*1),0)
    iosp_pj(i,2) = dim_avg_n(iosp_pjt(i,0+7*2:6+7*2),0)
    iosp_pj(i,3) = dim_avg_n(iosp_pjt(i,0+7*3:6+7*3),0)
  end do

  do i = 0, dimsizes(iofp_pjt(:,0))-1
    iofp_pj(i,0) = dim_avg_n(iofp_pjt(i,0+7*0:6+7*0),0)
    iofp_pj(i,1) = dim_avg_n(iofp_pjt(i,0+7*1:6+7*1),0)
    iofp_pj(i,2) = dim_avg_n(iofp_pjt(i,0+7*2:6+7*2),0)
    iofp_pj(i,3) = dim_avg_n(iofp_pjt(i,0+7*3:6+7*3),0)
  end do

  iosp_pij = array_append_record(iosp_pi,iosp_pj,0)
  iofp_pij = array_append_record(iofp_pi,iofp_pj,0)
;all phase.
  iop = array_append_record(iosp_pij,iofp_pij,0)
  apc_diff = (dim_avg_n(iosp_pij,0)-dim_avg_n(iofp_pij,0))/(dim_avg_n(iop,0))*100

  return apc_diff

end
;--------------------------------------------------------------------------
ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/hindcast/"
;ecmwf.
model_name = "ecmwf"
sdir = new((/8,8/),string)
fdir = new((/8,8/),string)
do iyear = 2016,2023
  do i = 1,8
    sdir(iyear-2016,i-1) = ddir+model_name+"/"+iyear+"/apc/apc_IOsP"+"_p"+i+"_a1_ens.nc"
    fdir(iyear-2016,i-1) = ddir+model_name+"/"+iyear+"/apc/apc_IOfP"+"_p"+i+"_a1_ens.nc"
  end do
end do
;----------------------------------------------------------------------------
;main code.
;----------------------------------------------------------------------------
  i = 1
  ecmwf_iosp_p1 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p1 = re_apc_values(fdir(:,i-1))

  i = 2
  ecmwf_iosp_p2 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p2 = re_apc_values(fdir(:,i-1))

  i = 3
  ecmwf_iosp_p3 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p3 = re_apc_values(fdir(:,i-1))

  i = 4
  ecmwf_iosp_p4 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p4 = re_apc_values(fdir(:,i-1))

  i = 5
  ecmwf_iosp_p5 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p5 = re_apc_values(fdir(:,i-1))

  i = 6
  ecmwf_iosp_p6 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p6 = re_apc_values(fdir(:,i-1))

  i = 7
  ecmwf_iosp_p7 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p7 = re_apc_values(fdir(:,i-1))

  i = 8
  ecmwf_iosp_p8 = re_apc_values(sdir(:,i-1))
  ecmwf_iofp_p8 = re_apc_values(fdir(:,i-1))


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

  print(dimsizes(ecmwf_iosp_p))
  print(dimsizes(ecmwf_iofp_p))
;ecmwf two events.
  ecmwf_apc_two_diff = new((/8/),float)
  ecmwf_apc_two_diff_w1 = new((/8/),float)
  ecmwf_apc_two_diff_w2 = new((/8/),float)
  ecmwf_apc_two_diff_w4 = new((/8/),float)
;P12
  ect = calc_weekly_apc_diff(ecmwf_iosp_p1, ecmwf_iosp_p2,ecmwf_iofp_p1, ecmwf_iofp_p2)
  ecmwf_apc_two_diff(0) = ect(2)
  ecmwf_apc_two_diff_w1(0) = ect(0)
  ecmwf_apc_two_diff_w2(0) = ect(1)
  ecmwf_apc_two_diff_w4(0) = ect(3)
;P23
  ect = calc_weekly_apc_diff(ecmwf_iosp_p2, ecmwf_iosp_p3,ecmwf_iofp_p2, ecmwf_iofp_p3)
  ecmwf_apc_two_diff(1) = ect(2)
  ecmwf_apc_two_diff_w1(1) = ect(0)
  ecmwf_apc_two_diff_w2(1) = ect(1)
  ecmwf_apc_two_diff_w4(1) = ect(3)
;P34
  ect = calc_weekly_apc_diff(ecmwf_iosp_p3, ecmwf_iosp_p4,ecmwf_iofp_p3, ecmwf_iofp_p4)
  ecmwf_apc_two_diff(2) = ect(2)
  ecmwf_apc_two_diff_w1(2) = ect(0)
  ecmwf_apc_two_diff_w2(2) = ect(1)
  ecmwf_apc_two_diff_w4(2) = ect(3)
;P45
  ect = calc_weekly_apc_diff(ecmwf_iosp_p4, ecmwf_iosp_p5,ecmwf_iofp_p4, ecmwf_iofp_p5)
  ecmwf_apc_two_diff(3) = ect(2)
  ecmwf_apc_two_diff_w1(3) = ect(0)
  ecmwf_apc_two_diff_w2(3) = ect(1)
  ecmwf_apc_two_diff_w4(3) = ect(3)
;P56
  ect = calc_weekly_apc_diff(ecmwf_iosp_p5, ecmwf_iosp_p6,ecmwf_iofp_p5, ecmwf_iofp_p6)
  ecmwf_apc_two_diff(4) = ect(2)
  ecmwf_apc_two_diff_w1(4) = ect(0)
  ecmwf_apc_two_diff_w2(4) = ect(1)
  ecmwf_apc_two_diff_w4(4) = ect(3)
;P67
  ect = calc_weekly_apc_diff(ecmwf_iosp_p6, ecmwf_iosp_p7,ecmwf_iofp_p6, ecmwf_iofp_p7)
  ecmwf_apc_two_diff(5) = ect(2)
  ecmwf_apc_two_diff_w1(5) = ect(0)
  ecmwf_apc_two_diff_w2(5) = ect(1)
  ecmwf_apc_two_diff_w4(5) = ect(3)
;P78
  ect = calc_weekly_apc_diff(ecmwf_iosp_p7, ecmwf_iosp_p8,ecmwf_iofp_p7, ecmwf_iofp_p8)
  ecmwf_apc_two_diff(6) = ect(2)
  ecmwf_apc_two_diff_w1(6) = ect(0)
  ecmwf_apc_two_diff_w2(6) = ect(1)
  ecmwf_apc_two_diff_w4(6) = ect(3)
;P81
  ect = calc_weekly_apc_diff(ecmwf_iosp_p8, ecmwf_iosp_p1,ecmwf_iofp_p8, ecmwf_iofp_p1)
  ecmwf_apc_two_diff(7) = ect(2)
  ecmwf_apc_two_diff_w1(7) = ect(0)
  ecmwf_apc_two_diff_w2(7) = ect(1)
  ecmwf_apc_two_diff_w4(7) = ect(3)


  wks = gsn_open_wks("x11", "bar")
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
  res@tmYLLabelFontHeightF = 0.01
  res@tmXBMode = "Explicit" 
  res@tmXBValues = (/1,2,3,4,5,6,7,8/)
  res@tmXBLabels = (/"P1","P2","P3","P4","P5","P6","P7","P8"/)
  res@tmXBLabelFontHeightF = 0.01  

  res@tiYAxisString = "~F21~Z500 APC %DIFF~C~~F21~    (IOsP~F18~-~F21~IOfP)"
  res@tiYAxisFontHeightF = 0.012

  res@gsnXYBarChart = True
  res@gsnXYBarChartBarWidth = 0.45
  res@gsnXYBarChartOutlineThicknessF = 0.01
  res@gsnAboveYRefLineColor = "red"           
  res@gsnBelowYRefLineColor = "blue"
  res@xyLineColors = (/"red","blue"/)          

  res@trYMaxF =  60
  res@trYMinF = -60
  res@trXMinF = 0.3
  res@trXMaxF = 8.7
  res@gsnYRefLine  = 0

  res@vpWidthF  = 0.4
  res@vpHeightF = 0.2

;week1
  res@vpXF = 0.1
  res@vpYF = 0.8
  res@gsnLeftString = "~F21~(a) week1"
  plot1 = gsn_csm_xy(wks,(/1,2,3,4,5,6,7,8/),ecmwf_apc_two_diff_w1,res)

;week2
  res@vpXF = 0.1
  res@vpYF = 0.5
  res@gsnLeftString = "~F21~(b) week2"
  plot2 = gsn_csm_xy(wks,(/1,2,3,4,5,6,7,8/),ecmwf_apc_two_diff_w2,res)

;week3
  res@trYMaxF =  200
  res@trYMinF = -200
  res@vpXF = 0.57
  res@vpYF = 0.8
  res@gsnLeftString = "~F21~(c) week3"
  res@tiYAxisString = ""
  res@tmYLOn = True
  plot3 = gsn_csm_xy(wks,(/1,2,3,4,5,6,7,8/),ecmwf_apc_two_diff,res)

;week4
  res@vpXF = 0.57
  res@vpYF = 0.5
  res@gsnLeftString = "~F21~(d) week4"
  plot4 = gsn_csm_xy(wks,(/1,2,3,4,5,6,7,8/),ecmwf_apc_two_diff_w4,res)


;stats.
  nrand = 5000
;week3
  fs = addfile(ddir+"boostrap_apc_sep_models_sp_w3.nc","r")
  sb = fs->apc
  do i = 0,nrand-1
    sb(:,i) = sb(:,i) - ecmwf_apc_two_diff
  end do
  sb_sort = dim_pqsort_n(sb,2,1)
  nup = toint(nrand*0.975-1)
  ndn = toint(nrand*0.025-1)
  sb_up = sb(:,nup)
  sb_dn = sb(:,ndn)

  xf := (/1,2,3,4,5,6,7,8/)
  do i = 0,7
    ErrorBar(wks,plot3,xf(i),ecmwf_apc_two_diff(i),\
             ecmwf_apc_two_diff(i)-sb_dn(i),ecmwf_apc_two_diff(i)-sb_up(i),"black") 
  end do
  ;do i = 0,7
  ;  ErrorBar(wks,plot3,xf(i),ecmwf_apc_two_diff(i),sb_up(i),sb_dn(i),"black") 
  ;end do
  delete([/fs,sb,sb_sort,nup,ndn,sb_up,sb_dn/])

;week1
  fs = addfile(ddir+"boostrap_apc_sep_models_sp_w1.nc","r")
  sb = fs->apc
  do i = 0,nrand-1
    sb(:,i) = sb(:,i) - ecmwf_apc_two_diff_w1
  end do
  sb_sort = dim_pqsort_n(sb,2,1)
  nup = toint(nrand*0.975-1)
  ndn = toint(nrand*0.025-1)
  sb_up = sb(:,nup)
  sb_dn = sb(:,ndn)

  do i = 0,7
    ErrorBar(wks,plot1,xf(i),ecmwf_apc_two_diff_w1(i),\
             ecmwf_apc_two_diff_w1(i)-sb_dn(i),ecmwf_apc_two_diff_w1(i)-sb_up(i),"black") 
  end do

  ;do i = 0,7
  ;  ErrorBar(wks,plot1,xf(i),ecmwf_apc_two_diff_w1(i),sb_up(i),sb_dn(i),"black") 
  ;end do

  delete([/fs,sb,sb_sort,nup,ndn,sb_up,sb_dn/])

;week2
  fs = addfile(ddir+"boostrap_apc_sep_models_sp_w2.nc","r")
  sb = fs->apc
  do i = 0,nrand-1
    sb(:,i) = sb(:,i) - ecmwf_apc_two_diff_w2
  end do
  sb_sort = dim_pqsort_n(sb,2,1)
  nup = toint(nrand*0.975-1)
  ndn = toint(nrand*0.025-1)
  sb_up = sb(:,nup)
  sb_dn = sb(:,ndn)

  do i = 0,7
    ErrorBar(wks,plot2,xf(i),ecmwf_apc_two_diff_w2(i),\
             ecmwf_apc_two_diff_w2(i)-sb_dn(i),ecmwf_apc_two_diff_w2(i)-sb_up(i),"black") 
  end do

  ;do i = 0,7
  ;  ErrorBar(wks,plot2,xf(i),ecmwf_apc_two_diff_w2(i),sb_up(i),sb_dn(i),"black") 
  ;end do
  delete([/fs,sb,sb_sort,nup,ndn,sb_up,sb_dn/])

;week4
  fs = addfile(ddir+"boostrap_apc_sep_models_sp_w4.nc","r")
  sb = fs->apc
  do i = 0,nrand-1
    sb(:,i) = sb(:,i) - ecmwf_apc_two_diff_w4
  end do
  sb_sort = dim_pqsort_n(sb,2,1)
  nup = toint(nrand*0.975-1)
  ndn = toint(nrand*0.025-1)
  sb_up = sb(:,nup)
  sb_dn = sb(:,ndn)

  do i = 0,7
    ErrorBar(wks,plot4,xf(i),ecmwf_apc_two_diff_w4(i),\
             ecmwf_apc_two_diff_w4(i)-sb_dn(i),ecmwf_apc_two_diff_w4(i)-sb_up(i),"black") 
  end do

  ;do i = 0,7
  ;  ErrorBar(wks,plot4,xf(i),ecmwf_apc_two_diff_w4(i),sb_up(i),sb_dn(i),"black") 
  ;end do
  delete([/fs,sb,sb_sort,nup,ndn,sb_up,sb_dn/])

  draw(plot1)
  draw(plot2)
  draw(plot3)
  draw(plot4)
          
  i = 0
  do while(i.eq.0)
    i = 0
  end do