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
;ecmwf all events.
  ecmwf_apc_diff = calc_weekly_apc_diff(ecmwf_iosp_p, ecmwf_iosp_p7,\
                                   ecmwf_iofp_p, ecmwf_iofp_p7)
  print(ecmwf_apc_diff)
;ecmwf two events.
  ecmwf_apc_two_diff = new((/8/),float)
;P12
  ect = calc_weekly_apc_diff(ecmwf_iosp_p1, ecmwf_iosp_p2,ecmwf_iofp_p1, ecmwf_iofp_p2)
  ecmwf_apc_two_diff(0) = ect(3)
;P23
  ect = calc_weekly_apc_diff(ecmwf_iosp_p2, ecmwf_iosp_p3,ecmwf_iofp_p2, ecmwf_iofp_p3)
  ecmwf_apc_two_diff(1) = ect(3)
;P34
  ect = calc_weekly_apc_diff(ecmwf_iosp_p3, ecmwf_iosp_p4,ecmwf_iofp_p3, ecmwf_iofp_p4)
  ecmwf_apc_two_diff(2) = ect(3)
;P45
  ect = calc_weekly_apc_diff(ecmwf_iosp_p4, ecmwf_iosp_p5,ecmwf_iofp_p4, ecmwf_iofp_p5)
  ecmwf_apc_two_diff(3) = ect(3)
;P56
  ect = calc_weekly_apc_diff(ecmwf_iosp_p5, ecmwf_iosp_p6,ecmwf_iofp_p5, ecmwf_iofp_p6)
  ecmwf_apc_two_diff(4) = ect(3)
;P67
  ect = calc_weekly_apc_diff(ecmwf_iosp_p6, ecmwf_iosp_p7,ecmwf_iofp_p6, ecmwf_iofp_p7)
  ecmwf_apc_two_diff(5) = ect(3)
;P78
  ect = calc_weekly_apc_diff(ecmwf_iosp_p7, ecmwf_iosp_p8,ecmwf_iofp_p7, ecmwf_iofp_p8)
  ecmwf_apc_two_diff(6) = ect(3)
;P81
  ect = calc_weekly_apc_diff(ecmwf_iosp_p8, ecmwf_iosp_p1,ecmwf_iofp_p8, ecmwf_iofp_p1)
  ecmwf_apc_two_diff(7) = ect(3)

  print(ecmwf_apc_two_diff)
  undef("apc_diff")
  undef("std")
  undef("Boot")
  fb = addfile(ddir+"boostrap_apc.nc","r")
  Boot = fb->apc

  nrand = dimsizes(Boot(:,0))
  ;rs_re = reshape(rs, (/2,nrand/10,10/))
  ;rs_me = dim_avg_n(rs_re,2)
  rs_me = Boot
  do i = 0,nrand - 1
    rs_me(i,0) = rs_me(i,0) - ecmwf_apc_diff(0)
    rs_me(i,1) = rs_me(i,1) - ecmwf_apc_diff(1)
    rs_me(i,2) = rs_me(i,2) - ecmwf_apc_diff(2)
    rs_me(i,3) = rs_me(i,3) - ecmwf_apc_diff(3)
  end do

  rs_sorted = dim_pqsort_n(rs_me, 2, 0)

  ind_up = toint(nrand*0.975-1)
  ind_dn = toint(nrand*0.025-1)

  sdf_up = rs_me(ind_up,:)
  sdf_dn = rs_me(ind_dn,:)
  res = True
  res@gsnDraw  = False
  res@gsnFrame = False

  res@gsnRightString = ""

  res@gsnRightStringFontHeightF = 0.018
  res@gsnLeftStringFontHeightF  = 0.018
  res@gsnRightStringOrthogonalPosF = 0.0
  res@gsnLeftStringOrthogonalPosF  = 0.0

  res@tmYROn = False
  res@tmXTOn = False
  res@tmYLLabelFontHeightF = 0.015
  res@tmXBMode = "Explicit" 
  res@tmXBValues = (/1,2,3,4/)
  res@tmXBLabels = (/"week1","week2","week3","week4"/)
  res@tmXBLabelFontHeightF = 0.015

  res@tiYAxisString = "~F21~Z500 APC %DIFF~C~    ~F21~(IOsP~F18~-~F21~IOfP)"
  res@tiYAxisFontHeightF = 0.015

  res@xyMonoLineColor   = True      
  res@xyLineColor       = "Red"
  res@xyLineThicknessF  = 4.  
  res@xyMonoDashPattern = 0

  wks = gsn_open_wks("x11", "bar")
  xf = ispan(1, 4, 1)
;ecmwf
  res@vpXF = 0.25
  res@vpYF = 0.85
  res@vpWidthF  = 0.5
  res@vpHeightF = 0.25
  res@trYMaxF =  50
  res@trYMinF = -50
  res@trXMinF = 1
  res@trXMaxF = 4
  res@gsnYRefLine  = 0
  res@gsnLeftString = "(a) All phases" 
  apc_diff = ecmwf_apc_diff
  plot4 = gsn_csm_xy(wks,xf,apc_diff,res)

  res@gsnXYFillColors = "LightPink"
  res@xyLineColor     = -1 
  mid_value = new((/2,4/),float)
  mid_value(0,:) = apc_diff-sdf_up
  mid_value(1,:) = apc_diff-sdf_dn
  mid_plot  = gsn_csm_xy(wks,xf,mid_value,res)
  overlay(plot4,mid_plot)
  draw(plot4)

;;real apc.
  undef("res")
  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@gsnRightString = ""

  res@gsnRightStringFontHeightF = 0.018
  res@gsnLeftStringFontHeightF  = 0.018
  res@gsnRightStringOrthogonalPosF = 0.0
  res@gsnLeftStringOrthogonalPosF  = 0.0

  res@tmYROn = False
  res@tmXTOn = False
  res@tmYLLabelFontHeightF = 0.015
  res@tmXBMode = "Explicit" 
  res@tmXBValues = (/1,2,3,4,5,6,7,8/)
  res@tmXBLabels = (/"P12","P23","P34","P45","P56","P67","P78","P81"/)
  res@tmXBLabelFontHeightF = 0.015  

  res@tiYAxisString = "~F21~Z500 APC %DIFF~C~~F21~    (IOsP~F18~-~F21~IOfP)"
  res@tiYAxisFontHeightF = 0.015 

  res@gsnXYBarChart = True
  res@gsnXYBarChartBarWidth = 0.45
  res@gsnXYBarChartOutlineThicknessF = 0.01
  res@gsnAboveYRefLineColor = "red"           
  res@gsnBelowYRefLineColor = "blue"
  res@xyLineColors = (/"red","blue"/)          

  res@vpWidthF  = 0.5
  res@vpHeightF = 0.309
  res@vpXF = 0.25
  res@vpYF = 0.5  
  res@trYMaxF =  120
  res@trYMinF = -90
  res@trXMinF = 0.3
  res@trXMaxF = 8.7
  res@gsnYRefLine  = 0
  res@gsnLeftString = "~F21~(b) week3"
  plot = gsn_csm_xy(wks,(/1,2,3,4,5,6,7,8/),ecmwf_apc_two_diff,res)


;stats.
  fs = addfile(ddir+"boostrap_apc_sep_models.nc","r")
  sb = fs->apc
  sbc = sb
  do i = 0,5000-1
    sbc(:,i) = sb(:,i) - ecmwf_apc_two_diff
  end do

  sb_sort = dim_pqsort_n(sbc,2,1)
  sb_s = dim_pqsort_n(sb,2,1)

  nup = toint(5000*0.975)-1
  ndn = toint(5000*0.025)-1
  printVarSummary(sbc)
  print(nup)
  sb_up = sbc(:,nup)
  sb_dn = sbc(:,ndn)

  print(sb_up)
  print(sb_dn)

  xf := (/1,2,3,4,5,6,7,8/)
  do i = 0,7
    ErrorBar(wks,plot,xf(i),ecmwf_apc_two_diff(i),\
             ecmwf_apc_two_diff(i)-sb_dn(i),ecmwf_apc_two_diff(i)-sb_up(i),"black") ;ecmwf_apc_two_diff(i)-
  end do

  draw(plot)
          
  ;plot = gsn_csm_xy(wks,xf,apc_real,res)      

; Create a plot with the area between both curves filled in blue.
  ;delete(res@xyLineColors)
  ;res@gsnXYFillColors = "LightBlue"
  ;res@xyLineColor     = -1                           
  ;bot_plot  = gsn_csm_xy (wks,time,mnmx(2:3,:),res)  

; Create a plot with the area between both curves filled in pink.
  ;res@gsnXYFillColors = "LightPink"
  ;res@xyLineColor     = -1                           
  ;mid_plot  = gsn_csm_xy (wks,time,mnmx(0:1,:),res)  

;  overlay(bot_plot,mid_plot)
;  overlay(bot_plot,top_plot)
;  draw(plot)
  i = 0
  do while(i.eq.0)
    i = 0
  end do
