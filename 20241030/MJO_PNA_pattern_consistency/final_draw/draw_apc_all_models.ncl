
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
;---------------------------------------------------------
;Function re_apc_values_single
;取单个模式的APC值。
;---------------------------------------------------------
function re_apc_values_single(sdir)
local f,sdir,x,x_dims,xout
begin
  f = addfile(sdir,"r")
  x = f->apc(:,:,0:29)
  x_dims = dimsizes(x)
  xout = reshape(x,(/x_dims(0)*x_dims(1),x_dims(2)/))

  return xout
end
;------------------------------------------------------------
;function
;取所有模式的APC值，连接到一起。
;------------------------------------------------------------
function re_apc_values(Idir,Ndir,Cdir,Edir)
local f,Idir,Ndir,Cdir,Edir,xi,xi_dims,xiout,xe,xn,xc
begin
  fi = addfile(Idir,"r")
  xi = fi->apc(:,:,0:29)
  xi_dims = dimsizes(xi)
  xiout = reshape(xi,(/xi_dims(0)*xi_dims(1),xi_dims(2)/))
  ;ncep
  fn = addfile(Ndir,"r")
  xn = fn->apc(:,:,0:29)
  xn_dims = dimsizes(xn)
  xnout = reshape(xn,(/xn_dims(0)*xn_dims(1),xn_dims(2)/))
  ;cma
  fc = addfile(Cdir,"r")
  xc = fc->apc(:,:,0:29)
  xc_dims = dimsizes(xc)
  xcout = reshape(xc,(/xc_dims(0)*xc_dims(1),xc_dims(2)/))
  ;ecmwf
  fe = addfile(Edir,"r")
  xe = fe->apc(:,:,0:29)
  xe_dims = dimsizes(xe)
  xeout = reshape(xe,(/xe_dims(0)*xe_dims(1),xe_dims(2)/))

  xout  = array_append_record(xiout,xnout,0)
  xout := array_append_record(xout,xcout,0)
  xout := array_append_record(xout,xeout,0)

  return xout
end
;------------------------------------------------------------------------
;function
;计算每一周的APC的DIFF。
;------------------------------------------------------------------------
;;calc phase values.
function calc_weekly_apc(iosp_pit,iosp_pjt,iofp_pit,iofp_pjt)
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

  apc_diff = new((/3,4/),float)

  iosp_pij = array_append_record(iosp_pi,iosp_pj,0)
  iofp_pij = array_append_record(iofp_pi,iofp_pj,0)
;single phase.
  iop  = array_append_record(iosp_pi,iofp_pi,0)
  apc_diff(0,:) = (dim_avg_n(iosp_pi,0)-dim_avg_n(iofp_pi,0))/(dim_avg_n(iop,0))*100

  iop := array_append_record(iosp_pi,iofp_pi,0)
  apc_diff(1,:) = (dim_avg_n(iosp_pj,0)-dim_avg_n(iofp_pj,0))/(dim_avg_n(iop,0))*100
;all phase.
  iop := array_append_record(iosp_pij,iofp_pij,0)
  apc_diff(2,:) = (dim_avg_n(iosp_pij,0)-dim_avg_n(iofp_pij,0))/(dim_avg_n(iop,0))*100

  return apc_diff

end
;--------------------------------------------------------------------------
ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/hindcast/"
;ecmwf.
model_name = "ecmwf"
edir = new((/2,8/),string)
do i = 1,8
  edir(0,i-1) = ddir+model_name+"/2020/apc/apc_IOfP_"+model_name+"_p"+i+".nc"
  edir(1,i-1) = ddir+model_name+"/2020/apc/apc_IOsP_"+model_name+"_p"+i+".nc"
end do
;----------------------------------------------------------------------------
;main code.
;----------------------------------------------------------------------------
;-------------------------------------------------------
;单个模式的APC。
;-------------------------------------------------------
  i = 1
  ecmwf_iosp_p1 = re_apc_values_single(edir(1,i-1))
  ecmwf_iofp_p1 = re_apc_values_single(edir(0,i-1))

  i = 2
  ecmwf_iosp_p2 = re_apc_values_single(edir(1,i-1))
  ecmwf_iofp_p2 = re_apc_values_single(edir(0,i-1))

  i = 3
  ecmwf_iosp_p3 = re_apc_values_single(edir(1,i-1))
  ecmwf_iofp_p3 = re_apc_values_single(edir(0,i-1))

  i = 4
  ecmwf_iosp_p4 = re_apc_values_single(edir(1,i-1))
  ecmwf_iofp_p4 = re_apc_values_single(edir(0,i-1))

  i = 5
  ecmwf_iosp_p5 = re_apc_values_single(edir(1,i-1))
  ecmwf_iofp_p5 = re_apc_values_single(edir(0,i-1))

  i = 6
  ecmwf_iosp_p6 = re_apc_values_single(edir(1,i-1))
  ecmwf_iofp_p6 = re_apc_values_single(edir(0,i-1))

  i = 7
  ecmwf_iosp_p7 = re_apc_values_single(edir(1,i-1))
  ecmwf_iofp_p7 = re_apc_values_single(edir(0,i-1))

  i = 8
  ecmwf_iosp_p8 = re_apc_values_single(edir(1,i-1))
  ecmwf_iofp_p8 = re_apc_values_single(edir(0,i-1))


  ecmwf_iofp_p6 := array_append_record(ecmwf_iofp_p6,ecmwf_iofp_p2,0)
  ecmwf_iofp_p6 := array_append_record(ecmwf_iofp_p6,ecmwf_iofp_p1,0)
  ecmwf_iofp_p6 := array_append_record(ecmwf_iofp_p6,ecmwf_iofp_p3,0)
  ecmwf_iofp_p6 := array_append_record(ecmwf_iofp_p6,ecmwf_iofp_p4,0)
  ecmwf_iofp_p6 := array_append_record(ecmwf_iofp_p6,ecmwf_iofp_p5,0)
  ecmwf_iofp_p6 := array_append_record(ecmwf_iofp_p6,ecmwf_iofp_p8,0)

;;iosp
  ecmwf_iosp_p6 := array_append_record(ecmwf_iosp_p6,ecmwf_iosp_p1,0)
  ecmwf_iosp_p6 := array_append_record(ecmwf_iosp_p6,ecmwf_iosp_p2,0)
  ecmwf_iosp_p6 := array_append_record(ecmwf_iosp_p6,ecmwf_iosp_p3,0)
  ecmwf_iosp_p6 := array_append_record(ecmwf_iosp_p6,ecmwf_iosp_p4,0)
  ecmwf_iosp_p6 := array_append_record(ecmwf_iosp_p6,ecmwf_iosp_p5,0)
  ecmwf_iosp_p6 := array_append_record(ecmwf_iosp_p6,ecmwf_iosp_p8,0)

;ecmwf
  ecmwf_apc_diff = calc_weekly_apc(ecmwf_iosp_p6, ecmwf_iosp_p7,\
                                   ecmwf_iofp_p6, ecmwf_iofp_p7)

  undef("apc_diff")
  undef("std")
  undef("Boot")
  fb = addfile(ddir+"boostrap_apc_sep_models.nc","r")
  Boot = fb->apc
  std = dim_stddev_n(Boot, 1)
  print(std)
  res@vpWidthF  = 0.35
  res@vpHeightF = 0.2
  res@gsnRightString = ""

;ecmwf
  res@vpXF = 0.12
  res@vpYF = 0.62  
  res@trYMaxF =  60
  res@trYMinF = -60
  res@gsnLeftString = "(b) ECMWF" 
  apc_diff = ecmwf_apc_diff(2,:) 
  Std = std(3,:)
  plot4 = gsn_csm_xy(wks,xf,apc_diff,res)   
  do i = 0,3
    ErrorBar(wks,plot4,xf(i),apc_diff(i),apc_diff(i)+Std(i),apc_diff(i)-Std(i),"black")
  end do
  draw(plot4)

  i = 0
  do while(i.eq.0)
    i = 0
  end do
