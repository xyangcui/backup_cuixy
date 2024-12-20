;--------------------------------------------
function re_apc_values(Idir,Ndir,Cdir,Edir)
local f,Idir,Ndir,Cdir,Edir,xi,xi_dims,xiout,xe,xn,xc
begin
  fi = addfile(Idir,"r")
  xi = fi->apc(:,:,0:29)
  xi_dims = dimsizes(xi)
  xiout = reshape(xi,(/xi_dims(0)*xi_dims(1),xi_dims(2)/))
  print("isac "+xi_dims(1))
  ;ncep
  fn = addfile(Ndir,"r")
  xn = fn->apc(:,:,0:29)
  xn_dims = dimsizes(xn)
  xnout = reshape(xn,(/xn_dims(0)*xn_dims(1),xn_dims(2)/))
  print("ncep "+xn_dims(1))
  ;cma
  fc = addfile(Cdir,"r")
  xc = fc->apc(:,:,0:29)
  xc_dims = dimsizes(xc)
  xcout = reshape(xc,(/xc_dims(0)*xc_dims(1),xc_dims(2)/))
  print("cma "+xc_dims(1))
  ;ecmwf
  fe = addfile(Edir,"r")
  xe = fe->apc(:,:,0:29)
  xe_dims = dimsizes(xe)
  xeout = reshape(xe,(/xe_dims(0)*xe_dims(1),xe_dims(2)/))
  print("ecmwf "+xe_dims(1))

  xout  = array_append_record(xiout,xnout,0)
  xout := array_append_record(xout,xcout,0)
  xout := array_append_record(xout,xeout,0)

  return xout
end

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
;-----------------------------------------------------------------------------------------

ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/hindcast/"

;isac.
model_name = "isac"
idir = new((/2,8/),string)
do i = 1,8
  idir(0,i-1) = ddir+model_name+"/apc_IOfP_"+model_name+"_p"+i+"_a1.nc"
  idir(1,i-1) = ddir+model_name+"/apc_IOsP_"+model_name+"_p"+i+"_a1.nc"
end do
;ncep.
model_name = "ncep"
ndir = new((/2,8/),string)
do i = 1,8
  ndir(0,i-1) = ddir+model_name+"/apc_IOfP_"+model_name+"_p"+i+"_a1.nc"
  ndir(1,i-1) = ddir+model_name+"/apc_IOsP_"+model_name+"_p"+i+"_a1.nc"
end do
;cma.
model_name = "cma"
cdir = new((/2,8/),string)
do i = 1,8
  cdir(0,i-1) = ddir+model_name+"/apc_IOfP_"+model_name+"_p"+i+"_a1.nc"
  cdir(1,i-1) = ddir+model_name+"/apc_IOsP_"+model_name+"_p"+i+"_a1.nc"
end do
;ecmwf.
model_name = "ecmwf"
edir = new((/2,8/),string)
do i = 1,8
  edir(0,i-1) = ddir+model_name+"/apc_IOfP_"+model_name+"_p"+i+"_a1.nc"
  edir(1,i-1) = ddir+model_name+"/apc_IOsP_"+model_name+"_p"+i+"_a1.nc"
end do

;------------------------------------------------------
;main code.
;------------------------------------------------------
  print("phase1")
  i = 1
  iosp_p1  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p1  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))
  print("")
  print("phase2")
  i = 2
  iosp_p2  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p2  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))
  print("")
  print("phase3")
  i = 3
  iosp_p3  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p3  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))
  print("")
  print("phase4")
  i = 4
  iosp_p4  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p4  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))
  print("")
  print("phase5")
  i = 5
  iosp_p5  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p5  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))
  print("")
  print("phase6")
  i = 6
  iosp_p6  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p6  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))
  print("")
  print("phase7")
  i = 7
  iosp_p7  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p7  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))
  print("")
  print("phase8")
  i = 8
  iosp_p8  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p8  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))

;;calc phase values.
  apc_diff = calc_weekly_apc(iosp_p6, iosp_p7, iofp_p6, iofp_p7)
  print(apc_diff)
;;boostrap
  fb = addfile(ddir+"boostrap_apc.nc","r")
  Boot = fb->apc

  std = dim_stddev_n(Boot, 0)
  
  nBoot = dimsizes(Boot(:,0))
  n025 = toint(0.025*nBoot)   ; indices for sorted array
  n975 = toint(0.975*nBoot)
  xBoot_025 = Boot(n025,:) 
  xBoot_975 = Boot(n975,:)     
  print(xBoot_975)
  print(xBoot_025)
  print(std)
