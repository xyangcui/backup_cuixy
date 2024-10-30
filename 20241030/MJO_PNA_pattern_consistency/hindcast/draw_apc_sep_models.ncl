;--------------------------------------------
function re_apc_values(sdir)
local f,sdir,x,x_dims,xout
begin
  f = addfile(sdir,"r")
  x = f->apc(:,:,0:29)
  x_dims = dimsizes(x)
  xout = reshape(x,(/x_dims(0)*x_dims(1),x_dims(2)/))

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
;--------------------------------------------------------------------------------------------

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
  i = 1
  isac_iosp_p1  = re_apc_values(idir(1,i-1))
  ncep_iosp_p1  = re_apc_values(ndir(1,i-1))
  cma_iosp_p1   = re_apc_values(cdir(1,i-1))
  ecmwf_iosp_p1 = re_apc_values(edir(1,i-1))

  isac_iofp_p1  = re_apc_values(idir(0,i-1))
  ncep_iofp_p1  = re_apc_values(ndir(0,i-1))
  cma_iofp_p1   = re_apc_values(cdir(0,i-1))
  ecmwf_iofp_p1 = re_apc_values(edir(0,i-1))

  i = 2
  isac_iosp_p2  = re_apc_values(idir(1,i-1))
  ncep_iosp_p2  = re_apc_values(ndir(1,i-1))
  cma_iosp_p2   = re_apc_values(cdir(1,i-1))
  ecmwf_iosp_p2 = re_apc_values(edir(1,i-1))

  isac_iofp_p2  = re_apc_values(idir(0,i-1))
  ncep_iofp_p2  = re_apc_values(ndir(0,i-1))
  cma_iofp_p2   = re_apc_values(cdir(0,i-1))
  ecmwf_iofp_p2 = re_apc_values(edir(0,i-1))

  i = 3
  isac_iosp_p3  = re_apc_values(idir(1,i-1))
  ncep_iosp_p3  = re_apc_values(ndir(1,i-1))
  cma_iosp_p3   = re_apc_values(cdir(1,i-1))
  ecmwf_iosp_p3 = re_apc_values(edir(1,i-1))

  isac_iofp_p3  = re_apc_values(idir(0,i-1))
  ncep_iofp_p3  = re_apc_values(ndir(0,i-1))
  cma_iofp_p3   = re_apc_values(cdir(0,i-1))
  ecmwf_iofp_p3 = re_apc_values(edir(0,i-1))

  i = 4
  isac_iosp_p4  = re_apc_values(idir(1,i-1))
  ncep_iosp_p4  = re_apc_values(ndir(1,i-1))
  cma_iosp_p4   = re_apc_values(cdir(1,i-1))
  ecmwf_iosp_p4 = re_apc_values(edir(1,i-1))

  isac_iofp_p4  = re_apc_values(idir(0,i-1))
  ncep_iofp_p4  = re_apc_values(ndir(0,i-1))
  cma_iofp_p4   = re_apc_values(cdir(0,i-1))
  ecmwf_iofp_p4 = re_apc_values(edir(0,i-1))

  i = 5
  isac_iosp_p5  = re_apc_values(idir(1,i-1))
  ncep_iosp_p5  = re_apc_values(ndir(1,i-1))
  cma_iosp_p5   = re_apc_values(cdir(1,i-1))
  ecmwf_iosp_p5 = re_apc_values(edir(1,i-1))

  isac_iofp_p5  = re_apc_values(idir(0,i-1))
  ncep_iofp_p5  = re_apc_values(ndir(0,i-1))
  cma_iofp_p5   = re_apc_values(cdir(0,i-1))
  ecmwf_iofp_p5 = re_apc_values(edir(0,i-1))

  i = 6
  isac_iosp_p6  = re_apc_values(idir(1,i-1))
  ncep_iosp_p6  = re_apc_values(ndir(1,i-1))
  cma_iosp_p6   = re_apc_values(cdir(1,i-1))
  ecmwf_iosp_p6 = re_apc_values(edir(1,i-1))

  isac_iofp_p6  = re_apc_values(idir(0,i-1))
  ncep_iofp_p6  = re_apc_values(ndir(0,i-1))
  cma_iofp_p6   = re_apc_values(cdir(0,i-1))
  ecmwf_iofp_p6 = re_apc_values(edir(0,i-1))

  i = 7
  isac_iosp_p7  = re_apc_values(idir(1,i-1))
  ncep_iosp_p7  = re_apc_values(ndir(1,i-1))
  cma_iosp_p7   = re_apc_values(cdir(1,i-1))
  ecmwf_iosp_p7 = re_apc_values(edir(1,i-1))

  isac_iofp_p7  = re_apc_values(idir(0,i-1))
  ncep_iofp_p7  = re_apc_values(ndir(0,i-1))
  cma_iofp_p7   = re_apc_values(cdir(0,i-1))
  ecmwf_iofp_p7 = re_apc_values(edir(0,i-1))

  i = 8
  isac_iosp_p8  = re_apc_values(idir(1,i-1))
  ncep_iosp_p8  = re_apc_values(ndir(1,i-1))
  cma_iosp_p8   = re_apc_values(cdir(1,i-1))
  ecmwf_iosp_p8 = re_apc_values(edir(1,i-1))

  isac_iofp_p8  = re_apc_values(idir(0,i-1))
  ncep_iofp_p8  = re_apc_values(ndir(0,i-1))
  cma_iofp_p8   = re_apc_values(cdir(0,i-1))
  ecmwf_iofp_p8 = re_apc_values(edir(0,i-1))


;isac
  isac_apc_diff = calc_weekly_apc(isac_iosp_p6, isac_iosp_p7,\
                                  isac_iofp_p6, isac_iofp_p7)
;ncep
  ncep_apc_diff = calc_weekly_apc(ncep_iosp_p6, ncep_iosp_p7,\
                                  ncep_iofp_p6, ncep_iofp_p7)
;cma
  cma_apc_diff =  calc_weekly_apc(cma_iosp_p6, cma_iosp_p7,\
                                  cma_iofp_p6, cma_iofp_p7)
;ecmwf
  ecmwf_apc_diff = calc_weekly_apc(ecmwf_iosp_p2, ecmwf_iosp_p3,\
                                   ecmwf_iofp_p2, ecmwf_iofp_p3)

  print(isac_apc_diff)
  print(ncep_apc_diff)
  print(cma_apc_diff)
  print(ecmwf_apc_diff)

  fb = addfile(ddir+"boostrap_apc_sep_models.nc","r")
  Boot = fb->apc

  std = dim_stddev_n(Boot, 1)

  print(std)