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

  iosp_pi = dim_avg_n(iosp_pit(i,0+7*2:6+7*2),0)
  iofp_pi = dim_avg_n(iofp_pit(i,0+7*2:6+7*2),0)

  iosp_pj = dim_avg_n(iosp_pjt(i,0+7*2:6+7*2),0)
  iofp_pj = dim_avg_n(iofp_pjt(i,0+7*2:6+7*2),0)

  apc_diff = new((/3,4/),float)

  iosp_pij = array_append_record(iosp_pi,iosp_pj,0)
  iofp_pij = array_append_record(iofp_pi,iofp_pj,0)
;all phase.
  iop := array_append_record(iosp_pij,iofp_pij,0)
  apc_diff = (dim_avg(iosp_pij,0)-dim_avg(iofp_pij,0))/(dim_avg(iop))*100

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


  fb = addfile(ddir+"boostrap_apc_sep_models.nc","r")
  Boot = fb->apc

  std = dim_stddev_n(Boot, 1)

  print(std)