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
function calc_weekly_apc(iosp_pi,iosp_pj,iofp_pi,iofp_pj,nB)
local iosp_p67,iofp_p67,boostrap_a,iosp_p67_week,iofp_p67_week,rseed1,rseed2,apc_sotred,nB
begin
  boostrap_a = new((/nB,4/),float)

;;calc week-mean.
  iosp_p67 = array_append_record(iosp_pi,iosp_pj,0)
  iofp_p67 = array_append_record(iofp_pi,iofp_pj,0)

  iosp_p67_week = new((/dimsizes(iosp_p67(:,0)),4/),float)
  iofp_p67_week = new((/dimsizes(iofp_p67(:,0)),4/),float)
  do i = 0, dimsizes(iosp_p67(:,0))-1
    iosp_p67_week(i,0) = dim_avg_n(iosp_p67(i,0+7*0:6+7*0),0)
    iosp_p67_week(i,1) = dim_avg_n(iosp_p67(i,0+7*1:6+7*1),0)
    iosp_p67_week(i,2) = dim_avg_n(iosp_p67(i,0+7*2:6+7*2),0)
    iosp_p67_week(i,3) = dim_avg_n(iosp_p67(i,0+7*3:6+7*3),0)
  end do

  do i = 0, dimsizes(iofp_p67(:,0))-1
    iofp_p67_week(i,0) = dim_avg_n(iofp_p67(i,0+7*0:6+7*0),0)
    iofp_p67_week(i,1) = dim_avg_n(iofp_p67(i,0+7*1:6+7*1),0)
    iofp_p67_week(i,2) = dim_avg_n(iofp_p67(i,0+7*2:6+7*2),0)
    iofp_p67_week(i,3) = dim_avg_n(iofp_p67(i,0+7*3:6+7*3),0)
  end do

  iop_p67 = array_append_record(iosp_p67_week,iofp_p67_week,0)

  rseed1 = toint(systemfunc(" date +%s"))
  rseed2 = toint((12345l*rseed1)%2147483398l)
  random_setallseed(rseed1, rseed2)              
  do i = 0, nB-1
    ;random seed.
    rseed1 = toint(systemfunc(" date +%s"))
    rseed2 = toint((12345l*rseed1)%2147483398l)+toint(random_normal(100, 10, 1))
    random_setallseed(rseed1, rseed2)

    iw = generate_sample_indices(dimsizes(iop_p67(:,0))-1,0)
    iwsf = iw(0:dimsizes(iosp_p67_week(:,0))-1)
    iwff = iw(dimsizes(iosp_p67_week(:,0)):)

    iosp_re = iop_p67(iwsf,:)
    iofp_re = iop_p67(iwff,:)
    iop = array_append_record(iosp_re,iofp_re,0)

    boostrap_a(i,:) = (dim_avg_n(iosp_re,0)-dim_avg_n(iofp_re,0))/(dim_avg_n(iop,0))*100

    delete([/iw,iwsf,iwff,iosp_re,iofp_re,iop/])
  end do

  apc_sotred = dim_pqsort_n(boostrap_a, 2, 1)

  return boostrap_a
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


  nBoot = 10000
  boostrap_acp = new((/4,nBoot,4/),float) 
;isac
  boostrap_acp(0,:,:) = calc_weekly_apc(isac_iosp_p6, isac_iosp_p7,\
                                  isac_iofp_p6, isac_iofp_p7,nBoot)
;ncep
  boostrap_acp(1,:,:) = calc_weekly_apc(ncep_iosp_p6, ncep_iosp_p7,\
                                  ncep_iofp_p6, ncep_iofp_p7,nBoot)
;cma
  boostrap_acp(2,:,:) =  calc_weekly_apc(cma_iosp_p6, cma_iosp_p7,\
                                  cma_iofp_p6, cma_iofp_p7,nBoot)
;ecmwf
  boostrap_acp(3,:,:) = calc_weekly_apc(ecmwf_iosp_p6, ecmwf_iosp_p7,\
                                   ecmwf_iofp_p6, ecmwf_iofp_p7,nBoot)

  boostrap_acp@dims = "(models,times,weeks)"
  boostrap_acp@model = "(isac,ncep,cma,ecmwf)"
  boostrap_acp@weeks = "(week1,week2,week3,week4)"

  system("/bin/rm -f "+ddir+"/boostrap_apc_sep_models.nc") 
  fs = addfile(ddir+"/boostrap_apc_sep_models.nc","c")  
  fs->apc = boostrap_acp  
  print("done")