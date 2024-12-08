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
  ecmwf_iofp_p := array_append_record(ecmwf_iofp_p,ecmwf_iofp_p7,0)
;;iosp
  ecmwf_iosp_p  = array_append_record(ecmwf_iosp_p1,ecmwf_iosp_p2,0)
  ecmwf_iosp_p := array_append_record(ecmwf_iosp_p,ecmwf_iosp_p3,0)
  ecmwf_iosp_p := array_append_record(ecmwf_iosp_p,ecmwf_iosp_p4,0)
  ecmwf_iosp_p := array_append_record(ecmwf_iosp_p,ecmwf_iosp_p5,0)
  ecmwf_iosp_p := array_append_record(ecmwf_iosp_p,ecmwf_iosp_p6,0)
  ecmwf_iosp_p := array_append_record(ecmwf_iosp_p,ecmwf_iosp_p8,0)
  ecmwf_iosp_p := array_append_record(ecmwf_iosp_p,ecmwf_iosp_p7,0)

  iosp_p = new((/dimsizes(ecmwf_iosp_p(:,0)),4/),float)
  iofp_p = new((/dimsizes(ecmwf_iofp_p(:,0)),4/),float)

  iosp_p(:,0) = dim_avg_n(ecmwf_iosp_p(:,0+7*0:6+7*0),1)
  iosp_p(:,1) = dim_avg_n(ecmwf_iosp_p(:,0+7*1:6+7*1),1)
  iosp_p(:,2) = dim_avg_n(ecmwf_iosp_p(:,0+7*2:6+7*2),1)
  iosp_p(:,3) = dim_avg_n(ecmwf_iosp_p(:,0+7*3:6+7*3),1)

  iofp_p(:,0) = dim_avg_n(ecmwf_iofp_p(:,0+7*0:6+7*0),1)
  iofp_p(:,1) = dim_avg_n(ecmwf_iofp_p(:,0+7*1:6+7*1),1)
  iofp_p(:,2) = dim_avg_n(ecmwf_iofp_p(:,0+7*2:6+7*2),1)
  iofp_p(:,3) = dim_avg_n(ecmwf_iofp_p(:,0+7*3:6+7*3),1) 

  nBoot = 10000
  boostrap_acp = new((/nBoot,4/),float)    
  rseed1 = toint(systemfunc(" date +%s"))           
  do i = 0, nBoot-1
    ;random seed.
    do while ( rseed1.eq.toint(systemfunc(" date +%s")) ) 
      sleep(2)
    end do
    rseed1=toint(systemfunc(" date +%s"))
    rseed2=toint((12345l*rseed1)%2147483398l)+toint(random_normal(100, 10, 1))
    random_setallseed(rseed1, rseed2)

    iws = generate_sample_indices(dimsizes(iosp_p(:,0)),1)
    iwf = generate_sample_indices(dimsizes(iofp_p(:,0)),1)

    iosp_re = iosp_p(iws,:)
    iofp_re = iofp_p(iwf,:)
    iop = array_append_record(iosp_re,iofp_re,0)
    boostrap_acp(i,:) = (dim_avg_n(iosp_re,0)-dim_avg_n(iofp_re,0))/(dim_avg_n(iop,0))*100
    delete([/iws,iwf,iosp_re,iofp_re,iop/])
  end do

  boostrap_acp@dims = "(times,weeks)"
  boostrap_acp@weeks = "(week1,week2,week3,week4)"

  system("/bin/rm -f "+ddir+"/boostrap_apc.nc") 
  fs = addfile(ddir+"/boostrap_apc.nc","c")  
  fs->apc = boostrap_acp  
  print("done")