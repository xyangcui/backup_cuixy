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
;;calc phase values.
function calc_weekly_apc(iosp_pit,iosp_pjt,iofp_pit,iofp_pjt,nran)
local iosp_pit,iosp_pjt,iofp_pit,iofp_pjt,iosp_pi,iofp_pi,i,boostrap_a,\
  iop,iosp_pij,iofp_pij,rseed1,rseed2,iws,iwf,iosp_re,iofp_re
begin

  iosp_pi = dim_avg_n(iosp_pit(:,0+7*2:6+7*2),1)
  iofp_pi = dim_avg_n(iofp_pit(:,0+7*2:6+7*2),1)

  iosp_pj = dim_avg_n(iosp_pjt(:,0+7*2:6+7*2),1)
  iofp_pj = dim_avg_n(iofp_pjt(:,0+7*2:6+7*2),1)

  iosp_pij = array_append_record(iosp_pi,iosp_pj,0)
  iofp_pij = array_append_record(iofp_pi,iofp_pj,0)

  boostrap_a = new(nran,float)
  rseed1 = toint(systemfunc(" date +%s"))             
  do i = 0, nran-1
    ;random seed.
    do while ( rseed1.eq.toint(systemfunc(" date +%s")) ) 
      sleep(2)
    end do
    rseed1=toint(systemfunc(" date +%s"))
    rseed2=toint((12345l*rseed1)%2147483398l)+toint(random_normal(100, 10, 1))
    random_setallseed(rseed1, rseed2)

    iws = generate_sample_indices(dimsizes(iosp_pij),1)
    iwf = generate_sample_indices(dimsizes(iofp_pij),1)
    iosp_re = iosp_pij(iws)
    iofp_re = iofp_pij(iwf)
    iop = array_append_record(iosp_re,iofp_re,0)
    boostrap_a(i) = (dim_avg(iosp_re)-dim_avg(iofp_re))/(dim_avg(iop))*100
    delete([/iws,iwf,iosp_re,iofp_re,iop/])
  end do

  return boostrap_a

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
  iosp_p1 = re_apc_values(sdir(:,i-1))
  iofp_p1 = re_apc_values(fdir(:,i-1))

  i = 2
  iosp_p2 = re_apc_values(sdir(:,i-1))
  iofp_p2 = re_apc_values(fdir(:,i-1))

  i = 3
  iosp_p3 = re_apc_values(sdir(:,i-1))
  iofp_p3 = re_apc_values(fdir(:,i-1))

  i = 4
  iosp_p4 = re_apc_values(sdir(:,i-1))
  iofp_p4 = re_apc_values(fdir(:,i-1))

  i = 5
  iosp_p5 = re_apc_values(sdir(:,i-1))
  iofp_p5 = re_apc_values(fdir(:,i-1))

  i = 6
  iosp_p6 = re_apc_values(sdir(:,i-1))
  iofp_p6 = re_apc_values(fdir(:,i-1))

  i = 7
  iosp_p7 = re_apc_values(sdir(:,i-1))
  iofp_p7 = re_apc_values(fdir(:,i-1))

  i = 8
  iosp_p8 = re_apc_values(sdir(:,i-1))
  iofp_p8 = re_apc_values(fdir(:,i-1))
;------------------------------------------------------------------------
  nBoot = 5000
  boostrap_acp = new((/8,nBoot/),float) 
;P12
  boostrap_acp(0,:) = calc_weekly_apc(iosp_p1,iosp_p2,iofp_p1,iofp_p2,nBoot)
;P23
  boostrap_acp(1,:) = calc_weekly_apc(iosp_p2,iosp_p3,iofp_p2,iofp_p3,nBoot)
;P34
  boostrap_acp(2,:) = calc_weekly_apc(iosp_p3,iosp_p4,iofp_p3,iofp_p4,nBoot)
;P45
  boostrap_acp(3,:) = calc_weekly_apc(iosp_p4,iosp_p5,iofp_p4,iofp_p5,nBoot)
;P56
  boostrap_acp(4,:) = calc_weekly_apc(iosp_p5,iosp_p6,iofp_p5,iofp_p6,nBoot)
;P67
  boostrap_acp(5,:) = calc_weekly_apc(iosp_p6,iosp_p7,iofp_p6,iofp_p7,nBoot)
;P78
  boostrap_acp(6,:) = calc_weekly_apc(iosp_p7,iosp_p8,iofp_p7,iofp_p8,nBoot)
;P81
  boostrap_acp(7,:) = calc_weekly_apc(iosp_p8,iosp_p1,iofp_p8,iofp_p1,nBoot)

  system("/bin/rm -f "+ddir+"/boostrap_apc_sep_models.nc") 
  fs = addfile(ddir+"/boostrap_apc_sep_models.nc","c")  
  fs->apc = boostrap_acp  
  print("done"+ systemfunc("date"))