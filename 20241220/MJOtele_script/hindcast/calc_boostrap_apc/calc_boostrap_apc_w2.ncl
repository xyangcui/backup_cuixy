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
function calc_weekly_apc(iosp_pit,nran)
local iosp_pit,iosp_pi,i
begin
  iosp_pi = new((/dimsizes(iosp_pit(:,0))/),float)
  do i = 0, dimsizes(iosp_pit(:,0))-1
    iosp_pi(i) = dim_avg_n(iosp_pit(i,0+7*1:6+7*1),0)
  end do

  boostrap_a = new(nran,float)
  rseed1 = toint(systemfunc(" date +%s"))             
  do i = 0, nran-1
    ;random seed.
    do while ( rseed1.eq.toint(systemfunc(" date +%s")) ) 
      sleep(2)
    end do
    rseed1=toint(systemfunc(" date +%s"))
    rseed2=toint((12345l*rseed1)%2147483398l)
    random_setallseed(rseed1, rseed2)

    iw = generate_sample_indices(dimsizes(iosp_pi),1)
    iosp_re = iosp_pi(iw)
    boostrap_a(i) = dim_avg(iosp_re)
    delete([/iw,iosp_re/])

  end do

  return boostrap_a

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


  print("week2 begin. "+ systemfunc("date"))
  nBoot = 5000
  boostrap_acp = new((/nBoot,3,8/),float) 

;iosp
;P1
  boostrap_acp(:,0,0) = calc_weekly_apc(ecmwf_iosp_p1,nBoot)
;P2
  boostrap_acp(:,0,1) = calc_weekly_apc(ecmwf_iosp_p2,nBoot)
;P3
  boostrap_acp(:,0,2) = calc_weekly_apc(ecmwf_iosp_p3,nBoot)
;P4
  boostrap_acp(:,0,3) = calc_weekly_apc(ecmwf_iosp_p4,nBoot)
;P5
  boostrap_acp(:,0,4) = calc_weekly_apc(ecmwf_iosp_p5,nBoot)
;P6
  boostrap_acp(:,0,5) = calc_weekly_apc(ecmwf_iosp_p6,nBoot)
;P7
  boostrap_acp(:,0,6) = calc_weekly_apc(ecmwf_iosp_p7,nBoot)
;P8
  boostrap_acp(:,0,7) = calc_weekly_apc(ecmwf_iosp_p8,nBoot)

;io
;P1
  boostrap_acp(:,1,0) = calc_weekly_apc(ecmwf_io_p1,nBoot)
;P2
  boostrap_acp(:,1,1) = calc_weekly_apc(ecmwf_io_p2,nBoot)
;P3
  boostrap_acp(:,1,2) = calc_weekly_apc(ecmwf_io_p3,nBoot)
;P4
  boostrap_acp(:,1,3) = calc_weekly_apc(ecmwf_io_p4,nBoot)
;P5
  boostrap_acp(:,1,4) = calc_weekly_apc(ecmwf_io_p5,nBoot)
;P6
  boostrap_acp(:,1,5) = calc_weekly_apc(ecmwf_io_p6,nBoot)
;P7
  boostrap_acp(:,1,6) = calc_weekly_apc(ecmwf_io_p7,nBoot)
;P8
  boostrap_acp(:,1,7) = calc_weekly_apc(ecmwf_io_p8,nBoot)

;iofp
;P1
  boostrap_acp(:,2,0) = calc_weekly_apc(ecmwf_iofp_p1,nBoot)
;P2
  boostrap_acp(:,2,1) = calc_weekly_apc(ecmwf_iofp_p2,nBoot)
;P3
  boostrap_acp(:,2,2) = calc_weekly_apc(ecmwf_iofp_p3,nBoot)
;P4
  boostrap_acp(:,2,3) = calc_weekly_apc(ecmwf_iofp_p4,nBoot)
;P5
  boostrap_acp(:,2,4) = calc_weekly_apc(ecmwf_iofp_p5,nBoot)
;P6
  boostrap_acp(:,2,5) = calc_weekly_apc(ecmwf_iofp_p6,nBoot)
;P7
  boostrap_acp(:,2,6) = calc_weekly_apc(ecmwf_iofp_p7,nBoot)
;P8
  boostrap_acp(:,2,7) = calc_weekly_apc(ecmwf_iofp_p8,nBoot)


  system("/bin/rm -f "+ddir+"/boostrap_apc_w2.nc") 
  fs = addfile(ddir+"/boostrap_apc_w2.nc","c")  
  fs->apc = boostrap_acp  

  print("week2 done. "+ systemfunc("date"))