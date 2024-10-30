;--------------------------------------------
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
  iosp_p1  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p1  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))

  i = 2
  iosp_p2  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p2  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))

  i = 3
  iosp_p3  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p3  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))

  i = 4
  iosp_p4  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p4  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))

  i = 5
  iosp_p5  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p5  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))

  i = 6
  iosp_p6  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p6  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))

  i = 7
  iosp_p7  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p7  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))

  i = 8
  iosp_p8  = re_apc_values(idir(1,i-1),ndir(1,i-1),cdir(1,i-1),edir(1,i-1))
  iofp_p8  = re_apc_values(idir(0,i-1),ndir(0,i-1),cdir(0,i-1),edir(0,i-1))

;;calc phase values.
;P23
  p2 = new((/2,30/),float)
  p2(0,:) = dim_avg_n(iofp_p2,0)
  p2(1,:) = dim_avg_n(iosp_p2,0)

  p3 = new((/2,30/),float)
  p3(0,:) = dim_avg_n(iofp_p3,0)
  p3(1,:) = dim_avg_n(iosp_p3,0)

  p23 = new((/2,30/),float)
  iosp_p23 = array_append_record(iosp_p2,iosp_p3,0)
  iofp_p23 = array_append_record(iofp_p2,iofp_p3,0)
  p23(0,:) = dim_avg_n(iofp_p23,0)
  p23(1,:) = dim_avg_n(iosp_p23,0)

;P67
  p6 = new((/2,30/),float)
  p6(0,:) = dim_avg_n(iofp_p6,0)
  p6(1,:) = dim_avg_n(iosp_p6,0)

  p7 = new((/2,30/),float)
  p7(0,:) = dim_avg_n(iofp_p7,0)
  p7(1,:) = dim_avg_n(iosp_p7,0)

  p67 = new((/2,30/),float)
  iosp_p67 = array_append_record(iosp_p6,iosp_p7,0)
  iofp_p67 = array_append_record(iofp_p6,iofp_p7,0)
  p67(0,:) = dim_avg_n(iofp_p67,0)
  p67(1,:) = dim_avg_n(iosp_p67,0)

;P81
  p8 = new((/2,30/),float)
  p8(0,:) = dim_avg_n(iofp_p8,0)
  p8(1,:) = dim_avg_n(iosp_p8,0)

  p1 = new((/2,30/),float)
  p1(0,:) = dim_avg_n(iofp_p1,0)
  p1(1,:) = dim_avg_n(iosp_p1,0)

  p81 = new((/2,30/),float)
  iosp_p81 = array_append_record(iosp_p1,iosp_p8,0)
  iofp_p81 = array_append_record(iofp_p1,iofp_p8,0)
  p81(0,:) = dim_avg_n(iofp_p81,0)
  p81(1,:) = dim_avg_n(iosp_p81,0)
;;calc week-mean.
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

  nBoot = 10000
  boostrap_acp = new((/nBoot,4/),float) 
  rseed1 = toint(systemfunc(" date +%s"))
  rseed2 = toint((12345l*rseed1)%2147483398l)
  random_setallseed(rseed1, rseed2)              
  do i = 0, nBoot-1
    ;random seed.
    if(toint(systemfunc(" date +%s")).eq.rseed1) then
      rseed1 = toint(systemfunc(" date +%s"))
    else
      rseed1 = toint(systemfunc(" date +%s"))
    end if
    rseed2 = toint((12345l*rseed1)%2147483398l)+toint(random_normal(100, 10, 1))
    random_setallseed(rseed1, rseed2)

    ;iws = generate_sample_indices(dimsizes(iosp_p67(:,0))-1,1)
    ;iwf = generate_sample_indices(dimsizes(iofp_p67(:,0))-1,1)

    iw = generate_sample_indices(dimsizes(iop_p67(:,0))-1,0)
    iwsf = iw(0:dimsizes(iosp_p67_week(:,0))-1)
    iwff = iw(dimsizes(iosp_p67_week(:,0)):)

    ;iws25 = toint(dimsizes(iws)*0.25)
    ;iwf25 = toint(dimsizes(iwf)*0.25)
    ;iwsf = iws(1:iws25)
    ;iwff = iwf(1:iwf25)

    if i.lt.5
      print("time "+i)
      print(iwsf(0:5))
    end if
    ;iosp_re = iosp_p67_week(iwsf,:)
    ;iofp_re = iofp_p67_week(iwff,:)
    iosp_re = iop_p67(iwsf,:)
    iofp_re = iop_p67(iwff,:)
    iop = array_append_record(iosp_re,iofp_re,0)

    boostrap_acp(i,:) = (dim_avg_n(iosp_re,0)-dim_avg_n(iofp_re,0))/(dim_avg_n(iop,0))*100

    ;delete([/iws,iwf,iws25,iwf25,iwsf,iwff,iosp_re,iofp_re,iop/])
    delete([/iw,iwsf,iwff,iosp_re,iofp_re,iop/])
  end do

  apc_sotred = dim_pqsort_n(boostrap_acp, 2, 0)

  boostrap_acp@dims = "(times,weeks)"
  boostrap_acp@weeks = "(week1,week2,week3,week4)"

  system("/bin/rm -f "+ddir+"/boostrap_apc.nc") 
  fs = addfile(ddir+"/boostrap_apc.nc","c")  
  fs->apc = boostrap_acp  
  print("done")