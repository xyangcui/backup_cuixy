;------------------------------------------------
;Function ACC  single point.
;------------------------------------------------
function calc_ACC(mvalue,ovalue,itime,ltime)
local mvalue,ovalue,i,acc
begin
 
  acc = new((/dimsizes(ltime)/),float)
  do i = 0,dimsizes(ltime)-1
    mv = mvalue(:,i)
    ov = ovalue(:,i)
    acc1 = dim_sum(mv*ov)
    acc2 = dim_sum(mv*mv)
    acc3 = dim_sum(ov*ov)
    if acc1.gt.0 then
      acc(i) = acc1/sqrt(acc2*acc3)
    else
      acc(i) = 0.
    end if

    delete([/mv,ov,acc1,acc2,acc3/])
  end do

  return acc
end

lats = 10
latn = 80
lonw = 140
lone = 270

pi = atan(1)*4.
model_name = "ecmwf"

ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/hindcast/"
sdir = new((/8/),string)
fdir = new((/8/),string)
osdir = new((/8/),string)
ofdir = new((/8/),string)
do iyear = 2016,2023
  sdir(iyear-2016) = ddir+model_name+"/"+iyear+"/"+"IOsP_p"+i+"_a1.nc"
  fdir(iyear-2016) = ddir+model_name+"/"+iyear+"/"+"IOfP_p"+i+"_a1.nc"
  osdir(iyear-2016) = ddir+model_name+"/"+iyear+"/"+"IOsP_ob_p"+i+"_a1.nc"
  ofdir(iyear-2016) = ddir+model_name+"/"+iyear+"/"+"IOfP_ob_p"+i+"_a1.nc"
end do

  fs  = addfile(sdir(0),"r")
  fos = addfile(osdir(0),"r")

  ff  = addfile(fdir(0),"r")
  fof = addfile(ofdir(0),"r")

  pf_z_s = fs->gh(:,:,:,{lats:latn},{lonw:lone})
  ob_z_s = fos->z(:,:,{lats:latn},{lonw:lone})

  pf_z_s := pf_z_s(initial_time|:,number|:,lead_time|:,lat|:,lon|:)

  pf_z_f = ff->gh(:,:,:,{lats:latn},{lonw:lone})
  ob_z_f = fof->z(:,:,{lats:latn},{lonw:lone})

  pf_z_f := pf_z_f(initial_time|:,number|:,lead_time|:,lat|:,lon|:)

  do j = 1 ,7

    fs  = addfile(sdir(j),"r")
    fos = addfile(osdir(j),"r")

    ff  = addfile(fdir(j),"r")
    fof = addfile(ofdir(j),"r")

    pf_s = fs->gh(:,:,:,{lats:latn},{lonw:lone})
    ob_s = fos->z(:,:,{lats:latn},{lonw:lone})

    pf_f = ff->gh(:,:,:,{lats:latn},{lonw:lone})
    ob_f = fof->z(:,:,{lats:latn},{lonw:lone})

    pf_s := pf_s(initial_time|:,number|:,lead_time|:,lat|:,lon|:)
    pf_f := pf_f(initial_time|:,number|:,lead_time|:,lat|:,lon|:)

    pf_z_s := array_append_record(pf_z_s,pf_s,0)
    pf_z_f := array_append_record(pf_z_f,pf_f,0)

    ob_z_s := array_append_record(ob_z_s,ob_s,0)
    ob_z_f := array_append_record(ob_z_f,ob_f,0)

    delete([/fs,fos,ff,fof,pf_s,ob_s,pf_f,ob_f,pf_s/])

  end do

  dims = dimsizes(pf_z_s)
  print("IOsP "+" P"+i+": "+dims(0))  
  dimf = dimsizes(pf_z_f)
  print("IOfP "+" P"+i+": "+dimf(0)) 

;ens mean
  pf_z_st := dim_avg_n_Wrap(pf_z_s, 1)
  pf_z_ft := dim_avg_n_Wrap(pf_z_f, 1)

;weeks mean.
  undef("pf_z_s")
  undef("pf_z_f")
  pf_z_s = new((/dims(0),4,dims(3),dims(4)/),float)
  pf_z_f = new((/dimf(0),4,dimf(3),dimf(4)/),float)

  ob_z_st = new((/dims(0),4,dims(3),dims(4)/),float)
  ob_z_ft = new((/dimf(0),4,dimf(3),dimf(4)/),float)

  do j = 0, dims(0)-1

    pf_z_s(j,0,:,:) = dim_avg_n(pf_z_st(j,0+7*0:6+7*0,:,:),0)
    pf_z_s(j,1,:,:) = dim_avg_n(pf_z_st(j,0+7*1:6+7*1,:,:),0)
    pf_z_s(j,2,:,:) = dim_avg_n(pf_z_st(j,0+7*2:6+7*2,:,:),0)
    pf_z_s(j,3,:,:) = dim_avg_n(pf_z_st(j,0+7*3:6+7*3,:,:),0)

    ob_z_st(j,0,:,:) = dim_avg_n(ob_z_s(j,0+7*0:6+7*0,:,:),0)
    ob_z_st(j,1,:,:) = dim_avg_n(ob_z_s(j,0+7*1:6+7*1,:,:),0)
    ob_z_st(j,2,:,:) = dim_avg_n(ob_z_s(j,0+7*2:6+7*2,:,:),0)
    ob_z_st(j,3,:,:) = dim_avg_n(ob_z_s(j,0+7*3:6+7*3,:,:),0)

  end do  

  do j = 0, dimf(0)-1

    pf_z_f(j,0,:,:) = dim_avg_n(pf_z_ft(j,0+7*0:6+7*0,:,:),0)
    pf_z_f(j,1,:,:) = dim_avg_n(pf_z_ft(j,0+7*1:6+7*1,:,:),0)
    pf_z_f(j,2,:,:) = dim_avg_n(pf_z_ft(j,0+7*2:6+7*2,:,:),0)
    pf_z_f(j,3,:,:) = dim_avg_n(pf_z_ft(j,0+7*3:6+7*3,:,:),0)

    ob_z_ft(j,0,:,:) = dim_avg_n(ob_z_f(j,0+7*0:6+7*0,:,:),0)
    ob_z_ft(j,1,:,:) = dim_avg_n(ob_z_f(j,0+7*1:6+7*1,:,:),0)
    ob_z_ft(j,2,:,:) = dim_avg_n(ob_z_f(j,0+7*2:6+7*2,:,:),0)
    ob_z_ft(j,3,:,:) = dim_avg_n(ob_z_f(j,0+7*3:6+7*3,:,:),0)

  end do

  copy_VarCoords(pf_z_st(:,0,:,:),pf_z_s(:,0,:,:))
  copy_VarCoords(pf_z_ft(:,0,:,:),pf_z_f(:,0,:,:))
  pf_z_s!1 = "lead_time"
  pf_z_f!1 = "lead_time"
  pf_z_s&lead_time = ispan(1, 4, 1)
  pf_z_f&lead_time = ispan(1, 4, 1)

;calc ACC single point
  nlat = dimsizes(pf_z_s&lat)
  nlon = dimsizes(pf_z_s&lon)

  acc_s = new((/dimsizes(pf_z_s&lead_time),nlat,nlon/),float)
  acc_f = new((/dimsizes(pf_z_f&lead_time),nlat,nlon/),float)

  do k = 0, nlat-1
    do j = 0, nlon-1
      acc_s(:,k,j) = calc_ACC(pf_z_s(:,:,k,j),ob_z_st(:,:,k,j),pf_z_s&initial_time,pf_z_s&lead_time)
      acc_f(:,k,j) = calc_ACC(pf_z_f(:,:,k,j),ob_z_ft(:,:,k,j),pf_z_f&initial_time,pf_z_f&lead_time)
    end do
  end do

  copy_VarCoords(pf_z_s(0,:,:,:), acc_s)
  copy_VarCoords(pf_z_f(0,:,:,:), acc_f)

  acc_d = acc_s - acc_f
  copy_VarCoords(acc_s, acc_d)

  system("/bin/rm -f "+ddir+model_name+"/"+"acc_IOfP_"+"p"+i+"_a1.nc")  
  system("/bin/rm -f "+ddir+model_name+"/"+"acc_IOsP_"+"p"+i+"_a1.nc") 
  system("/bin/rm -f "+ddir+model_name+"/"+"acc_diff_"+"p"+i+"_a1.nc") 

  ff = addfile(ddir+model_name+"/"+"acc_IOfP"+"_p"+i+"_a1.nc" ,"c")  
  ff->acc = acc_f

  fd = addfile(ddir+model_name+"/"+"acc_diff"+"_p"+i+"_a1.nc" ,"c")  
  fd->acc = acc_d

  fs = addfile(ddir+model_name+"/"+"acc_IOsP"+"_p"+i+"_a1.nc" ,"c")  
  fs->acc = acc_s