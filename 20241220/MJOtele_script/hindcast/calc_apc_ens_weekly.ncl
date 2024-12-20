;------------------------------------------------
function calc_APC(mvalue,ovalue,cLat,itime,ltime)
local mvalue,ovalue,cLat,cLaT,apc,apc1,apc2,apc3,i,j
begin
  cLaT = conform_dims(dimsizes(mvalue), cLat, 2)
  mvalue = mvalue * cLaT
  ovalue = ovalue * cLaT

  apc = new((/dimsizes(itime),dimsizes(ltime)/),float)
  do i = 0,dimsizes(itime)-1
    do j = 0, dimsizes(ltime)-1
      apc1 = mvalue(i,j,:,:)*ovalue(i,j,:,:)
      apc2 = mvalue(i,j,:,:)*mvalue(i,j,:,:)
      apc3 = ovalue(i,j,:,:)*ovalue(i,j,:,:)
      as1  = dim_sum(apc1)
      as1 := dim_sum(as1)
      as2  = dim_sum(apc2)
      as2 := dim_sum(as2)
      as3  = dim_sum(apc3)
      as3 := dim_sum(as3)
      if as1.gt.0 then
        apc(i,j) = as1/sqrt(as2*as3)
      else
        apc(i,j) = 0.
      end if
      
      delete([/as1,as2,as3/])
    end do
  end do

  return apc

end

lats = 20
latn = 70
lonw = 150
lone = 240

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
  dimf = dimsizes(pf_z_f)

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

  clat = tofloat(cos(pi*pf_z_s&lat/180)) 

  apc_s = calc_APC(pf_z_s,ob_z_st,clat,pf_z_s&initial_time,pf_z_s&lead_time)
  apc_f = calc_APC(pf_z_f,ob_z_ft,clat,pf_z_f&initial_time,pf_z_f&lead_time)

  copy_VarCoords(pf_z_s(:,:,0,0), apc_s)
  copy_VarCoords(pf_z_f(:,:,0,0), apc_f)

  system("/bin/rm -f "+ddir+model_name+"/"+"apc_IOfP_"+"p"+i+"_a1_ens_weekly.nc")  
  system("/bin/rm -f "+ddir+model_name+"/"+"apc_IOsP_"+"p"+i+"_a1_ens_weekly.nc") 

  ff = addfile(ddir+model_name+"/"+"apc_IOfP"+"_p"+i+"_a1_ens_weekly.nc" ,"c")  
  ff->apc = apc_f

  fs = addfile(ddir+model_name+"/"+"apc_IOsP"+"_p"+i+"_a1_ens_weekly.nc" ,"c")  
  fs->apc = apc_s