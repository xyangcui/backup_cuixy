;------------------------------------------------
; 每种事件有n次预报，每次预报有m个member。
;循环就是n次预报第一层，m个member第二层。
;最终结果为一条线。
;------------------------------------------------
function calc_APC(mvalue,ovalue,cLat,ltime)
local mvalue,ovalue,cLat,cLaT,apc,apc1,apc2,apc3,i
begin
  ;mvalue (lead,lat,lon) ovalue(lead,lat,lon)
  cLaT = conform_dims(dimsizes(mvalue), cLat, 1)
  mvalue = mvalue * cLaT
  ovalue = ovalue * cLaT
  
  apc = new((/dimsizes(ltime)/),float)
  do i = 0, dimsizes(ltime)-1
    apc1 = mvalue(i,:,:)*ovalue(i,:,:)
    apc2 = mvalue(i,:,:)*mvalue(i,:,:)
    apc3 = ovalue(i,:,:)*ovalue(i,:,:)
    as1  = dim_sum(apc1)
    as1 := dim_sum(as1)
    as2  = dim_sum(apc2)
    as2 := dim_sum(as2)
    as3  = dim_sum(apc3)
    as3 := dim_sum(as3)

    apc(i) = as1/sqrt(as2*as3)
    delete([/as1,as2,as3/])
  end do
  return apc
end
;-----------------------------------------------------------------
;function calc predictability.
;------------------------------------------------------------------
function calc_predictbility(pf_z_s,pf_z_f)
local ens,pf_z_s,pf_z_f,clat,pre_iosp,i,apc_s,indice,ens_mean,apc_s,\
    pre_iofp,apc_f,predict
begin
  ;pf (number,initial_time,leadtime,lat,lon)
  ens = pf_z_s&number
  clat = tofloat(cos(pi*pf_z_s&lat/180)) 

  pred_iosp = new((/dimsizes(pf_z_s&initial_time),dimsizes(pf_z_s&lead_time)/),float)

  do i = 0,dimsizes(pf_z_s&initial_time)-1   
    apc_s = new((/dimsizes(pf_z_s&number),dimsizes(pf_z_s&lead_time)/),float) 
    do j = 0,dimsizes(pf_z_s&number)-1
      indice = ind(ens.ne.i)
      ens_mean = dim_avg_n_Wrap(pf_z_s(indice,i,:,:,:), 0)
      apc_s(j,:) = calc_APC(pf_z_s(j,i,:,:,:),ens_mean,clat,pf_z_s&lead_time)
      delete([/indice,ens_mean/])
    end do
    pre_iosp(i,:) = dim_avg_n(apcs, 1) 
    delete(apc_s)
  end do

  pre_iofp = new((/dimsizes(pf_z_f&initial_time),dimsizes(pf_z_f&lead_time)/),float)

  do i = 0,dimsizes(pf_z_f&initial_time)-1   
    apc_f = new((/dimsizes(pf_z_f&number),dimsizes(pf_z_f&lead_time)/),float) 
    do j = 0,dimsizes(pf_z_f&number)-1
      indice = ind(ens.ne.i)
      ens_mean = dim_avg_n_Wrap(pf_z_f(indice,i,:,:,:), 0)
      apc_f(j,:) = calc_APC(pf_z_f(j,i,:,:,:),ens_mean,clat,pf_z_f&lead_time)
      delete([/indice,ens_mean/])
    end do
    pre_iofp(i,:) = dim_avg_n(apc_f, 1) 
    delete(apc_f)
  end do

  predict = new((/2,dimsizes(pf_z_f&lead_time)/),float)
  predict(0,:) = dim_avg_n(pre_iosp,0)
  predict(1,:) = dim_avg_n(pre_iofp,0)

  return predict
end


lats = 20
latn = 70
lonw = 150
lone = 240

pi = atan(1)*4.
model_name = "ecmwf"

ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/hindcast/"

fsdir = ddir+model_name+"/"+Year+"/"+"IOsP_p"+i+".nc"
obsdir = ddir+model_name+"/"+Year+"/"+"IOsP_ob_p"+i+".nc"

ffdir = ddir+model_name+"/"+Year+"/"+"IOfP_p"+i+".nc"
obfdir = ddir+model_name+"/"+Year+"/"+"IOfP_ob_p"+i+".nc"

Sdir = new((/8,8/),float)
Fdir = new((/8,8/),float)

do iy = 2016,2023
  do i = 1,8
    Sdir(iy-2016,i-1) = ddir+model_name+"/"+Year+"/"+"IOsP_p"+i+".nc"
    Fdir(iy-2016,i-1) = ddir+model_name+"/"+Year+"/"+"IOfP_p"+i+".nc"
  end do
end do

;读入列表之中。 按位相和事件 则共16个列表。
  

  fs = addfile(fsdir,"r")
  fos = addfile(obsdir,"r")

  ff = addfile(ffdir,"r")
  fof = addfile(obfdir,"r")

  pf_z_s = fs->gh(:,:,:,{lats:latn},{lonw:lone})
  ob_z_s = fos->z(:,:,{lats:latn},{lonw:lone})

  pf_z_f = ff->gh(:,:,:,{lats:latn},{lonw:lone})
  ob_z_f = fof->z(:,:,{lats:latn},{lonw:lone})

  predic = new((/8,2,45/))

  do iy = 0,dimsizes(Sdir)-1
    fs = addfile(Sdir(iy),"r")
    f_z_s = fs->gh(:,:,:,{lats:latn},{lonw:lone})
    ff = addfile(Fdir(iy),"r")
    f_z_f = ff->gh(:,:,:,{lats:latn},{lonw:lone})

    predic(iy,:,:) = calc_predictbility(f_z_s,f_z_f)

    delete([/fs,f_z_s,ff,f_z_f/])
  end do

  predic!0 = "year"
  predic&year = ispan(2016,2023,1)
  predic!1 = "events"
  predic&events = ispan(0,1,1)
  predic!2 = "lead_time"
  predic&lead_time = ispan(1,45,1)


  system("/bin/rm -f "+ddir+model_name+"/predic_all.nc")  


  ff = addfile(ddir+model_name+"/predic_all.nc" ,"c")  
  ff->predic = predic
