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

      apc(i,j) = as1/sqrt(as2*as3)
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

ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/hindcast/"

cfsdir = ddir+model_name+"/"+"IOsP_cf_p"+i+"_a0.8.nc"
pfsdir = ddir+model_name+"/"+"IOsP_pf_p"+i+"_a0.8.nc"
obsdir = ddir+model_name+"/"+"IOsP_ob_p"+i+"_a0.8.nc"

cffdir = ddir+model_name+"/"+"IOfP_cf_p"+i+"_a0.8.nc"
pffdir = ddir+model_name+"/"+"IOfP_pf_p"+i+"_a0.8.nc"
obfdir = ddir+model_name+"/"+"IOfP_ob_p"+i+"_a0.8.nc"

  
  fcs = addfile(cfsdir,"r")
  fps = addfile(pfsdir,"r")
  fos = addfile(obsdir,"r")

  fcf = addfile(cffdir,"r")
  fpf = addfile(pffdir,"r")
  fof = addfile(obfdir,"r")

  cf_z_s = fcs->gh(:,:,{lats:latn},{lonw:lone})
  pf_z_s = fps->gh(:,:,:,{lats:latn},{lonw:lone})
  ob_z_s = fos->z(:,:,{lats:latn},{lonw:lone})

  cf_z_f = fcf->gh(:,:,{lats:latn},{lonw:lone})
  pf_z_f = fpf->gh(:,:,:,{lats:latn},{lonw:lone})
  ob_z_f = fof->z(:,:,{lats:latn},{lonw:lone})

  printVarSummary(cf_z_s)

  apc_s = new((/dimsizes(pf_z_s&number)+1,dimsizes(cf_z_s&initial_time),dimsizes(cf_z_s&lead_time)/),float)
  clat = tofloat(cos(pi*ob_z_s&lat/180)) 
  apc_s(0,:,:) = calc_APC(cf_z_s,ob_z_s,clat,cf_z_s&initial_time,cf_z_s&lead_time)
  do j = 1,dimsizes(pf_z_s&number)
    apc_s(j,:,:) = calc_APC(pf_z_s(j-1,:,:,:,:),ob_z_s,clat,cf_z_s&initial_time,cf_z_s&lead_time)
  end do

  apc_f = new((/dimsizes(pf_z_f&number)+1,dimsizes(cf_z_f&initial_time),dimsizes(cf_z_f&lead_time)/),float)
  apc_f(0,:,:) = calc_APC(cf_z_f,ob_z_f,clat,cf_z_f&initial_time,cf_z_f&lead_time)
  do j = 1,dimsizes(pf_z_f&number)
    apc_f(j,:,:) = calc_APC(pf_z_f(j-1,:,:,:,:),ob_z_f,clat,cf_z_f&initial_time,cf_z_f&lead_time)
  end do

  copy_VarCoords(cf_z_s(:,:,0,0), apc_s(0,:,:))
  apc_s!0 = "number"
  apc_s&number = ispan(0,dimsizes(pf_z_s&number),1)

  copy_VarCoords(cf_z_f(:,:,0,0), apc_f(0,:,:))
  apc_f!0 = "number"
  apc_f&number = ispan(0,dimsizes(pf_z_f&number),1)

  system("/bin/rm -f "+ddir+model_name+"/apc_IOfP_"+model_name+"_p"+i+"_a0.8.nc")  
  system("/bin/rm -f "+ddir+model_name+"/apc_IOsP_"+model_name+"_p"+i+"_a0.8.nc") 

  ph = i+1
  ff = addfile(ddir+model_name+"/apc_IOfP_"+model_name+"_p"+i+"_a0.8.nc" ,"c")  
  ff->apc = apc_f

  fs = addfile(ddir+model_name+"/apc_IOsP_"+model_name+"_p"+i+"_a0.8.nc" ,"c")  
  fs->apc = apc_s