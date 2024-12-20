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
model_name = "ecmwf"

ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/hindcast/"

fsdir  = ddir+model_name+"/"+Year+"/"+"IO_p"+i+"_a1.nc"
obsdir = ddir+model_name+"/"+Year+"/"+"IO_ob_p"+i+"_a1.nc"

  fs = addfile(fsdir,"r")
  fos = addfile(obsdir,"r")

  pf_z_s = fs->gh(:,:,:,{lats:latn},{lonw:lone})
  ob_z_s = fos->z(:,:,{lats:latn},{lonw:lone})

  pf_z_s := dim_avg_n_Wrap(pf_z_s, 0)

  clat = tofloat(cos(pi*ob_z_s&lat/180)) 
  apc_s = calc_APC(pf_z_s,ob_z_s,clat,pf_z_s&initial_time,pf_z_s&lead_time)

  copy_VarCoords(pf_z_s(:,:,0,0), apc_s)

  system("/bin/rm -f "+ddir+model_name+"/"+Year+"/apc/apc_IO_"+"p"+i+"_a1_ens.nc")  

  fs = addfile(ddir+model_name+"/"+Year+"/apc/apc_IO"+"_p"+i+"_a1_ens.nc" ,"c")  
  fs->apc = apc_s