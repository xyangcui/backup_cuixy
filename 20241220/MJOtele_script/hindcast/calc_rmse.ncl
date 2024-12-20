;------------------------------------------------
;Function RMSE  single point.
;------------------------------------------------
function calc_RMSE(mvalue,ovalue,itime,ltime)
local mvalue,ovalue,cLat,cLaT,apc,apc1,apc2,apc3,i,j
begin

  rmse = new((/dimsizes(ltime)/),float)
  do i = 0,dimsizes(ltime)-1
    mv = dim_standardize(mvalue(:,i),0)
    ov = dim_standardize(ovalue(:,i),0)
    ;mv = mvalue(:,i)
    ;ov = ovalue(:,i)
    rmse1 = mv-ov
    rmse2 = dim_avg(rmse1*rmse1)
    rmse(i) = sqrt(rmse2)
    delete([/mv,ov,rmse1,rmse2/])
  end do

  return rmse
end

lats = 10
latn = 80
lonw = 140
lone = 270

pi = atan(1)*4.
model_name = "ecmwf"

ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/hindcast/"

fsdir = ddir+model_name+"/"+Year+"/"+"IOsP_p"+i+"_a1.nc"
obsdir = ddir+model_name+"/"+Year+"/"+"IOsP_ob_p"+i+"_a1.nc"

ffdir = ddir+model_name+"/"+Year+"/"+"IOfP_p"+i+"_a1.nc"
obfdir = ddir+model_name+"/"+Year+"/"+"IOfP_ob_p"+i+"_a1.nc"

  fs = addfile(fsdir,"r")
  fos = addfile(obsdir,"r")

  ff = addfile(ffdir,"r")
  fof = addfile(obfdir,"r")

  pf_z_s = fs->gh(:,:,:,{lats:latn},{lonw:lone})
  ob_z_s = fos->z(:,:,{lats:latn},{lonw:lone})

  pf_z_f = ff->gh(:,:,:,{lats:latn},{lonw:lone})
  ob_z_f = fof->z(:,:,{lats:latn},{lonw:lone})

;ens mean
  pf_z_s := dim_avg_n_Wrap(pf_z_s, 0)
  pf_z_f := dim_avg_n_Wrap(pf_z_f, 0)

;calc ACC single point
  nlat = dimsizes(pf_z_s&lat)
  nlon = dimsizes(pf_z_s&lon)

  rmse_s = new((/dimsizes(pf_z_s&lead_time),nlat,nlon/),float)
  rmse_f = new((/dimsizes(pf_z_s&lead_time),nlat,nlon/),float)

  do k = 0, nlat-1
    do j = 0, nlon-1
      rmse_s(:,k,j) = calc_RMSE(pf_z_s(:,:,k,j),ob_z_s(:,:,k,j),pf_z_s&initial_time,pf_z_s&lead_time)
      rmse_f(:,k,j) = calc_RMSE(pf_z_f(:,:,k,j),ob_z_f(:,:,k,j),pf_z_f&initial_time,pf_z_f&lead_time)
    end do
  end do

  copy_VarCoords(pf_z_s(0,:,:,:), rmse_s)
  copy_VarCoords(pf_z_f(0,:,:,:), rmse_f)

  system("/bin/rm -f "+ddir+model_name+"/"+Year+"/rmse/rmse_IOfP_"+"p"+i+"_a1.nc")  
  system("/bin/rm -f "+ddir+model_name+"/"+Year+"/rmse/rmse_IOsP_"+"p"+i+"_a1.nc") 

  ff = addfile(ddir+model_name+"/"+Year+"/rmse/rmse_IOfP"+"_p"+i+"_a1.nc" ,"c")  
  ff->rmse = rmse_f

  fs = addfile(ddir+model_name+"/"+Year+"/rmse/rmse_IOsP"+"_p"+i+"_a1.nc" ,"c")  
  fs->rmse = rmse_s