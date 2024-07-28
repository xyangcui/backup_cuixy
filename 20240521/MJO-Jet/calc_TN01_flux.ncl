load "/home/sunming/data5/cuixy/script/ncl_func/cal_WFsf.ncl"
;----------------------------------------------
;; calc TN-01 flux.
;----------------------------------------------
tday = 24.*60.*60.
pi   = atan(1.0)*4.
;------------------------------------------------------------------
undef("create_nc_file")

procedure create_nc_file(x:float,y:float,TIME,LAT,LON,DIR:string,VAR1:string,VAR2:string)
local x,LON,LAT,TIME,DIR,VAR1,VAR2,fout,y
begin
  system("/bin/rm -f "+DIR)

  fout = addfile(DIR,"c")

  dimNames = (/"time", "lat", "lon"/)  
  dimSizes = (/ -1   ,  dimsizes(LAT),  dimsizes(LON)/) 
  dimUnlim = (/ True , False, False/)   
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout, "time" ,typeof(TIME),getvardims(TIME))                           
  filevardef(fout, "lat"  ,typeof(LAT),getvardims(LAT))                          
  filevardef(fout, "lon"  ,typeof(LON),getvardims(LON))                          
  filevardef(fout, VAR1    ,typeof(x)  ,getvardims(x))   
  filevardef(fout, VAR2    ,typeof(y)  ,getvardims(y))   
                                                                                             
  filevarattdef(fout,"time" ,TIME)                    
  filevarattdef(fout,"lat"  ,LAT)                   
  filevarattdef(fout,"lon"  ,LON)                     
  filevarattdef(fout,VAR1   ,x) 
  filevarattdef(fout,VAR2   ,y)   

  setfileoption(fout,"DefineMode",False)

  fout->time         = (/TIME/)     
  fout->latitude     = (/LAT/)
  fout->longitude    = (/LON/) 
  fout->$VAR1$       = (/x/)
  fout->$VAR2$       = (/y/)
  
end
;------------------------------------------------------------
ddir = "/home/sunming/data5/cuixy/global_jets/data"
zdir = ddir + "/ERA5_daily_anom_phi250_1979-2022_2.5x2.5.nc"
udir = ddir+"/"+"ERA5_daily_u250_1979-2022_r2.5.nc"
vdir = ddir+"/"+"ERA5_daily_v250_1979-2022_r2.5.nc"
outdir = ddir+"/"+"ERA5_daily_waf_1979-2022_r2.5.nc"
;----------------------------------------------
; main code.
;----------------------------------------------
fz = addfile(zdir,"r")
fu = addfile(udir,"r")
fv = addfile(vdir,"r")

z = fz->phi(:,{90:-35},:)
u = fu->u(:,{250},{90:-35},:)
v = fv->v(:,{250},{90:-35},:)

time = z&time
lat  = z&lat
lon  = z&lon

TIME = cd_calendar(time,0)
year = toint(TIME(:,0))
mon  = toint(TIME(:,1))
day  = toint(TIME(:,2))

ddd     = day_of_year(year, mon, day) 
yyyyddd = year*1000 + ddd  

uclim = clmDayTLL(u, yyyyddd)
vclim = clmDayTLL(v,yyyyddd)

uclim2 = uclim
vclim2 = vclim
do i = 59,364
  uclim2(i,:,:) = uclim(i+1,:,:)
  vclim2(i,:,:) = vclim(i+1,:,:)
end do
;; convert to the same bnds.
um = new(dimsizes(u),typeof(u))
vm = new(dimsizes(v),typeof(v))
um = 0.
vm = 0.

do iyear = 1979,2022
  nt = ind(year.eq.iyear)
  if (isleapyear(iyear)) then
    um(nt,:,:) = uclim
    vm(nt,:,:) = vclim
  else 
    um(nt,:,:) = uclim2(0:364,:,:)
    vm(nt,:,:) = vclim2(0:364,:,:)
  end if
  delete(nt)
end do

print("clim U and V done. "+systemfunc("date") )

WAFu = new(dimsizes(z),typeof(z),z@_FillValue)
WAFv = new(dimsizes(z),typeof(z),z@_FillValue)
WAFu = 0.
WAFv = 0.
;cori parameter.
prem  = 25000./100000.
nlats = ind(lat.eq.-10)
nlatn = ind(lat.eq.10)
nlat  = toint(dimsizes(lat))
nlon  = toint(dimsizes(lon))

dlat = doubletofloat(abs(lat(1)-lat(0)))
dlon = doubletofloat(abs(lon(1)-lon(0)))

;t = 0
;do i = 0, dimsizes(time)-1
;  waf = new((/nlat,nlon,2/),float)
;  waf = 0.
;  tnf_xy(z(i,:,:),um(i,:,:),vm(i,:,:),dlat,dlon,nlat,nlon,prem,nlats,nlatn,waf)
;  WAFu(i,:,:) = waf(:,:,0)
;  WAFv(i,:,:) = waf(:,:,1)
;  delete(waf)
;  t = t + 1
;end do
nlati = ind(lat.gt.-10.and.lat.lt.10)

do i = 0, dimsizes(time)-1
  WAFuv = cal_WFsf(z(i,:,:),um(i,:,:),vm(i,:,:),lat,lon)
  WAFuv(:,nlati,:) = 0.
  WAFu(i,:,:) = WAFuv(0,:,:)*prem
  WAFv(i,:,:) = WAFuv(1,:,:)*prem
  delete(WAFuv)
end do
WAFu@long_name = "zonal wave flux"
WAFu@unit = "m^2/s^2"
WAFv@long_name = "meridional wave flux"
WAFv@unit = "m^2/s^2"
copy_VarCoords(u, WAFu)
copy_VarCoords(v, WAFv)
delete([/z,um,vm/])
create_nc_file(WAFu,WAFv,time,lat,lon,outdir,"u","v")
print("calc done. "+systemfunc("date"))