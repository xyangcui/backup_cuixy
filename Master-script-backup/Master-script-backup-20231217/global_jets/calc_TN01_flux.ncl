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

  dimNames = (/"time", "latitude", "longitude"/)  
  dimSizes = (/ -1   ,  dimsizes(LAT),  dimsizes(LON)/) 
  dimUnlim = (/ True , False, False/)   
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout, "time" ,typeof(TIME),getvardims(TIME))                           
  filevardef(fout, "latitude"  ,typeof(LAT),getvardims(LAT))                          
  filevardef(fout, "longitude"  ,typeof(LON),getvardims(LON))                          
  filevardef(fout, VAR1    ,typeof(x)  ,getvardims(x))   
  filevardef(fout, VAR2    ,typeof(y)  ,getvardims(y))   
                                                                                             
  filevarattdef(fout,"time" ,TIME)                    
  filevarattdef(fout,"latitude"  ,LAT)                   
  filevarattdef(fout,"longitude"  ,LON)                     
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
zdir = ddir + "/ERA5_daily_phi250_1979-2022_all.nc"
udir = ddir+"/"+"ERA5_daily_u250_1979-2022_all.nc"
vdir = ddir+"/"+"ERA5_daily_v250_1979-2022_all.nc"
outdir = ddir+"/"+"ERA5_daily_waf_1979-2022.nc"
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
lat  = z&latitude
lon  = z&longitude

TIME = cd_calendar(time,0)
year = toint(TIME(:,0))
mon  = toint(TIME(:,1))
day  = toint(TIME(:,2))

ddd     = day_of_year(year, mon, day) 
yyyyddd = year*1000 + ddd  

hClmDay = clmDayTLL(z, yyyyddd)  
psi = calcDayAnomTLL(z, yyyyddd, hClmDay)                                       
                              
be = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
en = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)
u := u(be:en,:,:)
v := v(be:en,:,:)
time1 = u&time
TIME1  = cd_calendar(time1, 0)
mon1   = toint(TIME1(:,1))
ntDJF  = ind(mon1.eq.1.or.mon1.eq.2.or.mon1.eq.12)
u := u(ntDJF,:,:)
v := v(ntDJF,:,:)

um = dim_avg_n_Wrap(u, 0)
vm = dim_avg_n_Wrap(u, 0)
print("DJF mean done. "+systemfunc("date") )

WAFu = new(dimsizes(z),typeof(z),z@_FillValue)
WAFv = new(dimsizes(z),typeof(z),z@_FillValue)
WAFu = 0.
WAFv = 0.
;cori parameter.
prem  = 25000./100000.
nlats = ind(lat.eq.-15)
nlatn = ind(lat.eq.15)
nlat  = dimsizes(lat)
nlon  = dimsizes(lon)

dlat = abs(lat(1)-lat(0))
dlon = abs(lon(1)-lon(0))

t = 0
do i = 0, dimsizes(time)-1
  waf = new((/nlat,nlon,2/),float)
  waf = 0.
  tnf_xy(psi(i,:,:),um,vm,dlat,dlon,nlat,nlon,prem,nlats,nlatn,waf)
  WAFu(i,:,:) = waf(:,:,0)
  WAFv(i,:,:) = waf(:,:,1)
  delete(waf)
  t = t + 1
end do
print(t)
;nlat = ind(lat.gt.-10.and.lat.lt.10)
;do i = 0, dimsizes(time)-1
    ;WAFuv = cal_WFsf(psi(i,:,:),um,vm)
    ;WAFuv(:,nlat,:) = 0.
    ;do j=0,dimsizes(lat)-1
    ;    do k = 0,dimsizes(lon)-1
    ;        if(um(j,k).lt.1.)then
    ;            WAFuv(:,j,k) = 0.
    ;        end if
    ;    end do
    ;end do
    ;WAFu(i,:,:) = WAFuv(0,:,:)*prem
    ;WAFv(i,:,:) = WAFuv(1,:,:)*prem
    ;delete(WAFuv)
;end do
WAFu@long_name = "zonal wave flux"
WAFu@unit = "m^2/s^2"
WAFv@long_name = "meridional wave flux"
WAFv@unit = "m^2/s^2"
copy_VarCoords(z, WAFu)
copy_VarCoords(z, WAFv)
delete([/z,um,vm/])
create_nc_file(WAFu,WAFv,time,lat,lon,outdir,"u","v")
print("calc done. "+systemfunc("date"))