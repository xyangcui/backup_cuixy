import xarray, sys, datetime

model_name = str(sys.argv[1])

ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/hindcast"

f1 = xarray.open_dataset(ddir+"/"+model_name+"_cf_z500_mean_r2.5.nc")  #control forecast.
f2 = xarray.open_dataset(ddir+"/"+model_name+"_pf_z500_mean_r2.5.nc")  #perturbed forecast.

xc = f1['gh']
xc = xc.isel(lead_time=range(1, len(xc.lead_time)))
xp = f2['gh']
xp = xp.isel(lead_time=range(1, len(xp.lead_time)))

xc = xc.assign_coords(number=0)

print(xp)

clim = xarray.open_dataarray(ddir+"/"+model_name+"_z500_daily_anom_1979-2022_r2.5.nc")
clim['time'] = clim.indexes['time'].normalize()

clim /= 9.8


x_clim = xarray.concat([xc,xp], dim="number")
x_clim = x_clim.mean("number")
print(x_clim)