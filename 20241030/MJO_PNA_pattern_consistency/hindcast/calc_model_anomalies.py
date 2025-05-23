import xarray, sys, datetime

model_name = str(sys.argv[1])

ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/hindcast"

f1 = xarray.open_dataset(ddir+"/"+model_name+"_cf_z500_mean_r2.5.nc")  #control forecast.
f2 = xarray.open_dataset(ddir+"/"+model_name+"_pf_z500_mean_r2.5.nc")  #perturbed forecast.

xc = f1['gh']
xc = xc.isel(lead_time=range(1, len(xc.lead_time)))
xp = f2['gh']
xp = xp.isel(lead_time=range(1, len(xp.lead_time)))

xc_clim = xc.assign_coords(dayofyear=xc.initial_time.dt.strftime("%m-%d")).groupby('dayofyear').mean()
xp_clim = xp.assign_coords(dayofyear=xp.initial_time.dt.strftime("%m-%d")).groupby('dayofyear').mean()

xc_clim = xc_clim.assign_coords(number=0)
x_clim = xarray.concat([xc_clim,xp_clim], dim="number")
x_clim = x_clim.mean("number")
print(model_name+"data clim done.")

xc_anom = xc.assign_coords(dayofyear=xc.initial_time.dt.strftime("%m-%d")).groupby('dayofyear')  - x_clim
xp_anom = xp.assign_coords(dayofyear=xp.initial_time.dt.strftime("%m-%d")).groupby('dayofyear')  - x_clim
print(model_name+"data anomaly done.")

xc_anom_final = xc_anom
xp_anom_final = xp_anom
#remove previous 120 days.
clim = xarray.open_dataarray(ddir+"/"+model_name+"_z500_daily_anom_1979-2022_r2.5.nc")
clim['time'] = clim.indexes['time'].normalize()
clim /= 9.8   #unit: gpm
#remove previous 120 days.
for start_time in range(0,len(xc_anom.initial_time)):
    for indice in range(0,len(xc_anom.lead_time)):
        forecast_time = datetime.datetime.strptime(str(xc_anom.initial_time[start_time].dt.strftime("%Y-%m-%d").to_numpy()),"%Y-%m-%d")+datetime.timedelta(days=indice)
        clim_t = clim.sel(time=[(forecast_time+datetime.timedelta(days=ftime)) for ftime in range(-119,1)])
        if (indice > 0):
            clim_t[120-indice:,:,:] = xc_anom[start_time,:indice,:,:]
        xc_anom_final[start_time,indice,:,:] = xc_anom[start_time,indice,:,:] - clim_t.mean('time')
        for ens in range(0,len(xp_anom.number)):
            if (indice > 0):
                clim_t[120-indice:,:,:] = xp_anom[ens,start_time,:indice,:,:]
            xp_anom_final[ens,start_time,indice,:,:] = xp_anom[ens,start_time,indice, :, :] - clim_t.mean('time')
print("remove previous 120 days done.")

xc_anom_final.to_netcdf(ddir+"/"+model_name+"_cf_z500_anom_r2.5.nc",engine="netcdf4",format="NETCDF4")
xp_anom_final.to_netcdf(ddir+"/"+model_name+"_pf_z500_anom_r2.5.nc",engine="netcdf4",format="NETCDF4")
print("store anomalies done.")