import xarray, sys, datetime

Year = str(sys.argv[1])

ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/hindcast/ecmwf"

f1 = xarray.open_dataset(ddir+"/"+Year+"/data/ecmwf_cf_z500_mean_r2.5.nc")  #control forecast.
f2 = xarray.open_dataset(ddir+"/"+Year+"/data/ecmwf_pf_z500_mean_r2.5.nc")  #perturbed forecast.
clim = xarray.open_dataarray(ddir+"/"+Year+"/data/ecmwf_z500_daily_anom_1979-2022_r2.5.nc")


xc_anom = f1['gh']
xc_anom = xc_anom.isel(lead_time=range(1, len(xc_anom.lead_time)))
xp_anom = f2['gh']
xp_anom = xp_anom.isel(lead_time=range(1, len(xp_anom.lead_time)))

print(datetime.datetime.now())
for start_time in range(0,len(xc_anom.initial_time)):
    for indice in range(0,len(xc_anom.lead_time)):
        forecast_time = datetime.datetime.strptime(str(xc_anom.initial_time[start_time].dt.strftime("%Y-%m-%d").to_numpy()),"%Y-%m-%d")+datetime.timedelta(days=indice)
        if (xc_anom.initial_time[start_time].dt.year==2022 and xc_anom.initial_time[start_time].dt.month>=10):
            i = 0
        else:
        	i = 1