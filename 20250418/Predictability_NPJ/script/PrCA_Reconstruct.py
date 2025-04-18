import numpy, xarray, pickle, datetime
from sklearn.preprocessing import StandardScaler

'''
    reconstruct u with several MPMs. 
    Both model and obs.
    like X = p1v1 + p2v2 + ... + pnvn
'''

lats = 10
latn = 80
lonw = 100
lone = 240
ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data/"

#read in MPMs. [mode lat lon]
n = 4
P = xarray.open_dataarray(ddir+'PrCA_patterns.nc').loc[:n,:,:]
Pt = P.stack(space=(['lat','lon'])).values  #[mode,space]

#reconstruct for obs. [mode,lead_time,time]
print('reconstruct obs begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
v = xarray.open_dataarray(ddir+'/PrCA_obs_variates.nc').loc[:n,:,:]
vt = v.stack(ec=('lead_time','time')).values
#with open(ddir+'/PrCA_obs_variates_std.pkl','rb') as f:
#    transfer = pickle.load(f)
#v_std = transfer.inverse_transform(vt)
RC_obs = vt.T @ Pt                

# coordinate arrays for the output data
lat  = P['lat']
lon  = P['lon']
time = v['time']
lead = v['lead_time']

# Create an empty DataArray with the specified dimensions and coordinates
RC = xarray.DataArray(
        numpy.zeros(( numpy.size(time), numpy.size(lead), numpy.size(lat), numpy.size(lon) )),  
        dims=["time", "lead_time", "lat", "lon"],
        coords={"time": time, "lead_time": lead, "lat": lat, "lon": lon},
    )

# arrange the patterns in an xarray
RC_stacked = 0.0 * RC.rename('u').stack(ec=(['lead_time','time']),space=(['lat','lon']))
RC_stacked.values = RC_obs
RC = RC_stacked.unstack().transpose('time','lead_time','lat','lon')
RC.to_netcdf(ddir+'ecmwf/ecmwf_obs_u200_rc.nc')

#reconstruct for model. [lead_time,number,time,mode]
print('reconstruct model begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
v = xarray.open_dataarray(ddir+'/PrCA_variates.nc').loc[:,:,:,:n].mean('number')
vt = v.stack(ec=('lead_time','time')).transpose('mode','ec').values
RC_model = vt.T @ Pt

# Create an empty DataArray with the specified dimensions and coordinates
RC = xarray.DataArray(
        numpy.zeros(( numpy.size(time), numpy.size(lead), numpy.size(lat), numpy.size(lon) )),  
        dims=["time", "lead_time", "lat", "lon"],
        coords={"time": time, "lead_time": lead, "lat": lat, "lon": lon},
    )

# arrange the patterns in an xarray
RC_stacked = 0.0 * RC.rename('u').stack(ec=(['lead_time','time']),space=(['lat','lon']))
RC_stacked.values = RC_model
RC = RC_stacked.unstack().transpose('time','lead_time','lat','lon')
RC.to_netcdf(ddir+'ecmwf/ecmwf_pf_u200_rc.nc')