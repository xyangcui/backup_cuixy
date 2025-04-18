import numpy,xarray,pickle,datetime,os
from sklearn.preprocessing import StandardScaler
import numpy as np
ddir = '/home/sunming/data5/cuixy/Subpre_NPJ/data/'

lats = 10
latn = 80
lonw = 100
lone = 240
# read in variables.
print('read in u begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
f = xarray.open_dataset(ddir+'/ecmwf/ecmwf_obs_anom_u200_v2.nc')
u = f['u'].loc['1982-12-01':'2022-03-31',:,lats:latn,lonw:lone]
u = u.sel(time=u.time.dt.month.isin([12, 1, 2, 3]))

# read in predictable patterns.
print('read in Pr patterns begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
x = xarray.open_dataarray(ddir+'/PrCA_patterns.nc')

# get variates in obs.
y_ec = u.stack(space=(['lat','lon']),ec=('lead_time','time'))

# move eofs to same to dimensions [space x mode]
E = x.stack(space=(['lat','lon'])).transpose('space','mode')

# preallocate xarray with dimensions we want: space reduced to Neofs
y_ec_reduced = 0.0 * y_ec[:np.size(E.mode),:]

# now reduce dimensionality by projecting EOFs onto the original data ['mode','time']
y_ec_reduced.values = np.linalg.inv(E.values.transpose() @ E.values) @ E.values.transpose() @ y_ec.values


# standardize it.
if os.path.exists(ddir+'/PrCA_obs_variates_std.pkl'):
    with open(ddir+'/PrCA_obs_variates_std.pkl','rb') as f:
        transfer = pickle.load(f)
    y_ec_std = transfer.transform(y_ec_reduced.values)

else:
    transfer = StandardScaler()
    u_transfer = transfer.fit(y_ec_reduced.values)
    with open(ddir+'/PrCA_obs_variates_std.pkl','wb') as f:
        pickle.dump(u_transfer, f)
    y_ec_std = u_transfer.transform(y_ec_reduced.values)

# store data.
print('store data begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
mode = np.arange(1,np.size(y_ec_reduced.space) + 1,1); time = u['time']; lead = u['lead_time']

y_ec_ori = xarray.DataArray(
        np.zeros(( np.size(mode), np.size(time), np.size(lead) )),  # Fill with zeros initially
        dims=['mode', 'time', 'lead_time'],
        coords={'mode': mode, 'time': time, 'lead_time': lead},
    )

y_ec_ori_stacked = 0.0 * y_ec_ori.rename('v').stack(ec=('lead_time','time')).transpose('mode','ec')
y_ec_ori_stacked.values = y_ec_std
y_ec_ori = y_ec_ori_stacked.unstack()

y_ec_ori.to_netcdf(ddir+'/PrCA_obs_variates.nc')