from eofs.standard import Eof
from sklearn.decomposition import PCA
import numpy, xarray
import numpy as np

def PCA_analysis(model_data,neofs,weights=None,lead_dim_name='lead',ens_dim_name='member',time_dim_name='time',space_dim_names=['lat','lon']):

    # move model data to [space x ensemble-lead-time] dimensions
    y_ec = model_data.stack(space=(space_dim_names),ec=(lead_dim_name,ens_dim_name,time_dim_name)).values.transpose('space','ec')

    pca = PCA(n_components=neofs)
    y_ec_pca = pca.fit(y_ec.T)  #transpose so have dims [ec,space]
    eof = y_ec_pca.components_
    var = y_ec_pca.explained_variance_ratio_
    eva = y_ec_pca.explained_variance_

    env = np.arange(1,neofs + 1,1)
    lat = model_data[space_dim_names[0]]
    lon = model_data[space_dim_names[1]]

    eofs = xarray.DataArray(
        np.zeros(( np.size(env), np.size(lat), np.size(lon) )),  # Fill with zeros initially
        dims=["env", space_dim_names[0], space_dim_names[1]],
        coords={"env": env, space_dim_names[0]: lat, space_dim_names[1]: lon},
    )

    eofs_stacked = 0.0 * eofs.rename('eof').stack(space=(space_dim_names))
    eofs_stacked.values = eof
    eofs = eofs_stacked.unstack()
    eofs.attrs['eval']  = eva
    eofs.attrs['pcvar'] = var

    return eofs

lats = 10
latn = 80
lonw = 120
lone = 240
ddir = '/home/sunming/data5/cuixy/Subpre_NPJ/data/'

#read in variables.
f = xarray.open_dataset(ddir+'/ecmwf/ecmwf_pf_anom_u200_r1.5.nc')
u = f['u'].loc[:,:,:,lats:latn,lonw:lone]

u = f['u'].loc[:,'1982-12-01':'2022-03-31',:,lats:latn,lonw:lone]
u = u.sel(time=u.time.dt.month.isin([12, 1, 2, 3]))
print(u.shape)
#u = u[:2,:2,:2,:,:]  #test

#add cosin weights.
coslat = numpy.cos(numpy.deg2rad(u.lat.values))[numpy.newaxis,numpy.newaxis,numpy.newaxis,:,numpy.newaxis]
U = u*numpy.sqrt(coslat)

#EOF
eof = PCA_analysis(U,neofs=20,weights=None,lead_dim_name="lead_time",ens_dim_name="number")

#store eof values.
eof.to_netcdf(ddir+'/ecmwf/eof_ecmwf_u200.nc')