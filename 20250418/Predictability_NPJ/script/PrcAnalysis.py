import numpy,xarray,pickle,datetime,os
from apt import APT
'''
    Predictable component Analysis for U200 (NPJ).
    main features of this script.
    variables: ECMWF U200 in NPJ domain with [number,time,lead_time,lat,lon]
    eofs: select first 20 eofs explaining about 82% variance.
    
    Main function APT is copied from https://doi.org/10.5281/zenodo.13647.

'''

ddir = '/home/sunming/data5/cuixy/Subpre_NPJ/data/'

lats = 10
latn = 80
lonw = 100
lone = 240
#read in variables.
print('read in u begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
if os.path.exists(ddir+'/ecmwf/ecmwf_pf_anom_u200_r1.5_weighted.nc'):
    #weighted file.
    f = xarray.open_dataset(ddir+'/ecmwf/ecmwf_pf_anom_u200_r1.5_weighted.nc')
    U = f['u'].loc[:,'1982-12-01':'2022-03-31',:,lats:latn,lonw:lone]
    U = U.sel(time=U.time.dt.month.isin([12, 1, 2, 3]))

else:
    #unweighted file and need to be weighted.
    f = xarray.open_dataset(ddir+'/ecmwf/ecmwf_pf_anom_u200_r1.5.nc')
    u = f['u'].loc[:,'1982-12-01':'2022-03-31',:,lats:latn,lonw:lone]
    u = u.sel(time=u.time.dt.month.isin([12, 1, 2, 3]))
    coslat = numpy.cos(numpy.deg2rad(u.lat.values))[numpy.newaxis,numpy.newaxis,numpy.newaxis,:,numpy.newaxis]
    U = u*numpy.sqrt(coslat)
    U.to_netcdf(ddir+'/ecmwf/ecmwf_pf_anom_u200_r1.5_weighted.nc')



#read in eofs.
print('read in eofs begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
n = 20
eofs = xarray.open_dataarray(ddir+'/eof_u200.nc')
eof = eofs.loc[:n,lats:latn,lonw:lone]

#APT
print('PrCAnalysis begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
weights = numpy.sqrt(numpy.cos(numpy.deg2rad(u.lat.values))[numpy.newaxis,:,numpy.newaxis])

patterns,predictable_variates,apt,q_patterns,P_Espace = APT(U,eof,weights=weights,lead_dim_name='lead_time',ens_dim_name='number')

#store data.
print('store data begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
patterns.to_netcdf(ddir+'/PrCA_patterns.nc')
predictable_variates.to_netcdf(ddir+'/PrCA_variates.nc')
apt.to_netcdf(ddir+'/PrCA_apt.nc')
q_patterns.to_netcdf(ddir+'/PrCA_q_patterns.nc')

with open(ddir+'/PrCA_PEspace', 'wb') as f:
    pickle.dump(P_Espace,f)