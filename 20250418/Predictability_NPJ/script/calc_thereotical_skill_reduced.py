import numpy,xarray, pickle, os,datetime
from math import sqrt

def SNR(ens,x,E):
    '''
     input
      ens [time,lat,lon]
      x [number,time,lat,lon]
    '''
    # move eofs to same to dimensions [space x mode]
    E = eofs.stack(space=(['lat','lon'])).transpose('space','evn')

    # first move model data to [space x ensembe-lead-time] dimensions
    y_ec = ens.stack(space=(['lat','lon'])).transpose('space','time')
    # preallocate xarray with dimensions we want: space reduced to Neofs
    ens_re = 0.0 * y_ec[:numpy.size(E.evn),:]
    # now reduce dimensionality by projecting EOFs onto the original data [evn,time]
    ens_re = numpy.linalg.inv(E.values.transpose() @ E.values) @ E.values.transpose() @ y_ec.values

    #signal
    S = numpy.dot(E,numpy.dot(numpy.dot(ens_re,ens_re.T),E.T)).diagonal()/(ens.shape[0]-1)

    # first move model data to [space x ensembe-lead-time] dimensions
    y_ec = x.stack(space=(['lat','lon']),ec=(['time','number'])).transpose('space','ec')
    # preallocate xarray with dimensions we want: space reduced to Neofs
    x_re = 0.0 * y_ec[:numpy.size(E.evn),:]
    # now reduce dimensionality by projecting EOFs onto the original data [evn,timexnumber]
    x_re.values = numpy.linalg.inv(E.values.transpose() @ E.values) @ E.values.transpose() @ y_ec.values

    #noise
    x_re = x_re.unstack('ec').transpose('number','time','space').values
    Ntmp = numpy.mean([numpy.cov(x_re[:,i,:].transpose()) for i in range(x_re.shape[1])], axis=0)
    N = numpy.dot(E,numpy.dot(Ntmp,E.T)).diagonal()

    #get skill
    acc = numpy.zeros_like(E[:,0])
    for i in range(numpy.size(acc)):
        snr = S[i]/N[i]
        acc[i] = snr/sqrt((snr+1)*snr)

    return acc    

lats = 10
latn = 80
lonw = 100
lone = 240
ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data/"

#read in data.
print('read in u begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
if os.path.exists(ddir+'/ecmwf/ecmwf_pf_anom_u200_r1.5_weighted.nc'):
    #weighted file.
    f = xarray.open_dataset(ddir+'/ecmwf/ecmwf_pf_anom_u200_r1.5_weighted.nc')
    u = f['u'].loc[:,'1982-12-01':'2022-03-31',:,lats:latn,lonw:lone]
    u = u.sel(time=u.time.dt.month.isin([12, 1, 2, 3]))

else:
    #unweighted file and need to be weighted.
    f = xarray.open_dataset(ddir+'/ecmwf/ecmwf_pf_anom_u200_r1.5.nc')
    u = f['u'].loc[:,'1982-12-01':'2022-03-31',:,lats:latn,lonw:lone]
    u = u.sel(time=u.time.dt.month.isin([12, 1, 2, 3]))
    coslat = numpy.cos(numpy.deg2rad(u.lat.values))[numpy.newaxis,numpy.newaxis,numpy.newaxis,:,numpy.newaxis]
    u = u*numpy.sqrt(coslat)
    u.to_netcdf(ddir+'/ecmwf/ecmwf_pf_anom_u200_r1.5_weighted.nc')


#get ens mean.
if os.path.exists(ddir+'/ecmwf/ecmwf_ens_mean_u200_r1.5.nc'):
    f = xarray.open_dataset(ddir+'/ecmwf/ecmwf_ens_mean_u200_r1.5.nc')
    ens = f['u'].loc['1982-12-01':'2022-03-31',:,lats:latn,lonw:lone]
    ens = ens.sel(time=ens.time.dt.month.isin([12, 1, 2, 3]))

else:
    ens = u.mean('number')
    ens.to_netcdf(ddir+'/ecmwf/ecmwf_ens_mean_u200_r1.5.nc')

# weight.
if os.path.exists(ddir+'/ecmwf/ecmwf_ens_mean_u200_r1.5_weighted.nc'):
    #weighted file.
    f = xarray.open_dataset(ddir+'/ecmwf/ecmwf_ens_mean_u200_r1.5_weighted.nc')
    ens = f['u'].loc['1982-12-01':'2022-03-31',:,lats:latn,lonw:lone]
    ens = ens.sel(time=u.time.dt.month.isin([12, 1, 2, 3]))

else:
    #unweighted file and need to be weighted.
    coslat = numpy.cos(numpy.deg2rad(ens.lat.values))[numpy.newaxis,numpy.newaxis,:,numpy.newaxis]
    ens.values = ens.values*numpy.sqrt(coslat)
    ens.to_netcdf(ddir+'/ecmwf/ecmwf_ens_mean_u200_r1.5_weighted.nc')

#read in eofs.
print('read in eofs begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
n = 50
eofs = xarray.open_dataarray(ddir+'/eof_u200.nc')
eof = eofs.loc[:n,lats:latn,lonw:lone]

# acc begin.
print('begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
lat = u.lat; lon = u.lon; lead_time = u.lead_time
acc = numpy.zeros([numpy.size(lead_time.values),numpy.size(lat.values)*numpy.size(lon.values)])

for i in numpy.arange(0,numpy.size(lead_time.values),1):
    print(i)
    acc[i,:] = SNR(ens[:,i,:,:],u[:,:,i,:,:],eof)
print('end. '+datetime.datetime.now().strftime('%m-%d %H:%M'))

#reshape it.

ac = acc.reshape(numpy.size(lead_time.values),numpy.size(lat.values),numpy.size(lon.values))
ds = xarray.Dataset({'acc': (('lead_time','lat','lon'), ac)},coords={'lead_time': lead_time,'lat': lat,'lon': lon})
ds.to_netcdf(ddir+"ecmwf_SNRskill.nc")