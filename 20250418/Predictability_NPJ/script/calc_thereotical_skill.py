import numpy,xarray, pickle, os,datetime
from math import sqrt, log, exp

def SNR(ens,x):
    '''
     input
      ens [time,lat,lon]
      x [number,time,lat,lon]
    '''
    # first move model data to [space x ensembe-lead-time] dimensions
    ens_re = ens.stack(space=(['lat','lon'])).transpose('space','time').values

    #signal
    S = numpy.dot(ens_re,ens_re.T).diagonal()/(ens.shape[0]-1)

    #noise
    x_re = x.stack(space=(['lat','lon'])).transpose('number','time','space').values
    N = numpy.zeros_like(S)
    for j in range(numpy.size(N)):
        N[j] = numpy.mean([( (x_re[:,i,j]-ens_re[j,i])[numpy.newaxis,:] @
                (x_re[:,i,j]-ens_re[j,i])[numpy.newaxis,:].T ) for i in range(x_re.shape[1])], axis=0)/(x_re.shape[0]-1)
        
    #get skill
    acc = numpy.zeros_like(ens_re[:,0])
    for i in range(numpy.size(acc)):
        snr = S[i]/N[i]
        acc[i] = snr/sqrt((snr+1)*snr)

    return acc    

def MI(ens,x):
    '''
        mutual-information based predictability:
    '''

    #signal
    print('signal begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
    ens_re = ens.stack(space=(['lat','lon'])).transpose('time','space').values
    F = numpy.dot(ens_re.T,ens_re).diagonal()/(ens_re.shape[0]-1)

    #noise
    print('noise begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
    x_re = x.stack(space=(['lat','lon'])).transpose('number','time','space').values
    #E = numpy.mean([numpy.cov(x_re[:,i,:].transpose()).diagonal() for i in range(x_re.shape[1])], axis=0)
    E = numpy.zeros_like(F)
    N = numpy.zeros_like(F)
    for j in range(numpy.size(E)):
        E[j] = numpy.mean([log( ( (x_re[:,i,j]-x_re[:,i,j].mean())[numpy.newaxis,:] @
                (x_re[:,i,j]-x_re[:,i,j].mean())[numpy.newaxis,:].T ) /(x_re.shape[0]-1) )  for i in range(x_re.shape[1])], axis=0)

        N[j] = numpy.mean([( (x_re[:,i,j]-x_re[:,i,j].mean())[numpy.newaxis,:] @
                (x_re[:,i,j]-x_re[:,i,j].mean())[numpy.newaxis,:].T ) for i in range(x_re.shape[1])], axis=0)/(x_re.shape[0]-1)

    T = N + F
    # get mutual information.
    acc = numpy.zeros_like(F)
    for i in range(numpy.size(F)):
        MI = (log(T[i])-E[i])*0.5
        acc[i] = sqrt(1-exp(-2.*MI))

    return acc 


lats = 10
latn = 80
lonw = 100
lone = 240
ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data/"

#read in data.
print('read in u begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
f = xarray.open_dataset(ddir+'/ecmwf/ecmwf_pf_anom_u200_r1.5.nc')
u = f['u'].loc[:,'1982-12-01':'2022-03-31',:,lats:latn,lonw:lone]
u = u.sel(time=u.time.dt.month.isin([12, 1, 2, 3]))

#get ens mean.
if os.path.exists(ddir+'/ecmwf/ecmwf_ens_mean_u200_r1.5.nc'):
    f = xarray.open_dataset(ddir+'/ecmwf/ecmwf_ens_mean_u200_r1.5.nc')
    ens = f['u'].loc['1982-12-01':'2022-03-31',:,lats:latn,lonw:lone]
    ens = ens.sel(time=ens.time.dt.month.isin([12, 1, 2, 3]))

else:
    ens = u.mean('number')
    ens.to_netcdf(ddir+'/ecmwf/ecmwf_ens_mean_u200_r1.5.nc')

# acc begin.
print('begin. '+datetime.datetime.now().strftime('%m-%d %H:%M'))
lat = u.lat; lon = u.lon; lead_time = u.lead_time
acc = numpy.zeros([numpy.size(lead_time.values),numpy.size(lat.values)*numpy.size(lon.values)])

for i in numpy.arange(0,numpy.size(lead_time.values),1):
    print(i)
    acc[i,:] = SNR(ens[:,i,:,:],u[:,:,i,:,:])
    #acc[i,:] = MI(ens[:,i,:,:],u[:,:,i,:,:])
print('end. '+datetime.datetime.now().strftime('%m-%d %H:%M'))

with open(ddir+'/ecmwf_SNRskill', 'wb') as f2:
    pickle.dump(acc,f2)

#reshape it.
ac = acc.reshape(numpy.size(lead_time.values),numpy.size(lat.values),numpy.size(lon.values))
ds = xarray.Dataset({'acc': (('lead_time','lat','lon'), ac)},coords={'lead_time': lead_time,'lat': lat,'lon': lon})
ds.to_netcdf(ddir+"ecmwf_SNRskill.nc")