import numpy,xarray, pickle
from math import sqrt
from math import log
from math import exp

def SNR(x,C0,d):
    '''
    input
    x: (times,dims)
    '''
    #calculate signal F=<x x.T>
    F = numpy.matmul(x.T,x)/(x.shape[0]-1)
    #ACC version
    skill = numpy.zeros(d)
    for i in range(d):
        S2 = F[i,i]/(C0[i,i]-F[i,i])
        if S2 <0:
            skill[i] = 1.
        else:
            skill[i] = S2/sqrt((S2+1)*S2)
    
    return skill

def MI():
    '''
        mutual-information based predictability:
    '''

    # first move model data to [space x ensembe-lead-time] dimensions
    y_ec = model_data.stack(space=(space_dim_names),ec=(lead_dim_name,ens_dim_name,time_dim_name))

    # get the climatological covariance matrix
    sigma_inf = numpy.cov(y_ec)

    # get the forecast covariance matrix
        # preallocate noise covariance matrix (one for each initialization)
    sigma_tau = np.zeros([np.size(model_data[lead_dim_name]),np.shape(sigma_inf)[0],np.shape(sigma_inf)[1]])

    # iterate to compute the noise covariance (time mean ensemble variance) at each lead time (weighted by dtau)
    for j in range(0,np.size(model_data[lead_dim_name])):
        sigma_tau[j,:,:] = np.mean([np.cov(y_ec[j,:,i,:].transpose()) for i in range(np.size(model_data[time_dim_name]))], axis=0) 

    # get mutual information.
    acc = np.zeros_like(y_ec[:,].values)
    for i in range(dim1):
        MI1 = log(sigma_inf[i,i])
        for j in range(dim2):
            MI = (MI1-log(sigma_tau[j,:,:]))*0.5
            acc[j,i] = sqrt(1-exp(-2.*MI))

    return acc

def PCC(x,C0):
    '''
    input
    x: (times,dims)
    '''
    #calculate signal F=<x x.T>
    F = numpy.matmul(x.T,x)/(x.shape[0]-1)
    E = C0 - F
    #PCC version
    S2 = numpy.trace(F)/numpy.trace(E)

    return S2/sqrt((S2+1)*S2)

lats = 10
latn = 80
lonw = 100
lone = 240
ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data/"


f = xarray.open_dataset("/home/sunming/data5/cuixy/Subpre_NPJ/data/ecmwf/ecmwf_pf_anom_u200.nc")
u = f['u'].loc['1982-12-01':'2022-03-31',:,lats:latn,lonw:lone]
print("sort months.")
u = u.sel(time=u.time.dt.month.isin([12, 1, 2, 3]))
print("sort months done.")
lat = u.lat; lon = u.lon; lead_time = u.lead_time

X = u.to_numpy().reshape(u.shape[0],u.shape[1],-1)
C0 = numpy.dot(X[:,0,:].T,X[:,0,:])/(X.shape[0]-1)

with open(ddir+'/covar_obs_u200', 'wb') as f2:
    pickle.dump(C0,f2)

#with open(ddir+'/ecmwf_obs_X', 'rb') as f1:
#    X = pickle.load(f1)

##covariance
#with open(ddir+'/covar_obs_u200', 'rb') as f2:
#    C0 = pickle.load(f2)


print("skill begin.")
acc = numpy.zeros([X.shape[1],X.shape[2]])
for i in numpy.arange(0,X.shape[1],1):
    print(i)
    acc[i,:] = SNR(X[:,i,:],C0,X.shape[2])

#reshape it.
ac = acc.reshape(u.shape[1],u.shape[2],u.shape[3])
ds = xarray.Dataset({'acc': (('lead_time','lat','lon'), ac)},coords={'lead_time': lead_time,'lat': lat,'lon': lon})
ds.to_netcdf(ddir+"ecmwf_theskill.nc")