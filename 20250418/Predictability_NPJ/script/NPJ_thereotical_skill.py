import pickle, xarray, numpy
from scipy.linalg import expm
from math import sqrt

def ACC(L,c0,C0,E,n,d,tau):
    '''
        input
         L:  LIM matrix
         C0: u200 covariance
         c0: state vector covariance
         E: EOF of U200 [k,d]
         n:  dims of state vector
         d:  dims of U200
         tau:LIM params  

        output
         skill: acc skill skill[d]      
    '''
    G   = expm(L.real*tau)
    tmp = numpy.dot(G,numpy.dot(c0,G.T))
    F = numpy.dot(E.T,numpy.dot(tmp[:n,:n],E))  #signal  

    skill = numpy.zeros(d)
    for i in range(d):
        S2 = F[i,i]/(C0[i,i]-F[i,i])
        if S2 <0:
            skill[i] = 1.
        else:
            skill[i] = S2/sqrt((S2+1)*S2)
    
    return skill 


n = 20   #EOFs of U200
lats = 10
latn = 80
lonw = 120
lone = 240

ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"

with open(ddir+'/LIM_params.pkl', 'rb') as file:
    data = pickle.load(file)

c0 = data['c0']   #lag-0 covariance matrix
L  = data['L']    #exp(L*tau)

#convert pcs skill into real skill. skill*EOFs
eofs = xarray.open_dataarray(ddir+'/eof_u200.nc')
eofs = eofs.loc[:n,lats:latn,lonw:lone]
eval = eofs.attrs.get('eval'); eval = numpy.sqrt(eval[:n])
eofs = eofs*eval[:, numpy.newaxis, numpy.newaxis]
eof  = eofs.to_numpy().reshape(eofs.shape[0],-1)

#covariance
print("covaraince begin.")
u = xarray.open_dataarray(ddir+"/LIM_obs_u200_v2.nc")
u = u.loc[:,:,lats:latn,lonw:lone]
lat = u.lat; lon = u.lon
X = u.sel(time=u.time.dt.month.isin([12, 1, 2, 3]))
X = u[:,0,:,:].to_numpy().reshape(u.shape[0],-1) #[t,d]
#C0 = numpy.cov(X,rowvar=False)   
C0 = numpy.matmul(X.T,X)/(X.shape[0]-1)
print(C0.shape) 
        
#with open(ddir+'/covar_obs_u200', 'wb') as f:
#    pickle.dump(C0,f)
#print("covaraince done.")

Tau = numpy.arange(0,47,1)
acc = numpy.zeros([Tau.shape[0],X.shape[1]])
for i,tau in enumerate(Tau):
    acc[i,:] = ACC(L,c0,C0,eof,n,X.shape[1],tau)

ac = acc.reshape(u.shape[1],u.shape[2],u.shape[3])
ds = xarray.Dataset({'acc': (('lead_time','lat','lon'), ac)},coords={'lead_time': Tau,'lat': lat,'lon': lon})
ds.to_netcdf(ddir+"/LIM_theskill.nc")