import numpy as np
import xarray as xr
import statsmodels.api as sm
from math import floor
from scipy import stats

# get effeicent freedom.
def compute_equation(N, X, Y):
    Neff = 1 / N  
    rho_xx = sm.tsa.acf(X, nlags=N-1)  
    rho_yy = sm.tsa.acf(Y, nlags=N-1)  
    
    for k in range(1, N+1):
        Neff += 2/N  *  ((N - k) / N) * rho_xx[k-1] * rho_yy[k-1]
        
    return floor(1/Neff)

def critical_correlation(alpha=0.05, df=100):
    t_critical = stats.t.ppf(1 - alpha, df)

    return t_critical / np.sqrt(t_critical**2 + df)

ddir = '/home/sunming/data5/cuixy/Subpre_NPJ/data/'

ac = xr.open_dataarray(ddir+'PrCA_ac.nc')
v_obs = xr.open_dataarray(ddir+'PrCA_obs_variates.nc')
v_mod = xr.open_dataarray(ddir+'PrCA_variates.nc')
v_mod = v_mod.mean('number')

v_obs = v_obs.transpose('mode','lead_time','time')
v_mod = v_mod.transpose('mode','lead_time','time')

# get freedom.
N  = v_obs.shape[2]

Nef = np.zeros_like(v_obs[:,:,0])
for i in range(Nef.shape[0]):
    Nef[i,:] = [compute_equation(N,v_obs[i,j,:],v_mod[i,j,:]) for j in range(Nef.shape[1])]

# get r critical.
r_cri = np.zeros_like(v_obs[:,:,0])
for i in range(r_cri.shape[0]):
    r_cri[i,:] = [critical_correlation(alpha=0.05, df=Nef[i,j]) for j in range(r_cri.shape[1])]

print(r_cri)

mode = v_mod.mode
lead = np.arange(0,46+1,1)

r = xr.DataArray(np.zeros(( np.size(mode), np.size(lead) )), 
    dims=["mode", "step"],coords={"mode": mode, "step": lead}).rename('ac')
r.values = r_cri
r.to_netcdf(ddir+'PrCA_ac_sig.nc')