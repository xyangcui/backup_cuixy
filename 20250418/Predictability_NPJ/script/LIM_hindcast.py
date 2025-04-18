import numpy, xarray, pickle
from numpy import linalg as LA
from LIM_utils import LIM, Nyquist_check, Q_test, forecast

tau0 = 1
ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"

x = xarray.open_dataarray(ddir+"/state_vectors.nc")
x = x.to_numpy()
ndim = x.shape[0]
ntim = x.shape[1]
print(ntim)
pcs = numpy.full([ndim,ntim,47],numpy.nan)

#ten-fold cross validation. 40years 10 ensemble 4years.
nyears= 40
ndays = 121
#create ensemble and verification.
for i in range(10):
	nsrt = 0 + i*ndays*4
	nend = ndays*4 + i*ndays*4
	xfcst  = x[:,nsrt:nend]
	xtrain = numpy.delete(x,numpy.s_[nsrt:nend],axis=1)
	print("")
	print("year= "+str(i*4+1982)+"-"+str(i*4+1985))
	output = LIM(xtrain,tau0)
	L = output['L']; Q = output['Q']
	Q_test(Q,Q_plot='no')
	#Nyquist_check(x,tau0)
	print("")
	#fcst.
	for tau in range(1,47):
		pcs_fcst = forecast(output['g'],output['normU'],output['v'],tau0,xfcst,tau)
		pcs[:,nsrt:nend,tau] = pcs_fcst
	pcs[:,nsrt:nend,0] = xfcst

#croods
var = range(ndim)
xt  = xarray.open_dataset(ddir+"/ecmwf/ecmwf_pf_anom_u200.nc")
lead_time = xt['u'].lead_time

xt = xarray.open_dataarray(ddir+"/ERA5_daily_anom_u200_1982-2022.nc")
xt = xt.loc['1982-12-01':'2022-03-31',:,:]
xt = xt.sel(time=xt.time.dt.month.isin([12, 1, 2, 3]))
time = xt.time

pcs = xarray.DataArray(pcs,dims=["var", "time", "lead_time"],coords={"var": var, "time": time, "lead_time": lead_time})
ds  = xarray.Dataset({'x': pcs}).to_netcdf(ddir+"/LIM_hindcast.nc")