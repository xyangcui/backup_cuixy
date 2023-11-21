import numpy  as np
import xarray as xr
import os
"""
This script is used to calculate stream function and velocity potential
by Libman method.

variables needed for calc:
convergence (time,lat,lon)
or
vorticity (time,lat,lon)
and
geopotential height (time,lat,lon)

variables for output:
streamfunction (time,lat,lon)
or
velocity potential (time,lat,lon)

defined variables in this script:
x: convergence that input for velocity potential
or 
vorticity that input for stream function
gph: geopotential height that input for boundary values of stream function.
y: stream function (phi)
or
velocity potential (chi)
time: time series
lat: latitude series
lon: longitude series
"""

def Libman(xl, yl, dxl, dyl, fil):

"""
this function use Libman method to solve 2-D.
Procedures:
Rn,i,j = (yn,i,j+1 - yn,i,j-1)/(deltax^2) + 
(yn,i+1,j - yn,i-1,j)/(deltax^2) -
2*(1/deltax^2 + 1/deltay^2)*yn,i,j + 
Xi,j

yn+1,i,j = yn,i,j + Rn,i,j / ( 2*(1/deltax^2 + 1/deltay^2) )

variables:
fil: critia
xl: vor or con
yl: inintialized sf or vp.
dxl: deltax
dyl: deltay
delt: to judge if procedure can stop.
R: residual.
"""
    delt = 0.0
    while delt > fil :
        delt = 0.0
        n1 = x1.shape[0]
        n2 = x1.shape[1]
        for i in range(1,n1-2,1):
            for j in range(1,n2-2,1):
                delta1 = 1.0 / (dxl[i,j]*dxl[i,j])
                delta2 = 1.0 / (dyl[i,j]*dyl[i,j])
                term1  = (yl[i,j+1] + yl[i,j-1]) * delta1
                term2  = (yl[i+1,j] + yl[i-1,j]) * delta2
                term3  = yl[i,j] * (2*delta1 + 2*delta2)
                R = term1 + term2 - term3 + xl[i,j]
                deltt  = R / (2*delta1 + 2*delta2)
                yl[i,j] = yl[i,j] + deltt
                deltt1 = np.abs(deltt)
                if deltt1 > delt :
                    delt = deltt
    return yl

### some useful parameters and definitions.
a  = 6371000.0              # earth's radius
pi = np.pi                  # pi
om = 2.0*pi/(24.0*60.0*60.)    # earth's angular momentum
fi = 100

fdirx   = os.getenv('dirx')        # file dir of vor or con.
fdiry   = os.getenv('diry')        # file dir of gph.
outdir  = os.getenv('outdir')      # file for output
xvar    = os.getenv('xvar')        # variable name for x.
yvar    = os.getenv('yvar')        # variable name for y.
ovar    = os.getenv('ovar')        # variable name for output.
### read in variables.
xf = xr.open_dataset(dirx)
yf = xr.open_dataset(diry)

x   = xf[xvar]
gph = yf[yvar]

x   = np.array(x)
gph = np.array(gph)

time1 = xf['time']
lon1  = xf['lon']
lat1  = xf['lat']

time  = np.array(time1)
lon   = np.array(lon1)
lat   = np.array(lat1)

### initiation for calc.
rad = pi/180.0 
f   = 2.0 * om * np.sin(rad*lat)     # coriori parameter.

## calc delta x and delta y.
deltax = np.abs(lon[0]-lon[1])
deltay = np.abs(lat[0]-lat[1])

dx = np.zeros(x[0,:,:], dtype = float)
dy = np.zeros(x[0,:,:], dtype = float)

dy = a * rad * deltay 
for i in range[lat.shape[0]]
    dx[i,:] = a * rad * deltax * np.cos(rad * lat[i]) 


y = np.zeros(x, dtype = float)

if ovar == 'sf':
    
    i = 0
    for j in range(lon.shape[0]):
        y[:,i,j] = gph[:,i,j] / f[i]
    i = lat.shape[0]
    for j in range(lon.shape[0]):
        y[:,i,j] = gph[:,i,j] / f[i]
    j = 0
    for i in range(lat.shape[0]):
        y[:,i,j] = gph[:,i,j] / f[i]
    j = lon.shape[0]
    for i in range(lat.shape[0]):
        y[:,i,j] = gph[:,i,j] / f[i]
    x = -1.0 * x


### formal calc by external function.
vout = np.zeros(x, dtype = float)
for k in range(time.shape[0]):
    vout[k,:,:] = Libman(x[k,:,:], y[k,:,:], dx, dy, fi)

ds = xr.Dataset({ovar: (('time','lat','lon'), vout)} coords={'time': time1,'lat': lat1,'lon': lon1})

ds.to_netcdf(outdir)

