import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps

ddir = '/home/sunming/data5/cuixy/Subpre_NPJ/data/'
ac = xr.open_dataarray(ddir+'PrCA_ac.nc')
ac_sig = xr.open_dataarray(ddir+'PrCA_ac_sig.nc')

x = ac.step.values
y = ac[:6,:].values
z = ac_sig[:6,:].values

plt.figure(figsize=(12, 6))
colors = colormaps['tab20']
for i in range(y.shape[0]):
    plt.plot(x, y[i,:], 
             color=colors(i),  # 从colormap获取颜色
             linewidth=2,
             label=f'mode {i+1}')

for i in range(z.shape[0]):
    plt.plot(x, z[i,:], 
             color=colors(i),  # 从colormap获取颜色
             linewidth=1,
             linestyle='--')

plt.legend(ncol=4, bbox_to_anchor=(1.05, 1), loc='best')

plt.xlim((0,46))
plt.ylim((-0.1,1))
plt.xlabel('lead day', fontsize=12)
plt.ylabel('AC', fontsize=12)
#plt.axhline(0.5,color='#000000',linestyle='--',linewidth=1)

plt.tight_layout()
plt.show()
