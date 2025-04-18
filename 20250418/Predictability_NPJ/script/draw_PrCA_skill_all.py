import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps

ddir = '/home/sunming/data5/cuixy/Subpre_NPJ/data/'
ac = xr.open_dataarray(ddir+'PrCA_ac.nc')
ac_sig = xr.open_dataarray(ddir+'PrCA_ac_sig.nc')

x = ac.step.values
y = ac.values
z = ac_sig.values

fig, axs = plt.subplots(5, 4, figsize=(20, 25))
colors = colormaps['tab20']

for i in range(5):     
    for j in range(4):
        n  = i*4 + j
        ax = axs[i, j]  
        ax.plot(x, y[n,:], 
             color=colors(n),  
             linewidth=2,
             label=f'mode {n+1}')

        ax.plot(x, z[n,:], 
             color=colors(n),  
             linewidth=1,
             linestyle='--')

        ax.set_title(f'mode {n + 1}', fontsize=10)
        ax.set_xlim(0, 46)
        ax.set_ylim(-0.1, 1)

        ax.set_xlabel('lead day', fontsize=9)

        #if i == 4:
        #    ax.set_xlabel('lead day', fontsize=9)
        #else:
        #    ax.set_xticklabels([])
        
        if j == 0:
            ax.set_ylabel('AC', fontsize=9)
        else:
            ax.set_yticklabels([])

plt.tight_layout()
fig.subplots_adjust(top=0.95, hspace=0.3, wspace=0.3)
plt.show()
