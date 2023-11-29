import sys
sys.path.insert(0, "/Users/helenarichie/GitHub/my_scripts/")
from hconfig import *
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

###################################################################################
cat = True
vlims_gas_small = (18.2, 19.6)
vlims_gas_large = (19.3, 22)
# vlims_gas_small = vlims_gas_large = (18.2, 21.8)
fnums_small = [0, 715, 1051]
fnums_large = [0, 196, 304]
pad = 0.1
spacing_small = 20  # pc, spacing of tick marks
spacing_large = 320  # pc, spacing of tick marks
fontsize = 28
labelpad = 12
tickwidth = 1.75
ticklength = 12
tmax_small = 2.1e6
tmax_large = 30.4e6
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': fontsize})
projdir_small = csvdir_small = "/Users/helenarichie/Desktop/cloud_survival/0426"
projdir_large = csvdir_large = "/Users/helenarichie/Desktop/cloud_survival/0726"
pngdir = os.path.join("/Users/helenarichie/Desktop/")
cmap = sns.color_palette("mako", as_cmap=True)
###################################################################################

fig, ax = plt.subplots(nrows=2, ncols=3, gridspec_kw={"width_ratios":[1, 2, 2.4]}, figsize=(24.1,10))
fnums = [fnums_small, fnums_large]
projdir = [projdir_small, projdir_large]
csvdir = [csvdir_small, csvdir_large]
vlims = [vlims_gas_small, vlims_gas_large]
spacing = [spacing_small, spacing_large]

for i, row in enumerate(ax):

    for j, col in enumerate(row):
        
        data = ReadHDF5(projdir[i], proj="xy", cat=cat, fnum=fnums[i][j])
        head = data.head
        conserved = data.conserved
        dx = head["dx"][0]
        nx, ny, nz = head["dims"]

        xlen = [nz, nx-2*nz, nx-2*nz]
        indices = [[[0, nz], [nz, 3*nz], [2*nz, 4*nz]], 
                   [[0, nz], [0+64+64, 2*nz+64+64], [nz, 3*nz]]]
        
        dx *= 1e3  # pc
        d_gas = conserved["density"]
        d_gas *= head["mass_unit"] / (head["length_unit"] ** 2)
        n_gas = d_gas[0]/(0.6*MP) # column density
        t = data.t_cgs() / yr_in_s  # yr
        im = ax[i][j].imshow(np.log10(n_gas[indices[i][j][0]:indices[i][j][1],:].T), origin="lower", vmin=vlims[i][0], vmax=vlims[i][1], extent=[0, xlen[j]*dx, 0, nz*dx], cmap=cmap)
        
        if j == 0:
            if i == 0:
                ax[i][j].hlines(0.12*dx*ny, 0.5*spacing[i], 1.5*spacing[i], color='white')
                ax[i][j].text(1.5*spacing[i]+2, 0.1*dx*ny, '20 pc', color='white', fontsize=fontsize)
            if i == 1:
                ax[i][j].hlines(0.12*dx*ny, 0.5*spacing[i], 1.5*spacing[i], color='white')
                ax[i][j].text(1.5*spacing[i]+50, 0.1*dx*ny, '320 pc', color='white', fontsize=fontsize)

        
        ax[i][j].text(0.5*spacing[i], 0.85*dx*nz, f'{round(t[0]/1e6, 1)} Myr', color='white', fontsize=fontsize) 

        ax[i][j].set_xticks(np.arange(0, xlen[j]*dx, spacing[i]))
        ax[i][j].set_yticks(np.arange(0, nz*dx, spacing[i]))
        ax[i][j].tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=ticklength, width=tickwidth)

        if (i == 0) and (j == 2):
            cbar = fig.colorbar(im, ax=ax[i][j], location='right', pad=0.02, shrink=0.912, aspect=10)
            cbar.set_label(r'$\mathrm{log}_{10}(N_{H, gas})$ [$\mathrm{cm}^{-2}$]', rotation=270, labelpad=40, fontsize=fontsize)
            cbar.ax.tick_params(length=ticklength, width=tickwidth, color="white", labelsize=fontsize-5)
            cbar.set_ticks(np.linspace(vlims[i][0], vlims[i][1], 4).round(1))
        if (i == 1) and (j == 2):
            cbar = fig.colorbar(im, ax=ax[i][j], location='right', pad=0.02, shrink=0.912, aspect=10)
            cbar.set_label(r'$\mathrm{log}_{10}(N_{H, gas})$ [$\mathrm{cm}^{-2}$]', rotation=270, labelpad=40, fontsize=fontsize)
            cbar.ax.tick_params(length=ticklength, width=tickwidth, color="white", labelsize=fontsize-5)
            cbar.set_ticks(np.linspace(vlims[i][0], vlims[i][1], 4).round(1))
plt.tight_layout()      
plt.savefig(pngdir + f"cloud_snapshots.png", dpi=300)
plt.close()
