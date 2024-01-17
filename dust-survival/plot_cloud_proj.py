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
vlims_gas_destroyed = (18.2, 19.6)
vlims_gas_disrupted = (19.3, 22)
vlims_gas_survived = (19.3, 22)
# vlims_gas_destroyed = vlims_gas_disrupted = (18.2, 21.8)
# fnums_destroyed = [0, 715, 1051]
fnums_destroyed = [0, 482, 865]
# fnums_disrupted = [0, 196, 304]
fnums_disrupted = [0, 143, 300]
# fnums_survived = [0, 102, 372]
fnums_survived = [0, 200, 425]
pad = 0.1
spacing_destroyed = 20  # pc, spacing of tick marks
spacing_disrupted = 320  # pc, spacing of tick marks
spacing_survived = 320  # pc, spacing of tick marks
fontsize = 28
labelpad = 12
tickwidth = 2
ticklength = 13
tmax_destroyed = 2.1e6
tmax_disrupted = 30.4e6
tmax_survived = 30.4e6
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': fontsize})
projdir_destroyed = "/Users/helenarichie/Desktop/dust-survival-snapshots/0426/proj_cl/"
projdir_disrupted = "/Users/helenarichie/Desktop/dust-survival-snapshots/0726/proj_cl/"
projdir_survived = "/Users/helenarichie/Desktop/dust-survival-snapshots/1107/proj_cl/"
pngdir = os.path.join("/Users/helenarichie/Desktop/")
cmap = sns.color_palette("mako", as_cmap=True)
###################################################################################

fig, ax = plt.subplots(nrows=3, ncols=3, gridspec_kw={"width_ratios":[1, 2, 2]}, figsize=(24.1,14.3), constrained_layout = True)
fnums = [fnums_destroyed, fnums_disrupted, fnums_survived]
projdir = [projdir_destroyed, projdir_disrupted, projdir_survived]
vlims = [vlims_gas_destroyed, vlims_gas_disrupted, vlims_gas_survived]
spacing = [spacing_destroyed, spacing_disrupted, spacing_survived]

for i, row in enumerate(ax):

    for j, col in enumerate(row):

        print(projdir[i])
        
        data = ReadHDF5(projdir[i], proj="xy", cat=cat, fnum=fnums[i][j])
        head = data.head
        conserved = data.conserved
        dx = head["dx"][0]
        nx, ny, nz = head["dims"]

        xlen = [nz, nx-2*nz, nx-2*nz]
        # [[0, nz], [nz, 3*nz], [2*nz, 4*nz]]
        # [[0, nz], [64+64, 2*nz+64+64], [nz, 3*nz]]
        indices = [[[0, nz], [int(0.65*nz), int(2.65*nz)], [int(1.5*nz), int(3.5*nz)]], 
                   [[0, nz], [int(0.15*nz), int(2.15*nz)], [int(0.75*nz), int(2.75*nz)]],
                   [[0, nz], [64-20, 2*nz+64-20], [int(1.6*nz), int(3.6*nz)]]]
        
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
            if i == 2:
                ax[i][j].hlines(0.12*dx*ny, 0.5*spacing[i], 1.5*spacing[i], color='white')
                ax[i][j].text(1.5*spacing[i]+50, 0.1*dx*ny, '320 pc', color='white', fontsize=fontsize)
        if j == 2:
            if i == 0:
                ax[i][j].text(6*spacing[i]+2, 0.85*dx*nz, 'destroyed', color='white', fontsize=fontsize)
            if i == 1:
                ax[i][j].text(7.5*spacing[i]+50, 0.85*dx*nz, 'disrupted', color='white', fontsize=fontsize)
            if i == 2:
                ax[i][j].text(7.5*spacing[i]+50, 0.85*dx*nz, 'survived', color='white', fontsize=fontsize)

        
        ax[i][j].text(0.5*spacing[i], 0.85*dx*nz, f'{round(t[0]/1e6, 1)} Myr', color='white', fontsize=fontsize) 

        ax[i][j].set_xticks(np.arange(0, xlen[j]*dx, spacing[i]))
        ax[i][j].set_yticks(np.arange(0, nz*dx, spacing[i]))
        ax[i][j].tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=ticklength, width=tickwidth)

        if (i == 0) and (j == 2):
            cbar = fig.colorbar(im, ax=ax[i][j], location='right', pad=0.02, shrink=0.99, aspect=10)
            cbar.set_label(r'$\mathrm{log}_{10}(N_{H, gas})$ [$\mathrm{cm}^{-2}$]', rotation=270, labelpad=40, fontsize=fontsize)
            cbar.ax.tick_params(length=ticklength, width=tickwidth, color="white", labelsize=fontsize-5)
            cbar.set_ticks(np.linspace(vlims[i][0], vlims[i][1], 4).round(1))
        if (i == 1) and (j == 2):
            cbar = fig.colorbar(im, ax=ax[i][j], location='right', pad=0.02, shrink=0.99, aspect=10)
            cbar.set_label(r'$\mathrm{log}_{10}(N_{H, gas})$ [$\mathrm{cm}^{-2}$]', rotation=270, labelpad=40, fontsize=fontsize)
            cbar.ax.tick_params(length=ticklength, width=tickwidth, color="white", labelsize=fontsize-5)
            cbar.set_ticks(np.linspace(vlims[i][0], vlims[i][1], 4).round(1))
        if (i == 2) and (j == 2):
            cbar = fig.colorbar(im, ax=ax[i][j], location='right', pad=0.02, shrink=0.99, aspect=10)
            cbar.set_label(r'$\mathrm{log}_{10}(N_{H, gas})$ [$\mathrm{cm}^{-2}$]', rotation=270, labelpad=40, fontsize=fontsize)
            cbar.ax.tick_params(length=ticklength, width=tickwidth, color="white", labelsize=fontsize-5)
            cbar.set_ticks(np.linspace(vlims[i][0], vlims[i][1], 4).round(1))
#plt.tight_layout()      
plt.savefig(pngdir + f"cloud_snapshots.png", dpi=300)
