import sys
sys.path.insert(0, "/Users/helenarichie/GitHub/my_scripts/")
from hconfig import *
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

###################################################################################
sim = "0726"
cat = True
vlims_gas_small = (18.2, 19.6)
vlims_gas_large = (19.3, 22)
vlims_dust_small = (-10, -6)
vlims_dust_large = (-11, -3.5)
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
projdir_small = "/Users/helenarichie/Desktop/cloud_survival/0426"
projdir_large = "/Users/helenarichie/Desktop/cloud_survival/0726"
pngdir = os.path.join("/Users/helenarichie/Desktop/")
cmap_gas = sns.color_palette("mako", as_cmap=True)
cmap_dust = sns.color_palette("rocket", as_cmap=True)
###################################################################################

fig, ax = plt.subplots(nrows=2, ncols=3, gridspec_kw={"width_ratios":[1, 2, 2.4]}, figsize=(24.1,10))
if sim == "0426":
    fnums = fnums_small
    projdir = projdir_small
    vlims_gas = vlims_gas_small
    vlims_dust = vlims_dust_small
    spacing = spacing_small
if sim == "0726":
    fnums = fnums_large
    projdir = projdir_large
    vlims_gas = vlims_gas_large
    vlims_dust = vlims_dust_large
    spacing = spacing_large

for i, row in enumerate(ax):

    for j, col in enumerate(row):
        
        data = ReadHDF5(projdir, proj="xy", cat=cat, dust=True, fnum=fnums[j])
        head = data.head
        conserved = data.conserved
        dx = head["dx"][0]
        nx, ny, nz = head["dims"]

        xlen = [nz, nx-2*nz, nx-2*nz]
        # indices = [[[0, nz], [nz, 3*nz], [2*nz, 4*nz]], 
        #            [[0, nz], [0+64+64, 2*nz+64+64], [nz, 3*nz]]]
        if sim == "0426":
            indices = [[0, nz], [nz, 3*nz], [2*nz, 4*nz]]
        if sim == "0726":
            indices = [[0, nz], [0+64+64, 2*nz+64+64], [nz, 3*nz]]
        
        dx *= 1e3  # pc
        d_gas = conserved["density"]
        d_gas *= head["mass_unit"] / (head["length_unit"] ** 2)
        n_gas = d_gas[0]/(0.6*MP) # column density
        d_dust = conserved["dust_density"][0]
        d_dust *= head["mass_unit"] / (head["length_unit"] ** 2)
        d_dust[d_dust<=0] = 1e-40
        t = data.t_cgs() / yr_in_s  # yr
        if i == 0:
            im = ax[i][j].imshow(np.log10(n_gas[indices[j][0]:indices[j][1],:].T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=[0, xlen[j]*dx, 0, nz*dx], cmap=cmap_gas)
        if i == 1:
            im = ax[i][j].imshow(np.log10(d_dust[indices[j][0]:indices[j][1],:].T), origin="lower", vmin=vlims_dust[0], vmax=vlims_dust[1], extent=[0, xlen[j]*dx, 0, nz*dx], cmap=cmap_dust)
        
        if j == 0:
            if i == 1:
                if sim == "0426":
                    ax[i][j].hlines(0.12*dx*ny, 0.5*spacing, 1.5*spacing, color='white')
                    ax[i][j].text(1.5*spacing+2, 0.1*dx*ny, '20 pc', color='white', fontsize=fontsize)
                if sim == "0726":
                    ax[i][j].hlines(0.12*dx*ny, 0.5*spacing, 1.5*spacing, color='white')
                    ax[i][j].text(1.5*spacing+50, 0.1*dx*ny, '320 pc', color='white', fontsize=fontsize)

        if i == 0:
            ax[i][j].text(0.5*spacing, 0.85*dx*nz, f'{round(t[0]/1e6, 1)} Myr', color='white', fontsize=fontsize) 


        ax[i][j].set_xticks(np.arange(0, xlen[j]*dx, spacing))
        ax[i][j].set_yticks(np.arange(0, nz*dx, spacing))
        ax[i][j].tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=ticklength, width=tickwidth)

        if (i == 0) and (j == 2):
            cbar = fig.colorbar(im, ax=ax[i][j], location='right', pad=0.02, shrink=0.912, aspect=10)
            cbar.set_label(r'$\mathrm{log}_{10}(N_{H, gas})$ [$\mathrm{cm}^{-2}$]', rotation=270, labelpad=40, fontsize=fontsize)
            cbar.ax.tick_params(length=ticklength, width=tickwidth, color="white", labelsize=fontsize-5)
            cbar.set_ticks(np.linspace(vlims_gas[0], vlims_gas[1], 4).round(1))
        if (i == 1) and (j == 2):
            cbar = fig.colorbar(im, ax=ax[i][j], location='right', pad=0.02, shrink=0.912, aspect=10)
            cbar.set_label(r'$\mathrm{log}_{10}(\Sigma_{dust})$ [$\mathrm{g}\,\mathrm{cm}^{-2}$]', rotation=270, labelpad=40, fontsize=fontsize)
            cbar.ax.tick_params(length=ticklength, width=tickwidth, color="white", labelsize=fontsize-5)
            cbar.set_ticks(np.linspace(vlims_dust[0], vlims_dust[1], 4).round(1))
plt.tight_layout()      
plt.savefig(pngdir + f"cloud_dust_snapshots_{sim}.png", dpi=300)
plt.close()
