import sys
sys.path.insert(0, "/Users/helenarichie/GitHub/my_scripts/")
from hconfig import *
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

###################################################################################
sim = "1107"
cat = True
vlims_gas_destroyed = (18.2, 19.6)
vlims_gas_disrupted = (19.3, 22)
vlims_gas_survived = (19.3, 22)
vlims_dust_destroyed = (-9.5, -6.2)
vlims_dust_survived = (-8, -4.2)
vlims_dust_disrupted = (-10.8, -4.1)
# vlims_gas_destroyed = vlims_gas_disrupted = (18.2, 21.8)
fnums_destroyed = [0, 482, 865]
fnums_disrupted = [0, 143, 300]
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
tmax_survived = 58.8e6
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': fontsize})
projdir_destroyed = "/Users/helenarichie/Desktop/dust-survival-snapshots/0426/proj_cl/"
projdir_disrupted = "/Users/helenarichie/Desktop/dust-survival-snapshots/0726/proj_cl/"
projdir_survived = "/Users/helenarichie/Desktop/dust-survival-snapshots/1107/proj_cl/"
pngdir = os.path.join("/Users/helenarichie/Desktop/")
cmap_gas = sns.color_palette("mako", as_cmap=True)
cmap_dust = sns.color_palette("rocket", as_cmap=True)
###################################################################################

fig, ax = plt.subplots(nrows=2, ncols=3, gridspec_kw={"width_ratios":[1, 2, 2]}, figsize=(24.1,9.4), constrained_layout = True)
if sim == "0426":
    fnums = fnums_destroyed
    projdir = projdir_destroyed
    vlims_gas = vlims_gas_destroyed
    vlims_dust = vlims_dust_destroyed
    spacing = spacing_destroyed
if sim == "0726":
    fnums = fnums_disrupted
    projdir = projdir_disrupted
    vlims_gas = vlims_gas_disrupted
    vlims_dust = vlims_dust_disrupted
    spacing = spacing_disrupted
if sim == "1107":
    fnums = fnums_survived
    projdir = projdir_survived
    vlims_gas = vlims_gas_survived
    vlims_dust = vlims_dust_survived
    spacing = spacing_survived

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
            indices = [[0, nz], [int(0.75*nz), int(2.75*nz)], [int(1.65*nz), int(3.65*nz)]]
        if sim == "0726":
            indices = [[0, nz], [int(0.35*nz), int(2.35*nz)], [int(0.95*nz), int(2.95*nz)]]
        if sim == "1107":
            indices = [[0, nz], [70, 2*nz+70], [int(1.6*nz), int(3.6*nz)]]
        
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
                    ax[i][j].hlines(0.88*dx*nz, 0.5*spacing, 1.5*spacing, color='white')
                    ax[i][j].text(1.5*spacing+2, 0.85*dx*nz, '20 pc', color='white', fontsize=fontsize)
                    ax[i][j].text(0.65*spacing, 0.125*dx*ny, 'dust', color='white', fontsize=fontsize)
                if sim == "0726":
                    ax[i][j].hlines(0.88*dx*nz, 0.5*spacing, 1.5*spacing, color='white')
                    ax[i][j].text(1.5*spacing+50, 0.85*dx*nz, '320 pc', color='white', fontsize=fontsize)
                    ax[i][j].text(0.59*spacing, 0.125*dx*ny, 'dust', color='white', fontsize=fontsize)
                if sim == "1107":
                    ax[i][j].hlines(0.88*dx*nz, 0.5*spacing, 1.5*spacing, color='white')
                    ax[i][j].text(1.5*spacing+50, 0.85*dx*nz, '320 pc', color='white', fontsize=fontsize)
                    ax[i][j].text(0.59*spacing, 0.125*dx*ny, 'dust', color='white', fontsize=fontsize)

        if i == 0:
            ax[i][j].text(0.5*spacing, 0.85*dx*nz, f'{round(t[0]/1e6, 1)} Myr', color='white', fontsize=fontsize)
            if j == 0:
                if sim == "0426":
                    ax[i][j].text(0.65*spacing, 0.125*dx*ny, 'gas', color='white', fontsize=fontsize)
                else:
                    ax[i][j].text(0.6*spacing, 0.125*dx*ny, 'gas', color='white', fontsize=fontsize)


        ax[i][j].set_xticks(np.arange(0, xlen[j]*dx, spacing))
        ax[i][j].set_yticks(np.arange(0, nz*dx, spacing))
        ax[i][j].tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=ticklength, width=tickwidth)

        if (i == 0) and (j == 2):
            cbar = fig.colorbar(im, ax=ax[i][j], location='right', pad=0.02, shrink=0.99, aspect=10)
            cbar.set_label(r'$\mathrm{log}_{10}(N_{H, gas})$ [$\mathrm{cm}^{-2}$]', rotation=270, labelpad=40, fontsize=fontsize)
            cbar.ax.tick_params(length=ticklength, width=tickwidth, color="white", labelsize=fontsize-5)
            cbar.set_ticks(np.linspace(vlims_gas[0], vlims_gas[1], 4).round(1))
        if (i == 1) and (j == 2):
            cbar = fig.colorbar(im, ax=ax[i][j], location='right', pad=0.02, shrink=0.99, aspect=10)
            cbar.set_label(r'$\mathrm{log}_{10}(\Sigma_{dust})$ [$\mathrm{g}\,\mathrm{cm}^{-2}$]', rotation=270, labelpad=40, fontsize=fontsize)
            cbar.ax.tick_params(length=ticklength, width=tickwidth, color="white", labelsize=fontsize-5)
            cbar.set_ticks(np.linspace(vlims_dust[0], vlims_dust[1], 4).round(1))

#plt.tight_layout()      
plt.savefig(pngdir + f"cloud_dust_snapshots_{sim}.png", dpi=300)
plt.close()
