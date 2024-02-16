import sys
sys.path.insert(0, "/ix/eschneider/helena/code/my_scripts/")
from hconfig import *
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3
plt.rcParams.update({'font.family': 'Helvetica'})

##############################################
date = "frontier/2024-02-09"
datestr = "0209"
# vlims_gas = (18.55, 20.5)    # destroyed
# vlims_dust = (-10.5, -5.2)   # destroyed
vlims_gas = (19, 22)    # disrupted
vlims_dust = (-10.5, -5.2)   # disrupted
# fnums = [0, 400, 800]     # survived
# spacing = 320             # survived
# fnums = [0, 300, 500]     # destroyed
# spacing = 20              # destroyed
fnums = [0, 240, 480]     # disrupted
# spacing = 20              # destroyed
spacing, unit = 320*1e-3, "kpc"
##############################################

##############################################
cat = True
pad = 0.1
fontsize = 28
labelpad = 12
tickwidth = 2
ticklength = 13
##############################################

##############################################
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
projdir = os.path.join(basedir, "hdf5/proj/")
pngdir = os.path.join(basedir, "png/")
##############################################

cmap_gas = sns.color_palette("mako", as_cmap=True)
cmap_dust = sns.color_palette("rocket", as_cmap=True)
plt.rcParams.update({'font.size': fontsize})

fig, ax = plt.subplots(nrows=2, ncols=3, gridspec_kw={"width_ratios":[1, 2, 2]}, figsize=(24.1,9.4), constrained_layout=True)

for i, row in enumerate(ax):
    for j, col in enumerate(row):
        data = ReadHDF5(projdir, proj="xy", cat=cat, dust=True, fnum=fnums[j])
        head = data.head
        conserved = data.conserved
        dx = head["dx"][0]
        nx, ny, nz = head["dims"]
        
        nz, ny = 400, 400

        # xlen = [nz, nx-2*nz, nx-2*nz]
        # indices = [[0, nz], [300, 2*nz+300], [int(1.2*nz), int(3.2*nz)]]  # survived
        # xlen = [400, 800, 800]
        # indices = [[0, 400], [250, 1050], [650, 1450]]  # destroyed
        # indices_y = [40, 440]
        # xlen = [nz, 2*nz, 2*nz]
        xlen = [nz, nx-2*nz, nx-2*nz]
        indices = [[0, nz], [300, 2*nz+300], [int(1.2*nz), int(3.2*nz)]]

        if unit == "pc":
            dx *= 1e3  # pc
        d_gas = conserved["density"]
        d_gas *= head["mass_unit"] / (head["length_unit"] ** 2)
        n_gas = d_gas[0]/(0.6*MP) # column density
        d_dust = conserved["dust_density"][0]
        d_dust[d_dust==0] = 1e-40
        d_dust *= head["mass_unit"] / (head["length_unit"] ** 2)
        t = data.t_cgs() / yr_in_s  # yr

        if i == 0:
            im = ax[i][j].imshow(np.log10(n_gas[indices[j][0]:indices[j][1],:].T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=[0, xlen[j]*dx, 0, nz*dx], cmap=cmap_gas)
            if datestr == "0205":
                im = ax[i][j].imshow(np.log10(n_gas[indices[j][0]:indices[j][1],indices_y[0]:indices_y[1]].T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=[0, xlen[j]*dx, 0, 400*dx], cmap=cmap_gas)
        if i == 1:
            im = ax[i][j].imshow(np.log10(d_dust[indices[j][0]:indices[j][1],:].T), origin="lower", vmin=vlims_dust[0], vmax=vlims_dust[1], extent=[0, xlen[j]*dx, 0, nz*dx], cmap=cmap_dust)
            if datestr == "0205":
                im = ax[i][j].imshow(np.log10(d_dust[indices[j][0]:indices[j][1],indices_y[0]:indices_y[1]].T), origin="lower", vmin=vlims_dust[0], vmax=vlims_dust[1], extent=[0, xlen[j]*dx, 0, 400*dx], cmap=cmap_dust)

        if j == 0:
            if i == 1:
                if datestr == "0205":
                    ax[i][j].hlines(0.88*dx*nz, 0.5*spacing, 1.5*spacing, color='white')
                    ax[i][j].text(1.5*spacing+2, 0.85*dx*nz, '20 pc', color='white', fontsize=fontsize)
                    ax[i][j].text(0.65*spacing, 0.125*dx*ny, 'dust', color='white', fontsize=fontsize)
                else:
                    ax[i][j].hlines(0.88*dx*nz, 0.5*spacing, 1.5*spacing, color='white')
                    ax[i][j].text(1.5*spacing+50, 0.85*dx*nz, '320 pc', color='white', fontsize=fontsize)
                    ax[i][j].text(0.59*spacing, 0.125*dx*ny, 'dust', color='white', fontsize=fontsize)

        
        if i == 0:
            ax[i][j].text(0.5*spacing, 0.85*dx*nz, f'{round(t[0]/1e6, 1)} Myr', color='white', fontsize=fontsize)
            if j == 0:
                if datestr == "0205":
                    ax[i][j].text(0.6*spacing, 0.125*dx*ny, 'gas', color='white', fontsize=fontsize)
                else:
                    ax[i][j].text(0.6*spacing, 0.125*dx*ny, 'gas', color='white', fontsize=fontsize)

        ax[i][j].set_xticks(np.arange(0, xlen[j]*dx, spacing))
        ax[i][j].set_yticks(np.arange(0, nz*dx, spacing))
        ax[i][j].tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=ticklength, width=tickwidth)

        if (i == 0) and (j == 2):
            cbar = fig.colorbar(im, ax=ax[i][j], location='right', pad=0.02, aspect=17)
            cbar.set_label(r'$\mathrm{log}_{10}(N_{H, gas})$ [$\mathrm{cm}^{-2}$]', rotation=270, labelpad=40, fontsize=fontsize)
            cbar.ax.tick_params(length=ticklength, width=tickwidth, color="white", labelsize=fontsize-5)
            cbar.set_ticks(np.linspace(vlims_gas[0], vlims_gas[1], 4).round(1))
        if (i == 1) and (j == 2):
            cbar = fig.colorbar(im, ax=ax[i][j], location='right', pad=0.02, aspect=17)
            cbar.set_label(r'$\mathrm{log}_{10}(\Sigma_{dust})$ [$\mathrm{g}\,\mathrm{cm}^{-2}$]', rotation=270, labelpad=40, fontsize=fontsize)
            cbar.ax.tick_params(length=ticklength, width=tickwidth, color="white", labelsize=fontsize-5)
            cbar.set_ticks(np.linspace(vlims_dust[0], vlims_dust[1], 4).round(1))

plt.savefig(pngdir + f"cloud_dust_snapshots_{datestr}.png", dpi=300)
plt.close()
