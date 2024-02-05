import h5py
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

###############################
date = "2024-01-29"
ns = 0
ne = 61
DE = True
SCALAR = True
###############################

###############################
crc = True
frontier = False
###############################

###############################
testing = True
cloud_wind = False
###############################

########### plotting #############
vlims = True
vlims_gas = (-23.0459 , -22.301) # g/cm^3
spacing, unit = 640*1e-3, "kpc" # spacing of tick marks in units and sets axes labels and units of dx (kpc or pc)
# spacing, unit = 40, "pc"
fontsize = 20
#################################

plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': fontsize})

if crc:
  if cloud_wind:
    basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
  if testing:
    basedir = f"/ix/eschneider/helena/data/testing/{date}/"
  datadir = os.path.join(basedir, "hdf5/edges/")
  pngdir = os.path.join(basedir, "png/edges/")
if frontier:
  basedir = f"/lustre/orion/ast181/scratch/helenarichie/{date}/"
  datadir = os.path.join(basedir, "hdf5/edges/")
  pngdir = os.path.join(basedir, "png/edges/")

n_faces = 6
n_start = 2
key_nums = range(0, n_faces+1)

# loop over the output times
for n in range(ns, ne+1):
    f = h5py.File(os.path.join(datadir, str(n)+"_edges.h5"))
    keys = list(f.keys())[n_start*n_faces:(n_start+1)*n_faces]
    head = f.attrs
    nx, ny, nz = head["dims"]
    dx, dy, dz = head["dx"]
    t = head["t"][0]

    extx = [ny*dx, nz*dx, nz*dx, ny*dx, nz*dx, nz*dx]
    exty = [nx*dx, nx*dx, ny*dx, nx*dx, nx*dx, ny*dx]

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(30, 20), gridspec_kw={"width_ratios": [nx, nx, nz]})

    ax = ax.flatten()

    for i, key in enumerate(keys):
        data = f[key][()]

        if vlims:
            im = ax[i].imshow(np.log10(data.T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=[0, extx[i], 0, exty[i]])
        else:
            im = ax[i].imshow(np.log10(data.T), origin="lower", extent=[0, extx[i], 0, exty[i]])

        ylabel = r'$\mathrm{log}_{10}(\rho_{gas})$ [$\mathrm{g}\mathrm{cm}^{-3}$]'
        divider = make_axes_locatable(ax[i])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(im, ax=ax[i], cax=cax)
        cbar.set_label(ylabel, fontsize=fontsize)
        ax[i].set_xticks(np.arange(0, nx*dx, spacing))
        ax[i].set_yticks(np.arange(0, ny*dx, spacing))
        ax[i].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=8)
        ax[i].set_title(key, fontsize=fontsize)
        ax[i].set_xlabel(r"$x~$[{}]".format(unit), fontsize=fontsize)
        ax[i].set_ylabel(r"$y~$[{}]".format(unit), fontsize=fontsize)
        ax[i].text(spacing, 0.1*dx*ny, f'{round(t/1e3, 2)} Myr', color='white', fontsize=fontsize)
    
    plt.tight_layout()
    plt.savefig(os.path.join(pngdir, str(n) + "_edges.png"), dpi=300)
    plt.close()

    print(f"Saving figure {n} of {ne}.\n")


