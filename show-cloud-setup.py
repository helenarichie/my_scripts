import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'font.family': 'sans serif'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.fontset':'stixsans'})
plt.rcParams.update({'axes.linewidth': 1.5})
plt.rcParams.update({'xtick.direction':'in'})
plt.rcParams.update({'xtick.major.size': 5})
plt.rcParams.update({'xtick.major.width': 1.25 })
plt.rcParams.update({'xtick.minor.size': 2.5})
plt.rcParams.update({'xtick.minor.width': 1.25 })
plt.rcParams.update({'ytick.direction':'in'})
plt.rcParams.update({'ytick.major.size': 5})
plt.rcParams.update({'ytick.major.width': 1.25 })
plt.rcParams.update({'ytick.minor.size': 2.5})
plt.rcParams.update({'ytick.minor.width': 1.25 })

# import warnings
# warnings.filterwarnings("ignore")

T_colors = ["lightskyblue", "violet", "forestgreen", "darkviolet", "thistle", "cadetblue", 
            "palegoldenrod", "darksalmon", "indigo"]

import h5py

import glob
import os
import sys
sys.path.insert(0, "/ihome/hrichie/her45/GitHub/analysis_scripts")

from cholla_py_utils import ChollaPyUtils
import cholla_plotter
from cholla_plotter import plotter

cholla_helper = ChollaPyUtils()

x_slice = 20
y_slice = 64
z_slice = 64
slices = (x_slice, y_slice, z_slice)

# directory with slices
datadir = "/ix/eschneider/helena/data/2022-11-02/hdf5/slices/"
outdir = datadir

# read in slice data
files = sorted(glob.glob(os.path.join(datadir, "0_slice.h5.0")))

n_out = len(files)

nx, ny, nz = None, None, None
dx, dy, dz = None, None, None

t_arr = []

d_gas, d_dust = [], []

d_gas_ch = []

vx = []
T, P = [], []

n = []

#i_out = np.arange(0, n_out)
i_out = [8]

units = set()

E_ch = []

print("Loading in slice data...\n")
for j, i in enumerate(i_out):
    f = h5py.File(datadir + str(i) + "_slice.h5.0", "r")
    head = f.attrs
    a_group_key = list(f.keys())
    # number of cells in each direction
    nx = head['dims'][0]
    ny = head['dims'][1]
    nz = head['dims'][2]

    dx = head['dx'][0]
    
    # time in years
    t_arr.append(np.array(cholla_helper.t_cgs(head["t"][0])/cholla_helper._YR_IN_S))

    # mass density in code units
    d_gas_ch = np.array(f["d_xy"])
    
    d_gas.append(np.array(cholla_helper.d_cgs(np.array(f["d_xy"]))))
    d_dust.append(np.array(cholla_helper.d_cgs(np.array(f["scalar_xy"]))))
    
    # x, y, and z momentum in code units
    mx_ch = np.array(f["mx_xy"])
    my_ch = np.array(f["my_xy"])
    mz_ch = np.array(f["mz_xy"])
    
    # x, y, and z velocity in code units
    vx_ch = mx_ch / d_gas_ch
    vy_ch = my_ch / d_gas_ch
    vz_ch = mz_ch / d_gas_ch

    # get x velocity in cgs units
    vx.append(vx_ch*cholla_helper._VELOCITY_UNIT)
    
    # get temperature in cgs units
    E_ch = np.array(f["E_xy"])
    T.append(np.array(cholla_helper.calc_T(E_ch, vx_ch, vy_ch, vz_ch, d_gas_ch)))
    
    # pressure in cgs units
    P.append(np.array(cholla_helper.calc_P_cgs(T[j], cholla_helper.calc_n_cgs(d_gas[j]))))
print("Data loaded.\n")

d_dust = np.array(d_dust)
wh_zero = np.where(d_dust<=0)
d_dust[wh_zero] = 1e-40
vlims_du = [np.log10(np.amin(d_dust.flatten())), np.log10(np.amax(d_dust.flatten()))]

d_gas = np.array(d_gas)
T = np.array(T)

vx = np.array(vx)

cloud = []
cloud_T = []

background = np.where(np.logical_and(T > 1e6, d_gas < 1e-24))
d_gas[background] = np.nan
T[background] = np.nan

cloud = np.array(cloud)
cloud_T = np.array(cloud_T)

for i, d in enumerate(d_gas):

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(20,7))
    
    # xy gas density projection
    im = axs[0][0].imshow(np.log10(d_gas[i].T), origin="lower", extent=[0, nx*dx, 0, nz*dx])
    ylabel = r'$\mathrm{log}_{10}(\rho)$ [$\mathrm{g}\mathrm{cm}^{-3}$]'
    divider = make_axes_locatable(axs[0][0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[0][0], cax=cax)
    cbar.set_label(ylabel)
    axs[0][0].set_xticks(nx*dx*np.arange(0, 1.25, 0.25))
    axs[0][0].set_yticks(nz*dx*np.arange(0, 1.25, 0.5))
    axs[0][0].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1)
    axs[0][0].set_title(r"Gas Density Slice")
    axs[0][0].set_xlabel(r"$x~$[kpc]")
    axs[0][0].set_ylabel(r"$z~$[kpc]")
    axs[0][0].text(10, 15, '$Time={:.1e}~yr$'.format(t_arr[i]), bbox={'facecolor': 'white', 'pad': 8}, fontsize=10, alpha=0.6)
    

    # xy density-weighted temperature projection
    im = axs[1][0].imshow(np.log10(d_dust[i].T), origin="lower", cmap="plasma", vmin=-28, vmax=-22, extent=[0, nx*dx, 0, nz*dx])
    ylabel = r'$\mathrm{log}_{10}(\rho)$ [$\mathrm{g}\mathrm{cm}^{-3}$]'
    divider = make_axes_locatable(axs[1][0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[1][0], cax=cax)
    cbar.set_label(ylabel)
    axs[1][0].set_xticks(nx*dx*np.arange(0, 1.25, 0.25))
    axs[1][0].set_yticks(nz*dx*np.arange(0, 1.25, 0.5))
    axs[1][0].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1)
    axs[1][0].set_title(r"Dust Density Slice")
    axs[1][0].set_xlabel(r"$x~$[kpc]")
    axs[1][0].set_ylabel(r"$z~$[kpc]")
    axs[1][0].text(10, 15, '$Time={:.1e}~yr$'.format(t_arr[i]), bbox={'facecolor': 'white', 'pad': 8}, fontsize=10, alpha=0.6)
    
    # xy temperature slice
    # im = axs[0][1].imshow(np.log10(T[i].T), origin="lower", cmap="inferno", vmin=2, vmax=vlims_T[1])

    im = axs[0][1].imshow(np.log10(T[i].T), origin="lower", cmap="inferno", extent=[0, nx*dx, 0, nz*dx])
    
    if np.isnan(np.log10(T[i].T).any()):
        print("there's a nan")

    ylabel = r'$\mathrm{log}_{10}(T) [\mathrm{K}]$'
    divider = make_axes_locatable(axs[0][1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[0][1], cax=cax)
    cbar.set_label(ylabel)
    axs[0][1].set_xticks(nx*dx*np.arange(0, 1.25, 0.25))
    axs[0][1].set_yticks(nz*dx*np.arange(0, 1.25, 0.5))
    axs[0][1].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1)
    axs[0][1].set_title(r"Temperature Slice")
    axs[0][1].set_xlabel(r"$x~$[kpc]")
    axs[0][1].set_ylabel(r"$z~$[kpc]")
    axs[0][1].text(10, 15, '$Time={:.1e}~yr$'.format(t_arr[i]), bbox={'facecolor': 'white', 'pad': 8}, fontsize=10, alpha=0.6)
    
    """
    # xz dust to gas ratio slice
    dust_to_gas = d_dust[i] / d_gas[i]
    im = axs[1][1].imshow(dust_to_gas.T, origin="lower", vmin=-1e-3, vmax=0.011, extent=[0, nx*dx, 0, nz*dx])
    ylabel = r'Dust to Gas Ratio'

    divider = make_axes_locatable(axs[1][1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[1][1], cax=cax)
    cbar.set_label(ylabel)
    axs[1][1].set_xticks(nx*dx*np.arange(0, 1.25, 0.25))
    axs[1][1].set_yticks(nz*dx*np.arange(0, 1.25, 0.5))
    axs[1][1].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1)
    axs[1][1].set_title(r"Dust to Gas Ratio Slice")
    axs[1][1].set_xlabel(r"$x~$[kpc]")
    axs[1][1].set_ylabel(r"$z~$[kpc]")
    axs[1][1].text(10, 15, '$Time={:.1e}~yr$'.format(t_arr[i]), bbox={'facecolor': 'white', 'pad': 8}, fontsize=10, alpha=0.6)
    """

    # xz velocity slice
    # im = axs[1][1].imshow(vx[i].T, origin="lower", vmin=-1e-3, vmax=0.011, extent=[0, nx*dx, 0, nz*dx])
    im = axs[1][1].imshow(vx[i].T*1e-5, origin="lower", extent=[0, nx*dx, 0, nz*dx])
    ylabel = r'x-velocity [km/s]'

    divider = make_axes_locatable(axs[1][1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[1][1], cax=cax)
    cbar.set_label(ylabel)
    axs[1][1].set_xticks(nx*dx*np.arange(0, 1.25, 0.25))
    axs[1][1].set_yticks(nz*dx*np.arange(0, 1.25, 0.5))
    axs[1][1].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1)
    axs[1][1].set_title(r"x-velocity")
    axs[1][1].set_xlabel(r"$x~$[kpc]")
    axs[1][1].set_ylabel(r"$z~$[kpc]")
    axs[1][1].text(10, 15, '$Time={:.1e}~yr$'.format(t_arr[i]), bbox={'facecolor': 'white', 'pad': 8}, fontsize=10, alpha=0.6)
    # plot and save
    save = True
    if save:
        plt.savefig(outdir + "init.png", dpi=300)
    plt.close()

    print(f"Saving figures {i+1} of {len(d_dust)}.\n")