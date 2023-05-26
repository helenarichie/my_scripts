from hconfig import *
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl

import numpy as np
import os
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 30})
from mpl_toolkits.axes_grid1 import make_axes_locatable

import seaborn as sns

hist_cmap = sns.cubehelix_palette(light=1, as_cmap=True, reverse=True)

############ hard-coded ###########
r_cl = 1 * 3.086e+18 # pc to cm
d_cl_init = 1e-24 # n = 1
r_grain = 0.1 # micron
gamma = 5/2
T_cool = 1e4 # K
chi = 1000
v_wind_i = 1000e3 # cm/s
tau_cc = (np.sqrt(chi)*r_cl/v_wind_i)/yr_in_s # cloud crushing time
T_max = 1e8
T_min = 1e2
n_min = 1e-5
n_max = 1e3
#############################################

# McKinnon et al. (2017)
def calc_tau_sp(n, T):
    YR_IN_S = 3.154e7;
    a1 = 1; # dust grain size in units of 0.1 micrometers
    d0 = n / (6e-4); # gas density in units of 10^-27 g/cm^3
    T_0 = 2e6; # K
    omega = 2.5;
    A = 0.17e9 * YR_IN_S; # Gyr in s

    tau_sp = A * (a1/d0) * ((T_0/T)**omega + 1) # s

    return tau_sp;

extent = np.log10([n_min, n_max, T_min, T_max])
n_bins = np.linspace(np.log10(n_min), np.log10(n_max), 250)
T_bins = np.linspace(np.log10(T_min), np.log10(T_max), 250)
T_arr = np.linspace(T_min, T_max, 1000)
n_arr = np.linspace(n_min, n_max, 1000)

nn, TT = np.meshgrid(n_arr, T_arr)

tau_sp = np.log10(calc_tau_sp(nn, TT)/3.154e7)

fig, ax = plt.subplots(figsize=(12, 10))
im = ax.contourf(np.log10(nn), np.log10(TT), tau_sp-6, cmap=sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True), alpha=0.5, extent=extent, origin="lower")
vis_contours = plt.contour(im, levels=[0.98], colors='k', origin="lower", linestyles="dashed")
ax.scatter(2, 4, marker="x", c="k", s=100, linewidths=2)
ax.text(1.63, 3.7, "cloud", fontsize=25)
ax.scatter(-2, 7.5, marker="x", c="k", s=100, linewidths=2)
ax.text(-2.91, 7.2, "starburst wind", fontsize=25)
ax.scatter(-2, 6.5, marker="x", c="k", s=100, linewidths=2)
ax.text(-2.35, 6.2, "wind", fontsize=25)

divider = make_axes_locatable(ax)

cax = divider.append_axes('right', size='5%', pad=0.1)
cbar = fig.colorbar(im, cax=cax, orientation='vertical')
# cbar.add_lines(vis_contours)
cbar.set_label(r'Sputtering Timescale $[\log_{10}(Myr)]$', rotation=270, labelpad=30)

cbar.ax.tick_params(length=9, width=2)
ax.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2)


ax.set_xlabel(r'$\mathrm{log}_{10}(T)$ [$\mathrm{cm}^{-3}$]')
ax.set_ylabel(r'$\mathrm{log}_{10}(T) \, [\mathrm{K}]$')

plt.tight_layout()

plt.savefig("/Users/helenarichie/Desktop/sputtering_contours.png", dpi=300, bbox_inches='tight', pad_inches=0)  