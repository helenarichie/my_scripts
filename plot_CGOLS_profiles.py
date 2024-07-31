import sys
sys.path.insert(0, "../github/my_scripts")
from hconfig import *
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl
from csv import writer

import numpy as np
import os
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 25})
from mpl_toolkits.axes_grid1 import make_axes_locatable

import seaborn as sns

############ hard-coded ###########
write = True
new = False
r_cl = 1 * 3.086e+18 # pc to cm
d_cl_init = 1e-24 # n = 1
r_grain = 0.1 # micron
gamma = 5/2
T_cool = 1e4
chi = 1000
v_wind_i = 1000e3 # cm/s
tau_cc = (np.sqrt(chi)*r_cl/v_wind_i)/yr_in_s # cloud crushing time
colors = sns.color_palette("cubehelix", 7)
mu = 0.6
plt.style.use('dark_background')
plt.rcParams.update({'font.family': 'Helvetica'})
##################################

bins, r_av, n_av, n_med, n_lo, n_hi, v_av, v_med, v_lo, v_hi, T_av, T_med, T_lo, T_hi, p_av, p_med, p_lo, p_hi, c_av, c_med, c_lo, c_hi, cs_av, cs_med, cs_lo, cs_hi, K_av, K_med, K_lo, K_hi, M_av, M_med, M_lo, M_hi = np.loadtxt('/Users/helenarichie/Downloads/CGOLS profiles/2048_central_35_hot_dweight.txt', unpack=True)

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(13,15))

ax[0].plot(r_av, np.log10(T_med), c="darkorange", linewidth=5)
ax[0].fill_between(r_av, np.log10(T_med-T_lo), np.log10(T_med+T_hi), alpha=0.2, color="darkorange")
ax[0].set_xlabel("$r~[kpc]$")
ax[0].set_ylabel(r"$T_{hot}~[K]$")
ax[0].set_xlim(0, 10)
ax[0].set_ylim(5.5, 8)
ax[0].tick_params(top=True, right=True)

ax[1].semilogy(r_av, n_med, c="blueviolet", linewidth=5)
ax[1].fill_between(r_av, n_med-n_lo, n_med+n_hi, alpha=0.2, color="blueviolet")
ax[1].set_xlabel("$r~[kpc]$")
ax[1].set_ylabel(r"$n_{hot}~[cm^{-3}]$")
ax[1].set_xlim(0, 10)
ax[1].set_ylim(1e-5, 1e-1)
ax[1].tick_params(top=True, right=True)

# plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

plt.savefig(f"/Users/helenarichie/Downloads/figs/CGOLS_hot.png")
plt.close()




