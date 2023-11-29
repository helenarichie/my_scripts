import sys
sys.path.insert(0, "../github/my_scripts")
from hconfig import *
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl

import numpy as np
import os
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 25})
from mpl_toolkits.axes_grid1 import make_axes_locatable

import seaborn as sns

############ hard-coded ###########
colors = sns.color_palette("cubehelix", 7)
mu = 0.6
##################################

# McKinnon et al. (2017)
def calc_tau_sp(n, T):
    YR_IN_S = 3.154e7;
    a1 = 1; # dust grain size in units of 0.1 micrometers
    d0 = n / (6e-4); # gas density in units of 10^-27 g/cm^3
    T_0 = 2e6; # K
    omega = 2.5;
    A = 0.17e9 * YR_IN_S; # Gyr in s

    tau_sp = A * (a1/d0) * ((T_0/T)**omega + 1); # s

    return tau_sp

# McKinnon et al. (2017)
def calc_dd_dt(d_dust, tau_sp):
    return -d_dust / (tau_sp/3)

a = np.array(range(2, 9)) # log K
T_range = 10**a # K

pairs = []

tmax = 1e9 # yr
h = 1e2 # yr
n_gas = 0.01808641 # cm^-3
d_dust_init = 0.01*n_gas*mu*MP # g cm^-3
t_arr = np.arange(0, tmax, h)

fig = plt.figure(figsize=(11,8))
plt.style.use('dark_background')
plt.rcParams.update({'font.family': 'Helvetica'})

plt.axvline(x=1e6, color='lightgrey', linestyle="--", zorder=0, linewidth=4, label="1 Myr")

for j, T in enumerate(T_range):
    tau_sp = calc_tau_sp(n_gas, T) / yr_in_s # yr
    d_dust = d_dust_init
    d_dust_solution = []

    for i, t in enumerate(t_arr):
        dt = h
        dd_dt = calc_dd_dt(d_dust, tau_sp)
        dd = dd_dt * dt

        while dd/d_dust > 0.01:
            dt_sub = 0.01 * d_dust / dd_dt
            d_dust += dt_sub * dd_dt
            dt -= dt_sub
            dd_dt = calc_dd_dt()
            dd = dt * dd_dt
        
        d_dust += dd

        d_dust_solution.append(d_dust)

    d_dust_solution = np.array(d_dust_solution)
    # d_dust_solution /= MP * mu

    wh_knee = (np.abs(t_arr - tau_sp)).argmin()

    plt.semilogx(t_arr, d_dust_solution/d_dust_init, lw=10, color='k')
    plt.semilogx(t_arr, d_dust_solution/d_dust_init, c=colors[j], linewidth=5, label=r'$10^{{{:d}}}$ K'.format(a[j]))
    #if wh_knee != np.argmax(t_arr):
    #    plt.scatter(t_arr[wh_knee], d_dust_solution[wh_knee], marker="x", color=colors[j], s=100, zorder=150, linewidths=5)
    #    plt.title(r"$\rho_{gas}=$" + f"{n_gas*mu*MP:.1e}" + r"$~[g\,cm^{-3}]$")

# plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
#plt.scatter(-10, -10, marker="x", s=150, label=r"$t_{sp}$", linewidths=5, c="k")
#plt.xlim(np.amin(t_arr)+0.01, np.amax(t_arr))
#plt.ylim(0 - d_dust_init/10, d_dust_init + d_dust_init/100)
plt.legend(fontsize=18, loc="lower left")
plt.xlabel("Time [yr]")
plt.ylabel(r"$\rho_{dust}/\rho_{dust,i}$")
plt.title("Dust sputtering in \rho=")
plt.savefig("/Users/helenarichie/Desktop/sputtering_curves.png")




