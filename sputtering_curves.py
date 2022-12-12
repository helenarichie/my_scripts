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


T_range = np.array(range(2, 9)) # log K
T_range = 10**T_range # K

pairs = []

tmax = 1e9 # yr
h = 1e2 # yr
n_gas = 0.1 # cm^-3
d_dust_init = 0.01*n_gas*mu*MP # g cm^-3
t_arr = np.arange(0, tmax, h)

fig = plt.figure(figsize=(11,8))

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

    plt.semilogx(t_arr, d_dust_solution, lw=10, color='k')
    plt.semilogx(t_arr, d_dust_solution, c=colors[j], linewidth=5, label=r"$T=$"+f"{T:.0e} K")
    plt.title(r"$\rho_{gas}=$" + f"{n_gas*mu*MP:.1e}" + r"$~[g\,cm^{-3}]$")

# plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
plt.legend(fontsize=18, loc="lower left")
plt.xlabel("Time [yr]")
plt.ylabel(r"$\rho_{dust}~[g\,cm^{-3}]$")
plt.savefig("figs/sputtering_curves.png")




