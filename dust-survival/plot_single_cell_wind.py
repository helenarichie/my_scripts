import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "/Users/helenarichie/GitHub/my_scripts/")
from hconfig import *

plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': 20})

csvdir = "/Users/helenarichie/Desktop/"
fontsize = 20

# get times of snapshots
t, r, rho_d = [], [], []
with open(os.path.join(csvdir, "single_cell_wind.csv")) as f:
    for line in f:
        line = line.split(",")
        t.append(float(line[0]))
        r.append(float(line[1]))
        rho_d.append(float(line[2]))

rho_d = np.array(rho_d)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 6.5))
ax.semilogx(t, rho_d/rho_d[0], linewidth=4, c="k")
ax.set_xlabel("Time [yr]", fontsize=fontsize)
ax.set_ylabel(r"$\rho_{dust}/\rho_{dust,i}$", fontsize=fontsize+2)
plt.tight_layout()
ax.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True)
plt.savefig("/Users/helenarichie/Desktop/single_cell_wind_t.png")


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 6.5))
ax.plot(np.log10(r), rho_d/rho_d[0], linewidth=4, c="k")
ax.set_xlim(np.amin(np.log10(r)), np.amax(np.log10(r)))
ax.set_xticks(np.linspace(np.amin(np.log10(r)), np.amax(np.log10(r)), 5).round(1))
ax.set_xlabel("$\log(r)$ [kpc]", fontsize=fontsize)
ax.set_ylabel(r"$\rho_{dust}/\rho_{dust,i}$", fontsize=fontsize+2)
ax.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True)
plt.tight_layout()
plt.savefig("/Users/helenarichie/Desktop/single_cell_wind_r.png")
