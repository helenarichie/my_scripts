uuimport sys
sys.path.insert(0, "/Users/helenarichie/GitHub/my_scripts/")
from hconfig import *
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': 20})

##################################################################
cat = True
pad = 0.1
fontsize = 20
labelpad = 12
tickwidth = 1
##################################################################

##################################################################
projdir = ["/Users/helenarichie/Desktop/cloud_survival/0426", "/Users/helenarichie/Desktop/cloud_survival/0726"]
csvdir = ["/Users/helenarichie/Desktop/cloud_survival/0426", "/Users/helenarichie/Desktop/cloud_survival/0726"]
pngdir = os.path.join("/Users/helenarichie/Desktop/")
##################################################################

##################################################################
tmax = [2.1e6, 30.4e6]
fnums = [[0, 715, 1051], [0, 196, 304]]
snapshot_times = [[1.4e6, 2.1e6], [19.6e6, 30.3e6]]
##################################################################

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6.5, 12.5))

for i, dir in enumerate(csvdir):
    # get value of dx to convert from density to mass
    data = ReadHDF5(projdir[i], proj="xy", cat=cat, fnum=fnums[0])
    head = data.head
    conserved = data.conserved
    dx = head["dx"][0]  # kpc

    # get times of snapshots
    t_arr = [] # yr
    with open(os.path.join(csvdir[i], "t_arr.csv")) as f:
        for line in f:
            line = line.split(",")
            t = line[0].strip("\n").strip("[").strip("]")
            t_arr.append(float(t))
    t_arr = np.array(t_arr)

    n_steps = len(t_arr)

    # get timestep
    dt = (t_arr[1] - t_arr[0])/1e3

    # get outflow rates for gas
    rate_cl = np.zeros((n_steps, 6))  # M_sun / yr
    with open(os.path.join(csvdir[i], "rate_cl.csv")) as f:
        for j, line in enumerate(f):
            line = line.split(",")
            rate_cl[j] = np.array(line, dtype=float)

    # get total cloud mass in simulation volume
    mass_cl = []
    with open(os.path.join(csvdir[i], "rho_cl_tot.csv")) as f:
        for line in f:
            line = line.split(",")
            density = float(line[0]) * density_conversion
            mass = density * dx**3
            mass_cl.append(mass)
    mass_cl = np.array(mass_cl)

    # use outflow rates to calculate how much cloud has exited the simulation volume
    mass_i = 0
    mass_out = []
    for rate in rate_cl:
        rate_i = np.sum(rate)
        mass_i += rate_i * dt
        mass_out.append(mass_i)
    mass_out = np.array(mass_out)

    ymin = np.amin([np.amin(mass_cl), np.amin(mass_out)]) - pad
    ymax = np.amax([np.amax(mass_cl), np.amax(mass_out)]) + pad
    xmin = np.amin(t_arr[t_arr<=tmax[i]]) - pad
    xmax = np.amax(t_arr[t_arr<=tmax[i]]) + pad

    mass_cl_init = mass_cl[0]
    print(f"initial mass: {mass_cl_init:e}")
        
    if i == 0:
        indices = [np.argmin(t_arr), np.argmin(np.abs(t_arr-snapshot_times[0][0])), np.argmin(np.abs(t_arr-snapshot_times[0][1]))]
    if i == 1:
        indices = [np.argmin(t_arr), np.argmin(np.abs(t_arr-snapshot_times[1][0])), np.argmin(np.abs(t_arr-snapshot_times[1][1]))]

    ax[i].plot(t_arr[t_arr<=tmax[i]]/1e6, mass_cl[t_arr<=tmax[i]], linewidth=4, c="#49b4ab", label="in box")
    ax[i].plot(t_arr[t_arr<=tmax[i]]/1e6, mass_out[t_arr<=tmax[i]], linewidth=4, linestyle="--", c="#49b4ab", label="exited box")
    ax[i].scatter(t_arr[indices]/1e6, mass_cl[indices], marker="o", c="#49b4ab", zorder=11, s=50, linewidths=1.5, edgecolors="k")

    ax[i].set_xlabel("Time [Myr]")
    ax[i].set_ylabel(r"Cloud Mass [$M_\odot$]")

    ax[i].set_xlim(xmin/1e6, xmax/1e6)
    ax[i].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True, labelsize="medium", zorder=10)
    ax[i].ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    if i == 0:
        ax[i].legend(fontsize=fontsize-5, loc="upper right")

plt.tight_layout()
plt.savefig(pngdir + f"cloud_mass.png", dpi=300)
plt.close()