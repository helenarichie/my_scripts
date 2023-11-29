uuimport sys
sys.path.insert(0, "/Users/helenarichie/GitHub/my_scripts/")
from hconfig import *
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

cat = True
vlims_gas_small = (18.2, 19.6)
vlims_gas_large = (19.3, 22)
# vlims_gas_small = vlims_gas_large = (18.2, 21.8)
fnums_small = [0, 715, 1051]
fnums_large = [0, 196, 304]
pad = 0.1
spacing_small = 20  # pc, spacing of tick marks
spacing_large = 320  # pc, spacing of tick marks
fontsize = 20
labelpad = 12
tickwidth = 1
tmax_small = 2.1e6
tmax_large = 30.4e6
plt.rcParams.update({'font.family': 'Helvetica'})
projdir_small = csvdir_small = "/Users/helenarichie/Desktop/cloud_survival/0426"
projdir_large = csvdir_large = "/Users/helenarichie/Desktop/cloud_survival/0726"
pngdir = os.path.join("/Users/helenarichie/Desktop/")

plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': fontsize})

fnums = [fnums_small, fnums_large]
projdir = [projdir_small, projdir_large]
csvdir = [csvdir_small, csvdir_large]
vlims = [vlims_gas_small, vlims_gas_large]
spacing = [spacing_small, spacing_large]
tmax = [tmax_small, tmax_large]

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6.5, 12.5))

for i, dir in enumerate(csvdir):
    data = ReadHDF5(projdir[i], proj="xy", cat=cat, fnum=fnums[0])
    head = data.head
    conserved = data.conserved
    dx = head["dx"][0] * 1e3  # pc

    t_arr = [] # yr
    with open(os.path.join(csvdir[i], "t_arr.csv")) as f:
        for line in f:
            line = line.split(",")
            t = line[0].strip("\n").strip("[").strip("]")
            t_arr.append(float(t))
    t_arr = np.array(t_arr)

    n_steps = len(t_arr)

    dt = (t_arr[1] - t_arr[0])/1e3

    rate_cl = np.zeros((n_steps, 6))  # M_sun / yr
    with open(os.path.join(csvdir[i], "rate_cl.csv")) as f:
        for j, line in enumerate(f):
            line = line.split(",")
            rate_cl[j] = np.array(line, dtype=float)

    mass_cl = []
    with open(os.path.join(csvdir[i], "rho_cl_tot.csv")) as f:
        for line in f:
            line = line.split(",")
            density = float(line[0]) * density_conversion
            mass = density * (dx/1e3)**3
            mass_cl.append(mass)
    mass_cl = np.array(mass_cl)

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
    
    if i == 0:
        indices = [np.argmin(t_arr), np.argmin(np.abs(t_arr-1.4e6)), np.argmin(np.abs(t_arr-2.1e6))]
    if i == 1:
        indices = [np.argmin(t_arr), np.argmin(np.abs(t_arr-19.6e6)), np.argmin(np.abs(t_arr-30.4e6))]

    ax[i].plot(t_arr[t_arr<=tmax[i]]/1e6, mass_cl[t_arr<=tmax[i]], linewidth=4, c="#49b4ab", label="in box")
    ax[i].scatter(t_arr[indices]/1e6, mass_cl[indices], marker="o", c="#49b4ab", zorder=1, s=50, linewidths=4)
    ax[i].plot(t_arr[t_arr<=tmax[i]]/1e6, mass_out[t_arr<=tmax[i]], linewidth=4, linestyle="--", c="#49b4ab", label="exited box")
    ax[i].set_xlim(xmin/1e6, xmax/1e6)
    ax[i].set_xlabel("Time [Myr]")
    ax[i].set_ylabel(r"Cloud Mass [$M_\odot$]")
    ax[i].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True, labelsize="medium", zorder=10)
    ax[i].legend(fontsize=fontsize-2, loc="center left")
    ax[i].ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.tight_layout()
plt.savefig(pngdir + f"cloud_mass.png", dpi=300)
plt.close()