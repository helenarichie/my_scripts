import sys
sys.path.insert(0, "/Users/helenarichie/GitHub/my_scripts/")
from hconfig import *
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

###################################################################################
cat = True
vlims_gas_small = (-10, -6)
vlims_gas_large = (-11, -3.5)
fnums_small = [0, 715, 1051]
fnums_large = [0, 196, 304]
pad = 0.1
spacing_small = 10  # pc, spacing of tick marks
spacing_large = 160  # pc, spacing of tick marks
fontsize = 20
labelpad = 12
tickwidth = 2
tmax_small = 2.1e6
tmax_large = 30.4e6
plt.rcParams.update({'font.family': 'Helvetica'})
projdir_small = csvdir_small = "/Users/helenarichie/Desktop/cloud_survival/0426"
projdir_large = csvdir_large = "/Users/helenarichie/Desktop/cloud_survival/0726"
pngdir = os.path.join("/Users/helenarichie/Desktop/")
cmap = sns.color_palette("rocket", as_cmap=True)
# cmap = "viridis"
###################################################################################

plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': fontsize})
# plt.style.use('dark_background')

fig, axs = plt.subplots(constrained_layout=True, figsize=(30,10))
(subfig_small, subfig_large) = fig.subfigures(nrows=2, ncols=1)

subfigsnest_small = subfig_small.subfigures(nrows=1, ncols=2, width_ratios=[5, 1])
axes_small_proj = subfigsnest_small[0].subplots(1, 3, gridspec_kw={'width_ratios': [1, 2, 2]})
axes_small_mass = subfigsnest_small[1].subplots(1, 1)

subfigsnest_large = subfig_large.subfigures(nrows=1, ncols=2, width_ratios=[5, 1])
axes_large_proj = subfigsnest_large[0].subplots(1, 3, gridspec_kw={'width_ratios': [1, 2, 2]})
axes_large_mass = subfigsnest_large[1].subplots(1, 1)

def plot_panel(projdir, csvdir, vlims, fnums, subfig, tmax, axes_proj, axes_mass, spacing, id):
    data = ReadHDF5(projdir, proj="xy", dust=True, cat=cat, fnum=fnums[0])
    head = data.head
    conserved = data.conserved
    dx = head["dx"][0]
    nx, ny, nz = head["dims"]

    xlen = [nz, nx - 2*nz, nx - 2*nz]
    if id == 0:
        indices = [[0, nz], [nz, 3*nz], [2*nz, 4*nz]]
    if id == 1:
        indices = [[0, nz], [0+64+64, 2*nz+64+64], [nz, 3*nz]]

    t_arr = [] # yr
    with open(os.path.join(csvdir, "t_arr.csv")) as f:
        for line in f:
            line = line.split(",")
            t = line[0].strip("\n").strip("[").strip("]")
            t_arr.append(float(t))
    t_arr = np.array(t_arr)

    n_steps = len(t_arr)

    dt = (t_arr[1] - t_arr[0])/1e3

    rate_cl = np.zeros((n_steps, 6))  # M_sun / yr
    with open(os.path.join(csvdir, "rate_dust.csv")) as f:
        for i, line in enumerate(f):
            line = line.split(",")
            rate_cl[i] = np.array(line, dtype=float)

    mass_cl = []
    with open(os.path.join(csvdir, "rho_d_tot.csv")) as f:
        for line in f:
            line = line.split(",")
            density = float(line[0]) * density_conversion
            mass = density * dx**3
            mass_cl.append(mass)
    mass_cl = np.array(mass_cl)

    mass_i = 0
    mass_out = []
    for i, rate in enumerate(rate_cl):
        rate_i = np.sum(rate)
        mass_i += rate_i * dt
        mass_out.append(mass_i)
    mass_out = np.array(mass_out)

    for i, axes in enumerate(axes_proj):
        data = ReadHDF5(projdir, proj="xy", dust=True, cat=cat, fnum=fnums[i])
        head = data.head
        conserved = data.conserved
        dx = head["dx"][0] * 1e3  # pc
        d_dust = conserved["dust_density"][0]
        d_dust *= head["mass_unit"] / (head["length_unit"] ** 2)
        t = data.t_cgs() / yr_in_s  # yr
        d_dust[indices[i][0]:indices[i][1],:][d_dust[indices[i][0]:indices[i][1],:]<=0] = 1e-40
        print(np.amax(d_dust[indices[i][0]:indices[i][1],:]))
        im = axes.imshow(np.log10(d_dust[indices[i][0]:indices[i][1],:].T), origin="lower", vmin=vlims[0], vmax=vlims[1], extent=[0, xlen[i]*dx, 0, nz*dx], cmap=cmap)
        #im = axes.imshow(np.log10(d_dust[indices[i][0]:indices[i][1],:].T), origin="lower", extent=[0, xlen[i]*dx, 0, nz*dx], cmap=cmap)

        if i == 0:
            axes.hlines(0.12*dx*ny, spacing, spacing+spacing, color='white')
            if id == 0:
                axes.text(spacing+spacing+2, 0.1*dx*ny, '10 pc', color='white', fontsize=fontsize)
            if id == 1:
                axes.text(spacing+spacing+50, 0.1*dx*ny, '160 pc', color='white', fontsize=fontsize)
        axes.text(spacing, 0.85*dx*nz, f'{round(t[0]/1e6, 1)} Myr', color='white', fontsize=fontsize) 

        axes.set_xticks(np.arange(0, xlen[i]*dx, spacing))
        axes.set_yticks(np.arange(0, nz*dx, spacing))
        axes.tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=9, width=tickwidth)

    cbar = subfig[0].colorbar(im, ax=axes_proj, location='right', pad=0.009, shrink=0.955, aspect=10)
    cbar.set_label(r'$\mathrm{log}_{10}(\Sigma_{dust})$ [$g\,\mathrm{cm}^{-2}$]', rotation=270, labelpad=20, fontsize=fontsize)
    cbar.ax.tick_params(length=9, width=tickwidth, color="white", labelsize="medium")
    cbar.set_ticks(np.linspace(vlims[0], vlims[1], 4).round(1))

    ymin = np.amin([np.amin(mass_cl), np.amin(mass_out)]) - pad
    ymax = np.amax([np.amax(mass_cl), np.amax(mass_out)]) + pad
    xmin = np.amin(t_arr[t_arr<=tmax]) - pad
    xmax = np.amax(t_arr[t_arr<=tmax]) + pad

    axes_mass.plot(t_arr[t_arr<=tmax]/1e6, mass_cl[t_arr<=tmax]/mass_cl[0], linewidth=4, c="k", label="in volume")
    axes_mass.plot(t_arr[t_arr<=tmax]/1e6, mass_out[t_arr<=tmax]/mass_cl[0], linewidth=4, linestyle="--", c="k", label="exited volume")
    axes_mass.set_xlim(xmin/1e6, xmax/1e6)
    axes_mass.set_xlabel("Time [Myr]")
    axes_mass.set_ylabel(r"Fraction of Initial Mass")
    axes_mass.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True, labelsize="medium")
    axes_mass.legend(fontsize=fontsize-2, loc="center")
    axes_mass.ticklabel_format(axis='y', style='sci', useOffset=True)

    print(mass_cl[0])

plot_panel(projdir_small, csvdir_small, vlims_gas_small, fnums_small, subfigsnest_small, tmax_small, axes_small_proj, axes_small_mass, spacing_small, 0)
plot_panel(projdir_large, csvdir_large, vlims_gas_large, fnums_large, subfigsnest_large, tmax_large, axes_large_proj, axes_large_mass, spacing_large, 1)

plt.savefig(pngdir + f"dust_survival.png", dpi=300)
plt.close()