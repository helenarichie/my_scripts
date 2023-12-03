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
vlims_gas = (19.3, 21.8)
vlims_dust = (-13.6, -3.2)
fnums = [304]
pad = 0.1
spacing = 160  # pc, spacing of tick marks
fontsize = 20
labelpad = 12
tickwidth = 2
tmax = 30.4e6
plt.rcParams.update({'font.family': 'Helvetica'})
projdir = csvdir = "/Users/helenarichie/Desktop/cloud_survival/0726"
pngdir = os.path.join("/Users/helenarichie/Desktop/")
cmap_dust = sns.color_palette("rocket", as_cmap=True)
cmap_gas = sns.color_palette("mako", as_cmap=True)
# cmap = "viridis"
###################################################################################

plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': fontsize})
# plt.style.use('dark_background')

fig, axs = plt.subplots(2, 2, gridspec_kw={'width_ratios': [1.1, 4]}, figsize=(20.5,10))

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
with open(os.path.join(csvdir, "rate_cl.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_cl[i] = np.array(line, dtype=float)

rate_dust = np.zeros((n_steps, 6))  # M_sun / yr
with open(os.path.join(csvdir, "rate_dust.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_dust[i] = np.array(line, dtype=float)

data = ReadHDF5(projdir, proj="xy", dust=True, cat=cat, fnum=fnums[0])
head = data.head
conserved = data.conserved
dx = head["dx"][0] * 1e3
nx, ny, nz = head["dims"]
xlen = nx - 1*nz
indices = [nz, xlen+nz]

mass_cl = []
with open(os.path.join(csvdir, "rho_cl_tot.csv")) as f:
    for line in f:
        line = line.split(",")
        density = float(line[0]) * density_conversion
        mass = density * (dx/1e3)**3
        mass_cl.append(mass)
mass_cl = np.array(mass_cl)

mass_i = 0
mass_out_cl = []
for i, rate in enumerate(rate_cl):
    rate_i = np.sum(rate)
    mass_i += rate_i * dt
    mass_out_cl.append(mass_i)
mass_out_cl = np.array(mass_out_cl)

mass_dust = []
with open(os.path.join(csvdir, "rho_d_tot.csv")) as f:
    for line in f:
        line = line.split(",")
        density = float(line[0]) * density_conversion
        mass = density * (dx/1e3)**3
        mass_dust.append(mass)
mass_dust = np.array(mass_dust)

mass_i = 0
mass_out_dust = []
mass_lost_dust = []
for i, rate in enumerate(rate_dust):
    rate_i = np.sum(rate)
    mass_lost_dust.append(rate_i * dt)
    mass_i += rate_i * dt
    mass_out_dust.append(mass_i)
mass_out_dust = np.array(mass_out_dust)
mass_lost_dust = np.array(mass_lost_dust)

mass_destroyed = np.zeros(len(mass_lost_dust))
mass_destroyed[0] = 0
mass_dest_cum = 0
for i, mass in enumerate(mass_dust):
    if (i+1 < len(mass_lost_dust)):
        mass_dest_cum += (mass - mass_dust[i+1] - mass_lost_dust[i+1])
        mass_destroyed[i+1] = mass_dest_cum

d_dust = conserved["dust_density"][0]
d_dust *= head["mass_unit"] / (head["length_unit"] ** 2)
d_gas = conserved["density"][0] * (head["mass_unit"] / (head["length_unit"] ** 2))
n_gas = d_gas /(0.6*MP)
t = data.t_cgs() / yr_in_s  # yr
d_dust[d_dust<=0] = 1e-40
#im_gas = axs[0][1].imshow(np.log10(n_gas[indices[0]:indices[1],:].T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=[0, xlen*dx, 0, nz*dx], cmap=cmap_gas)
#im_dust = axs[1][1].imshow(np.log10(d_dust[indices[0]:indices[1],:].T), origin="lower", vmin=vlims_dust[0], vmax=vlims_dust[1], extent=[0, xlen*dx, 0, nz*dx], cmap=cmap_dust)
im_gas = axs[0][1].imshow(np.log10(n_gas[indices[0]:indices[1],:].T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=[0, xlen*dx, 0, nz*dx], cmap=cmap_gas)
im_dust = axs[1][1].imshow(np.log10(d_dust[indices[0]:indices[1],:].T), origin="lower", vmin=vlims_dust[0], vmax=vlims_dust[1], extent=[0, xlen*dx, 0, nz*dx], cmap=cmap_dust)
#im = axes.imshow(np.log10(d_dust[indices[i][0]:indices[i][1],:].T), origin="lower", extent=[0, xlen[i]*dx, 0, nz*dx], cmap=cmap)

axs[0][1].hlines(0.12*dx*ny, spacing, spacing+spacing, color='white')
axs[0][1].text(spacing+spacing+50, 0.1*dx*ny, '160 pc', color='white', fontsize=fontsize)
axs[1][1].text(spacing, 0.85*dx*nz, f'{round(t[0]/1e6, 1)} Myr', color='white', fontsize=fontsize) 

axs[0][1].set_xticks(np.arange(0, xlen*dx, spacing))
axs[0][1].set_yticks(np.arange(0, ny*dx, spacing))
axs[0][1].tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=9, width=tickwidth)
axs[1][1].set_xticks(np.arange(0, xlen*dx, spacing))
axs[1][1].set_yticks(np.arange(0, ny*dx, spacing))
axs[1][1].tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=9, width=tickwidth)

divider = make_axes_locatable(axs[0][1])
cax = divider.append_axes("right", size="5%", pad=0.09)
cbar_gas = fig.colorbar(im_gas, ax=axs[0][1], cax=cax)
cbar_gas.set_label(r'$\mathrm{log}_{10}(N_{H, gas})$ [$\mathrm{cm}^{-2}$]', rotation=270, labelpad=30, fontsize=fontsize)
cbar_gas.ax.tick_params(length=9, width=tickwidth, color="white", labelsize="medium")
cbar_gas.set_ticks(np.linspace(vlims_gas[0], vlims_gas[1], 4).round(1))

divider = make_axes_locatable(axs[1][1])
cax = divider.append_axes("right", size="5%", pad=0.09)
cbar_dust = fig.colorbar(im_dust, ax=axs[1][1], cax=cax)
cbar_dust.set_label(r'$\mathrm{log}_{10}(\Sigma_{dust})$ [$\mathrm{g}\,\mathrm{cm}^{-2}$]', rotation=270, labelpad=25, fontsize=fontsize)
cbar_dust.ax.tick_params(length=9, width=tickwidth, color="white", labelsize="medium")
cbar_dust.set_ticks(np.linspace(vlims_dust[0], vlims_dust[1], 4).round(1))

ymin = np.amin([np.amin(mass_cl), np.amin(mass_out_cl)]) - pad
ymax = np.amax([np.amax(mass_cl), np.amax(mass_out_cl)]) + pad
xmin = np.amin(t_arr[t_arr<=tmax]) - pad
xmax = np.amax(t_arr[t_arr<=tmax]) + pad

axs[0][0].plot(t_arr[t_arr<=tmax]/1e6, mass_cl[t_arr<=tmax]/mass_cl[0], linewidth=4, c="k", label="in box")
axs[0][0].plot(t_arr[t_arr<=tmax]/1e6, mass_out_cl[t_arr<=tmax]/mass_cl[0], linewidth=4, linestyle="--", c="k", label="exited box")
axs[0][0].set_xlim(xmin/1e6, xmax/1e6)
axs[0][0].set_xlabel("Time [Myr]")
axs[0][0].set_ylabel(r"$m_{cloud}/m_{cloud,i}$")
axs[0][0].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True, labelsize="medium")
axs[0][0].legend(fontsize=fontsize-2, loc="center")
axs[0][0].ticklabel_format(axis='y', style='sci', useOffset=True)

ymin = np.amin([np.amin(mass_dust), np.amin(mass_out_dust)]) - pad
ymax = np.amax([np.amax(mass_dust), np.amax(mass_out_dust)]) + pad
xmin = np.amin(t_arr[t_arr<=tmax]) - pad
xmax = np.amax(t_arr[t_arr<=tmax]) + pad

axs[1][0].plot(t_arr[t_arr<=tmax]/1e6, mass_dust[t_arr<=tmax]/mass_dust[0], linewidth=4, c="k", label="in box")
axs[1][0].plot(t_arr[t_arr<=tmax]/1e6, mass_out_dust[t_arr<=tmax]/mass_dust[0], linewidth=4, linestyle="--", c="k", label="exited box")
axs[1][0].plot(t_arr[t_arr<=tmax]/1e6, mass_destroyed[t_arr<=tmax]/mass_dust[0], linewidth=4, c="r", zorder=10, label="sputtered")
axs[1][0].set_xlim(xmin/1e6, xmax/1e6)
axs[1][0].set_xlabel("Time [Myr]")
axs[1][0].set_ylabel(r"$m_{dust}/m_{dust,i}$")
axs[1][0].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True, labelsize="medium")
axs[1][0].legend(fontsize=fontsize-2, loc="center")
axs[1][0].ticklabel_format(axis='y', style='sci', useOffset=True)

print(mass_cl[0])

plt.tight_layout()
plt.savefig(pngdir + f"dust_survival.png", dpi=300)
plt.close()