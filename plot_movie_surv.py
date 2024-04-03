import sys
sys.path.insert(0, "/ccs/home/helenarichie/code/my_scripts/")
from hconfig import *
import seaborn as sns

##################################
date = "2024-02-07"
ns = 1525
ne = 1550
cat = True
##################################

########### location #############
crc = False
frontier = True
##################################

########## data type #############
debugging = False
cloud_wind = False
testing = True
##################################

########### plotting #############
dust = True
vlims_gas = (19.4, 22)
vlims_dust = (-8.5, -4.3)
spacing, unit = 400*1e-3, "kpc"
tmax = 77e3
##################################

##################################
pad = 0.1
fontsize = 30
labelpad = 12
tickwidth = 2.25
ticklength = 14
linewidth = 4
##################################

########## specify slice and png directories ############
if crc:
    if debugging:
        basedir = f"/ix/eschneider/helena/data/debugging/{date}/"
    if cloud_wind:
        basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
    if testing:
        basedir = f"/ix/eschneider/helena/data/testing/{date}/"

if frontier:
    basedir = f"/lustre/orion/ast181/scratch/helenarichie/{date}/"

if cat:
    datadir = os.path.join(basedir, "hdf5/proj/")
else:
    datadir = os.path.join(basedir, "hdf5/raw/")

pngdir = os.path.join(basedir, "png/movie")
csvdir = os.path.join(basedir, "csv/")
#########################################################

mass_cloud_unsorted = []
time_cloud_unsorted = []
with open(os.path.join(csvdir, "full_mass.csv"), "r") as f:
    for line in f:
        line = line.split(",")
        mass_cloud_unsorted.append(float(line[2]))
        time_cloud_unsorted.append(float(line[1]))

mass_cloud = [x for _, x in sorted(zip(time_cloud_unsorted, mass_cloud_unsorted))]
time_cloud = sorted(time_cloud_unsorted)
mass_cloud = np.array(mass_cloud)
time_cloud = np.array(time_cloud)

"""
mass_cloud = []
with open(os.path.join(csvdir, "mass_cloud.csv")) as f:
    for line in f:
        mass_cloud.append(float(line))
mass_cloud = np.array(mass_cloud)
"""

time = []
with open(os.path.join(csvdir, "time.csv")) as f:
    for line in f:
        time.append(float(line))
time = np.array(time)

mass_dust = []
with open(os.path.join(csvdir, "mass_dust.csv")) as f:
    for line in f:
        mass_dust.append(float(line))
mass_dust = np.array(mass_dust)

sputter_hot = []
with open(os.path.join(csvdir, "sputter_hot.csv")) as f:
    for line in f:
        sputter_hot.append(float(line))

sputter = []
with open(os.path.join(csvdir, "sputter.csv")) as f:
    for line in f:
        sputter.append(float(line))
"""
"""
time_output = []
with open(os.path.join(csvdir, "time_output.csv")) as f:
    for line in f:
        time_output.append(float(line))
time_output = np.array(time_output)

rate_dust = []
with open(os.path.join(csvdir, "rate_dust.csv")) as f:
    for line in f:
        line = line.split(",")
        rate_dust.append(np.array(line, dtype=float))

dt_out = time_output[2] - time_output[1]

rate_cloud = []
with open(os.path.join(csvdir, "rate_cloud.csv")) as f:
    for line in f:
        line = line.rstrip("\n").split(",")
        rate_cloud.append(np.array(line, dtype=float))

mass_out_dust = []
mass_dust_cumulative = 0
for i, rate in enumerate(rate_dust):
    rate = np.sum(rate)
    mass_dust_cumulative += rate * dt_out
    mass_out_dust.append(mass_dust_cumulative)
mass_out_dust = np.array(mass_out_dust)

mass_out_cloud = []
mass_cloud_cumulative = 0
for i, rate in enumerate(rate_cloud):
    rate = np.sum(rate)
    mass_cloud_cumulative += rate * dt_out
    mass_out_cloud.append(mass_cloud_cumulative)
mass_out_cloud = np.array(mass_out_cloud)
"""
"""

sputter_tot, sputter_tot_hot = [], []
mass_cumulative, mass_cumulative_hot = 0, 0
for i, mass in enumerate(sputter):
    mass_cumulative += mass
    mass_cumulative_hot += sputter_hot[i]
    sputter_tot.append(mass_cumulative)
    sputter_tot_hot.append(mass_cumulative_hot)
sputter_tot = np.array(sputter_tot)
sputter_tot_hot = np.array(sputter_tot_hot)

for i in range(ns, ne+1):
    cmap_gas = sns.color_palette("mako", as_cmap=True)
    cmap_dust = sns.color_palette("rocket", as_cmap=True)
    # read in data

    f = h5py.File(os.path.join(datadir, str(i)+"_proj.h5"), "r")
    
    head = f.attrs
    dx = head["dx"][0]
    t = head["t"][0]
    gamma = head["gamma"][0]
    nx, ny, nz = head["dims"][0], head["dims"][1], head["dims"][2]
    mu = 0.6

    xlen, ylen = nx, ny
    xs, ys = 0, 0

    density = np.array(f["d_xy"]) * head["mass_unit"] / (head["length_unit"] ** 2)
    d_dust = np.array(f["d_dust_xy"]) * head["mass_unit"] / (head["length_unit"] ** 2)
    d_dust[d_dust==0] = 1e-40

    n_gas = density/(0.6*MP)

    # plot
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(40,9), gridspec_kw={"width_ratios":[4, 30], 'wspace':0, 'hspace':0})

    im_g = ax[0][1].imshow(np.log10(n_gas[xs:(xs+xlen+1),ys:(ys+ylen+1)].T), origin="lower", extent=[0, xlen*dx, 0, ylen*dx], cmap=cmap_gas, vmin=vlims_gas[0], vmax=vlims_gas[1], aspect='auto')
    im_d = ax[1][1].imshow(np.log10(d_dust[xs:(xs+xlen+1),ys:(ys+ylen+1)].T), origin="lower", extent=[0, xlen*dx, 0, ylen*dx], cmap=cmap_dust, vmin=vlims_dust[0], vmax=vlims_dust[1], aspect='auto')

    #ax[0][1].set_yticks([spacing-0.25*spacing, 2*spacing-0.25*spacing])
    #ax[1][1].set_yticks([spacing-0.25*spacing, 2*spacing-0.25*spacing])
    ax[1][1].set_yticks(np.arange(0, ny*dx, spacing))
    ax[0][1].set_yticks(np.arange(0, ny*dx, spacing))
    ax[1][1].set_xticks(np.arange(0, nx*dx, spacing))
    ax[0][1].set_xticks(np.arange(0, nx*dx, spacing))
    ax[0][0].tick_params(axis='both', which='both', direction='in', color='black', labelbottom=0, top=1, right=1, length=ticklength, width=tickwidth)
    ax[0][1].tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, top=1, right=1, length=ticklength, width=tickwidth)
    ax[1][0].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=ticklength, width=tickwidth)
    ax[1][1].tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=ticklength, width=tickwidth)

    ax[0][1].text(0.5*spacing, 0.81*dx*ylen, f'{round(t/1e3, 1)} Myr', color='white', fontsize=fontsize)
    ax[0][1].text(0.6*spacing, 0.125*dx*ylen, 'gas', color='white', fontsize=fontsize)

    ax[1][1].hlines(0.86*dx*ylen, 0.5*spacing, 1.5*spacing, color='white')
    ax[1][1].text(1.5*spacing+0.05, 0.83*dx*ylen, '400 pc', color='white', fontsize=fontsize)
    ax[1][1].text(0.65*spacing, 0.125*dx*ylen, 'dust', color='white', fontsize=fontsize)

    ax[0][1].spines['bottom'].set_color('white')
    ax[0][1].spines['right'].set_color('white')
    ax[1][1].spines['top'].set_color('white')
    ax[1][1].spines['right'].set_color('white')

    f.close()

    cbar_height = 0.385
    cbar_ax = fig.add_axes([0.90, 0.495, 0.025, cbar_height])
    cbar = fig.colorbar(im_g, cax=cbar_ax, pad=0.1)
    cbar.set_label(r'$\mathrm{log}_{10}(N_{H, gas})$ [$\mathrm{cm}^{-2}$]', rotation=270, labelpad=50, fontsize=fontsize-4)
    cbar.ax.tick_params(length=ticklength, width=tickwidth, color="white", labelsize=fontsize-7)
    cbar_spacing = (vlims_gas[1] - vlims_gas[0])/8
    cbar.set_ticks(np.linspace(vlims_gas[0]+cbar_spacing, vlims_gas[1]-cbar_spacing, 4).round(1))
    #cbar.outline.set_linewidth(1.5)
    cbar.ax.zorder = -1

    cbar_ax = fig.add_axes([0.90, 0.11, 0.025, cbar_height])
    cbar = fig.colorbar(im_d, cax=cbar_ax)
    cbar.set_label(r'$\mathrm{log}_{10}(\Sigma_{dust})$ [$\mathrm{g}\,\mathrm{cm}^{-2}$]', rotation=270, labelpad=50, fontsize=fontsize-4)
    cbar.ax.tick_params(length=ticklength, width=tickwidth, color="white", labelsize=fontsize-7)
    cbar_spacing = (vlims_dust[1] - vlims_dust[0])/8
    cbar.set_ticks(np.linspace(vlims_dust[0]+cbar_spacing, vlims_dust[1]-cbar_spacing, 4).round(1))
    #cbar.outline.set_linewidth(1.5)
    cbar.ax.zorder = -1

    mass_dust_init = mass_dust[0]
    time_i = time[time<=t]
    mass_dust_i = mass_dust[time<=t]
    time_output_i = time_output[time_output<=t]
    mass_out_dust_i = mass_out_dust[time_output<=t]
    sputter_tot_i = sputter_tot[time<=t]
    sputter_tot_hot_i = sputter_tot_hot[time<=t]

    new = "#d43a4f"
    ax[1][0].plot(time_i/1e3, mass_dust_i/mass_dust_init, label="total", linewidth=linewidth, c=new, zorder=0)
    ax[1][0].plot(time_output_i/1e3, mass_out_dust_i/mass_dust_init, linestyle="--", linewidth=linewidth-1, c=new, zorder=0)
    ax[1][0].plot(time_i/1e3, (sputter_tot_i+sputter_tot_hot_i)/mass_dust_init, c="black", label=r"sputtered", linewidth=linewidth, zorder=1)
    ax[1][0].plot(time_i/1e3, sputter_tot_i/mass_dust_init, c="black", linestyle="--", linewidth=linewidth-2, zorder=1)
    ax[1][0].plot(time_i/1e3, sputter_tot_hot_i/mass_dust_init, c="black", linestyle="-.", linewidth=linewidth-2, zorder=1)

    mass_cl_init = mass_cloud[0]
    mass_cloud_i = mass_cloud[time_cloud<=t]
    mass_out_cloud_i = mass_out_cloud[time_output<=t]
    time_cloud_i = time_cloud[time_cloud<=t]

    ax[0][0].plot(time_cloud_i/1e3, mass_cloud_i/mass_cl_init, linewidth=linewidth, c="#49b4ab")
    ax[0][0].plot(time_output_i/1e3, mass_out_cloud_i/mass_cl_init, linewidth=linewidth-1, linestyle="--", c="#49b4ab")

    ax[0][0].set_xticks(np.linspace(0, tmax/1e3, 5).round(0))
    ax[1][0].set_xticks(np.linspace(0, tmax/1e3, 5).round(0))
    ax[0][0].set_yticks([0.0, 0.3, 0.6, 0.9, 1.2])
    ax[1][0].set_yticks(np.linspace(0.0, 1.0, 6).round(1))
    ax[1][0].set_ylabel(r"$m_{dust}/m_{dust,i}$", labelpad=10, fontsize=fontsize-2)
    ax[0][0].set_ylabel(r"$m_{cl}/m_{cl,i}$", labelpad=10, fontsize=fontsize-2)
    ax[1][0].set_xlabel("Time [Myr]", labelpad=10, fontsize=fontsize-4)
    ax[0][0].set_ylim(0-0.1, 1.3)
    ax[1][0].set_ylim(0-0.1, 1+0.1)
    ax[1][0].set_xlim(0, tmax/1e3)
    ax[0][0].set_xlim(0, tmax/1e3)

    plt.savefig(os.path.join(pngdir, f"{i}_movie.png"), dpi=300, bbox_inches="tight")
    plt.close()

