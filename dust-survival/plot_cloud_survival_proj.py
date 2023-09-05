import sys
sys.path.insert(0, "/ix/eschneider/helena/code/my_scripts")
from hconfig import *

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

##################################################################
date = "2023-04-26"
save = True
cat = True
dust = False
vlims_gas = (18.2, 19.8)
vlims_dust = (-13.6, -5)
nrow = 2
ncol = 4
fnums = [0, 71, 285]
fnum = 0
pad = 0.1
spacing = 10  # pc, spacing of tick marks
fontsize = 20
labelpad = 12
tickwidth = 2
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': 25})
##################################################################

##################################################################
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
projdir = os.path.join(basedir, "hdf5/proj/")
pngdir = os.path.join(basedir, "/ix/eschneider/helena/figs/dust_survival/")
csvdir = os.path.join(basedir, "csv/")
##################################################################

data = ReadHDF5(projdir, dust=dust, proj="xy", cat=cat, fnum=0)
head = data.head
conserved = data.conserved
dx = head["dx"][0]

t_arr = [] # yr
with open(os.path.join(csvdir, "t_arr.csv")) as f:
    for line in f:
        line = line.split(",")
        t = line[0].strip("\n").strip("[").strip("]")
        t_arr.append(float(t))
t_arr = np.array(t_arr)

n_steps = len(t_arr)

dt = (t_arr[1] - t_arr[0])/1e3

rate_dust = np.zeros((n_steps, 6))  # M_sun / yr
with open(os.path.join(csvdir, "rate_dust.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_dust[i] = np.array(line, dtype=float)

mass_dust_tot = []
with open(os.path.join(csvdir, "rho_d_tot.csv")) as f:
    for line in f:
        line = line.split(",")
        density = float(line[0]) * density_conversion
        mass = density * dx**3
        mass_dust_tot.append(mass)
mass_dust_tot = np.array(mass_dust_tot)

mass_i = 0
mass_out = []
for i, rate in enumerate(rate_dust):
    rate_i = np.sum(rate)
    mass_i += rate_i * dt
    mass_out.append(mass_i)
mass_out = np.array(mass_out)

ymin = np.amin([np.amin(mass_dust_tot), np.amin(mass_out)]) - pad
ymax = np.amax([np.amax(mass_dust_tot), np.amax(mass_out)]) + pad
xmin = np.amin(t_arr) - pad
xmax = np.amax(t_arr) + pad

dx = None
fig, axs = plt.subplots(nrows=nrow, ncols=ncol, figsize=(20, 10))

for i, fnum in enumerate(fnums):

    data = ReadHDF5(projdir, dust=dust, proj="xy", cat=cat, fnum=fnum)
    head = data.head
    conserved = data.conserved

    nx, ny, nz = head["dims"]
    dx = head["dx"][0] * 1e3  # pc
    d_gas = conserved["density"]
    d_gas *= head["mass_unit"] / (head["length_unit"] ** 2)

    gamma = head["gamma"]
    t = data.t_cgs() / yr_in_s  # yr

    d_dust = conserved["dust_density"]

    d_dust *= head["mass_unit"] / (head["length_unit"] ** 2)

    extent = [0, nz*dx, 0, nz*dx]

    # xy gas density projection
    n_gas = d_gas[0]/(0.6*MP) # column density
    im = axs[0][i].imshow(np.log10(n_gas[0:nz,:].T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=extent)
    im = axs[1][i].imshow(np.log10(n_gas[0:nz,:].T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=extent)
    #ylabel = r'$\mathrm{log}_{10}(N_{H, gas})$ [$\mathrm{cm}^{-2}$]'
    #divider = make_axes_locatable(axs[0])
    #cax = divider.append_axes("right", size="5%", pad=pad)
    #cbar = fig.colorbar(im, ax=axs[0], cax=cax, pad=pad)
    #cbar.ax.tick_params(length=9, width=tickwidth)
    #cbar.set_label(ylabel, fontsize=fontsize, labelpad=11)
    #cbar.set_ticks(np.linspace(vlims_gas[0], vlims_gas[1], 4).round(1))
    for j in range(0, ncol):
        for k in range(0, nrow):
            if ((j == 0) and (k == 0)) or ((j == 0) and (k == 1)):
                axs[k][j].hlines(0.12*dx*ny, spacing, spacing+spacing, color='white')
                axs[k][j].text(spacing+spacing+2, 0.1*dx*ny, '10 pc', color='white', fontsize=fontsize)
            if not (((j == 0) and (k == 3)) or ((j == 1) and (k == 3))):
                axs[k][j].set_xticks(np.arange(0, ny*dx, spacing))
                axs[k][j].set_yticks(np.arange(0, ny*dx, spacing))
                axs[k][j].tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=9, width=tickwidth)

axs[0][3].plot(t_arr/1e6, mass_dust_tot, linewidth=4, c="k")
#axs[0][3].set_ylim(ymin, ymax)
axs[0][3].set_xlim(xmin/1e6, xmax/1e6)
#axs[0][3].set_title("Total Mass")
axs[0][3].set_xlabel("Time [Myr]")
axs[0][3].set_ylabel(r"Dust Mass$~[M_\odot]$")
axs[0][3].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True, labelsize="small")

axs[0][3].plot(t_arr/1e6, mass_out, linewidth=4, linestyle="--", c="k")
#axs[1][3].set_ylim(ymin, ymax)
axs[0][3].set_xlim(xmin/1e6, xmax/1e6)
axs[0][3].set_xlabel("Time [Myr]", fontsize=fontsize)
axs[0][3].set_ylabel(r"Dust Mass$~[M_\odot]$", fontsize=fontsize)
axs[0][3].set_xticks(np.linspace(0, np.amax(t_arr/1e6), 5).round(1))
axs[0][3].set_yticks(np.linspace(0, np.amax(mass_dust_tot), 5).round(2))

fig.tight_layout()

save = True
if save:
    plt.savefig(pngdir + f"cloud_survival.png", dpi=300)
plt.close()
