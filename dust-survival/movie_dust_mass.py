import sys
sys.path.insert(0, "/ix/eschneider/helena/code/my_scripts")
from hconfig import *

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': 20})

##################################################################
date = "2023-04-26"
datestr = "0426"
cat = True
pad = 0.1
fontsize = 20
labelpad = 12
tickwidth = 1
tmax_1105 = 52.1e6
tmax_1103 = 52.1e6
tmax_0426 = 2.4e6
##################################################################

tmax = tmax_0426

##################################################################
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
projdir = os.path.join(basedir, "hdf5/proj/")
pngdir = os.path.join(basedir, "png/dust_mass/")
csvdir = os.path.join(basedir, "csv/")
##################################################################

# get value of dx to convert from density to mass
data = ReadHDF5(projdir, proj="xy", cat=cat, fnum=0)
head = data.head
conserved = data.conserved
dx = head["dx"][0]

# get times of snapshots
t_arr = [] # yr
with open(os.path.join(csvdir, "t_arr.csv")) as f:
    for line in f:
        line = line.split(",")
        t = line[0].strip("\n").strip("[").strip("]")
        t_arr.append(float(t))
t_arr = np.array(t_arr)

n_steps = len(t_arr)

# get timestep
dt = (t_arr[1] - t_arr[0])/1e3

# get outflow rates for dust
rate_dust = np.zeros((n_steps, 6))  # M_sun / yr
with open(os.path.join(csvdir, "rate_dust.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_dust[i] = np.array(line, dtype=float)

# get total dust masses in simulation volume
mass_dust_tot = []
with open(os.path.join(csvdir, "rho_d_tot.csv")) as f:
    for line in f:
        line = line.split(",")
        density = float(line[0]) * density_conversion
        mass = density * dx**3
        mass_dust_tot.append(mass)
mass_dust_tot = np.array(mass_dust_tot)

mass_dust_init = mass_dust_tot[0]
print(f"initial mass: {mass_dust_init:e}")

# use outflow rates to calculate how much dust has exited the simulation volume
mass_i = 0
mass_out = []
mass_lost = []
for i, rate in enumerate(rate_dust):
    rate_i = np.sum(rate)
    mass_lost.append(rate_i * dt)
    mass_i += rate_i * dt
    mass_out.append(mass_i)
mass_out = np.array(mass_out)
mass_lost = np.array(mass_lost)

ymin = np.amin([np.amin(mass_dust_tot), np.amin(mass_out)]) - pad
ymax = np.amax([np.amax(mass_dust_tot), np.amax(mass_out)]) + pad
xmin = np.amin(t_arr) - pad
xmax = np.amax(t_arr) + pad

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(7, 7))

# calculate how much dust was destroyed due to sputtering
mass_destroyed = np.zeros(len(mass_lost))
mass_destroyed[0] = 0
mass_dest_cum = 0
for i, mass in enumerate(mass_dust_tot):
    if (i+1 < len(mass_lost)):
        mass_dest_cum += (mass - mass_dust_tot[i+1] - mass_lost[i+1])
        mass_destroyed[i+1] = mass_dest_cum

for i, t in enumerate(t_arr):   
    plt.style.use('dark_background')

    t_arr_i = t_arr[t_arr<=t]
    mass_dust_tot_i = mass_dust_tot[t_arr<=t]
    mass_out_i = mass_out[t_arr<=t]
    mass_destroyed_i = mass_destroyed[t_arr<=t]

    axs.plot(t_arr_i/1e6, mass_dust_tot_i/mass_dust_init, linewidth=4, c="#d43a4f", label="in box")
    axs.plot(t_arr_i/1e6, mass_out_i/mass_dust_init, linewidth=4, linestyle="--", c="#d43a4f", label="exited box")
    axs.plot(t_arr_i/1e6, mass_destroyed_i/mass_dust_init, linewidth=4, c="white", zorder=10, label="sputtered")

    axs.set_xlabel("Time [Myr]", fontsize=fontsize)
    axs.set_ylabel(r"m_{dust}/m_{dust,i}$", fontsize=fontsize)

    axs.set_xlim(xmin/1e6, xmax/1e6)
    axs.tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1, length=9, width=2, reset=True)
    axs.set_xticks(np.linspace(0, np.amax(t_arr/1e6), 5).round(1))
    #axs.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    #axs.set_yticks(np.linspace(0, np.amax(mass_dust_tot[t_arr<=tmax]), 5).round(2))
    axs.legend(loc="upper right")    

    fig.tight_layout()
    
    print(f"Saving figure {i} of {len(t_arr)-1}.\n")
    
    plt.savefig(pngdir + f"{i}_dust_mass.png", dpi=300)

    plt.close()
