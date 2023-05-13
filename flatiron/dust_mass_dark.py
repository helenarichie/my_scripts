from hconfig import *

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

################# hard-coded, fill these in! ###################
date = "2023-03-22"
cat = True
rho_cl_init = 1e-24  # g cm^-3
r_cl_init = 5 * 3.086e+18  # cm, initial cloud radius
chi = 1e2  # background-wind density contrast
v_b = 1e2 * 1e5  # km/s * (cm/s / km/s), background wind speed
# xlims = (0, 1)
pad = 0.005
################################################################

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/" # crc
datadir = os.path.join(basedir, "hdf5/full/")
csvdir = os.path.join(basedir, "csv")
pngdir = os.path.join(basedir, "png/flatiron/")

data = ReadHDF5(datadir, fnum=0, cat=cat)
head = data.head
nx, ny, nz = head["dims"]
dx = head["dx"][0] # M_sun

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

for i, t in enumerate(t_arr):

    plt.style.use('dark_background')
    plt.rcParams.update({'font.family': 'Helvetica'})
    plt.rcParams.update({'font.size': 30})

    t_arr_i = t_arr[t_arr<=t]
    mass_dust_tot_i = mass_dust_tot[t_arr<=t]
    mass_out_i = mass_out[t_arr<=t]

    fig, ax = plt.subplots(1, 2, figsize=(25,10))
    ax[0].plot(t_arr_i/1e6, mass_dust_tot_i, linewidth=4)
    ax[0].set_ylim(ymin, ymax)
    ax[0].set_xlim(xmin/1e6, xmax/1e6)
    ax[0].set_title("Total Mass")
    ax[0].set_xlabel("Time [Myr]")
    ax[0].set_ylabel(r"Dust Mass$~[M_\odot]$")
    ax[0].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1, length=9, width=2)

    ax[1].plot(t_arr_i/1e6, mass_out_i, linewidth=4)
    ax[1].set_ylim(ymin, ymax)
    ax[1].set_xlim(xmin/1e6, xmax/1e6)
    ax[1].set_title("Mass Out")
    ax[1].set_xlabel("Time [Myr]")
    ax[1].set_ylabel(r"Dust Mass$~[M_\odot]$")
    ax[1].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1, length=9, width=2)
   
    plt.tight_layout()

    print(f"Saving figures {i} of {len(t_arr)-1}.\n")
    
    plt.savefig(pngdir + f"{i}_dust_mass.png", dpi=300)

    plt.close()