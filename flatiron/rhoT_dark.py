import sys
sys.path.insert(0, "/ix/eschneider/helena/code/my_scripts/")
from hconfig import *
from csv import writer
from matplotlib import colors
from labellines import labelLines
import seaborn as sns

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

################# hard-coded, fill these in! ###################
date = "2024-02-06"
rho_cl_i = 1e-23  # n = 1, needed to index cloud material
cutoff = rho_cl_i*density_conversion/3 # M_sun/kpc^3
cat = True
istart = 200
iend = 200
n_hydro = 1
################################################################

####################
gas = True
dust = False
a_grain = 1  # 0.1 micrometer
n_hydro = 1
tickwidth = 2
#####################

plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': 27})
# plt.style.use('dark_background')

hist_cmap = sns.cubehelix_palette(light=1, as_cmap=True, reverse=True)

d_min, d_max = None, None
vmin, vmax = None, None
if date == "2023-04-26":
    d_min, d_max = 8e-28, 1e-22
    T_min, T_max = 1e3, 1e8
    if gas:
        vmin, vmax = 5e-5, 5e2
    if dust:
        vmin, vmax = 5e-9, 0.05
if date == "2023-05-03":
    d_min, d_max = 8e-28, 1e-22
    T_min, T_max = 1e3, 1e8
    if gas:
        vmin, vmax = 5e-5, 5e2
    if dust:
        vmin, vmax = 5e-9, 0.05
if date == "2023-03-22":
    d_min, d_max = 5e-28, 1e-22
    T_min, T_max = 1e3, 7e7
    if gas:
        vmin, vmax = 5e-5, 5e2
    if dust:
        vmin, vmax = 5e-9, 0.05
if date == "2023-05-12":
    d_min, d_max = 5e-28, 1e-22
    T_min, T_max = 1e3, 5e7
    if gas:
        vmin, vmax = 5e-5, 5e2
    if dust:
        vmin, vmax = 5e-9, 0.05
if date == "2023-05-13":
    d_min, d_max = 5e-28, 1e-22
    T_min, T_max = 1e3, 7e7
    if gas:
        vmin, vmax = 5e-5, 5e2
    if dust:
        vmin, vmax = 5e-9, 0.05
if date == "2024-02-06":
    d_min, d_max = 5e-28, 1e-22
    T_min, T_max = 10, 7e7
    if gas:
        vmin, vmax = 5e-5, 5e2
    if dust:
        vmin, vmax = 5e-9, 0.05

def tau_sp_n(T, tau_sp):
    YR_IN_S = 3.154e7;
    a1 = a_grain; # dust grain size in units of 0.1 micrometers
    T_0 = 2e6; # K
    omega = 2.5;
    A = 0.17e9 * YR_IN_S; # Gyr in s

    tau_sp *= yr_in_s

    return A * 6e-4 * (a1/tau_sp) * ((T_0/T)**omega + 1)

a = np.arange(0, 20, 2)
a = np.array([0, 2, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18])
tau_sp = 10**a
T_sput_i = np.linspace(T_min, T_max, 1000)
T_sput = []
n_sput = []
tau_sps = []
for tau in tau_sp:
    tau_sps.append(np.linspace(tau, tau, 100))
    T_sput.append(T_sput_i)
    n_sput.append(tau_sp_n(T_sput_i, tau))
tau_sps = np.array(tau_sps)
n_sput = np.array(n_sput)

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/" # crc
datadir = os.path.join(basedir, "hdf5/full/")
csvdir = os.path.join(basedir, "csv/")
pngdir = os.path.join(basedir, "png/flatiron/")

data = ReadHDF5(datadir, fnum=200, dust=True, cat=cat)
head = data.head
nx, ny, nz = head["dims"]
dx = head["dx"] # kpc

# x, y, z
areas = [ny*nz*dx**2, nx*nz*dx**2, nx*ny*dx**2]
indices = [0, -1]
fluids = ["gas", "dust"]

if cat:
    files = glob.glob(os.path.join(datadir, "*.h5"))
else:
    files = glob.glob(os.path.join(datadir, "*.h5.0"))

if istart == 0:
    # create files to write quantities for cloud (dense gas only), total gas, and dust
    f = open(os.path.join(csvdir, "rate_cl.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "mass_cl.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "rate_gas.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "mass_gas.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "mass_dust.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "rate_dust.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "rho_d_tot.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "rho_cl_tot.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "T_cl_avg.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "t_arr.csv"), "w")
    f.close()

def calc_mass_loss_rate(rho, v, area):
    return np.sum(rho * np.abs(v) * area)

def get_rates(d, cutoff):
    rates = np.zeros(6)
    for i, area in enumerate(areas):
        for index in indices:
            if (index == 0) and (i == 0):
                d_i, v_i = d[0, :, :], conserved["momentum_x"][0][0, :, :] / conserved["density"][0][0, :, :]
                mask = np.logical_and(d_i >= cutoff, v_i < 0) # mask for high-density, outflowing gas
                rates[0] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == 0) and (i == 1):
                d_i, v_i = d[:, 0, :], conserved["momentum_y"][0][:, 0, :] / conserved["density"][0][:, 0, :]
                mask = np.logical_and(d_i >= cutoff, v_i < 0)
                rates[1] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == 0) and (i == 2):
                d_i, v_i = d[:, :, 0], conserved["momentum_z"][0][:, :, 0] / conserved["density"][0][:, :, 0]
                mask = np.logical_and(d_i >= cutoff, v_i < 0)
                rates[2] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == -1) and (i == 0):
                d_i, v_i = d[-1, :, :], conserved["momentum_x"][0][-1, :, :] / (conserved["density"][0][-1, :, :])
                mask = np.logical_and(d_i >= cutoff, v_i > 0)
                rates[3] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == -1) and (i == 1):
                d_i, v_i = d[:, -1, :], conserved["momentum_y"][0][:, -1, :] / conserved["density"][0][:, -1, :]
                mask = np.logical_and(d_i >= cutoff, v_i > 0)
                rates[4] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == -1) and (i == 2):
                d_i, v_i = d[:, :, -1], conserved["momentum_z"][0][:, :, -1] / conserved["density"][0][:, :, -1]
                mask = np.logical_and(d_i >= cutoff, v_i > 0)
                rates[5] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)
    
    return rates

def get_masses(d, cutoff):
    masses = np.zeros(6)
    for i, area in enumerate(areas):
        for index in indices:
            if (index == 0) and (i == 0):
                d_i = d[0, :, :]
                mask = (d_i >= cutoff)
                masses[0] = np.sum(d_i[mask]) * dx**3 # total density times total volume of cells

            if (index == 0) and (i == 1):
                d_i = d[:, 0, :]
                mask = (d_i >= cutoff)
                masses[1] = np.sum(d_i[mask]) * dx**3

            if (index == 0) and (i == 2):
                d_i = d[:, :, 0]
                mask = (d_i >= cutoff)
                masses[2] = np.sum(d_i[mask]) * dx**3

            if (index == -1) and (i == 0):
                d_i = d[-1, :, :]
                mask = (d_i >= cutoff)
                masses[3] = np.sum(d_i[mask]) * dx**3

            if (index == -1) and (i == 1):
                d_i = d[:, -1, :]
                mask = (d_i >= cutoff)
                masses[4] = np.sum(d_i[mask]) * dx**3

            if (index == -1) and (i == 2):
                d_i = d[:, :, -1]
                mask = (d_i >= cutoff)
                masses[5] = np.sum(d_i[mask]) * dx**3
    
    return masses

def write_csv(data):
    head = data.head
    conserved = data.conserved

    d_gas = data.d_cgs()
    d_dust = conserved["dust_density"] * head["density_unit"]
    T = data.T()
    t = data.t_cgs() / yr_in_s

    rho_d_tot_i = None
    rho_cl_tot_i = None
    T_cl_avg_i = None
    for i, dens in enumerate(d_dust):
        rho_d_tot_i = np.sum(d_dust[i][d_dust[i]>0])
        rho_cl_tot_i = np.sum(d_gas[i][d_gas[i]>=1/3*rho_cl_i])
        T_cl_avg_i = np.average(T[i][d_gas[i]>=1/3*rho_cl_i])

    with open(os.path.join(csvdir, "rho_d_tot.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow([rho_d_tot_i])
        f.close()

    with open(os.path.join(csvdir, "rho_cl_tot.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow([rho_cl_tot_i])
        f.close()

    with open(os.path.join(csvdir, "T_cl_avg.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow([T_cl_avg_i])
        f.close()

    with open(os.path.join(csvdir, "t_arr.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow([t[0]])
        f.close()

# loop through hdf5 files to write out rates and masses for cloud, gas, and dust
for i in range(istart, iend+1):
    i *= n_hydro
    # read in data
    data = ReadHDF5(datadir, fnum=i, dust=True, cat=cat)
    head = data.head
    conserved = data.conserved
    dx = head["dx"][0]

    """
    write_csv(data)

    # calculate and write rates and masses for cloud
    rates = get_rates(conserved["density"][0], cutoff)
    masses = get_masses(conserved["density"][0], cutoff)

    with open(os.path.join(csvdir, "rate_cl.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(rates)
        f.close()

    with open(os.path.join(csvdir, "mass_cl.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(masses)
        f.close()

    # calculate and write rates and masses for gas
    rates = get_rates(conserved["density"][0], 0)
    masses = get_masses(conserved["density"][0], 0)

    with open(os.path.join(csvdir, "rate_gas.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(rates)
        f.close()

    with open(os.path.join(csvdir, "mass_gas.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(masses)
        f.close()

    # calculate and write masses for dust
    masses = get_masses(conserved["dust_density"][0], 1)

    with open(os.path.join(csvdir, "mass_dust.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(masses)
        f.close()

    # calculate and write masses for dust
    rates = get_rates(conserved["dust_density"][0], 1)

    with open(os.path.join(csvdir, "rate_dust.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(rates)
        f.close()
    """

    dx = data.dx_cgs()[0]

    d = None
    d_weight = None
    cmap = None
    if gas:
        d = data.d_cgs()
        d_weight = d
        cmap = "viridis"
    if dust:
        d = data.d_cgs()
        d_weight = conserved["dust_density"] * head["density_unit"]
        wh_zero = np.where(d_weight<=0)
        d_weight[wh_zero] = 1e-40
        cmap = "plasma"

    T = data.T()
    t = data.t_cgs() / yr_in_s

    d_sput = n_sput * (MP * 0.6)
    weights = d_weight.flatten() * dx**3 * 5.02785e-34 # solar masses

    log_d = np.log10(d.flatten())
    log_T = np.log10(T.flatten())
    extent = np.log10([[d_min, d_max], [T_min, T_max]])

    fig, ax = plt.subplots(figsize=(20,15))
    hist = plt.hist2d(log_d, log_T, weights=weights,
                      bins=150, norm=colors.LogNorm(),
                      range=extent, cmap=hist_cmap, vmin=vmin, vmax=vmax)
    
    for j, tau in enumerate(tau_sp):
        plt.plot(np.log10(d_sput[j]), np.log10(T_sput[j]), linestyle="--", 
                 linewidth=3, color="white", label=r'$10^{{{:d}}}$ yr'.format(a[j]))

    labelLines(plt.gca().get_lines(), zorder=2.5)

    plt.xlabel(r"$\log (\rho_{gas}~[g\,cm^{-3}])$", fontsize=30)

    cbar_label = None
    if gas:
        cbar_label = "Gas Mass $M_\odot$"
    if dust:
        cbar_label = "Dust Mass $M_\odot$"
    plt.ylabel("$\log (T~[K])$", fontsize=30)
    plt.title(f"{round(t[0]/1e6, 1)} Myr", fontsize=32, pad=10)
    cbar = plt.colorbar()
    cbar.set_label(cbar_label, rotation=270, labelpad=30)

    cbar.ax.tick_params(length=9, width=tickwidth)
    ax.tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1, length=9, width=tickwidth)

    if gas:
        plt.savefig(pngdir + f"{i}_gas_rhoT.png", dpi=300)
    if dust:
        plt.savefig(pngdir + f"{i}_dust_rhoT.png", dpi=300)
    
    plt.close()