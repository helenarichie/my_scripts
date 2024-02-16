from hconfig import *
from csv import writer
from matplotlib import colors
# from labellines import labelLines
import seaborn as sns

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

################# hard-coded, fill these in! ###################
date = "2024-02-06"
rho_cl_i = 1e-23  # n = 1, needed to index cloud material
cutoff = rho_cl_i*density_conversion/3 # M_sun/kpc^3
cat = True
istart = 400
iend = 400
n_hydro = 1
mu = 0.6
################################################################

####################
sputtering_contours = False
a_grain = 1  # 0.1 micrometer
n_hydro = 1
tickwidth = 2
nbins = 100
crc = True
frontier = False
#####################

hist_cmap = sns.cubehelix_palette(light=1, as_cmap=True, reverse=True)

# T_min, T_max = 10, 7e7       # survived

# p_min, p_max = 3e2, 1e7   # disr
# T_min, T_max = 800, 8e7

p_min, p_max = 3e1, 1e6
T_min, T_max = 300, 7e6

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

if crc:
    basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/" # crc
if frontier:
    basedir = f"/lustre/orion/ast181/scratch/helenarichie/{date}/"

datadir = os.path.join(basedir, "hdf5/full/")
pngdir = os.path.join(basedir, "png/phase/")

# loop through hdf5 files to write out rates and masses for cloud, gas, and dust
for i in range(istart, iend+1):
    i *= n_hydro
    # read in data
    f = h5py.File(os.path.join(datadir, str(i)+".h5"), "r")
    head = f.attrs
    dx = head["dx"][0]
    t = head["t"][0]
    gamma = head["gamma"][0]
    mu = 0.6

    d = np.array(f["density"]).flatten()
    print("d")
    temp = (mu*(gamma-1.0)*1.15831413e14)*(np.array(f["GasEnergy"]).flatten()/d)
    print("temp")
    # pressure = (np.array(f["GasEnergy"]) - ((0.5 * ((np.array(f["momentum_x"]))**2+(np.array(f["momentum_y"]))**2+(np.array(f["momentum_z"])))**2)/np.array(f["density"]))).flatten()
    pressure = np.array(f["GasEnergy"]*head["energy_unit"]*(gamma-1)).flatten()
    print("pressure")
    weights = d * head["density_unit"] * (dx*head["length_unit"])**3 * 5.02785e-34 # solar masses
    print("weights")
    extent = np.log10([[p_min, p_max], [T_min, T_max]])

    fig, ax = plt.subplots(figsize=(10,7))

    hist = plt.hist2d(np.log10(pressure/KB), np.log10(temp), weights=weights,
                      bins=nbins, norm=colors.LogNorm(), range=extent, cmap=hist_cmap)
    
    if sputtering_contours:
        for j, tau in enumerate(tau_sp):
            plt.plot(np.log10(d_sput[j]), np.log10(T_sput[j]), linestyle="--", 
                    linewidth=2, color="lightgrey", label=r'$10^{{{:d}}}$ yr'.format(a[j]))

        labelLines(plt.gca().get_lines(), zorder=2.5, fontsize=15)

    plt.xlabel(r"$\log (P/k_B)$")
    plt.ylabel("$\log (T~[K])$")
    cbar_label = "Gas Mass $M_\odot$"
    plt.title(f"{round(t/1e3, 1)} Myr", pad=10)
    cbar = plt.colorbar()
    cbar.set_label(cbar_label, rotation=270, labelpad=30)

    cbar.ax.tick_params(length=9, width=tickwidth)
    ax.tick_params(axis='both', which='both', direction='in', color='k', top=1, right=1, length=9, width=tickwidth)

    plt.tight_layout()
    plt.savefig(pngdir + f"{i}_phase_tp.png", dpi=300, bbox_inches="tight")
    
    plt.close()
