from hconfig import *
from matplotlib import colors
from labellines import labelLines
# pip install matplotlib-label-lines

##########################################
date = "2023-05-09"
cat = True
gas = True
dust = False
r_cl = 5 * 3.086e+18 # pc to cm
chi = 100
v_wind_i = 100e3 # cm/s
a_grain = 0.1 # grain radius in 0.1 microns
##########################################

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
dnamein = os.path.join(basedir, "hdf5_cloudy/full_fine/")
dnameout = os.path.join(basedir, "png/phase_cloudy/")

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
if date == "2023-05-09":
    d_min, d_max = 8e-28, 1e-22
    T_min, T_max = 1e3, 1e8
    if gas:
        vmin, vmax = 5e-5, 5e2
    if dust:
        vmin, vmax = 5e-9, 0.05

tau_cc = (np.sqrt(chi)*r_cl/v_wind_i)/yr_in_s # cloud crushing time

if cat:
    files = glob.glob(os.path.join(dnamein, "*.h5"))
else:
    files = glob.glob(os.path.join(dnamein, "*.h5.0"))

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

def ideal_gas(n, P):
    return P / (n * KB)

P_init = (10**5.5)*KB
n_ideal = np.linspace(d_min, d_max, 100) / (0.6 * MP)
T_ideal = ideal_gas(n_ideal, (10**5.5)*KB)

for i in range(0, len(files)):
    data = ReadHDF5(dnamein, fnum=i, dust=True, cat=cat)
    head = data.head
    conserved = data.conserved
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

    fig, ax = plt.subplots(figsize=(25,20))
    hist = plt.hist2d(log_d, log_T, weights=weights,
                      bins=150, norm=colors.LogNorm(),
                      range=extent, cmap=cmap, vmin=vmin, vmax=vmax)
    
    for j, tau in enumerate(tau_sp):
        plt.plot(np.log10(d_sput[j]), np.log10(T_sput[j]), linestyle="--", 
                 linewidth=3, color="grey", label=r'$10^{{{:d}}}$ yr'.format(a[j]))

    labelLines(plt.gca().get_lines(), zorder=2.5)

    p_handle = plt.plot(np.log10(n_ideal * (MP * 0.6)), np.log10(T_ideal), c="k", linewidth=3, zorder=3, label="$\log_{10}(P/k_B)" + f"={round(np.log10(P_init/KB), 2)}$")
    plt.legend(handles=p_handle)

    plt.xlabel(r"$\log (\rho_{gas}~[g\,cm^{-3}])$")

    cbar_label = None
    if gas:
        cbar_label = "Gas Mass $M_\odot$"
    if dust:
        cbar_label = "Dust Mass $M_\odot$"
    plt.ylabel("$\log (T~[K])$")
    plt.title(f"Time={round(t[0]/1e6, 3)} Myr, " + r"$t/t_{cc}$=" + f"{round(t[0]/tau_cc, 3)}")
    cbar = plt.colorbar()
    cbar.set_label(cbar_label, rotation=270, labelpad=30)

    if gas:
        plt.savefig(dnameout + f"{i}_gas_rhoT.png", dpi=300)
    if dust:
        plt.savefig(dnameout + f"{i}_dust_rhoT.png", dpi=300)
    
    plt.close()