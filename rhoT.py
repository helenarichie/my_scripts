from hconfig import *
from matplotlib import colors
from labellines import labelLines
# pip install matplotlib-label-lines

date = "2023-03-07"
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
dnamein = os.path.join(basedir, "hdf5/full/")
dnameout = os.path.join(basedir, "png/phase/")
cat = False
gas = True
dust = False

T_min, T_max = 5e2, 5e7
d_min, d_max = None, None
if gas:
    d_min, d_max = 3e-29, 1.5e-23
if dust:
    d_min, d_max = 5e-31, 1.5e-25

r_cl = 5 * 3.086e+18 # pc to cm
chi = 1000
v_wind_i = 1000e3 # cm/s
tau_cc = (np.sqrt(chi)*r_cl/v_wind_i)/yr_in_s # cloud crushing time

if cat:
    files = glob.glob(os.path.join(dnamein, "*.h5"))
else:
    files = glob.glob(os.path.join(dnamein, "*.h5.0"))

def tau_sp_n(T, tau_sp):
    YR_IN_S = 3.154e7;
    a1 = 1; # dust grain size in units of 0.1 micrometers
    T_0 = 2e6; # K
    omega = 2.5;
    A = 0.17e9 * YR_IN_S; # Gyr in s

    tau_sp *= yr_in_s

    return A * 6e-4 * (a1/tau_sp) * ((T_0/T)**omega + 1)

a = np.arange(0, 20, 2)
tau_sp = 10**a
T_sput_i = np.linspace(T_min, T_max, 100)
T_sput = []
n_sput = []
tau_sps = []
for tau in tau_sp:
    tau_sps.append(np.linspace(tau, tau, 100))
    T_sput.append(T_sput_i)
    n_sput.append(tau_sp_n(T_sput_i, tau))
tau_sps = np.array(tau_sps)
n_sput = np.array(n_sput)

for i in range(0, len(files)):
    data = ReadHDF5(dnamein, fnum=i, nscalar=1, cat=cat)
    head = data.head
    conserved = data.conserved
    dx = data.dx_cgs()[0]

    d = None
    cmap = None
    if gas:
        d = data.d_cgs()
        cmap = "viridis"
    if dust:
        d = conserved["scalar0"] * head["density_unit"]
        wh_zero = np.where(d<=0)
        d[wh_zero] = 1e-40
        cmap = "plasma"

    T = data.T()
    t = data.t_cgs() / yr_in_s

    d_sput = n_sput * (MP * 0.6)
    weights = d.flatten() * dx**3 * 5.02785e-34 # solar masses

    log_d = np.log10(d.flatten())
    log_T = np.log10(T.flatten())
    extent = np.log10([[d_min, d_max], [T_min, T_max]])

    fig, ax = plt.subplots(figsize=(25,20))
    hist = plt.hist2d(log_d, log_T, weights=weights,
                      bins=150, norm=colors.LogNorm(),
                      range=extent, cmap=cmap)
    
    for j, tau in enumerate(tau_sp):
        plt.plot(np.log10(d_sput[j]), np.log10(T_sput[j]), linestyle="--", 
                 linewidth=3, color="grey", label=r'$10^{{{:d}}}$ yr'.format(a[j]))

    labelLines(plt.gca().get_lines(), zorder=2.5)

    if gas:
        plt.xlabel(r"$\log (\rho_{gas}~[g\,cm^{-3}])$")
    if dust:
        plt.xlabel(r"$\log (\rho_{dust}~[g\,cm^{-3}])$")
    plt.ylabel("$\log (T~[K])$")
    plt.title(f"Time={round(t[0]/1e6, 3)} Myr, " + r"$t/t_{cc}$=" + f"{round(t[0]/tau_cc, 3)}")
    cbar = plt.colorbar()
    cbar.set_label("$M_\odot$", rotation=270, labelpad=30)

    if gas:
        plt.savefig(dnameout + f"{i}_gas_rhoT.png", dpi=300)
    if dust:
        plt.savefig(dnameout + f"{i}_dust_rhoT.png", dpi=300)