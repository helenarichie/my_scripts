from hconfig import *
import seaborn as sns
from matplotlib import colors

date = "2023-04-03"
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
dnamein = os.path.join(basedir, "hdf5/full/")
dnameout = os.path.join(basedir, "png/phase/")
cat = True
T_min, T_max = 5e2, 1e8
n_min, n_max = 1e-5, 1.5e1

r_cl = 5 * 3.086e+18 # pc to cm
chi = 100
v_wind_i = 1000e3 # cm/s
tau_cc = (np.sqrt(chi)*r_cl/v_wind_i)/yr_in_s # cloud crushing time

def tau_sp_n(T, tau_sp):
    YR_IN_S = 3.154e7;
    a1 = 1; # dust grain size in units of 0.1 micrometers
    T_0 = 2e6; # K
    omega = 2.5;
    A = 0.17e9 * YR_IN_S; # Gyr in s

    tau_sp *= yr_in_s

    return A * 6e-4 * (a1/tau_sp) * ((T_0/T)**omega + 1)

tau_sp = np.linspace(2, 22, 20)
tau_sp = 10**tau_sp
T_sput = np.linspace(T_min, T_max, 100)
n_sput = []
for tau in tau_sp:
    n_sput.append(tau_sp_n(T_sput, tau))

if cat:
    files = glob.glob(os.path.join(dnamein, "*.h5"))
else:
    files = glob.glob(os.path.join(dnamein, "*.h5.0"))

for i in range(0, len(files)):
    data = ReadHDF5(dnamein, fnum=i, nscalar=1, cat=cat)
    head = data.head
    conserved = data.conserved

    d_gas = data.d_cgs()
    T = data.T()
    t_arr = data.t_cgs() / yr_in_s

    n_gas = d_gas / (MP * 0.6)

    log_n = np.log10(n_gas.flatten())
    log_T = np.log10(T.flatten())
    extent = np.log10([[n_min, n_max], [T_min, T_max]])

    hist_cmap = sns.cubehelix_palette(light=1, as_cmap=True, reverse=True)

    fig, ax = plt.subplots(figsize=(25,20))
    hist = plt.hist2d(log_n, log_T, 
                    bins=150, norm=colors.LogNorm(), 
                    cmap=hist_cmap,
                    range=extent)

    for j, n in enumerate(n_sput):
        plt.plot(np.log10(n), np.log10(T_sput), 
                c="k", linestyle="--", 
                linewidth=3, zorder=0, 
                alpha=0.1)

    plt.xlabel("$\log (n~[cm^{-3}])$")
    plt.ylabel("$\log (T~[K])$")
    plt.title(f"Time={round(t_arr[0]/1e6, 3)} Myr, " + r"$t/t_{cc}$=" + f"{round(t_arr[0]/tau_cc, 3)}")
    plt.colorbar(label="Number of Cells")

    plt.savefig(dnameout + f"{i}_nT.png", dpi=300)