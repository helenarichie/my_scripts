from hconfig import *
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl

from matplotlib import ticker

from matplotlib.colors import LogNorm

plt.rcParams.update({'font.size': 30})
from mpl_toolkits.axes_grid1 import make_axes_locatable

from labellines import labelLines

import seaborn as sns
from matplotlib.patches import Rectangle

hist_cmap = sns.cubehelix_palette(light=1, as_cmap=True, reverse=True)

############ hard-coded #######################################################
r_cl = 5 * 3.086e+18 # pc to cm
d_cl_init = 1e-24 # n = 1
a_grain = 1 # micron
gamma = 5/2
chi = 100
v_wind_i = 1000e3  # km/s to cm/s
tau_cc = (np.sqrt(chi)*r_cl/v_wind_i)/yr_in_s  # adiabatic cloud crushing time
T_min, T_max = 5e3, 1e8
d_min, d_max = 2e-28, 8e-24  # g cm^-3
###############################################################################

extent = np.log10([d_min, d_max, T_min, T_max])

def tau_sp_n(T, tau_sp):
    YR_IN_S = 3.154e7;
    a1 = a_grain; # dust grain size in units of 0.1 micrometers
    T_0 = 2e6; # K
    omega = 2.5;
    A = 0.17e9 * YR_IN_S; # Gyr in s

    tau_sp *= yr_in_s

    return A * 6e-4 * (a1/tau_sp) * ((T_0/T)**omega + 1)

# McKinnon et al. (2017)
def calc_tau_sp(n, T):
    YR_IN_S = 3.154e7;
    a1 = 1; # dust grain size in units of 0.1 micrometers
    d0 = n / (6e-4); # gas density in units of 10^-27 g/cm^3
    T_0 = 2e6; # K
    omega = 2.5;
    A = 0.17e9 * YR_IN_S; # Gyr in s

    tau_sp = A * (a1/d0) * ((T_0/T)**omega + 1) # s

    return tau_sp;

#print(f"{calc_tau_sp(1, 3e4)/yr_in_s:.3e}")

a = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
tau_sp = 10**a
T_sput_i = np.linspace(T_min, T_max, 1000)
T_sput = []
n_sput = []
tau_sps = []
for tau in tau_sp:
    tau_sps.append(np.linspace(tau, tau, 1000))
    T_sput.append(T_sput_i)
    n_sput.append(tau_sp_n(T_sput_i, tau))
tau_sps = np.array(tau_sps)
n_sput = np.array(n_sput)

tau_sp_sub = np.array([9.8e5])
T_sput_sub_i = np.linspace(T_min, T_max, 1000)
T_sput_sub = []
n_sput_sub = []
tau_sps_sub = []
for tau in tau_sp_sub:
    tau_sps_sub.append(np.linspace(tau, tau, 100))
    T_sput_sub.append(T_sput_sub_i)
    n_sput_sub.append(tau_sp_n(T_sput_sub_i, tau))
tau_sps_sub = np.array(tau_sps_sub[0], dtype=float)
n_sput_sub = np.array(n_sput_sub[0], dtype=float)
T_sput_sub = np.array(T_sput_sub[0], dtype=float)

"""
T_bins = np.linspace(np.log10(T_min), np.log10(T_max), 250)
T_arr = np.linspace(T_min, T_max, 1000)
n_arr = np.linspace(np.log10(d_min/(0.6*MP)), np.log10(d_max/(0.6*MP)), 1000)

# nn, TT = np.meshgrid(n_arr, T_arr)

tau_sp_wind = calc_tau_sp(nn, TT)/3.154e7
"""

d_sput = n_sput * (MP * 0.6)

fig, ax = plt.subplots(figsize=(11, 10))


ax.set_xlim(extent[0], extent[1])
ax.set_ylim(extent[2], extent[3])

plt.rcParams.update({'font.size': 23})
for j, tau in enumerate(tau_sp):
        ax.plot(np.log10(d_sput[j]), np.log10(T_sput[j]), linestyle="--", linewidth=1.5, color="black", label=r'$10^{{{:d}}}$'.format(a[j]), zorder=2)
        
labelLines(ax.get_lines(), zorder=2)

plt.rcParams.update({'font.size': 30})

#im = ax.contourf(np.log10(nn[tau_sp_wind>=9.8e5]), np.log10(TT[tau_sp_wind>=9.8e5]), tau_sp_wind[tau_sp_wind>=9.8e5], cmap=sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True), alpha=0.5, extent=extent, origin="lower")
ax.fill_between(np.log10(0.6*MP*n_sput_sub), np.log10(T_sput_sub), 9, color="lightgrey", zorder=1)
tau_sp_sub = np.array([9.8e6])
T_sput_sub_i = np.linspace(T_min, T_max, 1000)
T_sput_sub = []
n_sput_sub = []
tau_sps_sub = []
for tau in tau_sp_sub:
    tau_sps_sub.append(np.linspace(tau, tau, 100))
    T_sput_sub.append(T_sput_sub_i)
    n_sput_sub.append(tau_sp_n(T_sput_sub_i, tau))
tau_sps_sub = np.array(tau_sps_sub[0], dtype=float)
n_sput_sub = np.array(n_sput_sub[0], dtype=float)
T_sput_sub = np.array(T_sput_sub[0], dtype=float)
#ax.fill_between(np.log10(0.6*MP*n_sput_sub), np.log10(T_sput_sub), 9, color="lightgrey", zorder=0)


ax.scatter(np.log10((1)*0.6*MP), 4.48, marker="x", c="slateblue", s=100, linewidths=3.5, zorder=10)
t = ax.text(-24.33+0.12, 4.12+.48, "cloud", fontsize=25, zorder=10)
t.set_bbox(dict(facecolor='white', edgecolor="white"))
ax.scatter(np.log10((10**-2)*0.6*MP), 7.5, marker="x", c="steelblue", s=100, linewidths=3.5, zorder=10)
t = ax.text(-26.32, 7.06, "starburst\n   wind", fontsize=25, zorder=10)
t.set_bbox(dict(facecolor='white', edgecolor="white"))
#ax.add_patch(Rectangle([-26.5, 7.02], 0.519, 0.385, zorder=9, fill=True, color="white"))
#ax.add_patch(Rectangle([-26.5+0.519, 7.02], 0.509, 0.385, zorder=9, fill=True, color="lightgrey"))
ax.scatter(-27, 6.5, marker="x", c="seagreen", s=100, linewidths=3.5, zorder=10)
t = ax.text(-27.39+0.14, 6.05, "diffuse\n  wind", fontsize=25, zorder=10)
t.set_bbox(dict(facecolor='white', edgecolor="white"))

ax.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2)
ax.set_xticks([-27, -26, -25, -24])
ax.set_yticks([4, 5, 6, 7, 8])

ax.set_xlabel(r'$\mathrm{log}_{10}(\rho)$ [$\mathrm{g}\,\mathrm{cm}^{-3}$]')
ax.set_ylabel(r'$\mathrm{log}_{10}(T) \, [\mathrm{K}]$')

# add fill between for time spent in wind


plt.tight_layout()

plt.savefig("/Users/helenarichie/Desktop/sputtering_contours.png", dpi=300, bbox_inches='tight', pad_inches=0.1)  