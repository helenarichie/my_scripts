import sys
if __file__ == "/Users/helenarichie/GitHub/my_scripts/dust-survival/plot_sputtering_contours.py":
    print("hey")
sys.path.insert(0, "/Users/helenarichie/GitHub/my_scripts/")
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
a_grain = .1 # units of 0.1 micron
gamma = 5/2
chi = 100
v_wind_i = 1000e3  # km/s to cm/s
tau_cc = (np.sqrt(chi)*r_cl/v_wind_i)/yr_in_s  # adiabatic cloud crushing time
T_min, T_max = 3e3, 4e8
d_min, d_max = 1e-28, 5e-23  # g cm^-3
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

# McKinnon et al. (201
# 7)
def calc_tau_sp(n, T):
    YR_IN_S = 3.154e7;
    a1 = a_grain; # dust grain size in units of 0.1 micrometers
    d0 = n / (6e-4); # gas density in units of 10^-27 g/cm^3
    T_0 = 2e6; # K
    omega = 2.5;
    A = 0.17e9 * YR_IN_S; # Gyr in s

    tau_sp = A * (a1/d0) * ((T_0/T)**omega + 1) # s

    return tau_sp;

print(f"Base: {calc_tau_sp(2e-2, 3e7)/yr_in_s:.3e}")
print(f"5 kpc: {calc_tau_sp(1e-3, 6e6)/yr_in_s:.3e}")
print(f"10 kpc: {calc_tau_sp(2.5e-4, 3e6)/yr_in_s:.3e}")

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

tau_sp_sub = np.array([9785641.8061])
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


plt.style.use('dark_background')

fig, ax = plt.subplots(figsize=(11, 10))


ax.set_xlim(extent[0], extent[1])
ax.set_ylim(extent[2], extent[3])

plt.rcParams.update({'font.size': 23})
for j, tau in enumerate(tau_sp):
        if (a[j]%2 != 0):
            ax.plot(np.log10(d_sput[j]), np.log10(T_sput[j]), linestyle="--", linewidth=1.5, color="white", label=r'$10^{{{:d}}}$'.format(a[j]), zorder=2)
        
        else:
            ax.plot(np.log10(d_sput[j]), np.log10(T_sput[j]), linestyle="--", linewidth=1.5, color="white", zorder=2)
        
labelLines(ax.get_lines(), zorder=2)

plt.rcParams.update({'font.size': 30})

#im = ax.contourf(np.log10(nn[tau_sp_wind>=9.8e5]), np.log10(TT[tau_sp_wind>=9.8e5]), tau_sp_wind[tau_sp_wind>=9.8e5], cmap=sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True), alpha=0.5, extent=extent, origin="lower")
ax.fill_between(np.log10(0.6*MP*n_sput_sub), np.log10(T_sput_sub), 9, color="darkgrey", zorder=1)
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
#ax.fill_between(np.log10(0.6*MP*n_sput_sub), np.log10(T_sput_sub), 9, color="darkgrey", zorder=0)

n_cl = 1
log_T_cl = 4.48

##################################### Cloud ####################################
#ax.scatter(np.log10((n_cl)*0.6*MP), log_T_cl, marker="x", c="slateblue", s=100, linewidths=3.5, zorder=10)
#t = ax.text(-24.5, 4.65, "cool phase", fontsize=25, zorder=10)
#t.set_bbox(dict(facecolor='black', edgecolor="black"))
################################################################################

cgols_ns = np.array([1.5e-2, 2e-2, 9e-3, 4e-3, 2.1e-3, 1.5e-3, 1e-3, 7.5e-4, 5.5e-4, 4e-4, 3e-4, 2.5e-4])
log_cgols_Ts = np.array([7, 7.375, 7.25, 7, 6.875, 6.75, 6.7, 6.625, 6.575, 6.5, 6.4, 6.375])
# cgols = plt.scatter(np.log10(cgols_ns*0.6*MP), log_cgols_Ts, marker="x", s=100, linewidths=3.5)

distances = [0, 0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
t_sps = calc_tau_sp(cgols_ns, 10**log_cgols_Ts)

cgols_base  = [cgols_ns[1], log_cgols_Ts[1]]
cgols_2kpc  = [cgols_ns[3], log_cgols_Ts[3]]
cgols_5kpc  = [cgols_ns[6], log_cgols_Ts[6]]
cgols_10kpc = [cgols_ns[11], log_cgols_Ts[11]]

cgols_base_cool  = [2e1, 3.75]
cgols_2kpc_cool  = [5e-1, 3.9]
cgols_5kpc_cool  = [1e-1, 4]
cgols_10kpc_cool = [2e-2, 4]

##################################### hot ####################################
cgols = ax.scatter(np.log10((cgols_base[0])*0.6*MP), cgols_base[1], marker="x", c="coral", s=150, linewidths=4.5, zorder=10)
t = ax.text(-25.95, 7.1, "0.1 kpc", fontsize=25, zorder=10)
#t.set_bbox(dict(facecolor='black', edgecolor="black"))
ax.scatter(np.log10((cgols_2kpc[0])*0.6*MP), cgols_2kpc[1], marker="x", c="coral", s=150, linewidths=4.5, zorder=10)
t = ax.text(-26.65, 6.7, "2 kpc", fontsize=25, zorder=10)
#t.set_bbox(dict(facecolor='black', edgecolor="black"))
ax.scatter(np.log10((cgols_5kpc[0])*0.6*MP), cgols_5kpc[1], marker="x", c="coral", s=150, linewidths=4.5, zorder=10)
t = ax.text(-27.25, 6.35, "5 kpc", fontsize=25, zorder=10)
#t.set_bbox(dict(facecolor='black', edgecolor="black"))
ax.scatter(np.log10((cgols_10kpc[0])*0.6*MP), cgols_10kpc[1], marker="x", c="coral", s=150, linewidths=4.5, zorder=10)
t = ax.text(-27.9, 6.05, "10 kpc", fontsize=25, zorder=10)
t.set_bbox(dict(facecolor='black', edgecolor="black"))
################################################################################

##################################### hot profile ###################################
#ax.scatter(np.log10((cgols_ns[1:])*0.6*MP), log_cgols_Ts[1:], marker="x", c="coral", s=100, linewidths=3.5, zorder=10)
#ax.plot(np.log10((cgols_ns[1:])*0.6*MP), log_cgols_Ts[1:], c="coral", linewidth=3.5, zorder=10)
#t = ax.text(-28-0.03, 6.95, "r = 10 kpc", fontsize=25, zorder=10)
#t.set_bbox(dict(facecolor='black', edgecolor="black"))
################################################################################

##################################### cool ###################################
cgols_cool = ax.scatter(np.log10((cgols_base_cool[0])*0.6*MP), cgols_base_cool[1], marker="x", c="mediumslateblue", s=150, linewidths=4.5, zorder=10, label="CGOLS cool phase")
ax.scatter(np.log10((cgols_2kpc_cool[0])*0.6*MP), cgols_2kpc_cool[1], marker="x", c="mediumslateblue", s=150, linewidths=4.5, zorder=10)
ax.scatter(np.log10((cgols_5kpc_cool[0])*0.6*MP), cgols_5kpc_cool[1], marker="x", c="mediumslateblue", s=150, linewidths=4.5, zorder=10)
ax.scatter(np.log10((cgols_10kpc_cool[0])*0.6*MP), cgols_10kpc_cool[1], marker="x", c="mediumslateblue", s=150, linewidths=4.5, zorder=10)
t = ax.text(-23.1, 3.9, "0.1 kpc", fontsize=25, zorder=10)
t.set_bbox(dict(facecolor='black', edgecolor="black"))
###############################################################################


##################################### mixed ####################################
n_mix_base, T_mix_base = np.sqrt(cgols_base[0]*cgols_base_cool[0]), np.sqrt(10**cgols_base[1]*10**cgols_base_cool[1])
n_mix_2kpc, T_mix_2kpc = np.sqrt(cgols_2kpc[0]*cgols_2kpc_cool[0]), np.sqrt(10**cgols_2kpc[1]*10**cgols_2kpc_cool[1])
n_mix_5kpc, T_mix_5kpc = np.sqrt(cgols_5kpc[0]*cgols_5kpc_cool[0]), np.sqrt(10**cgols_5kpc[1]*10**cgols_5kpc_cool[1])
n_mix_10kpc, T_mix_10kpc = np.sqrt(cgols_10kpc[0]*cgols_10kpc_cool[0]), np.sqrt(10**cgols_10kpc[1]*10**cgols_10kpc_cool[1])

mixed = ax.scatter(np.log10((n_mix_base)*0.6*MP), np.log10(T_mix_base), marker="x", c="mediumseagreen", s=150, linewidths=4.5, zorder=10)
ax.scatter(np.log10((n_mix_2kpc)*0.6*MP), np.log10(T_mix_2kpc), marker="x", c="mediumseagreen", s=150, linewidths=4.5, zorder=10)
ax.scatter(np.log10((n_mix_5kpc)*0.6*MP), np.log10(T_mix_5kpc), marker="x", c="mediumseagreen", s=150, linewidths=4.5, zorder=10)
ax.scatter(np.log10((n_mix_10kpc)*0.6*MP), np.log10(T_mix_10kpc), marker="x", c="mediumseagreen", s=150, linewidths=4.5, zorder=10)
t = ax.text(-23.1, 3.9, "0.1 kpc", fontsize=25, zorder=10)
t.set_bbox(dict(facecolor='black', edgecolor="black"))
t = ax.text(-24.55, 4.05, "2 kpc", fontsize=25, zorder=10)
t.set_bbox(dict(facecolor='black', edgecolor="black"))
t = ax.text(-25.25, 4.15, "5 kpc", fontsize=25, zorder=10)
t.set_bbox(dict(facecolor='black', edgecolor="black"))
t = ax.text(-26.05, 4.15, "10 kpc", fontsize=25, zorder=10)
t.set_bbox(dict(facecolor='black', edgecolor="black"))
################################################################################

ax.add_patch(Rectangle([-23, 4.01], 0.685, 0.125, zorder=3, fill=True, color="black"))

#ax.add_patch(Rectangle([-26.5, 7.02], 0.519, 0.385, zorder=9, fill=True, color="black"))
#ax.add_patch(Rectangle([-26.5+0.519, 7.02], 0.509, 0.385, zorder=9, fill=True, color="darkgrey"))

################################################################################
#ax.scatter(-27, 6.5, marker="x", c="seagreen", s=100, linewidths=3.5, zorder=10)
#t = ax.text(-27.39+0.14, 6.05, "diffuse\n  wind", fontsize=25, zorder=10)
#t.set_bbox(dict(facecolor='white', edgecolor="white"))
################################################################################

ax.tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1, length=9, width=2)
ax.set_xticks([-28, -27, -26, -25, -24, -23])
ax.set_yticks([4, 5, 6, 7, 8])

ax.legend([cgols, cgols_cool, mixed], ["CGOLS hot phase", "CGOLS cool phase", "mixed phase"], loc="upper right", fontsize=23)
#ax.legend([cgols, cgols_cool], ["CGOLS hot phase", "CGOLS cool phase"], loc="upper right", fontsize=23)

# add fill between for time spent in wind

ax.set_xlabel(r'$\mathrm{log}_{10}(\rho)$ [$\mathrm{g}\,\mathrm{cm}^{-3}$]')
ax.set_ylabel(r'$\mathrm{log}_{10}(T) \, [\mathrm{K}]$')

plt.tight_layout()

plt.savefig("/Users/helenarichie/Desktop/sputtering_contours.png", dpi=300, bbox_inches='tight', pad_inches=0.1)  

fig, ax = plt.subplots(figsize=(11, 10))
plt.plot(distances, (t_sps/yr_in_s)*1e-6, linewidth=3.5, c="teal")

plt.xlabel(r"$r~[kpc]$")
plt.ylabel(r"$t_{sp}~[Myr]$")
plt.title(r"CGOLS Hot Phase $t_{sp}$")
plt.savefig("/Users/helenarichie/Desktop/cgols_t_sp.png", dpi=300, bbox_inches='tight', pad_inches=0.1)