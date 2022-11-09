from hconfig import *


######### hard-coded values! #########
date = "2022-10-10"
datadir = f"/ix/eschneider/helena/data/cloud_wind/{date}/hdf5/r100pc/full-vol/"
r_cl = 100 * 3.086e+18 # pc to cm
d_cl_init = 1e-24 # n = 1
# cloud_exit = 3.5e6 # yr
cloud_exit = 7e6
r_grain = 0.1 # micron
######################################

data = ReadHDF5(datadir, nscalar=1)
head = data.head
conserved = data.conserved

start = time.time()
print("\nloading gas densities...")
d_gas = data.d_cgs()
print(f"Finished --- {round(time.time() - start, 3)} s")

start = time.time()
print("\nloading dust densities...")
d_dust = conserved["scalar0"] * head["density_unit"]
print(f"Finished --- {round(time.time() - start, 3)} s")

start = time.time()
print("\nloading velocities...")
vx = data.vx_cgs()
print(f"Finished --- {round(time.time() - start, 3)} s")
gamma = head["gamma"]
t_arr = data.t_cgs() / yr_in_s


start = time.time()
print("\nloading temperatures...")
T = data.T()
print(f"Finished --- {round(time.time() - start, 3)} s")
T_hot = T[0][0][0][0] # initial wind temperature

dust_avg = []
dust_tot = []
cloud_tot = []
cloud_avg = []
cloud_avg_T = []

start = time.time()
print("\nloop starting...")
for i, dens in enumerate(d_dust):
    print(i)
    dust_tot.append(np.sum(d_dust[i][d_dust[i]>0]))
    dust_avg.append(np.average(d_dust[i][d_dust[i]>0]))
    cloud_tot.append(np.sum(d_gas[i][d_gas[i]>=1/3*d_cl_init]))
    cloud_avg.append(np.average(d_gas[i][d_gas[i]>=1/3*d_cl_init]))
    cloud_avg_T.append(np.average(T[i][d_gas[i]>=1/3*d_cl_init]))

print(f"Finished --- {round(time.time() - start, 3)} s")

dust_avg = np.array(dust_avg)
dust_tot = np.array(dust_tot)
cloud_avg = np.array(cloud_avg)
cloud_tot = np.array(cloud_tot)
cloud_avg_T = np.array(cloud_avg_T)

np.savetxt("dust_tot.csv", dust_tot, delimiter=",")
np.savetxt("cloud_tot.csv", cloud_tot, delimiter=",")

m_init = dust_tot[0]
m_cl_init = cloud_tot[0]

d_dust_i = dust_avg[0] # initial dust density
d_gas_cl_i = cloud_avg[0] # intial cloud gas mass density
d_gas_bg_i = d_gas[0][0][0][0] # intial wind mass density
T_cool = cloud_avg_T[0] # initial cloud temperature
vx_bg_i = vx[0][0][0][0] # initial background x-velocity

dtg_cl = d_dust_i / d_gas_cl_i
n_gas_cl_i = d_gas_cl_i / (0.6 * 1.672622e-24) # intial cloud gas number density
n_gas_bg_i = d_gas_bg_i / (0.6 * 1.672622e-24) # intial wind number density
chi = d_gas_cl_i / d_gas_bg_i
P0 = n_gas_bg_i*1.4e-16*T_hot # nkT, wind pressure
cs = np.sqrt(gamma*P0/d_gas_bg_i) # speed of sound in wind
mach = vx_bg_i / cs # mach number of the wind
t_cc = chi**(1/2)*r_cl/vx_bg_i # cloud crushing time

T_mix = np.sqrt(T_cool * T_hot)
n_mix = np.sqrt(n_gas_cl_i * n_gas_bg_i)

t_sp = 0.17*3.1536e16*(10*r_grain/(n_mix/6e-4))*(2e6/T_mix)**2.5 # sputtering time

others = np.array([d_gas_bg_i, d_gas_cl_i, dtg_cl, chi, T_hot, T_cool, mach, t_cc, t_sp])
np.savetxt("others.csv", others, delimiter=",")

print("plotting...")

fig = plt.figure(figsize=(15,12))

gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], width_ratios=[2, 1])
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[:, 1])

ax1.semilogx(t_arr, cloud_tot/m_cl_init, linewidth=3, color=colors[5])
ax1.set_ylabel(r'$m_g/m_{g,init}$', fontsize=14)
ax1.set_xlabel(r'Time [yr]', fontsize=14)
ax1.set_title(r"Fraction of Gas Mass Remaining")
####### hard-coded ############
#ax2.set_ylim(0, 1.01)
###############################
ax1.vlines(cloud_exit, 0, 1, transform=ax1.get_xaxis_transform(), colors="k", alpha=0.5, label="cloud exiting simulation volume")
ax1.legend()

ax2.semilogx(t_arr, dust_tot/m_init, linewidth=3, color=colors[5])
ax2.set_ylabel(r'$m_d/m_{d,init}$', fontsize=14)
ax2.set_xlabel(r'Time [yr]', fontsize=14)
ax2.set_title(r"Fraction of Dust Mass Remaining")
####### hard-coded ############
ax2.set_ylim(0.8, 1.01)
###############################
ax2.vlines(cloud_exit, 0, 1, transform=ax2.get_xaxis_transform(), colors="k", alpha=0.5, label="cloud exiting simulation volume")
ax2.legend()


ax3.axis("off")
ax3.text(0, 1.0, f"Initial background density: {d_gas_bg_i:.2e}" + "$~cm^{-3}$", fontsize=16)
ax3.text(0, 0.9, f"Initial cloud density: {d_gas_cl_i:.2e}" + "$~cm^{-3}$", fontsize=16)
ax3.text(0, 0.8, f"Cloud dust-to-gas ratio: {round(dtg_cl, 4)}", fontsize=16)
ax3.text(0, 0.7, f"$\chi$: {round(chi, 0)}", fontsize=16)
ax3.text(0, 0.6, "$T_{hot}$: " + f"{T_hot:.2e}$~K$", fontsize=16)
ax3.text(0, 0.5, "$T_{cool}$: " + f"{T_cool:.2e}$~K$", fontsize=16)
ax3.text(0, 0.4, f"Mach number: {mach:.2}", fontsize=16)
ax3.text(0, 0.3, "$t_{cc}$: " + f"{t_cc/3.154e+7:.2e}$~yr$", fontsize=16)
ax3.text(0, 0.2, "$t_{sp}$: " + f"{t_sp/3.154e+7:.2e}$~yr$", fontsize=16)
ax3.text(0, 0.1, "Cloud radius: " + f"{round(r_cl/3.086e+18, 0)}$~pc$", fontsize=16)
ax3.text(0, 0.0, "Dust grain size: " + f"{r_grain:.2}$~\mu m$", fontsize=16)


plt.savefig(os.path.join(datadir) + "cloud_mass_" + date + ".png", dpi=300)
plt.savefig("cloud_mass_" + date + "medium-long.png", dpi=300)


print("Plotting completed.\n")


print(f"Dust to gas ratio: {dtg_cl}")
print(f"Initial gas density: {d_gas_cl_i}")
print(f"Chi: {chi}")
print(f"Mach number: {mach}")
print(f"Hot phase temperature: {T_hot}")

