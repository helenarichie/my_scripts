from hconfig import *

date = "2024-01-18"

basedir = f"/ix/eschneider/helena/data/testing/{date}/"
pngdir = os.path.join(basedir, "png/")

cloud_velocities = []
cloud_masses = []
times = []
i = 0
with open(os.path.join(basedir, "output.log")) as f:
    for line in f:
        if line.startswith("Average cloud velocity"):
            line = line.lstrip("Average cloud velocity = ").rstrip(" km/s\n")
            cloud_velocities.append(float(line))
        if line.startswith("Mass"):
            line = line.lstrip("Mass = ").rstrip(" M_sun\n")
            cloud_masses.append(float(line))
        if line.startswith("n_step"):
            line = line.split(" ")
            line = line[1].lstrip(" sim time: ")
            line = float(line)
            times.append(line)
        i += 1

cloud_velocities = np.array(cloud_velocities, dtype=float)
cloud_masses = np.array(cloud_masses, dtype=float)
times = np.array(times, dtype=float)

t_arr = np.linspace(np.amin(times), np.amax(times), len(cloud_velocities))

fontsize = 20
tickwidth = 1
pad = 0.1
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': 20})
xmin = np.amin(t_arr)/1e3
xmax = np.amax(t_arr)/1e3


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6.5))
ax.plot(t_arr/1e3, cloud_velocities, linewidth=2, color="k")
ax.set_xlabel("Time [Myr]", fontsize=fontsize)
ax.set_ylabel(r"Cloud velocity [km/s]", fontsize=fontsize)
ax.set_xlim(xmin, xmax)
ax.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True)
ax.set_xticks(np.linspace(0, np.amax(t_arr)/1e3, 5).round(1))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.savefig(os.path.join(pngdir, "cloud_velocity.png"))


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 6.5))
ax.plot(t_arr/1e3, cloud_masses, linewidth=4, color="k")
ax.set_xlabel("Time [Myr]", fontsize=fontsize)
ax.set_ylabel(r"Cloud mass $[M_\odot]$", fontsize=fontsize)
ax.set_xlim(xmin, xmax)
ax.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True)
ax.set_xticks(np.linspace(0, np.amax(t_arr)/1e3, 5).round(1))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.tight_layout()
plt.savefig(os.path.join(pngdir, f"cloud_mass_{date}.png"))