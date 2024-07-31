from hconfig import *
from csv import writer

date = "frontier/2024-02-19"
datestr = "0219"

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
pngdir = os.path.join(basedir, "png/")
csvdir = os.path.join(basedir, "csv/")

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

velocity_cumulative = 500
with open(os.path.join(csvdir, "cloud_velocities.csv"), "w") as f:
    writer_obj = writer(f)
    for i, velocity in enumerate(cloud_velocities):
        velocity_cumulative -= velocity
        writer_obj.writerow([times[i], cloud_masses[i], cloud_velocities[i]])
    f.close()
print(velocity_cumulative)

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
plt.savefig(os.path.join(pngdir, f"cloud_mass_{datestr}.png"))