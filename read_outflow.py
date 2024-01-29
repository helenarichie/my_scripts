from hconfig import *

date = "2024-01-29"

basedir = f"/ix/eschneider/helena/data/testing/{date}/"
pngdir = os.path.join(basedir, "png/")

mass_cloud, rate_cloud, mass_cloud_bndry, mass_dust, rate_dust, mass_dust_bndry = [], [], [], [], [], []

total = [mass_cloud, rate_cloud, mass_cloud_bndry, mass_dust, rate_dust, mass_dust_bndry]
times, dts = [], []
i = 0
strings = [" Cloud mass: ", "Cloud rate: ", "Cloud mass bndry: ", "Dust mass: ", 
           "Dust rate: ", "Dust mass bndry: "]
with open(os.path.join(basedir, "output.log")) as f:
    for line in f:
        if line.startswith("@@ "):
            line = line.lstrip("@@").rstrip(" \n")
            line = line.split("  ")
            for j, entry in enumerate(line):
                entry = entry.lstrip(strings[j])
                total[j].append(float(entry))

with open(os.path.join(basedir, "output.log")) as f:
    for line in f:
        if line.startswith("n_step"):
            line = line.split(" ")
            times.append(float(line[6]))
            dts.append(float(line[11]))
        i += 1

times = np.array(times, dtype=float)

dts = []
for k, t in enumerate(times):
    if (k+1) < len(times):
        dts.append(times[k+1]-times[k])

last = dts[-1]
dts.append(last)

mass_dust_i = 0
mass_dust_out = []
mass_dust_lost = []
for i, rate in enumerate(rate_dust):
    rate_i = np.sum(rate)
    mass_dust_lost.append(rate_i * dts[i])
    mass_dust_i += rate_i * dts[i]
    mass_dust_out.append(mass_dust_i)
mass_dust_out = np.array(mass_dust_out)
mass_dust_lost = np.array(mass_dust_lost)

mass_destroyed = np.zeros(len(mass_dust_lost))
mass_destroyed[0] = 0
mass_dest_cum = 0
for i, mass in enumerate(mass_dust):
    if (i+1 < len(mass_dust_lost)):
        mass_dest_cum += (mass - mass_dust[i+1] - mass_dust_lost[i+1])
        mass_destroyed[i+1] = mass_dest_cum

mass_cloud_i = 0
mass_cloud_out = []
mass_cloud_lost = []
for i, rate in enumerate(rate_cloud):
    rate_i = np.sum(rate)
    mass_cloud_lost.append(rate_i * dts[i])
    mass_cloud_i += rate_i * dts[i]
    mass_cloud_out.append(mass_cloud_i)
mass_cloud_out = np.array(mass_cloud_out)
mass_cloud_lost = np.array(mass_cloud_lost)

fontsize = 20
tickwidth = 1
pad = 0.1
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': 20})
xmin = np.amin(times)/1e3
xmax = np.amax(times)/1e3

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 6.5))
ax.plot(times/1e3, mass_cloud, linewidth=4, color="#49b4ab")
ax.plot(times/1e3, mass_cloud_out, linewidth=4, linestyle="--", c="#49b4ab", label="exited box")
ax.set_xlabel("Time [Myr]", fontsize=fontsize)
ax.set_ylabel(r"Cloud mass [M_sun]", fontsize=fontsize)
ax.set_xlim(xmin, xmax)
ax.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True)
ax.set_xticks(np.linspace(0, np.amax(times)/1e3, 5).round(1))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax.set_title("runtime")
plt.savefig(os.path.join(pngdir, "cloud_mass.png"))

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 6.5))
ax.plot(times/1e3, mass_dust, linewidth=4, color="#d43a4f")
ax.plot(times/1e3, mass_dust_out, linewidth=4, linestyle="--", c="#d43a4f", label="exited box")
ax.plot(times/1e3, mass_destroyed, linewidth=4, c="k", zorder=10, label="sputtered")
ax.set_xlabel("Time [Myr]", fontsize=fontsize)
ax.set_ylabel(r"dust mass [M_sun]", fontsize=fontsize)
ax.set_xlim(xmin, xmax)
ax.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True)
ax.set_xticks(np.linspace(0, np.amax(times)/1e3, 5).round(1))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax.set_title("runtime")
plt.savefig(os.path.join(pngdir, "dust_mass.png"))