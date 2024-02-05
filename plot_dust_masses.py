from hconfig import *
from csv import writer

################################################################
date = "2024-02-04"
pad = 0.1
fontsize = 20
linewidth = 4
################################################################

plt.rcParams.update({'font.size': fontsize})

basedir = f"/ix/eschneider/helena/data/testing/{date}/"
pngdir = os.path.join(basedir, "png/")
csvdir = os.path.join(basedir, "csv/")

time = []
with open(os.path.join(csvdir, "time.csv")) as f:
    for line in f:
        time.append(float(line))
time = np.array(time)

time_output = []
with open(os.path.join(csvdir, "time_output.csv")) as f:
    for line in f:
        time_output.append(float(line))
time_output = np.array(time_output)

mass_dust = []
with open(os.path.join(csvdir, "mass_dust.csv")) as f:
    for line in f:
        mass_dust.append(float(line))

sputter_hot = []
with open(os.path.join(csvdir, "sputter_hot.csv")) as f:
    for line in f:
        sputter_hot.append(float(line))

sputter = []
with open(os.path.join(csvdir, "sputter.csv")) as f:
    for line in f:
        sputter.append(float(line))

rate_dust = []
with open(os.path.join(csvdir, "rate_dust.csv")) as f:
    for line in f:
        line = line.split(",")
        rate_dust.append(np.array(line, dtype=float))

dt_out = time_output[1] - time_output[0]

mass_out = []
mass_cumulative = 0
for rate in rate_dust:
    rate = np.sum(rate)
    mass_cumulative += rate * dt_out
    mass_out.append(mass_cumulative)

sputter_tot = []
mass_cumulative = 0
for mass in sputter:
    mass_cumulative += mass
    sputter_tot.append(mass_cumulative)

sputter_tot = np.array(sputter_tot)

print(f"Total sputtered: {sputter_tot[-1]:e}")
print(f"Total out: {mass_out[-1]:e}")
print(f"Sputtered + out: {mass_out[-1]+sputter_tot[-1]:e}")
print(f"Initial mass: {mass_dust[0]:e}")

ymin = np.amin([np.amin(mass_dust), np.amin(mass_out)]) - pad
ymax = np.amax([np.amax(mass_dust), np.amax(mass_out)]) + pad
xmin = np.amin(time_output) - pad
xmax = np.amax(time_output) + pad

plt.plot(figsize=(6.5, 6.5))
plt.plot(time/1e3, mass_dust, label="in box", linewidth=linewidth)
plt.plot(time_output/1e3, mass_out, label="exited box", linewidth=linewidth)
plt.plot(time/1e3, sputter_tot, label="sputtered", linewidth=linewidth)
plt.legend(fontsize=fontsize-5)
plt.xlabel("Time [Myr]")
plt.ylabel(r"Mass $[M_\odot]$")
plt.xlim(xmin/1e3, xmax/1e3)
plt.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True)
plt.xticks(np.linspace(0, np.amax(time_output/1e3), 5).round(1))
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.tight_layout()
plt.savefig(os.path.join(pngdir, "dust_mass.png"))
