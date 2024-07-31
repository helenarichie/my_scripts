from hconfig import *

################################################################
date = "frontier/2024-02-19"
pad = 0.1
fontsize = 20
linewidth = 4
edges = True
tmax = None
title = r"surv, $r_{cl}/16$, $a=0.1~{\mu m}$"
name = "surv_16"
################################################################

plt.rcParams.update({'font.size': fontsize})

###############################
crc = True
frontier = False
###############################

########## data type ############
debugging = False
cloud_wind = True
testing = False
#################################

if crc:
  if debugging:
      basedir = f"/ix/eschneider/helena/data/debugging/{date}/"
  if cloud_wind:
      basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
  if testing:
      basedir = f"/ix/eschneider/helena/data/testing/{date}/"
if frontier:
  basedir = f"/lustre/orion/ast181/scratch/helenarichie/{date}/"

pngdir = os.path.join(basedir, "png/")
csvdir = os.path.join(basedir, "csv/")

time = []
with open(os.path.join(csvdir, "time.csv")) as f:
    for line in f:
        time.append(float(line))
time = np.array(time)

mass_dust = []
with open(os.path.join(csvdir, "mass_dust.csv")) as f:
    for line in f:
        mass_dust.append(float(line))
mass_dust = np.array(mass_dust)

sputter_hot = []
with open(os.path.join(csvdir, "sputter_hot.csv")) as f:
    for line in f:
        sputter_hot.append(float(line))

sputter = []
with open(os.path.join(csvdir, "sputter.csv")) as f:
    for line in f:
        sputter.append(float(line))

if edges:
    time_output = []
    with open(os.path.join(csvdir, "time_output.csv")) as f:
        for line in f:
            time_output.append(float(line))
    time_output = np.array(time_output)

    rate_dust = []
    with open(os.path.join(csvdir, "rate_dust.csv")) as f:
        for line in f:
            line = line.split(",")
            rate_dust.append(np.array(line, dtype=float))

    dt_out = time_output[1] - time_output[0]
    dt_out_end = time_output[-1] - time_output[-2]

    mass_out = []
    mass_cumulative = 0
    for k, rate in enumerate(rate_dust):
        rate = np.sum(rate)
        mass_cumulative += rate * dt_out
        mass_out.append(mass_cumulative)
    mass_out = np.array(mass_out)

sputter_tot, sputter_tot_hot = [], []
mass_cumulative, mass_cumulative_hot = 0, 0
for i, mass in enumerate(sputter):
    mass_cumulative += mass
    mass_cumulative_hot += sputter_hot[i]
    sputter_tot.append(mass_cumulative)
    sputter_tot_hot.append(mass_cumulative_hot)
sputter_tot = np.array(sputter_tot)
sputter_tot_hot = np.array(sputter_tot_hot)


if tmax == None:
    tmax = np.amax(time)

print(f"Total sputtered: {sputter_tot[-1]+sputter_tot_hot[-1]:e}")
if edges:
    print(f"Total out: {mass_out[-1]:e}")
    print(f"Sputtered + out: {mass_out[-1]+sputter_tot[-1]+sputter_tot_hot[-1]:e}")
print(f"Initial mass: {mass_dust[0]:e}")
print(f"Sputtered/initial: {(sputter_tot[time<=tmax][-1]+sputter_tot_hot[time<=tmax][-1])/mass_dust[0]}")

# ymin = np.amin([np.amin(mass_dust), np.amin(mass_out)]) - pad
# ymax = np.amax([np.amax(mass_dust), np.amax(mass_out)]) + pad
ymin = np.amin(mass_dust) - pad
ymax = np.amax(mass_dust) + pad
if edges:
    xmin = np.amin(time_output[time_output<=tmax]) - pad
    xmax = np.amax(time_output[time_output<=tmax]) + pad
else:
    xmin = np.amin(time[time<=tmax]) - pad
    # xmax = np.amax(time[time<=tmax]) + pad
    xmax = tmax + pad

plt.figure(figsize=(6, 5.5))
plt.plot(time[time<=tmax]/1e3, mass_dust[time<=tmax], label="in box", linewidth=linewidth, c="#1f77b4")
if edges:
    plt.plot(time_output[time_output<=tmax]/1e3, mass_out[time_output<=tmax], label="exited box", linewidth=linewidth, c="tab:orange")
plt.plot(time[time<=tmax]/1e3, sputter_tot[time<=tmax], label=r"sputtered ($T<10^6~K$)", linewidth=linewidth, c="#2ca02c")
plt.plot(time[time<=tmax]/1e3, sputter_tot_hot[time<=tmax], label=r"sputtered ($T>10^6~K$)", linewidth=linewidth, c="#d62728")
# plt.plot(time[time<=tmax]/1e3, sputter_tot[time<=tmax]+sputter_tot_hot[time<=tmax], label=r"sputtered (tot)", c="k", linewidth=linewidth)
plt.legend(fontsize=fontsize-10, loc="center left")
plt.xlabel("Time [Myr]")
plt.ylabel(r"Mass $[M_\odot]$")
plt.xlim(xmin/1e3, xmax/1e3)
plt.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True)
plt.xticks(np.linspace(0, xmax/1e3-pad, 5).round(0))
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.title(title, pad=pad, fontsize=fontsize-2)
plt.tight_layout()
plt.savefig(os.path.join(pngdir, f"dust_mass_{name}.png"))
