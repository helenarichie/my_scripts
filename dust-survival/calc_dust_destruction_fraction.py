import sys
sys.path.insert(0, "/ix/eschneider/helena/code/my_scripts")
from hconfig import *

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

##################################################################
date = "2024-01-12"
cat = True
tmax_1105 = 52.1e6
tmax_1103 = 52.1e6
tmax_1106 = 58.9e6
tmax_1107 = 58.9e6
tmax_1108 = 72.2e6
tmax_0426 = 2.4e6
tmax_1208 = 1.9e6
tmax_1206 = 48.5e6
tmax_0503 = 2.4e6
tmax_0513 = 1.4e6
tmax_1212 = 60e6
tmax_0111 = tmax_0112 = 55.1e6
tmax = tmax_0112
##################################################################

##################################################################
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
projdir = os.path.join(basedir, "hdf5/full/")
pngdir = os.path.join(basedir, "png/")
csvdir = os.path.join(basedir, "csv/")
##################################################################

data = ReadHDF5(projdir, cat=cat, fnum=0)
head = data.head
conserved = data.conserved
dx = head["dx"][0]

t_arr = [] # yr
with open(os.path.join(csvdir, "t_arr.csv")) as f:
    for line in f:
        line = line.split(",")
        t = line[0].strip("\n").strip("[").strip("]")
        t_arr.append(float(t))
t_arr = np.array(t_arr)

n_steps = len(t_arr)

dt = (t_arr[1] - t_arr[0])/1e3

rate_dust = np.zeros((n_steps, 6))  # M_sun / yr
with open(os.path.join(csvdir, "rate_dust.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_dust[i] = np.array(line, dtype=float)

mass_dust_tot = []
with open(os.path.join(csvdir, "rho_d_tot.csv")) as f:
    for line in f:
        line = line.split(",")
        density = float(line[0]) * density_conversion
        mass = density * dx**3
        mass_dust_tot.append(mass)
mass_dust_tot = np.array(mass_dust_tot)

mass_i = 0
mass_out = []
mass_lost = []
for i, rate in enumerate(rate_dust):
    rate_i = np.sum(rate)
    mass_lost.append(rate_i * dt)
    mass_i += rate_i * dt
    mass_out.append(mass_i)
mass_out = np.array(mass_out)
mass_lost = np.array(mass_lost)

mass_destroyed = np.zeros(len(mass_lost))
mass_destroyed[0] = 0
mass_dest_cum = 0
for i, mass in enumerate(mass_dust_tot):
    if (i+1 < len(mass_lost)):
        mass_dest_cum += (mass - mass_dust_tot[i+1] - mass_lost[i+1])
        mass_destroyed[i+1] = mass_dest_cum

print(f"Inital mass: {mass_dust_tot[0]} M_sun")
print(f"Fraction sputtered: {mass_destroyed[t_arr<=tmax][-1]/mass_dust_tot[0]}")
print(f"Out + sputtered fraction: {mass_destroyed[t_arr<=tmax][-1]/mass_dust_tot[0] + mass_out[t_arr<=tmax][-1]/mass_dust_tot[0]}")
print(f"t_max: {tmax/1e6} Myr")
