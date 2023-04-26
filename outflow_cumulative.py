from hconfig import *

################# hard-coded, fill these in! ###################
date = "2023-03-22"
cat = True
rho_cl_init = 1e-24  # g cm^-3
r_cl_init = 5 * 3.086e+18  # cm, initial cloud radius
chi = 1e2  # background-wind density contrast
v_b = 1e3 * 1e5  # km/s * (cm/s / km/s), background wind speed
# xlims = (0, 1)
################################################################

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/" # crc
datadir = os.path.join(basedir, "hdf5/full/")
csvdir = os.path.join(basedir, "csv")
pngdir = os.path.join(basedir, "png/outflow/")

data = ReadHDF5(datadir, fnum=0, nscalar=1, cat=cat)
head = data.head
nx, ny, nz = head["dims"]
dx = head["dx"][0] # M_sun

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

mass_i = 0
mass_out = []
for rate in rate_dust:
    rate_i = rate[3]
    mass_i += rate_i * dt
    mass_out.append(mass_i)

fig = plt.figure(figsize=(14, 10))
plt.plot(t_arr/1e6, mass_out, linewidth=4)
plt.title("+x")
plt.xlabel("Time [Myr]")
plt.ylabel(r"Dust Mass$~[M_\odot]$")
plt.tight_layout()
plt.savefig(os.path.join(pngdir, "mass_out.png"))