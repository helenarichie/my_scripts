from hconfig import *

######### hard-coded values! #########
date = "2023-03-10"
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5/full/")
pngdir = os.path.join(basedir, "png/")
csvdir = os.path.join(basedir, "csv/")
cat = True
cutoff = 0
######################################

# Constants
M_sun = 5.02785e-34 # solar mass per gram

data = ReadHDF5(datadir, fnum=0, nscalar=1, cat=cat)
head = data.head
nx, ny, nz = head["dims"]
dx = data.dx_cgs()[0]
gamma = head["gamma"]

area_x = ny*nz*dx**2
area_y = nx*nz*dx**2
area_z = nx*ny*dx**2

t_arr = []
with open(os.path.join(csvdir, "t_arr.csv")) as f:
    for line in f:
        line = line.split(",")
        t_arr.append(float(line[0]))
t_arr = np.array(t_arr)

n_steps = len(t_arr)

rate_cloud = np.zeros((n_steps, 6))
with open(os.path.join(csvdir, "outflow_cloud.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_cloud[i] = np.array(line, dtype=float)

flux_cloud = np.zeros((n_steps, 6))
with open(os.path.join(csvdir, "flux_cloud.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        flux_cloud[i] = np.array(line, dtype=float)

rate_dust = np.zeros((n_steps, 6))
with open(os.path.join(csvdir, "outflow_dust.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_dust[i] = np.array(line, dtype=float)

flux_dust = np.zeros((n_steps, 6))
with open(os.path.join(csvdir, "flux_dust.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        flux_dust[i] = np.array(line, dtype=float)

faces = [["-x", "-y", "-z"],
         ["+x", "+y", "+z"]]

def plot_outflow(rate, faces, fig_name):
    fig, axs = plt.subplots(3, 2, figsize=(30,30))
    for i, sign in enumerate(faces):
        for j, face in enumerate(sign):
            print(j, i)
            axs[j][i].plot(t_arr/1e6, rate[:, i+j], linewidth="5")
            axs[j][i].set_xlabel(r"Time$~[Myr]$")
            axs[j][i].set_ylabel(r'Outflow Rate $[g\,s^-1]$')
            axs[j][i].set_title(face)
    plt.suptitle(fig_name)
    plt.savefig(pngdir + fig_name)


def plot_flux(flux, faces, fig_name):
    fig, axs = plt.subplots(3, 2, figsize=(30,30))
    for i, sign in enumerate(faces):
        for j, face in enumerate(sign):
            axs[j][i].plot(t_arr/1e6, flux[:, i+j], linewidth="5", color="orange")
            axs[j][i].set_xlabel(r"Time$~[Myr]$")
            axs[j][i].set_ylabel(r'Mass At Boundary $[g]$')
            axs[j][i].set_title(face)
    plt.suptitle(fig_name)
    plt.savefig(pngdir + fig_name)


plot_outflow(rate_cloud, faces, "cloud_rate")
plot_outflow(rate_dust, faces, "dust_rate")
plot_flux(flux_cloud, faces, "cloud_flux")
plot_flux(flux_dust, faces, "dust_flux")