from hconfig import *

######### hard-coded values! #########
date = "2023-03-14"
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5/full/")
pngdir = os.path.join(basedir, "png/outflow/")
csvdir = os.path.join(basedir, "csv/")
cat = False
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
    cumulative_flux = np.zeros(6)
    for i, line in enumerate(f):
        line = np.array(line.split(","), dtype=float)
        cumulative_flux += line
        flux_cloud[i] = cumulative_flux

rate_dust = np.zeros((n_steps, 6))
with open(os.path.join(csvdir, "outflow_dust.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_dust[i] = np.array(line, dtype=float)

flux_dust = np.zeros((n_steps, 6))
with open(os.path.join(csvdir, "flux_dust.csv")) as f:
    cumulative_flux = np.zeros(6)
    for i, line in enumerate(f):
        line = np.array(line.split(","), dtype=float)
        cumulative_flux += line
        flux_dust[i] = cumulative_flux

faces = [["-x", "-y", "-z"],
         ["+x", "+y", "+z"]]

def plot_outflow(rate, faces, fig_name):
    fig, axs = plt.subplots(2, 3, figsize=(40,20))
    count = 0
    for i, sign in enumerate(faces):
        for j, face in enumerate(sign):
            axs[i][j].plot(t_arr/1e6, rate[:, count], linewidth="5")
            axs[i][j].set_xlabel(r"Time$~[Myr]$")
            axs[i][j].set_ylabel(r'Outflow Rate $[g\,s^-1]$')
            axs[i][j].set_title(face)
            count += 1
    #plt.tight_layout()
    plt.suptitle(fig_name)
    plt.savefig(pngdir + fig_name)


def plot_flux(flux, faces, fig_name):
    fig, axs = plt.subplots(2, 3, figsize=(40,20))
    count = 0
    for i, sign in enumerate(faces):
        for j, face in enumerate(sign):
            axs[i][j].plot(t_arr/1e6, flux[:, count], linewidth="5", color="orange")
            axs[i][j].set_xlabel(r"Time$~[Myr]$")
            axs[i][j].set_ylabel(r'Mass At Boundary $[g]$')
            axs[i][j].set_title(face)
            count += 1
    #plt.tight_layout()
    plt.subplots_adjust(wspace=0.3)
    plt.suptitle(fig_name)
    plt.savefig(pngdir + fig_name)


plot_outflow(rate_cloud, faces, "cloud_rate")
plot_outflow(rate_dust, faces, "dust_rate")
plot_flux(flux_cloud, faces, "cloud_flux")
plot_flux(flux_dust, faces, "dust_flux")