from hconfig import *

######### hard-coded values! #########
date = "2023-03-20"
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5/full/")
pngdir = os.path.join(basedir, "png/outflow/")
csvdir = os.path.join(basedir, "csv/")
cat = True
cgs = True
rho_cl_init = 1e-24 # g cm^-3
r_cl_init = 5 * 3.086e+18 # cm
######################################

# Constants
M_sun = 5.02785e-34 # solar mass per gram

data = ReadHDF5(datadir, fnum=0, nscalar=1, cat=cat)
head = data.head
nx, ny, nz = head["dims"]
dx = data.dx_cgs()[0]
gamma = head["gamma"]

m_cl_init = rho_cl_init * dx**3
n_cells_init = 4/3 * np.pi * r_cl_init**3 / dx**3
m_cl_tot_init = n_cells_init * m_cl_init
m_dust_tot_init = 1e-2 * n_cells_init * m_cl_init

area_x = ny*nz*dx**2
area_y = nx*nz*dx**2
area_z = nx*ny*dx**2

t_arr = []
with open(os.path.join(csvdir, "t_arr.csv")) as f:
    for line in f:
        line = line.split(",")
        t = line[0].strip("\n").strip("[").strip("]")
        t_arr.append(float(t))
t_arr = np.array(t_arr)

n_steps = len(t_arr)

rate_cloud = np.zeros((n_steps, 6))  # g/s
with open(os.path.join(csvdir, "outflow_cloud.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_cloud[i] = np.array(line, dtype=float)

flux_cloud = np.zeros((n_steps, 6))  # g
with open(os.path.join(csvdir, "flux_cloud.csv")) as f:
    cumulative_flux = np.zeros(6)
    for i, line in enumerate(f):
        line = np.array(line.split(","), dtype=float)
        cumulative_flux += line
        flux_cloud[i] = cumulative_flux

rate_dust = np.zeros((n_steps, 6))  # g/s
with open(os.path.join(csvdir, "outflow_dust.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_dust[i] = np.array(line, dtype=float)

flux_dust = np.zeros((n_steps, 6))  # g
with open(os.path.join(csvdir, "flux_dust.csv")) as f:
    cumulative_flux = np.zeros(6)
    for i, line in enumerate(f):
        line = np.array(line.split(","), dtype=float)
        cumulative_flux += line
        flux_dust[i] = cumulative_flux

cloud_flux_min = np.zeros(6)
cloud_flux_max = np.zeros(6)
dust_flux_min = np.zeros(6)
dust_flux_max = np.zeros(6)
for i in range(0, 6):
    cloud_flux_min[i] = np.amin(flux_cloud[:, i])
    cloud_flux_max[i] = np.amax(flux_cloud[:, i])
    dust_flux_min[i] = np.amin(flux_dust[:, i])
    dust_flux_max[i] = np.amax(flux_dust[:, i])

ax2_lims_cloud = np.column_stack((cloud_flux_min/m_cl_tot_init, cloud_flux_max/m_cl_tot_init))
ax2_lims_dust = np.column_stack((dust_flux_min/m_dust_tot_init, dust_flux_max/m_dust_tot_init))

faces = [["-x", "-y", "-z"],
         ["+x", "+y", "+z"]]

def plot_outflow(rate, faces, fig_name, cgs):
    fig, axs = plt.subplots(2, 3, figsize=(40,20))
    count = 0
    for i, sign in enumerate(faces):
        for j, face in enumerate(sign):
            axs[i][j].plot(t_arr/1e6, rate[:, count], linewidth="5")
            axs[i][j].set_xlabel(r"Time$~[Myr]$")
            axs[i][j].set_ylabel(r'Outflow Rate $[g\,s^-1]$')
            axs[i][j].set_title(face)
            ax2 = axs[i][j].twinx()
            ax2.set_ylim(np.amin(rate[:, count]*M_sun/yr_in_s), np.amax(rate[:, count]*M_sun/yr_in_s))
            ax2.set_ylabel(r"Outflow Rate $[M_\odot\,yr^{-1}]$", rotation=270, labelpad=55)
            count += 1
    #plt.tight_layout()
    plt.subplots_adjust(wspace=0.5)
    plt.suptitle(fig_name)
    plt.savefig(pngdir + fig_name)


def plot_flux(flux, faces, fig_name, cgs, ax2_lims):
    fig, axs = plt.subplots(2, 3, figsize=(40,20))
    count = 0
    for i, sign in enumerate(faces):
        for j, face in enumerate(sign):
            if cgs:
                axs[i][j].plot(t_arr/1e6, flux[:, count], linewidth="5", color="orange")
                axs[i][j].set_ylabel(r'Mass At Boundary $[g]$')
            else:
                axs[i][j].plot(t_arr/1e6, flux[:, count]*M_sun, linewidth="5", color="orange")
                axs[i][j].set_ylabel(r'Mass At Boundary $[M_\odot]$')
            ax2 = axs[i][j].twinx()
            ax2.set_ylim(ax2_lims[count])
            ax2.set_ylabel("Fraction of Initial Mass", rotation=270, labelpad=35)
            axs[i][j].set_xlabel(r"Time$~[Myr]$")
            axs[i][j].set_title(face)
            count += 1
    #plt.tight_layout()
    plt.subplots_adjust(wspace=0.5)
    plt.suptitle(fig_name)
    plt.savefig(pngdir + fig_name)


plot_outflow(rate_cloud, faces, "cloud_rate", cgs)
plot_outflow(rate_dust, faces, "dust_rate", cgs)
plot_flux(flux_cloud, faces, "cloud_flux", cgs, ax2_lims_cloud)
plot_flux(flux_dust, faces, "dust_flux", cgs, ax2_lims_dust)