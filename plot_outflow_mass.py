from hconfig import *

######### hard-coded values! ###########
date = "2023-03-07"
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5/full/")
pngdir = os.path.join(basedir, "png/outflow/")
csvdir = os.path.join(basedir, "csv/")
cat = False
cgs = True
rho_cl_init = 1e-24  # g cm^-3
r_cl_init = 5 * 3.086e+18  # cm, initial cloud radius
chi = 1e2  # background-wind density contrast
v_b = 1e3 * 1e5  # km/s * (cm/s / km/s), background wind speed
#xlims = (0, 1)
########################################

# Constants
M_sun = 5.02785e-34  # solar mass per gram

# Klein et al. (1994) cloud crushing time
tau_cc = chi**(1/2) * r_cl_init / v_b  # s

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

rate_tot_cloud = np.zeros((n_steps, 6))  # g/s
with open(os.path.join(csvdir, "outflow_tot_cl.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_tot_cloud[i] = np.array(line, dtype=float)

rate_dense_cloud = np.zeros((n_steps, 6))  # g/s
with open(os.path.join(csvdir, "outflow_dense_cl.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_dense_cloud[i] = np.array(line, dtype=float)

masses_cloud = np.zeros((n_steps, 6))  # g
with open(os.path.join(csvdir, "masses_cl.csv")) as f:
    cumulative_flux = np.zeros(6)
    for i, line in enumerate(f):
        line = np.array(line.split(","), dtype=float)
        cumulative_flux += line
        masses_cloud[i] = cumulative_flux

rate_tot_dust = np.zeros((n_steps, 6))  # g/s
with open(os.path.join(csvdir, "outflow_tot_dust.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_tot_dust[i] = np.array(line, dtype=float)

rate_dense_dust = np.zeros((n_steps, 6))  # g/s
with open(os.path.join(csvdir, "outflow_dense_dust.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_dense_dust[i] = np.array(line, dtype=float)

masses_dust = np.zeros((n_steps, 6))  # g
with open(os.path.join(csvdir, "masses_dust.csv")) as f:
    cumulative_flux = np.zeros(6)
    for i, line in enumerate(f):
        line = np.array(line.split(","), dtype=float)
        cumulative_flux += line
        masses_dust[i] = cumulative_flux

cloud_mass_min = np.zeros(6)
cloud_mass_max = np.zeros(6)
dust_mass_min = np.zeros(6)
dust_mass_max = np.zeros(6)
for i in range(0, 6):
    cloud_mass_min[i] = np.amin(masses_cloud[:, i])
    cloud_mass_max[i] = np.amax(masses_cloud[:, i])
    dust_mass_min[i] = np.amin(masses_dust[:, i])
    dust_mass_max[i] = np.amax(masses_dust[:, i])

ax2_lims_cloud = np.column_stack((cloud_mass_min/m_cl_tot_init, cloud_mass_max/m_cl_tot_init))
ax2_lims_dust = np.column_stack((dust_mass_min/m_dust_tot_init, dust_mass_max/m_dust_tot_init))

faces = [["-x", "-y", "-z"],
         ["+x", "+y", "+z"]]

def plot_outflow(rate, rate_dense, faces, fig_name, cgs):
    fig, axs = plt.subplots(2, 3, figsize=(40,20))
    count = 0
    for i, sign in enumerate(faces):
        for j, face in enumerate(sign):
            #axs[i][j].plot(t_arr/1e6, rate[:, count], linewidth="5", label="tot")
            axs[i][j].plot(t_arr/1e6, rate_dense[:, count], linewidth="5", label="dense")
            #axs[i][j].set_xlim(xlims)
            axs[i][j].set_xlabel(r"Time$~[Myr]$")
            axs[i][j].set_ylabel(r'Outflow Rate $[g\,s^-1]$')
            axs[i][j].set_title(face, fontsize=40, y=1.1)
            #x2 = axs[i][j].twiny()
            #t_cc_range = np.arange(np.amin(t_arr/(tau_cc/yr_in_s)), np.amax(t_arr/(tau_cc/yr_in_s)), step=4)
            #x2.set_xlim(np.amin(t_cc_range), np.amax(t_cc_range))
            #x2.set_xlabel("$t/t_{cc}$")
            #x2.set_xticks(t_cc_range[1:])
            y2 = axs[i][j].twinx()
            y2.set_ylim(np.amin(rate[:, count]*M_sun/yr_in_s), np.amax(rate[:, count]*M_sun/yr_in_s))
            y2.set_ylabel(r"Outflow Rate $[M_\odot\,yr^{-1}]$", rotation=270, labelpad=55)
            count += 1
    #plt.tight_layout()
    axs[0][0].legend()
    plt.subplots_adjust(hspace=0.5, wspace=0.7)
    plt.suptitle(fig_name, fontsize=40)
    plt.savefig(pngdir + fig_name)


def plot_m_out(mass, faces, fig_name, cgs, ax2_lims):
    fig, axs = plt.subplots(2, 3, figsize=(40,20))
    count = 0
    for i, sign in enumerate(faces):
        for j, face in enumerate(sign):
            if cgs:
                axs[i][j].plot(t_arr/1e6, mass[:, count], linewidth="5", color="orange")
                axs[i][j].set_ylabel(r'Mass At Boundary $[g]$')
            else:
                axs[i][j].plot(t_arr/1e6, mass[:, count]*M_sun, linewidth="5", color="orange")
                axs[i][j].set_ylabel(r'Mass At Boundary $[M_\odot]$')
            axs[i][j].set_xlabel(r"Time$~[Myr]$")
            axs[i][j].set_title(face, fontsize=40, y=1.1)
            #axs[i][j].set_xlim(xlims)
            #x2 = axs[i][j].twiny()
            #t_cc_range = np.arange(np.amin(t_arr/(tau_cc/yr_in_s)), np.amax(t_arr/(tau_cc/yr_in_s)), step=4)
            #x2.set_xlim(np.amin(t_cc_range), np.amax(t_cc_range))
            #x2.set_xlabel("$t/t_{cc}$")
            #x2.set_xticks(t_cc_range[1:])
            y2 = axs[i][j].twinx()
            y2.set_ylim(np.amin(mass[:, count]*M_sun), np.amax(mass[:, count]*M_sun))
            y2.set_ylabel(r"Mass At Boundary $[M_\odot\,yr^{-1}]$", rotation=270, labelpad=55)
            count += 1
    #plt.tight_layout()
    plt.subplots_adjust(hspace=0.5, wspace=0.7)
    plt.suptitle(fig_name, fontsize=40)
    plt.savefig(pngdir + fig_name)


plot_outflow(rate_tot_cloud, rate_dense_cloud, faces, "cloud_rate", cgs)
plot_outflow(rate_tot_dust, rate_dense_dust, faces, "dust_rate", cgs)
plot_m_out(masses_cloud, faces, "cloud_m_out", cgs, ax2_lims_cloud)
plot_m_out(masses_dust, faces, "dust_m_out", cgs, ax2_lims_dust)