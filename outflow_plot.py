from hconfig import *

################# hard-coded, fill these in! ###################
date = "2023-04-26"
cat = True
rho_cl_init = 1e-24  # g cm^-3
r_cl_init = 5 * 3.086e+18  # cm, initial cloud radius
chi = 1e2  # background-wind density contrast
v_b = 1e3 * 1e5  # km/s * (cm/s / km/s), background wind speed
# xlims = (0, 1)
################################################################

# Klein et al. (1994) cloud crushing time
tau_cc = chi**(1/2) * r_cl_init / v_b / yr_in_s  # yr

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/" # crc
datadir = os.path.join(basedir, "hdf5/full/")
csvdir = os.path.join(basedir, "csv")
pngdir = os.path.join(basedir, "png/outflow/")

data = ReadHDF5(datadir, fnum=0, dust=True, cat=cat)
head = data.head
nx, ny, nz = head["dims"]
dx = head["dx"] # M_sun

m_cl_init = rho_cl_init * dx**3
n_cells_init = 4/3 * np.pi * r_cl_init**3 / dx**3
m_cl_tot_init = n_cells_init * m_cl_init
m_dust_tot_init = 1e-2 * n_cells_init * m_cl_init

# kpc^2
area_x = ny*nz*dx**2
area_y = nx*nz*dx**2
area_z = nx*ny*dx**2

faces = [["-x", "-y", "-z"],
         ["+x", "+y", "+z"]]

t_arr = [] # yr
with open(os.path.join(csvdir, "t_arr.csv")) as f:
    for line in f:
        line = line.split(",")
        t = line[0].strip("\n").strip("[").strip("]")
        t_arr.append(float(t))
t_arr = np.array(t_arr)

n_steps = len(t_arr)

rate_cl = np.zeros((n_steps, 6))  # M_sun / yr
with open(os.path.join(csvdir, "rate_cl.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_cl[i] = np.array(line, dtype=float)

mass_cl = np.zeros((n_steps, 6))  # M_sun / yr
with open(os.path.join(csvdir, "mass_cl.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        mass_cl[i] = np.array(line, dtype=float)

rate_gas = np.zeros((n_steps, 6))  # M_sun / yr
with open(os.path.join(csvdir, "rate_gas.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_gas[i] = np.array(line, dtype=float)

mass_gas = np.zeros((n_steps, 6))  # M_sun / yr
with open(os.path.join(csvdir, "mass_gas.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        mass_gas[i] = np.array(line, dtype=float)

mass_dust = np.zeros((n_steps, 6))  # M_sun / yr
with open(os.path.join(csvdir, "mass_dust.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        mass_dust[i] = np.array(line, dtype=float)

rate_dust = np.zeros((n_steps, 6))  # M_sun / yr
with open(os.path.join(csvdir, "rate_dust.csv")) as f:
    for i, line in enumerate(f):
        line = line.split(",")
        rate_dust[i] = np.array(line, dtype=float)

def plot_rate(rate, faces, fig_name, scale):
    fig, axs = plt.subplots(2, 3, figsize=(40,20))
    count = 0
    for i, sign in enumerate(faces):
        for j, face in enumerate(sign):
            if scale == "log":
                axs[i][j].semilogy(t_arr/1e6, rate[:, count], linewidth="5")
            if scale == "linear":
                axs[i][j].plot(t_arr/1e6, rate[:, count], linewidth="5")
            # axs[i][j].set_xlim(xlims)
            axs[i][j].set_xlabel(r"Time$~[Myr]$")
            axs[i][j].set_ylabel(r'Outflow Rate $[M_\odot\,kyr^-1]$')
            axs[i][j].set_title(face, fontsize=40, y=1.15)
            x2 = axs[i][j].twiny()
            t_cc_range = np.arange(np.amin(t_arr/(tau_cc)), np.amax(t_arr/(tau_cc)), step=4)
            x2.set_xlim(np.amin(t_cc_range), np.amax(t_cc_range))
            x2.set_xlabel("$t/t_{cc}$")
            x2.set_xticks(t_cc_range[1:])
            count += 1
    #plt.tight_layout()
    plt.subplots_adjust(hspace=0.5, wspace=0.7)
    plt.suptitle(fig_name, fontsize=40)
    plt.savefig(pngdir + fig_name)


def plot_mass(mass, faces, fig_name, scale):
    fig, axs = plt.subplots(2, 3, figsize=(40,20))
    count = 0
    for i, sign in enumerate(faces):
        for j, face in enumerate(sign):
            if scale == "log":
                axs[i][j].semilogy(t_arr/1e6, mass[:, count], linewidth="5", color="orange")
            if scale == "linear":
                axs[i][j].plot(t_arr/1e6, mass[:, count], linewidth="5", color="orange")
            axs[i][j].set_ylabel(r'Mass At Boundary $[M_\odot]$')
            axs[i][j].set_xlabel(r"Time$~[Myr]$")
            axs[i][j].set_title(face, fontsize=40, y=1.15)
            # axs[i][j].set_xlim(xlims)
            x2 = axs[i][j].twiny()
            t_cc_range = np.arange(np.amin(t_arr/(tau_cc)), np.amax(t_arr/(tau_cc)), step=4)
            x2.set_xlim(np.amin(t_cc_range), np.amax(t_cc_range))
            x2.set_xlabel("$t/t_{cc}$")
            x2.set_xticks(t_cc_range[1:])
            count += 1
    #plt.tight_layout()
    plt.subplots_adjust(hspace=0.5, wspace=0.7)
    plt.suptitle(fig_name, fontsize=40)
    plt.savefig(pngdir + fig_name)


plot_rate(rate_cl, faces, "rate_cl", "linear")
plot_mass(mass_cl, faces, "mass_cl", "linear")
plot_rate(rate_gas, faces, "rate_gas", "log")
plot_mass(mass_gas, faces, "mass_gas", "log")
plot_mass(mass_dust, faces, "mass_dust", "log")
plot_rate(rate_dust, faces, "rate_dust", "log")