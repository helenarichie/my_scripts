from hconfig import *

################# hard-coded, fill these in! ###################
date = "2023-03-08"
################################################################

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5",)
csvdir = os.path.join(basedir, "csv")
pngdir = os.path.join(basedir, "png")

rho_d_tot = []  # total dust density in sim volume for each time
rho_cl_tot = []  # total cloud density in sim volume for each time
T_cl_avg = []  # average T for each time
t_arr = []

print("\nReading in data...")

with open(os.path.join(csvdir, "rho_d_tot.csv")) as f:
    for line in f:
        line = line.split(",")
        rho_d_tot.append(float(line[0]))

with open(os.path.join(csvdir, "rho_cl_tot.csv")) as f:
    for line in f:
        line = line.split(",")
        rho_cl_tot.append(float(line[0]))

with open(os.path.join(csvdir, "t_arr.csv")) as f:
    for line in f:
        line = line.strip().split(",")
        line = line[0].strip("]").strip("[")
        t_arr.append(float(line))

print("\nDone.")

rho_d_tot = np.array(rho_d_tot) # why was I multiplying this by 256?
rho_cl_tot = np.array(rho_cl_tot)
t_arr = np.array(t_arr)

rho_cl_tot_i = rho_cl_tot[0]
rho_d_tot_i = rho_d_tot[0]

print("\nPlotting...")

fig, ax = plt.subplots(2, 1, figsize=(15,22))

ax[0].semilogx(t_arr, rho_cl_tot/rho_cl_tot_i, linewidth=7, color="red")
ax[0].set_xlim(1, 1.2e6)
ax[0].set_ylabel(r'$\rho_{cl}/\rho_{cl,init}$')
ax[0].set_xlabel(r'Time [yr]')
ax[0].set_title(r"Cloud")

ax[1].semilogx(t_arr, rho_d_tot/rho_d_tot_i, linewidth=7, color="blue")
ax[1].set_xlim(1, 1.2e6)
ax[1].set_ylabel(r'$\rho_{d}/\rho_{d,init}$')
ax[1].set_xlabel(r'Time [yr]')
ax[1].set_title(r"Dust")

plt.savefig(os.path.join(pngdir, "cloud_evolution.png"), dpi=300)

print("\nDone.\n")