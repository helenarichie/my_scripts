from hconfig import *

################# hard-coded, fill these in! ###################
date1 = ""
date2 = ""
title = "$n=10^{-2}$ Wind"
xlim = (0, 1)
comp = True
################################################################

date = date1

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5",)
csvdir = os.path.join(basedir, "csv")
pngdir = os.path.join(basedir, "png")

rho_d_tot_16 = []  # total dust density in sim volume for each time
rho_cl_tot_16 = []  # total cloud density in sim volume for each time
T_cl_avg_16 = []  # average T for each time
t_arr_16 = []

print("\nReading in data...")

with open(os.path.join(csvdir, "rho_d_tot.csv")) as f:
    for line in f:
        line = line.split(",")
        rho_d_tot_16.append(float(line[0]))

with open(os.path.join(csvdir, "rho_cl_tot.csv")) as f:
    for line in f:
        line = line.split(",")
        rho_cl_tot_16.append(float(line[0]))

with open(os.path.join(csvdir, "t_arr.csv")) as f:
    for line in f:
        line = line.strip().split(",")
        line = line[0].strip("]").strip("[")
        t_arr_16.append(float(line))

print("\nDone.")

rho_d_tot_16 = np.array(rho_d_tot_16) # why was I multiplying this by 256?
rho_cl_tot_16 = np.array(rho_cl_tot_16)
t_arr_16 = np.array(t_arr_16)

rho_cl_tot_i_16 = rho_cl_tot_16[0]
rho_d_tot_i_16 = rho_d_tot_16[0]

rho_d_tot_32 = []  # total dust density in sim volume for each time
rho_cl_tot_32 = []  # total cloud density in sim volume for each time
T_cl_avg_32 = []  # average T for each time
t_arr_32 = []

if comp:
    date = date2

    basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
    datadir = os.path.join(basedir, "hdf5",)
    csvdir = os.path.join(basedir, "csv")
    pngdir = os.path.join(basedir, "png")

    print("\nReading in data...")

    with open(os.path.join(csvdir, "rho_d_tot.csv")) as f:
        for line in f:
            line = line.split(",")
            rho_d_tot_32.append(float(line[0]))

    with open(os.path.join(csvdir, "rho_cl_tot.csv")) as f:
        for line in f:
            line = line.split(",")
            rho_cl_tot_32.append(float(line[0]))

    with open(os.path.join(csvdir, "t_arr.csv")) as f:
        for line in f:
            line = line.strip().split(",")
            line = line[0].strip("]").strip("[")
            t_arr_32.append(float(line))

    print("\nDone.")

    rho_d_tot_32 = np.array(rho_d_tot_32) # why was I multiplying this by 256?
    rho_cl_tot_32 = np.array(rho_cl_tot_32)
    t_arr_32 = np.array(t_arr_32)

    rho_cl_tot_i_32 = rho_cl_tot_32[0]
    rho_d_tot_i_32 = rho_d_tot_32[0]

print("\nPlotting...")

fig, ax = plt.subplots(2, 1, figsize=(15,22))

ax[0].scatter(t_arr_16/1e6, rho_cl_tot_16/rho_cl_tot_i_16, linewidth=7, color="purple", marker="o", label="16")
if comp:
    ax[0].scatter(t_arr_32/1e6, rho_cl_tot_32/rho_cl_tot_i_32, linewidth=7, color="mediumorchid", marker="o", label="32")
#ax[0].semilogx(t_arr_16, rho_cl_tot_16/rho_cl_tot_i_16, linewidth=7, color="red", linestyle="--", label="16")
#ax[0].semilogx(t_arr_32, rho_cl_tot_32/rho_cl_tot_i_32, linewidth=7, color="red", label="32")
ax[0].set_xlim(xlim)
ax[0].set_ylabel(r'$\rho_{cl}/\rho_{cl,init}$')
ax[0].set_xlabel(r'Time [Myr]')
ax[0].set_title(r"Cloud")
ax[0].legend()

ax[1].scatter(t_arr_16/1e6, rho_d_tot_16/rho_d_tot_i_16, linewidth=7, color="blue", marker="o", label="16")
if comp:
    ax[1].scatter(t_arr_32/1e6, rho_d_tot_32/rho_d_tot_i_32, linewidth=7, color="cornflowerblue", marker="o", label="32")
#ax[1].semilogx(t_arr_16, rho_d_tot_16/rho_d_tot_i_16, linewidth=7, color="blue", linestyle="--", label="16")
#ax[1].semilogx(t_arr_32, rho_d_tot_32/rho_d_tot_i_32, linewidth=7, color="blue", label="32")
ax[1].set_xlim(xlim)
ax[1].set_ylabel(r'$\rho_{d}/\rho_{d,init}$')
ax[1].set_xlabel(r'Time [Myr]')
ax[1].set_title(r"Dust")
ax[1].legend()

plt.suptitle(title)

plt.savefig(os.path.join(pngdir, "cloud_evolution.png"), dpi=300)

print("\nDone.\n")