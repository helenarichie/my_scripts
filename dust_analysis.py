from hconfig import *

################# hard-coded, fill these in! ###################
date = "2023-03-09"
write_csv = True  # do you need to read data from HDF5 files?
rho_cl_i = 1e-24  # n = 1, needed to index cloud material
################################################################

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5/full",)
csvdir = os.path.join(basedir, "csv")
pngdir = os.path.join(basedir, "png")

rho_d_tot = []  # total dust density in sim volume for each time
rho_cl_tot = []  # total cloud density in sim volume for each time
T_cl_avg = []  # average T for each time
t_arr = []

# read data in from HDF5 files and write it out to csv files
if write_csv:
    # read in full-volume data from HDF5 files with one scalar field
    data = ReadHDF5(datadir, nscalar=1, cat=False)
    head = data.head
    conserved = data.conserved

    start = time.time()
    print("Loading gas densities...")
    d_gas = data.d_cgs()
    print(f"Finished --- {round(time.time() - start, 3)} s")

    start = time.time()
    print("\nLoading dust densities...")
    d_dust = conserved["scalar0"] * head["density_unit"]
    print(f"Finished --- {round(time.time() - start, 3)} s")

    start = time.time()
    print("\nLoading temperatures...")
    T = data.T()
    print(f"Finished --- {round(time.time() - start, 3)} s")

    gamma = head["gamma"]
    t_arr = data.t_cgs() / yr_in_s

    start = time.time()
    print("\nGetting cloud and dust densities...")
    for i, dens in enumerate(d_dust):
        rho_d_tot.append(np.sum(d_dust[i][d_dust[i]>0]))
        rho_cl_tot.append(np.sum(d_gas[i][d_gas[i]>=1/3*rho_cl_i]))
        T_cl_avg.append(np.average(T[i][d_gas[i]>=1/3*rho_cl_i]))
    print(f"Finished --- {round(time.time() - start, 3)} s")

    np.savetxt(os.path.join(csvdir, "rho_d_tot.csv"), rho_d_tot, delimiter=",")
    np.savetxt(os.path.join(csvdir, "rho_cl_tot.csv"), rho_cl_tot, delimiter=",")
    np.savetxt(os.path.join(csvdir, "T_cl_avg.csv"), T_cl_avg, delimiter=",")
    np.savetxt(os.path.join(csvdir, "t_arr.csv"), t_arr, delimiter=",")

# if the csv files were already saved
else: 
    with open(os.path.join(csvdir, "rho_d_tot.csv")) as f:
        for line in f:
            line = line.split(",")
            rho_d_tot.append(float(line[0]))

    with open(os.path.join(csvdir, "rho_cl_tot.csv")) as f:
        for line in f:
            line = line.split(",")
            rho_cl_tot.append(float(line[0]))
 
    with open(os.path.join(csvdir, "T_cl_avg.csv")) as f:
        for line in f:
            line = line.split(",")
            T_cl_avg.append(float(line[0]))

    with open(os.path.join(csvdir, "t_arr.csv")) as f:
        for line in f:
            line = line.split(",")
            t_arr.append(float(line[0]))

rho_d_tot = np.array(rho_d_tot) # why was I multiplying this by 256?
rho_cl_tot = np.array(rho_cl_tot)
T_cl_avg = np.array(T_cl_avg)
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