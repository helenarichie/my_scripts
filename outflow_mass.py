from hconfig import *

from csv import writer

######### hard-coded values! #########
date = "2023-03-07"
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5/full/")
pngdir = os.path.join(basedir, "png/")
csvdir = os.path.join(basedir, "csv/")
cat = False
cutoff = 0.05*1e-24  # 5% of initial density
######################################

data = ReadHDF5(datadir, fnum=0, nscalar=1, cat=cat)
head = data.head
nx, ny, nz = head["dims"]
dx = data.dx_cgs()[0]

# x, y, z
areas = [ny*nz*dx**2, nx*nz*dx**2, nx*ny*dx**2]
indices = ["0", "-1"]
fluids = ["gas", "dust"]

def calc_mass_loss_rate(rho, v, area):
        return np.sum(rho * np.abs(v)) * area

if cat:
    files = glob.glob(os.path.join(datadir, "*.h5"))
else:
    files = glob.glob(os.path.join(datadir, "*.h5.0"))

f = open(os.path.join(csvdir, "outflow_tot_cl.csv"), "w")
f.close()
f = open(os.path.join(csvdir, "outflow_dense_cl.csv"), "w")
f.close()
f = open(os.path.join(csvdir, "masses_cl.csv"), "w")
f.close()
f = open(os.path.join(csvdir, "outflow_tot_dust.csv"), "w")
f.close()
f = open(os.path.join(csvdir, "outflow_dense_dust.csv"), "w")
f.close()
f = open(os.path.join(csvdir, "masses_dust.csv"), "w")
f.close()

def get_rates(d, cutoff):
    rates_tot, rates_dense, masses = np.zeros(6), np.zeros(6), np.zeros(6)
    for i, area in enumerate(areas):
        for index in indices:
            if (index == "0") and (i == 0):
                d_i, v_i = d[0, :, :], data.vx_cgs()[0][0, :, :]
                mask = np.logical_and(d_i >= cutoff, v_i < 0) # mask for high-density, outflowing gas
                rates_tot[0] = calc_mass_loss_rate(d_i, v_i, area)
                rates_dense[0] = calc_mass_loss_rate(d_i[mask], v_i[mask], area)
                masses[0] = np.sum(d_i[mask]) * d_i[mask].size * dx**2 # total density times volume of face

            if (index == "0") and (i == 1):
                d_i, v_i = d[:, 0, :], data.vy_cgs()[0][:, 0, :]
                mask = np.logical_and(d_i >= cutoff, v_i < 0)
                rates_tot[1] = calc_mass_loss_rate(d_i, v_i, area)
                rates_dense[1] = calc_mass_loss_rate(d_i[mask], v_i[mask], area)
                masses[1] = np.sum(d_i[mask]) * d_i[mask].size * dx**2

            if (index == "0") and (i == 2):
                d_i, v_i = d[:, :, 0], data.vz_cgs()[0][:, :, 0]
                mask = np.logical_and(d_i >= cutoff, v_i < 0)
                rates_tot[2] = calc_mass_loss_rate(d_i, v_i, area)
                rates_dense[2] = calc_mass_loss_rate(d_i[mask], v_i[mask], area)
                masses[2] = np.sum(d_i[mask]) * d_i[mask].size * dx**2

            if (index == "-1") and (i == 0):
                d_i, v_i = d[-1, :, :], data.vx_cgs()[0][-1, :, :]
                mask = np.logical_and(d_i >= cutoff, v_i > 0)
                rates_tot[3] = calc_mass_loss_rate(d_i, v_i, area)
                rates_dense[3] = calc_mass_loss_rate(d_i[mask], v_i[mask], area)
                masses[3] = np.sum(d_i[mask]) * d_i[mask].size * dx**2

            if (index == "-1") and (i == 1):
                d_i, v_i = d[:, -1, :], data.vy_cgs()[0][:, -1, :]
                mask = np.logical_and(d_i >= cutoff, v_i > 0)
                rates_tot[4] = calc_mass_loss_rate(d_i, v_i, area)
                rates_dense[4] = calc_mass_loss_rate(d_i[mask], v_i[mask], area)
                masses[4] = np.sum(d_i[mask]) * d_i[mask].size * dx**2

            if (index == "-1") and (i == 2):
                d_i, v_i = d[:, :, -1], data.vz_cgs()[0][:, :, -1]
                mask = np.logical_and(d_i >= cutoff, v_i > 0)
                rates_tot[5] = calc_mass_loss_rate(d_i, v_i, area)
                rates_dense[5] = calc_mass_loss_rate(d_i[mask], v_i[mask], area)
                masses[5] = np.sum(d_i[mask]) * d_i[mask].size * dx**2
    
    return rates_tot, rates_dense, masses

for i in range(0, len(files)):
    data = ReadHDF5(datadir, fnum=i, nscalar=1, cat=cat)
    head = data.head
    conserved = data.conserved
    nx, ny, nz = head["dims"]
    gamma = head["gamma"]
    dx = data.dx_cgs()[0]

    rates_tot, rates_dense, masses = get_rates(data.d_cgs()[0], cutoff)

    with open(os.path.join(csvdir, "outflow_tot_cl.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(rates_tot)
        f.close()

    with open(os.path.join(csvdir, "outflow_dense_cl.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(rates_dense)
        f.close()

    with open(os.path.join(csvdir, "masses_cl.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(masses)
        f.close()

    rates_tot, rates_dense, masses = get_rates(conserved["scalar0"][0] * head["density_unit"], 1e-2*cutoff)

    with open(os.path.join(csvdir, "outflow_tot_dust.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(rates_tot)
        f.close()

    with open(os.path.join(csvdir, "outflow_dense_dust.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(rates_dense)
        f.close()

    with open(os.path.join(csvdir, "masses_dust.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(masses)
        f.close()