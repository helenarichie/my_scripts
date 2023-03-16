from hconfig import *

from csv import writer

######### hard-coded values! #########
date = "2023-03-14"
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5/full/")
pngdir = os.path.join(basedir, "png/")
csvdir = os.path.join(basedir, "csv/")
cat = False
cutoff = 1/3*1e-24
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
        return np.sum(rho * v) * area

if cat:
    files = glob.glob(os.path.join(datadir, "*.h5"))
else:
    files = glob.glob(os.path.join(datadir, "*.h5.0"))

f = open(os.path.join(csvdir, "outflow_cloud.csv"), "w")
f.close()
f = open(os.path.join(csvdir, "flux_cloud.csv"), "w")
f.close()
f = open(os.path.join(csvdir, "outflow_dust.csv"), "w")
f.close()
f = open(os.path.join(csvdir, "flux_dust.csv"), "w")
f.close()

def get_rates(d, cutoff):
    rates, fluxes = np.zeros(6), np.zeros(6)
    for i, area in enumerate(areas):
        for index in indices:
            if (index == "0") and (i == 0):
                rates[0] = calc_mass_loss_rate(d[0, :, :], data.vx_cgs()[0][0, :, :], area)
                flux = d[0, :, :][d[0, :, :]>=cutoff]
                fluxes[0] = np.sum(flux) * flux.size * dx**2
            if (index == "0") and (i == 1):
                rates[1] = calc_mass_loss_rate(d[:, 0, :], data.vy_cgs()[0][:, 0, :], area)
                flux = d[:, 0, :][d[:, 0, :]>=cutoff]
                fluxes[1] = np.sum(flux) * flux.size * dx**2
            if (index == "0") and (i == 2):
                rates[2] = calc_mass_loss_rate(d[:, :, 0], data.vz_cgs()[0][:, :, 0], area)
                flux = d[:, :, 0][d[:, :, 0]>=cutoff]
                fluxes[2] = np.sum(flux) * flux.size * dx**2
            if (index == "-1") and (i == 0):
                rates[3] = calc_mass_loss_rate(d[-1, :, :], data.vx_cgs()[0][-1, :, :], area)
                flux = d[-1, :, :][d[-1, :, :]>=cutoff]
                fluxes[3] = np.sum(flux) * flux.size * dx**2
            if (index == "-1") and (i == 1):
                rates[4] = calc_mass_loss_rate(d[:, -1, :], data.vy_cgs()[0][:, -1, :], area)
                flux = d[:, -1, :][d[:, -1, :]>=cutoff]
                fluxes[4] = np.sum(flux) * flux.size * dx**2
            if (index == "-1") and (i == 2):
                rates[5] = calc_mass_loss_rate(d[:, :, -1], data.vz_cgs()[0][:, :, -1], area)
                flux = d[:, :, -1][d[:, :, -1]>=cutoff]
                fluxes[5] = np.sum(flux) * flux.size * dx**2
    
    return rates, fluxes


for i in range(0, len(files)):
    data = ReadHDF5(datadir, fnum=i, nscalar=1, cat=cat)
    head = data.head
    conserved = data.conserved
    nx, ny, nz = head["dims"]
    gamma = head["gamma"]
    dx = data.dx_cgs()[0]

    outflow_rates, fluxes = get_rates(data.d_cgs()[0], cutoff)

    with open(os.path.join(csvdir, "outflow_cloud.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(outflow_rates)
        f.close()

    with open(os.path.join(csvdir, "flux_cloud.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(fluxes)
        f.close()

    outflow_rates, fluxes = get_rates(conserved["scalar0"][0] * head["density_unit"], 1e-2*cutoff)

    with open(os.path.join(csvdir, "outflow_dust.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(outflow_rates)
        f.close()

    with open(os.path.join(csvdir, "flux_dust.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(fluxes)
        f.close()