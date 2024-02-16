from hconfig import *
from csv import writer

def calc_mass_loss_rate(rho, v, area):
    return np.sum(rho * np.abs(v) * area)

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

################################################################
date = "2024-02-11"
rho_cl_i = 1e-23  # needed to index cloud material
cutoff = rho_cl_i*density_conversion/3 # M_sun/kpc^3
ns = 0
ne = 1200
DE = True
SCALAR = True
################################################################

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5/edges/")
csvdir = os.path.join(basedir, "csv/")

if ns == 0:
    f = open(os.path.join(csvdir, "rate_cloud.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "rate_dust.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "time_output.csv"), "w")
    f.close()

f = h5py.File(os.path.join(datadir, "0_edges.h5"))
head = f.attrs
nx, ny, nz = head["dims"]
dx, dy, dz = head["dx"]

files = glob.glob(os.path.join(datadir, "*.h5"))

n_faces = 6
ns_gas = 2
ns_dust = 6
ms = ["mz_minus_xy", "my_minus_xz", "mx_minus_yz", "mz_plus_xy", "my_plus_xz", "mx_plus_yz"]

# loop through hdf5 files to write out rates and masses for cloud and dust
for n in range(ns, len(files)):
    # read in data
    fh5 = h5py.File(os.path.join(datadir, str(n)+"_edges.h5"))
    keys_gas = list(fh5.keys())[ns_gas*n_faces:(ns_gas+1)*n_faces]
    keys_dust = list(fh5.keys())[ns_dust*n_faces:(ns_dust+1)*n_faces]
    head = fh5.attrs
    nx, ny, nz = head["dims"]
    dx, dy, dz = head["dx"]
    t = head["t"]

    with open(os.path.join(csvdir, "t_arr.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow([t])
        f.close()

    rates_cloud, rates_dust = [], []
    for i, key in enumerate(keys_gas):
        velocity, mask = None, None

        gas = fh5[keys_gas[i]][()]
        dust = fh5[keys_dust[i]][()]

        if i <= 2:
            velocity = fh5[ms[i]] / gas
            mask = np.logical_and(gas >= cutoff, velocity < 0)
        if i > 2:
            velocity = fh5[ms[i]] / gas
            mask = np.logical_and(gas >= cutoff, velocity > 0)  

        rates_cloud.append(calc_mass_loss_rate(gas[mask], velocity[mask], dx**2))

        if i <= 2:
            velocity = fh5[ms[i]] / gas
            mask = np.logical_and(dust >= 0, velocity < 0)
        if i > 2:
            velocity = fh5[ms[i]] / gas
            mask = np.logical_and(dust >= 0, velocity > 0)
        
        rates_dust.append(calc_mass_loss_rate(dust[mask], velocity[mask], dx**2))

    with open(os.path.join(csvdir, "rate_cloud.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(rates_cloud)
        f.close()

    with open(os.path.join(csvdir, "rate_dust.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(rates_dust)
        f.close()

    with open(os.path.join(csvdir, "time_output.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(t)
        f.close()
