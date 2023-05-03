from hconfig import *
from csv import writer

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

################# hard-coded, fill these in! ###################
date = "2023-04-26"
cutoff = 0.05*1e-24*density_conversion # 5% of initial density, M_sun/kpc^3
cat = True
istart = 301
################################################################

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/" # crc
datadir = os.path.join(basedir, "hdf5/full/")
csvdir = os.path.join(basedir, "csv/")

data = ReadHDF5(datadir, fnum=0, dust=True, cat=cat)
head = data.head
nx, ny, nz = head["dims"]
dx = head["dx"] # kpc

# x, y, z
areas = [ny*nz*dx**2, nx*nz*dx**2, nx*ny*dx**2]
indices = [0, -1]
fluids = ["gas", "dust"]

if cat:
    files = glob.glob(os.path.join(datadir, "*.h5"))
else:
    files = glob.glob(os.path.join(datadir, "*.h5.0"))

if istart == 0:
    # create files to write quantities for cloud (dense gas only), total gas, and dust
    f = open(os.path.join(csvdir, "rate_cl.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "mass_cl.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "rate_gas.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "mass_gas.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "mass_dust.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "rate_dust.csv"), "w")
    f.close()

def calc_mass_loss_rate(rho, v, area):
    return np.sum(rho * np.abs(v) * area)

def get_rates(d, cutoff):
    rates = np.zeros(6)
    for i, area in enumerate(areas):
        for index in indices:
            if (index == 0) and (i == 0):
                d_i, v_i = d[0, :, :], conserved["momentum_x"][0][0, :, :] / conserved["density"][0][0, :, :]
                mask = np.logical_and(d_i >= cutoff, v_i < 0) # mask for high-density, outflowing gas
                rates[0] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == 0) and (i == 1):
                d_i, v_i = d[:, 0, :], conserved["momentum_y"][0][:, 0, :] / conserved["density"][0][:, 0, :]
                mask = np.logical_and(d_i >= cutoff, v_i < 0)
                rates[1] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == 0) and (i == 2):
                d_i, v_i = d[:, :, 0], conserved["momentum_z"][0][:, :, 0] / conserved["density"][0][:, :, 0]
                mask = np.logical_and(d_i >= cutoff, v_i < 0)
                rates[2] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == -1) and (i == 0):
                d_i, v_i = d[-1, :, :], conserved["momentum_x"][0][-1, :, :] / (conserved["density"][0][-1, :, :])
                mask = np.logical_and(d_i >= cutoff, v_i > 0)
                rates[3] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == -1) and (i == 1):
                d_i, v_i = d[:, -1, :], conserved["momentum_y"][0][:, -1, :] / conserved["density"][0][:, -1, :]
                mask = np.logical_and(d_i >= cutoff, v_i > 0)
                rates[4] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == -1) and (i == 2):
                d_i, v_i = d[:, :, -1], conserved["momentum_z"][0][:, :, -1] / conserved["density"][0][:, :, -1]
                mask = np.logical_and(d_i >= cutoff, v_i > 0)
                rates[5] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)
    
    return rates

def get_masses(d, cutoff):
    masses = np.zeros(6)
    for i, area in enumerate(areas):
        for index in indices:
            if (index == 0) and (i == 0):
                d_i = d[0, :, :]
                mask = (d_i >= cutoff)
                masses[0] = np.sum(d_i[mask]) * dx**3 # total density times total volume of cells

            if (index == 0) and (i == 1):
                d_i = d[:, 0, :]
                mask = (d_i >= cutoff)
                masses[1] = np.sum(d_i[mask]) * dx**3

            if (index == 0) and (i == 2):
                d_i = d[:, :, 0]
                mask = (d_i >= cutoff)
                masses[2] = np.sum(d_i[mask]) * dx**3

            if (index == -1) and (i == 0):
                d_i = d[-1, :, :]
                mask = (d_i >= cutoff)
                masses[3] = np.sum(d_i[mask]) * dx**3

            if (index == -1) and (i == 1):
                d_i = d[:, -1, :]
                mask = (d_i >= cutoff)
                masses[4] = np.sum(d_i[mask]) * dx**3

            if (index == -1) and (i == 2):
                d_i = d[:, :, -1]
                mask = (d_i >= cutoff)
                masses[5] = np.sum(d_i[mask]) * dx**3
    
    return masses

# loop through hdf5 files to write out rates and masses for cloud, gas, and dust
for i in range(istart, len(files)):
    # read in data
    data = ReadHDF5(datadir, fnum=i, dust=True, cat=cat)
    head = data.head
    conserved = data.conserved
    dx = head["dx"][0]

    # calculate and write rates and masses for cloud
    rates = get_rates(conserved["density"][0], cutoff)
    masses = get_masses(conserved["density"][0], cutoff)

    with open(os.path.join(csvdir, "rate_cl.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(rates)
        f.close()

    with open(os.path.join(csvdir, "mass_cl.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(masses)
        f.close()

    # calculate and write rates and masses for gas
    rates = get_rates(conserved["density"][0], 0)
    masses = get_masses(conserved["density"][0], 0)

    with open(os.path.join(csvdir, "rate_gas.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(rates)
        f.close()

    with open(os.path.join(csvdir, "mass_gas.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(masses)
        f.close()

    # calculate and write masses for dust
    masses = get_masses(conserved["dust_density"][0], 1)

    with open(os.path.join(csvdir, "mass_dust.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(masses)
        f.close()

    # calculate and write masses for dust
    rates = get_rates(conserved["dust_density"][0], 1)

    with open(os.path.join(csvdir, "rate_dust.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow(rates)
        f.close()