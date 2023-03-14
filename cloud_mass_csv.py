from hconfig import *
from csv import writer

################# hard-coded, fill these in! ###################
date = "2023-03-09"
write_csv = True  # do you need to read data from HDF5 files?
rho_cl_i = 1e-24  # n = 1, needed to index cloud material
cat = False
################################################################

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5/raw/full/",)
csvdir = os.path.join(basedir, "csv")
pngdir = os.path.join(basedir, "png")

t_arr = []
wind_init = []

def write_csv(path, fnum=None, nscalar=0, cat=True):
    data = ReadHDF5(path, fnum=fnum, nscalar=nscalar, cat=cat)
    head = data.head
    conserved = data.conserved

    d_gas = data.d_cgs()
    d_dust = conserved["scalar0"] * head["density_unit"]
    T = data.T()
    t = data.t_cgs() / yr_in_s

    rho_d_tot_i = None
    rho_cl_tot_i = None
    T_cl_avg_i = None
    for i, dens in enumerate(d_dust):
        rho_d_tot_i = np.sum(d_dust[i][d_dust[i]>0])
        rho_cl_tot_i = np.sum(d_gas[i][d_gas[i]>=1/3*rho_cl_i])
        T_cl_avg_i = np.average(T[i][d_gas[i]>=1/3*rho_cl_i])

    with open(os.path.join(csvdir, "rho_d_tot.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow([rho_d_tot_i])
        f.close()

    with open(os.path.join(csvdir, "rho_cl_tot.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow([rho_cl_tot_i])
        f.close()

    with open(os.path.join(csvdir, "T_cl_avg.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow([T_cl_avg_i])
        f.close()

    with open(os.path.join(csvdir, "t_arr.csv"), "a") as f:
        writer_obj = writer(f)
        writer_obj.writerow([t[0]])
        f.close()

for i in range(0, 201):
    write_csv(datadir, fnum=i, nscalar=1, cat=cat)