from hconfig import *
from csv import writer

################# hard-coded, fill these in! ###################
date = "2023-03-20"
cutoff = 0.05*1e-24  # 5% of initial density
cat = True
################################################################

data = ReadHDF5(datadir, fnum=0, nscalar=1, cat=cat)
head = data.head
nx, ny, nz = head["dims"]
dx = data.dx_cgs()[0]

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/" # crc
datadir = os.path.join(basedir, "hdf5/full/")
csvdir = os.path.join(basedir, "csv")

t_arr = []

f = open(os.path.join(csvdir, "rho_d_tot.csv"), "w")
f.close()
f = open(os.path.join(csvdir, "rho_g_tot.csv"), "w")
f.close()
f = open(os.path.join(csvdir, "t_arr.csv"), "w")
f.close()