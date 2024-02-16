from hconfig import *
from csv import writer

#################################
date = "2024-02-07"
density_cl_init = 1e-23
write, plot = True, True
ns = 0
ne = 4000
cat = True
#################################

########### location ############
crc = True
frontier = False
#################################

########## data type ############
debugging = False
cloud_wind = True
testing = False
#################################

########## specify slice and png directories ############
if crc:
    if debugging:
        basedir = f"/ix/eschneider/helena/data/debugging/{date}/"
    if cloud_wind:
        basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
    if testing:
        basedir = f"/ix/eschneider/helena/data/testing/{date}/"

if frontier:
    basedir = f"/lustre/orion/ast181/scratch/helenarichie/{date}/"

if cat:
    datadir = os.path.join(basedir, "hdf5/slice/")
else:
    datadir = os.path.join(basedir, "hdf5/raw/")

csvdir = os.path.join(basedir, "csv/")
pngdir = os.path.join(basedir, "png/")
#########################################################

if write:

    if ns == 0:
        f = open(os.path.join(csvdir, "avg_v.csv"), "w")
        f.close()

    for i in range(ns, ne+1):
        # read in data
        data = ReadHDF5(datadir, fnum=i, slice="xy", cat=cat)

        head = data.head
        conserved = data.conserved

        wh_cl = np.where(data.d_cgs() >= density_cl_init/3)
        vx = np.average(data.vx_cgs()[wh_cl])
        vy = np.average(data.vy_cgs()[wh_cl])
        vz = np.average(data.vz_cgs()[wh_cl])
        t = head["t"][0]

        with open(os.path.join(csvdir, "avg_v.csv"), "a") as f:
            writer_obj = writer(f)
            writer_obj.writerow([t, vx, vy, vz])
            f.close()

if plot:
    t, vx, vy, vz = [], [], [], []
    with open(os.path.join(csvdir, "avg_v.csv"), "r") as f:
        for line in f:
            line = line.split(",")
            t.append(float(line[0]))
            vx.append(float(line[1]))
            vy.append(float(line[2]))
            vz.append(float(line[3]))
    t, vx, vy, vz = np.array(t), np.array(vx), np.array(vy), np.array(vz)

    plt.plot(t*1e-3, vx*1e-5, linewidth=3, c=colors[0])
    # plt.plot(t, vy*1e-5)
    # plt.plot(t, vz*1e-5)
    plt.xlabel("Time [Myr]")
    plt.ylabel("$v_{cl,x}~[km/s]$")

    plt.tight_layout()
    plt.savefig(os.path.join(pngdir, "v_avg.png"))
