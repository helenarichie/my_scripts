"""
/ix/eschneider/helena/code/my_scripts/
/ccs/home/hrichie/code/my_scripts/
"""

import sys
sys.path.insert(0, "")
import plot_generic
import os
import hconfig

date = ""
ftype = ""
basedir = f""
cat = False
if cat:
    datadir = os.path.join(basedir, "hdf5", ftype)
else:
    datadir = os.path.join(basedir, "hdf5", "raw")
ns = 0
ne = 1
pngdir = os.path.join(basedir, "png/" + str(ftype) + "/")
fields = [""]
vlims = [[]]

plot_generic.plot_generic(datadir, pngdir, ns, ne, cat, ftype, fields, vlims)


"""
dust --- r"$\mathrm{log}_{10}(\rho_{dust})$ [$\mathrm{M}_\odot\mathrm{kpc}^{-3}$]"
"""