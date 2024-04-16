import sys
sys.path.insert(0, "")
import plot_generic
import os
import hconfig

date = ""
basedir = f""
datadir = os.path.join(basedir, "hdf5/raw/")
ns = 0
ne = 1
cat = False
ftype = ""
pngdir = os.path.join(basedir, "png/" + str(ftype) + "/")
fields = [""]
shapes = [""]  # xy, xz, yz
xlabels = ["$x~$[kpc]"]
ylabels = ["$y~$[kpc]"]
clabels = [""]
vlims = [[]]

plot_generic.plot_generic(datadir, pngdir, ns, ne, cat, ftype, fields, shapes, xlabels, ylabels, clabels, vlims, print_keys=False)


"""
dust --- r"$\mathrm{log}_{10}(\rho_{dust})$ [$\mathrm{M}_\odot\mathrm{kpc}^{-3}$]"
"""