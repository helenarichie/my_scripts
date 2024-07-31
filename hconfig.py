
import sys
if __file__ == "/ix/eschneider/helena/code/my_scripts/hconfig.py":
    sys.path.insert(0, "/ix/eschneider/helena/code/analysis_scripts")
if __file__ == "/Users/helenarichie/GitHub/my_scripts/hconfig.py":
    sys.path.insert(0, "/Users/helenarichie/GitHub/analysis_scripts/")
if __file__ == "/ccs/home/helenarichie/code/my_scripts/hconfig.py":
    sys.path.insert(0, "/ccs/home/helenarichie/code/analysis_scripts")
from read_hdf5 import ReadHDF5
import os
import glob
import h5py

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.family': 'sans serif'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.fontset':'stixsans'})
plt.rcParams.update({'axes.linewidth': 1.5})
plt.rcParams.update({'xtick.direction':'in'})
plt.rcParams.update({'xtick.major.size': 7})
plt.rcParams.update({'xtick.major.width': 1.5})
plt.rcParams.update({'xtick.minor.size': 5})
plt.rcParams.update({'xtick.minor.width': 2 })
plt.rcParams.update({'ytick.direction':'in'})
plt.rcParams.update({'ytick.major.size': 7})
plt.rcParams.update({'ytick.major.width': 1.5})
plt.rcParams.update({'ytick.minor.size': 5})
plt.rcParams.update({'ytick.minor.width': 2 })
plt.rcParams.update({'lines.linewidth': 4})

colors = ["lightskyblue", "violet", "forestgreen", "darkviolet", "thistle", "cadetblue", 
            "palegoldenrod", "darksalmon", "indigo"]
yr_in_s = 3.154e7 # s
MP = 1.672622e-24 # g
KB = 1.3807e-16 # cm^2 g s^-2 K^-1

# import warnings
# warnings.filterwarnings("ignore")
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time

plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': 10})
# plt.style.use('dark_background')
