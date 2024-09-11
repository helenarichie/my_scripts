import sys
sys.path.insert(0, "../github/my_scripts")
from hconfig import *
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl

import numpy as np
import os
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 25})
from mpl_toolkits.axes_grid1 import make_axes_locatable

# McKinnon et al. (2017)
def calc_tau_sp(n, T, a):
    YR_IN_S = 3.154e7;
    a1 = a/0.1; # dust grain size in units of 0.1 micrometers
    d0 = n / (6e-4); # gas density in units of 10^-27 g/cm^3
    T_0 = 2e6; # K
    omega = 2.5;
    A = 0.17e9 * YR_IN_S; # Gyr in s

    tau_sp = A * (a1/d0) * ((T_0/T)**omega + 1); # s

    return tau_sp

molecular = [100, 15]
CNM = [20, 75]
WNM = [3000, 0.5]
phases = [molecular, CNM, WNM]
labels = ["molecular", "CNM", "WNM"]

a = np.array([0.001, 0.01, 0.1, 1])

fig = plt.figure(figsize=(7, 6.5))
for i, phase in enumerate(phases):
    plt.semilogy(a*0.1, calc_tau_sp(phase[0], phase[1], a)/(yr_in_s*1e6), label=labels[i])
plt.yticks(np.logspace(8, 13, 5).round(1))
plt.legend(fontsize=20)
plt.tight_layout()
plt.xlabel("grain size [micron]")
plt.ylabel("Sputtering Time [Myr]")
plt.savefig("sputtering_radius")
