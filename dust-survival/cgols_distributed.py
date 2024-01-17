import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("/Users/helenarichie/GitHub/")
import galacticwindcloud.plot.code_plot as cp
import galacticwindcloud.utils as utils

# Locate the cloud catalog data
import galacticwindcloud.data as data_module 
data_dir = data_module.__path__[0]

def load(filename):
    return utils.z_select_dx(np.load(filename), 2048, 10/2048., 1.0) # 0.5 instead of 1.0?
    # return utils.cone_select(np.load(filename), 2048)

def make_histogram(array,res,color,label):
    radii = cp.radius_pc(res, array[:,1])s
    sum_density = array[:,2]
    mass_m_sun = sum_density * cp.cell_vol_kpc(res)    
    log_radii = np.log10(radii)
    x = plt.hist(log_radii,weights=mass_m_sun,color=color,histtype='step',label=label,bins=20, range=(0.5,2.5))
    # x = plt.hist(log_radii,color=color,histtype='step',label=label,bins=20,range=(0.5,2.5))
    plt.xlabel('$log_{10}(radius)$ [pc]')
    plt.ylabel('Mass in bin [$\mathrm{M}_\odot$]') 
    # plt.ylabel('$\mathrm{N}_\mathrm{clouds}$')
    return x

def load(filename):
    return utils.z_select_dx(np.load(filename), 2048, 10/2048., 0.5) # 0.5 instead of 1.0?

fd = '{}/2048{}/30.{}.2048.merged.npy'.format(data_dir,'dist','dist')

mass, radii, dontcare = make_histogram(load(fd), 2048, 'k', 'Distributed')

radii = 10**radii

r_mid = []
for i, r in enumerate(radii):
    if (i+1) < len(radii):
        r_mid.append((radii[i]+radii[i+1])/2)
r_mid = np.array(r_mid)

r_crit = 58
r_smooth = 245

wh_dest = np.where(r_mid <= r_crit)
wh_disr = np.where(np.logical_and(r_mid > r_crit, r_mid < r_smooth))
wh_surv = np.where(r_mid >= r_smooth)
print(f"Total cloud mass:     {np.sum(mass):.5e}")
print(f"Destroyed cloud mass: {np.sum(mass[wh_dest]):.5e}")
print(f"Disrupted cloud mass: {np.sum(mass[wh_disr]):.5e}")
print(f"Survived cloud mass:  {np.sum(mass[wh_surv]):.5e}")

def critical_radius(p):
    M = 4  # Mach number
    chi3 = 1  # density contrast in units of 1e3
    Tcl4 = 3  # cloud temperature in units of 1e4
    Tmincool = 3e4  # minimum cool phase temperature
    p3 = 1e3  # pressure in units of p/kB/(1e3K cm−3),
    lambd = 0.41  # cooling rate in units of 1e−21.4
    return 618*M*chi3**1.5*Tcl4**1.5*Tmincool/((p/1e3)*lambd*10**4.25)
