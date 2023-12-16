import sys
sys.path.insert(0, "/ix/eschneider/helena/code/my_scripts")
from hconfig import *

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

##################################################################
date = "2023-11-03"
datestr = "1103"
cat = True
fontsize = 20
labelpad = 12
tickwidth = 1
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': fontsize})
fnum_1103 = 304
rho_cl_init = 10*0.6*MP * density_conversion
##################################################################

fnum = fnum_1103
fmax = 600

##################################################################
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
fulldir = os.path.join(basedir, "hdf5/full/")
pngdir = os.path.join(basedir, "png/")
##################################################################

for fnum in range(fnum, fnum+1):
    data = ReadHDF5(fulldir, cat=cat, fnum=fnum, dust=True)
    head = data.head
    conserved = data.conserved
    dx = head["dx"][0]
    nx, ny, nz = head["dims"]

    d_gas = conserved["density"][0]
    d_dust = conserved["dust_density"][0]

    gas_sum = []
    dust_sum = []
    for i, d in enumerate(d_gas):
        # only get values for gas and dust that are in the cloud
        gas_sum.append(np.sum(d_gas[i,:,:][d_gas[i,:,:]>=rho_cl_init/3]))
        dust_sum.append(np.sum(d_dust[i,:,:][d_gas[i,:,:]>=rho_cl_init/3]))

    gas_sum = np.array(gas_sum)
    dust_sum = np.array(dust_sum)

    dust_to_gas = dust_sum/gas_sum

    x_arr = np.zeros(nx)
    for i, x in enumerate(x_arr):
        if i > 0:
            x_arr[i] = x_arr[i-1] + dx
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 7))

    ax.plot(x_arr, dust_to_gas, linewidth=4, c="#f16948", label="in box")
    ax.set_xlabel(r"Distance$~[kpc]$", fontsize=fontsize)
    ax.set_ylabel(r"Dust-to-gas ratio", fontsize=fontsize)
    #ax.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True)
    #ax.set_xlim(xmin/1e6, xmax/1e6)
    #ax.set_xticks(np.linspace(0, np.amax(t_arr[t_arr<=tmax]/1e6), 5).round(1))
    #ax.set_yticks(np.linspace(0, np.amax(mass_dust_tot[t_arr<=tmax]), 5).round(2))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.savefig(pngdir + f"{fnum}_dust_to_gas_{datestr}.png", dpi=300)