import sys
sys.path.insert(0, "/ix/eschneider/helena/code/my_scripts")
from hconfig import *

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

##################################################################
date_1107 = "2023-11-07"
date_1103 = "2023-11-03"
datestr = "1107"
cat = True
fontsize = 20
labelpad = 12
tickwidth = 1
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': fontsize})
fnum_1103 = 304
fnum_1107 = 304
rho_cl_init = 10*0.6*MP * density_conversion
##################################################################

fnum = fnum_1103
fnum = fnum_1107
fmax = 600

fnums = [fnum_1103, fnum_1107]

##################################################################
basedir_1103 = f"/ix/eschneider/helena/data/cloud_wind/{date_1103}/"
basedir_1107 = f"/ix/eschneider/helena/data/cloud_wind/{date_1107}/"
pngdir = os.path.join(basedir_1107, "png/dtg/")
##################################################################

fulldirs = [os.path.join(basedir_1103, "hdf5/full/"), os.path.join(basedir_1107, "hdf5/full/")]
labels = ["disrupted", "survived"]
colors = ["#2e5e79", "#b63253"]

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 7))

for i, fnum in enumerate(fnums):
    print(i)
    data = ReadHDF5(fulldirs[i], cat=cat, fnum=fnum, dust=True)
    head = data.head
    conserved = data.conserved
    dx = head["dx"][0]
    nx, ny, nz = head["dims"]

    d_gas = conserved["density"][0]
    d_dust = conserved["dust_density"][0]

    gas_sum = []
    dust_sum = []
    for k, d in enumerate(d_gas):
        # only get values for gas and dust that are in the cloud
        gas_sum.append(np.sum(d_gas[k,:,:][d_gas[k,:,:]>=rho_cl_init/10]))
        dust_sum.append(np.sum(d_dust[k,:,:][d_gas[k,:,:]>=rho_cl_init/10]))

    gas_sum = np.array(gas_sum)
    dust_sum = np.array(dust_sum)

    dust_to_gas = dust_sum/gas_sum

    x_arr = np.zeros(nx)
    for j, x in enumerate(x_arr):
        if j > 0:
            x_arr[j] = x_arr[j-1] + dx

    wh_real = np.isnan(dust_to_gas, where=False)

    # , c="#f16948"
    ax.plot(x_arr, dust_to_gas, linewidth=4, label=labels[i], color=colors[i])
    ax.set_xlabel(r"r$~[kpc]$", fontsize=fontsize)
    ax.set_ylabel(r"Dust-to-gas ratio", fontsize=fontsize)
    ax.set_xlim(1.7, np.amax(x_arr[wh_real]))
    #ax.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True)
    ax.set_xticks(np.linspace(1.7, np.amax(x_arr[wh_real]), 5).round(1))
    #ax.set_yticks(np.linspace(0, np.amax(mass_dust_tot[t_arr<=tmax]), 5).round(2))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax.legend(loc="upper right", fontsize=fontsize-5)

plt.savefig(pngdir + f"{fnum}_dust_to_gas_{datestr}.png", dpi=300)