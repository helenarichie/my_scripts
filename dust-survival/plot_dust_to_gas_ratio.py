import sys
sys.path.insert(0, "/ix/eschneider/helena/code/my_scripts")
from hconfig import *

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

##################################################################
date_0206 = "2024-02-06"
date_0207 = "2024-02-07"
datestr = "0206"
cat = True
fontsize = 20
labelpad = 12
tickwidth = 1
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': fontsize})
fnum_0207 = 3000
fnum_0206 = 300
rho_cl_init = 10*0.6*MP * density_conversion
xmin = 2.1
##################################################################

fnum = fnum_0207
fnum = fnum_0206

fnums = [fnum_0207, fnum_0206]

##################################################################
basedir_0207 = f"/ix/eschneider/helena/data/cloud_wind/{date_0207}/"
basedir_0206 = f"/ix/eschneider/helena/data/cloud_wind/{date_0206}/"
pngdir = os.path.join(basedir_0206, "png/")
##################################################################

fulldirs = [os.path.join(basedir_0207, "hdf5/full/"), os.path.join(basedir_0206, "hdf5/full/")]
labels = ["disrupted", "survived"]
colors = ["plum", "steelblue"]

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
    ax.set_ylabel(r"dust-to-gas ratio", fontsize=fontsize, labelpad=5)
    ax.set_xlim(xmin, np.amax(x_arr[wh_real]))
    #ax.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True)
    ax.set_xticks(np.linspace(xmin, np.amax(x_arr[wh_real]), 5).round(1))
    #ax.set_yticks(np.linspace(0, np.amax(mass_dust_tot[t_arr<=tmax]), 5).round(2))
    #ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax.legend(loc="upper right", fontsize=fontsize-5)

plt.savefig(pngdir + f"{fnum}_dust_to_gas_{datestr}.png", dpi=300, bbox_inches="tight")