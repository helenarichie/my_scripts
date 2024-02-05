import sys
sys.path.insert(0, "/ix/eschneider/helena/code/my_scripts")
from hconfig import *
from csv import writer

#################################
date = "2024-01-29"
ns = 0
ne = 5
n_procs = 2
#################################

########### plotting ############
dust = True
cat = True
vlims = False
unit = "kpc"  # kpc or pc
spacing = 640*1e-3  # in kpc or pc
# spacing = 20
fontsize = 28
vlims_gas = (19.5, 25.5)  # g/cm^3
vlims_dust = (-14, -3)  # g/cm^3
pad = 0.1
labelpad = 12
tickwidth = 2
ticklen = 9
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': fontsize})
#################################

########### location ############
crc = True
frontier = False
#################################

########## data type ############
debugging = False
cloud_wind = False
testing = True
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
    projdir = os.path.join(basedir, "hdf5/proj/")
else:
    projdir = os.path.join(basedir, "hdf5/raw/")

csvdir = os.path.join(basedir, "csv/")
pngdir = os.path.join(basedir, "png/proj/")
#########################################################

data = ReadHDF5(projdir, dust=dust, proj="xy", cat=cat)
head = data.head
conserved = data.conserved

nx, ny, nz = head["dims"]
dx = head["dx"][0]
if unit == "pc":
    dx *= 1e3

d_gas = conserved["density"][0]
d_gas *= head["mass_unit"] / (head["length_unit"] ** 2)

gamma = head["gamma"]
t_arr = data.t_cgs() / yr_in_s  # yr

d_dust = conserved["dust_density"][0]

d_dust *= head["mass_unit"] / (head["length_unit"] ** 2)
##################################################################
# create files to write quantities for cloud (dense gas only), total gas, and dust
#f = open(os.path.join(csvdir, "d_dust_neg.csv"), "w")
#f.close()

for i in range(ns, ne):

    # read in data
    if dust:
        data = ReadHDF5(projdir, fnum=i, dust=dust, proj="xy", cat=cat)
    else:
        data = ReadHDF5(projdir, fnum=i, proj="xy", cat=cat)

    """with open(os.path.join(csvdir, "d_dust_neg.csv"), "a") as f:
        d_dust_ch = d_dust[i] / (head["mass_unit"] / (head["length_unit"] ** 2))
        writer_obj = writer(f)
        if len(d_dust_ch[d_dust_ch<0]) > 0:
            writer_obj.writerow([np.amin(d_dust_ch[d_dust_ch<0]), np.amax(d_dust_ch[d_dust_ch<0]), len(d_dust_ch[d_dust_ch<0])])
        else:
            writer_obj.writerow(["None"])
        f.close()"""


    # plt.style.use('dark_background')

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(17,10))

    # xy gas density projection
    n_gas = d_gas/(0.6*MP) # column density
    if vlims:
        im = axs[0].imshow(np.log10(n_gas.T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=[0, nx*dx, 0, nz*dx])
    else:
        im = axs[0].imshow(np.log10(n_gas.T), origin="lower", extent=[0, nx*dx, 0, nz*dx])
    ylabel = r'$\mathrm{log}_{10}(N_{H, gas})$ [$\mathrm{cm}^{-2}$]'
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="5%", pad=pad)
    cbar = fig.colorbar(im, ax=axs[0], cax=cax, pad=pad)
    cbar.ax.tick_params(length=9, width=tickwidth)
    cbar.set_label(ylabel, fontsize=fontsize, labelpad=11)
    cbar.set_ticks(np.linspace(vlims_gas[0], vlims_gas[1], 4).round(1))
    axs[0].set_xticks(np.arange(0, nx*dx, spacing))
    axs[0].set_yticks(np.arange(0, ny*dx, spacing))
    axs[0].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1, length=ticklen, width=tickwidth)
    # axs[0].set_title(r"Gas Density Projection", fontsize=fontsize)
    axs[0].set_xlabel(r"$x~$[{}]".format(unit), fontsize=fontsize)
    axs[0].set_ylabel(r"$y~$[{}]".format(unit), fontsize=fontsize)
    axs[0].text(spacing, 0.1*dx*ny, f'{round(t_arr[i]/1e6, 2)} Myr', color='white', fontsize=fontsize)

    negative = d_dust[i][d_dust[i]<0]
    # d_dust[i][d_dust<0] = 1e-40
    # d_dust[i][d_dust==0] = 1e-40

    if len(negative) > 0:
        print(len(negative))
        print("Max: ", np.amax(negative))
        print("Min: ", np.amin(negative))
    if np.any(np.isnan(d_dust[i])):
        print(d_dust[i][np.isnan(d_dust[i])])

    # xy dust density projection
    im = axs[1].imshow(np.log10(d_dust.T), origin="lower", cmap="plasma", vmin=vlims_dust[0], vmax=vlims_dust[1], extent=[0, nx*dx, 0, nz*dx])
    ylabel = r'$\mathrm{log}_{10}(\Sigma_{dust})$ [$\mathrm{g}\,\mathrm{cm}^{-2}$]'
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("right", size="5%", pad=pad)
    cbar = fig.colorbar(im, ax=axs[1], cax=cax)
    cbar.ax.tick_params(length=9, width=tickwidth)
    cbar.set_label(ylabel, fontsize=fontsize, labelpad=11)
    cbar.set_ticks(np.linspace(vlims_dust[0], vlims_dust[1], 4).round(1))
    axs[1].set_xticks(np.arange(0, nx*dx, spacing))
    axs[1].set_yticks(np.arange(0, ny*dx, spacing))
    axs[1].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1, length=ticklen, width=tickwidth)
    # axs[1].set_title(r"Dust Density Slice", fontsize=fontsize)
    axs[1].set_xlabel(r"$x~$[{}]".format(unit), fontsize=fontsize)
    axs[1].set_ylabel(r"$y~$[{}]".format(unit), fontsize=fontsize)

    fig.tight_layout()

    # plot and save
    plt.savefig(pngdir + f"{i}_proj.png", dpi=300)
    plt.close()

    print(f"Saving figures {i} of {len(os.listdir(projdir))-1}.\n")