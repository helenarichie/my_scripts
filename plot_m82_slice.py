from hconfig import *

#################################
date = "2024-04-05"
ns = 0
ne = 30
cat = True
#################################

########### location ############
crc = True
frontier = False
#################################

########## data type ############
debugging = False
m82 = True
testing = False
#################################

########### plotting #############
dust = True
basic_scalar = True
vlims = True
vlims_gas = (32.7 , 38.5) # g/cm^3
vlims_dust = (30, 36.7) # g/cm^3
vlims_p = (3, 6) # P/k_b (K/cm^3)
vlims_T = (5, 7) # K
vlims_v = (0, 150)
spacing, unit = 1000*1e-3, "kpc" # spacing of tick marks in units and sets axes labels and units of dx (kpc or pc)
# spacing, unit = 40, "pc"
fontsize = 20
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': 20})
#################################

########## specify slice and png directories ############
if crc:
    if debugging:
        basedir = f"/ix/eschneider/helena/data/debugging/{date}/"
    if m82:
        basedir = f"/ix/eschneider/helena/data/m82/{date}/"
    if testing:
        basedir = f"/ix/eschneider/helena/data/testing/{date}/"

if frontier:
    basedir = f"/lustre/orion/ast181/scratch/helenarichie/{date}/"

if cat:
    datadir = os.path.join(basedir, "hdf5/slice/")
else:
    datadir = os.path.join(basedir, "hdf5/raw/")

pngdir = os.path.join(basedir, "png/slice/")
#########################################################

for i in range(ns, ne+1):

    # read in data
    f = h5py.File(os.path.join(datadir, str(i)+"_slice.h5"))
    head = f.attrs
    dx = head["dx"][0]
    t = head["t"][0]*1e3
    gamma = head["gamma"][0]
    mu = 0.6

    nx, ny, nz = head["dims"][0], head["dims"][1], head["dims"][2]
    dx = head["dx"][0]

    d_gas = np.array(f["d_xz"]) * head["density_unit"] * (dx*head["length_unit"])**3
    if dust:
        d_dust = np.array(f["d_dust_xz"]) * head["density_unit"] * (dx*head["length_unit"])**3
        print("Dust min: ", np.amin(d_dust))
        print("Dust max: ", np.amax(d_dust))
    if basic_scalar:
        scalar = np.array(f["basic_scalar_xz"])
        print("Scalar min: ", np.amin(scalar))
        print("Scalar max: ", np.amax(scalar))

    # p_gas = (data.energy_cgs() - 0.5*d_gas*((vx)**2 + (data.vy_cgs())**2 + (data.vz_cgs())**2)) * (head["gamma"] - 1.0)

    # plot
    if basic_scalar:
        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(23, 13))
    else:
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15,13))
    if vlims:
        im = ax[0].imshow(np.log10(d_gas.T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=[0, nx*dx, 0, nz*dx])
    else:
        im = ax[0].imshow(np.log10(d_gas.T), origin="lower", extent=[0, nx*dx, 0, nz*dx])
    
    ylabel = r'$\mathrm{log}_{10}(\rho_{gas})$ [$\mathrm{g}\mathrm{cm}^{-3}$]'
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=ax[0], cax=cax)
    cbar.set_label(ylabel, fontsize=fontsize)
    ax[0].set_xticks(np.arange(0, nx*dx, spacing))
    ax[0].set_yticks(np.arange(0, nz*dx, spacing))
    ax[0].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=8)
    ax[0].set_title(r"Gas Density Slice", fontsize=fontsize)
    ax[0].set_xlabel(r"$x~$[{}]".format(unit), fontsize=fontsize)
    ax[0].set_ylabel(r"$z~$[{}]".format(unit), fontsize=fontsize)
    ax[0].text(spacing, 0.1*dx*nz, f'{round(t/1e6, 2)} Myr', color='white', fontsize=fontsize)
    
    # xy dust density
    if dust:
        wh_zero = np.where(d_dust==0)
        d_dust[wh_zero] = 1e-40
        wh_neg = np.where(d_dust<0)
        d_dust[wh_neg] = np.nan

        if vlims:
            im = ax[1].imshow(np.log10(d_dust.T), origin="lower", cmap="plasma", vmin=vlims_dust[0], vmax=vlims_dust[1], extent=[0, nx*dx, 0, nz*dx])
        else:
            im = ax[1].imshow(np.log10(d_dust.T), origin="lower", cmap="plasma", extent=[0, nx*dx, 0, nz*dx], vmin=30)
        ylabel = r'$\mathrm{log}_{10}(\rho_{dust})$ [$\mathrm{g}\mathrm{cm}^{-3}$]'
        divider = make_axes_locatable(ax[1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(im, ax=ax[1], cax=cax)
        cbar.set_label(ylabel, fontsize=fontsize)
        ax[1].set_xticks(np.arange(0, nx*dx, spacing))
        ax[1].set_yticks(np.arange(0, nz*dx, spacing))
        ax[1].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=8)
        ax[1].set_title(r"Dust Density Slice", fontsize=fontsize)
        ax[1].set_xlabel(r"$x~$[{}]".format(unit), fontsize=fontsize)
        ax[1].set_ylabel(r"$z~$[{}]".format(unit), fontsize=fontsize)
    
    if basic_scalar:
        if vlims:
            im = ax[2].imshow(np.log10(scalar.T), origin="lower", cmap="jet", extent=[0, nx*dx, 0, nz*dx], vmin=2, vmax=10)
        else: 
            im = ax[2].imshow(np.log10(scalar.T), origin="lower", cmap="jet", extent=[0, nx*dx, 0, nz*dx])
        ylabel = "scalar"
        divider = make_axes_locatable(ax[2])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(im, ax=ax[2], cax=cax)
        cbar.set_label(ylabel, fontsize=fontsize)
        ax[2].set_xticks(np.arange(0, nx*dx, spacing))
        ax[2].set_yticks(np.arange(0, nz*dx, spacing))
        ax[2].tick_params(axis="both", which="both", direction="in", color="black", top=1, right=1, length=8)
        ax[2].set_title("scalar", fontsize=fontsize)
        ax[2].set_xlabel(r"$x~$[{}]".format(unit), fontsize=fontsize)
        ax[2].set_ylabel(r"$z~$[{}]".format(unit), fontsize=fontsize)

    fig.tight_layout()
    
    # plot and save
    plt.savefig(pngdir + f"{i}_m82_slice.png", dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saving figure {i} of {ne}.\n")
