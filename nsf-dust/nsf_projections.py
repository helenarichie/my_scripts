from hconfig import *

# date = input("\nDate: ")
date = "2022-11-05"

# directory with slices
datadir = f"/ix/eschneider/helena/data/{date}/hdf5/"
outdir = f"/ix/eschneider/helena/data/{date}/png/nsf/"

fnum = 34

data = ReadHDF5(datadir, proj="xy", fnum=fnum)
head = data.head
conserved = data.conserved

######################################################
nx, ny, nz = head["dims"]
dx = head["dx"][0] * head["length_unit"]
d_gas = conserved["density"]
d_gas *= head["mass_unit"] / (head["length_unit"] ** 2)

gamma = head["gamma"]
t_arr = data.t_cgs() / yr_in_s
######################################################


datadir = f"/ix/eschneider/helena/data/{date}/hdf5/"
outdir = f"/ix/eschneider/helena/data/{date}/png/nsf/"
data = ReadHDF5(datadir, nscalar=1, fnum=fnum)
head = data.head
conserved = data.conserved

d_dust = conserved["scalar0"] * head["density_unit"]

wh_zero = np.where(d_dust<=0)
d_dust[wh_zero] = 1e-40
vlims_du = [np.log10(np.amin(d_dust.flatten())), np.log10(np.amax(d_dust.flatten()))]

d_dust = np.sum(d_dust*dx, axis=3)


for i, d in enumerate(d_gas):

    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(13,5))
    
    # xy gas density projection
    # im = axs[0].imshow(np.log10(d_gas[i].T), origin="lower", vmin=-28, vmax=-23, extent=[0, nx*dx, 0, nz*dx])
    #im = axs.imshow(np.log10(d_gas[i].T), origin="lower", vmin=-26.25, vmax=-24.25)
    im = axs.imshow(np.log10(d_gas[i].T), origin="lower")
    ylabel = r'$\mathrm{log}_{10}(\Sigma_{gas})$ [$\mathrm{g}\,\mathrm{cm}^{-2}$]'
    divider = make_axes_locatable(axs)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs, cax=cax)
    cbar.set_label(ylabel)
    axs.set_xticks(nx*np.arange(0, 1.25, 1/15))
    axs.set_yticks(nz*np.arange(0, 1.25, 0.25))
    axs.tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=7)
    axs.hlines(0.12*ny, 50, (50+51.2), color='white')
    axs.text(105, 0.1*ny, '150 pc', color='white')
    #axs[0].set_title(r"Gas Density Projection")
    #axs[0].plot([], [], " ", label="${:.1e}~yr$".format(t_arr[i]))
    #axs[0].legend()
    fig.tight_layout()

        # plot and save
    save = True
    if save:
        plt.savefig(outdir + f"{fnum}_gas_proj.png", dpi=300)
        plt.savefig(f"{fnum}_gas_proj.png", dpi=300)
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(13,5))
    # im = axs[1].imshow(np.log10(T[i].T), origin="lower", cmap="inferno", vmin=2, vmax=vlims_T[1], extent=[0, nx*dx, 0, nz*dx])
    #im = axs.imshow(np.log10(d_dust[i].T), origin="lower", cmap="plasma", vmin=-28.75, vmax=-26.75)
    im = axs.imshow(np.log10(d_dust[i].T), origin="lower", cmap="plasma", vmin=-5.25, vmax=-7.25)
    
    if np.isnan(np.log10(d_dust[i].T).any()):
        print("there's a nan")

    ylabel = r'$\mathrm{log}_{10}(\Sigma_{dust})$ [$\mathrm{g}\,\mathrm{cm}^{-2}$]'
    divider = make_axes_locatable(axs)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs, cax=cax)
    cbar.set_label(ylabel)
    axs.set_xticks(nx*np.arange(0, 1.25, 1/15))
    axs.set_yticks(nz*np.arange(0, 1.25, 0.25))
    axs.tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=7)
    axs.hlines(0.12*ny, 50, (50+51.2), color='white')
    axs.text(105, 0.1*ny, '150 pc', color='white')
    #axs[1].set_title(r"Dust Density Projection")
    #axs[1].set_xlabel(r"$x~$[kpc]")
    #axs[1].set_ylabel(r"$z~$[kpc]")
    fig.tight_layout()
    
    # plot and save
    save = True
    if save:
        plt.savefig(outdir + f"{fnum}_dust_proj.png", dpi=300)
        plt.savefig(f"{fnum}_dust_proj.png", dpi=300)
    plt.close()

    print(f"Saving figure {i+1} of {len(d_gas)}.\n")