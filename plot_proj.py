from hconfig import *

# date = input("\nDate: ")
date = "2022-11-16"

# directory with slices
datadir = f"/ix/eschneider/helena/data/cloud_wind/{date}/chi-1e3/hdf5/proj/"
outdir = f"/ix/eschneider/helena/data/cloud_wind/{date}/chi-1e3/png/proj/"

data = ReadHDF5(datadir, proj="xy")
head = data.head
conserved = data.conserved

######################################################
nx = head["dims"][0]
nz = head["dims"][-1]
dx = head["dx"][0]
dx_cgs = head["dx"][0] * head["length_unit"]
d_gas = conserved["density"]
gamma = head["gamma"]
t_arr = data.t_cgs() / yr_in_s
T = conserved["temperature"]
######################################################

d_gas *= head["mass_unit"] / (head["length_unit"] ** 2) # surface density
# d_dust = np.sum(d_dust*dx_cgs, axis=3) # calculate surface density from volume density

wh_zero = np.where(T<=0)
T[wh_zero] = 1e-40
# vlims_T = [np.log10(np.amin(T.flatten())), np.log10(np.amax(T.flatten()))]
vlims_T = [1, 7]

for i, d in enumerate(d_gas):

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(15,10))
    
    # xy gas density projection
    # im = axs[0].imshow(np.log10(d_gas[i].T), origin="lower", vmin=-28, vmax=-23, extent=[0, nx*dx, 0, nz*dx])
    im = axs[0].imshow(np.log10(d_gas[i].T), origin="lower", extent=[0, nx*dx, 0, nz*dx])
    ylabel = r'$\mathrm{log}_{10}(\Sigma)$ [$\mathrm{g}\mathrm{cm}^{-2}$]'
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[0], cax=cax)
    cbar.set_label(ylabel)
    axs[0].set_xticks(nx*dx*np.arange(0, 1.25, 0.25))
    axs[0].set_yticks(nz*dx*np.arange(0, 1.25, 0.5))
    axs[0].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1)
    axs[0].set_title(r"Gas Density Projection")
    axs[0].set_xlabel(r"$x~$[kpc]")
    axs[0].plot([], [], " ", label="${:.1e}~yr$".format(t_arr[i]))
    axs[0].legend()
    

    im = axs[1].imshow(np.log10(T[i].T), origin="lower", cmap="inferno", extent=[0, nx*dx, 0, nz*dx])
    
    if np.isnan(np.log10(T[i].T).any()):
        print("there's a nan")

    ylabel = r'$\mathrm{log}_{10}(T) [\mathrm{K}]$'
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[1], cax=cax)
    cbar.set_label(ylabel)
    axs[1].set_xticks(nx*dx*np.arange(0, 1.25, 0.25))
    axs[1].set_yticks(nz*dx*np.arange(0, 1.25, 0.5))
    axs[1].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1)
    axs[1].set_title(r"Density-Weighted Temperature Projection")
    axs[1].set_xlabel(r"$x~$[kpc]")
    axs[1].set_ylabel(r"$z~$[kpc]")
    axs[1].plot([], [], " ", label="${:.1e}~yr$".format(t_arr[i]))
    axs[1].legend()

    fig.tight_layout()
    
    # plot and save
    save = True
    if save:
        plt.savefig(outdir + f"{i}_proj.png", dpi=300)
    plt.close()

    print(f"Saving figure {i+1} of {len(d_gas)}.\n")