from hconfig import *

date = input("\nDate: ")
nfile = input("\nFile number: ")

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5")

data = ReadHDF5(datadir, nscalar=1, slice="xy", fnum=nfile)
head = data.head
conserved = data.conserved

######################################################
nx = head["dims"][0]
nz = head["dims"][-1]
dx = head["dx"][0]
d_gas = data.d_cgs()
d_dust = conserved["scalar0"] * head["density_unit"]
vx = data.vx_cgs()
gamma = head["gamma"]
t_arr = data.t_cgs() / yr_in_s
T = data.T()
######################################################

wh_zero = np.where(d_dust<=0)
d_dust[wh_zero] = 1e-40
vlims_du = [np.log10(np.amin(d_dust.flatten())), np.log10(np.amax(d_dust.flatten()))]

wh_zero = np.where(T<=0)
T[wh_zero] = 1e-40
#vlims_T = [np.log10(np.amin(T.flatten())), np.log10(np.amax(T.flatten()))]
vlims_T = [2, 6]

for i, d in enumerate(d_gas):

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(25,9))
    
    # xy gas density projection
    #im = axs[0][0].imshow(np.log10(d_gas[i].T), origin="lower", vmin=-28, vmax=-23, extent=[0, nx*dx, 0, nz*dx])
    im = axs[0][0].imshow(np.log10(d_gas[i].T), origin="lower", extent=[0, nx*dx, 0, nz*dx])
    ylabel = r'$\mathrm{log}_{10}(\rho)$ [$\mathrm{g}\mathrm{cm}^{-3}$]'
    divider = make_axes_locatable(axs[0][0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[0][0], cax=cax)
    cbar.set_label(ylabel)
    axs[0][0].set_xticks(nx*dx*np.arange(0, 1.25, 0.25))
    axs[0][0].set_yticks(nz*dx*np.arange(0, 1.25, 0.5))
    axs[0][0].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1)
    axs[0][0].set_title(r"Gas Density Slice")
    axs[0][0].set_xlabel(r"$x~$[kpc]")
    axs[0][0].plot([], [], " ", label="${:.1e}~yr$".format(t_arr[i]))
    axs[0][0].legend()
    

    # xy density-weighted temperature projection
    im = axs[1][0].imshow(np.log10(d_dust[i].T), origin="lower", cmap="plasma", vmin=-28, vmax=-24, extent=[0, nx*dx, 0, nz*dx])
    ylabel = r'$\mathrm{log}_{10}(\rho)$ [$\mathrm{g}\mathrm{cm}^{-3}$]'
    divider = make_axes_locatable(axs[1][0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[1][0], cax=cax)
    cbar.set_label(ylabel)
    axs[1][0].set_xticks(nx*dx*np.arange(0, 1.25, 0.25))
    axs[1][0].set_yticks(nz*dx*np.arange(0, 1.25, 0.5))
    axs[1][0].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1)
    axs[1][0].set_title(r"Dust Density Slice")
    axs[1][0].set_xlabel(r"$x~$[kpc]")
    axs[1][0].set_ylabel(r"$z~$[kpc]")
    axs[1][0].plot([], [], " ", label="${:.1e}~yr$".format(t_arr[i]))
    axs[1][0].legend()
    

    #im = axs[0][1].imshow(np.log10(T[i].T), origin="lower", cmap="inferno", vmin=0, vmax=vlims_T[1], extent=[0, nx*dx, 0, nz*dx])
    im = axs[0][1].imshow(np.log10(T[i].T), origin="lower", cmap="inferno", extent=[0, nx*dx, 0, nz*dx])
    
    if np.isnan(np.log10(T[i].T).any()):
        print("there's a nan")

    ylabel = r'$\mathrm{log}_{10}(T) [\mathrm{K}]$'
    divider = make_axes_locatable(axs[0][1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[0][1], cax=cax)
    cbar.set_label(ylabel)
    axs[0][1].set_xticks(nx*dx*np.arange(0, 1.25, 0.25))
    axs[0][1].set_yticks(nz*dx*np.arange(0, 1.25, 0.5))
    axs[0][1].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1)
    axs[0][1].set_title(r"Temperature Slice")
    axs[0][1].set_xlabel(r"$x~$[kpc]")
    axs[0][1].set_ylabel(r"$z~$[kpc]")
    axs[0][1].plot([], [], " ", label="${:.1e}~yr$".format(t_arr[i]))
    axs[0][1].legend()

    # xz velocity slice
    # im = axs[1][1].imshow(vx[i].T, origin="lower", vmin=-1e-3, vmax=0.011, extent=[0, nx*dx, 0, nz*dx])
    # im = axs[1][1].imshow(vx[i].T*1e-5, origin="lower", extent=[0, nx*dx, 0, nz*dx], vmin=0, vmax=200)
    im = axs[1][1].imshow(vx[i].T*1e-5, origin="lower", extent=[0, nx*dx, 0, nz*dx])
    ylabel = r'x-velocity [km/s]'

    divider = make_axes_locatable(axs[1][1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[1][1], cax=cax)
    cbar.set_label(ylabel)
    axs[1][1].set_xticks(nx*dx*np.arange(0, 1.25, 0.25))
    axs[1][1].set_yticks(nz*dx*np.arange(0, 1.25, 0.5))
    axs[1][1].tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1)
    axs[1][1].set_title(r"x-velocity")
    axs[1][1].set_xlabel(r"$x~$[kpc]")
    axs[1][1].set_ylabel(r"$z~$[kpc]")
    axs[1][1].plot([], [], " ", label="${:.1e}~yr$".format(t_arr[i]))
    axs[1][1].legend()

    fig.tight_layout()
    
    # plot and save
    save = True
    if save:
        plt.savefig(basedir + "snapshot.png", dpi=300)
    plt.close()

    print(f"Saving figure {i+1} of {len(d_dust)}.\n")