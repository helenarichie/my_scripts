from hconfig import *

# date = input("\nDate: ")
#################################
date = "2023-03-22"
cat = True
vlims_gas = (-4.75, -3.0)
vlims_dust = (-8.75, -7.25)
save = True
gas = True
dust = True
#################################

# directory with slices
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
proj_data = os.path.join(basedir, "hdf5/proj")
full_data = os.path.join(basedir, "hdf5/full")
pngdir = os.path.join(basedir, "png/x/")

for i in range(0, len(os.listdir(full_data))):

    data = ReadHDF5(full_data, fnum=i, nscalar=1, cat=cat)
    head = data.head
    conserved = data.conserved

    dx = head["dx"][0]
    nx, ny, nz = head["dims"]

    t = data.t_cgs() / yr_in_s

    d_dust = conserved["scalar0"][0]

    wh_zero = np.where(d_dust<1e-10)
    d_dust[wh_zero] = 1e-10
    vlims_du = [np.log10(np.amin(d_dust.flatten())), np.log10(np.amax(d_dust.flatten()))]

    fig, axs = plt.subplots(nrows=1, ncols=1)

    # xz dust density projection
    #im = axs.imshow(np.log10(d_gas[i].T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=[0, nx*dx, 0, nz*dx])
    im = axs.imshow(np.log10(d_dust[-1, :, :].T), origin="lower", cmap="plasma")

    if np.isnan(np.log10(d_dust.T).any()):
        print("there's a nan")

    ylabel = r'$\mathrm{log}_{10}(\rho_{gas})$ [$\mathrm{M_\odot}\,\mathrm{kpc}^{-3}$]'
    divider = make_axes_locatable(axs)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs, cax=cax)
    cbar.set_label(ylabel, fontsize=20)
    cbar.ax.tick_params(labelsize=20)
    axs.xaxis.set_tick_params(labelbottom=False)
    axs.yaxis.set_tick_params(labelleft=False)
    #axs.set_xticks(nx*dx*np.arange(0, 1.25, 0.25))
    #axs.set_yticks(nz*dx*np.arange(0, 1.25, 0.5))
    #axs.set_xlabel(r"$x~$[kpc]")
    #axs.set_ylabel(r"$z~$[kpc]")
    #axs.tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1, length=7)
    axs.text(2*dx*nz, 2*dx*nz, f'{round(t[0]/1e6, 2)} Myr', color='white', fontsize=20)  
    axs.set_title(r"+x boundary slice", fontsize=20)

    # plot and save
    save = True
    if save:
        plt.tight_layout()
        plt.savefig(pngdir + f"{i}.png", dpi=300)
    plt.close()