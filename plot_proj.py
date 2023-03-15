from hconfig import *

# date = input("\nDate: ")
#################################
date = "2023-03-14"
cat = False
vlims_gas = (-4.75, -3.0)
vlims_dust = (-8, -7.25)
save = True
gas = True
dust = True
#################################

# directory with slices
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
proj_data = os.path.join(basedir, "hdf5/")
full_data = os.path.join(basedir, "hdf5/")
pngdir = os.path.join(basedir, "png/proj/")


if gas:

    data = ReadHDF5(proj_data, nscalar=1, proj="xy", cat=cat)
    head = data.head
    conserved = data.conserved

    nx, ny, nz = head["dims"]
    dx = head["dx"][0]
    d_gas = conserved["density"]
    d_gas *= head["mass_unit"] / (head["length_unit"] ** 2)

    gamma = head["gamma"]
    t_arr = data.t_cgs() / yr_in_s

    for i, d in enumerate(d_gas):

        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(15,5))
        
        # xz gas density projection
        #im = axs.imshow(np.log10(d_gas[i].T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=[0, nx*dx, 0, nz*dx])
        im = axs.imshow(np.log10(d_gas[i].T), origin="lower", extent=[0, nx*dx, 0, nz*dx])
        ylabel = r'$\mathrm{log}_{10}(\Sigma_{gas})$ [$\mathrm{g}\,\mathrm{cm}^{-2}$]'
        divider = make_axes_locatable(axs)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(im, ax=axs, cax=cax)
        cbar.set_label(ylabel)
        axs.set_xticks(nx*dx*np.arange(0, 1.25, 0.25))
        axs.set_yticks(nz*dx*np.arange(0, 1.25, 0.5))
        axs.set_xlabel(r"$x~$[kpc]")
        axs.set_ylabel(r"$z~$[kpc]")
        axs.tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1, length=7)
        axs.text(0.05*dx*nx, 0.1*dx*nz, f'{round(t_arr[i]/1e6, 2)} Myr', color='white', fontsize=20)    
        axs.set_title(r"Gas Density Projection")

        # plot and save
        save = True
        if save:
            plt.savefig(pngdir + f"{i}_gas_proj.png", dpi=300)
        plt.close()

        print(f"Saving figure {i+1} of {len(d_gas)}.\n")

    del d_gas

if dust:
    for i in range(0, len(os.listdir(full_data))):
        data = ReadHDF5(full_data, fnum=i, nscalar=1, cat=cat)
        head = data.head
        conserved = data.conserved

        dx_cgs = head["dx"][0] * head["length_unit"]
        nx, ny, nz = head["dims"]
        dx = head["dx"][0]

        t = data.t_cgs() / yr_in_s

        d_dust = conserved["scalar0"] * head["density_unit"]

        wh_zero = np.where(d_dust<=0)
        d_dust[wh_zero] = 1e-40
        vlims_du = [np.log10(np.amin(d_dust.flatten())), np.log10(np.amax(d_dust.flatten()))]

        d_dust = np.sum(d_dust*dx_cgs, axis=3)

        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(13,5))

        # xz dust density projection
        #im = axs.imshow(np.log10(d_gas[i].T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=[0, nx*dx, 0, nz*dx])
        im = axs.imshow(np.log10(d_dust[0].T), origin="lower", cmap="plasma", vmin=vlims_dust[0], extent=[0, nx*dx, 0, nz*dx])

        if np.isnan(np.log10(d_dust.T).any()):
            print("there's a nan")

        ylabel = r'$\mathrm{log}_{10}(\Sigma_{gas})$ [$\mathrm{g}\,\mathrm{cm}^{-2}$]'
        divider = make_axes_locatable(axs)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(im, ax=axs, cax=cax)
        cbar.set_label(ylabel)
        axs.set_xticks(nx*dx*np.arange(0, 1.25, 0.25))
        axs.set_yticks(nz*dx*np.arange(0, 1.25, 0.5))
        axs.set_xlabel(r"$x~$[kpc]")
        axs.set_ylabel(r"$z~$[kpc]")
        axs.tick_params(axis='both', which='both', direction='in', color='white', top=1, right=1, length=7)
        axs.text(0.05*dx*nx, 0.1*dx*nz, f'{round(t[0]/1e6, 2)} Myr', color='white', fontsize=20)  
        axs.set_title(r"Dust Density Projection")

        # plot and save
        save = True
        if save:
            plt.savefig(pngdir + f"{i}_dust_proj.png", dpi=300)
        plt.close()

        print(f"Saving figure {i} of {len(os.listdir(full_data))}.\n")