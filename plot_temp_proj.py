from hconfig import *

# date = input("\nDate: ")
#################################
date = "2023-03-22"
cat = True
save = True
#################################

# directory with slices
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
proj_data = os.path.join(basedir, "hdf5/proj")
pngdir = os.path.join(basedir, "png/proj/")

data = ReadHDF5(proj_data, proj="xy", cat=cat)
head = data.head
conserved = data.conserved

nx, ny, nz = head["dims"]
dx = head["dx"][0]
T = conserved["temperature"]

gamma = head["gamma"]
t_arr = data.t_cgs() / yr_in_s

for i, d in enumerate(T):

    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(15,5))
    
    # xz gas density projection
    #im = axs.imshow(np.log10(d_gas[i].T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=[0, nx*dx, 0, nz*dx])
    im = axs.imshow(np.log10(T[i].T), origin="lower", extent=[0, nx*dx, 0, nz*dx])
    ylabel = r'$\mathrm{log}_{10}(T)$ [$K$]'
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
    axs.set_title(r"Gas Temperature Projection")

    # plot and save
    save = True
    if save:
        plt.savefig(pngdir + f"{i}_temp_proj.png", dpi=300)
    plt.close()

    print(f"Saving figure {i+1} of {len(T)}.\n")

del T