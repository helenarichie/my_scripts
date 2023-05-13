from hconfig import *

##################################################################
date = "2023-03-22"
save = True
cat = True
dust = False
vlims_gas = (-6.0005, -2.5)
vlims_dust = (-10.5, -4.5)
##################################################################

##################################################################
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
projdir = os.path.join(basedir, "hdf5/proj/")
pngdir = os.path.join(basedir, "png/flatiron/")

data = ReadHDF5(projdir, dust=dust, proj="xy", cat=cat)
head = data.head
conserved = data.conserved

nx, ny, nz = head["dims"]
dx = head["dx"][0] * 1e3  # pc
d_gas = conserved["density"]
d_gas *= head["mass_unit"] / (head["length_unit"] ** 2)

gamma = head["gamma"]
t_arr = data.t_cgs() / yr_in_s  # yr

d_dust = conserved["dust_density"]

d_dust *= head["mass_unit"] / (head["length_unit"] ** 2)
##################################################################

for i, d in enumerate(d_gas):

    plt.style.use('dark_background')

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(20, 10))

    spacing = 20  # pc, spacing of tick marks
    fontsize = 28
    pad = 0.1
    labelpad = 12
    tickwidth = 2
    plt.rcParams.update({'font.family': 'Helvetica'})
    plt.rcParams.update({'font.size': 25})

    # xy gas density projection
    im = axs[0].imshow(np.log10(d_gas[i].T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=[0, nx*dx, 0, nz*dx])
    ylabel = r'$\mathrm{log}_{10}(\Sigma_{gas})$ [$\mathrm{g}\,\mathrm{cm}^{-2}$]'
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="5%", pad=pad)
    cbar = fig.colorbar(im, ax=axs[0], cax=cax, pad=pad)
    cbar.set_label(ylabel, fontsize=fontsize, labelpad=11)
    axs[0].hlines(0.13*dx*ny, spacing, spacing+spacing, color='white')
    axs[0].text(spacing+spacing+2, 0.1*dx*ny, '20 pc', color='white', fontsize=fontsize)
    axs[0].set_xticks(np.arange(0, nx*dx, spacing))
    axs[0].set_yticks(np.arange(0, ny*dx, spacing))
    axs[0].tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=9, width=tickwidth)
    #axs[0].text(spacing, 0.1*dx*nz, f'{round(t_arr[i]/1e6, 2)} Myr', color='white', fontsize=fontsize)    
    #axs[0].set_title(r"Gas Density Projection")

    d_dust[i][d_dust[i]<=0] = 1e-40

    # xy dust density projection
    im = axs[1].imshow(np.log10(d_dust[i].T), origin="lower", cmap="plasma", vmin=vlims_dust[0], vmax=vlims_dust[1], extent=[0, nx*dx, 0, nz*dx])
    if np.isnan(np.log10(d_dust[i].T).any()):
        print("there's a nan")
    ylabel = r'$\mathrm{log}_{10}(\Sigma_{dust})$ [$\mathrm{g}\,\mathrm{cm}^{-2}$]'
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("right", size="5%", pad=pad)
    cbar = fig.colorbar(im, ax=axs[1], cax=cax)
    cbar.set_label(ylabel, fontsize=fontsize, labelpad=2)
    axs[1].set_xticks(np.arange(0, nx*dx, spacing))
    axs[1].set_yticks(np.arange(0, ny*dx, spacing))
    axs[1].tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=9, width=tickwidth)
    axs[1].text(spacing, 0.1*dx*nz, f'{round(t_arr[i]/1e6, 1)} Myr', color='white', fontsize=fontsize)  
    #axs[1].set_title(r"Dust Density Projection")

    fig.tight_layout()

    # plot and save
    save = True
    if save:
        plt.savefig(pngdir + f"{i}_proj.png", dpi=300)
    plt.close()

    print(f"Saving figures {i} of {len(os.listdir(projdir))-1}.\n")
