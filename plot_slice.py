from hconfig import *

#################################
date = "2023-07-25"
cat = True
dust = True
pressure = False
vlims = True
vlims_gas = (-26, -21) # g/cm^3
vlims_dust = (-30, -25.5) # g/cm^3
vlims_p = (2, 7) # P/k_b (K/cm^3)
vlims_T = (2, 8) # K
vlims_v = (-200, 1050)
# vlims_v = (-50, 150)
spacing = 640*1e-3 # spacing of tick marks in units
fontsize = 20
unit = "kpc" # sets axes labels and units of dx
fnum = None
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': 20})
#################################

# directory with slices
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5/slice/")
pngdir = os.path.join(basedir, "png/slice/")

data = ReadHDF5(datadir, nscalar=1, fnum=fnum, slice="xy", cat=cat)
head = data.head
conserved = data.conserved

nx = head["dims"][0]
ny = head["dims"][-1]
dx = head["dx"][0]
if unit == "pc":
    dx *= 1e3
d_gas = data.d_cgs()
d_dust = None
p_gas = None
if dust:
    d_dust = conserved["scalar0"] * head["density_unit"]
if pressure:
    p_gas = (data.energy_cgs() - 0.5*d_gas*((data.vx_cgs())**2 + (data.vy_cgs())**2 + (data.vz_cgs())**2)) * (head["gamma"] - 1.0) 
vx = data.vx_cgs()
t_arr = data.t_cgs() / yr_in_s
T = data.T()

for i, d in enumerate(d_gas):

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(25,9))
    
    # xy gas density slice
    if vlims:
        im = axs[0][0].imshow(np.log10(d_gas[i].T), origin="lower", vmin=vlims_gas[0], vmax=vlims_gas[1], extent=[0, nx*dx, 0, ny*dx])
    else:
        im = axs[0][0].imshow(np.log10(d_gas[i].T), origin="lower", extent=[0, nx*dx, 0, ny*dx])
    ylabel = r'$\mathrm{log}_{10}(\rho_{gas})$ [$\mathrm{g}\mathrm{cm}^{-3}$]'
    divider = make_axes_locatable(axs[0][0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[0][0], cax=cax)
    cbar.set_label(ylabel, fontsize=fontsize)
    axs[0][0].set_xticks(np.arange(0, nx*dx, spacing))
    axs[0][0].set_yticks(np.arange(0, ny*dx, spacing))
    axs[0][0].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=8)
    axs[0][0].set_title(r"Gas Density Slice", fontsize=fontsize)
    axs[0][0].set_xlabel(r"$x~$[{}]".format(unit), fontsize=fontsize)
    axs[0][0].set_ylabel(r"$y~$[{}]".format(unit), fontsize=fontsize)
    axs[0][0].text(spacing, 0.1*dx*ny, f'{round(t_arr[i]/1e6, 2)} Myr', color='white', fontsize=fontsize)
    

    # xy dust density
    if dust:
        wh_zero = np.where(d_dust[i]<=0)
        d_dust[i][wh_zero] = 1e-40

        if vlims:
            im = axs[1][0].imshow(np.log10(d_dust[i].T), origin="lower", cmap="plasma", vmin=vlims_dust[0], vmax=vlims_dust[1], extent=[0, nx*dx, 0, ny*dx])
        else:
            im = axs[1][0].imshow(np.log10(d_dust[i].T), origin="lower", cmap="plasma", vmin=vlims_dust[0], extent=[0, nx*dx, 0, ny*dx])
        ylabel = r'$\mathrm{log}_{10}(\rho_{dust})$ [$\mathrm{g}\mathrm{cm}^{-3}$]'
        divider = make_axes_locatable(axs[1][0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(im, ax=axs[1][0], cax=cax)
        cbar.set_label(ylabel, fontsize=fontsize)
        axs[1][0].set_xticks(np.arange(0, nx*dx, spacing))
        axs[1][0].set_yticks(np.arange(0, ny*dx, spacing))
        axs[1][0].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=8)
        axs[1][0].set_title(r"Dust Density Slice", fontsize=fontsize)
        axs[1][0].set_xlabel(r"$x~$[{}]".format(unit), fontsize=fontsize)
        axs[1][0].set_ylabel(r"$y~$[{}]".format(unit), fontsize=fontsize)

    # xy pressure
    if pressure:

        if vlims:
            im = axs[1][0].imshow(np.log10(p_gas[i].T/KB), origin="lower", cmap="magma", vmin=vlims_p[0], vmax=vlims_p[1], extent=[0, nx*dx, 0, ny*dx])
        else:
            im = axs[1][0].imshow(np.log10(p_gas[i].T/KB), origin="lower", cmap="magma", extent=[0, nx*dx, 0, ny*dx])
        ylabel = r'$\mathrm{log}_{10}(P_{gas}/k_B)$ [$K\,(cm^{-3})$]'
        divider = make_axes_locatable(axs[1][0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(im, ax=axs[1][0], cax=cax)
        cbar.set_label(ylabel, fontsize=fontsize)
        axs[1][0].set_xticks(np.arange(0, nx*dx, spacing))
        axs[1][0].set_yticks(np.arange(0, ny*dx, spacing))
        axs[1][0].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=8)
        axs[1][0].set_title(r"Gas Pressure Slice", fontsize=fontsize)
        axs[1][0].set_xlabel(r"$x~$[{}]".format(unit), fontsize=fontsize)
        axs[1][0].set_ylabel(r"$y~$[{}]".format(unit), fontsize=fontsize)
    

    # xy temperature slice
    if vlims:
        im = axs[0][1].imshow(np.log10(T[i].T), origin="lower", cmap="inferno", vmin=vlims_T[0], vmax=vlims_T[1], extent=[0, nx*dx, 0, ny*dx])
    else:
        im = axs[0][1].imshow(np.log10(T[i].T), origin="lower", cmap="inferno", extent=[0, nx*dx, 0, ny*dx])
    ylabel = r'$\mathrm{log}_{10}(T_{gas}) [\mathrm{K}]$'
    divider = make_axes_locatable(axs[0][1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[0][1], cax=cax)
    cbar.set_label(ylabel, fontsize=fontsize)
    axs[0][1].set_xticks(np.arange(0, nx*dx, spacing))
    axs[0][1].set_yticks(np.arange(0, ny*dx, spacing))
    axs[0][1].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=8)
    axs[0][1].set_title(r"Temperature Slice", fontsize=fontsize)
    axs[0][1].set_xlabel(r"$x~$[{}]".format(unit), fontsize=fontsize)
    axs[0][1].set_ylabel(r"$y~$[{}]".format(unit), fontsize=fontsize)

    # xy velocity slice
    if vlims:
        im = axs[1][1].imshow(vx[i].T*1e-5, origin="lower", vmin=vlims_v[0], vmax=vlims_v[1], extent=[0, nx*dx, 0, ny*dx])
    else:
        im = axs[1][1].imshow(vx[i].T*1e-5, origin="lower", extent=[0, nx*dx, 0, ny*dx])
    ylabel = r'$v_x$ [km/s]'
    divider = make_axes_locatable(axs[1][1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[1][1], cax=cax)
    cbar.set_label(ylabel, fontsize=fontsize)
    axs[1][1].set_xticks(np.arange(0, nx*dx, spacing))
    axs[1][1].set_yticks(np.arange(0, ny*dx, spacing))
    axs[1][1].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=8)
    axs[1][1].set_title(r"x-velocity slice", fontsize=fontsize)
    axs[1][1].set_xlabel(r"$x~$[{}]".format(unit), fontsize=fontsize)
    axs[1][1].set_ylabel(r"$y~$[{}]".format(unit), fontsize=fontsize)

    fig.tight_layout()
    
    # plot and save
    save = True
    if save:
        plt.savefig(pngdir + f"{i}_slice.png", dpi=300)
    plt.close()

    print(f"Saving figure {i} of {len(d_gas)}.\n")
