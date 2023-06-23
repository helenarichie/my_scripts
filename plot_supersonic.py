from hconfig import *

#################################
date = "2023-05-13"
cat = True
dust = False
pressure = True
vlims = True
vlims_vx = (0, 29)
vlims_vy = (0, 2)
vlims_vz = (0, 2)
spacing = 40 # spacing of tick marks, pc
fontsize = 20
fnum = None
unit = "pc"
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': 20})
#################################

# directory with slices
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5/slice/")
pngdir = os.path.join(basedir, "png/mach/")

data = ReadHDF5(datadir, fnum=fnum, nscalar=1, slice="xz", cat=cat)
head = data.head
conserved = data.conserved

nx = head["dims"][0]
ny = head["dims"][-1]
dx = head["dx"][0]
gamma = head["gamma"]
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
vy = data.vy_cgs()
vz = data.vz_cgs()
t_arr = data.t_cgs() / yr_in_s
T = data.T()

for i, d in enumerate(d_gas):

    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10,10))

    cs = np.sqrt(gamma*p_gas/d_gas)

    if vlims:
        im = axs[0].imshow((np.abs(vx[i])/cs[i]).T, origin="lower", vmin=vlims_vx[0], vmax=vlims_vx[1], extent=[0, nx*dx, 0, ny*dx])
    else:
        im = axs[0].imshow((np.abs(vx[i])/cs[i]).T, origin="lower", extent=[0, nx*dx, 0, ny*dx])
    ylabel = r'$\mathcal{M}_x$'
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[0], cax=cax)
    cbar.set_label(ylabel, fontsize=fontsize)
    # cbar.set_ticks(np.linspace(np.nanmin(np.abs(vx[i])/cs[i]), np.nanmax(np.abs(vx[i])/cs[i]), 5).round(1))
    cbar.set_ticks(np.linspace(vlims_vx[0], vlims_vx[1], 5).round(1))
    axs[0].set_xticks(np.arange(0, nx*dx, spacing))
    axs[0].set_yticks(np.arange(0, ny*dx, spacing))
    axs[0].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=8)
    #axs[0].set_title(r"x-velocity slice", fontsize=fontsize)
    axs[0].set_xlabel(r"$x~$[{}]".format(unit), fontsize=fontsize)
    axs[0].set_ylabel(r"$z~$[{}]".format(unit), fontsize=fontsize)

    cs = np.sqrt(gamma*p_gas/d_gas)

    if vlims:
        im = axs[1].imshow((np.abs(vy[i])/cs[i]).T, origin="lower", vmin=vlims_vy[0], vmax=vlims_vy[1], extent=[0, nx*dx, 0, ny*dx])
    else:
        im = axs[1].imshow((np.abs(vy[i])/cs[i]).T, origin="lower", extent=[0, nx*dx, 0, ny*dx])
    ylabel = r'$\mathcal{M}_y$'
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[1], cax=cax)
    cbar.set_label(ylabel, fontsize=fontsize)
    cbar.set_ticks(np.linspace(vlims_vy[0], vlims_vy[1], 5).round(1))
    #cbar.set_ticks(np.linspace(np.nanmin(np.abs(vy[i])/cs[i]), np.nanmax(np.abs(vy[i])/cs[i]), 5).round(1))
    axs[1].set_xticks(np.arange(0, nx*dx, spacing))
    axs[1].set_yticks(np.arange(0, ny*dx, spacing))
    axs[1].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=8)
    #axs[1].set_title(r"", fontsize=fontsize)
    axs[1].set_xlabel(r"$x~$[{}]".format(unit), fontsize=fontsize)
    axs[1].set_ylabel(r"$z~$[{}]".format(unit), fontsize=fontsize)

    # xy velocity slice

    cs = np.sqrt(gamma*p_gas/d_gas)

    if vlims:
        im = axs[2].imshow((np.abs(vz[i])/cs[i]).T, origin="lower", vmin=vlims_vz[0], vmax=vlims_vz[1], extent=[0, nx*dx, 0, ny*dx])
    else:
        im = axs[2].imshow((np.abs(vz[i])/cs[i]).T, origin="lower", extent=[0, nx*dx, 0, ny*dx])
    ylabel = r'$\mathcal{M}_z$'
    divider = make_axes_locatable(axs[2])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, ax=axs[2], cax=cax)
    cbar.set_label(ylabel, fontsize=fontsize)
    cbar.set_ticks(np.linspace(vlims_vz[0], vlims_vz[1], 5).round(1))
    #cbar.set_ticks(np.linspace(np.nanmin(np.abs(vz[i])/cs[i]), np.nanmax(np.abs(vz[i])/cs[i]), 5).round(1))
    axs[2].set_xticks(np.arange(0, nx*dx, spacing))
    axs[2].set_yticks(np.arange(0, ny*dx, spacing))
    axs[2].tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=8)
    #axs[2].set_title(r"x-velocity slice", fontsize=fontsize)
    axs[2].set_xlabel(r"$x~$[{}]".format(unit), fontsize=fontsize)
    axs[2].set_ylabel(r"$z~$[{}]".format(unit), fontsize=fontsize)

    fig.tight_layout()
    
    # plot and save
    save = True
    if save:
        plt.savefig(pngdir + f"{i}_mach_slice.png", dpi=300)
    plt.close()

    print(f"Saving figure {i+1} of {len(d_gas)}.\n")
