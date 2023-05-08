from hconfig import *

date = input("\nDate: ")
nfile = int(input("\nFile number: "))
cat_in = input("\ncat?: ")
cat = None
n_procs = None
if cat_in.lower() == "true":
    cat = True
    n_procs = int(input("\nNumber of processes: "))
else:
    cat = False
    n_procs = 1

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
datadir = os.path.join(basedir, "hdf5/")
dnameout = dnamein = datadir

ns = nfile
ne = nfile

DE = False # set to True if Dual Energy flag was used
SCALAR = True # set to True if Scalar was used

# loop over the output times
for n in range(ns, ne+1):

  # open the output file for writing
  fileout = h5py.File(dnameout+str(n)+'_slice.h5', 'w')

  # loop over files for a given output time
  for i in range(0, n_procs):

    # open the input file for reading
    filein = h5py.File(dnamein+str(n)+'_slice.h5.'+str(i), 'r')
    # read in the header data from the input file
    head = filein.attrs

    # if it's the first input file, write the header attributes
    # and create the datasets in the output file
    if (i == 0):
      gamma = head['gamma']
      t = head['t']
      dt = head['dt']
      dx = head['dx']
      n_step = head['n_step']
      nx = head['dims'][0]
      ny = head['dims'][1]
      nz = head['dims'][2]
      length_unit = head["length_unit"]
      mass_unit = head["mass_unit"]
      time_unit = head["time_unit"]
      density_unit = head["density_unit"]
      velocity_unit = head["velocity_unit"]
      energy_unit = head["energy_unit"]

      fileout.attrs['gamma'] = gamma
      fileout.attrs['t'] = t
      fileout.attrs['dt'] = dt
      fileout.attrs['dx'] = dx
      fileout.attrs['n_step'] = n_step
      fileout.attrs['dims'] = [nx, ny, nz]
      fileout.attrs['length_unit'] = length_unit
      fileout.attrs['time_unit'] = time_unit
      fileout.attrs['mass_unit'] = mass_unit
      fileout.attrs['density_unit'] = density_unit
      fileout.attrs['velocity_unit'] = velocity_unit
      fileout.attrs['energy_unit'] = energy_unit

      d_xy = np.zeros((nx,ny))
      d_xz = np.zeros((nx,nz))
      d_yz = np.zeros((ny,nz))
      mx_xy = np.zeros((nx,ny))
      mx_xz = np.zeros((nx,nz))
      mx_yz = np.zeros((ny,nz))
      my_xy = np.zeros((nx,ny))
      my_xz = np.zeros((nx,nz))
      my_yz = np.zeros((ny,nz))
      mz_xy = np.zeros((nx,ny))
      mz_xz = np.zeros((nx,nz))
      mz_yz = np.zeros((ny,nz))
      E_xy = np.zeros((nx,ny))
      E_xz = np.zeros((nx,nz))
      E_yz = np.zeros((ny,nz))
      if DE:
       GE_xy = np.zeros((nx,ny))
       GE_xz = np.zeros((nx,nz))
       GE_yz = np.zeros((ny,nz))
      if SCALAR:
       scalar_xy = np.zeros((nx,ny))
       scalar_xz = np.zeros((nx,nz))
       scalar_yz = np.zeros((ny,nz))

    # write data from individual processor file to
    # correct location in concatenated file
    nxl = head['dims_local'][0]
    nyl = head['dims_local'][1]
    nzl = head['dims_local'][2]
    xs = head['offset'][0]
    ys = head['offset'][1]
    zs = head['offset'][2]

    d_xy[xs:xs+nxl,ys:ys+nyl] += filein['d_xy']
    d_xz[xs:xs+nxl,zs:zs+nzl] += filein['d_xz']
    d_yz[ys:ys+nyl,zs:zs+nzl] += filein['d_yz']
    mx_xy[xs:xs+nxl,ys:ys+nyl] += filein['mx_xy']
    mx_xz[xs:xs+nxl,zs:zs+nzl] += filein['mx_xz']
    mx_yz[ys:ys+nyl,zs:zs+nzl] += filein['mx_yz']
    my_xy[xs:xs+nxl,ys:ys+nyl] += filein['my_xy']
    my_xz[xs:xs+nxl,zs:zs+nzl] += filein['my_xz']
    my_yz[ys:ys+nyl,zs:zs+nzl] += filein['my_yz']
    mz_xy[xs:xs+nxl,ys:ys+nyl] += filein['mz_xy']
    mz_xz[xs:xs+nxl,zs:zs+nzl] += filein['mz_xz']
    mz_yz[ys:ys+nyl,zs:zs+nzl] += filein['mz_yz']
    E_xy[xs:xs+nxl,ys:ys+nyl] += filein['E_xy']
    E_xz[xs:xs+nxl,zs:zs+nzl] += filein['E_xz']
    E_yz[ys:ys+nyl,zs:zs+nzl] += filein['E_yz']
    if DE:
      GE_xy[xs:xs+nxl,ys:ys+nyl] += filein['GE_xy']
      GE_xz[xs:xs+nxl,zs:zs+nzl] += filein['GE_xz']
      GE_yz[ys:ys+nyl,zs:zs+nzl] += filein['GE_yz']
    if SCALAR:
      scalar_xy[xs:xs+nxl,ys:ys+nyl] += filein['scalar_xy']
      scalar_xz[xs:xs+nxl,zs:zs+nzl] += filein['scalar_xz']
      scalar_yz[ys:ys+nyl,zs:zs+nzl] += filein['scalar_yz']

    filein.close()

  # wrte out the new datasets
  fileout.create_dataset('d_xy', data=d_xy)
  fileout.create_dataset('d_xz', data=d_xz)
  fileout.create_dataset('d_yz', data=d_yz)
  fileout.create_dataset('mx_xy', data=mx_xy)
  fileout.create_dataset('mx_xz', data=mx_xz)
  fileout.create_dataset('mx_yz', data=mx_yz)
  fileout.create_dataset('my_xy', data=my_xy)
  fileout.create_dataset('my_xz', data=my_xz)
  fileout.create_dataset('my_yz', data=my_yz)
  fileout.create_dataset('mz_xy', data=mz_xy)
  fileout.create_dataset('mz_xz', data=mz_xz)
  fileout.create_dataset('mz_yz', data=mz_yz)
  fileout.create_dataset('E_xy', data=E_xy)
  fileout.create_dataset('E_xz', data=E_xz)
  fileout.create_dataset('E_yz', data=E_yz)
  if DE:
    fileout.create_dataset('GE_xy', data=GE_xy)
    fileout.create_dataset('GE_xz', data=GE_xz)
    fileout.create_dataset('GE_yz', data=GE_yz)
  if SCALAR:
    fileout.create_dataset('scalar_xy', data=scalar_xy)
    fileout.create_dataset('scalar_xz', data=scalar_xz)
    fileout.create_dataset('scalar_yz', data=scalar_yz)

  fileout.close()

data = ReadHDF5(datadir, nscalar=1, slice="xy", fnum=nfile, cat=cat)
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
