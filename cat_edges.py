import h5py
import numpy as np
import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

###############################
date = "2024-02-11"
ns = 0
ne = 1200
n_procs = 6
DE = True
SCALAR = True
###############################

###############################
crc = True
frontier = False
###############################

###############################
testing = False
cloud_wind = True
###############################

if crc:
  if cloud_wind:
    basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
  if testing:
    basedir = f"/ix/eschneider/helena/data/testing/{date}/"
  dnamein = os.path.join(basedir, "hdf5/raw/")
  dnameout = os.path.join(basedir, "hdf5/edges/")
if frontier:
  basedir = f"/lustre/orion/ast181/scratch/helenarichie/{date}/"
  dnamein = os.path.join(basedir, "hdf5/raw")
  dnameout = os.path.join(basedir, "hdf5/edges/")

# loop over the output times
for n in range(ns, ne+1):

  print(f"Concatenating boundaries {n} of {ne}.\n")

  # open the output file for writing
  fileout = h5py.File(dnameout+str(n)+'_edges.h5', 'w')

  # loop over files for a given output time
  for i in range(0, n_procs):

    # open the input file for reading
    filein = h5py.File(dnamein+str(n)+'_edges.h5.'+str(i), 'r')
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

      d_minus_xy = np.zeros((nx,ny))
      d_minus_xz = np.zeros((nx,nz))
      d_minus_yz = np.zeros((ny,nz))
      d_plus_xy = np.zeros((nx,ny))
      d_plus_xz = np.zeros((nx,nz))
      d_plus_yz = np.zeros((ny,nz))
      mx_minus_xy = np.zeros((nx,ny))
      mx_minus_xz = np.zeros((nx,nz))
      mx_minus_yz = np.zeros((ny,nz))
      mx_plus_xy = np.zeros((nx,ny))
      mx_plus_xz = np.zeros((nx,nz))
      mx_plus_yz = np.zeros((ny,nz))
      my_minus_xy = np.zeros((nx,ny))
      my_minus_xz = np.zeros((nx,nz))
      my_minus_yz = np.zeros((ny,nz))
      my_plus_xy = np.zeros((nx,ny))
      my_plus_xz = np.zeros((nx,nz))
      my_plus_yz = np.zeros((ny,nz))
      mz_minus_xy = np.zeros((nx,ny))
      mz_minus_xz = np.zeros((nx,nz))
      mz_minus_yz = np.zeros((ny,nz))
      mz_plus_xy = np.zeros((nx,ny))
      mz_plus_xz = np.zeros((nx,nz))
      mz_plus_yz = np.zeros((ny,nz))
      E_minus_xy = np.zeros((nx,ny))
      E_minus_xz = np.zeros((nx,nz))
      E_minus_yz = np.zeros((ny,nz))
      E_plus_xy = np.zeros((nx,ny))
      E_plus_xz = np.zeros((nx,nz))
      E_plus_yz = np.zeros((ny,nz))
      if DE:
       GE_minus_xy = np.zeros((nx,ny))
       GE_minus_xz = np.zeros((nx,nz))
       GE_minus_yz = np.zeros((ny,nz))
       GE_plus_xy = np.zeros((nx,ny))
       GE_plus_xz = np.zeros((nx,nz))
       GE_plus_yz = np.zeros((ny,nz))
      if SCALAR:
       scalar_minus_xy = np.zeros((nx,ny))
       scalar_minus_xz = np.zeros((nx,nz))
       scalar_minus_yz = np.zeros((ny,nz))
       scalar_plus_xy = np.zeros((nx,ny))
       scalar_plus_xz = np.zeros((nx,nz))
       scalar_plus_yz = np.zeros((ny,nz))

    # write data from individual processor file to
    # correct location in concatenated file
    nxl = head['dims_local'][0]
    nyl = head['dims_local'][1]
    nzl = head['dims_local'][2]
    xs = head['offset'][0]
    ys = head['offset'][1]
    zs = head['offset'][2]

    d_minus_xy[xs:xs+nxl,ys:ys+nyl] += filein['d_minus_xy']
    d_minus_xz[xs:xs+nxl,zs:zs+nzl] += filein['d_minus_xz']
    d_minus_yz[ys:ys+nyl,zs:zs+nzl] += filein['d_minus_yz']
    d_plus_xy[xs:xs+nxl,ys:ys+nyl] += filein['d_plus_xy']
    d_plus_xz[xs:xs+nxl,zs:zs+nzl] += filein['d_plus_xz']
    d_plus_yz[ys:ys+nyl,zs:zs+nzl] += filein['d_plus_yz']
    mx_minus_xy[xs:xs+nxl,ys:ys+nyl] += filein['mx_minus_xy']
    mx_minus_xz[xs:xs+nxl,zs:zs+nzl] += filein['mx_minus_xz']
    mx_minus_yz[ys:ys+nyl,zs:zs+nzl] += filein['mx_minus_yz']
    mx_plus_xy[xs:xs+nxl,ys:ys+nyl] += filein['mx_plus_xy']
    mx_plus_xz[xs:xs+nxl,zs:zs+nzl] += filein['mx_plus_xz']
    mx_plus_yz[ys:ys+nyl,zs:zs+nzl] += filein['mx_plus_yz']
    my_minus_xy[xs:xs+nxl,ys:ys+nyl] += filein['my_minus_xy']
    my_minus_xz[xs:xs+nxl,zs:zs+nzl] += filein['my_minus_xz']
    my_minus_yz[ys:ys+nyl,zs:zs+nzl] += filein['my_minus_yz']
    my_plus_xy[xs:xs+nxl,ys:ys+nyl] += filein['my_plus_xy']
    my_plus_xz[xs:xs+nxl,zs:zs+nzl] += filein['my_plus_xz']
    my_plus_yz[ys:ys+nyl,zs:zs+nzl] += filein['my_plus_yz']
    mz_minus_xy[xs:xs+nxl,ys:ys+nyl] += filein['mz_minus_xy']
    mz_minus_xz[xs:xs+nxl,zs:zs+nzl] += filein['mz_minus_xz']
    mz_minus_yz[ys:ys+nyl,zs:zs+nzl] += filein['mz_minus_yz']
    mz_plus_xy[xs:xs+nxl,ys:ys+nyl] += filein['mz_plus_xy']
    mz_plus_xz[xs:xs+nxl,zs:zs+nzl] += filein['mz_plus_xz']
    mz_plus_yz[ys:ys+nyl,zs:zs+nzl] += filein['mz_plus_yz']
    E_minus_xy[xs:xs+nxl,ys:ys+nyl] += filein['E_minus_xy']
    E_minus_xz[xs:xs+nxl,zs:zs+nzl] += filein['E_minus_xz']
    E_minus_yz[ys:ys+nyl,zs:zs+nzl] += filein['E_minus_yz']
    E_plus_xy[xs:xs+nxl,ys:ys+nyl] += filein['E_plus_xy']
    E_plus_xz[xs:xs+nxl,zs:zs+nzl] += filein['E_plus_xz']
    E_plus_yz[ys:ys+nyl,zs:zs+nzl] += filein['E_plus_yz']
    if DE:
      GE_minus_xy[xs:xs+nxl,ys:ys+nyl] += filein['GE_minus_xy']
      GE_minus_xz[xs:xs+nxl,zs:zs+nzl] += filein['GE_minus_xz']
      GE_minus_yz[ys:ys+nyl,zs:zs+nzl] += filein['GE_minus_yz']
      GE_plus_xy[xs:xs+nxl,ys:ys+nyl] += filein['GE_plus_xy']
      GE_plus_xz[xs:xs+nxl,zs:zs+nzl] += filein['GE_plus_xz']
      GE_plus_yz[ys:ys+nyl,zs:zs+nzl] += filein['GE_plus_yz']
    if SCALAR:
      scalar_minus_xy[xs:xs+nxl,ys:ys+nyl] += filein['scalar_minus_xy']
      scalar_minus_xz[xs:xs+nxl,zs:zs+nzl] += filein['scalar_minus_xz']
      scalar_minus_yz[ys:ys+nyl,zs:zs+nzl] += filein['scalar_minus_yz']
      scalar_plus_xy[xs:xs+nxl,ys:ys+nyl] += filein['scalar_plus_xy']
      scalar_plus_xz[xs:xs+nxl,zs:zs+nzl] += filein['scalar_plus_xz']
      scalar_plus_yz[ys:ys+nyl,zs:zs+nzl] += filein['scalar_plus_yz']

    filein.close()

  # wrte out the new datasets
  fileout.create_dataset('d_minus_xy', data=d_minus_xy)
  fileout.create_dataset('d_minus_xz', data=d_minus_xz)
  fileout.create_dataset('d_minus_yz', data=d_minus_yz)
  fileout.create_dataset('d_plus_xy', data=d_plus_xy)
  fileout.create_dataset('d_plus_xz', data=d_plus_xz)
  fileout.create_dataset('d_plus_yz', data=d_plus_yz)
  fileout.create_dataset('mx_minus_xy', data=mx_minus_xy)
  fileout.create_dataset('mx_minus_xz', data=mx_minus_xz)
  fileout.create_dataset('mx_minus_yz', data=mx_minus_yz)
  fileout.create_dataset('mx_plus_xy', data=mx_plus_xy)
  fileout.create_dataset('mx_plus_xz', data=mx_plus_xz)
  fileout.create_dataset('mx_plus_yz', data=mx_plus_yz)
  fileout.create_dataset('my_minus_xy', data=my_minus_xy)
  fileout.create_dataset('my_minus_xz', data=my_minus_xz)
  fileout.create_dataset('my_minus_yz', data=my_minus_yz)
  fileout.create_dataset('my_plus_xy', data=my_plus_xy)
  fileout.create_dataset('my_plus_xz', data=my_plus_xz)
  fileout.create_dataset('my_plus_yz', data=my_plus_yz)
  fileout.create_dataset('mz_minus_xy', data=mz_minus_xy)
  fileout.create_dataset('mz_minus_xz', data=mz_minus_xz)
  fileout.create_dataset('mz_minus_yz', data=mz_minus_yz)
  fileout.create_dataset('mz_plus_xy', data=mz_plus_xy)
  fileout.create_dataset('mz_plus_xz', data=mz_plus_xz)
  fileout.create_dataset('mz_plus_yz', data=mz_plus_yz)
  fileout.create_dataset('E_minus_xy', data=E_minus_xy)
  fileout.create_dataset('E_minus_xz', data=E_minus_xz)
  fileout.create_dataset('E_minus_yz', data=E_minus_yz)
  fileout.create_dataset('E_plus_xy', data=E_plus_xy)
  fileout.create_dataset('E_plus_xz', data=E_plus_xz)
  fileout.create_dataset('E_plus_yz', data=E_plus_yz)
  if DE:
    fileout.create_dataset('GE_minus_xy', data=GE_minus_xy)
    fileout.create_dataset('GE_minus_xz', data=GE_minus_xz)
    fileout.create_dataset('GE_minus_yz', data=GE_minus_yz)
    fileout.create_dataset('GE_plus_xy', data=GE_plus_xy)
    fileout.create_dataset('GE_plus_xz', data=GE_plus_xz)
    fileout.create_dataset('GE_plus_yz', data=GE_plus_yz)
  if SCALAR:
    fileout.create_dataset('scalar_minus_xy', data=scalar_minus_xy)
    fileout.create_dataset('scalar_minus_xz', data=scalar_minus_xz)
    fileout.create_dataset('scalar_minus_yz', data=scalar_minus_yz)
    fileout.create_dataset('scalar_plus_xy', data=scalar_plus_xy)
    fileout.create_dataset('scalar_plus_xz', data=scalar_plus_xz)
    fileout.create_dataset('scalar_plus_yz', data=scalar_plus_yz)

  fileout.close()
