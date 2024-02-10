import h5py
import numpy as np
import os

###############################
date = "2024-01-29"
ns = 0
ne = 5
n_procs = 2
dust = True
###############################

###############################
crc = True
frontier = False
###############################

########## data type ############
debugging = False
cloud_wind = False
testing = True
#################################

if crc:
  if debugging:
      basedir = f"/ix/eschneider/helena/data/debugging/{date}/"
  if cloud_wind:
      basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
  if testing:
      basedir = f"/ix/eschneider/helena/data/testing/{date}/"
  dnamein = os.path.join(basedir, "hdf5/raw/")
  dnameout = os.path.join(basedir, "hdf5/proj/")
if frontier:
  basedir = f"/lustre/orion/ast181/scratch/helenarichie/{date}/"
  dnamein = os.path.join(basedir, "hdf5/raw/")
  dnameout = os.path.join(basedir, "hdf5/proj/")

# loop over the output times
for n in range(ns, ne+1):

  # open the output file for writing
  fileout = h5py.File(dnameout+str(n)+'_proj.h5', 'w')

  # loop over files for a given output time
  for i in range(0, n_procs):

    # open the input file for reading
    filein = h5py.File(dnamein+str(n)+'_proj.h5.'+str(i), 'r')
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
      T_xy = np.zeros((nx,ny))
      T_xz = np.zeros((nx,nz))
      if dust:
        dust_xy = np.zeros((nx,ny))
        dust_xz = np.zeros((nx,nz))

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
    T_xy[xs:xs+nxl,ys:ys+nyl] += filein['T_xy']
    # T_xz[xs:xs+nxl,zs:zs+nzl] += filein['T_xz'] this wasn't working for some reason?
    if dust:
      dust_xy[xs:xs+nxl,ys:ys+nyl] += filein['d_dust_xy']
      #dust_xz[xs:xs+nxl,zs:zs+nzl] += filein['d_dust_xz']

    filein.close()

  # wrte out the new datasets
  fileout.create_dataset('d_xy', data=d_xy)
  fileout.create_dataset('d_xz', data=d_xz)
  fileout.create_dataset('T_xy', data=T_xy)
  fileout.create_dataset('T_xz', data=T_xz)
  if dust:
    fileout.create_dataset('d_dust_xy', data=dust_xy)
    fileout.create_dataset('d_dust_xz', data=dust_xz)

  fileout.close()
