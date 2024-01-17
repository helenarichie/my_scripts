#!/usr/bin/env python3
# Example file for concatenating 3D hdf5 datasets
import sys
sys.path.insert(0, "/ix/eschneider/helena/code/my_scripts/")
import h5py
import numpy as np
import os
from hconfig import *
from csv import writer

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

#######################
ns = 0
ne = 1104
n_hydro = 1
n_proc = 8 # number of processors that did the calculations
scalar = False
dust = True
date = "2024-01-11"
istart = 0
cat = True
rho_cl_i = 1e-23  # needed to index cloud material
cutoff = 0.05*rho_cl_i*density_conversion # 5% of initial density, M_sun/kpc^3
#######################

istart = 0*n_proc
iend = 1*n_proc

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
dnamein = os.path.join(basedir, "hdf5/raw/")
dnameout = fulldir = os.path.join(basedir, "hdf5/full/")
projdir = os.path.join(basedir, "hdf5/proj/")
csvdir = os.path.join(basedir, "csv/")

data = ReadHDF5(projdir, fnum=0, proj="xy", dust=True, cat=cat)
head = data.head
nx, ny, nz = head["dims"]
dx = head["dx"] # kpc

# x, y, z
areas = [ny*nz*dx**2, nx*nz*dx**2, nx*ny*dx**2]
indices = [0, -1]
fluids = ["gas", "dust"]

if istart == 0:
    # create files to write quantities for cloud (dense gas only), total gas, and dust
    f = open(os.path.join(csvdir, "rate_cl.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "mass_cl.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "rate_gas.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "mass_gas.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "mass_dust.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "rate_dust.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "rho_d_tot.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "rho_cl_tot.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "T_cl_avg.csv"), "w")
    f.close()
    f = open(os.path.join(csvdir, "t_arr.csv"), "w")

# loop over outputs
for n in range(ns, ne+1):
  n *= n_hydro
  # loop over files for a given output
  for i in range(istart, iend):

    # open the output file for writing (don't overwrite if exists)
    fileout = h5py.File(dnameout+str(n)+'.h5', 'a')
    # open the input file for reading
    filein = h5py.File(dnamein+str(n)+'.h5.'+str(i), 'r')
    # read in the header data from the input file
    head = filein.attrs

    # if it's the first input file, write the header attributes
    # and create the datasets in the output file
    if (i == 0):
      nx = head['dims'][0]
      ny = head['dims'][1]
      nz = head['dims'][2]
      fileout.attrs['dims'] = [nx, ny, nz]
      fileout.attrs['dx'] = [head['dx'][0]]
      fileout.attrs['gamma'] = [head['gamma'][0]]
      fileout.attrs['t'] = [head['t'][0]]
      fileout.attrs['dt'] = [head['dt'][0]]
      fileout.attrs['n_step'] = [head['n_step'][0]]

      units = ['time_unit', 'mass_unit', 'length_unit', 'energy_unit', 'velocity_unit', 'density_unit']
      for unit in units:
        fileout.attrs[unit] = [head[unit][0]]

      d  = fileout.create_dataset("density", (nx, ny, nz), chunks=True, dtype=filein['density'].dtype)
      mx = fileout.create_dataset("momentum_x", (nx, ny, nz), chunks=True, dtype=filein['momentum_x'].dtype)
      my = fileout.create_dataset("momentum_y", (nx, ny, nz), chunks=True, dtype=filein['momentum_y'].dtype)
      mz = fileout.create_dataset("momentum_z", (nx, ny, nz), chunks=True, dtype=filein['momentum_z'].dtype)
      E  = fileout.create_dataset("Energy", (nx, ny, nz), chunks=True, dtype=filein['Energy'].dtype)
      if scalar:
        scalar0 = fileout.create_dataset('scalar0', (nx, ny, nz), chunks=True, dtype=filein['scalar0'].dtype)
      if dust:
        dust_density = fileout.create_dataset('dust_density', (nx, ny, nz), chunks=True, dtype=filein['dust_density'].dtype)
      try:
        GE = fileout.create_dataset("GasEnergy", (nx, ny, nz), chunks=True, dtype=filein['GasEnergy'].dtype)
      except KeyError:
        print('No Dual energy data present');
      try:
        [nx_mag, ny_mag, nz_mag] = head['magnetic_field_dims']
        bx = fileout.create_dataset("magnetic_x", (nx_mag, ny_mag, nz_mag), chunks=True, dtype=filein['magnetic_x'].dtype)
        by = fileout.create_dataset("magnetic_y", (nx_mag, ny_mag, nz_mag), chunks=True, dtype=filein['magnetic_y'].dtype)
        bz = fileout.create_dataset("magnetic_z", (nx_mag, ny_mag, nz_mag), chunks=True, dtype=filein['magnetic_z'].dtype)
      except KeyError:
        print('No magnetic field data present');

    # write data from individual processor file to
    # correct location in concatenated file
    nxl = head['dims_local'][0]
    nyl = head['dims_local'][1]
    nzl = head['dims_local'][2]
    xs = head['offset'][0]
    ys = head['offset'][1]
    zs = head['offset'][2]
    fileout['density'][xs:xs+nxl,ys:ys+nyl,zs:zs+nzl]  = filein['density']
    fileout['momentum_x'][xs:xs+nxl,ys:ys+nyl,zs:zs+nzl] = filein['momentum_x']
    fileout['momentum_y'][xs:xs+nxl,ys:ys+nyl,zs:zs+nzl] = filein['momentum_y']
    fileout['momentum_z'][xs:xs+nxl,ys:ys+nyl,zs:zs+nzl] = filein['momentum_z']
    fileout['Energy'][xs:xs+nxl,ys:ys+nyl,zs:zs+nzl]  = filein['Energy']
    if scalar:
      fileout['scalar0'][xs:xs+nxl,ys:ys+nyl,zs:zs+nzl] = filein['scalar0']
    if dust:
      fileout['dust_density'][xs:xs+nxl,ys:ys+nyl,zs:zs+nzl] = filein['dust_density']
    try:
      fileout['GasEnergy'][xs:xs+nxl,ys:ys+nyl,zs:zs+nzl] = filein['GasEnergy']
    except KeyError:
        print('No Dual energy data present');
    try:
      [nxl_mag, nyl_mag, nzl_mag] = head['magnetic_field_dims_local']
      fileout['magnetic_x'][xs:xs+nxl_mag,ys:ys+nyl_mag,zs:zs+nzl_mag] = filein['magnetic_x']
      fileout['magnetic_y'][xs:xs+nxl_mag,ys:ys+nyl_mag,zs:zs+nzl_mag] = filein['magnetic_y']
      fileout['magnetic_z'][xs:xs+nxl_mag,ys:ys+nyl_mag,zs:zs+nzl_mag] = filein['magnetic_z']
    except KeyError:
        print('No magnetic field data present');

    filein.close()

  d_gas = d * head["density_unit"]
  d_dust = dust_density * head["density_unit"]

  with open(os.path.join(csvdir, "t_arr.csv"), "a") as f:
    writer_obj = writer(f)
    writer_obj.writerow([head['t'][0]])
    f.close()

  # calculate and write rates and masses for cloud
  rates = get_rates(d_gas, cutoff)
  masses = get_masses(d_gas, cutoff)

  with open(os.path.join(csvdir, "rate_cl.csv"), "a") as f:
    writer_obj = writer(f)
    writer_obj.writerow(rates)
    f.close()

  with open(os.path.join(csvdir, "mass_cl.csv"), "a") as f:
    writer_obj = writer(f)
    writer_obj.writerow(masses)
    f.close()

    # calculate and write masses for dust
  masses = get_masses(d_dust, 1)

  with open(os.path.join(csvdir, "mass_dust.csv"), "a") as f:
    writer_obj = writer(f)
    writer_obj.writerow(masses)
    f.close()

  # calculate and write masses for dust
  rates = get_rates(d_dust, 1)

  with open(os.path.join(csvdir, "rate_dust.csv"), "a") as f:
    writer_obj = writer(f)
    writer_obj.writerow(rates)
    f.close()

  fileout.close()

def calc_mass_loss_rate(rho, v, area):
    return np.sum(rho * np.abs(v) * area)

def get_rates(dens, cutoff):
    rates = np.zeros(6)
    for i, area in enumerate(areas):
        for index in indices:
            if (index == 0) and (i == 0):
                d_i, v_i = dens[0, :, :], mx[0, :, :] / d_gas[0, :, :]
                mask = np.logical_and(d_i >= cutoff, v_i < 0) # mask for high-density, outflowing gas
                rates[0] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == 0) and (i == 1):
                d_i, v_i = dens[:, 0, :], my[:, 0, :] / d_gas[:, 0, :]
                mask = np.logical_and(d_i >= cutoff, v_i < 0)
                rates[1] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == 0) and (i == 2):
                d_i, v_i = dens[:, :, 0], mz[:, :, 0] / d_gas[:, :, 0]
                mask = np.logical_and(d_i >= cutoff, v_i < 0)
                rates[2] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == -1) and (i == 0):
                d_i, v_i = dens[-1, :, :], mx[-1, :, :] / (d_gas[-1, :, :])
                mask = np.logical_and(d_i >= cutoff, v_i > 0)
                rates[3] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == -1) and (i == 1):
                d_i, v_i = dens[:, -1, :], my[:, -1, :] / d_gas[:, -1, :]
                mask = np.logical_and(d_i >= cutoff, v_i > 0)
                rates[4] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)

            if (index == -1) and (i == 2):
                d_i, v_i = dens[:, :, -1], mz[:, :, -1] / d_gas[:, :, -1]
                mask = np.logical_and(d_i >= cutoff, v_i > 0)
                rates[5] = calc_mass_loss_rate(d_i[mask], v_i[mask], dx**2)
    
    return rates

def get_masses(dens, cutoff):
    masses = np.zeros(6)
    for i, area in enumerate(areas):
        for index in indices:
            if (index == 0) and (i == 0):
                d_i = dens[0, :, :]
                mask = (d_i >= cutoff)
                masses[0] = np.sum(d_i[mask]) * dx**3 # total density times total volume of cells

            if (index == 0) and (i == 1):
                d_i = dens[:, 0, :]
                mask = (d_i >= cutoff)
                masses[1] = np.sum(d_i[mask]) * dx**3

            if (index == 0) and (i == 2):
                d_i = dens[:, :, 0]
                mask = (d_i >= cutoff)
                masses[2] = np.sum(d_i[mask]) * dx**3

            if (index == -1) and (i == 0):
                d_i = dens[-1, :, :]
                mask = (d_i >= cutoff)
                masses[3] = np.sum(d_i[mask]) * dx**3

            if (index == -1) and (i == 1):
                d_i = dens[:, -1, :]
                mask = (d_i >= cutoff)
                masses[4] = np.sum(d_i[mask]) * dx**3

            if (index == -1) and (i == 2):
                d_i = dens[:, :, -1]
                mask = (d_i >= cutoff)
                masses[5] = np.sum(d_i[mask]) * dx**3
    
    return masses