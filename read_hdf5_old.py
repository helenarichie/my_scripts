import h5py
import numpy as np
import glob
import os
import re
import time

class ReadHDF5:
    """
    
    Reads in Cholla HDF5 file output.

    """

    def __init__(self, path, DE=False, dust=False, nscalar=0, slice=None, proj=None, fnum=None, cat=False):
        """
        Parameters
        ----------
        path : str
            Relative or global path to Cholla full-volume HDF5 file output.
        DE : Bool, optional
            Specifies whether dual energy was used. Default is false.
        nscalar : int, optional
            Specifies number of passive scalars Cholla was run with. Defualt 
            is 0.
        """
        self.mu = 0.6 # mean molecular weight
        self.MP = 1.672622e-24 # proton mass, g
        self.KB = 1.380658e-16 # Boltzmann constant, cm^2 g s^-2 K^-1

        self.path = path
        self.slice = slice
        self.proj = proj

        if slice != None:
            self.head, self.conserved = self.read_hdf5_slice(DE, dust, nscalar, slice, fnum, cat)
        elif proj != None:
            self.head, self.conserved = self.read_hdf5_proj(proj, fnum, cat)
        else:
            self.head, self.conserved = self.read_hdf5_file(DE, dust, nscalar, fnum, cat)

        print("\nData loaded.\n")


    def read_hdf5_file(self, DE=False, dust=False, nscalar=0, fnum=None, cat=False):
        files_unsorted = []
        if fnum == None:
            if cat == False:
                files_unsorted = glob.glob(os.path.join(self.path, "*.h5.0"))
            else:
                files_unsorted = glob.glob(os.path.join(self.path, "*.h5"))
        else:
            if cat == False:
                files_unsorted = glob.glob(os.path.join(self.path, f"{fnum}.h5.0"))
            else:
                files_unsorted = glob.glob(os.path.join(self.path, f"{fnum}.h5"))
        
        if len(files_unsorted) == 0:
            raise Exception("No HDF5 files found in this directory.")
        
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        files = sorted(files_unsorted, key=alphanum_key)
    
        print("Reading in full-volume data...\n")

        # get values of header attributes that don't change
        f0 = h5py.File(files[0])
        head0 = f0.attrs

        dims = head0["dims"]
        dx = head0["dx"]
        gamma = head0["gamma"][0]
        length_unit = head0["length_unit"][0]
        time_unit = head0["time_unit"][0]
        mass_unit = head0["mass_unit"][0]
        density_unit = head0["density_unit"][0]
        velocity_unit = head0["velocity_unit"][0]
        energy_unit = head0["energy_unit"][0]

        f0.close()

        # initialize numpy arrays to store data
        array_dims = np.insert(dims, 0, len(files))

        # header attributes
        t = np.zeros(len(files))
        n_step = np.zeros(len(files))

        density = np.zeros(array_dims)
        momentum_x = np.zeros(array_dims)
        momentum_y = np.zeros(array_dims)
        momentum_z = np.zeros(array_dims)
        energy = np.zeros(array_dims)
        gas_energy, dust_density, scalar = 0, 0, 0
        if DE:
            gas_energy = np.zeros(array_dims)
        if dust:
            dust_density = np.zeros(array_dims)
        if nscalar > 0:
            scalar = np.zeros(np.insert(array_dims, 0, nscalar))
        
        for i, o_file in enumerate(files):
            start = time.time()
            
            f = h5py.File(os.path.join(self.path, o_file))
            head = f.attrs

            # header attributes
            t[i] = head["t"]
            n_step[i] = head["n_step"]

            # conserved variables
            density[i] = f["density"]
            momentum_x[i] = f["momentum_x"]
            momentum_y[i] = f["momentum_y"]
            momentum_z[i] = f["momentum_z"]
            energy[i] = f["Energy"]
            if DE:
                gas_energy[i] = f["GasEnergy"]
            if nscalar > 0:
                for j in range(nscalar):
                    scalar[j][i] = f["scalar" + str(j)]
            if dust:
                dust_density[i] = f["dust_density"]
            
            f.close()
            print(os.path.basename(o_file), f"--- {round(time.time() - start, 3)} s")

        conserved = {
            "density" : density,
            "momentum_x" : momentum_x,
            "momentum_y" : momentum_y,
            "momentum_z" : momentum_z,
            "energy" : energy,
            "gas_energy" : gas_energy,
            "dust_density" : dust_density,
        }
        for i in range(nscalar):
            conserved["scalar" + str(i)] = scalar[i]

        head = {
            "t" : t,
            "n_step" : n_step,
            "dims" : dims,
            "dx" : dx,
            "gamma" : gamma,
            "time_unit" : time_unit,
            "mass_unit" : mass_unit,
            "length_unit" : length_unit,
            "density_unit" : density_unit,
            "velocity_unit" : velocity_unit,
            "energy_unit" : energy_unit,
        }

        return head, conserved

    
    def read_hdf5_slice(self, DE=False, dust=False, nscalar=0, slice=None, fnum=None, cat=False):
        if not ((slice == "xy") or (slice == "xz") or (slice == "yz")):
            raise Exception('Slices must be specified as "xy", "xz", or "yz".')

        files_unsorted = []
        if fnum == None:
            if cat == False:
                files_unsorted = glob.glob(os.path.join(self.path, "*_slice.h5.0"))
            else:
                files_unsorted = glob.glob(os.path.join(self.path, "*_slice.h5"))
        else:
            if cat == False:
                files_unsorted = glob.glob(os.path.join(self.path, f"{fnum}_slice.h5.0"))
            else:
                files_unsorted = glob.glob(os.path.join(self.path, f"{fnum}_slice.h5"))
        if len(files_unsorted) == 0:
            raise Exception("No HDF5 files found in this directory.")
        
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        files = sorted(files_unsorted, key=alphanum_key)
    
        print("Reading in slice data...\n")

        # get values of header attributes that don't change
        f0 = h5py.File(files[0])
        head0 = f0.attrs

        print(list(f0.keys()))

        (nx, ny, nz) = head0["dims"]
        dx = head0["dx"]
        gamma = head0["gamma"][0]
        length_unit = head0["length_unit"][0]
        time_unit = head0["time_unit"][0]
        mass_unit = head0["mass_unit"][0]
        density_unit = head0["density_unit"][0]
        velocity_unit = head0["velocity_unit"][0]
        energy_unit = head0["energy_unit"][0]

        f0.close()

        # initialize numpy arrays to store data
        array_dims = None
        if slice == "xy":
            array_dims = (len(files), nx, ny)
        if slice == "xz":
            array_dims = (len(files), nx, nz)
        if slice == "yz":
            array_dims = (len(files), ny, nz)

        # header attributes
        t = np.zeros(len(files))
        n_step = np.zeros(len(files))

        density = np.zeros(array_dims)
        momentum_x = np.zeros(array_dims)
        momentum_y = np.zeros(array_dims)
        momentum_z = np.zeros(array_dims)
        energy = np.zeros(array_dims)
        gas_energy, dust_density, scalar = 0, 0, 0
        if DE:
            gas_energy = np.zeros(array_dims)
        if dust:
            dust_density = np.zeros(array_dims)
        if nscalar > 0:
            scalar = np.zeros(np.insert(array_dims, 0, nscalar))
        
        for i, o_file in enumerate(files):
            start = time.time()
            
            f = h5py.File(os.path.join(self.path, o_file))
            head = f.attrs

            # header attributes
            t[i] = head["t"]
            n_step[i] = head["n_step"]

            # conserved variables
            density[i] = f[f"d_{slice}"]
            momentum_x[i] = f[f"mx_{slice}"]
            momentum_y[i] = f[f"my_{slice}"]
            momentum_z[i] = f[f"mz_{slice}"]
            energy[i] = f[f"E_{slice}"]
            if DE:
                gas_energy[i] = f[f"GE_{slice}"]
            if dust:
                dust_density[i] = f[f"d_dust_{slice}"]
            if nscalar > 0:
                for j in range(nscalar):
                    scalar[j][i] = f[f"scalar_{slice}"]
            
            f.close()
            print(os.path.basename(o_file), f"--- {round(time.time() - start, 3)} s")

        conserved = {
            "density" : density,
            "momentum_x" : momentum_x,
            "momentum_y" : momentum_y,
            "momentum_z" : momentum_z,
            "energy" : energy,
            "gas_energy" : gas_energy,
            "dust_density" : dust_density,
        }
        for i in range(nscalar):
            conserved["scalar" + str(i)] = scalar[i]

        head = {
            "t" : t,
            "n_step" : n_step,
            "dims" : (nx, ny, nz),
            "dx" : dx,
            "gamma" : gamma,
            "time_unit" : time_unit,
            "mass_unit" : mass_unit,
            "length_unit" : length_unit,
            "density_unit" : density_unit,
            "velocity_unit" : velocity_unit,
            "energy_unit" : energy_unit,
        }

        return head, conserved


    def read_hdf5_proj(self, proj=None, fnum=None, cat=False):
        files_unsorted = []
        if fnum == None:
            if cat == False:
                files_unsorted = glob.glob(os.path.join(self.path, "*_proj.h5.0"))
            else:
                files_unsorted = glob.glob(os.path.join(self.path, "*_proj.h5"))
        else:
            if cat == False:
                files_unsorted = glob.glob(os.path.join(self.path, f"{fnum}_proj.h5.0"))
            else:
                files_unsorted = glob.glob(os.path.join(self.path, f"{fnum}_proj.h5"))
        if len(files_unsorted) == 0:
            raise Exception("No HDF5 files found in this directory.")
        
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        files = sorted(files_unsorted, key=alphanum_key)

        print("Reading in projection data...\n")

        # get values of header attributes that don't change
        f0 = h5py.File(files[0])
        head0 = f0.attrs

        (nx, ny, nz) = head0["dims"]
        dx = head0["dx"]
        gamma = head0["gamma"][0]
        length_unit = head0["length_unit"][0]
        time_unit = head0["time_unit"][0]
        mass_unit = head0["mass_unit"][0]
        density_unit = head0["density_unit"][0]
        velocity_unit = head0["velocity_unit"][0]
        energy_unit = head0["energy_unit"][0]

        f0.close()

        # initialize numpy arrays to store data
        array_dims = None
        if proj == "xy":
            array_dims = (len(files), nx, ny)
        if proj == "xz":
            array_dims = (len(files), nx, nz)

        # header attributes
        t = np.zeros(len(files))
        n_step = np.zeros(len(files))

        density = np.zeros(array_dims)
        temperature = np.zeros(array_dims)
        
        for i, o_file in enumerate(files):
            start = time.time()
            
            f = h5py.File(os.path.join(self.path, o_file))
            head = f.attrs

            # header attributes
            t[i] = head["t"]
            n_step[i] = head["n_step"]

            # conserved variables
            density[i] = f[f"d_{proj}"]
            temperature[i] = f[f"T_{proj}"]
            
            f.close()
            print(os.path.basename(o_file), f"--- {round(time.time() - start, 3)} s")

        conserved = {
            "density" : density,
            "temperature" : temperature,
        }

        head = {
            "t" : t,
            "n_step" : n_step,
            "dims" : (nx, ny, nz),
            "dx" : dx,
            "gamma" : gamma,
            "time_unit" : time_unit,
            "mass_unit" : mass_unit,
            "length_unit" : length_unit,
            "density_unit" : density_unit,
            "velocity_unit" : velocity_unit,
            "energy_unit" : energy_unit,
        }

        return head, conserved

    
    # density
    def d_cgs(self):
        return self.conserved["density"] * self.head["density_unit"]

    # x, y, and z momentum
    def mx_cgs(self):
        return self.conserved["momentum_x"] * self.head["mass_unit"] * self.head["velocity_unit"]

    def my_cgs(self):
        return self.conserved["momentum_y"] * self.head["mass_unit"] * self.head["velocity_unit"]
    
    def mz_cgs(self):
        return self.conserved["momentum_z"] * self.head["mass_unit"] * self.head["velocity_unit"]

    # energy and gas energy
    def energy_cgs(self):
        return self.conserved["energy"] * self.head["energy_unit"]
        
    def gas_energy_cgs(self, DE):
        if DE:
            return self.conserved["GasEnergy"] * self.head["energy_unit"]
        else:
            print("No gas energy field found")
            return None

    # x, y, and z velocity
    def vx_cgs(self):
        return self.conserved["momentum_x"] / self.conserved["density"] * self.head["velocity_unit"]
    
    def vy_cgs(self):
        return self.conserved["momentum_y"] / self.conserved["density"] * self.head["velocity_unit"]

    def vz_cgs(self):
        return self.conserved["momentum_z"] / self.conserved["density"] * self.head["velocity_unit"]

    # temperature
    def T(self):
        return ((self.head["gamma"] - 1) / (((self.conserved["density"] * self.head["density_unit"]) / (self.mu*self.MP))*self.KB)) * ((self.conserved["energy"] * self.head["energy_unit"]) - (1/2 * (self.conserved["density"] * self.head["density_unit"]) * ((self.conserved["momentum_x"] / self.conserved["density"] * self.head["velocity_unit"])**2 + (self.conserved["momentum_y"] / self.conserved["density"] * self.head["velocity_unit"])**2 + (self.conserved["momentum_z"] / self.conserved["density"] * self.head["velocity_unit"])**2)))
    
    # number density
    def n_cgs(self):
        return self.conserved["density"] * self.head["density_unit"] / (self.mu*self.MP)

    # time
    def t_cgs(self):
        return self.head["t"] * self.head["time_unit"]
    
    # cell dimensions
    def dx_cgs(self):
        return self.head["dx"] * self.head["length_unit"]

    
