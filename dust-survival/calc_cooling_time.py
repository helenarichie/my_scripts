import matplotlib.pyplot as plt
import numpy as np
import h5py

kb = 1.3807e-16

# f = h5py.File('/Users/helenarichie/Documents/Grad School/research/data/cooling/z_0.000.hdf5','r')
f = h5py.File("/ix/eschneider/helena/data/cooling/z_0.000.hdf5")
T = np.array(f['Solar/Temperature_bins/'])
L = np.array(f['Solar/Net_cooling/'])
f.close()

H = 10**np.arange(0, 4, 0.01)
T1 = 10**np.arange(0, 4, 0.01)
T2 = np.arange(4, 9, 0.01)
cool = np.zeros(np.size(T2))
cool1 = 2e-26 * (1e7 * np.exp(-1.148e5 / (T1+1000)) + 1.4e-2 * np.sqrt(T1) * np.exp(-92.0/T1))

i = 0
for temp in T2:
#  if (temp > 4.0 and temp < 4.43):
#    cool[i] = np.power(10.0, (-15.0 * (temp - 4.30) * (temp - 4.30) - 21.85))
  if (temp > 4.0 and temp < 5.9):
    cool[i] = np.power(10.0, (-1.3 * (temp - 5.25) * (temp - 5.25) - 21.25))
  if (temp > 5.9 and temp < 7.4):
    cool[i] = np.power(10.0, (0.7 * (temp - 7.1) * (temp - 7.1) - 22.8))
  if (temp > 7.4):
    cool[i] = np.power(10.0, (0.45*temp - 26.065))
  i += 1

LK1 = 1e-24
al1 = 2.3
Tk1 = 8e3
c1 = LK1*np.power(10**T2 / Tk1, al1)
LK2 = 6e-22
al2 = -0.65
Tk2 = 1e5
c2 = LK2*np.power(10**T2 / Tk2, al2)
LK3 = 1.75e-23
al3 = 0.45
Tk3 = 2e7
c3 = LK3*np.power(10**T2 / Tk3, al3)
Tref = 1e9
Lref = 1e-24
al = -0.75
c = Lref*np.power((10**T2/Tref), al)

def cooling_time(argn, argT):
    wh_T = np.argmin(abs(T[101:]-argT))
    # lambd = L[101:,80][wh_T]
    lambd = 4.102041029866097e-23
    return 3 * argn * kb * argT / (2 * lambd * argn ** 2)

print(np.sqrt(10*1e-2))
print(np.sqrt(3e3*3e6))

print(f"Cooling time: {cooling_time(np.sqrt(1e-2), np.sqrt(3e6))/3.154e7*1e-3} kyr")
