# Script to calculate non-thermal sputtering timescales and drag timescales
# Use plasmapy to calculate Coulomb logarithm
from hconfig import *
e = 4.8032e-10  # cm^3/2 g^1/2 s^-1
gamma = 5/2
mu = 0.6
rho = 1e-24  # g/cm^3
P = 1e6 * KB
v_rel = 1e3 * 1e5  # cm/s
ln_coulomb = 8
z = 1  # ionization fraction, q/e
Z_gr = 36
a = 0.1  # micrometers
T = 1e7  # K


def sound_speed(P, rho):
    return np.sqrt(gamma * P / rho)

def sonic_parameter(v_rel, c_s):
    return v_rel / (np.sqrt(2) * c_s)

def G_s(s, phi, z, ln_coulomb):
    return 1.5 * (1.0 + 0.44*s**2) ** 0.5 + z**2 * phi**2 * ln_coulomb * (1.3 + s**3)**(-1)

def Coulomb_potential(Z_gr, a, T):
    return Z_gr*e**2 / (a*KB*T)

c_s = sound_speed(P, rho)
s = sonic_parameter(v_rel, c_s)
phi = Coulomb_potential(Z_gr, a, T)
G = G_s(s, phi, z, ln_coulomb)

print("Sound speed: " + str(c_s*1e-5) + " km/s")
print("Sonic parameter, s: ", s)
print("Coulomb potential parameter: ", phi)
print("G(s): ", G)

def drag_timescale(a, n, T, G):
    return 0.59 * a * n**(-1) * (T/1e6)**(-1/2) * G**(-1)

t_drag = drag_timescale(a, rho/(mu*MP), T, G)

print("Drag timescale: ", str(t_drag), " Myr")