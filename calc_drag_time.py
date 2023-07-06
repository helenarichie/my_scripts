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


# Hu et al. (2019) polynomial coefficients for non-thermal sputtering erosion rates
a0_C, a1_C, a2_C, a3_C, a4_C, a5_C = -48.7, 55.7, -28.8, 7.44, -0.965, 0.005
a0_Si, a1_Si, a2_Si, a3_Si, a4_Si, a5_Si = -32.4, 28.6, -11.6, 2.23, -0.208, 0.00736

v_rel_arr = np.logspace(1, 5, 1000)

yield_C = []
yield_Si = []
for v in np.log10(v_rel_arr):
    yield_C.append(10**np.sum(a0_C + a1_C*v + a2_C*v**2 + a3_C*v**3 + a4_C*v**4 + a5_C*v**5))
    yield_Si.append(10**np.sum(a0_Si + a1_Si*v + a2_Si*v**2 + a3_Si*v**3 + a4_Si*v**4 + a5_Si*v**5))

yield_Si = np.array(yield_Si)

fig, ax = plt.subplots(figsize=(12,8))
plt.loglog(v_rel_arr, yield_C, label="C", color="k", linewidth=3)
plt.loglog(v_rel_arr, yield_Si, linestyle="--", label="Si", color="k", linewidth=3)
plt.xlabel(r"$v_{rel}~[km\,s^{-1}]$")
plt.ylabel(r"$Y~[\mu m\,yr^{-1}\,cm^3]$")
plt.ylim(1e-9, 2e-5)
plt.xlim(1e1, 1e5)
plt.legend()
# ax.tick_params(bottom=True, top=True, left=True, right=True)
plt.savefig("/Users/helenarichie/Desktop/nonth_sputtering_yield.png")

wh_Y = np.argmin(np.abs(v_rel - v_rel_arr))
Y_Si = yield_Si[wh_Y] 

def sputtering_timescale(a, n, Y):
    """Non-thermal sputtering timescale"""
    return 0.33 * a * (n**(-1)) * (Y/1e-6)**(-1)

t_sp_nonth = sputtering_timescale(a, rho/(mu*MP), Y_Si)

print("Nonth. sputtering timescale (Si): ", str(t_sp_nonth), " Myr")




fig, ax = plt.subplots(figsize=(12,8))
plt.semilogx(np.linspace(1e-9, np.amax(yield_Si), 10000), sputtering_timescale(a, rho/(mu*MP), np.linspace(1e-9, np.amax(yield_Si), 10000)), color="k", linewidth=3)
plt.xlabel(r"$Y~[\mu m\,yr^{-1}\,cm^3]$")
plt.ylabel(r"$t_{sp}~[Myr]$")
plt.xlim(1e-9, np.amax(yield_Si))
plt.ylim(-1, 30)
# ax.tick_params(bottom=True, top=True, left=True, right=True)
plt.savefig("/Users/helenarichie/Desktop/nonth_sputtering_time.png")


c_s = sound_speed(P, rho)
s = sonic_parameter(v_rel_arr, c_s)
phi = Coulomb_potential(Z_gr, a, T)
G = G_s(s, phi, z, ln_coulomb)
drags = drag_timescale(a, rho/(mu*MP), T, G)

fig, ax = plt.subplots(figsize=(12,8))
ax.ticklabel_format(axis="y", style="sci")
plt.semilogx(v_rel_arr, drags, color="k", linewidth=3)
plt.xlabel(r"$v_{rel}~[km\,s^{-1}]$")
plt.ylabel(r"$t_{drag}~[Myr]$")
#plt.xlim(1e1, np.amax(yield_Si))
#plt.ylim(np.amin(drags), np.amax(drags))
# ax.tick_params(bottom=True, top=True, left=True, right=True)
plt.savefig("/Users/helenarichie/Desktop/drag_time_v.png")

T_arr = np.logspace(1, 8, 1000)
c_s = sound_speed(P, rho)
s = sonic_parameter(v_rel, c_s)
phi = Coulomb_potential(Z_gr, a, T_arr)
G = G_s(s, phi, z, ln_coulomb)
drags = drag_timescale(a, rho/(mu*MP), T_arr, G)

fig, ax = plt.subplots(figsize=(12,8))
ax.ticklabel_format(axis="y", style="sci")
plt.semilogx(v_rel_arr, drags, color="k", linewidth=3)
plt.xlabel(r"$T~[K]$")
plt.ylabel(r"$t_{drag}~[Myr]$")
#plt.xlim(1e1, np.amax(yield_Si))
#plt.ylim(np.amin(drags), np.amax(drags))
# ax.tick_params(bottom=True, top=True, left=True, right=True)
plt.savefig("/Users/helenarichie/Desktop/drag_time_T.png")

n_arr = np.logspace(-2.5, 1, 1000)
c_s = sound_speed(P, n_arr*mu*MP)
s = sonic_parameter(v_rel, c_s)
phi = Coulomb_potential(Z_gr, a, T_arr)
G = G_s(s, phi, z, ln_coulomb)
drags = drag_timescale(a, n_arr, T_arr, G)

fig, ax = plt.subplots(figsize=(12,8))
ax.ticklabel_format(axis="y", style="sci")
plt.semilogx(n_arr, drags, color="k", linewidth=3)
plt.xlabel(r"$n~[cm^{-3}]$")
plt.ylabel(r"$t_{drag}~[Myr]$")
#plt.xlim(1e1, np.amax(yield_Si))
#plt.ylim(np.amin(drags), np.amax(drags))
# ax.tick_params(bottom=True, top=True, left=True, right=True)
plt.savefig("/Users/helenarichie/Desktop/drag_time_n.png")