import numpy as np
import matplotlib.pyplot as plt

kpc_in_cm = 3.086e21
yr_in_s = 3.154e7

############ hard-coded ###########
v_wind = 1000e3  # cm/s
rho_d_init = 1e-26  # g/cm^3
tmax = 1e9 * yr_in_s  # yr to s
h = 1e2 * yr_in_s  # yr to s
x0 = 0.320 * kpc_in_cm # kpc to cm
##################################

bins, r_av, n_av, n_med, n_lo, n_hi, v_av, v_med, v_lo, v_hi, T_av, T_med, T_lo, T_hi, p_av, p_med, p_lo, p_hi, c_av, c_med, c_lo, c_hi, cs_av, cs_med, cs_lo, cs_hi, K_av, K_med, K_lo, K_hi, M_av, M_med, M_lo, M_hi = np.loadtxt('/Users/helenarichie/Documents/Grad School/research/data/CGOLS profiles/2048_central_35_hot_dweight.txt', unpack=True)


# McKinnon et al. (2017)
def calc_tau_sp(n, T):
    YR_IN_S = 3.154e7;
    a1 = 1; # dust grain size in units of 0.1 micrometers
    d0 = n / (6e-4); # gas density in units of 10^-27 g/cm^3
    T_0 = 2e6; # K
    omega = 2.5;
    A = 0.17e9 * YR_IN_S; # Gyr in s

    tau_sp = A * (a1/d0) * ((T_0/T)**omega + 1); # s

    return tau_sp

# McKinnon et al. (2017)
def calc_dd_dt(d_dust, tau_sp):
    return -d_dust / (tau_sp/3)

t_arr = np.arange(0, tmax, h)
d_dust_solution = []
x_distance = []

d_dust = rho_d_init
x_distance_i = x0
for i, t in enumerate(t_arr):
    dt = h
    x_distance_i += v_wind * dt
    x_distance.append(x_distance_i)

    n_i, T_i = n_med[np.argmin(abs(r_av*kpc_in_cm-x_distance_i))], T_med[np.argmin(abs(r_av*kpc_in_cm-x_distance_i))]

    tau_sp = calc_tau_sp(n_i, T_i)
    print(i)
    print(f"t_sp = {round(tau_sp/yr_in_s*1e-6, 4)} Myr")

    dd_dt = calc_dd_dt(d_dust, tau_sp)
    dd = dd_dt * dt

    while dd/d_dust > 0.01:
        dt_sub = 0.01 * d_dust / dd_dt
        d_dust += dt_sub * dd_dt
        dt -= dt_sub
        dd_dt = calc_dd_dt(d_dust, tau_sp)
        dd = dt * dd_dt
    
    d_dust += dd

    d_dust_solution.append(d_dust)

d_dust_solution = np.array(d_dust_solution)
x_distance = np.array(x_distance)
# d_dust_solution /= MP * mu

t_arr /= yr_in_s
x_distance /= kpc_in_cm

plt.rcParams.update({'font.size': 16})

fig = plt.figure()
plt.semilogx(t_arr, d_dust_solution, linewidth=2)
plt.xlabel("Time [yr]")
plt.ylabel(r"$\rho_{dust}~[g\,cm^{-3}]$")
plt.savefig("/Users/helenarichie/Desktop/analytic_hot_phase_time.png")

fig = plt.figure()
plt.plot(x_distance, d_dust_solution, linewidth=2)
plt.xlabel("Distance [kpc]")
plt.ylabel(r"$\rho_{dust}~[g\,cm^{-3}]$")
plt.savefig("/Users/helenarichie/Desktop/analytic_hot_phase_distance.png")

