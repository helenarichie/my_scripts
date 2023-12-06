from hconfig import *

# McKinnon et al. (2017)
def calc_tau_sp(n, T):
    YR_IN_S = 3.154e7;
    a1 = 0.1; # dust grain size in units of 0.1 micrometers
    d0 = n / (6e-4); # gas density in units of 10^-27 g/cm^3
    T_0 = 2e6; # K
    omega = 2.5;
    A = 0.17e9 * YR_IN_S; # Gyr in s

    tau_sp = A * (a1/d0) * ((T_0/T)**omega + 1); # s
    
    return tau_sp / YR_IN_S

# T = 10**np.linspace(1, 8, 8)
T = [3e7]
n = [0.01]

# n = [1e-3, 1e-2, 1e-1, 1]
for n_i in n:
    print(n_i)
    for T_i in T:
        print("{:e}".format(T_i) + " K:")
        print("{:e}".format(calc_tau_sp(n_i, T_i)) + " yr")
