import sys
sys.path.insert(0, "/Users/helenarichie/GitHub/my_scripts/")
from hconfig import *

# McKinnon et al. (2017)
def calc_tau_sp(n, T):
    YR_IN_S = 3.154e7;
    a1 = .1; # dust grain size in units of 0.1 micrometers
    d0 = n / (6e-4); # gas density in units of 10^-27 g/cm^3
    T_0 = 2e6; # K
    omega = 2.5;
    A = 0.17e9 * YR_IN_S; # Gyr in s

    tau_sp = A * (a1/d0) * ((T_0/T)**omega + 1); # s
    
    return tau_sp / YR_IN_S

# T = 10**np.linspace(1, 8, 8)
# T = [np.sqrt(3e7*3e4)]
# n = [np.sqrt(10*1e-2)]

T = [np.sqrt(3e6*3e3)]
n = [np.sqrt(10*0.01)]

n_cool_base = 8.624272
T_cool_base = 6388.117
n_hot_base = 0.01808641
T_hot_base = 24279940.0

# n = [np.sqrt(n_cool_base * n_hot_base)]
# T = [np.sqrt(T_cool_base * T_hot_base)]
T = [np.sqrt(3e6*3e4)]
n = [np.sqrt(10*0.01)]

# n = [1e-3, 1e-2, 1e-1, 1]
for n_i in n:
    #print(n_i)
    for T_i in T:
        #print("{:e}".format(T_i) + " K:")
        print("{:f}".format(calc_tau_sp(n_i, T_i)*1e-6) + " Myr")
