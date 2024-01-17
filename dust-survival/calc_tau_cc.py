def calc_tau_cc(chi, r_cl, v_w):
    r_cl *= 3.086e+13  # pc to km
    return chi**(1/2)*r_cl / v_w  # s

print(f"t_cc: {calc_tau_cc(1e2, 5, 1e3)/3.154e+7*1e-6} Myr")