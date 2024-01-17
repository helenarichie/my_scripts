import numpy as np

def shear_time(v_w, r_cl):
    v_w *= 100000
    r_cl *= 3.086e+18
    return r_cl / v_w

print(f"Shear time: {shear_time(1000, 100)/3.154e7/1e3} kyr")