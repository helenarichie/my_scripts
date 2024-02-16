from hconfig import *

fontsize = 20
linewidth = 4
pad = 0.1
mu = 0.6

plt.rcParams.update({'font.size': fontsize})

a_cool, n_cool, T_cool, = 0.1, 10, 3e3
p_cool = n_cool * T_cool * KB

a_mild, n_mild, T_mild, = 0.1, 1e-2, 3e6
p_mild = n_mild * T_mild * KB

a_hot, n_hot, T_hot, = 0.1, 1e-2, 3e7
p_hot = n_hot * T_hot * KB

a_mix, n_mix, T_mix, = 0.1, np.sqrt(n_cool*n_hot), np.sqrt(T_cool*T_hot)
p_mix = n_mix * T_mix * KB

cool = [a_cool, n_cool, T_cool, p_cool]
mild = [a_mild, n_mild, T_mild, p_mild]
mix = [a_mix, n_mix, T_mix, p_mix]
hot = [a_hot, n_hot, T_hot, p_hot]

def G(v, p, rho):
    s = v*100000 / (np.sqrt(2)*np.sqrt(5/2*(p/rho)))
    return s

def drag_timescale(a, n, T, v, p):
    return 0.59 * a * n**(-1) * (T/1e6)**(-1/2) / G(v, p, n*mu*MP)

def nonth_timescale(a, n, Y):
    return 0.33 * a * n**(-1) * (Y/1e-6)**(-1)

def nonth_yield_C(v):
    return 0.005*v**5 - 0.965*v**4 + 7.44*v**3 - 28.8*v**2 + 55.7*v - 48.7

def nonth_yield_Si(v):
    return 0.00736*v**5 - 0.208*v**4 + 2.23*v**3 - 11.6*v**2 + 28.6*v - 32.4

v_range = np.logspace(1, 3.5, 1000)

t_surv, vx_surv, vy_surv, vz_surv = [], [], [], []
survdir = "/ix/eschneider/helena/data/cloud_wind/2024-02-06/csv/"
with open(os.path.join(survdir, "avg_v.csv"), "r") as f:
    for line in f:
        line = line.split(",")
        t_surv.append(float(line[0]))
        vx_surv.append(float(line[1]))
        vy_surv.append(float(line[2]))
        vz_surv.append(float(line[3]))
t_surv, vx_surv, vy_surv, vz_surv = np.array(t_surv), np.array(vx_surv), np.array(vy_surv), np.array(vz_surv)

t_disr, vx_disr, vy_disr, vz_disr = [], [], [], []
disrdir = "/ix/eschneider/helena/data/cloud_wind/2024-02-07/csv/"
with open(os.path.join(disrdir, "avg_v.csv"), "r") as f:
    for line in f:
        line = line.split(",")
        t_disr.append(float(line[0]))
        vx_disr.append(float(line[1]))
        vy_disr.append(float(line[2]))
        vz_disr.append(float(line[3]))
t_disr, vx_disr, vy_disr, vz_disr = np.array(t_disr), np.array(vx_disr), np.array(vy_disr), np.array(vz_disr)

plt.plot(t_surv*1e-3, vx_surv*1e-5, linewidth=3, c=colors[0])

v_range_yield = np.logspace(1, 5, 1000)
plt.figure()
plt.loglog(v_range_yield, 10**nonth_yield_C(np.log10(v_range_yield)), c="k",label="C", linewidth=linewidth)
plt.loglog(v_range_yield, 10**nonth_yield_Si(np.log10(v_range_yield)), c="k", linestyle="--", label="Si", linewidth=linewidth)
plt.ylim(1e-9, 3e-5)
plt.xlim(1e1, 1e5)
plt.xticks(np.logspace(1, 5, 5))
plt.xlabel("Time [Myr]")
plt.ylabel(r"$v_{cl,x}$ [km/s]")
plt.legend(loc="upper left")
plt.savefig("yield.png", bbox_inches="tight")

plt.figure(figsize=(6, 5.5))
plt.loglog(v_range, drag_timescale(cool[0], cool[1], cool[2], v_range, cool[3]), linestyle="--", c=colors[5], linewidth=linewidth-0.5, zorder=2, alpha=0.7)
plt.loglog(v_range, nonth_timescale(cool[0], cool[1], 10**nonth_yield_Si(np.log10(v_range))), c=colors[5], linewidth=linewidth, zorder=2)
#plt.loglog(v_range, drag_timescale(mild[0], mild[1], mild[2], v_range, mild[3]), linestyle="--", color="#fa9fb5", linewidth=linewidth)
#plt.loglog(v_range, nonth_timescale(mild[0], mild[1], 10**nonth_yield_Si(np.log10(v_range))), c="#fa9fb5", linewidth=linewidth)

plt.loglog(v_range, drag_timescale(hot[0], hot[1], hot[2], v_range, hot[3]), linestyle="--", c="palevioletred", linewidth=linewidth-0.5, zorder=0, alpha=0.7)
plt.loglog(v_range, nonth_timescale(hot[0], hot[1], 10**nonth_yield_Si(np.log10(v_range))), c="palevioletred", linewidth=linewidth, zorder=0)

plt.loglog(v_range, drag_timescale(mix[0], mix[1], hot[2], v_range, mix[3]), c=colors[7], linestyle="--", linewidth=linewidth-0.5, zorder=1, alpha=0.7)
plt.loglog(v_range, nonth_timescale(mix[0], mix[1], 10**nonth_yield_Si(np.log10(v_range))), c=colors[7], linewidth=linewidth, zorder=1)

plt.plot(0, 0, c="k", linestyle="--", label=r"$t_{drag}$", linewidth=linewidth)
plt.plot(0, 0, c="k", label=r"$t_{sp,nonth}$", linewidth=linewidth)

plt.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True)
plt.ylabel("Time [Myr]")
plt.xlabel(r"$v_{rel}$ [km/s]")
plt.xlim(np.amin(v_range) - pad, np.amax(v_range) + pad)
#plt.xticks(np.logspace(np.log10(np.amin(v_range)), np.log10(np.amax(v_range)), 5))
plt.legend(fontsize=fontsize-5)
plt.tight_layout()
plt.savefig("drag_nonth.png", bbox_inches="tight")

plt.figure(figsize=(6, 5.5))

plt.plot(t_surv*1e-3, vx_surv*1e-5, linewidth=linewidth, c="steelblue", label="survived")
plt.plot(t_disr*1e-3, vx_disr*1e-5, linewidth=linewidth, c="mediumaquamarine", label="disrupted")

plt.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=9, width=2, reset=True)
plt.xlabel("Time [Myr]")
plt.ylabel(r"$v_{cl,x}$ [km/s]")
plt.xlim(np.amin(t_surv*1e-3) - pad, 56 + pad)
plt.xticks(np.linspace(np.amin(t_surv*1e-3), 56, 5))
plt.legend(fontsize=fontsize-5)
plt.tight_layout()
plt.savefig("v_cl.png", bbox_inches="tight")
