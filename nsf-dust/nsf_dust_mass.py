import sys
sys.path.insert(0, "../")
from hconfig import *
from matplotlib.ticker import FormatStrFormatter

plt.rcParams.update({'font.size': 35})

t_exit = 4.8e7

dust = []
with open("../csv/dust_tot.csv") as f:
    for line in f:
        line = line.split(",")
        dust.append(float(line[0]))
dust = np.array(dust) * 256

gas = []
with open("../csv/cloud_tot.csv") as f:
    for line in f:
        line = line.split(",")
        gas.append(float(line[0]))
gas = np.array(gas)

t_arr = []
with open("../csv/t_arr.csv") as f:
    for line in f:
        line = line.split(",")
        t_arr.append(float(line[0]))
t_arr = np.array(t_arr)

others = []
with open("../csv/others.csv") as f:
    for line in f:
        line = line.split(",")
        others.append(float(line[0]))

gas_i = gas[0]
dust_i = dust[0]

dust_sub = dust[t_arr<=t_exit]
gas_sub = gas[t_arr<=t_exit]

t_arr = t_arr[t_arr<=t_exit]

# gas
fig = plt.figure(figsize=(12,10))

def set_size(w, h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)

plt.plot(t_arr/1e6, gas_sub/gas_i, linewidth=5, color="purple")

plt.ylabel(r'$m_{cloud}/m_{cloud,i}$', fontsize=40)
plt.xlabel(r'Time [Myr]', fontsize=40)

plt.ylim(0.93, np.amax(gas_sub/gas_i))
plt.xlim(1, 5e1)

plt.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=10)
plt.yticks(np.linspace(1, 2.15, 5))
plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#plt.xticks(np.log10([1e1, 1e2, 1e3, 1e5, 1e5, 1e6, 1e7]))
#plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.e'))
#plt.tight_layout()
plt.savefig("../png/gas_mass.png", dpi=300, bbox_inches="tight")

# dust
fig = plt.figure(figsize=(12, 10))

plt.plot(t_arr/1e6, dust_sub/dust_i, linewidth=5, color="blue")

plt.ylabel(r'$m_{dust}/m_{dust,i}$', fontsize=40)
plt.xlabel(r'Time [Myr]', fontsize=40)

plt.xlim(1, 5e1)
plt.ylim(0.9941, 1.0004)

plt.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1, length=10)
plt.yticks(np.linspace(0.9945, 1.000, 5))
plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

#plt.tight_layout()
plt.savefig("../png/dust_mass.png", dpi=300, bbox_inches="tight")