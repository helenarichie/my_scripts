import sys
sys.path.insert(0, "/ix/eschneider/helena/code/my_scripts")
from hconfig import *
import seaborn as sns

density_conversion = 5.028e-34/(3.24e-22)**3 # g/cm^3 to M_sun/kpc^3

##################################################################
date = "frontier/2024-02-07"
spacing= 400*1e-3 # spacing of tick marks in units and sets axes labels and units of dx (kpc or pc)
# spacing = 40
cat = True
pad = 0.1
fontsize = 28
labelpad = 12
tickwidth = 2
ticklength = 13
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams.update({'font.size': fontsize})
fnum = 700
##################################################################

##################################################################
basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/"
pngdir = os.path.join(basedir, "png/")
slicedir = os.path.join(basedir, "hdf5/slice/")
##################################################################ss

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 5))

data = ReadHDF5(slicedir, nscalar=1, fnum=fnum, slice="xy", cat=cat)
head = data.head
conserved = data.conserved
dx = head["dx"][0]
nx, ny, nz = head["dims"]

left_offset = 900
right_offset = 2800
xlen = nx-left_offset-right_offset

d_gas = data.d_cgs()[0][left_offset:(nx-right_offset),:]
d_dust = conserved["scalar0"][0][left_offset:(nx-right_offset),:] * head["density_unit"]
t_arr = data.t_cgs() / yr_in_s

im = ax.imshow((d_dust/d_gas).T, origin="lower", extent=[0, xlen*dx, 0, ny*dx], cmap=sns.color_palette("cubehelix", as_cmap=True), vmin=0, vmax=0.01)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="4%", pad=0.08)
cbar = plt.colorbar(im, cax=cax)
cbar.set_label("dust-to-gas ratio", fontsize=fontsize, rotation=270, labelpad=40)
cbar.ax.tick_params(length=ticklength, width=tickwidth, color="white", labelsize=fontsize-5)
ax.tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1, length=ticklength, width=tickwidth)
ax.hlines(0.75*dx*nz, spacing, 2*spacing, color='white', linewidth=2)
ax.text(2.2*spacing, 0.725*dx*nz, '400 pc', color='white', fontsize=fontsize)
ax.set_xticks(np.arange(0, xlen*dx, spacing))
ax.set_yticks(np.arange(0, ny*dx, spacing))
ax.text(xlen*dx-2.5*spacing, 0.725*dx*nz, f'{round(t_arr[0]/1e6, 1)} Myr', color='white', fontsize=fontsize)

plt.tight_layout()
plt.savefig(pngdir + f"dtg_slice.png", dpi=300, bbox_inches="tight")