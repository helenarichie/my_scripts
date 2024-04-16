import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_generic(datadir, pngdir, ns, ne, cat, ftype, fields, shapes, xlabels, ylabels, clabels, vlims, print_keys=False):

    for i in range(ns, ne+1):
        if ftype == "slice":
            if cat:
                f = h5py.File(os.path.join(datadir, str(i), str(i)+"_slice.h5"), "r")
            else:
                f = h5py.File(os.path.join(datadir, str(i), str(i)+"_slice.h5.0"), "r")
        if ftype == "proj":
            if cat:
                f = h5py.File(os.path.join(datadir, str(i), str(i)+"_proj.h5"), "r")
            else:
                f = h5py.File(os.path.join(datadir, str(i), str(i)+"_proj.h5.0"), "r")
        if ftype == "rot_proj":
            if cat:
                f = h5py.File(os.path.join(datadir, str(i), str(i)+"_rot_proj.h5"), "r")
            else:
                f = h5py.File(os.path.join(datadir, str(i), str(i)+"_rot_proj.h5.0"), "r")
        if print_keys:
            print(list(f.keys()))
        head = f.attrs
        dx = head["dx"][0]
        t = head["t"][0]
        gamma = head["gamma"][0]
        nx, ny, nz = head["dims"][0], head["dims"][1], head["dims"][2]
        dx = head["dx"][0]
        for j, field_name in enumerate(fields):

                if shapes[j] == "xy":
                    xlen = nx*dx
                    ylen = ny*dx

                field = np.array(f[field_name])

                def plot(vlim=False):
                    fig, ax = plt.subplots()
                    if vlim:
                        im = ax.imshow(np.log10(field.T), origin="lower", vmin=vlims[j][0], vmax=vlims[j][1], extent=[0, xlen, 0, ylen])
                    else:
                        im = ax.imshow(np.log10(field.T), origin="lower", extent=[0, xlen, 0, ylen])
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="5%", pad=0.05)
                    cbar = fig.colorbar(im, ax=ax, cax=cax)
                    cbar.set_label(clabels[j])
                    cbar.ax.tick_params(axis="y", direction="in")
                    ax.set_xticks(np.linspace(0, xlen, 11).round(1))
                    ax.set_yticks(np.linspace(0, ylen, 11).round(1))
                    ax.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1)
                    ax.set_title(f"{field_name}")
                    ax.set_xlabel(xlabels[j])
                    ax.set_ylabel(ylabels[j])
                    fig.tight_layout()

                    return fig
            
                fig = plot(vlim=True)
                fig.savefig(pngdir + f"{i}_{field_name}_{ftype}.png", dpi=300, bbox_inches="tight")
            
                fig = plot()
                fig.savefig(pngdir + f"{i}_{field_name}_{ftype}_unscaled.png", dpi=300, bbox_inches="tight")

                plt.close()
                
        f.close()
        print(f"Saving figure {i} of {ne}.\n")