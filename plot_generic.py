import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_generic(datadir, pngdir, ns, ne, cat, ftype, fields, vlims, print_keys=False, floor=None, units="cholla", derive=False, n_xtick=11, n_ytick=11, figsize=None):
    for i in range(ns, ne+1):
        if ftype == "slice":
            if cat:
                f = h5py.File(os.path.join(datadir, str(i)+"_slice.h5"), "r")
            else:
                f = h5py.File(os.path.join(datadir, str(i), str(i)+"_slice.h5.0"), "r")
        if ftype == "proj":
            if cat:
                f = h5py.File(os.path.join(datadir, str(i)+"_proj.h5"), "r")
            else:
                f = h5py.File(os.path.join(datadir, str(i), str(i)+"_proj.h5.0"), "r")
        if ftype == "rot_proj":
            if cat:
                f = h5py.File(os.path.join(datadir, str(i)+"_rot_proj.h5"), "r")
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
        mu = 0.6 # mean molecular weight
        MP = 1.672622e-24 # proton mass, g
        KB = 1.380658e-16 # Boltzmann constant, cm^2 g s^-2 K^-1
        for j, field_name in enumerate(fields):

                if field_name[-2:] == "xy":
                    xlen = nx*dx
                    ylen = ny*dx
                    xlabel, ylabel = "$x~$[kpc]", "$y~$[kpc]"
                if field_name[-2:] == "xz":
                    xlen = nx*dx
                    ylen = nz*dx
                    xlabel, ylabel = "$x~$[kpc]", "$z~$[kpc]"

                if field_name.startswith("d_"):
                    if units == "cgs":
                        if ftype == "slice":
                            conversion = head["density_unit"]
                            clabel = r"$\mathrm{log}_{10}(\rho)$ [$\mathrm{g}\mathrm{cm}^{-3}$]"                           
                    if units == "cholla":
                        if ftype == "slice":
                            clabel = r"$\mathrm{log}_{10}(\rho)$ [$\mathrm{M}_\odot\mathrm{kpc}^{-3}$]"
                        if ftype == "proj":
                            clabel = r'$\mathrm{log}_{10}(\Sigma_{gas})$ [$\mathrm{M}_\odot\,\mathrm{kpc}^{-2}$]'
                if field_name.startswith("T_"):
                    if ftype == "slice":
                        clabel = r"$\mathrm{log}_{10}(T)$ [K]"
                if field_name.startswith("mx_"):
                    if units == "cgs":
                        if ftype == "slice":
                            conversion = head["mass_unit"] * head["velocity_unit"]
                            clabel = r"$\mathrm{log}_{10}(p_x)$ [$g\,cm\,s^{-1}$]"
                    if units == "cholla":
                        if ftype == "slice":
                            clabel = r"$\mathrm{log}_{10}(p_x)$ [$M_\odot\,kpc\,kyr^{-1}$]"
                if field_name.startswith("my_"):
                    if units == "cgs":
                        if ftype == "slice":
                            conversion = head["mass_unit"] * head["velocity_unit"]
                            clabel = r"$\mathrm{log}_{10}(p_y)$ [$g\,cm\,s^{-1}$]"
                    if units == "cholla":
                        if ftype == "slice":
                            clabel = r"$\mathrm{log}_{10}(p_y)$ [$M_\odot\,kpc\,kyr^{-1}$]"
                if field_name.startswith("mz_"):
                    if units == "cgs":
                        if ftype == "slice":
                            conversion = head["mass_unit"] * head["velocity_unit"]
                            clabel = r"$\mathrm{log}_{10}(p_z)$ [$g\,cm\,s^{-1}$]"
                    if units == "cholla":
                        if ftype == "slice":
                            clabel = r"$\mathrm{log}_{10}(p_z)$ [$M_\odot\,kpc\,kyr^{-1}$]"
                if field_name.startswith("scalar"):
                    conversion = 1
                    clabel = ""

                if units == "cholla":
                    conversion = 1

                if not derive:
                    field = np.array(f[field_name]) * conversion
                else:
                    if field_name.startswith("T"):
                        if units == "cgs":
                            density = np.array(f["d_xz"])
                            energy = np.array(f["E_xz"])
                            momentum_x = np.array(f["mx_xz"])
                            momentum_y = np.array(f["my_xz"])
                            momentum_z = np.array(f["mz_xz"])
                            field = ((head["gamma"] - 1) / (((density * head["density_unit"]) / (mu*MP))*KB)) * ((energy * head["energy_unit"]) - (1/2 * (density * head["density_unit"]) * ((momentum_x / density * head["velocity_unit"])**2 + (momentum_y / density * head["velocity_unit"])**2 + (momentum_z / density * head["velocity_unit"])**2)))
                

                if floor != None:
                    field[field<=floor] = floor

                def plot(vlim=False):
                    if figsize == None:
                        fig, ax = plt.subplots()
                    else:
                        fig, ax = plt.subplots(figsize=figsize)
                    if vlim:
                        im = ax.imshow(np.log10(field.T), origin="lower", vmin=vlims[j][0], vmax=vlims[j][1], extent=[0, xlen, 0, ylen])
                    else:
                        im = ax.imshow(np.log10(field.T), origin="lower", extent=[0, xlen, 0, ylen])
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="5%", pad=0.05)
                    cbar = fig.colorbar(im, ax=ax, cax=cax)
                    cbar.set_label(clabel)
                    cbar.ax.tick_params(axis="y", direction="in")
                    try:
                        ax.set_xticks(np.linspace(0, xlen, n_xtick).round(1))
                        ax.set_yticks(np.linspace(0, ylen, n_ytick).round(1))
                    except NameError:
                        print('field_name must end in "xy", "yz", or "xz".')
                    ax.tick_params(axis='both', which='both', direction='in', color='black', top=1, right=1)
                    ax.set_title(f"{field_name}, {t} kyr")
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel(ylabel)
                    fig.tight_layout()

                    return fig
            
                fig = plot(vlim=True)
                fig.savefig(pngdir + f"{i}_{field_name}_{ftype}.png", dpi=300, bbox_inches="tight")
                plt.close()

                fig = plot()
                fig.savefig(pngdir + f"{i}_{field_name}_{ftype}_unscaled.png", dpi=300, bbox_inches="tight")
                plt.close()
                
        f.close()
        print(f"Saving figure {i} of {ne}.\n")
