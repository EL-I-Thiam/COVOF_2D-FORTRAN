#!/usr/bin/env python3
# ==============================================================================
# COVOF 2D — Génération des Visualisations des Champs de Vitesse et Tourbillon
# Cas 3.1 (Cylindre dense, Figures 2 & 4) et Cas 3.2 (Nappe liquide, Éqs. 17-18)
# Référence : M. El Ouafa, S. Vincent, V. Le Chenadec, Fluids 2021, 6(1), 23
# ==============================================================================

import os
import re
import glob
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.ndimage import zoom, gaussian_filter

def read_vtr(filepath):
    """Lit un fichier VTK RectilinearGrid (.vtr) et extrait les coordonnées et champs."""
    tree = ET.parse(filepath)
    root = tree.getroot()
    rg = root.find('RectilinearGrid')
    ext = [int(x) for x in rg.attrib['WholeExtent'].split()]
    nx, ny = ext[1] - ext[0], ext[3] - ext[2]
    piece = rg.find('Piece')
    coords = piece.find('Coordinates')
    da_list = coords.findall('DataArray')
    x_edges = np.fromstring(da_list[0].text, sep=' ')
    y_edges = np.fromstring(da_list[1].text, sep=' ')
    xc = 0.5 * (x_edges[:-1] + x_edges[1:])
    yc = 0.5 * (y_edges[:-1] + y_edges[1:])
    celldata = piece.find('CellData')
    fields = {}
    for da in celldata.findall('DataArray'):
        name = da.attrib['Name']
        ncomp = int(da.attrib.get('NumberOfComponents', 1))
        txt = re.sub(r'([0-9])([+-][0-9]{2,3})', r'\1E\2', da.text)
        arr = np.fromstring(txt, sep=' ')
        if ncomp == 1:
            fields[name] = arr.reshape((ny, nx))
        else:
            fields[name] = arr.reshape((ny, nx, ncomp))
    return xc, yc, fields

def refine_field(xc, yc, field, factor=4, sigma=0.6):
    """Raffine un champ 2D par interpolation spline cubique pour un rendu type ParaView."""
    f_smooth = gaussian_filter(field, sigma=sigma) if sigma > 0 else field
    f_fine = zoom(f_smooth, factor, order=3)
    x_fine = np.linspace(xc[0], xc[-1], f_fine.shape[1])
    y_fine = np.linspace(yc[0], yc[-1], f_fine.shape[0])
    return x_fine, y_fine, f_fine

def generate_case3_1_figures(root_cpp, root_f90, out_dirs):
    """
    Reproduit fidèlement la Figure 2 (Champs U, V, |omega| à t = 0.144 s)
    et la Figure 4 (Enstrophie Ens1(t) dans le fluide 1) de Fluids 2021.
    """
    vtr_cpp = os.path.join(root_cpp, "output_free_fall_dense_cylinder", "free_fall_dense_cylinder_p0_002300.vtr")
    vtr_f90 = os.path.join(root_f90, "output_free_fall_dense_cylinder", "free_fall_dense_cylinder_p0_002300.vtr")

    xc_c, yc_c, f_c = read_vtr(vtr_cpp)
    xc_f, yc_f, f_f = read_vtr(vtr_f90)

    dx, dy = xc_c[1] - xc_c[0], yc_c[1] - yc_c[0]

    u_c, v_c = f_c['Velocity'][:, :, 0], f_c['Velocity'][:, :, 1]
    c_c = f_c['PhaseIndicator']
    omega_c = np.abs(np.gradient(v_c, dx, axis=1) - np.gradient(u_c, dy, axis=0))

    u_f, v_f = f_f['Velocity'][:, :, 0], f_f['Velocity'][:, :, 1]
    c_f = f_f['PhaseIndicator']
    omega_f = np.abs(np.gradient(v_f, dx, axis=1) - np.gradient(u_f, dy, axis=0))

    # Centre de masse du cylindre à t = 0.144 s
    X_grid, Y_grid = np.meshgrid(xc_c, yc_c)
    cyl_mask = 1.0 - c_c
    x_cyl = np.sum(X_grid * cyl_mask) / np.sum(cyl_mask)
    y_cyl = np.sum(Y_grid * cyl_mask) / np.sum(cyl_mask)
    r_cyl = 0.0125

    # Raffinement type ParaView pour affichage lisse des champs et isolignes
    xf, yf, u_c_fine = refine_field(xc_c, yc_c, u_c, factor=5, sigma=0.5)
    _, _, v_c_fine = refine_field(xc_c, yc_c, v_c, factor=5, sigma=0.5)
    _, _, w_c_fine = refine_field(xc_c, yc_c, omega_c, factor=5, sigma=0.7)

    _, _, u_f_fine = refine_field(xc_f, yc_f, u_f, factor=5, sigma=0.5)
    _, _, v_f_fine = refine_field(xc_f, yc_f, v_f, factor=5, sigma=0.5)
    _, _, w_f_fine = refine_field(xc_f, yc_f, omega_f, factor=5, sigma=0.7)

    XF, YF = np.meshgrid(xf, yf)
    mask_disk = ((XF - x_cyl)**2 + (YF - y_cyl)**2) <= (r_cyl * 1.02)**2

    u_c_plot = np.ma.masked_where(mask_disk, u_c_fine)
    v_c_plot = np.ma.masked_where(mask_disk, v_c_fine)
    w_c_plot = np.ma.masked_where(mask_disk, w_c_fine)

    u_f_plot = np.ma.masked_where(mask_disk, u_f_fine)
    v_f_plot = np.ma.masked_where(mask_disk, v_f_fine)
    w_f_plot = np.ma.masked_where(mask_disk, w_f_fine)

    # ==========================================================================
    # FIGURE 2 (Fluids 2021, Page 9) : Champs de vitesse U, V et Tourbillon |omega|
    # ==========================================================================
    fig, axes = plt.subplots(2, 3, figsize=(15, 11.5), dpi=300)
    plt.subplots_adjust(wspace=0.28, hspace=0.30)

    u_levels = np.linspace(-1.0, 1.0, 21)
    v_levels = np.linspace(-1.4, 1.2, 21)
    w_levels = np.linspace(10.0, 420.0, 18)

    configs = [
        (axes[0, 0], u_c_plot, u_f_plot, [-1.0, 1.0], u_levels, [0.0, 0.1], [0.0, 0.11],
         "Vitesse horizontale $U$ [m/s]\nC++17 (champ) & Fortran (isolignes)", "U [m/s]"),
        (axes[0, 1], v_c_plot, v_f_plot, [-1.4, 1.2], v_levels, [0.0, 0.1], [0.0, 0.10],
         "Vitesse verticale $V$ [m/s]\nC++17 (champ) & Fortran (isolignes)", "V [m/s]"),
        (axes[0, 2], w_c_plot, w_f_plot, [0.0, 420.0], w_levels, [0.02, 0.08], [0.030, 0.095],
         "Magnitude du Tourbillon $|\\omega|$ [$\\text{s}^{-1}$]\nC++17 (champ) & Fortran (isolignes)", "|ω| [s⁻¹]"),
        (axes[1, 0], u_f_plot, u_c_plot, [-1.0, 1.0], u_levels, [0.0, 0.1], [0.0, 0.11],
         "Vitesse horizontale $U$ [m/s]\nFortran 2018 (champ) & C++17 (isolignes)", "U [m/s]"),
        (axes[1, 1], v_f_plot, v_c_plot, [-1.4, 1.2], v_levels, [0.0, 0.1], [0.0, 0.10],
         "Vitesse verticale $V$ [m/s]\nFortran 2018 (champ) & C++17 (isolignes)", "V [m/s]"),
        (axes[1, 2], w_f_plot, w_c_plot, [0.0, 420.0], w_levels, [0.02, 0.08], [0.030, 0.095],
         "Magnitude du Tourbillon $|\\omega|$ [$\\text{s}^{-1}$]\nFortran 2018 (champ) & C++17 (isolignes)", "|ω| [s⁻¹]"),
    ]

    for ax, field_fill, field_iso, vlim, iso_lev, xlim, ylim, title, cbar_label in configs:
        im = ax.pcolormesh(XF, YF, field_fill, cmap='jet', vmin=vlim[0], vmax=vlim[1], shading='auto')
        ax.contour(XF, YF, field_iso, levels=iso_lev, colors='k', linewidths=0.65, linestyles='solid', alpha=0.85)
        circ = Circle((x_cyl, y_cyl), r_cyl, facecolor='white', edgecolor='black', linewidth=1.5, zorder=5)
        ax.add_patch(circ)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect('equal')
        ax.set_xlabel('$x$ [m]', fontsize=9.5, fontweight='bold')
        ax.set_ylabel('$y$ [m]', fontsize=9.5, fontweight='bold')
        ax.set_title(title, fontsize=10, fontweight='bold', pad=8)
        cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.14, fraction=0.048)
        cbar.set_label(cbar_label, fontsize=8.5, fontweight='bold')
        cbar.ax.tick_params(labelsize=8)

    fig.suptitle("Cas 3.1 (Fluids 2021, Fig. 2) : Chute libre d'un cylindre dense à $t = 0.144\\,\\text{s}$ ($\\rho_2/\\rho_1 = 10^5, \\mu_2/\\mu_1 = 10^{10}$)",
                 fontsize=13, fontweight='bold', y=0.97)

    for d in out_dirs:
        os.makedirs(d, exist_ok=True)
        out_p = os.path.join(d, "fig_cas3_1_velocity_vorticity_fields.png")
        fig.savefig(out_p, dpi=300, bbox_inches='tight')
    plt.close(fig)

    # ==========================================================================
    # FIGURE 4 (Fluids 2021, Page 10) + Norme de Vitesse & Lignes de Courant
    # ==========================================================================
    steps = range(100, 2400, 100)
    t_list, ens_cpp, ens_f90 = [], [], []
    for st in steps:
        fc_path = os.path.join(root_cpp, "output_free_fall_dense_cylinder", f"free_fall_dense_cylinder_p0_{st:06d}.vtr")
        ff_path = os.path.join(root_f90, "output_free_fall_dense_cylinder", f"free_fall_dense_cylinder_p0_{st:06d}.vtr")
        if os.path.exists(fc_path) and os.path.exists(ff_path):
            _, _, dc = read_vtr(fc_path)
            _, _, df = read_vtr(ff_path)
            uc_s, vc_s = dc['Velocity'][:, :, 0], dc['Velocity'][:, :, 1]
            cc_s = dc['PhaseIndicator']
            wc_s = np.gradient(vc_s, dx, axis=1) - np.gradient(uc_s, dy, axis=0)
            uf_s, vf_s = df['Velocity'][:, :, 0], df['Velocity'][:, :, 1]
            cf_s = df['PhaseIndicator']
            wf_s = np.gradient(vf_s, dx, axis=1) - np.gradient(uf_s, dy, axis=0)
            t_list.append(st * 6.25e-5)
            ens_cpp.append(np.sum(cc_s * 0.5 * wc_s**2) * dx * dy)
            ens_f90.append(np.sum(cf_s * 0.5 * wf_s**2) * dx * dy)

    fig2, (ax_s1, ax_s2, ax_ens) = plt.subplots(1, 3, figsize=(16, 5.4), dpi=300)
    plt.subplots_adjust(wspace=0.42)

    umag_c = np.sqrt(u_c_fine**2 + v_c_fine**2)
    umag_f = np.sqrt(u_f_fine**2 + v_f_fine**2)

    im1 = ax_s1.pcolormesh(XF, YF, umag_c, cmap='turbo', vmin=0.0, vmax=1.5, shading='auto')
    ax_s1.streamplot(xf, yf, u_c_fine, v_c_fine, color='white', linewidth=0.8, density=1.3, arrowsize=0.9)
    ax_s1.add_patch(Circle((x_cyl, y_cyl), r_cyl, facecolor='none', edgecolor='black', linewidth=2.0, zorder=6))
    ax_s1.set_xlim([0.0, 0.1])
    ax_s1.set_ylim([0.0, 0.2])
    ax_s1.set_aspect('equal')
    ax_s1.set_title("Norme $\\|\\mathbf{u}\\|$ & Lignes de courant\nCOVOF 2D C++17 ($t=0.144\\,\\text{s}$)", fontsize=10.5, fontweight='bold')
    ax_s1.set_xlabel('$x$ [m]', fontweight='bold')
    ax_s1.set_ylabel('$y$ [m]', fontweight='bold')
    plt.colorbar(im1, ax=ax_s1, label='$\\|\\mathbf{u}\\|$ [m/s]', fraction=0.046, pad=0.04)

    im2 = ax_s2.pcolormesh(XF, YF, umag_f, cmap='turbo', vmin=0.0, vmax=1.5, shading='auto')
    ax_s2.streamplot(xf, yf, u_f_fine, v_f_fine, color='white', linewidth=0.8, density=1.3, arrowsize=0.9)
    ax_s2.add_patch(Circle((x_cyl, y_cyl), r_cyl, facecolor='none', edgecolor='black', linewidth=2.0, zorder=6))
    ax_s2.set_xlim([0.0, 0.1])
    ax_s2.set_ylim([0.0, 0.2])
    ax_s2.set_aspect('equal')
    ax_s2.set_title("Norme $\\|\\mathbf{u}\\|$ & Lignes de courant\nCOVOF 2D Fortran 2018 ($t=0.144\\,\\text{s}$)", fontsize=10.5, fontweight='bold')
    ax_s2.set_xlabel('$x$ [m]', fontweight='bold')
    ax_s2.set_ylabel('$y$ [m]', fontweight='bold')
    plt.colorbar(im2, ax=ax_s2, label='$\\|\\mathbf{u}\\|$ [m/s]', fraction=0.046, pad=0.04)

    ax_ens.loglog(t_list, ens_cpp, 's-', color='#d62728', label='COVOF 2D C++17 (AL)', markersize=5, linewidth=1.8)
    ax_ens.loglog(t_list, ens_f90, 'o--', color='#1f77b4', label='COVOF 2D Fortran 2018 (AL)', markersize=5, fillstyle='none', linewidth=1.8)
    ax_ens.set_xlabel('Temps $t$ [s] (échelle log)', fontsize=10.5, fontweight='bold')
    ax_ens.set_ylabel('Enstrophie fluide 1 $\\text{Ens}_1$ [$\\text{m}^3\\cdot\\text{s}^{-2}$]', fontsize=10.5, fontweight='bold')
    ax_ens.set_title("Cas 3.1 (Fluids 2021, Fig. 4) :\nÉvolution temporelle de l'Enstrophie $\\text{Ens}_1(t)$", fontsize=10.5, fontweight='bold')
    ax_ens.grid(True, which='both', linestyle='--', alpha=0.6)
    ax_ens.legend(loc='lower right', fontsize=9.5)

    for d in out_dirs:
        out_p2 = os.path.join(d, "fig_cas3_1_streamlines_enstrophy.png")
        fig2.savefig(out_p2, dpi=300, bbox_inches='tight')
    plt.close(fig2)

def generate_case3_2_figures(root_cpp, root_f90, out_dirs):
    """
    Génère les visualisations complètes des champs de vitesse (U, V à t=0 et t=tf, |omega|, profils)
    pour le Cas 3.2 : Atomisation d'une nappe liquide périodique (Éqs. 17-18, Fluids 2021).
    """
    vtr_cpp_0 = os.path.join(root_cpp, "output_liquid_sheet_ratio100", "liquid_sheet_ratio100_p0_000000.vtr")
    vtr_f90_0 = os.path.join(root_f90, "output_liquid_sheet_ratio100", "liquid_sheet_ratio100_p0_000000.vtr")
    vtr_cpp_f = sorted(glob.glob(os.path.join(root_cpp, "output_liquid_sheet_ratio100", "liquid_sheet_ratio100_p0_*.vtr")))[-1]
    vtr_f90_f = sorted(glob.glob(os.path.join(root_f90, "output_liquid_sheet_ratio100", "liquid_sheet_ratio100_p0_*.vtr")))[-1]

    xc_c, yc_c, f_c0 = read_vtr(vtr_cpp_0)
    _, _, f_f0 = read_vtr(vtr_f90_0)
    _, _, f_cf = read_vtr(vtr_cpp_f)
    _, _, f_ff = read_vtr(vtr_f90_f)

    u_c0, v_c0 = f_c0['Velocity'][:, :, 0], f_c0['Velocity'][:, :, 1]
    u_f0, v_f0 = f_f0['Velocity'][:, :, 0], f_f0['Velocity'][:, :, 1]
    c_c0 = f_c0['PhaseIndicator']

    u_cf, v_cf = f_cf['Velocity'][:, :, 0], f_cf['Velocity'][:, :, 1]
    u_ff, v_ff = f_ff['Velocity'][:, :, 0], f_ff['Velocity'][:, :, 1]
    c_cf = f_cf['PhaseIndicator']
    c_ff = f_ff['PhaseIndicator']

    x_mm = np.linspace(xc_c[0] * 1e3, xc_c[-1] * 1e3, len(xc_c))
    y_mm = np.linspace(yc_c[0] * 1e3, yc_c[-1] * 1e3, len(yc_c))
    X_mm, Y_mm = np.meshgrid(x_mm, y_mm)

    fig, axes = plt.subplots(2, 3, figsize=(16.5, 10.2), dpi=300)
    plt.subplots_adjust(wspace=0.34, hspace=0.32)

    # 1. Vitesse longitudinale U(x,y) C++17 à t = t_f + Lignes de courant
    ax = axes[0, 0]
    im = ax.pcolormesh(X_mm, Y_mm, u_cf, cmap='jet', vmin=0.0, vmax=31.0, shading='auto')
    ax.streamplot(x_mm, y_mm, u_cf, v_cf, color='white', linewidth=0.7, density=1.1, arrowsize=0.8)
    ax.contour(X_mm, Y_mm, c_cf, levels=[0.5], colors='white', linewidths=2.2)
    ax.contour(X_mm, Y_mm, c_ff, levels=[0.5], colors='black', linewidths=1.1, linestyles='--')
    ax.set_xlim([0.0, 3.0])
    ax.set_ylim([0.0, 3.0])
    ax.set_title("Vitesse longitudinale $u(x,y)$ [m/s]\nCOVOF 2D C++17 & Lignes de courant ($t=t_f$)", fontsize=10, fontweight='bold')
    ax.set_xlabel('$x$ [mm]', fontweight='bold')
    ax.set_ylabel('$y$ [mm]', fontweight='bold')
    ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='$u$ [m/s]', fraction=0.046, pad=0.04)

    # 2. Vitesse transversale initiale V(x,y, t=0) (Équation 17) : C++17 (champ) + Fortran (isolignes)
    ax = axes[0, 1]
    im = ax.pcolormesh(X_mm, Y_mm, v_c0 * 1e2, cmap='coolwarm', vmin=-4.0, vmax=4.0, shading='auto')
    v_iso = np.linspace(-3.6, 3.6, 13)
    ax.contour(X_mm, Y_mm, v_f0 * 1e2, levels=v_iso, colors='k', linewidths=0.75, alpha=0.8)
    ax.contour(X_mm, Y_mm, c_c0, levels=[0.5], colors='yellow', linewidths=1.8, linestyles='--')
    ax.set_title("Perturbation solénoïdale $v(x,y,t=0)$ [cm/s]\nC++17 (champ) & Fortran (isolignes, Éq. 17)", fontsize=10, fontweight='bold')
    ax.set_xlabel('$x$ [mm]', fontweight='bold')
    ax.set_ylabel('$y$ [mm]', fontweight='bold')
    ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='$v(t=0)$ [cm/s]', fraction=0.046, pad=0.04)

    # 3. Tourbillon de cisaillement interfacial |omega(x,y)| C++17
    dx, dy = xc_c[1] - xc_c[0], yc_c[1] - yc_c[0]
    omega_c = np.abs(np.gradient(v_cf, dx, axis=1) - np.gradient(u_cf, dy, axis=0)) * 1e-5
    ax = axes[0, 2]
    im = ax.pcolormesh(X_mm, Y_mm, omega_c, cmap='inferno', vmin=0.0, vmax=np.max(omega_c)*0.95, shading='auto')
    ax.contour(X_mm, Y_mm, c_cf, levels=[0.5], colors='cyan', linewidths=1.5, linestyles='--')
    ax.set_title("Tourbillon interfacial $|\\omega(x,y)|$ [$10^5\\,\\text{s}^{-1}$]\nCouches de cisaillement gaz-liquide (C++17)", fontsize=10, fontweight='bold')
    ax.set_xlabel('$x$ [mm]', fontweight='bold')
    ax.set_ylabel('$y$ [mm]', fontweight='bold')
    ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='$|\\omega|$ [$10^5$ s⁻¹]', fraction=0.046, pad=0.04)

    # 4. Vitesse longitudinale U(x,y) Fortran 2018 à t = t_f + Lignes de courant
    ax = axes[1, 0]
    im = ax.pcolormesh(X_mm, Y_mm, u_ff, cmap='jet', vmin=0.0, vmax=31.0, shading='auto')
    ax.streamplot(x_mm, y_mm, u_ff, v_ff, color='white', linewidth=0.7, density=1.1, arrowsize=0.8)
    ax.contour(X_mm, Y_mm, c_ff, levels=[0.5], colors='white', linewidths=2.2)
    ax.contour(X_mm, Y_mm, c_cf, levels=[0.5], colors='black', linewidths=1.1, linestyles='--')
    ax.set_xlim([0.0, 3.0])
    ax.set_ylim([0.0, 3.0])
    ax.set_title("Vitesse longitudinale $u(x,y)$ [m/s]\nCOVOF 2D Fortran 2018 & Lignes de courant ($t=t_f$)", fontsize=10, fontweight='bold')
    ax.set_xlabel('$x$ [mm]', fontweight='bold')
    ax.set_ylabel('$y$ [mm]', fontweight='bold')
    ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='$u$ [m/s]', fraction=0.046, pad=0.04)

    # 5. Vitesse transversale V(x,y, t=0) Fortran 2018 (champ) + C++17 (isolignes)
    ax = axes[1, 1]
    im = ax.pcolormesh(X_mm, Y_mm, v_f0 * 1e2, cmap='coolwarm', vmin=-4.0, vmax=4.0, shading='auto')
    ax.contour(X_mm, Y_mm, v_c0 * 1e2, levels=v_iso, colors='k', linewidths=0.75, alpha=0.8)
    ax.contour(X_mm, Y_mm, c_c0, levels=[0.5], colors='yellow', linewidths=1.8, linestyles='--')
    ax.set_title("Perturbation solénoïdale $v(x,y,t=0)$ [cm/s]\nFortran 2018 (champ) & C++17 (isolignes, Éq. 17)", fontsize=10, fontweight='bold')
    ax.set_xlabel('$x$ [mm]', fontweight='bold')
    ax.set_ylabel('$y$ [mm]', fontweight='bold')
    ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='$v(t=0)$ [cm/s]', fraction=0.046, pad=0.04)

    # 6. Profils analytiques vs numériques U(y) et V(x, yc) (Équations 17-18 de Fluids 2021)
    ax = axes[1, 2]
    mid_x = len(xc_c) // 2
    mid_y = len(yc_c) // 2
    ax.plot(y_mm, u_c0[:, mid_x], ':', color='#555555', linewidth=1.8, label='$u(y, t=0)$ Exact (Éq. 18)')
    ax.plot(y_mm, u_cf[:, mid_x], '-', color='#1f77b4', linewidth=2.2, label='$u(y, t_f)$ C++17 [m/s]')
    ax.plot(y_mm, u_ff[:, mid_x], '--', color='#2ca02c', linewidth=2.0, label='$u(y, t_f)$ Fortran 2018 [m/s]')
    ax2 = ax.twinx()
    L_val = 3.0e-3
    v_exact = 0.04 * np.sin(2.0 * np.pi * xc_c / L_val) * 1e2
    ax2.plot(x_mm, v_exact, '-', color='#d62728', linewidth=2.0, alpha=0.8, label='$v(x, y_c, t=0)$ Exact [cm/s]')
    ax2.plot(x_mm, v_c0[mid_y, :] * 1e2, 'o', color='#ff7f0e', markersize=3.5, markevery=3, label='$v(x, y_c, t=0)$ C++17')
    ax2.plot(x_mm, v_f0[mid_y, :] * 1e2, 'x', color='#9467bd', markersize=4.0, markevery=3, label='$v(x, y_c, t=0)$ Fortran')
    ax.set_xlabel('Coordonnée spatiale $y$ ou $x$ [mm]', fontweight='bold')
    ax.set_ylabel('Vitesse longitudinale $u$ [m/s]', color='#1f77b4', fontweight='bold')
    ax2.set_ylabel('Perturbation initiale $v(t=0)$ [cm/s]', color='#d62728', fontweight='bold')
    ax2.set_ylim([-5.5, 5.5])
    ax.set_title("Profils de Vitesse $u(y)$ et $v(x, y_c)$ (Éqs. 17–18)\nSaut de cisaillement & Mode sinusoïdal", fontsize=10, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='lower left', fontsize=7.5, framealpha=0.92)

    fig.suptitle(r"Cas 3.2 (Fluids 2021, Éqs. 17–18 & Fig. 5) : Champs de Vitesse de la Nappe Liquide sous Cisaillement ($\rho_2/\rho_1 = 100$)",
                 fontsize=13, fontweight='bold', y=0.97)

    for d in out_dirs:
        out_p = os.path.join(d, "fig_cas3_2_velocity_fields.png")
        fig.savefig(out_p, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("  [OK] Figures Cas 3.1 et Cas 3.2 régénérées avec succès.")

if __name__ == "__main__":
    root_cpp = "/Users/elhadjiibrahimathiam/Documents/PROJET_ANTIGRAVITY/COVOF_2D-C++"
    root_f90 = "/Users/elhadjiibrahimathiam/Documents/PROJET_ANTIGRAVITY/COVOF_2D-FORTRAN"
    brain_figs = "/Users/elhadjiibrahimathiam/.gemini/antigravity/brain/d891e39d-85a0-4365-9d86-f527fe9e206c/figures"
    out_dirs = [
        os.path.join(root_cpp, "docs", "figures"),
        os.path.join(root_f90, "docs", "figures"),
        brain_figs
    ]
    generate_case3_1_figures(root_cpp, root_f90, out_dirs)
    generate_case3_2_figures(root_cpp, root_f90, out_dirs)
