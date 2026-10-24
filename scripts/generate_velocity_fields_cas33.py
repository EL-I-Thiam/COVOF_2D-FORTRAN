#!/usr/bin/env python3
# ==============================================================================
# COVOF 2D — Visualisations des Champs de Vitesse, Lignes de Courant et Régimes
# Cas 3.3 : Flambement d'un Jet Visqueux (Figures 6, 7, 8 & Champs U, V, ||u||)
# Référence : M. El Ouafa, S. Vincent, V. Le Chenadec, Fluids 2021, 6(1), 23
# ==============================================================================

import os
import re
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
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
    xc = np.linspace(0.5 * (x_edges[0] + x_edges[1]), 0.5 * (x_edges[-2] + x_edges[-1]), nx)
    yc = np.linspace(0.5 * (y_edges[0] + y_edges[1]), 0.5 * (y_edges[-2] + y_edges[-1]), ny)
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

def refine_field(xc, yc, field, factor=4, sigma=0.5):
    """Raffine un champ 2D par spline cubique pour un rendu lisse type ParaView."""
    f_smooth = gaussian_filter(field, sigma=sigma) if sigma > 0 else field
    f_fine = zoom(f_smooth, factor, order=3)
    x_fine = np.linspace(xc[0], xc[-1], f_fine.shape[1])
    y_fine = np.linspace(yc[0], yc[-1], f_fine.shape[0])
    return x_fine, y_fine, f_fine

def generate_figure6_three_regimes(root_cpp, root_f90, out_dirs):
    """
    Reproduit la Figure 6 (Page 13 de Fluids 2021) :
    Évolution de l'interface du jet visqueux pour les 3 nombres de Reynolds du Tableau 4 :
      - Colonne 1 (Rouge)  : Cas 1 (Re = 0.1728, H/D = 12.5) à t = 0.25, 0.75, 1.25, 2.0 s
      - Colonne 2 (Bleu)   : Cas 2 (Re = 0.60,   H/D = 10.0) à t = 0.25, 0.75, 1.25, 2.0 s
      - Colonne 3 (Olive)  : Cas 3 (Re = 3.60,   H/D = 10.0) à t = 0.20, 0.30, 0.60, 0.9 s
    """
    cases_info = [
        ("output_viscous_jet_full",  "viscous_jet_full",  [250, 750, 1250, 2000], [0.25, 0.75, 1.25, 2.00], "#f01818", "Cas 1 : $Re = 0.1728$ ($H/D=12.5$)"),
        ("output_viscous_jet_case2", "viscous_jet_case2", [250, 750, 1250, 2000], [0.25, 0.75, 1.25, 2.00], "#0018f0", "Cas 2 : $Re = 0.60$ ($H/D=10.0$)"),
        ("output_viscous_jet_case3", "viscous_jet_case3", [800, 1200, 2400, 3600], [0.20, 0.30, 0.60, 0.90], "#889614", "Cas 3 : $Re = 3.60$ ($H/D=10.0$)"),
    ]

    fig, axes = plt.subplots(4, 3, figsize=(12.5, 15.5), dpi=300)
    plt.subplots_adjust(wspace=0.18, hspace=0.18)

    for col, (folder, prefix, steps, times, color_hex, col_title) in enumerate(cases_info):
        for row, (st, t_val) in enumerate(zip(steps, times)):
            ax = axes[row, col]
            f_cpp = os.path.join(root_cpp, folder, f"{prefix}_p0_{st:06d}.vtr")
            f_f90 = os.path.join(root_f90, folder, f"{prefix}_p0_{st:06d}.vtr")

            xc, yc, dc = read_vtr(f_cpp)
            _, _, df = read_vtr(f_f90)

            xf, yf, c_c_fine = refine_field(xc, yc, dc['PhaseIndicator'], factor=4, sigma=0.45)
            _, _, c_f_fine = refine_field(xc, yc, df['PhaseIndicator'], factor=4, sigma=0.45)
            XF, YF = np.meshgrid(xf, yf)

            # Liquide : C < 0.5
            ax.contourf(XF, YF, c_c_fine, levels=[-1.0, 0.5], colors=[color_hex], alpha=0.92)
            ax.contour(XF, YF, c_c_fine, levels=[0.5], colors=[color_hex], linewidths=1.2)
            ax.contour(XF, YF, c_f_fine, levels=[0.5], colors=['black'], linewidths=1.0, linestyles='--')

            ax.set_xlim([0.0, 1.0])
            ax.set_ylim([0.0, 1.0])
            ax.set_aspect('equal')
            if row == 0:
                ax.set_title(f"{col_title}\n$t = {t_val:.2f}\\,\\text{{s}}$", fontsize=10.5, fontweight='bold')
            else:
                ax.set_title(f"$t = {t_val:.2f}\\,\\text{{s}}$", fontsize=10, fontweight='bold')

            if col == 0:
                ax.set_ylabel('$y$ [m]', fontsize=9.5, fontweight='bold')
            else:
                ax.set_yticklabels([])
            if row == 3:
                ax.set_xlabel('$x$ [m]', fontsize=9.5, fontweight='bold')
            else:
                ax.set_xticklabels([])

    fig.suptitle(r"Cas 3.3 (Fluids 2021, Fig. 6 & Fig. 8) : Évolution du Jet Visqueux ($\rho_2/\rho_1 = 1529, \mu_2/\mu_1 = 2.7\times 10^5$)"
                 "\nRemplissage couleur : COVOF 2D C++17 | Contour noir tireté : COVOF 2D Fortran 2018",
                 fontsize=12.5, fontweight='bold', y=0.96)

    for d in out_dirs:
        out_p = os.path.join(d, "fig_cas3_3_figure6_three_regimes.png")
        fig.savefig(out_p, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("  [OK] Figure 6 (3 régimes de Reynolds) sauvegardée.")

def generate_figure7_streamlines_and_velocity(root_cpp, root_f90, out_dirs):
    """
    1) Reproduit la Figure 7 (Page 14 de Fluids 2021) :
       Lignes de courant et recirculations diphasiques autour du jet (Re = 0.1728)
       aux instants t = 0.25 s, 0.75 s, 1.25 s et 2.00 s (comparaison C++17 vs Fortran 2018).
    2) Génère la planche complète des champs de vitesse (||u||, V, U + isolignes) aux 4 instants.
    """
    steps = [250, 750, 1250, 2000]
    times = [0.25, 0.75, 1.25, 2.00]

    # ==========================================================================
    # FIGURE 7 (Page 14 de Fluids 2021) : Lignes de courant + Interface rouge
    # ==========================================================================
    fig7, axes7 = plt.subplots(2, 4, figsize=(17, 8.8), dpi=300)
    plt.subplots_adjust(wspace=0.18, hspace=0.25)

    for col, (st, t_val) in enumerate(zip(steps, times)):
        f_cpp = os.path.join(root_cpp, "output_viscous_jet_full", f"viscous_jet_full_p0_{st:06d}.vtr")
        f_f90 = os.path.join(root_f90, "output_viscous_jet_full", f"viscous_jet_full_p0_{st:06d}.vtr")
        xc, yc, dc = read_vtr(f_cpp)
        _, _, df = read_vtr(f_f90)

        xf, yf, c_c_fine = refine_field(xc, yc, dc['PhaseIndicator'], factor=4, sigma=0.45)
        _, _, c_f_fine = refine_field(xc, yc, df['PhaseIndicator'], factor=4, sigma=0.45)
        XF, YF = np.meshgrid(xf, yf)

        uc, vc = dc['Velocity'][:, :, 0], dc['Velocity'][:, :, 1]
        uf, vf = df['Velocity'][:, :, 0], df['Velocity'][:, :, 1]

        # Ligne 0 : C++17
        ax_c = axes7[0, col]
        ax_c.contourf(XF, YF, c_c_fine, levels=[-1.0, 0.5], colors=['#f01818'], alpha=0.92)
        ax_c.contour(XF, YF, c_c_fine, levels=[0.5], colors=['#b00000'], linewidths=1.2)
        ax_c.streamplot(xc, yc, uc, vc, color='#222222', linewidth=0.75, density=1.35, arrowsize=0.85)
        ax_c.set_xlim([0.0, 1.0])
        ax_c.set_ylim([0.0, 1.0])
        ax_c.set_aspect('equal')
        ax_c.set_title(f"COVOF 2D C++17 — $t = {t_val:.2f}\\,\\text{{s}}$", fontsize=10.5, fontweight='bold')
        if col == 0:
            ax_c.set_ylabel('$y$ [m]', fontweight='bold')
        ax_c.set_xlabel('$x$ [m]', fontweight='bold')

        # Ligne 1 : Fortran 2018
        ax_f = axes7[1, col]
        ax_f.contourf(XF, YF, c_f_fine, levels=[-1.0, 0.5], colors=['#f01818'], alpha=0.92)
        ax_f.contour(XF, YF, c_f_fine, levels=[0.5], colors=['#b00000'], linewidths=1.2)
        ax_f.streamplot(xc, yc, uf, vf, color='#222222', linewidth=0.75, density=1.35, arrowsize=0.85)
        ax_f.set_xlim([0.0, 1.0])
        ax_f.set_ylim([0.0, 1.0])
        ax_f.set_aspect('equal')
        ax_f.set_title(f"COVOF 2D Fortran 2018 — $t = {t_val:.2f}\\,\\text{{s}}$", fontsize=10.5, fontweight='bold')
        if col == 0:
            ax_f.set_ylabel('$y$ [m]', fontweight='bold')
        ax_f.set_xlabel('$x$ [m]', fontweight='bold')

    fig7.suptitle(r"Cas 3.3 (Fluids 2021, Fig. 7) : Lignes de Courant et Recirculations Tourbillonnaires ($Re = 0.1728, H/D = 12.5$)"
                  "\nInteraction Jet de Propergol (rouge) / Tourbillons de Gaz Ambiant aux instants $t = 0.25, 0.75, 1.25, 2.00\\,\\text{s}$",
                  fontsize=12.5, fontweight='bold', y=0.98)

    for d in out_dirs:
        out_p7 = os.path.join(d, "fig_cas3_3_figure7_streamlines.png")
        fig7.savefig(out_p7, dpi=300, bbox_inches='tight')
    plt.close(fig7)
    print("  [OK] Figure 7 (Lignes de courant et recirculations) sauvegardée.")

    # ==========================================================================
    # PLANCHE COMPLÈTE DES CHAMPS DE VITESSE (||u||, V, U) DU CAS 3.3
    # ==========================================================================
    fig_v, axes_v = plt.subplots(3, 4, figsize=(18.5, 13.8), dpi=300)
    plt.subplots_adjust(wspace=0.35, hspace=0.32)

    for col, (st, t_val) in enumerate(zip(steps, times)):
        f_cpp = os.path.join(root_cpp, "output_viscous_jet_full", f"viscous_jet_full_p0_{st:06d}.vtr")
        f_f90 = os.path.join(root_f90, "output_viscous_jet_full", f"viscous_jet_full_p0_{st:06d}.vtr")
        xc, yc, dc = read_vtr(f_cpp)
        _, _, df = read_vtr(f_f90)

        uc, vc = dc['Velocity'][:, :, 0], dc['Velocity'][:, :, 1]
        uf, vf = df['Velocity'][:, :, 0], df['Velocity'][:, :, 1]
        cc, cf = dc['PhaseIndicator'], df['PhaseIndicator']

        xf, yf, uc_fine = refine_field(xc, yc, uc, factor=4, sigma=0.45)
        _, _, vc_fine = refine_field(xc, yc, vc, factor=4, sigma=0.45)
        _, _, uf_fine = refine_field(xc, yc, uf, factor=4, sigma=0.45)
        _, _, vf_fine = refine_field(xc, yc, vf, factor=4, sigma=0.45)
        _, _, cc_fine = refine_field(xc, yc, cc, factor=4, sigma=0.45)
        _, _, cf_fine = refine_field(xc, yc, cf, factor=4, sigma=0.45)
        XF, YF = np.meshgrid(xf, yf)

        umag_c = np.sqrt(uc_fine**2 + vc_fine**2)

        # Ligne 0 : Norme de vitesse ||u|| + Streamlines + Interface C++17/Fortran
        ax0 = axes_v[0, col]
        im0 = ax0.pcolormesh(XF, YF, umag_c, cmap='turbo', vmin=0.0, vmax=1.25, shading='auto')
        ax0.streamplot(xc, yc, uc, vc, color='white', linewidth=0.65, density=1.1, arrowsize=0.75)
        ax0.contour(XF, YF, cc_fine, levels=[0.5], colors='white', linewidths=1.8)
        ax0.contour(XF, YF, cf_fine, levels=[0.5], colors='black', linewidths=1.0, linestyles='--')
        ax0.set_xlim([0.0, 1.0])
        ax0.set_ylim([0.0, 1.0])
        ax0.set_aspect('equal')
        ax0.set_title(f"Norme $\\|\\mathbf{{u}}\\|$ [m/s] — $t = {t_val:.2f}\\,\\text{{s}}$", fontsize=10, fontweight='bold')
        if col == 0:
            ax0.set_ylabel('Norme $\\|\\mathbf{u}\\|$\n$y$ [m]', fontweight='bold')
        plt.colorbar(im0, ax=ax0, fraction=0.046, pad=0.04)

        # Ligne 1 : Vitesse verticale V(x,y) C++17 (champ) + Fortran (isolignes)
        ax1 = axes_v[1, col]
        im1 = ax1.pcolormesh(XF, YF, vc_fine, cmap='jet', vmin=-1.2, vmax=0.3, shading='auto')
        v_lev = np.linspace(-1.1, 0.2, 14)
        ax1.contour(XF, YF, vf_fine, levels=v_lev, colors='k', linewidths=0.55, alpha=0.75)
        ax1.contour(XF, YF, cc_fine, levels=[0.5], colors='white', linewidths=1.6)
        ax1.set_xlim([0.0, 1.0])
        ax1.set_ylim([0.0, 1.0])
        ax1.set_aspect('equal')
        ax1.set_title(f"Vitesse verticale $V$ [m/s] — $t = {t_val:.2f}\\,\\text{{s}}$", fontsize=10, fontweight='bold')
        if col == 0:
            ax1.set_ylabel('Vitesse $V$ (C++/F90)\n$y$ [m]', fontweight='bold')
        plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)

        # Ligne 2 : Vitesse horizontale U(x,y) Fortran 2018 (champ) + C++17 (isolignes)
        ax2 = axes_v[2, col]
        im2 = ax2.pcolormesh(XF, YF, uf_fine, cmap='coolwarm', vmin=-0.35, vmax=0.35, shading='auto')
        u_lev = np.linspace(-0.30, 0.30, 13)
        ax2.contour(XF, YF, uc_fine, levels=u_lev, colors='k', linewidths=0.55, alpha=0.75)
        ax2.contour(XF, YF, cf_fine, levels=[0.5], colors='black', linewidths=1.6)
        ax2.set_xlim([0.0, 1.0])
        ax2.set_ylim([0.0, 1.0])
        ax2.set_aspect('equal')
        ax2.set_title(f"Vitesse horizontale $U$ [m/s] — $t = {t_val:.2f}\\,\\text{{s}}$", fontsize=10, fontweight='bold')
        if col == 0:
            ax2.set_ylabel('Vitesse $U$ (F90/C++)\n$y$ [m]', fontweight='bold')
        ax2.set_xlabel('$x$ [m]', fontweight='bold')
        plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)

    fig_v.suptitle(r"Cas 3.3 (Fluids 2021) : Champs de Vitesse Complets ($\|\mathbf{u}\|, V, U$) et Isolignes C++17 / Fortran 2018 ($Re = 0.1728$)",
                   fontsize=13, fontweight='bold', y=0.97)

    for d in out_dirs:
        out_pv = os.path.join(d, "fig_cas3_3_velocity_fields.png")
        fig_v.savefig(out_pv, dpi=300, bbox_inches='tight')
    plt.close(fig_v)
    print("  [OK] Planche complète des champs de vitesse Cas 3.3 sauvegardée.")

if __name__ == "__main__":
    root_cpp = "/Users/elhadjiibrahimathiam/Documents/PROJET_ANTIGRAVITY/COVOF_2D-C++"
    root_f90 = "/Users/elhadjiibrahimathiam/Documents/PROJET_ANTIGRAVITY/COVOF_2D-FORTRAN"
    brain_figs = "/Users/elhadjiibrahimathiam/.gemini/antigravity/brain/d891e39d-85a0-4365-9d86-f527fe9e206c/figures"
    out_dirs = [
        os.path.join(root_cpp, "docs", "figures"),
        os.path.join(root_f90, "docs", "figures"),
        brain_figs
    ]
    generate_figure6_three_regimes(root_cpp, root_f90, out_dirs)
    generate_figure7_streamlines_and_velocity(root_cpp, root_f90, out_dirs)
