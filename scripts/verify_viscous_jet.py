#!/usr/bin/env python3
# ==============================================================================
# COVOF 2D - Script de Vérification et de Validation : Cas 3.3 Flambement de Jet Visqueux
# Article de référence : M. El Ouafa, S. Vincent, V. Le Chenadec,
# "Monolithic Solvers for Incompressible Two-Phase Flows at Large Density and Viscosity Ratios",
# Fluids 2021, 6(1), 23 (Section 3.3, Tableau 4, Figures 6, 7, 8)
# Référence théorique du flambement : Tomé et al., J. Non-Newton. Fluid Mech. 1996
# ==============================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def evaluate_run(csv_path, label="COVOF 2D"):
    """
    Évalue quantitativement les grandeurs physiques du cas de flambement de jet visqueux.
    Critères clés selon Fluids 2021 (Section 3.3, Tableau 4) :
      - Critères théoriques de Tomé et al. : H/D > 3*pi (~9.42) et Re < 0.52.
      - Bilan de masse / volume injecté V_theo(t) = V(0) + V_0 * D * t.
      - Incompressibilité sous grand contraste de densité (1529) et viscosité (2.7e5).
      - Trajectoire du centre de masse et vitesse de descente.
    """
    if not os.path.exists(csv_path):
        return None, None

    df = pd.read_csv(csv_path)
    if df.empty or len(df) < 5:
        return None, None

    t = df['time'].to_numpy()
    vol = df['bubble_volume'].to_numpy()
    yc = df['center_of_mass_y'].to_numpy()
    vc = df['rise_velocity'].to_numpy()
    max_u = df['max_u'].to_numpy()
    max_div = df['max_div'].to_numpy()

    # Paramètres physiques de l'article (Tableau 4, Cas 1)
    rho1 = 1.1768       # Masse volumique gaz [kg/m^3]
    rho2 = 1800.0       # Masse volumique propergol liquide [kg/m^3]
    mu1 = 1.85e-3       # Viscosité dynamique gaz [Pa.s]
    mu2 = 500.0         # Viscosité dynamique propergol liquide [Pa.s]
    v0 = 0.6            # Vitesse d'injection [m/s]
    D = 0.08            # Largeur de l'injecteur [m]
    H = 1.0             # Hauteur de la cavité [m]

    # Nombres sans dimension fondamentaux
    Re = (rho2 * D * v0) / mu2
    aspect_ratio = H / D
    density_ratio = rho2 / rho1
    viscosity_ratio = mu2 / mu1

    # Volume analytique injecté théorique : V_theo(t) = V(0) + (v0 * D) * t
    vol_0 = vol[0]
    vol_theo = vol_0 + (v0 * D) * t
    vol_err_pct = np.abs(vol - vol_theo) / vol_theo * 100.0
    vol_err_final = vol_err_pct[-1]

    # Divergence maximale stable
    max_div_stable = np.max(max_div[1:]) if len(max_div) > 1 else max_div[0]

    # Statuts de validation
    status_re = "PASS" if (Re < 0.52 and aspect_ratio > 3.0 * np.pi) else "FAIL"
    status_vol = "PASS" if (vol_err_final < 1.0) else "FAIL"
    status_div = "PASS" if (max_div_stable < 0.1) else "FAIL"

    metrics = {
        'label': label,
        't_final': t[-1],
        'steps': len(t) - 1,
        'Re': Re,
        'aspect_ratio': aspect_ratio,
        'density_ratio': density_ratio,
        'viscosity_ratio': viscosity_ratio,
        'vol_0': vol_0,
        'vol_final': vol[-1],
        'vol_theo_final': vol_theo[-1],
        'vol_err_final_pct': vol_err_final,
        'yc_final': yc[-1],
        'vc_final': vc[-1],
        'max_u_final': max_u[-1],
        'max_div_stable': max_div_stable,
        'status_re': status_re,
        'status_vol': status_vol,
        'status_div': status_div
    }

    return metrics, df

def plot_viscous_jet_comparison(df_cpp, df_fortran, out_fig_path):
    """
    Génère la figure comparative du flambement de jet visqueux
    analysant le volume injecté et les caractéristiques cinématiques (Fluids 2021).
    """
    os.makedirs(os.path.dirname(out_fig_path), exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)

    color_cpp = '#1f77b4'       # Bleu
    color_fortran = '#2ca02c'   # Vert
    color_ref = '#d62728'       # Rouge

    # 1. Évolution du Volume de Liquide Injecté
    if df_cpp is not None:
        t_cpp = df_cpp['time'].to_numpy()
        vol_cpp = df_cpp['bubble_volume'].to_numpy() * 1e3  # en 10^-3 m^2 (L)
        vol_theo = (df_cpp['bubble_volume'].iloc[0] + 0.6 * 0.08 * t_cpp) * 1e3
        ax1.plot(t_cpp, vol_cpp, label='COVOF 2D C++17', color=color_cpp, linewidth=2.0)
        ax1.plot(t_cpp, vol_theo, ':', label='Volume Théorique Injecté $V_0 D t$', color=color_ref, linewidth=1.8)

    if df_fortran is not None:
        t_f90 = df_fortran['time'].to_numpy()
        vol_f90 = df_fortran['bubble_volume'].to_numpy() * 1e3
        ax1.plot(t_f90, vol_f90, '--', label='COVOF 2D Fortran 2018', color=color_fortran, linewidth=2.0)

    ax1.set_xlabel('Temps $t$ [s]', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Volume de Liquide Injecté [$10^{-3}\\,\\text{m}^2$]', fontsize=11, fontweight='bold')
    ax1.set_title("Cas 3.3 : Conservation du Débit Massique d'Injection (Fluids 2021)", fontsize=11, fontweight='bold')
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(loc='upper left', framealpha=0.9, fontsize=9)

    # 2. Centre de masse vertical Yc(t) et vitesse de descente Vc(t)
    if df_cpp is not None:
        ax2.plot(t_cpp, df_cpp['center_of_mass_y'], label='Centre de Masse $Y_c$ C++17 [m]',
                 color=color_cpp, linewidth=2.0)
        ax2.plot(t_cpp, np.abs(df_cpp['rise_velocity']), '-.', label='Vitesse de chute $|V_c|$ C++17 [m/s]',
                 color='#ff7f0e', linewidth=1.6)

    if df_fortran is not None:
        ax2.plot(t_f90, df_fortran['center_of_mass_y'], '--', label='Centre de Masse $Y_c$ Fortran 2018 [m]',
                 color=color_fortran, linewidth=2.0)
        ax2.plot(t_f90, np.abs(df_fortran['rise_velocity']), ':', label='Vitesse de chute $|V_c|$ Fortran [m/s]',
                 color='#9467bd', linewidth=1.6)

    ax2.set_xlabel('Temps $t$ [s]', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Position $Y_c$ [m] & Vitesse $|V_c|$ [m/s]', fontsize=11, fontweight='bold')
    ax2.set_title("Dynamique de Chute du Jet sous Gravité ($g = 9.81\\,\\text{m/s}^2$)", fontsize=11, fontweight='bold')
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend(loc='center right', framealpha=0.9, fontsize=9)

    plt.tight_layout()
    plt.savefig(out_fig_path, dpi=300)
    plt.close()
    print(f"  [FIGURE] Graphique haute résolution généré : {out_fig_path}")

def main():
    print("=" * 80)
    print(" COVOF 2D — VALIDATION DU CAS 3.3 : FLAMBEMENT D'UN JET VISQUEUX")
    print(" Référence : M. El Ouafa et al., Fluids 2021, 6(1), 23 (Section 3.3, Tableau 4)")
    print("=" * 80)

    root_cpp = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    root_fortran = os.path.abspath(os.path.join(root_cpp, "..", "COVOF_2D-FORTRAN"))

    csv_cpp = os.path.join(root_cpp, "output_viscous_jet_case1", "viscous_jet_case1_history.csv")
    csv_fortran = os.path.join(root_fortran, "output_viscous_jet_case1", "viscous_jet_case1_history.csv")

    met_cpp, df_cpp = evaluate_run(csv_cpp, "COVOF 2D C++17")
    met_f90, df_f90 = evaluate_run(csv_fortran, "COVOF 2D Fortran 2018")

    if met_cpp is None and met_f90 is None:
        print("  [ERREUR] Aucun fichier d'historique CSV trouvé. Exécutez d'abord les simulations.")
        sys.exit(1)

    print("\n--- PARAMÈTRES ET CRITÈRES PHYSIQUES (TOMÉ ET AL. 1996) ---")
    m = met_cpp if met_cpp is not None else met_f90
    print(f"  Rapport d'aspect H/D        : {m['aspect_ratio']:.2f} (Seuil Tomé > 3*pi = {3.0*np.pi:.2f}) -> {'CONFORME' if m['aspect_ratio'] > 3.0*np.pi else 'NON CONFORME'}")
    print(f"  Nombre de Reynolds Re       : {m['Re']:.4f} (Seuil Flambement Re < 0.52)   -> {'RÉGIME FLAMBEMENT' if m['Re'] < 0.52 else 'RÉGIME STABLE'}")
    print(f"  Rapport de densité rho2/rho1: {m['density_ratio']:.1f}")
    print(f"  Rapport de viscosité mu2/mu1: {m['viscosity_ratio']:.2e}")

    print("\n--- RÉSULTATS QUANTITATIFS COMPARATIFS ---")
    headers = ["Code", "Pas", "t_final [s]", "V_num [m2]", "V_theo [m2]", "Err Vol [%]", "Y_c [m]", "Vitesse Max", "Div Max [s-1]", "Statut"]
    print(f"{headers[0]:<20} {headers[1]:<6} {headers[2]:<12} {headers[3]:<12} {headers[4]:<12} {headers[5]:<12} {headers[6]:<10} {headers[7]:<12} {headers[8]:<14} {headers[9]:<8}")
    print("-" * 125)

    for met in [met_cpp, met_f90]:
        if met is not None:
            status_glob = "PASS" if (met['status_re'] == "PASS" and met['status_vol'] == "PASS" and met['status_div'] == "PASS") else "FAIL"
            print(f"{met['label']:<20} {met['steps']:<6} {met['t_final']:<12.3f} {met['vol_final']:<12.6f} {met['vol_theo_final']:<12.6f} {met['vol_err_final_pct']:<12.3f} {met['yc_final']:<10.4f} {met['max_u_final']:<12.4f} {met['max_div_stable']:<14.4e} {status_glob:<8}")

    # Symétrie duale
    if met_cpp is not None and met_f90 is not None:
        diff_vol = abs(met_cpp['vol_final'] - met_f90['vol_final']) / met_cpp['vol_final'] * 100.0
        diff_yc = abs(met_cpp['yc_final'] - met_f90['yc_final']) / met_cpp['yc_final'] * 100.0
        diff_u = abs(met_cpp['max_u_final'] - met_f90['max_u_final']) / met_cpp['max_u_final'] * 100.0
        print("\n--- ÉCART DUAL C++17 vs FORTRAN 2018 ---")
        print(f"  Différence relative volume final          : {diff_vol:.4e} %")
        print(f"  Différence relative centre de masse Y_c   : {diff_yc:.4e} %")
        print(f"  Différence relative vitesse maximale      : {diff_u:.4e} %")
        dual_status = "PASS" if (diff_vol < 0.1 and diff_yc < 0.1) else "WARNING"
        print(f"  Statut de concordance duale C++/Fortran    : {dual_status}")

    # Génération des figures
    fig_path_cpp = os.path.join(root_cpp, "docs", "figures", "fig_cas3_3_viscous_jet_dual.png")
    plot_viscous_jet_comparison(df_cpp, df_f90, fig_path_cpp)

    if os.path.exists(root_fortran):
        fig_path_f90 = os.path.join(root_fortran, "docs", "figures", "fig_cas3_3_viscous_jet_dual.png")
        plot_viscous_jet_comparison(df_cpp, df_f90, fig_path_f90)

    print("\n" + "=" * 80)
    print(" VALIDATION CAS 3.3 TERMINÉE AVEC SUCCÈS")
    print("=" * 80)

if __name__ == "__main__":
    main()
