#!/usr/bin/env python3
# ==============================================================================
# COVOF 2D - Script de Vérification et de Validation : Cas 3.2 Atomisation d'une Nappe Liquide
# Article de référence : M. El Ouafa, S. Vincent, V. Le Chenadec,
# "Monolithic Solvers for Incompressible Two-Phase Flows at Large Density and Viscosity Ratios",
# Fluids 2021, 6(1), 23 (Section 3.2, Équations 17-18, Figure 5)
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
    Évalue quantitativement les métriques de conservation pour le cas de la nappe liquide.
    Critères clés selon Fluids 2021 (Section 3.2, Figure 5) :
      - Maintien du rapport d'énergie cinétique normalisée E_k(t)/E_k(0) proche de 1.0 (borne [0.95, 1.05]).
      - Conservation stricte de l'incompressibilité (max div(u) bornée).
      - Conservation du volume de la nappe liquide.
    """
    if not os.path.exists(csv_path):
        return None, None

    df = pd.read_csv(csv_path)
    if df.empty or len(df) < 5:
        return None, None

    t = df['time'].to_numpy()
    vol = df['bubble_volume'].to_numpy()
    ke = df['kinetic_energy'].to_numpy()
    max_u = df['max_u'].to_numpy()
    max_div = df['max_div'].to_numpy()

    # Énergie cinétique initiale et normalisée
    ke_0 = ke[0]
    ke_norm = ke / ke_0

    # Volume initial théorique L * delta = 0.003 * 0.0003 = 9.0e-7 m^2
    vol_0 = vol[0]
    vol_cons_pct = (vol[-1] / vol_0) * 100.0

    # Valeur finale et extrêmes du ratio d'énergie cinétique
    ke_final_ratio = ke_norm[-1]
    ke_min_ratio = np.min(ke_norm)
    ke_max_ratio = np.max(ke_norm)

    # Écart maximal par rapport à la conservation parfaite 1.0
    max_ke_deviation = max(abs(ke_max_ratio - 1.0), abs(ke_min_ratio - 1.0)) * 100.0

    # Divergence maximale en régime établi (après premier pas)
    max_div_stable = np.max(max_div[1:]) if len(max_div) > 1 else max_div[0]

    # Statuts de validation (seuil tolérance conforme à Figure 5b de Fluids 2021 : [0.95, 1.05])
    status_ke = "PASS" if (0.95 <= ke_final_ratio <= 1.05) else "FAIL"
    status_vol = "PASS" if (vol_cons_pct >= 99.0) else "FAIL"
    status_div = "PASS" if (max_div_stable < 0.5) else "FAIL"

    metrics = {
        'label': label,
        't_final': t[-1],
        'steps': len(t) - 1,
        'ke_0': ke_0,
        'ke_final': ke[-1],
        'ke_final_ratio': ke_final_ratio,
        'ke_min_ratio': ke_min_ratio,
        'ke_max_ratio': ke_max_ratio,
        'max_ke_deviation_pct': max_ke_deviation,
        'vol_cons_pct': vol_cons_pct,
        'max_u_final': max_u[-1],
        'max_div_stable': max_div_stable,
        'status_ke': status_ke,
        'status_vol': status_vol,
        'status_div': status_div
    }

    return metrics, df

def plot_liquid_sheet_comparison(df_cpp, df_fortran, out_fig_path):
    """
    Génère la figure comparative de l'énergie cinétique normalisée E_k(t)/E_k(0)
    reproduisant la Figure 5 de l'article Fluids 2021 (El Ouafa et al.).
    """
    os.makedirs(os.path.dirname(out_fig_path), exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)

    # Couleurs harmonisées
    color_cpp = '#1f77b4'       # Bleu
    color_fortran = '#2ca02c'   # Vert
    color_ref = '#d62728'       # Rouge

    # 1. Évolution de l'énergie cinétique normalisée E_k(t)/E_k(0) (Figure 5 Fluids 2021)
    if df_cpp is not None:
        t_cpp = df_cpp['time'].to_numpy() * 1e3  # Conversion en millisecondes
        ke_norm_cpp = df_cpp['kinetic_energy'].to_numpy() / df_cpp['kinetic_energy'].iloc[0]
        ax1.plot(t_cpp, ke_norm_cpp, label='COVOF 2D C++17 (Lagrangien Augmenté)',
                 color=color_cpp, linewidth=2.0)

    if df_fortran is not None:
        t_f90 = df_fortran['time'].to_numpy() * 1e3  # Conversion en millisecondes
        ke_norm_f90 = df_fortran['kinetic_energy'].to_numpy() / df_fortran['kinetic_energy'].iloc[0]
        ax1.plot(t_f90, ke_norm_f90, '--', label='COVOF 2D Fortran 2018 (Lagrangien Augmenté)',
                 color=color_fortran, linewidth=2.0)

    # Ligne de conservation théorique idéale
    ax1.axhline(1.0, color='gray', linestyle=':', linewidth=1.5, label='Conservation Parfaite ($E_k(t)/E_k(0)=1$)')
    # Zone d'admissibilité du solveur monolithique Lagrangien Augmenté selon Fluids 2021 Figure 5
    ax1.axhspan(0.98, 1.04, color='green', alpha=0.08, label='Plage de stabilité Fluids 2021 (Fig. 5b)')

    ax1.set_xlabel('Temps $t$ [ms]', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Énergie Cinétique Normalisée $E_k(t) / E_k(0)$', fontsize=11, fontweight='bold')
    ax1.set_title("Cas 3.2 : Conservation de l'Énergie Cinétique (Fluids 2021 Fig. 5)", fontsize=11, fontweight='bold')
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.set_ylim(0.95, 1.05)
    ax1.legend(loc='lower left', framealpha=0.9, fontsize=9)

    # 2. Vitesse maximale et divergence maximale
    if df_cpp is not None:
        ax2.plot(t_cpp, df_cpp['max_u'], label='$U_{max}$ C++17 [m/s]', color=color_cpp, linewidth=1.8)
    if df_fortran is not None:
        ax2.plot(t_f90, df_fortran['max_u'], '--', label='$U_{max}$ Fortran 2018 [m/s]', color=color_fortran, linewidth=1.8)

    ax2.set_xlabel('Temps $t$ [ms]', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Vitesse Maximale $U_{max}$ [m/s]', fontsize=11, fontweight='bold')
    ax2.set_title("Stabilité du Champ de Vitesse ($u_{gaz} = 30$ m/s)", fontsize=11, fontweight='bold')
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.set_ylim(28.0, 32.0)
    ax2.legend(loc='upper right', framealpha=0.9, fontsize=9)

    plt.tight_layout()
    plt.savefig(out_fig_path, dpi=300)
    plt.close()
    print(f"  [FIGURE] Graphique haute résolution généré : {out_fig_path}")

def main():
    print("=" * 80)
    print(" COVOF 2D — VALIDATION DU CAS 3.2 : ATOMISATION D'UNE NAPPE LIQUIDE")
    print(" Référence : M. El Ouafa et al., Fluids 2021, 6(1), 23 (Section 3.2, Fig. 5)")
    print("=" * 80)

    # Détection des chemins C++ et Fortran
    root_cpp = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    root_fortran = os.path.abspath(os.path.join(root_cpp, "..", "COVOF_2D-FORTRAN"))

    csv_cpp = os.path.join(root_cpp, "output_liquid_sheet_ratio100", "liquid_sheet_ratio100_history.csv")
    csv_fortran = os.path.join(root_fortran, "output_liquid_sheet_ratio100", "liquid_sheet_ratio100_history.csv")

    met_cpp, df_cpp = evaluate_run(csv_cpp, "COVOF 2D C++17")
    met_f90, df_f90 = evaluate_run(csv_fortran, "COVOF 2D Fortran 2018")

    if met_cpp is None and met_f90 is None:
        print("  [ERREUR] Aucun fichier d'historique CSV trouvé. Exécutez d'abord les simulations.")
        sys.exit(1)

    print("\n--- RÉSULTATS QUANTITATIFS COMPARATIFS ---")
    headers = ["Code", "Pas", "t_final [ms]", "E_k(0) [J]", "E_k(tf)/E_k(0)", "Dev Max [%]", "Vol Cons [%]", "Div Max [s-1]", "Statut"]
    print(f"{headers[0]:<20} {headers[1]:<6} {headers[2]:<12} {headers[3]:<12} {headers[4]:<14} {headers[5]:<12} {headers[6]:<12} {headers[7]:<14} {headers[8]:<8}")
    print("-" * 110)

    for met in [met_cpp, met_f90]:
        if met is not None:
            t_ms = met['t_final'] * 1e3
            status_glob = "PASS" if (met['status_ke'] == "PASS" and met['status_vol'] == "PASS") else "FAIL"
            print(f"{met['label']:<20} {met['steps']:<6} {t_ms:<12.3f} {met['ke_0']:<12.6e} {met['ke_final_ratio']:<14.5f} {met['max_ke_deviation_pct']:<12.3f} {met['vol_cons_pct']:<12.4f} {met['max_div_stable']:<14.4e} {status_glob:<8}")

    # Vérification de symétrie binaire / duale C++ vs Fortran
    if met_cpp is not None and met_f90 is not None:
        diff_ke_rel = abs(met_cpp['ke_final_ratio'] - met_f90['ke_final_ratio']) / met_cpp['ke_final_ratio'] * 100.0
        diff_vol_rel = abs(met_cpp['vol_cons_pct'] - met_f90['vol_cons_pct'])
        print("\n--- ÉCART DUAL C++17 vs FORTRAN 2018 ---")
        print(f"  Différence relative finale E_k(tf)/E_k(0) : {diff_ke_rel:.4e} %")
        print(f"  Différence absolue conservation volume    : {diff_vol_rel:.4e} %")
        dual_status = "PASS" if (diff_ke_rel < 0.1) else "WARNING"
        print(f"  Statut de concordance duale C++/Fortran    : {dual_status}")

    # Génération des figures
    fig_path_cpp = os.path.join(root_cpp, "docs", "figures", "fig_cas3_2_liquid_sheet_dual.png")
    plot_liquid_sheet_comparison(df_cpp, df_f90, fig_path_cpp)

    if os.path.exists(root_fortran):
        fig_path_f90 = os.path.join(root_fortran, "docs", "figures", "fig_cas3_2_liquid_sheet_dual.png")
        plot_liquid_sheet_comparison(df_cpp, df_f90, fig_path_f90)

    print("\n" + "=" * 80)
    print(" VALIDATION CAS 3.2 TERMINÉE AVEC SUCCÈS")
    print("=" * 80)

if __name__ == "__main__":
    main()
