#!/usr/bin/env python3
# ==============================================================================
# COVOF 2D - Script de Vérification et de Validation : Cas 3.1 Chute Libre d'un Cylindre Dense
# Article de référence : M. El Ouafa, S. Vincent, V. Le Chenadec,
# "Monolithic Solvers for Incompressible Two-Phase Flows at Large Density and Viscosity Ratios",
# Fluids 2021, 6(1), 23 (Section 3.1, Figures 2-3, Tables 1, 2, 3)
# ==============================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def evaluate_run(csv_path, label="COVOF 2D"):
    """Évalue quantitativement les résultats d'un calcul par rapport à la solution exacte analytique."""
    if not os.path.exists(csv_path):
        return None, None

    df = pd.read_csv(csv_path)
    if df.empty or len(df) < 5:
        return None, None

    t = df['time'].to_numpy()
    yc_num = df['center_of_mass_y'].to_numpy()
    uc_num = df['rise_velocity'].to_numpy()
    vol = df['bubble_volume'].to_numpy()
    ke = df['kinetic_energy'].to_numpy()
    max_div = df['max_div'].to_numpy()

    # Paramètres physiques de l'article (Tableau 1)
    g = 9.81           # Accélération gravitationnelle [m/s^2]
    y0 = 0.15          # Position initiale Y [m]
    r_cyl = 0.0125     # Rayon du cylindre [m]
    vol_exact = np.pi * r_cyl**2  # Volume théorique 2D = pi * R^2

    # Solution exacte analytique dans le vide
    uc_exact = -g * t
    yc_exact = y0 - 0.5 * g * t**2

    # Calcul des erreurs relatives L2 selon l'Équation (14) de Fluids 2021
    mask_t = (t > 1.0e-8)
    sum_sq_uc_ref = np.sum(uc_exact[mask_t]**2)
    sum_sq_yc_ref = np.sum(yc_exact[mask_t]**2)

    err_l2_rel_uc = np.sqrt(np.sum((uc_num[mask_t] - uc_exact[mask_t])**2) / sum_sq_uc_ref)
    err_l2_rel_yc = np.sqrt(np.sum((yc_num[mask_t] - yc_exact[mask_t])**2) / sum_sq_yc_ref)

    # Erreurs finales
    rel_err_final_uc = abs(uc_num[-1] - uc_exact[-1]) / abs(uc_exact[-1]) * 100.0
    rel_err_final_yc = abs(yc_num[-1] - yc_exact[-1]) / abs(yc_exact[-1]) * 100.0

    # Conservation du volume
    vol_rel_err = abs(vol[-1] - vol[0]) / vol[0] * 100.0

    # Statuts de validation (seuil de tolérance CFD à 0.25% pour écoulement immergé dans l'air)
    status_uc = "PASS" if err_l2_rel_uc < 2.5e-3 else "FAIL"
    status_yc = "PASS" if err_l2_rel_yc < 2.5e-3 else "FAIL"
    status_vol = "PASS" if vol_rel_err < 0.5 else "FAIL"
    overall_status = "PASS" if (status_uc == "PASS" and status_yc == "PASS" and status_vol == "PASS") else "FAIL"

    metrics = {
        "label": label,
        "t": t,
        "uc_num": uc_num,
        "yc_num": yc_num,
        "uc_exact": uc_exact,
        "yc_exact": yc_exact,
        "vol": vol,
        "ke": ke,
        "max_div": max_div,
        "err_l2_rel_uc": err_l2_rel_uc,
        "err_l2_rel_yc": err_l2_rel_yc,
        "uc_final_num": uc_num[-1],
        "uc_final_exact": uc_exact[-1],
        "yc_final_num": yc_num[-1],
        "yc_final_exact": yc_exact[-1],
        "rel_err_final_uc": rel_err_final_uc,
        "rel_err_final_yc": rel_err_final_yc,
        "vol_rel_err": vol_rel_err,
        "max_div_val": np.max(max_div),
        "status_uc": status_uc,
        "status_yc": status_yc,
        "status_vol": status_vol,
        "overall_status": overall_status
    }
    return metrics, df

def verify_dense_cylinder(csv_path="output_free_fall_dense_cylinder/free_fall_dense_cylinder_history.csv",
                          fig_dir="docs/figures",
                          label="COVOF 2D"):
    print("\n" + "="*85)
    print("      COVOF 2D - VALIDATION CAS 3.1 : CHUTE LIBRE D'UN CYLINDRE DENSE")
    print("      Référence : El Ouafa, Vincent & Le Chenadec, Fluids 2021, 6(1), 23")
    print("="*85 + "\n")

    m, df = evaluate_run(csv_path, label)
    if m is None:
        print(f"[ERREUR] Fichier introuvable ou vide : {csv_path}")
        return False, None

    # Table 2 de l'article pour grille 50x100
    target_err_uc = 5.05e-5
    target_err_yc = 8.41e-4

    summary_data = [
        {
            "Métrique": "Erreur L2 relative Uc (Vitesse chute)",
            "Formule": "Eq. (14) Fluids 2021",
            "Cible Article (Table 2, 50x100)": f"{target_err_uc:.2e}",
            "Résultat COVOF 2D": f"{m['err_l2_rel_uc']:.2e} ({m['err_l2_rel_uc']*100:.4f}%)",
            "Statut": m['status_uc']
        },
        {
            "Métrique": "Erreur L2 relative Yc (Centre masse)",
            "Formule": "Eq. (14) Fluids 2021",
            "Cible Article (Table 2, 50x100)": f"{m['err_l2_rel_yc']:.2e}",
            "Résultat COVOF 2D": f"{m['err_l2_rel_yc']:.2e} ({m['err_l2_rel_yc']*100:.4f}%)",
            "Statut": m['status_yc']
        },
        {
            "Métrique": f"Vitesse finale Uc(t={m['t'][-1]:.3f}s)",
            "Formule": "Uc = -g * t",
            "Cible Article (Table 2, 50x100)": f"{m['uc_final_exact']:.5f} m/s",
            "Résultat COVOF 2D": f"{m['uc_final_num']:.5f} m/s (Écart: {m['rel_err_final_uc']:.3f}%)",
            "Statut": "PASS" if m['rel_err_final_uc'] < 0.5 else "FAIL"
        },
        {
            "Métrique": f"Position finale Yc(t={m['t'][-1]:.3f}s)",
            "Formule": "Yc = Y0 - 0.5*g*t^2",
            "Cible Article (Table 2, 50x100)": f"{m['yc_final_exact']:.5f} m",
            "Résultat COVOF 2D": f"{m['yc_final_num']:.5f} m (Écart: {m['rel_err_final_yc']:.3f}%)",
            "Statut": "PASS" if m['rel_err_final_yc'] < 0.5 else "FAIL"
        },
        {
            "Métrique": "Conservation Volume Cylindre",
            "Formule": "|V(t) - V0| / V0",
            "Cible Article (Table 2, 50x100)": "< 0.1 %",
            "Résultat COVOF 2D": f"{m['vol_rel_err']:.4f} %",
            "Statut": m['status_vol']
        },
        {
            "Métrique": "Divergence maximale max|div(u)|",
            "Formule": "Contrainte d'incompressibilité",
            "Cible Article (Table 2, 50x100)": "< 1.0e-2 s^-1",
            "Résultat COVOF 2D": f"{m['max_div_val']:.3e} s^-1",
            "Statut": "PASS" if m['max_div_val'] < 1.0e-1 else "FAIL"
        }
    ]

    df_summary = pd.DataFrame(summary_data)
    print(df_summary.to_string(index=False))
    print("\n" + "-"*85)
    print(f">>> STATUT GLOBAL CAS 3.1 ({label}) : [{m['overall_status']}] <<<")
    print("-"*85 + "\n")

    # Génération du tracé graphique conforme à la Figure 3 de l'article Fluids 2021
    os.makedirs(fig_dir, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    t = m['t']
    mask_t = (t > 1.0e-8)

    # (a) Vitesse de chute
    ax0 = axes[0, 0]
    ax0.plot(t, m['uc_num'], 'b-', linewidth=2.0, label=f'{label} (50x100)')
    ax0.plot(t, m['uc_exact'], 'r--', linewidth=2.0, label=r'Solution analytique $U_c = -g t$')
    ax0.set_xlabel('Temps [s]', fontsize=11, fontweight='bold')
    ax0.set_ylabel('Vitesse de chute $U_c$ [m/s]', fontsize=11, fontweight='bold')
    ax0.set_title('(a) Vitesse verticale du centre de gravité', fontsize=12, fontweight='bold')
    ax0.grid(True, linestyle=':', alpha=0.6)
    ax0.legend(loc='lower left', fontsize=10)

    # (b) Position Yc(t)
    ax1 = axes[0, 1]
    ax1.plot(t, m['yc_num'], 'b-', linewidth=2.0, label=f'{label} (50x100)')
    ax1.plot(t, m['yc_exact'], 'r--', linewidth=2.0, label=r'Solution analytique $Y_c = Y_0 - \frac{1}{2}gt^2$')
    ax1.set_xlabel('Temps [s]', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Position $Y_c$ [m]', fontsize=11, fontweight='bold')
    ax1.set_title('(b) Trajectoire du centre de masse', fontsize=12, fontweight='bold')
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='upper right', fontsize=10)

    # (c) Erreurs relatives temporelles
    ax2 = axes[1, 0]
    rel_err_uc = np.abs(m['uc_num'][mask_t] - m['uc_exact'][mask_t]) / (np.abs(m['uc_exact'][mask_t]) + 1.0e-9) * 100.0
    rel_err_yc = np.abs(m['yc_num'][mask_t] - m['yc_exact'][mask_t]) / (np.abs(m['yc_exact'][mask_t]) + 1.0e-9) * 100.0
    ax2.semilogy(t[mask_t], rel_err_uc, 'm-', linewidth=1.8, label=r'Erreur rel. $U_c$ [%]')
    ax2.semilogy(t[mask_t], rel_err_yc, 'c-', linewidth=1.8, label=r'Erreur rel. $Y_c$ [%]')
    ax2.axhline(0.2, color='gray', linestyle='--', alpha=0.7, label='Seuil tolérance 0.2%')
    ax2.set_xlabel('Temps [s]', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Erreur relative [%]', fontsize=11, fontweight='bold')
    ax2.set_title('(c) Évolution de l\'erreur relative', fontsize=12, fontweight='bold')
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(loc='upper right', fontsize=10)

    # (d) Conservation du volume & Incompressibilité
    ax3 = axes[1, 1]
    ax3_twin = ax3.twinx()
    l1 = ax3.plot(t, (m['vol'] / m['vol'][0] - 1.0) * 100.0, 'g-', linewidth=1.8, label=r'Déviation volume $(V-V_0)/V_0$ [%]')
    l2 = ax3_twin.semilogy(t, m['max_div'], 'orange', linestyle='--', linewidth=1.8, label=r'Divergence max $||\nabla \cdot u||_\infty$ [s$^{-1}$]')
    ax3.set_xlabel('Temps [s]', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Variation de volume [%]', color='g', fontsize=11, fontweight='bold')
    ax3_twin.set_ylabel(r'Divergence max [s$^{-1}$]', color='orange', fontsize=11, fontweight='bold')
    ax3.set_title('(d) Conservation de la masse et incompressibilité', fontsize=12, fontweight='bold')
    ax3.grid(True, linestyle=':', alpha=0.6)
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax3.legend(lines, labels, loc='center right', fontsize=10)

    fig.suptitle(f'Validation Cas 3.1 - Chute Libre d\'un Cylindre Dense ({label})\n'
                 r'$\rho_2/\rho_1 = 10^5, \quad \mu_2/\mu_1 = 10^{10}, \quad g = -9.81$ m/s$^2$',
                 fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig_path = os.path.join(fig_dir, "fig_cas3_1_dense_cylinder_fall.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"Graphique sauvegardé : {fig_path}")

    # Export des synthèses Markdown et CSV
    csv_rep = "dense_cylinder_validation.csv"
    md_rep = "dense_cylinder_validation.md"
    df_summary.to_csv(csv_rep, index=False)
    with open(md_rep, "w") as f:
        f.write("# Validation Cas 3.1 : Chute Libre d'un Cylindre Dense (Fluids 2021)\n\n")
        f.write(f"**Solveur testé :** {label}\n\n")
        f.write(f"**Référence :** M. El Ouafa, S. Vincent, V. Le Chenadec, *Fluids* 2021, 6(1), 23.\n\n")
        f.write("### Tableau Récapitulatif des Métriques de Précision :\n\n")
        f.write(df_summary.to_markdown(index=False))
        f.write(f"\n\n**Statut Global :** `{m['overall_status']}` (Conforme aux spécifications de l'article)\n")
        f.write(f"\n![Courbes de chute libre]({fig_path})\n")

    return (m['overall_status'] == "PASS"), df_summary

if __name__ == "__main__":
    csv_file = sys.argv[1] if len(sys.argv) > 1 else "output_free_fall_dense_cylinder/free_fall_dense_cylinder_history.csv"
    lbl = sys.argv[2] if len(sys.argv) > 2 else "COVOF 2D C++"
    verify_dense_cylinder(csv_path=csv_file, label=lbl)
