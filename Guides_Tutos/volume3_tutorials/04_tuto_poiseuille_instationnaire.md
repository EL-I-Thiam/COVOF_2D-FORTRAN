# Tutoriel 4 : Écoulement de Poiseuille Plan Instationnaire (Démarrage Impulsionnel)

Ce tutoriel aborde le régime transitoire d'établissement visqueux d'un écoulement de Poiseuille plan à partir du repos sous l'effet soudain d'un forçage gravitationnel axial constant. La solution analytique en série de Fourier de Szymanski permet de valider avec une haute précision les schémas d'intégration en temps (`EULER` et `GEAR` / BDF2) et la convergence spatio-temporelle.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre la physique transitoire de diffusion visqueuse de quantité de mouvement depuis le repos.
- Configurer la périodicité axiale (`PERIODIC`) et les schémas temporels `EULER` et `GEAR` dans un fichier `.json`.
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Mesurer dynamiquement les ordres de convergence spatial ($\mathcal{O}(\Delta y^2)$) et temporel ($\mathcal{O}(\Delta t)$ vs $\mathcal{O}(\Delta t^2)$) face à la série analytique.
- Analyser l'évolution temporelle et les profils de vitesse sous Python et ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `DirectSolver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation` (`Simulation::export_poiseuille_unsteady_results`).
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation` (`sim_export_poiseuille_unsteady_results`).

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas étudie l'établissement transitoire d'un écoulement de Poiseuille plan partant du repos sous l'action brusque d'une force volumique axiale constante. Il valide le couplage entre le schéma d'intégration temporelle (`EULER` d'ordre 1 ou `GEAR` / BDF2 d'ordre 2) et la discrétisation spatiale visqueuse d'ordre 2.

Dans un canal périodique en $x$ de dimensions $L_x = 0.16\,\mathrm{m}$ et $H = L_y = 0.16\,\mathrm{m}$, un fluide au repos initial $u(x,y,0) = 0, v(x,y,0) = 0$ de propriétés $\rho = 1.0\,\mathrm{kg/m}^3$ et $\mu = 0.01\,\mathrm{Pa\cdot s}$ ($\nu = \mu/\rho = 0.01\,\mathrm{m^2/s}$) est soumis pour $t > 0$ à un forçage axial constant :
$$k = \rho g_x = 0.0555\,\mathrm{N/m}^3 \qquad (g_x = 0.0555\,\mathrm{m/s}^2).$$
L'équation d'évolution de la vitesse longitudinale $u(y,t)$ s'écrit :
$$\rho \frac{\partial u}{\partial t} = k + \mu \frac{\partial^2 u}{\partial y^2}, \qquad u(0,t) = u(H,t) = 0, \qquad u(y,0) = 0.$$

### 2.2 Solution Analytique Exacte et Lois d'Erreur
Par séparation des variables, la solution analytique exacte superpose le profil parabolique asymptotique ($u_{\max,\infty} = \frac{k H^2}{8\mu} = 0.01776\,\mathrm{m/s}$) et une série de Fourier amortie exponentiellement :
$$u_{\mathrm{ex}}(y, t) = \frac{k}{2\mu} y(H - y) - \frac{4 k H^2}{\mu \pi^3} \sum_{m=0}^{\infty} \frac{1}{(2m+1)^3} \sin\left(\frac{(2m+1)\pi y}{H}\right) \exp\left(-\nu \left(\frac{(2m+1)\pi}{H}\right)^2 t\right).$$
Avec le schéma d'Euler implicite (`EULER`) d'ordre 1 en temps et le schéma centré d'ordre 2 en espace, l'erreur globale obéit à :
$$\epsilon_{\mathrm{EULER}}(\Delta t, \Delta y) = \alpha_1 \Delta t + \gamma \Delta y^2 + \mathcal{O}(\Delta t \Delta y),$$
tandis qu'avec le schéma de Gear (`GEAR` / BDF2) d'ordre 2 en temps, l'erreur suit :
$$\epsilon_{\mathrm{GEAR}}(\Delta t, \Delta y) = \alpha_2 \Delta t^2 + \gamma \Delta y^2.$$

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine | $[0, L_x] \times [0, L_y]$ | $L_x = 0.16\,\mathrm{m}$, $L_y = H = 0.16\,\mathrm{m}$ |
| Propriétés et forçage | $\rho, \mu, g_x$ | $\rho = 1.0\,\mathrm{kg/m}^3$, $\mu = 0.01\,\mathrm{Pa\cdot s}$, $g_x = 0.0555\,\mathrm{m/s}^2$, $g_y = 0.0$ |
| Conditions aux limites | Gauche/Droite ; Bas/Haut | `PERIODIC` en $x$ ; `DIRICHLET` ($u=v=0$) en $y=0$ et $y=H$ |
| Campagne 1 (Matrice $N_y \times \Delta t$) | $t_{\mathrm{final}} = 1.0\,\mathrm{s}$ | $N_y \in \{16, 32, 64, 128, 256\} \times \Delta t \in \{10^{-1}, 10^{-2}, 10^{-3}, 10^{-4}\}\,\mathrm{s}$ ($N_x=16$) |
| Campagne 2 (`EULER` vs `GEAR`) | $t_{\mathrm{final}} = 1.0\,\mathrm{s}$, $N_y=256$ | $\Delta t \in \{0.1, 0.05, 0.025, 0.01\}\,\mathrm{s}$ |
| Schémas et Solveurs | Temps, Couplage, Krylov | `EULER`, `GEAR` ; `AUGMENTED_LAGRANGIAN`, `AAL`, `MONOLITHIC` avec `BICGSTAB2` + `BLOCK_TRIANGULAR` ou `DIRECT` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_4_unsteady_al_bicgstab2_ny64_dt0_01.json`
*(Pour Fortran 2018, utiliser une copie dans `case/fortran/` avec `"output_dir": "output/fortran/ny64_dt0_01"`. L'identifiant `"test_case": "CAS_3_8_POISEUILLE_UNSTEADY"` déclenche l'export de `poiseuille_unsteady_error.csv`.)*

```json
{
  "case_name": "cas_4_unsteady_al_bicgstab2_ny64_dt0_01",
  "mode": "NAVIER_STOKES",
  "test_case": "CAS_3_8_POISEUILLE_UNSTEADY",
  "mesh": {
    "nx": 16,
    "ny": 64,
    "lx": 0.16,
    "ly": 0.16
  },
  "physics": {
    "enable_two_phase": false,
    "enable_surface_tension": false,
    "enable_convection": false,
    "gx": 0.0555,
    "gy": 0.0,
    "fluid1": {
      "rho": 1.0,
      "mu": 0.01
    }
  },
  "time": {
    "time_scheme": "EULER",
    "dt": 0.01,
    "t_final": 1.0,
    "cfl_max": 0.5,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "PERIODIC",  "u": 0.0, "v": 0.0 },
    "right":  { "type": "PERIODIC",  "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-12,
    "max_iter": 5000,
    "r_al": 100.0,
    "uzawa_iter": 1,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "output": {
    "output_dir": "output/cpp/ny64_dt0_01",
    "case_name": "cas_4_unsteady_al_bicgstab2_ny64_dt0_01",
    "vtk_interval": 25,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

L'exécution parallèle utilise **MPI pur** (`mpirun -np 4`) pour les solveurs de Krylov, et `mpirun -np 1` si le solveur `DIRECT` (mise en cache LU dans `DirectSolver` lorsque $\Delta t$, $\rho$ et $\mu$ sont constants) ou `STANDARD_PROJECTION` est sélectionné.

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/ny64_dt0_01 output/fortran/ny64_dt0_01

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_4_unsteady_al_bicgstab2_ny64_dt0_01.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_4_unsteady_al_bicgstab2_ny64_dt0_01.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Convergence spatiale d'ordre 2 ($\Delta t \ll \Delta y^2/\nu$) :**
   Pour $\Delta t = 10^{-4}\,\mathrm{s}$, extraire $E_{L_1}, E_{L_2}, E_{L_\infty}$ depuis `poiseuille_unsteady_error.csv` sur $N_y \in \{16, 32, 64, 128, 256\}$ et vérifier que l'ordre spatial $p_{\mathrm{espace}} = \ln(E(2\Delta y)/E(\Delta y)) / \ln(2) \in [1.95, 2.05]$.
2. **Convergence temporelle d'ordre 1 (`EULER`) vs ordre 2 (`GEAR`) :**
   Sur grille fine $N_y = 256$, calculer l'ordre temporel $p_{\mathrm{temps}} = \ln(E(2\Delta t)/E(\Delta t)) / \ln(2)$ entre $\Delta t \in \{0.1, 0.05, 0.025, 0.01\}\,\mathrm{s}$. Vérifier $p_{\mathrm{temps}} \approx 1.0$ pour `EULER` et $p_{\mathrm{temps}} \approx 2.0$ pour `GEAR`.
3. **Invariance des couplages et parité C++17 / Fortran 2018 :**
   Vérifier la cohérence entre `AUGMENTED_LAGRANGIAN`, `AAL` et `MONOLITHIC` ainsi que la parité C++17 / Fortran 2018 : $|E_{L_1}^{\mathrm{C++}} - E_{L_1}^{\mathrm{Fortran}}| < 10^{-13}$.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` (`9` colonnes en C++, `11` colonnes en Fortran) et `<case_name>_<step:06d>.pvtr` / `<case_name>_p<rank>_<step:06d>.vtr`.
- `poiseuille_unsteady_error.csv` : `nx,ny,dx,dy,dt,t_final,iters,err_l1,err_l2,err_linf,u_max_num,u_max_exact,cpu_time`.
- `profile_u_y.csv` : `y,u_num,u_exact,error` à $t = t_{\mathrm{final}}$.
- `solution_field_2d.csv` (`i,j,x,y,u,v,p,u_exact,error_u`) et `solution_field.vtk`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 4 (Poiseuille instationnaire)."""
import csv
import math
from pathlib import Path

def read_unsteady_error(out_dir: str) -> dict:
    with (Path(out_dir) / "poiseuille_unsteady_error.csv").open(newline="") as f:
        row = next(csv.DictReader(f))
    return {k: float(v) for k, v in row.items()}

if __name__ == "__main__":
    cases = ["output/cpp/ny16_dt1e-4", "output/cpp/ny32_dt1e-4",
             "output/cpp/ny64_dt1e-4", "output/cpp/ny128_dt1e-4"]
    prev = None
    for d in cases:
        if not (Path(d) / "poiseuille_unsteady_error.csv").exists():
            continue
        m = read_unsteady_error(d)
        p = math.log(prev["err_l2"] / m["err_l2"]) / math.log(prev["dy"] / m["dy"]) if prev else float("nan")
        print(f"ny={int(m['ny']):3d} | dt={m['dt']:.1e} | L2={m['err_l2']:.4e} "
              f"| u_max={m['u_max_num']:.6f} (ex={m['u_max_exact']:.6f}) | p={p:.3f}")
        prev = m
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série temporelle `<case_name>_..pvtr` dans ParaView (**File -> Open**) pour animer l'établissement transitoire de la couche limite pariétale vers l'axe central.
2. Appliquer le filtre **Plot Over Line** verticalement en $x = 0.08\,\mathrm{m}$ sur plusieurs pas de temps, ou comparer directement `u_num` et `u_exact` depuis `solution_field.vtk`.
