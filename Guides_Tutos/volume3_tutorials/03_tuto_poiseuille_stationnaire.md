# Tutoriel 3 : Écoulement de Poiseuille Plan Laminaire Stationnaire

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_3_poiseuille/` *(validation rapide `Run_Valide-CASES`, cas 3)*.

Ce tutoriel traite de l'écoulement laminaire cisaillé 1D stationnaire en canal plan entre deux parois immobiles, entraîné par un profil parabolique imposé en entrée et équilibré par un gradient de pression longitudinal. Il valide la discrétisation de Navier-Stokes incompressible sur maillage décalé MAC et l'extrapolation quadratique à 3 points dans les mailles fantômes pariétales.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les équations de Navier-Stokes incompressibles en canal plan laminaire et le couplage vitesse-pression.
- Configurer correctement les conditions aux limites `DIRICHLET` (`"profile": "PARABOLIC"`, `"u_mean"`) et `NEUMANN` sur grille MAC.
- Exécuter le cas test en séquentiel et en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Vérifier dynamiquement l'exactitude machine du profil parabolique et du gradient de pression axial.
- Visualiser les champs 2D et extraire les profils transverses sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `DirectSolver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide la résolution stationnaire des équations de Navier-Stokes incompressibles en canal plan laminaire et vérifie la propriété d'exactitude polynomiale d'ordre 2 du stencil visqueux couplé à l'extrapolation quadratique à 3 points dans les mailles fantômes de paroi.

Dans un canal plan $\Omega = [0, L_x] \times [0, H]$ de longueur $L_x = 1.0\,\mathrm{m}$ et de hauteur $H = L_y = 0.1\,\mathrm{m}$, un fluide newtonien de masse volumique $\rho = 1.0\,\mathrm{kg/m}^3$ et de viscosité $\mu = 1.0\times 10^{-3}\,\mathrm{Pa\cdot s}$ entre par la section $x=0$ avec une vitesse débitante moyenne $V_0 = u_{\mathrm{mean}} = 1.0\times 10^{-3}\,\mathrm{m/s}$, correspondant à une vitesse maximale sur l'axe $u_{\max} = \frac{3}{2} V_0 = 1.5\times 10^{-3}\,\mathrm{m/s}$ et un nombre de Reynolds $\mathrm{Re}_H = \rho V_0 H / \mu = 0.1$.

### 2.2 Solution Analytique Exacte et Extrapolation Pariétale
Le profil de vitesse longitudinale est purement parabolique en $y$ et le champ de pression décroît linéairement en $x$ :
$$u_{\mathrm{ex}}(y) = \frac{6 V_0}{H^2} y (H - y) = \frac{4 u_{\max}}{H^2} y (H - y), \qquad v_{\mathrm{ex}}(x,y) = 0,$$
$$G_{x,\mathrm{ex}} = -\frac{\mathrm{d}p_{\mathrm{ex}}}{\mathrm{d}x} = \frac{12 \mu V_0}{H^2} = \frac{8 \mu u_{\max}}{H^2} = 1.2\times 10^{-3}\,\mathrm{Pa/m}, \qquad p_{\mathrm{ex}}(x) = G_{x,\mathrm{ex}} (L_x - x).$$

Sur le maillage décalé MAC, les vitesses $u_{i,j}$ sont situées à distance $\Delta y / 2$ des parois $y = 0$ et $y = H$. Dans `MpiComm.cpp` et `module_mpi_comm.f90`, COVOF 2D applique l'extrapolation quadratique de Lagrange à 3 points :
$$u_{i,0} = \frac{8}{3} u_{\mathrm{wall}} - 2 u_{i,1} + \frac{1}{3} u_{i,2}, \qquad u_{i,N_y+1} = \frac{8}{3} u_{\mathrm{wall}} - 2 u_{i,N_y} + \frac{1}{3} u_{i,N_y-1}.$$
En injectant cette relation dans le stencil visqueux centré au premier nœud intérieur $j=1$, l'approximation discrète devient exacte pour tout polynôme $u(y) \in \mathbb{P}_2$, garantissant une erreur spatiale au niveau du zéro machine sur tout maillage.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine | $[0, L_x] \times [0, L_y]$ | $L_x = 1.0\,\mathrm{m}$, $L_y = H = 0.1\,\mathrm{m}$ |
| Propriétés du fluide | $\rho, \mu, V_0$ | $\rho = 1.0\,\mathrm{kg/m}^3$, $\mu = 10^{-3}\,\mathrm{Pa\cdot s}$, $V_0 = u_{\mathrm{mean}} = 10^{-3}\,\mathrm{m/s}$ |
| Séquence de grilles | $(N_x \times N_y)$ | $16\times 16$, $32\times 32$, $64\times 64$, $128\times 128$ |
| Conditions aux limites | Gauche, Droite, Bas/Haut | `DIRICHLET` (`"profile": "PARABOLIC"`, `"u_mean": 0.001`), `NEUMANN` (sortie libre), `DIRICHLET` ($u=v=0$) |
| Couplage vitesse-pression | `solver.type` | `AUGMENTED_LAGRANGIAN`, `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (`AAL`), `MONOLITHIC`, `STANDARD_PROJECTION` |
| Solveurs linéaires & Précond. | `linear_solver` / `precond_type` | `DIRECT`, `BICGSTAB2`, `BICGSTAB`, `GMRES`, `MINRES` avec `NONE`, `JACOBI`, `SSOR`, `ILU0`, `BLOCK_TRIANGULAR` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_3_poiseuille_grid64_al.json`
*(Pour Fortran 2018, utiliser une copie dans `case/fortran/` avec `"output_dir": "output/fortran/grid64_al"`. Noter que `"test_case": "CAS_3_7_POISEUILLE"` active l'export automatique des fichiers `poiseuille_error.csv`, `profile_u_y.csv` et `profile_p_x.csv`.)*

```json
{
  "case_name": "cas_3_poiseuille_grid64_al",
  "mode": "NAVIER_STOKES",
  "test_case": "CAS_3_7_POISEUILLE",
  "steady_state_tol": 1e-13,
  "mesh": {
    "nx": 64,
    "ny": 64,
    "lx": 1.0,
    "ly": 0.1
  },
  "physics": {
    "enable_two_phase": false,
    "enable_surface_tension": false,
    "enable_convection": true,
    "convective_scheme": "UPWIND",
    "convective_linearization": "PICARD",
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": {
      "rho": 1.0,
      "mu": 0.001
    }
  },
  "time": {
    "time_scheme": "EULER",
    "dt": 0.1,
    "t_final": 0.5,
    "cfl_max": 0.5,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u": 0.001, "u_mean": 0.001, "v": 0.0, "profile": "PARABOLIC" },
    "right":  { "type": "NEUMANN",   "u": 0.0,   "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0,   "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0,   "v": 0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "linear_solver": "BICGSTAB2",
    "precond_type": "BLOCK_TRIANGULAR",
    "use_precond": true,
    "tol": 1e-12,
    "max_iter": 1000,
    "r_al": 1.0,
    "uzawa_iter": 1
  },
  "output": {
    "output_dir": "output/cpp/grid64_al",
    "case_name": "cas_3_poiseuille_grid64_al",
    "vtk_interval": 5,
    "log_interval": 1,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur**. Avec les solveurs de Krylov (`BICGSTAB2`, `GMRES`, `BICGSTAB`, et `MINRES` uniquement sur systèmes symétriques sans convection ni préconditionneur asymétrique `BLOCK_TRIANGULAR`), le calcul s'exécute en parallèle sur plusieurs rangs MPI (`mpirun -np 4`) ; avec `DIRECT` ou `STANDARD_PROJECTION`, utiliser `mpirun -np 1`.

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/grid64_al output/fortran/grid64_al

# 2. Exécution parallèle MPI avec le moteur C++17 (solveurs Krylov)
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_3_poiseuille_grid64_al.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018 (solveurs Krylov)
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_3_poiseuille_grid64_al.json

# 4. Exécution séquentielle MPI (pour variantes linear_solver: "DIRECT" ou solver.type: "STANDARD_PROJECTION")
mpirun -np 1 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_3_poiseuille_grid64_direct.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Exactitude machine du profil de vitesse et du gradient de pression :**
   À partir de `poiseuille_error.csv`, vérifier sur toutes les grilles ($16\times 16$ à $128\times 128$) que $E_{L_1}, E_{L_2}, E_{L_\infty} < 10^{-11}\,\mathrm{m/s}$ et que la composante transversale satisfait $\|v\|_\infty < 10^{-13}\,\mathrm{m/s}$.
   Vérifier également que le gradient axial de pression numérique $G_{x,\mathrm{num}}$ restitué dans `profile_p_x.csv` coïncide avec $G_{x,\mathrm{ex}} = 12\mu V_0 / H^2$ (`err_grad_p` relatif $< 10^{-9}$).
2. **Invariance des schémas convectifs et comparaison des couplages :**
   Puisque $\partial u_{\mathrm{ex}}/\partial x = 0$ et $v_{\mathrm{ex}} = 0$, vérifier l'invariance vis-à-vis de `"convective_scheme"` (`UPWIND`, `CENTERED`, `TVD_SUPERBEE` — tous trois effectivement câblés dans `apply_system_matrix`) et `"convective_linearization"` (`PICARD`, `ADAMS_BASHFORTH`, `NEWTON`), ainsi que la suprématie du préconditionneur MPI-parallèle `BLOCK_TRIANGULAR` sur `JACOBI`, `SSOR` et `ILU0`.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier que $\|u^{\mathrm{C++}} - u^{\mathrm{Fortran}}\|_\infty < 10^{-12}\,\mathrm{m/s}$ et $\|p^{\mathrm{C++}} - p^{\mathrm{Fortran}}\|_\infty < 10^{-12}\,\mathrm{Pa}$.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` (`9` colonnes en C++, `11` colonnes en Fortran) et `<case_name>_<step:06d>.pvtr` / `<case_name>_p<rank>_<step:06d>.vtr`.
- `poiseuille_error.csv` : `nx,ny,dx,dy,iters,err_l1,err_l2,err_linf,err_grad_p,cpu_time`.
- `profile_u_y.csv` : `y,u_mid,u_three_quarter,u_ex,error` (coupes à $x=L_x/2$ et $x=3L_x/4$).
- `profile_p_x.csv` : `x,p_mid,p_ex,error` (coupe axiale à $y=H/2$).
- `solution_field_2d.csv` (`i,j,x,y,u,v,p,u_ex,p_ex,error_u`) et `solution_field.vtk` (`Velocity_U`, `Pressure`, `Error_U`).

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 3 (Poiseuille stationnaire)."""
import csv
from pathlib import Path

def verify_poiseuille(out_dir: str, mu: float = 1e-3, v0: float = 1e-3, ly: float = 0.1) -> None:
    d = Path(out_dir)
    with (d / "poiseuille_error.csv").open(newline="") as f:
        err = next(csv.DictReader(f))
    gx_ex = 12.0 * mu * v0 / (ly ** 2)
    rel_gp = float(err["err_grad_p"]) / gx_ex
    print(f"[{out_dir}] L1={float(err['err_l1']):.3e} | L2={float(err['err_l2']):.3e} "
          f"| Linf={float(err['err_linf']):.3e} | rel_err(dP/dx)={rel_gp:.3e}")

    with (d / "solution_field_2d.csv").open(newline="") as f:
        max_v = max(abs(float(r["v"])) for r in csv.DictReader(f))
    print(f"[{out_dir}] max|v| = {max_v:.3e} m/s")

if __name__ == "__main__":
    verify_poiseuille("output/cpp/grid64_al")
    verify_poiseuille("output/fortran/grid64_al")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir `output/cpp/grid64_al/cas_3_poiseuille_grid64_al_000005.pvtr` (ou `solution_field.vtk`) avec **File -> Open**.
2. Colorer le domaine par `Pressure` pour contrôler la décroissance linéaire amont-aval, puis appliquer **Plot Over Line** verticalement en $x = 0.5\,\mathrm{m}$ sur `Velocity` (composante $X$) afin de vérifier la superposition au profil parabolique exact.
