# Tutoriel 1 : Opérateur Laplacien avec Conditions de Dirichlet Homogènes

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_1_laplacian/` *(validation rapide `Run_Valide-CASES`, cas 1)*.

Ce premier tutoriel constitue le point de départ fondamental de la vérification de COVOF 2D. Il permet de valider la discrétisation spatiale à l'ordre 2 de l'opérateur Laplacien scalaire 2D (équation de Poisson) avec conditions aux limites de Dirichlet homogènes sur grille cartésienne structurée.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements mathématiques et la discrétisation aux différences/volumes finis à 5 points de l'équation de Poisson 2D.
- Configurer un fichier d'entrée `.json` conforme en mode `"LAPLACIAN"` pour les moteurs C++17 et Fortran 2018.
- Exécuter le cas test avec les binaires précompilés `covof2d_cpp` et `covof2d_fortran` via `mpirun`.
- Quantifier dynamiquement l'erreur discrète ($E_{L_1}$, $E_{L_2}$, $E_{L_\infty}$) et l'ordre de convergence de Richardson face à la solution analytique exacte.
- Exploiter les sorties `.csv` et `.vtk` sous Python et ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `CSRMatrix`, `Simulation` (`Simulation::run_laplacian`).
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_simulation` (`sim_run_laplacian`).

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce premier cas vérifie la consistance, la précision et l'ordre de convergence spatial de la discrétisation de l'opérateur Laplacien bidimensionnel stationnaire (équation de Poisson) avec conditions aux limites de Dirichlet homogènes.

Sur le domaine rectangulaire $\Omega = [0, a] \times [0, b]$ de dimensions $a = L_x = 1.0\,\mathrm{m}$ et $b = L_y = 0.5\,\mathrm{m}$, on résout :
$$\Delta \Psi = \frac{\partial^2 \Psi}{\partial x^2} + \frac{\partial^2 \Psi}{\partial y^2} = S \quad \text{dans } \Omega = ]0, a[ \times ]0, b[, \qquad \text{avec } S = 1.0\,\mathrm{m}^{-2},$$
soumis à la condition de Dirichlet homogène sur l'ensemble du bord $\partial\Omega$ :
$$\Psi(0, y) = \Psi(a, y) = \Psi(x, 0) = \Psi(x, b) = 0, \qquad \forall (x,y) \in \partial\Omega.$$

### 2.2 Solution Analytique Exacte et Discrétisation à 5 Points
Par décomposition sur la base propre orthogonale de l'opérateur Laplacien avec conditions de Dirichlet homogènes sur $[0,a]\times[0,b]$, seuls les modes impairs $(2k_1+1, 2k_2+1)$ contribuent au développement en série de Fourier double :
$$\Psi_{\mathrm{ex}}(x,y) = -\frac{16 S}{\pi^4} \sum_{k_1=0}^{\infty} \sum_{k_2=0}^{\infty} \frac{\sin\left(\frac{(2k_1+1)\pi x}{a}\right)\sin\left(\frac{(2k_2+1)\pi y}{b}\right)}{(2k_1+1)(2k_2+1)\left[\left(\frac{2k_1+1}{a}\right)^2 + \left(\frac{2k_2+1}{b}\right)^2\right]}.$$
Dans `Simulation::run_laplacian` et `sim_run_laplacian`, cette série est évaluée aux nœuds intérieurs $(x_i, y_j) = (i\Delta x, j\Delta y)$ pour $i \in \{1, \dots, N_x-1\}$ et $j \in \{1, \dots, N_y-1\}$ en tronquant la somme à $60 \times 60$ modes impairs ($k_1, k_2 \in \{0, \dots, 59\}$).

Sur un maillage uniforme à mailles carrées $\Delta x = a/N_x = \Delta y = b/N_y$, l'opérateur $-\Delta_h \Psi = -S$ est assemblé sous forme de matrice creuse symétrique définie positive (SPD) à 5 points d'ordre 2 :
$$\frac{-\Psi_{i-1,j} + 2\Psi_{i,j} - \Psi_{i+1,j}}{\Delta x^2} + \frac{-\Psi_{i,j-1} + 2\Psi_{i,j} - \Psi_{i,j+1}}{\Delta y^2} = -S,$$
où les valeurs nulles imposées aux nœuds de bord ($i \in \{0, N_x\}$ ou $j \in \{0, N_y\}$) s'éliminent directement du système intérieur de taille $(N_x - 1)(N_y - 1)$.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Géométrie | Domaine $[0, L_x] \times [0, L_y]$ | $L_x = 1.0\,\mathrm{m}$, $L_y = 0.5\,\mathrm{m}$ (rapport d'aspect $2:1$) |
| Séquence de grilles | $(N_x \times N_y)$ à mailles carrées | $16\times 8$, $32\times 16$, $64\times 32$, $128\times 64$, $256\times 128$ |
| Pas d'espace | $\Delta x = \Delta y$ | $L_x / N_x = L_y / N_y \in \{6.25\times 10^{-2}, \dots, 3.90625\times 10^{-3}\}\,\mathrm{m}$ |
| Terme source | $S$ (`source_term`) | $1.0\,\mathrm{m}^{-2}$ |
| Conditions aux limites | Gauche, Droite, Bas, Haut | `DIRICHLET` homogène ($\Psi = 0$) |
| Solveurs disponibles en mode `LAPLACIAN` | `solver.linear_solver` | `DIRECT` (Cholesky bande SPD de demi-largeur $N_x - 1$), `BICGSTAB2` (précond. Jacobi), `CG` (précond. Jacobi) |
| Tolérance et itérations | `solver.tol`, `solver.max_iter` | $\mathrm{tol} = 10^{-12}$, $\mathrm{max\_iter} = 5000$ |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype Solveur Direct (`case/cpp/cas_1_grid64_direct.json`)
*(Pour Fortran 2018, placer une copie dans `case/fortran/` en remplaçant `"output/cpp/grid64_direct"` par `"output/fortran/grid64_direct"`.)*

```json
{
  "case_name": "cas_1_grid64_direct",
  "mode": "LAPLACIAN",
  "test_case": "CAS_1_DIRICHLET",
  "source_term": 1.0,
  "mesh": {
    "nx": 64,
    "ny": 32,
    "lx": 1.0,
    "ly": 0.5
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "right":  { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "DIRECT",
    "linear_solver": "DIRECT",
    "precond_type": "NONE",
    "use_precond": false,
    "tol": 1e-12,
    "max_iter": 5000
  },
  "output": {
    "output_dir": "output/cpp/grid64_direct",
    "case_name": "cas_1_grid64_direct",
    "vtk_interval": 1,
    "log_interval": 1,
    "write_csv_history": true
  }
}
```

### 4.2 Variante Solveur de Krylov (`case/cpp/cas_1_grid64_cg_jacobi.json`)
Pour tester le Gradient Conjugué préconditionné par Jacobi (ou `BICGSTAB2`), remplacer le bloc `"solver"` par :

```json
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "linear_solver": "CG",
    "precond_type": "JACOBI",
    "use_precond": true,
    "tol": 1e-12,
    "max_iter": 5000
  }
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs reposent exclusivement sur **MPI pur**. Le mode `"LAPLACIAN"` assemblant un système global monolithique sur le domaine complet, l'exécution de référence s'effectue avec `mpirun -np 1` en s'assurant que le répertoire de sortie cible existe.

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/grid64_direct output/fortran/grid64_direct

# 2. Exécution avec le moteur C++17 (MPI séquentiel pour mode LAPLACIAN / DIRECT)
mpirun -np 1 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_1_grid64_direct.json

# 3. Exécution avec le moteur Fortran 2018 (MPI séquentiel pour mode LAPLACIAN / DIRECT)
mpirun -np 1 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_1_grid64_direct.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Convergence spatiale d'ordre 2 (extrapolation de Richardson) :**
   À partir des erreurs discrètes extraites dynamiquement depuis `laplacian_error.csv` pour deux maillages successifs de pas $h_1$ et $h_2 = h_1 / 2$, l'ordre observé est :
   $$p_\alpha = \frac{\ln\big(E_{L_\alpha}(h_1) / E_{L_\alpha}(h_2)\big)}{\ln(h_1 / h_2)}, \qquad \alpha \in \{1, 2, \infty\}.$$
   Critère de validation : $p_\alpha \in [1.95, 2.05]$ dans le régime asymptotique pour les deux moteurs.
2. **Cohérence entre solveurs `DIRECT`, `CG` et `BICGSTAB2` :**
   Vérifier sur une même grille que le solveur direct (Cholesky bande, `iters = 1`) et les solveurs itératifs convergent vers la même solution discrète à $\mathrm{tol} = 10^{-12}$ près :
   $$|E_{L_2}^{\mathrm{Krylov}} - E_{L_2}^{\mathrm{DIRECT}}| < 10^{-10}.$$
3. **Symétrie géométrique quatre-quadrants et parité C++17 / Fortran 2018 :**
   Vérifier à partir de `solution_field_2d.csv` que $|\Psi(x,y) - \Psi(a-x,y)| < 10^{-12}$, $|\Psi(x,y) - \Psi(x,b-y)| < 10^{-12}$, et que l'écart maximal entre les deux binaires satisfait $\|\Psi^{\mathrm{C++}} - \Psi^{\mathrm{Fortran}}\|_\infty < 10^{-11}$.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `laplacian_error.csv` : `nx,ny,dx,dy,iters,err_l1,err_l2,err_linf` (plus `cpu_time` en Fortran).
- `solver_history.csv` : historique des résidus `iter,residual` (en C++ pour `DIRECT`/`CG`/`BICGSTAB2`, et pour `DIRECT`/`CG` en Fortran).
- `solution_field_2d.csv` : champ complet `i,j,x,y,psi_num,psi_ex,error_abs`.
- `profile_x_mid.csv` et `profile_y_mid.csv` : coupes 1D selon $y = b/2$ et $x = a/2$ (`x,y,psi_num,psi_ex,error_abs`).
- `solution_field.vtk` : fichier VTK `STRUCTURED_POINTS` contenant `psi_numeric`, `psi_exact`, `error_abs`.
- `<case_name>_history.csv` : en-tête initialisé par `VTKWriter` (`9` colonnes en C++, `11` colonnes en Fortran).

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 1 (Laplacien Dirichlet) sur fichiers réels."""
import csv
import math
from pathlib import Path

def read_laplacian_error(output_dir: str) -> dict:
    path = Path(output_dir) / "laplacian_error.csv"
    with path.open(newline="") as f:
        row = next(csv.DictReader(f))
    return {k: float(v) for k, v in row.items()}

def check_symmetry_and_parity(dir_cpp: str, dir_ftn: str) -> None:
    def load_field(d: str):
        pts = {}
        with (Path(d) / "solution_field_2d.csv").open(newline="") as f:
            for r in csv.DictReader(f):
                pts[(int(r["i"]), int(r["j"]))] = float(r["psi_num"])
        return pts

    cpp = load_field(dir_cpp)
    ftn = load_field(dir_ftn)
    nx = max(i for i, _ in cpp)
    ny = max(j for _, j in cpp)
    sym_x = max(abs(cpp[(i, j)] - cpp[(nx - i, j)]) for i, j in cpp)
    sym_y = max(abs(cpp[(i, j)] - cpp[(i, ny - j)]) for i, j in cpp)
    parity = max(abs(cpp[k] - ftn[k]) for k in cpp)
    print(f"Symétrie X={sym_x:.3e}, Symétrie Y={sym_y:.3e}, Parité C++/Fortran={parity:.3e}")

if __name__ == "__main__":
    grids = [16, 32, 64, 128, 256]
    prev = None
    for nx in grids:
        d = f"output/cpp/grid{nx}_direct"
        if not (Path(d) / "laplacian_error.csv").exists():
            continue
        m = read_laplacian_error(d)
        order = math.log(prev["err_l2"] / m["err_l2"]) / math.log(prev["dx"] / m["dx"]) if prev else float("nan")
        print(f"nx={int(m['nx']):3d} | dx={m['dx']:.5e} | L2={m['err_l2']:.5e} | Linf={m['err_linf']:.5e} | p_L2={order:.3f}")
        prev = m
```

### 6.4 Visualisation sous ParaView
1. Ouvrir `output/cpp/grid64_direct/solution_field.vtk` avec **File -> Open** puis cliquer sur **Apply**.
2. Sélectionner le champ scalaire `psi_numeric` (ou `error_abs`) dans la barre d'outils et appliquer le filtre **Warp By Scalar** (facteur d'échelle adapté) pour visualiser la surface $\Psi(x,y)$ en 3D.
3. Utiliser le filtre **Plot Over Line** le long de $y = 0.25\,\mathrm{m}$ et $x = 0.50\,\mathrm{m}$ pour superposer `psi_numeric` et `psi_exact`.
