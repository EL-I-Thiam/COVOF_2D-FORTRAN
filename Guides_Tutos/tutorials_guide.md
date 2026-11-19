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
# Tutoriel 2 : Opérateur Laplacien avec Conditions de Neumann Non Homogènes et Jauge de Pression

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_2_neumann/` *(validation rapide `Run_Valide-CASES`, cas 2)*.

Ce tutoriel aborde la résolution d'une équation de Laplace sous conditions aux limites de Neumann pur. Il met en lumière deux aspects cruciaux du calcul CFD : le traitement rigoureux des flux pariétaux sur grille structurée et la levée du mode singulier à constante (jauge de pression) pour garantir l'unicité de la solution et la convergence du solveur.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre la compatibilité de flux d'un problème de Neumann pur et l'ancrage de la constante de jauge $\Psi(0,0) = 0$.
- Configurer un fichier `.json` conforme en mode `"LAPLACIAN"` avec `"test_case": "CAS_3_2_NEUMANN"`.
- Exécuter le cas test avec les binaires précompilés `covof2d_cpp` et `covof2d_fortran` via `mpirun`.
- Vérifier dynamiquement la propriété d'exactitude polynomiale machine sur $\Psi_{\mathrm{ex}}(x,y) = x^2 - y^2$.
- Exploiter les fichiers `.csv` et `.vtk` sous Python et ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `CSRMatrix`, `Simulation` (`Simulation::run_laplacian_cas3_2`).
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_simulation` (`sim_run_laplacian_cas3_2`).

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas vérifie le traitement des conditions aux limites de Neumann non homogènes (flux normaux imposés sur $\partial\Omega$) et la régularisation de jauge d'une matrice Laplacienne singulière dont le noyau est engendré par les fonctions constantes.

Sur le domaine carré unitaire $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$, on résout l'équation de Laplace homogène ($S = 0$) :
$$\Delta \Psi = \frac{\partial^2 \Psi}{\partial x^2} + \frac{\partial^2 \Psi}{\partial y^2} = 0 \quad \text{dans } \Omega,$$
avec les conditions aux limites de Neumann suivantes sur les quatre frontières :
$$\left.\frac{\partial \Psi}{\partial n}\right|_{x=0} = -\left.\frac{\partial \Psi}{\partial x}\right|_{x=0} = 0, \qquad \left.\frac{\partial \Psi}{\partial n}\right|_{x=1} = +\left.\frac{\partial \Psi}{\partial x}\right|_{x=1} = +2,$$
$$\left.\frac{\partial \Psi}{\partial n}\right|_{y=0} = -\left.\frac{\partial \Psi}{\partial y}\right|_{y=0} = 0, \qquad \left.\frac{\partial \Psi}{\partial n}\right|_{y=1} = +\left.\frac{\partial \Psi}{\partial y}\right|_{y=1} = -2.$$
La condition de compatibilité $\oint_{\partial\Omega} \frac{\partial\Psi}{\partial n}\,\mathrm{d}\gamma = 0 + 2(1) + 0 - 2(1) = 0$ est exactement satisfaite. En fixant la constante de jauge par la condition d'ancrage $\Psi(0,0) = 0$ (ligne $0$ de la matrice fixée à l'identité dans `run_laplacian_cas3_2` / `sim_run_laplacian_cas3_2`), la solution analytique exacte est le polynôme harmonique quadratique :
$$\Psi_{\mathrm{ex}}(x,y) = x^2 - y^2.$$

### 2.2 Propriété Théorique d'Exactitude Machine
Le développement de Taylor du schéma centré à 5 points montre que l'erreur de troncature locale s'écrit :
$$\tau_{i,j} = \frac{\Delta x^2}{12}\frac{\partial^4 \Psi_{\mathrm{ex}}}{\partial x^4} + \frac{\Delta y^2}{12}\frac{\partial^4 \Psi_{\mathrm{ex}}}{\partial y^4} + \mathcal{O}(\Delta x^4) \equiv 0,$$
car toutes les dérivées d'ordre supérieur ou égal à $3$ de $\Psi_{\mathrm{ex}}(x,y) = x^2 - y^2$ sont identiquement nulles, et le traitement par miroir d'ordre 2 aux nœuds de Neumann ($2/\Delta x^2$ et second membre $\pm 4/\Delta x$) est également exact pour les polynômes quadratiques. L'erreur discrète $E_{L_1}, E_{L_2}, E_{L_\infty}$ doit donc atteindre le zéro machine ($\sim 10^{-14} - 10^{-12}$) quelle que soit la finesse du maillage.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Géométrie | Domaine $[0, L_x] \times [0, L_y]$ | $L_x = 1.0\,\mathrm{m}$, $L_y = 1.0\,\mathrm{m}$ |
| Séquence de grilles | $(N_x \times N_y)$ | $16\times 16$, $32\times 32$, $64\times 64$, $128\times 128$, $256\times 256$ |
| Terme source | $S$ (`source_term`) | $0.0\,\mathrm{m}^{-2}$ |
| Dispatching cas analytique | `test_case` | `"CAS_3_2_NEUMANN"` (active les flux de Neumann exacts et l'ancrage $\Psi(0,0)=0$) |
| Conditions aux limites | `boundary_conditions` | `NEUMANN` sur les 4 bords |
| Solveurs disponibles | `solver.linear_solver` | `DIRECT` (factorisation LU bande de largeur $2N_x + 3$) et `BICGSTAB2` (préconditionneur ligne tridiagonal Thomas) |
| Tolérance et itérations | `solver.tol`, `solver.max_iter` | $\mathrm{tol} = 10^{-14}$, $\mathrm{max\_iter} = 10000$ |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype Solveur Direct (`case/cpp/cas_2_neumann_grid64_direct.json`)
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/grid64_direct"`.)*

```json
{
  "case_name": "cas_2_neumann_grid64_direct",
  "mode": "LAPLACIAN",
  "test_case": "CAS_3_2_NEUMANN",
  "source_term": 0.0,
  "mesh": {
    "nx": 64,
    "ny": 64,
    "lx": 1.0,
    "ly": 1.0
  },
  "boundary_conditions": {
    "left":   { "type": "NEUMANN", "u": 0.0, "v": 0.0 },
    "right":  { "type": "NEUMANN", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "NEUMANN", "u": 0.0, "v": 0.0 },
    "top":    { "type": "NEUMANN", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "DIRECT",
    "linear_solver": "DIRECT",
    "precond_type": "NONE",
    "use_precond": false,
    "tol": 1e-14,
    "max_iter": 10000
  },
  "output": {
    "output_dir": "output/cpp/grid64_direct",
    "case_name": "cas_2_neumann_grid64_direct",
    "vtk_interval": 1,
    "log_interval": 1,
    "write_csv_history": true
  }
}
```

### 4.2 Variante Solveur Itératif `BICGSTAB2` (`case/cpp/cas_2_neumann_grid64_bicgstab2.json`)
Pour activer le solveur itératif `BiCGSTAB(2)` couplé au préconditionneur ligne tridiagonal (algorithme de Thomas ligne par ligne), remplacer le bloc `"solver"` par :

```json
  "solver": {
    "type": "MONOLITHIC",
    "linear_solver": "BICGSTAB2",
    "precond_type": "BLOCK_TRIANGULAR",
    "use_precond": true,
    "tol": 1e-14,
    "max_iter": 10000
  }
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur**. Comme pour le Tutoriel 1, le mode `"LAPLACIAN"` assemble la matrice globale complète sans décomposition par sous-domaines MPI et s'exécute avec `mpirun -np 1`.

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/grid64_direct output/fortran/grid64_direct

# 2. Exécution avec le moteur C++17 (MPI séquentiel)
mpirun -np 1 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_2_neumann_grid64_direct.json

# 3. Exécution avec le moteur Fortran 2018 (MPI séquentiel)
mpirun -np 1 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_2_neumann_grid64_direct.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Exactitude machine sur polynôme quadratique $\Psi_{\mathrm{ex}}(x,y) = x^2 - y^2$ :**
   À partir de `laplacian_error.csv`, vérifier que pour `DIRECT` les erreurs discrètes $E_{L_1}, E_{L_2}, E_{L_\infty}$ restent au niveau de la précision machine ($< 10^{-12}$) sur toutes les grilles testées ($16^2$ à $256^2$).
2. **Convergence du solveur itératif `BICGSTAB2` et intégrité de la jauge :**
   Vérifier dans `solver_history.csv` que le résidu préconditionné décroît de manière monotone sous `tol = 1e-14`, que l'erreur finale satisfait $E_{L_2} \le 10^{-11}$ sans dérive de constante, et que $\Psi(0,0) = 0$.
3. **Antisymétrie diagonale et parité C++17 / Fortran 2018 :**
   Vérifier à partir de `solution_field_2d.csv` que l'antisymétrie $\Psi(y,x) = -\Psi(x,y)$ est respectée à $< 10^{-12}$ près et que l'écart maximal entre C++17 et Fortran 2018 vérifie $\|\Psi^{\mathrm{C++}} - \Psi^{\mathrm{Fortran}}\|_\infty < 10^{-12}$.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `laplacian_error.csv` : `nx,ny,dx,dy,iters,err_l1,err_l2,err_linf,cpu_time`.
- `solver_history.csv` : historique `iter,residual` des itérations de résolution.
- `solution_field_2d.csv` : grille nodale $(N_x+1)\times(N_y+1)$ avec colonnes `i,j,x,y,psi_num,psi_ex,error_abs`.
- `profile_x_mid.csv` et `profile_y_mid.csv` : profils 1D le long de $y = 0.5$ et $x = 0.5$.
- `solution_field.vtk` : fichier VTK `STRUCTURED_GRID` contenant les champs scalaires `psi_numeric`, `psi_exact`, `error_abs`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 2 (Laplacien Neumann) sur fichiers réels."""
import csv
from pathlib import Path

def verify_neumann_case(dir_cpp: str, dir_ftn: str | None = None) -> None:
    with (Path(dir_cpp) / "laplacian_error.csv").open(newline="") as f:
        err = next(csv.DictReader(f))
    print(f"[C++] nx={err['nx']}, iters={err['iters']}, L1={float(err['err_l1']):.3e}, "
          f"L2={float(err['err_l2']):.3e}, Linf={float(err['err_linf']):.3e}")

    pts = {}
    with (Path(dir_cpp) / "solution_field_2d.csv").open(newline="") as f:
        for r in csv.DictReader(f):
            pts[(int(r["i"]), int(r["j"]))] = float(r["psi_num"])
    antisym = max(abs(pts[(i, j)] + pts[(j, i)]) for i, j in pts if (j, i) in pts)
    print(f"[C++] Antisymétrie max |Psi(x,y)+Psi(y,x)| = {antisym:.3e} | Psi(0,0) = {pts[(0,0)]:.3e}")

    if dir_ftn and (Path(dir_ftn) / "solution_field_2d.csv").exists():
        pts_f = {}
        with (Path(dir_ftn) / "solution_field_2d.csv").open(newline="") as f:
            for r in csv.DictReader(f):
                pts_f[(int(r["i"]), int(r["j"]))] = float(r["psi_num"])
        parity = max(abs(pts[k] - pts_f[k]) for k in pts)
        print(f"[Parité C++/Fortran] Linf = {parity:.3e}")

if __name__ == "__main__":
    verify_neumann_case("output/cpp/grid64_direct", "output/fortran/grid64_direct")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir `output/cpp/grid64_direct/solution_field.vtk` avec **File -> Open** puis **Apply**.
2. Afficher `psi_numeric` en palette divergente (`Cool to Warm`) pour mettre en évidence la selle harmonique $x^2 - y^2$ (positif vers $(1,0)$, négatif vers $(0,1)$, nul sur la diagonale $x=y$).
3. Appliquer le filtre **Contour** sur `psi_numeric` avec des isovaleurs régulières et vérifier l'orthogonalité des isolignes vis-à-vis des frontières à flux nul ($x=0$ et $y=0$).
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
# Tutoriel 4 : Écoulement de Poiseuille Plan Instationnaire (Démarrage Impulsionnel)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_4_poiseuille_instat/` *(validation rapide `Run_Valide-CASES`, cas 4)*.

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
# Tutoriel 5 : Écoulement de Poiseuille Plan à Viscosité Variable Exponentielle

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_5_poiseuille_visc_variable/` *(validation rapide `Run_Valide-CASES`, cas 5)*.

Ce tutoriel explore l'effet d'une hétérogénéité spatiale continue de la viscosité dynamique $\mu(y) = \mu_0 \exp\big(a (y - e)^2 / e^2\big)$ sur un écoulement de Poiseuille. Il valide l'évaluation conservative des contraintes visqueuses aux nœuds de vorticité sur grille MAC ainsi que le comportement des solveurs et préconditionneurs à coefficients variables.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre l'influence d'une viscosité variable spatialement sur le profil d'équilibre de quantité de mouvement.
- Configurer le modèle `"viscosity_model": "EXPONENTIAL"` (`visc_mu0`, `visc_a`, `visc_e`) dans COVOF 2D.
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Quantifier dynamiquement l'ordre de convergence spatiale $\mathcal{O}(\Delta y^2)$ et la linéarité de la contrainte de cisaillement $\tau_{xz}(z) = -G_x z$.
- Visualiser les champs de viscosité, vitesse et cisaillement sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker` (`update_physical_properties`), `NavierStokesSolver`, `BiCGSTAB2Solver`, `DirectSolver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas vérifie la discrétisation conservative du terme de diffusion visqueuse $\nabla \cdot \big[\mu(\mathbf{x})(\nabla \mathbf{u} + \nabla \mathbf{u}^T)\big]$ lorsque la viscosité dynamique varie fortement dans l'espace.

Dans un canal plan $\Omega = [0, L_x] \times [0, L_y]$ avec $L_x = 0.05\,\mathrm{m}$ et $L_y = 2e = 0.02\,\mathrm{m}$ (demi-largeur $e = 0.01\,\mathrm{m}$, coordonnée transversale centrée $z = y - e \in [-e, +e]$), la viscosité dynamique suit une loi exponentielle symétrique :
$$\mu(z) = \mu_0 \exp\left(a \frac{z^2}{e^2}\right), \qquad \text{avec } \mu_0 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}, \quad a = 2.0, \quad e = 0.01\,\mathrm{m}.$$
Le rapport de viscosité entre les parois ($|z|=e$) et l'axe central ($z=0$) vaut $\mu(\pm e)/\mu_0 = \mathrm{e}^a = \mathrm{e}^2 \approx 7.3891$. Sous un forçage axial constant $G_x = \rho g_x = 0.294207317\,\mathrm{Pa/m}$ ($\rho = 1.0\,\mathrm{kg/m}^3$), l'équation d'équilibre stationnaire s'écrit :
$$\frac{\mathrm{d}}{\mathrm{d}z}\left(\mu(z)\frac{\mathrm{d}u}{\mathrm{d}z}\right) = -G_x, \qquad u(-e) = u(+e) = 0.$$

### 2.2 Solution Analytique Exacte
Une première intégration avec la condition de symétrie $\frac{\mathrm{d}u}{\mathrm{d}z}(0) = 0$ fournit la contrainte de cisaillement linéaire :
$$\tau_{xz}(z) = \mu(z)\frac{\mathrm{d}u}{\mathrm{d}z} = -G_x z.$$
En intégrant une seconde fois avec $u(\pm e) = 0$, on obtient le profil analytique gaussien tronqué :
$$u_{\mathrm{ex}}(z) = \frac{G_x e^2}{2 a \mu_0}\left[\exp\left(-a\frac{z^2}{e^2}\right) - \exp(-a)\right],$$
de vitesse maximale au centre $u_{\max} = u_{\mathrm{ex}}(0) = \frac{G_x e^2}{2 a \mu_0}(1 - \mathrm{e}^{-a}) \approx 0.343771\,\mathrm{m/s}$ pour $a = 2.0$. Sur le maillage décalé MAC, les flux de cisaillement aux faces Nord et Sud des volumes de contrôle de $u_{i,j}$ utilisent la viscosité évaluée aux nœuds de vorticité (`mu_node`), assurant une précision formelle d'ordre 2 ($\mathcal{O}(\Delta y^2)$).

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------------- | :----------------------- | :------------------------------------------------------------------ |
| Domaine | $[0, L_x] \times [0, L_y]$ | $L_x = 0.05\,\mathrm{m}$, $L_y = 0.02\,\mathrm{m}$ ($e = 0.01\,\mathrm{m}$) |
| Loi de viscosité | `"viscosity_model": "EXPONENTIAL"` | $\mu_0 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$, $a = 2.0$, $e = 0.01\,\mathrm{m}$, $\rho = 1.0\,\mathrm{kg/m}^3$ |
| Forçage axial | $g_x, g_y$ | $g_x = 0.294207317\,\mathrm{m/s}^2$, $g_y = 0.0$ |
| Séquence de grilles | $N_x \times N_y$ | $N_x = 16$, $N_y \in \{16, 32, 64, 128, 256, 512\}$ |
| Intégration temporelle | $\Delta t$, $t_{\mathrm{final}}$, `steady_state_tol` | $\Delta t = 0.5\,\mathrm{s}$, $t_{\mathrm{final}} = 50.0\,\mathrm{s}$, $\mathrm{steady\_state\_tol} = 10^{-12}$ |
| Solveurs testés | Formulations & Krylov | `AUGMENTED_LAGRANGIAN`, `AAL`, `MONOLITHIC` avec `BICGSTAB2` (`BLOCK_TRIANGULAR`, `SSOR`, `ILU0`) ou `DIRECT` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_5_var_visc_al_ny128.json`
*(Pour Fortran 2018, utiliser une copie dans `case/fortran/` avec `"output_dir": "output/fortran/al_ny128"`.)*

```json
{
  "case_name": "cas_5_var_visc_al_ny128",
  "mode": "NAVIER_STOKES",
  "test_case": "CAS_5_VAR_MU",
  "steady_state_tol": 1e-12,
  "mesh": {
    "nx": 16,
    "ny": 128,
    "lx": 0.05,
    "ly": 0.02
  },
  "physics": {
    "enable_two_phase": false,
    "enable_surface_tension": false,
    "enable_convection": false,
    "viscosity_model": "EXPONENTIAL",
    "visc_a": 2.0,
    "visc_e": 0.01,
    "visc_mu0": 1.85e-05,
    "gx": 0.294207317,
    "gy": 0.0,
    "fluid1": {
      "rho": 1.0,
      "mu": 1.85e-05
    }
  },
  "time": {
    "time_scheme": "EULER",
    "dt": 0.5,
    "t_final": 50.0,
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
    "precond_type": "BLOCK_TRIANGULAR",
    "steady_state_tol": 1e-12
  },
  "output": {
    "output_dir": "output/cpp/al_ny128",
    "case_name": "cas_5_var_visc_al_ny128",
    "vtk_interval": 200,
    "log_interval": 50,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur** (`mpirun -np 4`) pour les solveurs itératifs de Krylov, et `mpirun -np 1` lorsque `"linear_solver": "DIRECT"` est sélectionné.

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/al_ny128 output/fortran/al_ny128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_5_var_visc_al_ny128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_5_var_visc_al_ny128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Convergence spatiale d'ordre 2 face à $u_{\mathrm{ex}}(z)$ :**
   Extraire le profil stationnaire $u(y_j)$ depuis le dernier fichier `.pvtr`/`.vtr` et calculer dynamiquement $E_{L_1}$, $E_{L_2}$, $E_{L_\infty}$ et l'ordre de Richardson sur la séquence $N_y \in \{16, 32, 64, 128, 256, 512\}$. Vérifier $p_{L_1}, p_{L_2}, p_{L_\infty} \in [1.95, 2.05]$.
2. **Linéarité de la contrainte de cisaillement et sensibilité à $a$ (`visc_a`) :**
   Vérifier que $\tau_{xz}(y_{j+1/2}) = \mu(y_{j+1/2})\frac{u_{j+1}-u_j}{\Delta y}$ coïncide avec $-G_x(y_{j+1/2} - e)$. Pour $a \to 0$, vérifier que la solution rejoint le profil parabolique classique à viscosité constante.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier depuis `<case_name>_history.csv` et le dernier `.pvtr` que $|u_{\max}^{\mathrm{C++}} - u_{\max}^{\mathrm{Fortran}}| < 10^{-12}\,\mathrm{m/s}$ et $\|\mathbf{u}^{\mathrm{C++}} - \mathbf{u}^{\mathrm{Fortran}}\|_\infty < 10^{-12}\,\mathrm{m/s}$.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` : `step,time,dt,center_of_mass_y,rise_velocity,bubble_volume,kinetic_energy,max_u,max_div` (`9` colonnes C++ / `11` colonnes Fortran avec `avg_rho_gas,avg_p_gas`).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : contiennent `PhaseIndicator`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 5 (Viscosité variable exponentielle) depuis .vtr."""
import math
import xml.etree.ElementTree as ET
from pathlib import Path

def load_pvtr_velocity_x(pvtr_path: str):
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    ys_all, ux_all = [], []
    for piece in root.findall(".//Piece"):
        vtr = ET.parse(pvtr.parent / piece.attrib["Source"]).getroot()
        ext = [int(v) for v in piece.attrib["Extent"].split()]
        nx_l, ny_l = ext[1] - ext[0], ext[3] - ext[2]
        y_coords = [float(v) for v in vtr.find(".//DataArray[@Name='YCoordinates']").text.split()]
        vel = [float(v) for v in vtr.find(".//DataArray[@Name='Velocity']").text.split()]
        for j in range(ny_l):
            yc = 0.5 * (y_coords[j] + y_coords[j + 1])
            u_row = [vel[3 * (j * nx_l + i)] for i in range(nx_l)]
            ys_all.append(yc)
            ux_all.append(sum(u_row) / nx_l)
    return ys_all, ux_all

def check_exact(pvtr_path: str, gx: float = 0.294207317, mu0: float = 1.85e-5, a: float = 2.0, e: float = 0.01):
    ys, u_num = load_pvtr_velocity_x(pvtr_path)
    errs = []
    for y, un in zip(ys, u_num):
        z = y - e
        u_ex = (gx * e * e / (2.0 * a * mu0)) * (math.exp(-a * (z / e) ** 2) - math.exp(-a))
        errs.append(abs(un - u_ex))
    l2 = math.sqrt(sum(v * v for v in errs) / len(errs))
    print(f"[{pvtr_path}] N_samples={len(errs)} | L1={sum(errs)/len(errs):.4e} | L2={l2:.4e} | Linf={max(errs):.4e}")

if __name__ == "__main__":
    p = sorted(Path("output/cpp/al_ny128").glob("*.pvtr"))[-1]
    check_exact(str(p))
```

### 6.4 Visualisation sous ParaView
1. Ouvrir le dernier fichier `<case_name>_<step:06d>.pvtr` dans ParaView (**File -> Open**).
2. Inspecter le champ `Viscosity` pour visualiser le gradient exponentiel symétrique entre le cœur du canal ($\mu_0$) et les parois ($\mu_0 \mathrm{e}^a$).
3. Appliquer le filtre **Plot Over Line** transversalement en $x = 0.025\,\mathrm{m}$ sur `Velocity` (composante $X$) pour observer l'aplatissement au centre et le resserrement du profil gaussien tronqué près des parois.
# Tutoriel 6 : Goutte Statique en Équilibre Capillaire et Loi de Young-Laplace

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_6_laplace_droplet/` *(validation rapide `Run_Valide-CASES`, cas 6)*.

Ce tutoriel constitue le banc d'essai de référence de l'hydrostatique diphasique capillaire. Il valide le modèle de force superficielle continue CSF (*Continuum Surface Force*) à équilibre de forces (*Balanced-Force*), les lissages d'interface (diffusif et SVOF de Helmholtz), ainsi que l'option hybride par marqueurs lagrangiens CVI pour la maîtrise des courants parasites autour d'une goutte circulaire à l'équilibre.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre le modèle CSF et la loi hydrostatique de Young-Laplace sur maillage décalé MAC.
- Respecter la contrainte de pas de temps capillaire explicite de Brackbill lors du raffinement du maillage.
- Configurer un fichier `.json` conforme avec spécification explicite d'une goutte circulaire (`axis_a`, `axis_b`, `pert_amp`).
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Quantifier dynamiquement le saut de pression $\Delta P$, l'intensité des courants parasites $\|\mathbf{u}\|_\infty$ et la conservation du volume depuis les fichiers `.vtr`.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `InterfaceMarkers` (optionnel CVI), `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_markers` (optionnel CVI), `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide le couplage hydro-capillaire entre les équations de Navier-Stokes diphasiques et la force de tension superficielle CSF discrétisée aux faces des vitesses exactement aux mêmes points que $-\nabla p$.

Dans un domaine carré fermé $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$ en apesanteur ($\mathbf{g} = \mathbf{0}$) avec conditions de glissement libre (`FREE_SLIP`) sur les 4 parois, une goutte circulaire de rayon $R_0 = 0.25\,\mathrm{m}$ ($C = 1$) est centrée en $(x_c, y_c) = (0.5, 0.5)\,\mathrm{m}$ au repos initial $\mathbf{u}(\mathbf{x}, 0) = \mathbf{0}$. Les deux fluides ont des propriétés identiques $\rho_1 = \rho_2 = 1.0\,\mathrm{kg/m}^3$, $\mu_1 = \mu_2 = 1.0\,\mathrm{Pa\cdot s}$, et un coefficient de tension superficielle $\sigma = 1.0\,\mathrm{N/m}$.

### 2.2 Loi Analytique de Young-Laplace et Contrainte de Stabilité
À l'équilibre exact, la vitesse est identiquement nulle et la courbure géométrique vaut $\kappa_{\mathrm{ex}} = 1/R_0 = 4.0000\,\mathrm{m}^{-1}$. Le saut de pression à travers l'interface vaut :
$$\mathbf{u}_{\mathrm{ex}}(\mathbf{x}, t) \equiv \mathbf{0}, \qquad \Delta P_{\mathrm{ex}} = P_{\mathrm{int}} - P_{\mathrm{ext}} = \frac{\sigma}{R_0} = \frac{1.0}{0.25} = 4.0000\,\mathrm{Pa}.$$

Deux mécanismes numériques gouvernent la précision et la stabilité :
1. **Contrainte de pas de temps capillaire (Brackbill) :** La force capillaire $\mathbf{F}_\sigma = \sigma \kappa \nabla C$ étant traitée explicitement dans le second membre, le pas de temps doit respecter :
   $$\Delta t < \Delta t_{\mathrm{cap}} = \sqrt{\frac{(\rho_1 + \rho_2)\,\Delta x^3}{4\pi \sigma}}.$$
2. **Régularisation et marqueurs :** L'activation de `enable_diffusion_smoothing` (ou de `enable_svof_smoothing`) régularise le calcul de la courbure $\kappa = -\nabla \cdot (\nabla \widetilde{C}/\|\nabla \widetilde{C}\|)$, tandis que l'option `use_markers: true` (`CVI_RK2`) calcule une courbure lagrangienne régulière sans diffusion de l'interface.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et Goutte | $\Omega$, $(x_c, y_c)$, $R_0$ | $[0,1]\times[0,1]\,\mathrm{m}^2$, $(0.5, 0.5)\,\mathrm{m}$, $R_0 = 0.25\,\mathrm{m}$ (`axis_a = axis_b = 0.25`, `pert_amp = 0.0`) |
| Propriétés physiques | $\rho_1=\rho_2, \mu_1=\mu_2, \sigma, \mathbf{g}$ | $\rho = 1.0\,\mathrm{kg/m}^3$, $\mu = 1.0\,\mathrm{Pa\cdot s}$, $\sigma = 1.0\,\mathrm{N/m}$, $\mathbf{g} = \mathbf{0}$ |
| Séquence de grilles et $\Delta t \le 0.5\Delta t_{\mathrm{cap}}$ | $N \in \{32, 64, 128, 256\}$, $t_{\mathrm{final}}=0.02\,\mathrm{s}$ | $\Delta t \in \{10^{-3}, 3\times 10^{-4}, 10^{-4}, 4\times 10^{-5}\}\,\mathrm{s}$ |
| Schémas VOF et régularisation | `vof` & `physics` | `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC`, `VOF_TVD` ; `enable_diffusion_smoothing`, `enable_svof_smoothing`, `use_markers` |
| Couplage & Solveurs | `solver.type`, `linear_solver`, `precond_type` | `AUGMENTED_LAGRANGIAN`, `AAL`, `MONOLITHIC`, `STANDARD_PROJECTION` ; `BICGSTAB2`, `GMRES`, `MINRES` ; `BLOCK_TRIANGULAR`, `JACOBI`, `SSOR`, `ILU0` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_6_wy_plic_n128.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/wy_plic_n128"`.)*

```json
{
  "case_name": "cas_6_wy_plic_n128",
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 1.0,
    "ly": 1.0
  },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": true,
    "enable_convection": false,
    "enable_svof_smoothing": false,
    "svof_filter_width": 0.5,
    "enable_diffusion_smoothing": true,
    "diffusion_smoothing_steps": 3,
    "sigma": 1.0,
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": { "rho": 1.0, "mu": 1.0 },
    "fluid2": { "rho": 1.0, "mu": 1.0 }
  },
  "initial_condition": {
    "type": "DROPLET",
    "center_x": 0.5,
    "center_y": 0.5,
    "radius": 0.25,
    "axis_a": 0.25,
    "axis_b": 0.25,
    "pert_amp": 0.0
  },
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.0001,
    "t_final": 0.02,
    "cfl_max": 0.3,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "top":    { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 400,
    "r_al": 100.0,
    "uzawa_iter": 1,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.5
  },
  "output": {
    "output_dir": "output/cpp/wy_plic_n128",
    "case_name": "cas_6_wy_plic_n128",
    "vtk_interval": 20,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur** (`mpirun -np 4` avec les solveurs de Krylov, ou `mpirun -np 1` si `STANDARD_PROJECTION` ou `DIRECT` est choisi).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/wy_plic_n128 output/fortran/wy_plic_n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_6_wy_plic_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_6_wy_plic_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Saut de pression de Young-Laplace et convergence en maillage :**
   À $t = t_{\mathrm{final}}$, extraire dynamiquement depuis les fichiers `.vtr` la pression moyenne interne $P_{\mathrm{int}}$ sur $\{r < 0.15\,\mathrm{m},\, C > 0.99\}$ et externe $P_{\mathrm{ext}}$ sur $\{r > 0.38\,\mathrm{m},\, C < 0.01\}$, puis évaluer :
   $$\Delta P_{\mathrm{num}} = P_{\mathrm{int}} - P_{\mathrm{ext}}, \qquad E_{\mathrm{rel}}(\Delta P) = \frac{|\Delta P_{\mathrm{num}} - \sigma/R_0|}{\sigma/R_0}.$$
   Vérifier la décroissance de $E_{\mathrm{rel}}(\Delta P)$ avec le raffinement du maillage ($N = 32 \to 256$).
2. **Maîtrise des courants parasites et conservation du volume :**
   Calculer à $t_{\mathrm{final}}$ le nombre capillaire parasite $\mathrm{Ca}_{\mathrm{parasite}} = \mu \|\mathbf{u}\|_\infty / \sigma$ ainsi que la variation de volume $|\sum C_{i,j}(t_{\mathrm{final}}) - \sum C_{i,j}(0)| / \sum C_{i,j}(0)$ directement à partir de `PhaseIndicator` dans les `.vtr`. Comparer les variantes sans lissage, avec lissage diffusif, avec lissage SVOF et avec marqueurs CVI (`marker_curvature: true`).
3. **Invariance des solveurs, préconditionnement MPI et parité C++17 / Fortran 2018 :**
   Vérifier l'absence de dérive de pression moyenne (projection globale à moyenne nulle dans `BlockTriangularPreconditioner` en MPI pur) et comparer $\Delta P_{\mathrm{num}}$ ainsi que $\|\mathbf{u}\|_\infty$ entre `covof2d_cpp` et `covof2d_fortran` *(note : utiliser `BICGSTAB2` ou `GMRES` avec `BLOCK_TRIANGULAR`, et réserver `MINRES` aux formulations strictement symétriques)*.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` : `step,time,dt,center_of_mass_y,rise_velocity,bubble_volume,kinetic_energy,max_u,max_div` (`9` colonnes C++ / `11` colonnes Fortran avec `avg_rho_gas,avg_p_gas`).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `PhaseIndicator_Smoothed` (si SVOF activé), `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Curvature`, `Velocity`.
- `<case_name>_markers_<step:06d>.vtk` (si `use_markers: true` et `write_markers_vtk: true`).

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 6 (Young-Laplace) sur fichiers .pvtr/.vtr réels."""
import math
import xml.etree.ElementTree as ET
from pathlib import Path

def analyze_young_laplace(pvtr_path: str, sigma: float = 1.0, r0: float = 0.25) -> None:
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    p_in, p_out, u_max, vol_drop = [], [], 0.0, 0.0

    for piece in root.findall(".//Piece"):
        vtr = ET.parse(pvtr.parent / piece.attrib["Source"]).getroot()
        ext = [int(v) for v in piece.attrib["Extent"].split()]
        nx_l, ny_l = ext[1] - ext[0], ext[3] - ext[2]
        xc = [float(v) for v in vtr.find(".//DataArray[@Name='XCoordinates']").text.split()]
        yc = [float(v) for v in vtr.find(".//DataArray[@Name='YCoordinates']").text.split()]
        C = [float(v) for v in vtr.find(".//DataArray[@Name='PhaseIndicator']").text.split()]
        P = [float(v) for v in vtr.find(".//DataArray[@Name='Pressure']").text.split()]
        V = [float(v) for v in vtr.find(".//DataArray[@Name='Velocity']").text.split()]
        for j in range(ny_l):
            y = 0.5 * (yc[j] + yc[j + 1])
            dy = yc[j + 1] - yc[j]
            for i in range(nx_l):
                x = 0.5 * (xc[i] + xc[i + 1])
                dx = xc[i + 1] - xc[i]
                k = j * nx_l + i
                vol_drop += C[k] * dx * dy
                r = math.hypot(x - 0.5, y - 0.5)
                if r < 0.15 and C[k] > 0.99:
                    p_in.append(P[k])
                elif r > 0.38 and C[k] < 0.01:
                    p_out.append(P[k])
                spd = math.hypot(V[3 * k], V[3 * k + 1])
                if spd > u_max:
                    u_max = spd

    dp_num = (sum(p_in) / len(p_in)) - (sum(p_out) / len(p_out))
    dp_ex = sigma / r0
    print(f"[{pvtr.name}] dP_num={dp_num:.6f} Pa (ex={dp_ex:.6f} Pa, err={abs(dp_num-dp_ex)/dp_ex:.3e}) "
          f"| max|u|={u_max:.4e} m/s | V_drop={vol_drop:.6f} m^2")

if __name__ == "__main__":
    last_pvtr = sorted(Path("output/cpp/wy_plic_n128").glob("*.pvtr"))[-1]
    analyze_young_laplace(str(last_pvtr))
```

### 6.4 Visualisation sous ParaView
1. Ouvrir le dernier fichier `<case_name>_<step:06d>.pvtr` dans ParaView (**File -> Open**).
2. Tracer l'isocontour `PhaseIndicator = 0.5` avec le filtre **Contour**, puis afficher `Pressure` en fond pour observer le saut de pression uniforme à l'intérieur de la goutte.
3. Appliquer le filtre **Plot Over Line** le long de $y = 0.5\,\mathrm{m}$ sur `Pressure` pour vérifier la raideur de la marche de Young-Laplace et la planéité des plateaux interne et externe.
4. Superposer des vecteurs **Glyph** sur `Velocity` pour localiser les éventuels courants parasites résiduels près de l'interface.
# Tutoriel 7 : Advection d'une Tâche Circulaire dans un Tourbillon Cisaillé

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_7_vortex_shear/` *(validation rapide `Run_Valide-CASES`, cas 7)*.

Ce tutoriel isole le transport cinématique pur de l'interface fluide VOF sans résolution hydrodynamique (mode `"VORTEX_ADVECTION"`). Une tâche circulaire est soumise à un champ de vitesse cisaillé instationnaire réversible qui l'étire en une fine spirale avant d'inverser le flux à $t = T/2$ pour reformer le disque initial à $t = T$.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les schémas d'advection VOF géométriques (`WEYMOUTH_YUE`, `CIAM`) et algébriques (`VOF_TVD`), ainsi que la règle de cohérence avec les méthodes de reconstruction et les limiteurs TVD.
- Configurer le mode cinématique `"VORTEX_ADVECTION"` et l'option marqueurs lagrangiens CVI.
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Quantifier dynamiquement l'erreur géométrique de retour $E_{L_1}$, $E_{L_2}$, $E_{L_\infty}$ et la conservation de masse depuis `vortex_error.csv`.
- Visualiser l'étirement maximal à $t = T/2$ et le retour à $t = T$ sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `InterfaceMarkers` (optionnel CVI), `VTKWriter`, `Simulation` (`Simulation::run_vortex_advection`).
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_markers` (optionnel CVI), `module_vtk_writer`, `module_simulation` (`sim_run_vortex_advection`).

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas cinématique (benchmark de Rider & Kothe) valide les schémas de transport et de reconstruction d'interface Volume-of-Fluid (VOF) sous étirement tourbillonnaire sévère et réversibilité temporelle exacte.

Dans la cavité carrée unitaire $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$, un disque fluide initial ($C = 1$) de rayon $R_0 = 0.125\,\mathrm{m}$ est centré en $(x_c, y_c) = (0.5, 0.75)\,\mathrm{m}$ à $t = 0$. L'indicatrice de phase $C(\mathbf{x}, t) \in [0, 1]$ obéit à l'équation d'advection pure sans terme source :
$$\frac{\partial C}{\partial t} + \nabla \cdot (C \mathbf{u}) = 0.$$
Le champ de vitesse solénoïdal ($\nabla \cdot \mathbf{u} \equiv 0$) dérive de la fonction de courant instationnaire modulée par un cosinus temporel de période $T = 6.0\,\mathrm{s}$ :
$$\Psi(x, y, t) = \frac{1}{\pi} \sin^2(\pi x) \sin^2(\pi y) \cos\left(\frac{\pi t}{T}\right),$$
évalué au demi-pas de temps $t + \Delta t/2$ aux faces des mailles décalées MAC pour assurer une précision temporelle d'ordre 2 :
$$u(x, y, t) = -\frac{\partial \Psi}{\partial y} = -\sin^2(\pi x) \sin(2\pi y) \cos\left(\frac{\pi t}{T}\right),$$
$$v(x, y, t) = +\frac{\partial \Psi}{\partial x} = +\sin(2\pi x) \sin^2(\pi y) \cos\left(\frac{\pi t}{T}\right).$$

### 2.2 Réversibilité Exacte et Règle de Cohérence VOF
De $t = 0$ à $t = T/2 = 3.0\,\mathrm{s}$, le disque est étiré en une fine spirale ; à $t = 3.0\,\mathrm{s}$, le champ de vitesse change de signe et réenroule le filament jusqu'à $t = T = 6.0\,\mathrm{s}$ où la solution exacte coïncide identiquement avec l'état initial :
$$C_{\mathrm{ex}}(\mathbf{x}, T) \equiv C(\mathbf{x}, 0), \qquad M(T) = \int_\Omega C(\mathbf{x}, T)\,\mathrm{d}\Omega = M(0).$$
Dans `Config.cpp` et `module_config.f90`, une règle stricte de cohérence s'applique :
- Si `advection_method` vaut `"VOF_TVD"`, `reconstruction_method` est forcé à `"NONE"` et un limiteur (`"SUPERBEE"`, `"VAN_LEER"`, `"MUSCL"` ou `"HRIC"`) est utilisé.
- Si `advection_method` vaut `"WEYMOUTH_YUE"` ou `"CIAM"`, `tvd_limiter` est forcé à `"NONE"` et la reconstruction géométrique (`"PLIC"`, `"ELVIRA"`, `"LVIRA"` ou `"SLIC"`) est activée.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et Disque initial | $\Omega$, $(x_c, y_c)$, $R_0$ | $[0,1]\times[0,1]\,\mathrm{m}^2$, $(0.5, 0.75)\,\mathrm{m}$, $R_0 = 0.125\,\mathrm{m}$ (`axis_a = axis_b = 0.125`, `pert_amp = 0.0`) |
| Période cinématique | $T = t_{\mathrm{final}}$ | $6.0\,\mathrm{s}$ (étirement maximal à $T/2 = 3.0\,\mathrm{s}$) |
| Séquence de grilles et $\Delta t$ | $(N\times N, \Delta t)$ | $(64^2, 0.005\,\mathrm{s})$, $(128^2, 0.0025\,\mathrm{s})$, $(256^2, 0.00125\,\mathrm{s})$, $(512^2, 0.0005\,\mathrm{s})$ |
| Méthodes d'advection VOF | `vof.advection_method` | `"WEYMOUTH_YUE"`, `"CIAM"`, `"VOF_TVD"` |
| Reconstructions et limiteurs | `reconstruction_method` / `tvd_limiter` | `"PLIC"`, `"ELVIRA"`, `"LVIRA"`, `"SLIC"` ; `"SUPERBEE"`, `"VAN_LEER"`, `"MUSCL"`, `"HRIC"` |
| Option marqueurs lagrangiens | `vof.use_markers` | `false` (VOF pur) ou `true` (`"CVI_RK2"`, `markers_per_cell: 4`, `write_markers_vtk: true`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_7_weymouth_yue_plic_n128.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/wy_plic_n128"`.)*

```json
{
  "case_name": "cas_7_weymouth_yue_plic_n128",
  "mode": "VORTEX_ADVECTION",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 1.0,
    "ly": 1.0
  },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": false,
    "enable_convection": true,
    "sigma": 0.0,
    "fluid1": { "rho": 1000.0, "mu": 0.001 },
    "fluid2": { "rho": 1.0,    "mu": 1.0e-05 }
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.5
  },
  "initial_condition": {
    "type": "VORTEX_SHEAR",
    "center_x": 0.5,
    "center_y": 0.75,
    "radius": 0.125,
    "axis_a": 0.125,
    "axis_b": 0.125,
    "pert_amp": 0.0
  },
  "time": {
    "dt": 0.0025,
    "t_final": 6.0,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "right":  { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-12,
    "max_iter": 500,
    "r_al": 100.0,
    "uzawa_iter": 1,
    "use_precond": true,
    "precond_type": "JACOBI"
  },
  "output": {
    "output_dir": "output/cpp/wy_plic_n128",
    "case_name": "cas_7_weymouth_yue_plic_n128",
    "vtk_interval": 240,
    "log_interval": 200,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Le mode `"VORTEX_ADVECTION"` est entièrement parallélisé via **MPI pur** (échange des halos VOF et réductions globales des erreurs $L_1, L_2, L_\infty$ et des masses).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/wy_plic_n128 output/fortran/wy_plic_n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_7_weymouth_yue_plic_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_7_weymouth_yue_plic_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Erreur géométrique de retour à $t = T = 6.0\,\mathrm{s}$ :**
   Extraire directement depuis `vortex_error.csv` les normes calculées par les binaires :
   $$E_{L_1} = \frac{1}{N_x N_y}\sum_{i=1}^{N_x}\sum_{j=1}^{N_y} |C_{i,j}(T) - C_{i,j}(0)|, \qquad E_{L_2} = \sqrt{\frac{1}{N_x N_y}\sum_{i=1}^{N_x}\sum_{j=1}^{N_y} |C_{i,j}(T) - C_{i,j}(0)|^2}.$$
   Calculer l'ordre de convergence spatiale $p = \ln(E_{L_1}(2h)/E_{L_1}(h)) / \ln(2)$ entre grilles successives ($p \in [1.0, 2.0]$ attendu selon la résolution de la queue du filament).
2. **Conservation de la masse et bornitude :**
   Vérifier dans `vortex_error.csv` que `mass_conservation = mass_final / mass_init` reste à $1.00000 \pm 10^{-5}$ pour `WEYMOUTH_YUE` et `VOF_TVD`, et que $0 \le C_{i,j}(t) \le 1$ dans chaque fichier `.vtr`.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier que $|E_{L_1}^{\mathrm{C++}} - E_{L_1}^{\mathrm{Fortran}}| < 10^{-12}$ et que les erreurs $E_{L_2}$ et `mass_conservation` concordent.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `vortex_error.csv` : `nx,ny,dt,T,advection,reconstruction,limiter,err_l1,err_l2,err_linf,mass_init,mass_final,mass_conservation`.
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` (exportés à $t=0$, à $t \approx T/2$, tous les `vtk_interval` pas, et à $t = T$).
- `<case_name>_markers_<step:06d>.vtk` (si `use_markers: true` et `write_markers_vtk: true`).

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 7 (Vortex cisaillé) sur vortex_error.csv."""
import csv
import math
from pathlib import Path

def read_vortex_error(out_dir: str) -> dict:
    with (Path(out_dir) / "vortex_error.csv").open(newline="") as f:
        return next(csv.DictReader(f))

if __name__ == "__main__":
    prev = None
    for n in [64, 128, 256, 512]:
        d = f"output/cpp/wy_plic_n{n}"
        if not (Path(d) / "vortex_error.csv").exists():
            continue
        r = read_vortex_error(d)
        l1 = float(r["err_l1"])
        mc = float(r["mass_conservation"]) * 100.0
        order = math.log(float(prev["err_l1"]) / l1) / math.log(n / int(prev["nx"])) if prev else float("nan")
        print(f"N={n:3d} | {r['advection']}+{r['reconstruction']} | L1={l1:.4e} "
              f"| L2={float(r['err_l2']):.4e} | Masse={mc:.5f}% | p={order:.2f}")
        prev = r
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série temporelle `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Afficher `PhaseIndicator` et appliquer **Contour** sur la valeur $0.5$ aux instants clés $t = 0$, $t = T/2 = 3.0\,\mathrm{s}$ (étirement maximal en spirale) et $t = T = 6.0\,\mathrm{s}$.
3. Comparer visuellement la superposition de l'isocontour final $t = T$ avec l'isocontour initial $t = 0$ (et, si `use_markers: true`, charger les fichiers `*_markers_*.vtk` avec rendu `Point Gaussian`).
# Tutoriel 8 : Écoulement de la Cavité Carrée Entraînée 2D à $Re = 1000$

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_8_cavite_entrainee/` *(validation rapide `Run_Valide-CASES`, cas 8)*.

Ce tutoriel valide la précision spatiale d'ordre 2 des opérateurs convectifs et visqueux incompressibles de COVOF 2D sur le benchmark classique de la cavité carrée entraînée (*Lid-Driven Cavity*) à $Re = 1000$. Il permet de comparer les schémas spatiaux (`CENTERED`, `TVD_SUPERBEE`, `UPWIND`), les linéarisations (`NEWTON`, `ADAMS_BASHFORTH`, `PICARD`) et les couplages vitesse-pression (`AL`, `AAL`, `MONOLITHIC`).

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Configurer un écoulement tourbillonnaire confiné avec cisaillement pariétal entraîné sur la paroi supérieure.
- Comparer quantitativement la dissipation numérique des schémas convectifs `CENTERED`, `TVD_SUPERBEE` et `UPWIND` câblés dans `apply_system_matrix`.
- Évaluer dynamiquement l'ordre de convergence spatiale de Richardson sur la fonction de courant $\psi(x, y)$, la vorticité $\omega_z$ et les profils médians $u(0.5, y)$, $v(x, 0.5)$.
- Exécuter le cas en parallèle MPI pur (`covof2d_cpp` et `covof2d_fortran`) et exploiter les sorties `.pvtr`/`.vtr`.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Références Bibliographiques

Sur le domaine carré $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$, un fluide newtonien incompressible ($\rho = 1.0\,\mathrm{kg/m}^3$, $\mu = 10^{-3}\,\mathrm{Pa\cdot s}$) est entraîné par la paroi supérieure mobile à vitesse tangentielle $u_{\mathrm{top}} = -1.0\,\mathrm{m/s}$, les parois gauche, droite et inférieure étant fixes et adhérentes ($\mathbf{u} = \mathbf{0}$) :
$$Re = \frac{\rho |u_{\mathrm{top}}| L}{\mu} = 1000.$$

La littérature de référence (Botella & Peyret 1998 ; Ghia et al. 1982, adaptée au sens $u_{\mathrm{top}} = -1$) fournit les repères suivants :
- **Tourbillon primaire central :** $\psi_{\max} \approx 0.1189366$, $\omega_z \approx 2.067753$ en $(0.4692, 0.5652)$.
- **Tourbillon secondaire inférieur gauche (aval) :** $\psi_{\min,1} \approx -1.729717\times 10^{-3}$, $\omega_z \approx -1.109789$ en $(0.1360, 0.1118)$.
- **Tourbillon secondaire inférieur droit (amont) :** $\psi_{\min,2} \approx -2.334528\times 10^{-4}$, $\omega_z \approx -0.35214$ en $(0.9167, 0.0781)$.
- **Extrema des profils médians :** $\max_y u(0.5, y) \approx 0.3885698$ vers $y \approx 0.1717$ ; sur $v(x, 0.5)$, $v_{\max} \approx 0.5270771$ vers $x \approx 0.0909$ et $v_{\min} \approx -0.3769447$ vers $x \approx 0.8422$.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et maillages | `mesh.{nx,ny,lx,ly}` | $L_x = L_y = 1.0\,\mathrm{m}$, grilles $32^2, 64^2, 128^2, 256^2, 512^2$ |
| Fluide et convection | `physics` | `enable_two_phase: false`, $\rho=1.0$, $\mu=10^{-3}$, `CENTERED` / `TVD_SUPERBEE` / `UPWIND`, `NEWTON` / `AB2` / `PICARD` |
| Temps | `time.{time_scheme,dt,t_final}` | `"GEAR"`, $\Delta t = 0.05\,\mathrm{s}$ ($0.02\,\mathrm{s}$ sur $256^2, 512^2$), $t_{\mathrm{final}} = 50.0\,\mathrm{s}$ |
| Conditions aux limites | `boundary_conditions` | `top`: `DIRICHLET` ($u = -1.0, v = 0$) ; `left`, `right`, `bottom`: `DIRICHLET` ($u=v=0$) |
| Solveurs et précond. | `solver` | `AL` ($r_{\mathrm{al}}=100$), `AAL`, `MONOLITHIC` ; `BICGSTAB2` / `GMRES` ; `BLOCK_TRIANGULAR` ($\mathrm{tol}=10^{-9}$, `max_iter=800`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_8_cavite_n128.json`
*(Pour Fortran 2018, utiliser une copie dans `case/fortran/` avec `"output_dir": "output/fortran/n128"`.)*

```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 1.0,
    "ly": 1.0
  },
  "physics": {
    "enable_two_phase": false,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "NEWTON",
    "fluid1": { "rho": 1.0, "mu": 0.001 },
    "gx": 0.0,
    "gy": 0.0
  },
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.05,
    "t_final": 50.0,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u":  0.0, "v": 0.0 },
    "right":  { "type": "DIRICHLET", "u":  0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u":  0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": -1.0, "v": 0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "r_al": 100.0,
    "uzawa_iter": 3,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-09,
    "max_iter": 800,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "output": {
    "output_dir": "output/cpp/n128",
    "case_name": "cas_8_cavite_n128",
    "vtk_interval": 200,
    "log_interval": 20,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Exécuter impérativement avec `cwd` placé dans le répertoire du cas de validation et répertoires de sortie relatifs courts pré-créés :

```bash
# 1. Créer les répertoires de sortie
mkdir -p output/cpp/n128 output/fortran/n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_8_cavite_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_8_cavite_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification Dynamique
1. **Fonction de courant et vorticité :** Intégrer $\partial_y \psi = u$ (avec $\psi = 0$ sur la paroi inférieure) et calculer $\omega_z = \partial_x v - \partial_y u$ à partir de `Velocity` lu dans le dernier `.pvtr`/`.vtr`.
2. **Convergence spatiale d'ordre 2 :** Calculer les ordres de Richardson $p$ entre $32^2 \to 64^2 \to 128^2 \to 256^2 \to 512^2$ sur les extrema $\psi_{\max}$, $\psi_{\min,1}$ et les profils médians.
3. **Invariance et comparaison des schémas :** Quantifier l'écart entre `CENTERED`, `TVD_SUPERBEE` et `UPWIND` (notamment sur la capture des tourbillons secondaires de coin) et la parité C++17 / Fortran 2018.

### 6.2 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 8 (Cavité entraînée Re=1000) depuis .pvtr/.vtr."""
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

def read_pvtr_velocity(pvtr_path: str):
    p = Path(pvtr_path)
    root = ET.parse(p).getroot().find("PRectilinearGrid")
    ext = [int(x) for x in root.attrib["WholeExtent"].split()]
    nx, ny = ext[1] - ext[0], ext[3] - ext[2]
    u = np.zeros((ny, nx))
    v = np.zeros((ny, nx))
    for piece in root.findall("Piece"):
        pe = [int(x) for x in piece.attrib["Extent"].split()]
        vtr = ET.parse(p.parent / piece.attrib["Source"]).getroot().find("RectilinearGrid").find("Piece")
        for da in vtr.find("CellData").findall("DataArray"):
            if da.attrib.get("Name") == "Velocity":
                arr = np.fromstring(da.text, sep=" ").reshape((pe[3]-pe[2], pe[1]-pe[0], 3))
                u[pe[2]:pe[3], pe[0]:pe[1]] = arr[:, :, 0]
                v[pe[2]:pe[3], pe[0]:pe[1]] = arr[:, :, 1]
    return u, v

if __name__ == "__main__":
    u, v = read_pvtr_velocity("output/cpp/n128/cas_8_cavite_n128_001000.pvtr")
    ny, nx = u.shape
    dy = 1.0 / ny
    psi = np.cumsum(u, axis=0) * dy
    print(f"psi_max = {psi.max():.6f} | psi_min = {psi.min():.6e} | max|u(0.5,y)| = {np.max(np.abs(u[:, nx//2])):.6f}")
```

### 6.3 Visualisation sous ParaView
1. Ouvrir `output/cpp/n128/cas_8_cavite_n128_001000.pvtr` avec **File -> Open**.
2. Tracer les lignes de courant avec le filtre **Stream Tracer** (semences diagonales) et calculer la vorticité via **Compute Derivatives** pour observer le tourbillon central et les deux tourbillons secondaires de coin.
# Tutoriel 9 : Écoulement Laminaire derrière une Marche Descendante ($Re = 100$ à $1000$)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_9_marche_descendante/` *(validation rapide `Run_Valide-CASES`, cas 9)*.

Ce tutoriel valide la pénalisation volumique de Darcy-Brinkman pour obstacle solide rectangulaire (`type: "BOX"`), la condition d'entrée parabolique décalée `POISEUILLE_STEP`, la sortie ouverte `NEUMANN` et l'évolution des zones de recirculation en fonction du nombre de Reynolds.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Modéliser une géométrie à marche descendante (*Backward-Facing Step*) sur maillage cartésien par pénalisation de Brinkman (`obstacles` avec $K = 10^{-14}\,\mathrm{m}^2$).
- Configurer le profil d'entrée dédié `"profile": "POISEUILLE_STEP"` (`y_step`, `u_mean`).
- Mesurer dynamiquement la longueur de recollement primaire $X_r/H$ sur la paroi inférieure et l'apparition du décollement secondaire supérieur $(X_s/H, X_{rs}/H)$ pour $Re \ge 400$.
- Comparer les schémas convectifs `CENTERED`, `TVD_SUPERBEE` et `UPWIND` en parallèle MPI pur.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `NavierStokesSolver` (`is_*_in_obstacle`), `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Grandeurs de Référence

Dans le canal $\Omega = [0, 40.0] \times [0, 2.0612244898]\,\mathrm{m}^2$, la marche solide occupe $\Omega_{\mathrm{step}} = [0, 4.0] \times [0, 1.0]\,\mathrm{m}^2$ (hauteur de marche $H = 1.0\,\mathrm{m}$, rapport d'expansion $L_y / (L_y - H) \approx 1.942$). En entrée ($x=0$), le fluide ($\rho = 1.0\,\mathrm{kg/m}^3$, $\bar{U} = 1.0\,\mathrm{m/s}$) entre avec le profil établi sur $y \in [H, L_y]$ :
$$u(0, y) = 6\,\bar{U}\,\frac{(y - H)(L_y - y)}{(L_y - H)^2} \quad (y \in [H, L_y]), \qquad u(0, y) = 0 \quad (y \in [0, H]).$$

Le nombre de Reynolds est défini par $Re = \frac{2 \rho \bar{U} H}{\mu} \in [100, 1000]$, soit $\mu = \frac{2}{Re}\,\mathrm{Pa\cdot s}$. Les grandeurs extraites par changement de signe du cisaillement pariétal $\tau_w(x)$ sont comparées aux résultats de référence d'Armaly et al. (1983), Gartling (1990) et Erturk (2008).

---

## 3. Paramètres Physiques, Numériques et Discrétisation

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et maillages | `mesh.{nx,ny,lx,ly}` | $L_x = 40.0\,\mathrm{m}$, $L_y = 2.0612244898\,\mathrm{m}$ ; $128\times 32$, $256\times 64$, $512\times 128$ |
| Marche immergée | `obstacles[0]` | `{"type": "BOX", "x0": 0.0, "y0": 0.0, "x1": 4.0, "y1": 1.0, "permeability": 1e-14}` |
| Fluide et convection | `physics` | `enable_two_phase: false`, $\rho=1.0$, $\mu = 2/Re$, `CENTERED`/`TVD_SUPERBEE`/`UPWIND`, `NEWTON`/`AB2`/`PICARD` |
| Temps | `time.{time_scheme,dt,t_final}` | `"GEAR"`, $\Delta t = 0.05\,\mathrm{s}$, $t_{\mathrm{final}} = 100.0\,\mathrm{s}$ |
| Conditions aux limites | `boundary_conditions` | `left`: `DIRICHLET` (`POISEUILLE_STEP`, `y_step: 1.0`, `u_mean: 1.0`) ; `right`: `NEUMANN` ; `bottom`/`top`: `DIRICHLET (0,0)` |
| Solveurs et précond. | `solver` | `AL` ($r_{\mathrm{al}}=10$), `AAL`, `MONOLITHIC`, `DIRECT` ; `BICGSTAB2` / `GMRES` ; `BLOCK_TRIANGULAR` ($\mathrm{tol}=10^{-8}$, `max_iter=800`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_9_marche_re400_256x64.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/re400_256x64"`.)*

```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 256,
    "ny": 64,
    "lx": 40.0,
    "ly": 2.0612244898
  },
  "physics": {
    "enable_two_phase": false,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "NEWTON",
    "fluid1": { "rho": 1.0, "mu": 0.005 },
    "gx": 0.0,
    "gy": 0.0
  },
  "obstacles": [
    {
      "type": "BOX",
      "x0": 0.0,
      "y0": 0.0,
      "x1": 4.0,
      "y1": 1.0,
      "permeability": 1.0e-14
    }
  ],
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.05,
    "t_final": 100.0,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "profile": "POISEUILLE_STEP", "y_step": 1.0, "u_mean": 1.0, "u": 1.0, "v": 0.0 },
    "right":  { "type": "NEUMANN",   "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "r_al": 10.0,
    "uzawa_iter": 3,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 800,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "output": {
    "output_dir": "output/cpp/re400_256x64",
    "case_name": "cas_9_marche_re400_256x64",
    "vtk_interval": 200,
    "log_interval": 20,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/re400_256x64 output/fortran/re400_256x64

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_9_marche_re400_256x64.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_9_marche_re400_256x64.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification Dynamique
1. **Longueurs de recollement :** Déterminer par interpolation linéaire sur la première couche de mailles près de $y=0$ (pour $x > 4.0$) le zéro de $u(x, \Delta y/2)$ donnant $X_r/H = (x_{\mathrm{recol}} - 4.0)/H$, ainsi que les points de séparation/recollement sur $y = L_y$ pour $Re \ge 400$.
2. **Conservation globale du débit et rigidité de l'obstacle :** Contrôler $\int u(0, y)\,\mathrm{d}y \approx \int u(L_x, y)\,\mathrm{d}y$ et vérifier que la vitesse résiduelle dans $\Omega_{\mathrm{step}}$ est $\mathcal{O}(K) \ll 10^{-8}\,\mathrm{m/s}$.
3. **Parité C++17 / Fortran 2018 et solveurs :** Comparer $X_r/H$, les itérations Krylov et la divergence maximale entre les deux moteurs.

### 6.2 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Extraction dynamique de Xr/H pour le Tutoriel 9 (Marche descendante)."""
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

def extract_reattachment(pvtr_path: str, lx: float = 40.0, x_step: float = 4.0) -> float:
    p = Path(pvtr_path)
    root = ET.parse(p).getroot().find("PRectilinearGrid")
    ext = [int(x) for x in root.attrib["WholeExtent"].split()]
    nx = ext[1] - ext[0]
    u_bot = np.zeros(nx)
    for piece in root.findall("Piece"):
        pe = [int(x) for x in piece.attrib["Extent"].split()]
        if pe[2] == 0:
            vtr = ET.parse(p.parent / piece.attrib["Source"]).getroot().find("RectilinearGrid").find("Piece")
            for da in vtr.find("CellData").findall("DataArray"):
                if da.attrib.get("Name") == "Velocity":
                    arr = np.fromstring(da.text, sep=" ").reshape((pe[3]-pe[2], pe[1]-pe[0], 3))
                    u_bot[pe[0]:pe[1]] = arr[0, :, 0]
    dx = lx / nx
    xc = (np.arange(nx) + 0.5) * dx
    mask = xc > x_step + 0.25
    sub_x, sub_u = xc[mask], u_bot[mask]
    idx = np.where((sub_u[:-1] <= 0.0) & (sub_u[1:] > 0.0))[0]
    if len(idx) == 0:
        return float("nan")
    i = idx[0]
    xr = sub_x[i] - sub_u[i] * (sub_x[i+1] - sub_x[i]) / (sub_u[i+1] - sub_u[i])
    return (xr - x_step)

if __name__ == "__main__":
    xr = extract_reattachment("output/cpp/re400_256x64/cas_9_marche_re400_256x64_002000.pvtr")
    print(f"Xr/H = {xr:.4f}")
```

### 6.3 Visualisation sous ParaView
1. Ouvrir `output/cpp/re400_256x64/cas_9_marche_re400_256x64_002000.pvtr` dans ParaView.
2. Filtrer la marche solide par **Threshold** ou **Clip** et tracer les lignes de courant (**Stream Tracer**) derrière l'arête $(4.0, 1.0)$ pour visualiser la bulle primaire inférieure et la zone de séparation supérieure.
# Tutoriel 10 : Sillage Laminaire Stationnaire autour d'un Cylindre Immergé ($Re = 20$)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_10_sillage_cylindre/` *(validation rapide `Run_Valide-CASES`, cas 10)*.

Ce tutoriel valide la représentation d'obstacles courbes immergés sur grille cartésienne par pénalisation volumique de Darcy-Brinkman (`"type": "CYLINDER"` dans `"obstacles"`) sur l'écoulement laminaire stationnaire autour d'un cylindre circulaire à $Re = 20$.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Distinguer un obstacle rigide fixe de Darcy-Brinkman (`obstacles` avec `type: "CYLINDER"`) d'une phase fluide mobile diphasique (`initial_condition.type: "CYLINDER"`).
- Calculer les efforts hydrodynamiques par intégration volumique de la force de Brinkman : $F_x = \int_{\Omega_{\mathrm{cyl}}} (\mu / \mathcal{K}) u\,\mathrm{d}\Omega$.
- Mesurer dynamiquement la longueur de recirculation symétrique $L_r/D$, l'angle de décollement $\theta_s$ et le coefficient de traînée $C_d$.
- Valider la symétrie transverse et la convergence spatiale entre les deux moteurs MPI `covof2d_cpp` et `covof2d_fortran`.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Grandeurs de Référence

Dans le domaine $\Omega = [0, 30] \times [0, 15]\,\mathrm{m}^2$, un cylindre circulaire fixe de diamètre $D = 1.0\,\mathrm{m}$ ($R = 0.5\,\mathrm{m}$) est placé en $(x_c, y_c) = (7.5, 7.5)\,\mathrm{m}$ (`permeability: 1e-14`). Une vitesse uniforme $U_\infty = 1.0\,\mathrm{m/s}$ est prescrite à l'entrée $x=0$, avec sortie `NEUMANN` à $x=30$ et glissement libre `FREE_SLIP` sur $y=0$ et $y=15$. Avec $\rho = 1.0\,\mathrm{kg/m}^3$ et $\mu = 0.05\,\mathrm{Pa\cdot s}$ :
$$Re = \frac{\rho U_\infty D}{\mu} = 20, \qquad C_d = \frac{2 F_x}{\rho U_\infty^2 D}, \qquad C_l = \frac{2 F_y}{\rho U_\infty^2 D}.$$

Références externes de la littérature (Dennis & Chang 1970 ; Fornberg 1980 ; Linnick & Fasel 2005) : $L_r/D \approx 0.90\text{--}0.93$, $C_d \approx 2.04\text{--}2.09$ et $C_l = 0$ par symétrie géométrique.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et maillages | `mesh.{nx,ny,lx,ly}` | $L_x = 30.0\,\mathrm{m}$, $L_y = 15.0\,\mathrm{m}$ ; $200\times 100$, $400\times 200$, $800\times 400$ |
| Cylindre immergé | `obstacles[0]` | `{"type": "CYLINDER", "xc": 7.5, "yc": 7.5, "radius": 0.5, "permeability": 1e-14}` |
| Fluide et convection | `physics` | `enable_two_phase: false`, $\rho=1.0$, $\mu=0.05$, `CENTERED` vs `TVD_SUPERBEE`, `NEWTON` vs `AB2` |
| Temps | `time.{time_scheme,dt,t_final}` | `"GEAR"`, $\Delta t = 0.05\,\mathrm{s}$, $t_{\mathrm{final}} = 50.0\,\mathrm{s}$ |
| Conditions aux limites | `boundary_conditions` | `left`: `DIRICHLET (1,0)` ; `right`: `NEUMANN` ; `bottom`/`top`: `FREE_SLIP` |
| Solveurs et précond. | `solver` | `AL` ($r_{\mathrm{al}}=10$), `AAL`, `MONOLITHIC`, `SP` ; `BICGSTAB2` ; `BLOCK_TRIANGULAR` ($\mathrm{tol}=10^{-8}$, `max_iter=800`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_10_sillage_re20_400x200.json`
*(Pour Fortran 2018, utiliser une copie dans `case/fortran/` avec `"output_dir": "output/fortran/grid_400x200"`.)*

```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 400,
    "ny": 200,
    "lx": 30.0,
    "ly": 15.0
  },
  "physics": {
    "enable_two_phase": false,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "NEWTON",
    "fluid1": { "rho": 1.0, "mu": 0.05 },
    "gx": 0.0,
    "gy": 0.0
  },
  "obstacles": [
    {
      "type": "CYLINDER",
      "xc": 7.5,
      "yc": 7.5,
      "radius": 0.5,
      "permeability": 1.0e-14
    }
  ],
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.05,
    "t_final": 50.0,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u": 1.0, "v": 0.0 },
    "right":  { "type": "NEUMANN",   "u": 0.0, "v": 0.0 },
    "bottom": { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "top":    { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "r_al": 10.0,
    "uzawa_iter": 3,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 800,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "output": {
    "output_dir": "output/cpp/grid_400x200",
    "case_name": "cas_10_sillage_re20_400x200",
    "vtk_interval": 200,
    "log_interval": 20,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/grid_400x200 output/fortran/grid_400x200

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_10_sillage_re20_400x200.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_10_sillage_re20_400x200.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification Dynamique
1. **Longueur du sillage $L_r/D$ :** Identifier sur la ligne médiane $y = y_c = 7.5\,\mathrm{m}$ en aval du bord arrière ($x > 8.0\,\mathrm{m}$) l'abscisse $x_{\mathrm{stag}}$ où $u(x, y_c)$ s'annule et calculer $L_r/D = (x_{\mathrm{stag}} - 8.0) / 1.0$.
2. **Traînée de Brinkman $C_d$ et symétrie $|C_l|$ :** Intégrer $(\mu/\mathcal{K}) u\,\Delta x\Delta y$ et $(\mu/\mathcal{K}) v\,\Delta x\Delta y$ sur les cellules contenues dans $(x-x_c)^2+(y-y_c)^2 \le R^2$ et contrôler la symétrie transverse $|u(x, y_c+\eta) - u(x, y_c-\eta)|$.
3. **Convergence et parité C++17 / Fortran 2018 :** Étudier la convergence sur $200\times 100 \to 400\times 200 \to 800\times 400$ et vérifier la cohérence inter-codes et inter-solveurs.

### 6.2 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Extraction dynamique de Lr/D et Cd de Brinkman pour le Tutoriel 10 (Sillage cylindre Re=20)."""
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

def analyze_wake(pvtr_path: str, lx: float = 30.0, ly: float = 15.0,
                 mu: float = 0.05, perm: float = 1e-14):
    p = Path(pvtr_path)
    root = ET.parse(p).getroot().find("PRectilinearGrid")
    ext = [int(x) for x in root.attrib["WholeExtent"].split()]
    nx, ny = ext[1] - ext[0], ext[3] - ext[2]
    u, v = np.zeros((ny, nx)), np.zeros((ny, nx))
    for piece in root.findall("Piece"):
        pe = [int(x) for x in piece.attrib["Extent"].split()]
        vtr = ET.parse(p.parent / piece.attrib["Source"]).getroot().find("RectilinearGrid").find("Piece")
        for da in vtr.find("CellData").findall("DataArray"):
            if da.attrib.get("Name") == "Velocity":
                arr = np.fromstring(da.text, sep=" ").reshape((pe[3]-pe[2], pe[1]-pe[0], 3))
                u[pe[2]:pe[3], pe[0]:pe[1]] = arr[:, :, 0]
                v[pe[2]:pe[3], pe[0]:pe[1]] = arr[:, :, 1]
    dx, dy = lx / nx, ly / ny
    xc = (np.arange(nx) + 0.5) * dx
    yc = (np.arange(ny) + 0.5) * dy
    X, Y = np.meshgrid(xc, yc)
    inside = (X - 7.5)**2 + (Y - 7.5)**2 <= 0.5**2
    fx = float(np.sum((mu / perm) * u[inside]) * dx * dy)
    cd = 2.0 * fx / (1.0 * 1.0 * 1.0)
    u_mid = 0.5 * (u[ny//2 - 1, :] + u[ny//2, :])
    mask = xc > 8.0
    sub_x, sub_u = xc[mask], u_mid[mask]
    idx = np.where((sub_u[:-1] <= 0.0) & (sub_u[1:] > 0.0))[0]
    lr = float("nan") if len(idx) == 0 else (sub_x[idx[0]] - sub_u[idx[0]]*(sub_x[idx[0]+1]-sub_x[idx[0]])/(sub_u[idx[0]+1]-sub_u[idx[0]]) - 8.0)
    print(f"[{pvtr_path}] Cd = {cd:.4f} | Lr/D = {lr:.4f}")

if __name__ == "__main__":
    analyze_wake("output/cpp/grid_400x200/cas_10_sillage_re20_400x200_001000.pvtr")
```

### 6.3 Visualisation sous ParaView
1. Ouvrir `output/cpp/grid_400x200/cas_10_sillage_re20_400x200_001000.pvtr` dans ParaView.
2. Superposer le champ `Pressure` et les lignes de courant (**Stream Tracer**) autour de $(7.5, 7.5)$ pour observer les deux lobes tourbillonnaires symétriques attachés en aval du cylindre.
# Tutoriel 11 : Convection Naturelle Laminaire en Cavité Différentiellement Chauffée ($Pr = 0.71, Ra = 10^3$ à $10^6$)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_11_cavite_chauffage_differentiel/` *(validation rapide `Run_Valide-CASES`, cas 11)*.

Ce tutoriel valide le solveur thermique implicite (`physics.enable_energy = true`), les conditions aux limites thermiques mixtes (Dirichlet isothermes et Neumann adiabatiques) et le terme de flottabilité de Boussinesq sur le benchmark de De Vahl Davis (1983) / Le Quéré (1991).

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Renseigner avec précision les propriétés thermiques requises (`Cp`, `k_thermal`, `beta` dans `physics.fluid1`, et `T_ref` avec `T` majuscule dans `physics`).
- Paramétrer indépendamment les conditions aux limites hydrodynamiques (`DIRICHLET` adhérent) et thermiques (`T_type: "DIRICHLET"` ou `"NEUMANN"`).
- Extraire dynamiquement le nombre de Nusselt moyen $\overline{Nu}(x=0)$ et $\overline{Nu}(x=H)$ depuis `Temperature` dans `.vtr`/`.pvtr`.
- Vérifier l'ordre de convergence spatiale 2 et la centro-symétrie exacte du champ thermo-convectif pour $Ra \in \{10^3, 10^4, 10^5, 10^6\}$.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `NavierStokesSolver` (`solve_energy`, flottabilité), `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm` (`apply_thermal_bc`), `module_navier_stokes` (`ns_solve_energy`), `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Références Bibliographiques

Sur la cavité carrée $\Omega = [0, H] \times [0, H]$ ($H = 1.0\,\mathrm{m}$) à parois rigides adhérentes ($\mathbf{u} = \mathbf{0}$), le système couple Navier-Stokes incompressible avec poussée de Boussinesq $-\rho_0 \beta (T - T_{\mathrm{ref}})\mathbf{g}$ et l'équation de l'énergie :
$$\rho_0 C_p \left(\frac{\partial T}{\partial t} + \mathbf{u} \cdot \nabla T\right) = k_{\mathrm{th}} \Delta T.$$

Avec $\rho_0 = 1.0$, $\mu = 1.0$, $C_p = 0.71$, $k_{\mathrm{th}} = 1.0$, $\beta = 1.0$, $T_{\mathrm{ref}} = 0$, $T_G(0,y) = +0.5$, $T_D(H,y) = -0.5$ ($\Delta T = 1.0$) et $\partial_y T = 0$ en $y=0, H$, le nombre de Rayleigh est fixé par $g_y = -Ra / Pr$ ($Pr = \mu C_p / k_{\mathrm{th}} = 0.71$). Les valeurs de référence (De Vahl Davis 1983 ; Le Quéré 1991) sont :
$$\overline{Nu}_{\mathrm{r\acute{e}f}}(10^3) = 1.11780, \qquad \overline{Nu}_{\mathrm{r\acute{e}f}}(10^4) = 2.24475, \qquad \overline{Nu}_{\mathrm{r\acute{e}f}}(10^5) = 4.52164, \qquad \overline{Nu}_{\mathrm{r\acute{e}f}}(10^6) = 8.82513.$$

---

## 3. Paramètres Physiques, Numériques et Discrétisation

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et maillages | `mesh.{nx,ny,lx,ly}` | $H = 1.0\,\mathrm{m}$, $N \in \{16, 32, 64, 128, 256, 512\}$ |
| Thermophysique | `physics` | `enable_energy: true`, `T_ref: 0.0`, $\rho_0=1.0, \mu=1.0, C_p=0.71, k_{\mathrm{th}}=1.0, \beta=1.0$, $g_y = -Ra/0.71$ |
| Convection et temps | `physics` & `time` | `CENTERED` vs `UPWIND`, `PICARD`/`AB2`, `"GEAR"`, $\Delta t = 0.002\,\mathrm{s}$ (pour $Ra=10^4$), $t_{\mathrm{final}} = 2.0\,\mathrm{s}$ |
| Conditions aux limites | `boundary_conditions` | `u=v=0` partout ; `left`: `T_type: "DIRICHLET", T: 0.5` ; `right`: `T_type: "DIRICHLET", T: -0.5` ; `bottom`/`top`: `T_type: "NEUMANN"` |
| Solveurs et précond. | `solver` | `AL` ($r_{\mathrm{al}}=100$), `AAL`, `MONOLITHIC` ; `BICGSTAB2` / `GMRES` ; `BLOCK_TRIANGULAR` ($\mathrm{tol}=10^{-9}$, `max_iter=800`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_11_thermique_ra1e4_n128.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/ra1e4_n128"`.)*

```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 1.0,
    "ly": 1.0
  },
  "physics": {
    "enable_two_phase": false,
    "enable_energy": true,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "PICARD",
    "T_ref": 0.0,
    "gx": 0.0,
    "gy": -14084.507042253521,
    "fluid1": {
      "rho": 1.0,
      "mu": 1.0,
      "Cp": 0.71,
      "k_thermal": 1.0,
      "beta": 1.0
    }
  },
  "initial_condition": {
    "t_init": 0.0,
    "u0": 0.0,
    "v0": 0.0
  },
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.002,
    "t_final": 2.0,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u": 0.0, "v": 0.0, "T_type": "DIRICHLET", "T":  0.5 },
    "right":  { "type": "DIRICHLET", "u": 0.0, "v": 0.0, "T_type": "DIRICHLET", "T": -0.5 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0, "T_type": "NEUMANN",   "T":  0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0, "T_type": "NEUMANN",   "T":  0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "r_al": 100.0,
    "uzawa_iter": 3,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-09,
    "max_iter": 800,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "output": {
    "output_dir": "output/cpp/ra1e4_n128",
    "case_name": "cas_11_thermique_ra1e4_n128",
    "vtk_interval": 100,
    "log_interval": 20,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/ra1e4_n128 output/fortran/ra1e4_n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_11_thermique_ra1e4_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_11_thermique_ra1e4_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification Dynamique
1. **Nombre de Nusselt moyen :** Calculer $\overline{Nu}(x=0) = \frac{1}{H}\int_0^H -\partial_x T(0,y)\,\mathrm{d}y$ à l'ordre 2 sur la paroi chaude (`T_0 = 0.5`) et comparer avec $\overline{Nu}(x=H)$ et la référence de De Vahl Davis / Le Quéré.
2. **Centro-symétrie et bilan énergétique :** Vérifier $\max |T(H-x, H-y) + T(x, y)| \ll 10^{-4}$ et $|\overline{Nu}_{x=0} - \overline{Nu}_{x=H}| \ll 10^{-3}$.
3. **Convergence de Richardson et parité C++17 / Fortran 2018 :** Vérifier l'ordre 2 sur $\overline{Nu}(N)$ pour $N \in \{16, 32, 64, 128, 256, 512\}$.

### 6.2 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Calcul dynamique du Nusselt moyen et de la centro-symétrie (Tutoriel 11)."""
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

def check_thermal_cavity(pvtr_path: str):
    p = Path(pvtr_path)
    root = ET.parse(p).getroot().find("PRectilinearGrid")
    ext = [int(x) for x in root.attrib["WholeExtent"].split()]
    nx, ny = ext[1] - ext[0], ext[3] - ext[2]
    T = np.zeros((ny, nx))
    for piece in root.findall("Piece"):
        pe = [int(x) for x in piece.attrib["Extent"].split()]
        vtr = ET.parse(p.parent / piece.attrib["Source"]).getroot().find("RectilinearGrid").find("Piece")
        for da in vtr.find("CellData").findall("DataArray"):
            if da.attrib.get("Name") == "Temperature":
                T[pe[2]:pe[3], pe[0]:pe[1]] = np.fromstring(da.text, sep=" ").reshape((pe[3]-pe[2], pe[1]-pe[0]))
    dx = 1.0 / nx
    # Dérivée quadratique à la paroi gauche T_wall = 0.5 (centres en dx/2 et 3dx/2)
    nu_hot = float(np.mean((8.0 * 0.5 - 9.0 * T[:, 0] + T[:, 1]) / (3.0 * dx)))
    sym_err = float(np.max(np.abs(T + T[::-1, ::-1])))
    print(f"[{pvtr_path}] Nu_hot = {nu_hot:.6f} | sym_err(T) = {sym_err:.3e}")

if __name__ == "__main__":
    check_thermal_cavity("output/cpp/ra1e4_n128/cas_11_thermique_ra1e4_n128_001000.pvtr")
```

### 6.3 Visualisation sous ParaView
1. Ouvrir `output/cpp/ra1e4_n128/cas_11_thermique_ra1e4_n128_001000.pvtr` dans ParaView.
2. Appliquer le filtre **Contour** sur `Temperature` (10 isothermes entre `-0.5` et `+0.5`) et **Stream Tracer** sur `Velocity` pour observer la stratification thermique centrale et les couches limites thermiques pariétales.
# Tutoriel 12 : Convection Naturelle en Cavité Poreuse de Darcy-Brinkman ($Ra_f = 10^7, Da = 10^{-6}, Pr = 1$)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_12_convection_poreux_darcy_brinkman/` *(validation rapide `Run_Valide-CASES`, cas 12)*.

Ce tutoriel valide le couplage entre la convection naturelle de Boussinesq et la perméabilité volumique finie de Darcy-Brinkman (`obstacles` avec $\mathcal{K} > 10^{-12}\,\mathrm{m}^2$ sur l'ensemble du domaine), permettant d'explorer la transition du régime de Darcy pur ($Da \ll 1$) vers le régime de Navier-Stokes libre.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Configurer un milieu poreux isotrope saturé par une boîte poreuse recouvrant $\Omega$ (`"permeability": 1.0e-8`, $Da = \mathcal{K}/H^2 = 10^{-6}$).
- Vérifier l'activation simultanée de la résistance de Darcy $-(\mu/\mathcal{K})\mathbf{u}$, de la diffusion visqueuse de Brinkman $\mu \Delta \mathbf{u}$ et de l'équation d'énergie thermique.
- Mesurer l'ordre de convergence spatiale 2 sur $\overline{Nu}(N)$ ainsi que l'effet du nombre de Darcy $Da \in \{10^{-6}, 10^{-4}, 10^{-2}\}$.
- Valider la parité C++17 / Fortran 2018 en parallèle MPI pur (`covof2d_cpp` et `covof2d_fortran`).

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `NavierStokesSolver` (`get_*_obstacle_perm`, `solve_energy`), `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Grandeurs de Référence

Dans la cavité carrée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ ($H = 0.1\,\mathrm{m}$), la matrice poreuse uniforme possède une perméabilité $\mathcal{K} = 10^{-8}\,\mathrm{m}^2$. Avec $\rho = 1.0\,\mathrm{kg/m}^3$, $\mu = 10^{-5}\,\mathrm{Pa\cdot s}$, $C_p = 1000\,\mathrm{J/(kg\cdot K)}$, $k_{\mathrm{th}} = 0.01\,\mathrm{W/(m\cdot K)}$, $\beta = 10^{-3}\,\mathrm{K}^{-1}$, $g_y = -10.0\,\mathrm{m/s}^2$, $T_G = +0.5$ et $T_D = -0.5$ ($\Delta T = 1.0\,\mathrm{K}$, $T_{\mathrm{ref}} = 0$) :
$$Pr = \frac{\mu C_p}{k_{\mathrm{th}}} = 1.0, \qquad Da = \frac{\mathcal{K}}{H^2} = 10^{-6}, \qquad Ra_f = \frac{\rho^2 g \beta \Delta T H^3 C_p}{\mu k_{\mathrm{th}}} = 10^7, \qquad Ra^* = Ra_f \, Da = 10.$$

Référence de la littérature (Lauriat & Prasad 1989 ; Nithiarasu et al. 1997) : en régime de Darcy pur ($Da \to 0$) à $Ra^* = 10$, $\overline{Nu}_{\mathrm{th}} \approx 1.079$, valeur légèrement réduite par l'effet visqueux de paroi de Brinkman à $Da = 10^{-6}$.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et maillages | `mesh.{nx,ny,lx,ly}` | $H = 0.1\,\mathrm{m}$, grilles $32^2, 64^2, 128^2, 256^2, 512^2$ |
| Matrice poreuse | `obstacles[0]` | `{"type": "BOX", "x0": 0.0, "y0": 0.0, "x1": 0.1, "y1": 0.1, "permeability": 1e-8}` ($Da \in \{10^{-6}, 10^{-4}, 10^{-2}\}$) |
| Thermophysique | `physics` | `enable_energy: true`, `T_ref: 0.0`, $\rho=1.0, \mu=10^{-5}, C_p=1000, k_{\mathrm{th}}=0.01, \beta=10^{-3}, g_y=-10.0$ |
| Temps | `time.{time_scheme,dt,t_final}` | `"GEAR"`, $\Delta t = 0.05\,\mathrm{s}$, $t_{\mathrm{final}} = 10.0\,\mathrm{s}$ |
| Conditions aux limites | `boundary_conditions` | `u=v=0` (`DIRICHLET`) ; `left`: $T=+0.5$ ; `right`: $T=-0.5$ ; `bottom`/`top`: `T_type: "NEUMANN"` |
| Solveurs et précond. | `solver` | `AL` ($r_{\mathrm{al}}=100$), `MONOLITHIC`, `DIRECT` ; `BICGSTAB2` ; `BLOCK_TRIANGULAR` ($\mathrm{tol}=10^{-10}$, `max_iter=600`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_12_poreux_n128.json`
*(Pour Fortran 2018, utiliser une copie dans `case/fortran/` avec `"output_dir": "output/fortran/n128"`.)*

```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 0.1,
    "ly": 0.1
  },
  "physics": {
    "enable_two_phase": false,
    "enable_energy": true,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "PICARD",
    "T_ref": 0.0,
    "gx": 0.0,
    "gy": -10.0,
    "fluid1": {
      "rho": 1.0,
      "mu": 1.0e-05,
      "Cp": 1000.0,
      "k_thermal": 0.01,
      "beta": 0.001
    }
  },
  "obstacles": [
    {
      "type": "BOX",
      "x0": 0.0,
      "y0": 0.0,
      "x1": 0.1,
      "y1": 0.1,
      "permeability": 1.0e-08
    }
  ],
  "initial_condition": {
    "t_init": 0.0,
    "u0": 0.0,
    "v0": 0.0
  },
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.05,
    "t_final": 10.0,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u": 0.0, "v": 0.0, "T_type": "DIRICHLET", "T":  0.5 },
    "right":  { "type": "DIRICHLET", "u": 0.0, "v": 0.0, "T_type": "DIRICHLET", "T": -0.5 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0, "T_type": "NEUMANN",   "T":  0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0, "T_type": "NEUMANN",   "T":  0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "r_al": 100.0,
    "uzawa_iter": 3,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-10,
    "max_iter": 600,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "output": {
    "output_dir": "output/cpp/n128",
    "case_name": "cas_12_poreux_n128",
    "vtk_interval": 50,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/n128 output/fortran/n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_12_poreux_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_12_poreux_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification Dynamique
1. **Convergence spatiale d'ordre 2 sur $\overline{Nu}$ :** Extraire $\overline{Nu}(N)$ de manière identique pour C++ et Fortran depuis `Temperature` dans `.vtr`/`.pvtr`, et calculer les ordres de Richardson $p$ sur $32^2 \to 512^2$.
2. **Transition Darcy $\to$ Brinkman :** Varier `permeability` ($\mathcal{K} \in \{10^{-8}, 10^{-6}, 10^{-4}\}\,\mathrm{m}^2$, soit $Da \in \{10^{-6}, 10^{-4}, 10^{-2}\}$) et observer l'épaississement de la couche limite de Brinkman et l'augmentation de $\overline{Nu}$.
3. **Parité C++17 / Fortran 2018 :** Vérifier $|\overline{Nu}_{\mathrm{C++}} - \overline{Nu}_{\mathrm{F90}}| \ll 10^{-6}$.

### 6.2 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique de Nu moyen pour le Tutoriel 12 (Cavité poreuse Darcy-Brinkman)."""
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

def porous_nusselt(pvtr_path: str, H: float = 0.1) -> float:
    p = Path(pvtr_path)
    root = ET.parse(p).getroot().find("PRectilinearGrid")
    ext = [int(x) for x in root.attrib["WholeExtent"].split()]
    nx, ny = ext[1] - ext[0], ext[3] - ext[2]
    T = np.zeros((ny, nx))
    for piece in root.findall("Piece"):
        pe = [int(x) for x in piece.attrib["Extent"].split()]
        vtr = ET.parse(p.parent / piece.attrib["Source"]).getroot().find("RectilinearGrid").find("Piece")
        for da in vtr.find("CellData").findall("DataArray"):
            if da.attrib.get("Name") == "Temperature":
                T[pe[2]:pe[3], pe[0]:pe[1]] = np.fromstring(da.text, sep=" ").reshape((pe[3]-pe[2], pe[1]-pe[0]))
    dx = H / nx
    nu_hot = float(np.mean(H * (8.0 * 0.5 - 9.0 * T[:, 0] + T[:, 1]) / (3.0 * dx)))
    print(f"[{pvtr_path}] Nu_hot = {nu_hot:.6f}")
    return nu_hot

if __name__ == "__main__":
    porous_nusselt("output/cpp/n128/cas_12_poreux_n128_000200.pvtr")
```

### 6.3 Visualisation sous ParaView
1. Ouvrir `output/cpp/n128/cas_12_poreux_n128_000200.pvtr` dans ParaView.
2. Tracer `Temperature` (isothermes quasiment conductives à $Ra^*=10$ puis fortement convectives lorsque $Da$ augmente) et `Velocity` le long des parois verticales.
# Tutoriel 13 : Écoulement Laminaire autour d'un Obstacle Rectangulaire Fin en Canal ($Re = 82.5$)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_1_obstacle_fin/` *(campagne avancée `Run_CASES`, cas 1)*.

Ce tutoriel met en œuvre la pénalisation rigide de Darcy-Brinkman sur un obstacle mince vertical fixe en canal (`BOX`, $e = 2\,\mathrm{mm}$, $S/H = 0.514$) soumis à une entrée parabolique de Poiseuille (`profile: "PARABOLIC"`).

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Configurer un obstacle mince rigide de Brinkman (`obstacles` avec `permeability: 1e-14`) et analyser la sensibilité à $\mathcal{K} \in \{10^{-10}, 10^{-12}, 10^{-14}, 10^{-16}\}\,\mathrm{m}^2$.
- Paramétrer le profil d'entrée parabolique plein canal (`"profile": "PARABOLIC"`, `u_mean`, `u`) couplé à une sortie `NEUMANN`.
- Extraire dynamiquement les longueurs de recirculation aval $x_1/S$ (paroi inférieure) et $(x_2/S, x_3/S)$ (paroi supérieure), la perte de charge $\Delta P$ et les 9 profils axiaux $u(x^*, y/H)$.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Grandeurs de Référence

Dans le canal $\Omega = [0, 0.50] \times [0, 0.025]\,\mathrm{m}^2$ ($H = 0.025\,\mathrm{m}$), un obstacle mince de hauteur $S = 0.01285\,\mathrm{m}$ et d'épaisseur $e = 0.002\,\mathrm{m}$ est posé en $x_{\mathrm{obs}} = 0.15\,\mathrm{m}$ (`[0.149, 0.151] x [0.0, 0.01285]`). Le fluide ($\rho = 1.0\,\mathrm{kg/m}^3$, $\mu = 1.557576\times 10^{-4}\,\mathrm{Pa\cdot s}$) entre avec $u_{\mathrm{mean}} = 1.0\,\mathrm{m/s}$ ($u_{\max} = 1.5\,\mathrm{m/s}$), donnant $Re = \rho \bar{U} S / \mu = 82.5$. La littérature de référence rapporte un tourbillon primaire aval en paroi inférieure ($x_1/S \approx 4.99$) et un décollement supérieur entre $x_2/S \approx 4.01$ et $x_3/S \approx 9.19$.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et maillages | `mesh.{nx,ny,lx,ly}` | $L_x = 0.50\,\mathrm{m}$, $L_y = 0.025\,\mathrm{m}$ ; $250\times 50$, $500\times 100$, $1000\times 200$ |
| Obstacle mince | `obstacles[0]` | `{"type": "BOX", "x0": 0.149, "y0": 0.0, "x1": 0.151, "y1": 0.01285, "permeability": 1e-14}` |
| Fluide et convection | `physics` | `enable_two_phase: false`, $\rho=1.0$, $\mu = 1.557576\times 10^{-4}$, `CENTERED` vs `TVD_SUPERBEE` vs `UPWIND`, `NEWTON` vs `AB2` |
| Temps | `time.{time_scheme,dt,t_final}` | `"GEAR"`, $\Delta t = 5\times 10^{-4}\,\mathrm{s}$, $t_{\mathrm{final}} = 2.0\,\mathrm{s}$ |
| Conditions aux limites | `boundary_conditions` | `left`: `DIRICHLET` (`PARABOLIC`, `u_mean: 1.0`, `u: 1.5`) ; `right`: `NEUMANN` ; `bottom`/`top`: `DIRICHLET (0,0)` |
| Solveurs et précond. | `solver` | `AL` ($r_{\mathrm{al}}=10$), `AAL`, `MONOLITHIC` ; `BICGSTAB2` ; `BLOCK_TRIANGULAR` ($\mathrm{tol}=10^{-8}$, `max_iter=800`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_1_obstacle_fin_500x100.json`
*(Pour Fortran 2018, utiliser `"output_dir": "output/fortran/grid_500x100"`.)*

```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 500,
    "ny": 100,
    "lx": 0.50,
    "ly": 0.025
  },
  "physics": {
    "enable_two_phase": false,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "NEWTON",
    "fluid1": { "rho": 1.0, "mu": 0.0001557576 },
    "gx": 0.0,
    "gy": 0.0
  },
  "obstacles": [
    {
      "type": "BOX",
      "x0": 0.149,
      "y0": 0.0,
      "x1": 0.151,
      "y1": 0.01285,
      "permeability": 1.0e-14
    }
  ],
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.0005,
    "t_final": 2.0,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "profile": "PARABOLIC", "u_mean": 1.0, "u": 1.5, "v": 0.0 },
    "right":  { "type": "NEUMANN",   "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "r_al": 10.0,
    "uzawa_iter": 3,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 800,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "output": {
    "output_dir": "output/cpp/grid_500x100",
    "case_name": "cas_1_obstacle_fin_500x100",
    "vtk_interval": 200,
    "log_interval": 20,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

```bash
mkdir -p output/cpp/grid_500x100 output/fortran/grid_500x100
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_1_obstacle_fin_500x100.json
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_1_obstacle_fin_500x100.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

1. **Longueurs $x_1/S, x_2/S, x_3/S$ :** Détecter dynamiquement l'inversion de $u(x, \Delta y/2)$ et $u(x, H - \Delta y/2)$ pour $x > 0.151\,\mathrm{m}$ à partir du dernier `.pvtr`.
2. **Rigidité Brinkman et perte de charge :** Vérifier $\max |u|$ dans l'obstacle et $\Delta P = \bar{p}(0.05) - \bar{p}(0.45)$.

```python
#!/usr/bin/env python3
"""Extraction dynamique de x1/S sur le Tutoriel 13 (Obstacle fin en canal)."""
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

def check_thin_obstacle(pvtr_path: str, lx: float = 0.50, x_obs: float = 0.151, S: float = 0.01285):
    p = Path(pvtr_path)
    root = ET.parse(p).getroot().find("PRectilinearGrid")
    ext = [int(x) for x in root.attrib["WholeExtent"].split()]
    nx, ny = ext[1] - ext[0], ext[3] - ext[2]
    u = np.zeros((ny, nx))
    for piece in root.findall("Piece"):
        pe = [int(x) for x in piece.attrib["Extent"].split()]
        vtr = ET.parse(p.parent / piece.attrib["Source"]).getroot().find("RectilinearGrid").find("Piece")
        for da in vtr.find("CellData").findall("DataArray"):
            if da.attrib.get("Name") == "Velocity":
                u[pe[2]:pe[3], pe[0]:pe[1]] = np.fromstring(da.text, sep=" ").reshape((pe[3]-pe[2], pe[1]-pe[0], 3))[:, :, 0]
    xc = (np.arange(nx) + 0.5) * (lx / nx)
    m = xc > x_obs + 0.002
    sub_x, sub_u = xc[m], u[0, m]
    idx = np.where((sub_u[:-1] <= 0.0) & (sub_u[1:] > 0.0))[0]
    x1_s = float("nan") if len(idx) == 0 else (sub_x[idx[0]] - sub_u[idx[0]]*(sub_x[idx[0]+1]-sub_x[idx[0]])/(sub_u[idx[0]+1]-sub_u[idx[0]]) - 0.15) / S
    print(f"[{pvtr_path}] x1/S = {x1_s:.4f}")

if __name__ == "__main__":
    check_thin_obstacle("output/cpp/grid_500x100/cas_1_obstacle_fin_500x100_004000.pvtr")
```
# Tutoriel 14 : Campagne Dual-Code Newton, AAL et SVOF sur Cylindre Dense ($50 \times 100$)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_2_cylindre_dense_newton_aal/` *(campagne avancée `Run_CASES`, cas 2)*.

Ce tutoriel constitue une campagne paramétrique rapide ($50\times 100$, $\Delta x = \Delta y = 0.002\,\mathrm{m}$) dédiée à l'analyse croisée de la linéarisation convective de Newton (`NEWTON` vs `ADAMS_BASHFORTH` vs `PICARD`), du Lagrangien Augmenté Adaptatif (`AAL` vs `MONOLITHIC`), du schéma temporel (`GEAR` vs `EULER`) et du lissage SVOF ($\sigma = 0$) entre `covof2d_cpp` et `covof2d_fortran`.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comparer les trois linéarisations convectives (`NEWTON`, `ADAMS_BASHFORTH`, `PICARD`) et les trois schémas spatiaux (`CENTERED`, `UPWIND`, `TVD_SUPERBEE`) sur un saut $r_\rho \approx 85\,000$, $r_\mu \approx 5.4\times 10^9$.
- Quantifier le comportement de `AAL` ($K_{\mathrm{aal}} \in \{10, 100\}$, `uzawa_iter` $\in \{1, 3, 5\}$) et `MONOLITHIC` avec le préconditionneur `BLOCK_TRIANGULAR` 100 % parallèle MPI.
- Mesurer l'impact rhéologique pur du lissage SVOF de Helmholtz (`enable_svof_smoothing`) lorsque $\sigma = 0$ (`enable_diffusion_smoothing: false`).
- Vérifier la parité stricte C++17 $\leftrightarrow$ Fortran 2018 sur la trajectoire galiléenne $U_c(t) = -gt$, $Y_c(t) = Y_0 - \frac{1}{2}gt^2$.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Matrice des Paramètres

Sur $\Omega = [0, 0.1] \times [0, 0.2]\,\mathrm{m}^2$ discrétisé en $50\times 100$ mailles ($\Delta t = 6.25\times 10^{-5}\,\mathrm{s}$, $t_{\mathrm{final}} = 0.144\,\mathrm{s}$), le cylindre de rayon $R = 0.0125\,\mathrm{m}$ centré en $(0.05, 0.15)\,\mathrm{m}$ (`axis_a = axis_b = 0.0125`, `pert_amp = 0.0`) chute librement dans l'air ($\rho_1 = 1.1768$, $\mu_1 = 1.85\times 10^{-5}$, $\rho_2 = 10^5$, $\mu_2 = 10^5$, $\sigma = 0$). La solution analytique à $t_{\mathrm{final}} = 0.144\,\mathrm{s}$ vaut $U_c^{\mathrm{ex}} = -1.41264\,\mathrm{m/s}$ et $Y_c^{\mathrm{ex}} = 0.048288\,\mathrm{m}$.

| Catégorie | Clé JSON | Plage / Options étudiées |
| :--- | :--- | :--- |
| Convection | `physics.convective_scheme` / `convective_linearization` | `CENTERED`, `UPWIND`, `TVD_SUPERBEE` ; `NEWTON`, `ADAMS_BASHFORTH`, `PICARD` |
| Temps & SVOF | `time.time_scheme` ; `physics.enable_svof_smoothing` | `GEAR` vs `EULER` ; `true` ($k = 1.0$) vs `false` |
| Solveurs & Précond. | `solver` | `AAL` ($K_{\mathrm{aal}} \in \{10, 100\}$, `uzawa_iter` $\in \{1, 3, 5\}$), `MONOLITHIC` ; `BLOCK_TRIANGULAR`, `JACOBI`, `SSOR`, `ILU0` ($\mathrm{tol}=10^{-6}$, `max_iter=500`) |

---

## 3. Fichier de Configuration JSON Pas à Pas (`case/cpp/cas_2_aal_gear_newton_svof.json`)
*(Pour Fortran 2018, utiliser `"output_dir": "output/fortran/aal_gear_newton_svof"`.)*

```json
{
  "mode": "NAVIER_STOKES",
  "mesh": { "nx": 50, "ny": 100, "lx": 0.1, "ly": 0.2 },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": false,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "NEWTON",
    "enable_svof_smoothing": true,
    "svof_filter_width": 1.0,
    "enable_diffusion_smoothing": false,
    "sigma": 0.0,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 1.1768,   "mu": 1.85e-05 },
    "fluid2": { "rho": 100000.0, "mu": 100000.0 }
  },
  "initial_condition": {
    "type": "CYLINDER",
    "center_x": 0.05,
    "center_y": 0.15,
    "radius": 0.0125,
    "axis_a": 0.0125,
    "axis_b": 0.0125,
    "pert_amp": 0.0
  },
  "time": { "time_scheme": "GEAR", "dt": 6.25e-05, "t_final": 0.144, "cfl_max": 0.5, "adaptive_dt": false },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "ADAPTATIVE_AUGMENTED_LAGRANGIAN",
    "adaptative_al": true,
    "k_aal": 100.0,
    "uzawa_tol": 1e-06,
    "uzawa_iter": 3,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-06,
    "max_iter": 500,
    "r_al": 0.0,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": { "advection_method": "WEYMOUTH_YUE", "reconstruction_method": "PLIC", "tvd_limiter": "NONE", "cfl_vof": 0.5, "use_markers": false },
  "output": { "output_dir": "output/cpp/aal_gear_newton_svof", "case_name": "cas_2_aal_gear_newton_svof", "vtk_interval": 100, "log_interval": 100, "write_csv_history": true }
}
```

---

## 4. Procédure d'Exécution et Vérification Dynamique Dual-Code

```bash
mkdir -p output/cpp/aal_gear_newton_svof output/fortran/aal_gear_newton_svof
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_2_aal_gear_newton_aal.json
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_2_aal_gear_newton_aal.json
```

```python
#!/usr/bin/env python3
"""Vérification de parité C++/Fortran pour le Tutoriel 14 (Cylindre dense 50x100)."""
import csv, math
from pathlib import Path

def eval_csv(path: str, y0: float = 0.15, gy: float = -9.81):
    with Path(path).open(newline="") as f:
        rows = [r for r in csv.DictReader(f) if float(r["time"]) > 0.0]
    t = [float(r["time"]) for r in rows]
    uc = [float(r["rise_velocity"]) for r in rows]
    yc = [float(r["center_of_mass_y"]) for r in rows]
    err_u = math.sqrt(sum((a - gy*ti)**2 for a, ti in zip(uc, t)) / sum((gy*ti)**2 for ti in t))
    err_y = math.sqrt(sum((a - (y0 + 0.5*gy*ti*ti))**2 for a, ti in zip(yc, t)) / sum((y0 + 0.5*gy*ti*ti)**2 for ti in t))
    print(f"[{path}] Err_L2(U)={err_u:.3e} | Err_L2(Y)={err_y:.3e} | Uc_end={uc[-1]:.6f}")

if __name__ == "__main__":
    eval_csv("output/cpp/aal_gear_newton_svof/cas_2_aal_gear_newton_svof_history.csv")
    eval_csv("output/fortran/aal_gear_newton_svof/cas_2_aal_gear_newton_svof_history.csv")
```
# Tutoriel 15 : Chute Libre d'un Cylindre Dense dans le Vide (Forts Contrastes)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_3_dense_cylinder/` *(campagne avancée `Run_CASES`)*.

Ce tutoriel valide le couplage fluide-solide par approche 1-fluide à fort contraste de masse volumique et de viscosité ($\rho_2 / \rho_1 \approx 85\,000$, $\mu_2 / \mu_1 \approx 5.4\times 10^9$). Il démontre que l'inclusion cylindrique dense suit rigoureusement la cinématique galiléenne de chute libre dans le vide sans dispersion ni amortissement numérique parasite.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Maîtriser le comportement du couplage 1-fluide Navier-Stokes/VOF face à des sauts extrêmes de densité et de viscosité.
- Configurer le schéma temporel `GEAR` (BDF2), la linéarisation `ADAMS_BASHFORTH` et le lissage SVOF de Helmholtz (`enable_svof_smoothing`).
- Comparer les formulations de couplage `MONOLITHIC` (FC), `AUGMENTED_LAGRANGIAN` (AL), `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (AAL) et `STANDARD_PROJECTION` (SP).
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Évaluer dynamiquement les erreurs $L_2$ temporelles sur $U_c(t)$ et $Y_c(t)$ depuis `*_history.csv`.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas diphasique sévère simule la chute libre sous gravité d'un cylindre circulaire ultra-dense et ultra-visqueux dans l'air au repos. Il met à l'épreuve les solveurs de couplage vitesse-pression et le préconditionneur `BLOCK_TRIANGULAR` à densité variable face à des sauts $r_\rho = \rho_2/\rho_1 \approx 84\,976$ et $r_\mu = \mu_2/\mu_1 \approx 5.405\times 10^9$.

Dans une cavité rectangulaire $\Omega = [0, 0.1] \times [0, 0.2]\,\mathrm{m}^2$ soumise à la pesanteur $\mathbf{g} = (0, -9.81)\,\mathrm{m/s}^2$, un cylindre de rayon $R = 0.0125\,\mathrm{m}$ est initialement centré en $(x_c, Y_0) = (0.05, 0.15)\,\mathrm{m}$ sans vitesse initiale (`"type": "CYLINDER"`, qui affecte $C = 0$ à l'intérieur du cylindre et $C = 1$ à l'extérieur) :
- **Fluide 1 (Air ambiant, $C=1$)** : $\rho_1 = 1.1768\,\mathrm{kg/m}^3$, $\mu_1 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$.
- **Fluide 2 (Cylindre dense, $C=0$)** : $\rho_2 = 100\,000\,\mathrm{kg/m}^3$, $\mu_2 = 100\,000\,\mathrm{Pa\cdot s}$, $\sigma = 0$.

### 2.2 Solution Analytique Galiléenne et Formulations Numériques
Compte tenu du rapport d'inertie $r_\rho \gg 1$, la traînée aérodynamique et la poussée d'Archimède exercées par l'air sont négligeables devant le poids du cylindre. La vitesse verticale $U_c(t)$ et l'ordonnée du centre de masse $Y_c(t)$ obéissent à la chute libre galiléenne :
$$U_c^{\mathrm{ex}}(t) = g_y t = -9.81\,t, \qquad Y_c^{\mathrm{ex}}(t) = Y_0 + \frac{1}{2} g_y t^2 = 0.15 - 4.905\,t^2.$$
À l'instant final $t_{\mathrm{final}} = 0.144\,\mathrm{s}$, la cinématique théorique vaut :
$$U_c^{\mathrm{ex}}(0.144\,\mathrm{s}) = -1.41264\,\mathrm{m/s}, \qquad Y_c^{\mathrm{ex}}(0.144\,\mathrm{s}) = 0.048288\,\mathrm{m}.$$

Les briques numériques clés activées pour ce cas sont :
1. **Schéma temporel Gear / BDF2** (`"time_scheme": "GEAR"`) et linéarisation convective d'Adams-Bashforth (`"convective_linearization": "ADAMS_BASHFORTH"`).
2. **Lissage SVOF de Helmholtz** (`(I - k\Delta x\Delta y \Delta)C_s = C`) malgré $\sigma = 0$ : régularise le saut extrême $r_\rho \approx 85\,000, r_\mu \approx 5.4\times 10^9$ avec moyenne harmonique nodale $\mu_{\mathrm{node}}$ (à comparer impérativement avec `enable_svof_smoothing: false`).
3. **Schéma convectif `CENTERED` et budget Krylov :** `convective_scheme` est câblé dans l'opérateur complet (`CENTERED` recommandé selon *Fluids* 2021, à comparer à `UPWIND` et `TVD_SUPERBEE`), et `BLOCK_TRIANGULAR` (100 % parallèle MPI via CSR `Fuu`/`Fvv`/`SchurLaplace`) requiert `max_iter >= 500` sur $100\times 200$ (`800` sur $200\times 400$).
4. **Couplages comparés :** `MONOLITHIC` (FC), `AUGMENTED_LAGRANGIAN` (AL), `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (AAL avec facteur `k_aal` et critère `uzawa_tol`) et `STANDARD_PROJECTION` (SP).

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et Cylindre | $\Omega$, $(x_c, Y_0)$, $R$ | $[0, 0.1]\times[0, 0.2]\,\mathrm{m}^2$, $(0.05, 0.15)\,\mathrm{m}$, $R = 0.0125\,\mathrm{m}$ (`axis_a = axis_b = 0.0125`, `pert_amp = 0.0`) |
| Maillages et Pas de temps | $(N_x \times N_y)$, $\Delta t$, $t_{\mathrm{final}}$ | $100\times 200$ et $200\times 400$, $\Delta t = 6.25\times 10^{-5}\,\mathrm{s}$, $t_{\mathrm{final}} = 0.144\,\mathrm{s}$ (2304 pas) |
| Conditions aux limites | Gauche/Droite ; Bas/Haut | `FREE_SLIP` sur $x=0, 0.1$ ; `DIRICHLET` ($u=v=0$) sur $y=0, 0.2$ |
| Solveurs et Préconditionneurs | Couplage, Krylov, Précond. | `MONOLITHIC` (FC), `AUGMENTED_LAGRANGIAN` (AL), `AAL` (`k_aal` $\in \{10, 100, 1000\}$, réf. `10` / `100`, `max_iter = 500`/`800`), `STANDARD_PROJECTION` ; `BICGSTAB2`, `GMRES`, `BICGSTAB` ; `BLOCK_TRIANGULAR`, `JACOBI`, `SSOR`, `ILU0` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_8_fc_100x200.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/fc_100x200"`.)*

```json
{
  "case_name": "cas_8_fc_100x200",
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 100,
    "ny": 200,
    "lx": 0.1,
    "ly": 0.2
  },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": false,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_svof_smoothing": true,
    "svof_filter_width": 1.0,
    "enable_diffusion_smoothing": false,
    "sigma": 0.0,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 1.1768,   "mu": 1.85e-05 },
    "fluid2": { "rho": 100000.0, "mu": 100000.0 }
  },
  "initial_condition": {
    "type": "CYLINDER",
    "center_x": 0.05,
    "center_y": 0.15,
    "radius": 0.0125,
    "axis_a": 0.0125,
    "axis_b": 0.0125,
    "pert_amp": 0.0
  },
  "time": {
    "time_scheme": "GEAR",
    "dt": 6.25e-05,
    "t_final": 0.144,
    "cfl_max": 0.5,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "MONOLITHIC",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-06,
    "max_iter": 500,
    "r_al": 0.0,
    "uzawa_iter": 1,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.5,
    "use_markers": false
  },
  "output": {
    "output_dir": "output/cpp/fc_100x200",
    "case_name": "cas_8_fc_100x200",
    "vtk_interval": 100,
    "log_interval": 100,
    "write_csv_history": true
  }
}
```

Pour activer le Lagrangien Augmenté Adaptatif (`AAL`), remplacer le bloc `"solver"` par :

```json
  "solver": {
    "type": "ADAPTATIVE_AUGMENTED_LAGRANGIAN",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-06,
    "max_iter": 500,
    "r_al": 0.0,
    "adaptative_al": true,
    "k_aal": 10.0,
    "uzawa_iter": 5,
    "uzawa_tol": 1e-06,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  }
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

L'exécution utilise **MPI pur** (`mpirun -np 4` pour `MONOLITHIC`, `AL` et `AAL` avec solveurs Krylov ; `mpirun -np 1` si `STANDARD_PROJECTION` est sélectionné).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/fc_100x200 output/fortran/fc_100x200

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_8_fc_100x200.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_8_fc_100x200.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Erreurs relatives $L_2$ temporelles sur $U_c(t)$ et $Y_c(t)$ :**
   À partir des colonnes `time`, `center_of_mass_y` ($Y_c$) et `rise_velocity` ($U_c$) de `<case_name>_history.csv` (car pour `IC_CYLINDER`, $1 - C$ représente le cylindre dense), calculer dynamiquement :
   $$\mathrm{Err}_{\mathrm{rel},2}(U_c) = \sqrt{\frac{\sum_{n=1}^{N_t} \big(U_c(t^n) - U_c^{\mathrm{ex}}(t^n)\big)^2}{\sum_{n=1}^{N_t} \big(U_c^{\mathrm{ex}}(t^n)\big)^2}}, \qquad \mathrm{Err}_{\mathrm{rel},2}(Y_c) = \sqrt{\frac{\sum_{n=1}^{N_t} \big(Y_c(t^n) - Y_c^{\mathrm{ex}}(t^n)\big)^2}{\sum_{n=1}^{N_t} \big(Y_c^{\mathrm{ex}}(t^n)\big)^2}}.$$
   Critères cibles : $\mathrm{Err}_{\mathrm{rel},2}(U_c) < 10^{-2}$, $\mathrm{Err}_{\mathrm{rel},2}(Y_c) < 10^{-2}$, et variation relative de volume $|V(t_{\mathrm{final}}) - V(0)| / V(0) < 10^{-4}$ mesurée depuis `PhaseIndicator` dans les `.vtr`.
2. **Comparaison des couplages (`FC`, `AL`, `AAL`, `SP`) et du filtre SVOF :**
   Comparer la préservation de la circularité du cylindre, l'absence de traînée parasite et le nombre moyen d'itérations Krylov relevé dans les journaux d'exécution.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier la superposition des trajectoires $Y_c(t)$ et des vitesses $U_c(t)$ entre `covof2d_cpp` et `covof2d_fortran` ainsi que le maintien de $\|\nabla\cdot\mathbf{u}\|_\infty$ (`max_div`).

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` : `step,time,dt,center_of_mass_y,rise_velocity,bubble_volume,kinetic_energy,max_u,max_div` (`9` colonnes C++ / `11` colonnes Fortran avec `avg_rho_gas,avg_p_gas`).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `PhaseIndicator_Smoothed`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 15 (Chute libre cylindre dense) depuis *_history.csv."""
import csv
import math
from pathlib import Path

def verify_free_fall(csv_path: str, y0: float = 0.15, gy: float = -9.81) -> None:
    times, yc, uc = [], [], []
    with Path(csv_path).open(newline="") as f:
        for r in csv.DictReader(f):
            t = float(r["time"])
            if t > 0.0:
                times.append(t)
                yc.append(float(r["center_of_mass_y"]))
                uc.append(float(r["rise_velocity"]))

    u_ex = [gy * t for t in times]
    y_ex = [y0 + 0.5 * gy * t * t for t in times]
    err_u = math.sqrt(sum((a - b) ** 2 for a, b in zip(uc, u_ex)) / sum(b ** 2 for b in u_ex))
    err_y = math.sqrt(sum((a - b) ** 2 for a, b in zip(yc, y_ex)) / sum(b ** 2 for b in y_ex))
    print(f"[{csv_path}] t_end={times[-1]:.4f} s | Uc={uc[-1]:.6f} (ex={u_ex[-1]:.6f}) "
          f"| Yc={yc[-1]:.6f} (ex={y_ex[-1]:.6f}) | Err_L2(U)={err_u:.3e} | Err_L2(Y)={err_y:.3e}")

if __name__ == "__main__":
    verify_free_fall("output/cpp/fc_100x200/cas_8_fc_100x200_history.csv")
    verify_free_fall("output/fortran/fc_100x200/cas_8_fc_100x200_history.csv")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série temporelle `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Afficher le champ `Density` en échelle logarithmique pour distinguer le cylindre dense ($\rho_2 = 10^5\,\mathrm{kg/m}^3$) de l'air ambiant ($\rho_1 = 1.1768\,\mathrm{kg/m}^3$), et tracer l'isocontour `PhaseIndicator = 0.5` afin de vérifier l'absence de déformation du solide rigide durant sa chute.
3. Utiliser le filtre **Glyph** sur `Velocity` pour visualiser l'écoulement d'air déplacé autour du cylindre.
# Tutoriel 16 : Sédimentation Confinée 2D d'une Inclusion Cylindrique et Loi de Faxén

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_4_sedimentation_faxen/` *(campagne avancée `Run_CASES`)*.

Ce tutoriel modélise la sédimentation gravitaire d'un cylindre solide en régime rampant de Stokes ($\mathrm{Re}_p \ll 1$) entre deux parois planes verticales. La vitesse terminale quasi-stationnaire obtenue numériquement est confrontée à la loi analytique de Faxén intégrant l'effet de confinement pariétal.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre la physique du régime de Stokes confiné et la formule analytique de correction pariétale de Faxén.
- Configurer la pénalisation de viscosité ($\mu_p \gg \mu_f$) avec le solveur `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (`AAL`) et le lissage SVOF.
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Extraire dynamiquement la vitesse verticale du centre de masse de la particule depuis les sorties `.vtr` et la comparer à la valeur analytique de Faxén.
- Visualiser le champ de vitesse autour du cylindre sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide le couplage entre la pénalisation visqueuse d'une particule rigide ($\mu_p / \mu_f \approx 6530$) et la méthode du Lagrangien Augmenté Adaptatif (`ADAPTATIVE_AUGMENTED_LAGRANGIAN`) pour la sédimentation lente d'un cylindre confiné entre deux parois verticales ($\mathrm{Re}_p \ll 1$).

Dans une cavité carrée fermée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ ($2l = 0.1\,\mathrm{m}$, demi-largeur $l = 0.05\,\mathrm{m}$) à parois adhérentes (`DIRICHLET`, $\mathbf{u} = \mathbf{0}$), un cylindre de rayon $a = 0.00625\,\mathrm{m}$ (rapport de confinement $k_x = a/l = 1/8 = 0.125$) est centré en $(0.05, 0.05)\,\mathrm{m}$ avec `"type": "DROPLET"` ($C = 1$ dans la particule). Les propriétés des deux phases sont :
- **Fluide 1 (Particule solide, $C=1$)** : $\rho_p = 1500\,\mathrm{kg/m}^3$, $\mu_p = 1000\,\mathrm{Pa\cdot s}$.
- **Fluide 2 (Fluide porteur, $C=0$)** : $\rho_f = 1000\,\mathrm{kg/m}^3$, $\mu_f = 0.153125\,\mathrm{Pa\cdot s}$, sous gravité $g_y = -9.81\,\mathrm{m/s}^2$.

### 2.2 Loi Analytique de Vitesse Terminale de Faxén
En régime de Stokes 2D confiné ($k_x = a/l = 0.125$), le développement asymptotique de Faxén donne la vitesse terminale analytique du cylindre ($g = |g_y| = 9.81\,\mathrm{m/s}^2$) :
$$U_{c,\mathrm{th}} = \frac{(\rho_f - \rho_p)\,g\,a^2}{4\,\mu_f}\,\left(-0.9157 - \ln(k_x) + 1.7244\,k_x^2 - 1.7302\,k_x^4\right) \approx -0.05060\,\mathrm{m/s}.$$
*(Note technique sur les diagnostics : dans `VTKWriter::write_diagnostics`, les colonnes `center_of_mass_y` et `rise_velocity` de `*_history.csv` intègrent l'indicatrice $1 - C$, qui correspond au fluide porteur lorsque la particule est définie par `"type": "DROPLET"`. La vitesse verticale de la particule $V_p(t)$ s'extrait donc directement par intégration de $C \cdot v$ sur `PhaseIndicator` et `Velocity` dans les fichiers `.vtr`.)*

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et confinement $k_x$ | $[0, 2l]^2$, $a$, $k_x = a/l$ | $2l = 0.1\,\mathrm{m}$ ($l=0.05\,\mathrm{m}$), $a = 0.00625\,\mathrm{m}$, $k_x = 0.125$, centre $(0.05, 0.05)\,\mathrm{m}$ (`axis_a = axis_b = 0.00625`, `pert_amp = 0.0`) |
| Particule ($C=1$) et Fluide ($C=0$) | $\rho_1, \mu_1$ ; $\rho_2, \mu_2$ | $\rho_1 = \rho_p = 1500\,\mathrm{kg/m}^3, \mu_1 = 1000\,\mathrm{Pa\cdot s}$ ; $\rho_2 = \rho_f = 1000\,\mathrm{kg/m}^3, \mu_2 = 0.153125\,\mathrm{Pa\cdot s}$ |
| Séquence de grilles et $\Delta t$ | $64\times 128$, $128\times 128$, $256\times 256$ | $\Delta t = 5\times 10^{-4}\,\mathrm{s}$, $t_{\mathrm{final}} = 0.15\,\mathrm{s}$ |
| Solveurs et filtrage SVOF | `solver` & `physics` | `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (`k_aal=100`), `AUGMENTED_LAGRANGIAN`, `MONOLITHIC` ; `enable_svof_smoothing: true` (`svof_filter_width: 0.5`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_9_sedimentation_faxen_n128.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/n128"`.)*

```json
{
  "case_name": "cas_9_sedimentation_faxen_n128",
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 0.1,
    "ly": 0.1
  },
  "physics": {
    "enable_two_phase": true,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "PICARD",
    "enable_surface_tension": false,
    "enable_svof_smoothing": true,
    "svof_filter_width": 0.5,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 1500.0, "mu": 1000.0 },
    "fluid2": { "rho": 1000.0, "mu": 0.153125 }
  },
  "initial_condition": {
    "type": "DROPLET",
    "center_x": 0.05,
    "center_y": 0.05,
    "radius": 0.00625,
    "axis_a": 0.00625,
    "axis_b": 0.00625,
    "pert_amp": 0.0
  },
  "time": {
    "dt": 0.0005,
    "t_final": 0.15,
    "time_scheme": "GEAR",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "right":  { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "ADAPTATIVE_AUGMENTED_LAGRANGIAN",
    "adaptative_al": true,
    "k_aal": 100.0,
    "uzawa_tol": 1e-7,
    "uzawa_iter": 6,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 800,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.3
  },
  "output": {
    "output_dir": "output/cpp/n128",
    "case_name": "cas_9_sedimentation_faxen_n128",
    "vtk_interval": 30,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs fonctionnent en **MPI pur** (`mpirun -np 4`).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/n128 output/fortran/n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_9_sedimentation_faxen_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_9_sedimentation_faxen_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Précision sur la vitesse terminale de Faxén :**
   À partir des fichiers `.pvtr`/`.vtr` sur le plateau quasi-stationnaire ($t \in [0.05, 0.10]\,\mathrm{s}$, avant l'approche de la paroi inférieure), calculer la vitesse moyenne verticale du cylindre :
   $$V_p(t) = \frac{\int_\Omega C(\mathbf{x}, t)\,v(\mathbf{x}, t)\,\mathrm{d}\Omega}{\int_\Omega C(\mathbf{x}, t)\,\mathrm{d}\Omega}, \qquad \mathrm{Err}_{\mathrm{rel}}(V_p) = \frac{|V_p - U_{c,\mathrm{th}}|}{|U_{c,\mathrm{th}}|}.$$
   Vérifier la convergence de $V_p$ vers $U_{c,\mathrm{th}} \approx -0.05060\,\mathrm{m/s}$ à mesure que le maillage est raffiné ($64\times 128 \to 256\times 256$).
2. **Rigidité de la particule et énergie cinétique :**
   Vérifier la conservation de la forme circulaire et du volume $V_p(t) = \int_\Omega C\,\mathrm{d}\Omega$ ($|\Delta V_p / V_p(0)| < 10^{-4}$) ainsi que la stabilisation de `kinetic_energy` dans `*_history.csv`.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier la concordance des courbes $V_p(t)$ et `kinetic_energy(t)` entre les deux exécutables.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` : `step,time,dt,center_of_mass_y,rise_velocity,bubble_volume,kinetic_energy,max_u,max_div` (`9` colonnes C++ / `11` colonnes Fortran avec `avg_rho_gas,avg_p_gas`).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `PhaseIndicator_Smoothed`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 16 (Sédimentation de Faxén) depuis les .pvtr."""
import math
import xml.etree.ElementTree as ET
from pathlib import Path

def particle_velocity_from_pvtr(pvtr_path: str):
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    num, den = 0.0, 0.0
    for piece in root.findall(".//Piece"):
        vtr = ET.parse(pvtr.parent / piece.attrib["Source"]).getroot()
        ext = [int(v) for v in piece.attrib["Extent"].split()]
        nx_l, ny_l = ext[1] - ext[0], ext[3] - ext[2]
        xc = [float(v) for v in vtr.find(".//DataArray[@Name='XCoordinates']").text.split()]
        yc = [float(v) for v in vtr.find(".//DataArray[@Name='YCoordinates']").text.split()]
        C = [float(v) for v in vtr.find(".//DataArray[@Name='PhaseIndicator']").text.split()]
        V = [float(v) for v in vtr.find(".//DataArray[@Name='Velocity']").text.split()]
        for j in range(ny_l):
            dy = yc[j + 1] - yc[j]
            for i in range(nx_l):
                dx = xc[i + 1] - xc[i]
                k = j * nx_l + i
                w = C[k] * dx * dy
                num += w * V[3 * k + 1]
                den += w
    return num / den, den

if __name__ == "__main__":
    kx = 0.00625 / 0.05
    u_faxen = ((1000.0 - 1500.0) * 9.81 * (0.00625 ** 2) / (4.0 * 0.153125)) * (
        -0.9157 - math.log(kx) + 1.7244 * kx**2 - 1.7302 * kx**4
    )
    for p in sorted(Path("output/cpp/n128").glob("*.pvtr")):
        vp, area = particle_velocity_from_pvtr(str(p))
        print(f"{p.name} | Vp={vp:.6f} m/s (Faxen={u_faxen:.6f}, err={abs(vp-u_faxen)/abs(u_faxen):.2%}) | Aire={area:.6e}")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Tracer l'isocontour `PhaseIndicator = 0.5` pour suivre le cylindre solide et appliquer le filtre **Stream Tracer** (ou **Glyph**) sur `Velocity` pour visualiser les deux tourbillons de recirculation générés entre la particule et les parois verticales.
# Tutoriel 17 : Rupture de Barrage 2D sur Fond Sec (Solution de Saint-Venant)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_5_rupture_barrage/` *(campagne avancée `Run_CASES`)*.

Ce tutoriel aborde l'effondrement gravitationnel instantané d'une colonne d'eau sur fond sec (*Dam Break*). La cinématique d'avancée du front d'onde et le profil parabolique de la surface libre sont confrontés à la solution analytique des caractéristiques de Saint-Venant (solution de Ritter).

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre la physique de l'effondrement d'une colonne liquide à grand rapport de densité ($\rho_1/\rho_2 \approx 850$) et le modèle analytique en eau peu profonde de Saint-Venant.
- Configurer la condition initiale rectangulaire `"DAM_BREAK"` (`box_x0`, `box_y0`, `box_x1`, `box_y1`) et les schémas `GEAR` + `ADAMS_BASHFORTH` + `TVD_SUPERBEE`.
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Extraire dynamiquement la hauteur d'eau intégrée $h_{\mathrm{num}}(x,t)$ et l'abscisse du front $x_{\mathrm{front}}(t)$ depuis les fichiers `.vtr`.
- Visualiser la propagation de la vague et les tourbillons d'air sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide le couplage Navier-Stokes diphasique à surface libre ($\rho_1/\rho_2 \approx 850$, $\mu_1/\mu_2 = 100$) sous gravité sur le problème classique de rupture de barrage 2D sur fond sec.

Dans un canal $\Omega = [0, 1.2] \times [0, 0.2]\,\mathrm{m}^2$ à parois glissantes (`FREE_SLIP`), une colonne d'eau occupe initialement $[0, 0.6] \times [0, 0.1]\,\mathrm{m}^2$ ($x_0 = 0.6\,\mathrm{m}, h_g = 0.1\,\mathrm{m}$) avec $\rho_1 = 1000\,\mathrm{kg/m}^3, \mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$, $\rho_2 = 1.1768\,\mathrm{kg/m}^3, \mu_2 = 10^{-5}\,\mathrm{Pa\cdot s}$, $\sigma = 0.075\,\mathrm{N/m}$ et $g_y = -9.81\,\mathrm{m/s}^2$.

### 2.2 Solution Analytique des Caractéristiques de Saint-Venant (Ritter)
Avec la célérité d'onde initiale $c_0 = \sqrt{g h_g} \approx 0.990454\,\mathrm{m/s}$ ($g = |g_y| = 9.81\,\mathrm{m/s}^2$), la hauteur d'eau analytique en théorie d'eau peu profonde (avant réflexion sur le mur droit $x = 1.2\,\mathrm{m}$ à $t_{\mathrm{refl}} = (L_x - x_0)/(2c_0) \approx 0.303\,\mathrm{s}$ pour le front et sans effet de bord amont) vaut :
$$h_{\mathrm{SV}}(x,t) = \begin{cases} h_g & \text{si } x \le x_0 - c_0 t, \\ \dfrac{1}{9g}\left(2c_0 - \dfrac{x - x_0}{t}\right)^2 & \text{si } x_0 - c_0 t < x < x_0 + 2c_0 t, \\ 0 & \text{si } x \ge x_0 + 2c_0 t. \end{cases}$$
L'abscisse théorique du front d'onde sec est $x_{\mathrm{front,SV}}(t) = x_0 + 2 c_0 t$. Numériquement, le profil de surface libre se calcule directement par intégration verticale de la fraction volumique d'eau :
$$h_{\mathrm{num}}(x_i, t) = \sum_{j=1}^{N_y} C_{i,j}(t)\,\Delta y.$$

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et retenue initiale | $[0, 1.2]\times[0, 0.2]\,\mathrm{m}^2$, $x_0, h_g$ | $x_0 = 0.6\,\mathrm{m}$, $h_g = 0.1\,\mathrm{m}$ (`DAM_BREAK`), $g_y = -9.81\,\mathrm{m/s}^2$ |
| Fluides 1 (Eau) et 2 (Air) | $\rho_1, \mu_1, \rho_2, \mu_2, \sigma$ | $\rho_1 = 1000\,\mathrm{kg/m}^3$, $\mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$ ; $\rho_2 = 1.1768\,\mathrm{kg/m}^3$, $\mu_2 = 10^{-5}\,\mathrm{Pa\cdot s}$ ; $\sigma = 0.075\,\mathrm{N/m}$ |
| Séquence de grilles et $\Delta t$ | $120\times 20$, $240\times 40$, $480\times 80$ | $\Delta t \in \{10^{-3}, 5\times 10^{-4}, 2.5\times 10^{-4}\}\,\mathrm{s}$, $t_{\mathrm{final}} = 0.50\,\mathrm{s}$ (contrôles à $0.24\,\mathrm{s}$ et $0.40\,\mathrm{s}$) |
| Schémas VOF et solveurs | `vof` & `solver` | `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC`, `VOF_TVD` (`SUPERBEE`) ; `ADAPTATIVE_AUGMENTED_LAGRANGIAN`, `AUGMENTED_LAGRANGIAN`, `MONOLITHIC` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_12_rupture_barrage_240x40.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/grid_240x40"`.)*

```json
{
  "case_name": "cas_12_rupture_barrage_240x40",
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 240,
    "ny": 40,
    "lx": 1.2,
    "ly": 0.2
  },
  "physics": {
    "enable_two_phase": true,
    "enable_convection": true,
    "convective_scheme": "TVD_SUPERBEE",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_surface_tension": true,
    "enable_diffusion_smoothing": true,
    "diffusion_smoothing_steps": 2,
    "sigma": 0.075,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 1000.0, "mu": 1.0e-3 },
    "fluid2": { "rho": 1.1768, "mu": 1.0e-5 }
  },
  "initial_condition": {
    "type": "DAM_BREAK",
    "box_x0": 0.0,
    "box_y0": 0.0,
    "box_x1": 0.6,
    "box_y1": 0.1
  },
  "time": {
    "dt": 0.0005,
    "t_final": 0.50,
    "time_scheme": "GEAR",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "top":    { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "ADAPTATIVE_AUGMENTED_LAGRANGIAN",
    "adaptative_al": true,
    "k_aal": 100.0,
    "uzawa_tol": 1e-6,
    "uzawa_iter": 5,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 800,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.3
  },
  "output": {
    "output_dir": "output/cpp/grid_240x40",
    "case_name": "cas_12_rupture_barrage_240x40",
    "vtk_interval": 40,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur** (`mpirun -np 4`).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/grid_240x40 output/fortran/grid_240x40

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_12_rupture_barrage_240x40.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_12_rupture_barrage_240x40.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Profil de surface libre $h_{\mathrm{num}}(x,t)$ et front d'onde $x_{\mathrm{front}}(t)$ :**
   À $t = 0.24\,\mathrm{s}$ (avant réflexion sur la paroi droite), extraire $h_{\mathrm{num}}(x_i, t) = \sum_{j=1}^{N_y} C_{i,j}(t)\Delta y$ depuis les `.vtr` et comparer à $h_{\mathrm{SV}}(x_i, t)$ en norme $L_2$. À $t = 0.40\,\mathrm{s}$, observer la remontée de vague sur la paroi droite.
2. **Conservation du volume d'eau et comparaison des schémas VOF :**
   Vérifier par sommation de `PhaseIndicator * dx * dy` que la variation relative de volume d'eau $|\Delta V / V_0| < 10^{-4}$ avec `WEYMOUTH_YUE+PLIC` et `VOF_TVD`.
3. **Performance de `AAL` vs `AL` et parité C++17 / Fortran 2018 :**
   Comparer la convergence Krylov et vérifier l'accord de $h_{\mathrm{num}}(x,t)$ et $E_k(t)$ entre `covof2d_cpp` et `covof2d_fortran`.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` (`9` colonnes en C++ / `11` colonnes en Fortran).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Curvature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 17 (Rupture de barrage / Saint-Venant)."""
import math
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

def h_saint_venant(x: float, t: float, x0: float = 0.6, hg: float = 0.1, g: float = 9.81) -> float:
    c0 = math.sqrt(g * hg)
    if x <= x0 - c0 * t:
        return hg
    if x < x0 + 2.0 * c0 * t:
        return ((2.0 * c0 - (x - x0) / t) ** 2) / (9.0 * g)
    return 0.0

def extract_height_profile(pvtr_path: str):
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    cols = defaultdict(float)
    vol_water = 0.0
    for piece in root.findall(".//Piece"):
        vtr = ET.parse(pvtr.parent / piece.attrib["Source"]).getroot()
        ext = [int(v) for v in piece.attrib["Extent"].split()]
        nx_l, ny_l = ext[1] - ext[0], ext[3] - ext[2]
        xc = [float(v) for v in vtr.find(".//DataArray[@Name='XCoordinates']").text.split()]
        yc = [float(v) for v in vtr.find(".//DataArray[@Name='YCoordinates']").text.split()]
        C = [float(v) for v in vtr.find(".//DataArray[@Name='PhaseIndicator']").text.split()]
        for j in range(ny_l):
            dy = yc[j + 1] - yc[j]
            for i in range(nx_l):
                dx = xc[i + 1] - xc[i]
                x = 0.5 * (xc[i] + xc[i + 1])
                cval = C[j * nx_l + i]
                cols[round(x, 8)] += cval * dy
                vol_water += cval * dx * dy
    return sorted(cols.items()), vol_water

if __name__ == "__main__":
    p = Path("output/cpp/grid_240x40/cas_12_rupture_barrage_240x40_000480.pvtr")
    if p.exists():
        prof, vol = extract_height_profile(str(p))
        l2 = math.sqrt(sum((h - h_saint_venant(x, 0.24)) ** 2 for x, h in prof) / len(prof))
        x_front = max((x for x, h in prof if h > 1e-3), default=0.6)
        print(f"t=0.24 s | L2(h)={l2:.4e} m | x_front={x_front:.4f} m (SV={0.6+2*math.sqrt(9.81*0.1)*0.24:.4f}) | V={vol:.6f}")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Tracer l'isocontour `PhaseIndicator = 0.5` sur fond coloré par la norme de `Velocity` afin de suivre le collapse gravitationnel, l'avancée du front sec et la réflexion contre la paroi aval.
3. Utiliser **Glyph** sur `Velocity` pour observer le contre-écoulement de l'air au-dessus de la nappe liquide.
# Tutoriel 18 : Oscillations Capillaires Amorties d'une Bulle de Gaz (Théorie de Lamb)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_6_oscillation_bulle_capillaire/` *(campagne avancée `Run_CASES`)*.

Ce tutoriel met en œuvre une bulle de gaz initialement déformée en ellipse selon le mode quadrupolaire $n = 2$ dans un métal liquide dense, sous l'action combinée de la tension superficielle et de la dissipation visqueuse. La fréquence propre d'oscillation et l'amortissement sont confrontés à la théorie linéaire de Rayleigh-Lamb.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre la dynamique inertio-capillaire d'une bulle oscillante à fort rapport de masse volumique ($\rho_1/\rho_2 \approx 5950$).
- Configurer une bulle elliptique explicite (`"axis_a"`, `"axis_b"`, `"pert_amp": 0.0`) et la régularisation CSF (`enable_diffusion_smoothing`).
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Extraire dynamiquement les demi-axes $a(t), b(t)$, la fréquence d'oscillation $f_{\mathrm{num}}$ et le bilan énergétique depuis les `.vtr` et `*_history.csv`.
- Visualiser les cycles d'oscillation et la pression capillaire sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide la dynamique instationnaire du couplage entre l'inertie diphasique à fort rapport de densité ($\rho_1/\rho_2 \approx 5950$), la dissipation visqueuse et la force de tension superficielle CSF lors des oscillations libres d'une bulle initialement elliptique.

Dans une cavité carrée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ sans gravité ($\mathbf{g} = \mathbf{0}$) avec parois à glissement libre (`FREE_SLIP`), une bulle d'air ($\rho_2 = 1.17683\,\mathrm{kg/m}^3$, $\mu_2 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$) est centrée en $(0.05, 0.05)\,\mathrm{m}$ dans un métal liquide dense ($\rho_1 = 7000\,\mathrm{kg/m}^3$, $\mu_1 = 4.0\times 10^{-3}\,\mathrm{Pa\cdot s}$) avec $\sigma = 1.8\,\mathrm{N/m}$. À $t=0$, la bulle de rayon équivalent $R_0 = 0.02\,\mathrm{m}$ est déformée sur le mode $n=2$ avec une excentricité $\varepsilon = 0.20$ :
$$a_0 = R_0 (1 + \varepsilon) = 0.024\,\mathrm{m}, \qquad b_0 = \frac{R_0}{1 + \varepsilon} \approx 0.0166666667\,\mathrm{m} \qquad (a_0 b_0 = R_0^2).$$

### 2.2 Relation de Dispersion Capillaire Linéaire
En milieu infini non visqueux, la pulsation propre $\omega_n$ et la fréquence $f_{\mathrm{th}}$ du mode $n = 2$ de Rayleigh-Lamb cylindrique valent :
$$\omega_n^2 = \frac{(n^3 - n)\,\sigma}{(\rho_1 + \rho_2)\,R_0^3} \implies \omega_2 = \sqrt{\frac{6\,\sigma}{(\rho_1 + \rho_2)\,R_0^3}}, \qquad f_{\mathrm{th}} = \frac{\omega_2}{2\pi} \approx 5.70632\,\mathrm{Hz} \quad (T_{\mathrm{th}} \approx 0.17524\,\mathrm{s}).$$
En cavité carrée confinée ($L / R_0 = 5$) et pour une amplitude finie $\varepsilon = 0.20$, l'effet de confinement pariétal abaisse légèrement la fréquence d'oscillation effective par rapport à la valeur en domaine infini, tandis que l'énergie cinétique $E_k(t)$ oscille à la fréquence double $2 f_{\mathrm{num}}$.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et bulle initiale | $[0, L]^2$, $(x_c, y_c)$, $a_0, b_0$ | $L = 0.1\,\mathrm{m}$, centre $(0.05, 0.05)\,\mathrm{m}$, $R_0 = 0.02\,\mathrm{m}$, $a_0 = 0.024\,\mathrm{m}$, $b_0 = 0.0166666667\,\mathrm{m}$ (`pert_amp = 0.0`) |
| Fluide 1 (Métal liquide) | $\rho_1, \mu_1$ | $\rho_1 = 7000.0\,\mathrm{kg/m}^3$, $\mu_1 = 4.0\times 10^{-3}\,\mathrm{Pa\cdot s}$ |
| Fluide 2 (Air) | $\rho_2, \mu_2$ | $\rho_2 = 1.17683\,\mathrm{kg/m}^3$, $\mu_2 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$ |
| Tension superficielle et lissage | $\sigma$, `enable_diffusion_smoothing` | $\sigma = 1.8\,\mathrm{N/m}$, `enable_diffusion_smoothing: true`, `diffusion_smoothing_steps: 3` |
| Séquence de grilles et $\Delta t$ | $N \in \{32, 64, 128, 256, 512\}$ | $\Delta t = 2.5\times 10^{-4}\,\mathrm{s}$ (ou adapté selon $\Delta t_{\mathrm{cap}}$ sur grille fine), $t_{\mathrm{final}} = 0.6\,\mathrm{s}$ |
| Solveurs et schémas VOF | `solver` & `vof` | `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (`AAL`), `AUGMENTED_LAGRANGIAN`, `MONOLITHIC` ; `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_13_oscillation_bulle_n128.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/n128"`.)*

```json
{
  "case_name": "cas_13_oscillation_bulle_n128",
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 0.1,
    "ly": 0.1
  },
  "physics": {
    "enable_two_phase": true,
    "compressible": false,
    "enable_convection": true,
    "enable_surface_tension": true,
    "enable_diffusion_smoothing": true,
    "diffusion_smoothing_steps": 3,
    "sigma": 1.8,
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": { "rho": 7000.0,  "mu": 4.0e-3 },
    "fluid2": { "rho": 1.17683, "mu": 1.85e-5 }
  },
  "initial_condition": {
    "type": "BUBBLE",
    "center_x": 0.05,
    "center_y": 0.05,
    "radius": 0.02,
    "axis_a": 0.024,
    "axis_b": 0.0166666667,
    "pert_amp": 0.0
  },
  "time": {
    "dt": 0.00025,
    "t_final": 0.6,
    "time_scheme": "GEAR",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "top":    { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "ADAPTATIVE_AUGMENTED_LAGRANGIAN",
    "adaptative_al": true,
    "k_aal": 100.0,
    "uzawa_tol": 1e-6,
    "uzawa_iter": 6,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 600,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.3
  },
  "output": {
    "output_dir": "output/cpp/n128",
    "case_name": "cas_13_oscillation_bulle_n128",
    "vtk_interval": 20,
    "log_interval": 5,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur** (`mpirun -np 4`).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/n128 output/fortran/n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_13_oscillation_bulle_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_13_oscillation_bulle_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Fréquence d'oscillation $f_{\mathrm{num}}$ et amortissement visqueux :**
   Extraire la série temporelle d'énergie cinétique $E_k(t)$ depuis `*_history.csv` et les demi-axes $a(t), b(t)$ depuis les fichiers `.vtr` (à partir des moments géométriques d'ordre 2 de $1 - C$). Déterminer la période d'oscillation $T_{\mathrm{num}}$ par détection des maxima de $a(t)$ (ou à partir de la demi-période des pics de $E_k(t)$) et vérifier la convergence de $f_{\mathrm{num}} = 1/T_{\mathrm{num}}$ avec le raffinement du maillage.
2. **Conservation du volume de la bulle et stabilité à haut rapport de densité ($\rho_1/\rho_2 \approx 5950$) :**
   Vérifier par sommation de $(1-C)\Delta x\Delta y$ dans les `.vtr` que $|V(t_{\mathrm{final}}) - V(0)| / V(0) < 10^{-4}$ avec `WEYMOUTH_YUE+PLIC`.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier la superposition des courbes $E_k(t)$, $a(t)$ et $b(t)$ entre `covof2d_cpp` et `covof2d_fortran`.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` (`9` colonnes en C++ / `11` colonnes en Fortran).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Curvature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 18 (Oscillation de bulle / Lamb)."""
import csv
import math
import xml.etree.ElementTree as ET
from pathlib import Path

def estimate_freq_from_ke(csv_path: str) -> float:
    t_vals, ke_vals = [], []
    with Path(csv_path).open(newline="") as f:
        for r in csv.DictReader(f):
            t_vals.append(float(r["time"]))
            ke_vals.append(float(r["kinetic_energy"]))
    peaks = [t_vals[i] for i in range(1, len(ke_vals) - 1)
             if ke_vals[i] > ke_vals[i - 1] and ke_vals[i] >= ke_vals[i + 1] and ke_vals[i] > 1e-8]
    if len(peaks) < 3:
        return float("nan")
    dt_peaks = [(peaks[i + 1] - peaks[i]) for i in range(len(peaks) - 1)]
    t_half = sum(dt_peaks) / len(dt_peaks)
    return 1.0 / (2.0 * t_half)  # E_k oscille à 2*f

def bubble_axes_from_pvtr(pvtr_path: str):
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    m0, mxx, myy = 0.0, 0.0, 0.0
    for piece in root.findall(".//Piece"):
        vtr = ET.parse(pvtr.parent / piece.attrib["Source"]).getroot()
        ext = [int(v) for v in piece.attrib["Extent"].split()]
        nx_l, ny_l = ext[1] - ext[0], ext[3] - ext[2]
        xc = [float(v) for v in vtr.find(".//DataArray[@Name='XCoordinates']").text.split()]
        yc = [float(v) for v in vtr.find(".//DataArray[@Name='YCoordinates']").text.split()]
        C = [float(v) for v in vtr.find(".//DataArray[@Name='PhaseIndicator']").text.split()]
        for j in range(ny_l):
            y = 0.5 * (yc[j] + yc[j + 1])
            dy = yc[j + 1] - yc[j]
            for i in range(nx_l):
                x = 0.5 * (xc[i] + xc[i + 1])
                dx = xc[i + 1] - xc[i]
                w = (1.0 - C[j * nx_l + i]) * dx * dy
                m0 += w
                mxx += w * (x - 0.05) ** 2
                myy += w * (y - 0.05) ** 2
    return 2.0 * math.sqrt(mxx / m0), 2.0 * math.sqrt(myy / m0), m0

if __name__ == "__main__":
    f_th = math.sqrt(6.0 * 1.8 / ((7000.0 + 1.17683) * 0.02 ** 3)) / (2.0 * math.pi)
    csv_f = Path("output/cpp/n128/cas_13_oscillation_bulle_n128_history.csv")
    if csv_f.exists():
        f_num = estimate_freq_from_ke(str(csv_f))
        print(f"Fréquence numérique f_num = {f_num:.4f} Hz (théorie milieu infini = {f_th:.4f} Hz)")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Appliquer le filtre **Contour** sur `PhaseIndicator = 0.5` et colorer le fond par `Pressure` pour observer l'inversion alternée de la surpression aux pôles horizontaux puis verticaux.
3. Afficher les vecteurs **Glyph** de `Velocity` pour visualiser l'écoulement quadripolaire interne et externe.
# Tutoriel 19 : Instabilité Interfaciale Diphasique de Rayleigh-Taylor (Dispersion Linéaire)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_7_instabilite_interface_gravitaire/` *(campagne avancée `Run_CASES`)*.

Ce dernier tutoriel simule l'instabilité hydrodynamique de Rayleigh-Taylor se développant à l'interface entre un fluide lourd superposé à un fluide léger soumis à la gravité. Le taux de croissance exponentiel en régime linéaire et le seuil de stabilisation capillaire sont confrontés à la relation de dispersion analytique hydro-capillaire.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre la compétition entre déstabilisation gravitationnelle (nombre d'Atwood $\mathrm{At}$) et stabilisation capillaire (seuil critique $\sigma_c$).
- Configurer la condition initiale `"RAYLEIGH_TAYLOR"` (`center_y`, `pert_amp`, `pert_k`) et les schémas d'ordre 2 (`GEAR`, `ADAMS_BASHFORTH`, `WEYMOUTH_YUE+PLIC`).
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Extraire dynamiquement l'amplitude de la perturbation $A(t)$ depuis les `.vtr` et estimer $\omega_{\mathrm{num}}$ par ajustement sur $A_0 \cosh(\omega t)$.
- Visualiser le développement du champignon non-linéaire sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide la compétition entre l'instabilité gravitationnelle de Rayleigh-Taylor (nombre d'Atwood $\mathrm{At} = (\rho_1 - \rho_2)/(\rho_1 + \rho_2) = 0.5$) et la stabilisation capillaire CSF, du régime linéaire exponentiel jusqu'à l'enroulement non-linéaire en champignon.

Dans un domaine $\Omega = [0, 1.0] \times [0, 4.0]\,\mathrm{m}^2$ avec parois latérales `FREE_SLIP` (ou `PERIODIC`) et parois haute/basse `DIRICHLET` adhérentes ($\mathbf{u} = \mathbf{0}$), le fluide lourd ($\rho_1 = 3.0\,\mathrm{kg/m}^3, \mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$, $C=1$ pour $y \ge y_{\mathrm{int}}$) surmonte initialement le fluide léger ($\rho_2 = 1.0\,\mathrm{kg/m}^3, \mu_2 = 10^{-3}\,\mathrm{Pa\cdot s}$, $C=0$ pour $y < y_{\mathrm{int}}$) sous $g_y = -9.81\,\mathrm{m/s}^2$. L'interface initiale est perturbée selon (`"type": "RAYLEIGH_TAYLOR"`, `"center_y": 2.0`, `"pert_amp": 0.05`, `"pert_k": 1.0`) :
$$y_{\mathrm{int}}(x, 0) = 2.0 + A_0 \cos\left(\frac{2\pi m x}{L_x}\right), \qquad k = \frac{2\pi m}{L_x} = 2\pi\,\mathrm{m}^{-1}, \quad A_0 = 0.05\,\mathrm{m}.$$

### 2.2 Relation de Dispersion Hydro-Capillaire Linéaire
En régime linéaire non visqueux, le taux de croissance théorique $\omega_{\mathrm{th}}(k)$ et le seuil de coupure capillaire $\sigma_c$ ($g = |g_y| = 9.81\,\mathrm{m/s}^2$) valent :
$$\omega_{\mathrm{th}}^2(k) = \frac{(\rho_1 - \rho_2)\,g\,k - \sigma\,k^3}{\rho_1 + \rho_2}, \qquad \sigma_c = \frac{(\rho_1 - \rho_2)\,g}{k^2} \approx 0.4970\,\mathrm{N/m}.$$
- Pour $\sigma = 0\,\mathrm{N/m}$ : régime instable pur avec $\omega_{\mathrm{th}}(0) \approx 5.5515\,\mathrm{s}^{-1}$.
- Pour $\sigma = 0.20\,\mathrm{N/m}$ : croissance ralentie par la capillarité avec $\omega_{\mathrm{th}}(0.2) \approx 4.2915\,\mathrm{s}^{-1}$.
- Pour $\sigma = 0.60\,\mathrm{N/m} > \sigma_c$ : $\omega_{\mathrm{th}}^2 < 0$, l'interface est stabilisée et oscille de façon dispersive sans croître.

En régime linéaire ($k A(t) \ll 1$), l'amplitude de la perturbation depuis le repos évolue selon $A(t) = A_0 \cosh(\omega t)$.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et perturbation initiale | $[0, 1]\times[0, 4]\,\mathrm{m}^2$, $y_c, A_0, m$ | $y_c = 2.0\,\mathrm{m}$, $A_0 = 0.05\,\mathrm{m}$, $m = 1$ ($k = 2\pi\,\mathrm{m}^{-1}$, `RAYLEIGH_TAYLOR`) |
| Fluides ($\mathrm{At} = 0.5$) et gravité | $\rho_1, \rho_2, \mu_1, \mu_2, g_y$ | $\rho_1 = 3.0\,\mathrm{kg/m}^3$, $\rho_2 = 1.0\,\mathrm{kg/m}^3$, $\mu_1 = \mu_2 = 10^{-3}\,\mathrm{Pa\cdot s}$, $g_y = -9.81\,\mathrm{m/s}^2$ |
| Trois régimes capillaires $\sigma$ | $\sigma \in \{0.0, 0.20, 0.60\}\,\mathrm{N/m}$ | Instable pur ($\sigma=0$), ralenti ($\sigma=0.20$), stabilisé ($\sigma=0.60 > \sigma_c$) |
| Séquence de grilles et $\Delta t$ | $64\times 256$, $128\times 512$, $256\times 1024$ | $\Delta t = 10^{-3}\,\mathrm{s}$, $t_{\mathrm{final}} = 1.0\,\mathrm{s}$ |
| Schémas VOF et solveurs | `vof` & `solver` | `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC`, `VOF_TVD` ; `ADAPTATIVE_AUGMENTED_LAGRANGIAN`, `AUGMENTED_LAGRANGIAN`, `MONOLITHIC` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_14_rayleigh_taylor_sigma0_128x512.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/sigma0_128x512"`. Pour les cas $\sigma > 0$, passer `"enable_surface_tension": true` et `"enable_diffusion_smoothing": true`.)*

```json
{
  "case_name": "cas_14_rayleigh_taylor_sigma0_128x512",
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 512,
    "lx": 1.0,
    "ly": 4.0
  },
  "physics": {
    "enable_two_phase": true,
    "enable_convection": true,
    "convective_scheme": "TVD_SUPERBEE",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_surface_tension": false,
    "sigma": 0.0,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 3.0, "mu": 1.0e-3 },
    "fluid2": { "rho": 1.0, "mu": 1.0e-3 }
  },
  "initial_condition": {
    "type": "RAYLEIGH_TAYLOR",
    "center_y": 2.0,
    "pert_amp": 0.05,
    "pert_k": 1.0
  },
  "time": {
    "dt": 0.001,
    "t_final": 1.0,
    "time_scheme": "GEAR",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "ADAPTATIVE_AUGMENTED_LAGRANGIAN",
    "adaptative_al": true,
    "k_aal": 100.0,
    "uzawa_tol": 1e-6,
    "uzawa_iter": 4,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 800,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.3
  },
  "output": {
    "output_dir": "output/cpp/sigma0_128x512",
    "case_name": "cas_14_rayleigh_taylor_sigma0_128x512",
    "vtk_interval": 50,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur** (`mpirun -np 4` ou plus sur maillage fin $256\times 1024$).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/sigma0_128x512 output/fortran/sigma0_128x512

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_14_rayleigh_taylor_sigma0_128x512.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_14_rayleigh_taylor_sigma0_128x512.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Taux de croissance linéaire $\omega_{\mathrm{num}}$ et seuil capillaire $\sigma_c$ :**
   À partir des fichiers `.pvtr`/`.vtr` sur la fenêtre linéaire $t \in [0.05, 0.30]\,\mathrm{s}$, calculer la demi-amplitude crête-à-creux $A(t) = \frac{1}{2}\big(y_{\mathrm{int},\max}(t) - y_{\mathrm{int},\min}(t)\big)$ de l'isocontour $C = 0.5$ et estimer $\omega_{\mathrm{num}}$ par ajustement sur $A(t) = A_0 \cosh(\omega_{\mathrm{num}} t)$. Comparer avec $\omega_{\mathrm{th}}(\sigma)$ pour $\sigma \in \{0.0, 0.20\}\,\mathrm{N/m}$, et vérifier la stabilisation oscillatoire pour $\sigma = 0.60 > \sigma_c$.
2. **Symétrie bilatérale et conservation de masse :**
   Vérifier la symétrie axiale $C(L_x - x, y, t) = C(x, y, t)$ tout au long de l'enroulement non-linéaire et vérifier que la variation relative de volume $|\Delta V / V_0| < 10^{-4}$ depuis `PhaseIndicator` dans les `.vtr`.
3. **Parité C++17 / Fortran 2018 :**
   Comparer $A(t)$, `kinetic_energy` (`*_history.csv`) et la morphologie des panaches entre `covof2d_cpp` et `covof2d_fortran`.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` (`9` colonnes en C++ / `11` colonnes en Fortran).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Curvature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 19 (Rayleigh-Taylor) depuis les .pvtr."""
import math
import xml.etree.ElementTree as ET
from pathlib import Path

def interface_amplitude_and_volume(pvtr_path: str):
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    y_mix, vol_heavy = [], 0.0
    for piece in root.findall(".//Piece"):
        vtr = ET.parse(pvtr.parent / piece.attrib["Source"]).getroot()
        ext = [int(v) for v in piece.attrib["Extent"].split()]
        nx_l, ny_l = ext[1] - ext[0], ext[3] - ext[2]
        xc = [float(v) for v in vtr.find(".//DataArray[@Name='XCoordinates']").text.split()]
        yc = [float(v) for v in vtr.find(".//DataArray[@Name='YCoordinates']").text.split()]
        C = [float(v) for v in vtr.find(".//DataArray[@Name='PhaseIndicator']").text.split()]
        for j in range(ny_l):
            y = 0.5 * (yc[j] + yc[j + 1])
            dy = yc[j + 1] - yc[j]
            for i in range(nx_l):
                dx = xc[i + 1] - xc[i]
                cval = C[j * nx_l + i]
                vol_heavy += cval * dx * dy
                if 0.05 < cval < 0.95:
                    y_mix.append(y)
    amp = 0.5 * (max(y_mix) - min(y_mix)) if y_mix else 0.0
    return amp, vol_heavy

if __name__ == "__main__":
    k = 2.0 * math.pi
    omega_th = math.sqrt((3.0 - 1.0) * 9.81 * k / (3.0 + 1.0))
    dt = 0.001
    for pvtr in sorted(Path("output/cpp/sigma0_128x512").glob("*.pvtr")):
        t = int(pvtr.stem.split("_")[-1]) * dt
        amp, vol = interface_amplitude_and_volume(str(pvtr))
        omega_est = math.acosh(max(1.0, amp / 0.05)) / t if 0.05 <= t <= 0.30 else float("nan")
        print(f"t={t:.3f} s | A={amp:.4f} m | omega_est={omega_est:.3f} s^-1 (th={omega_th:.3f}) | V1={vol:.6f}")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Colorer par `PhaseIndicator` et tracer l'isocontour `PhaseIndicator = 0.5` pour observer le développement symétrique de la pointe descendante (*spike*) et de la bulle ascendante (*bubble*) avec enroulement de Kelvin-Helmholtz aux bords.
3. Utiliser le filtre **Compute Derivatives** (vorticité) ou **Glyph** sur `Velocity` pour visualiser les deux tourbillons contrarotatifs qui alimentent le champignon.
# Tutoriel 20 : Compression d'Air par Injection d'Eau (Piston Compressible Bas Mach)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_8_piston_compressible/` *(campagne avancée `Run_CASES`)*.

Ce tutoriel valide la formulation compressible à bas nombre de Mach (Lagrangien Augmenté Compressible `C-AAL` / `MONOLITHIC` compressible) sur un cas d'injection d'un liquide quasi-incompressible comprimant une poche de gaz dans une enceinte étanche. La montée en pression et en masse volumique est confrontée à la loi analytique exacte issue de la conservation de la masse.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les équations 1-fluide faiblement compressibles à bas nombre de Mach et l'équation de continuité avec compressibilité isotherme $\chi_T(C)$.
- Configurer correctement les paramètres EOS (`"eos"`, `"chi_T"` avec majuscule) et l'état initial compressible (`"p_init"`, `"u"`, `"v"`).
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Quantifier dynamiquement $\overline{\rho}_{\mathrm{air}}(t)$ et $\overline{p}_{\mathrm{air}}(t)$ face à la loi analytique de compression.
- Visualiser le mouvement de l'interface plane et la montée en pression homogène sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `EOS`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_eos`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide l'extension du modèle diphasique 1-fluide aux écoulements faiblement compressibles à bas nombre de Mach (`"compressible": true`), dans lesquels un liquide quasi-incompressible comprime un gaz compressible au sein d'une cavité fermée.

Le système couple l'équation de continuité compressible isotherme $\chi_T(C) = C\,\chi_{T,1} + (1-C)\,\chi_{T,2}$, l'équation de quantité de mouvement, l'advection de la fonction couleur avec terme de dilatation $C\,\nabla\cdot\mathbf{u}$ et la loi d'état linéaire (`"eos": "LINEAR"`) :
$$\chi_T(C)\frac{\partial p}{\partial t} + \nabla \cdot \mathbf{u} = 0, \qquad \rho(C)\left(\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u}\right) = -\nabla p + \nabla\cdot\left[\mu(C)\left(\nabla\mathbf{u} + \nabla\mathbf{u}^T - \frac{2}{3}(\nabla\cdot\mathbf{u})\mathbf{I}\right)\right],$$
$$\frac{\partial C}{\partial t} + \nabla\cdot(\mathbf{u}C) = C\,\nabla\cdot\mathbf{u}, \qquad \rho_i(p) = \rho_{i,0}\left[1 + \chi_{T,i}\,(p - p_0)\right].$$

### 2.2 Configuration Géométrique et Loi Analytique Exacte
Dans une cavité carrée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ sans gravité ($\mathbf{g}=\mathbf{0}$) ni capillarité ($\sigma=0$), l'eau ($C=1$) occupe initialement $y \in [0, 0.025]\,\mathrm{m}$, surmontée d'une hauteur $H_0 = 0.075\,\mathrm{m}$ d'air ($C=0$) à $p_0 = 101\,325\,\mathrm{Pa}$. De l'eau est injectée par la face inférieure $y=0$ à la vitesse verticale constante $V_0 = 0.1\,\mathrm{m/s}$, les trois autres parois étant à glissement libre (`FREE_SLIP`). L'interface plane monte selon $y_{\mathrm{int}}(t) = 0.025 + V_0 t$, et la masse volumique ainsi que la pression de l'air vérifient la loi exacte :
$$\rho_{\mathrm{air}}(t) = \frac{\rho_{\mathrm{air},0}}{1 - \dfrac{V_0\,t}{H_0}}, \qquad p_{\mathrm{isotherme}}(t) = \frac{p_0}{1 - \dfrac{V_0\,t}{H_0}}.$$
*(Pour la loi d'état linéaire $\rho_2(p) = \rho_{2,0}[1 + \chi_{T,2}(p - p_0)]$ avec $\chi_{T,2} = 1/p_0 = 9.8692327\times 10^{-6}\,\mathrm{Pa}^{-1}$, on retrouve exactement $\rho_{\mathrm{air}}(t)$ et $p(t) = p_0 + \frac{1}{\chi_{T,2}}\frac{V_0 t / H_0}{1 - V_0 t / H_0} = \frac{p_0}{1 - V_0 t / H_0}$.)*

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et interface initiale | $[0, L_x]\times[0, L_y]$, $h_0$, $H_0$ | $L_x = L_y = 0.1\,\mathrm{m}$, eau sur $y\in[0, 0.025]\,\mathrm{m}$, air sur $H_0 = 0.075\,\mathrm{m}$ |
| Fluide 1 (Eau liquide) | $\rho_{1,0}, \mu_1, \chi_{T,1}$, `eos` | $\rho_{1,0} = 1000.0\,\mathrm{kg/m}^3$, $\mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$, $\chi_{T,1} = 4.44\times 10^{-10}\,\mathrm{Pa}^{-1}$, `"LINEAR"` |
| Fluide 2 (Air compressible) | $\rho_{2,0}, \mu_2, \chi_{T,2}$, `eos` | $\rho_{2,0} = 1.176829268\,\mathrm{kg/m}^3$, $\mu_2 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$, $\chi_{T,2} = 9.8692327\times 10^{-6}\,\mathrm{Pa}^{-1}$, `"LINEAR"` |
| Conditions aux limites | Bas ($y=0$), Haut/Gauche/Droite | Injection `DIRICHLET` $(u=0, v=V_0=0.1\,\mathrm{m/s})$ en bas ; `FREE_SLIP` sur les 3 autres bords |
| Séquence de grilles et $\Delta t$ | $N \in \{32, 64, 128, 256, 512\}$ | $\Delta t = 10^{-3}\,\mathrm{s}$, $t_{\mathrm{final}} = 0.5\,\mathrm{s}$ (500 pas ; contrôles à $t \in \{0.1, 0.2, 0.3, 0.4, 0.5\}\,\mathrm{s}$) |
| Solveurs et schémas VOF | `solver` & `vof` | `MONOLITHIC` ($r_{\mathrm{al}}=0$) ou `C-AAL` avec `BICGSTAB2` / `GMRES` / `DIRECT` ; `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC`, `VOF_TVD` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_10_piston_compressible_n128.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/n128"`. Respecter strictement la casse `"chi_T"`.)*

```json
{
  "case_name": "cas_10_piston_compressible_n128",
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 0.1,
    "ly": 0.1
  },
  "physics": {
    "enable_two_phase": true,
    "compressible": true,
    "flow_regime": "COMPRESSIBLE",
    "enable_convection": true,
    "enable_surface_tension": false,
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": { "rho": 1000.0,      "mu": 1.0e-3, "eos": "LINEAR", "chi_T": 4.44e-10 },
    "fluid2": { "rho": 1.176829268, "mu": 1.85e-5, "eos": "LINEAR", "chi_T": 9.8692327e-6 }
  },
  "initial_condition": {
    "type": "DAM_BREAK",
    "box_x0": 0.0,
    "box_y0": 0.0,
    "box_x1": 0.1,
    "box_y1": 0.025,
    "p_init": 101325.0,
    "u": 0.0,
    "v": 0.0
  },
  "time": {
    "dt": 0.001,
    "t_final": 0.5,
    "time_scheme": "EULER",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.1 },
    "top":    { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "MONOLITHIC",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-10,
    "max_iter": 1000,
    "r_al": 0.0,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.5
  },
  "output": {
    "output_dir": "output/cpp/n128",
    "case_name": "cas_10_piston_compressible_n128",
    "vtk_interval": 50,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur** (`mpirun -np 4` avec solveurs itératifs Krylov ; `mpirun -np 1` si `DIRECT` est sélectionné).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/n128 output/fortran/n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_10_piston_compressible_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_10_piston_compressible_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Précision sur $\overline{\rho}_{\mathrm{air}}(t)$ et $\overline{p}_{\mathrm{air}}(t)$ :**
   À partir des fichiers `.pvtr`/`.vtr` aux instants $t \in \{0.1, 0.2, 0.3, 0.4, 0.5\}\,\mathrm{s}$, extraire la masse volumique et la pression moyennes dans la zone gazeuse pure $\{C < 0.01\}$ et calculer l'écart relatif par rapport à $\rho_{\mathrm{th}}(t)$ et $p_{\mathrm{th}}(t)$ :
   $$E_{\mathrm{rel}}(\rho, t) = \frac{|\overline{\rho}_{\mathrm{air}}(t) - \rho_{\mathrm{th}}(t)|}{\rho_{\mathrm{th}}(t)}, \qquad E_{\mathrm{rel}}(p, t) = \frac{|\overline{p}_{\mathrm{air}}(t) - p_{\mathrm{th}}(t)|}{p_{\mathrm{th}}(t)}.$$
   Vérifier la convergence de ces erreurs vers zéro lors du raffinement du maillage ($32^2 \to 512^2$).
2. **Conservation de la masse totale de gaz et planéité de l'interface :**
   Vérifier que la masse totale d'air $M_{\mathrm{air}}(t) = \int_\Omega (1 - C)\,\rho\,\mathrm{d}\Omega$ reste constante ($|M_{\mathrm{air}}(t) - M_{\mathrm{air}}(0)| / M_{\mathrm{air}}(0) < 10^{-3}$) et que le front liquide-gaz reste horizontal à $y_{\mathrm{int}}(t) = 0.025 + V_0 t$.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier la concordance des courbes $\overline{\rho}_{\mathrm{air}}(t)$ et $\overline{p}_{\mathrm{air}}(t)$ entre `covof2d_cpp` et `covof2d_fortran` ($|\overline{p}^{\mathrm{C++}} - \overline{p}^{\mathrm{Fortran}}| / p_{\mathrm{th}} < 10^{-6}$).

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` (`9` colonnes en C++ ; `11` colonnes en Fortran incluant `avg_rho_gas,avg_p_gas`).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 20 (Piston compressible) sur fichiers .pvtr."""
import xml.etree.ElementTree as ET
from pathlib import Path

def eval_gas_state(pvtr_path: str):
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    sum_rho, sum_p, sum_vol, sum_mass_air = 0.0, 0.0, 0.0, 0.0
    for piece in root.findall(".//Piece"):
        vtr = ET.parse(pvtr.parent / piece.attrib["Source"]).getroot()
        ext = [int(v) for v in piece.attrib["Extent"].split()]
        nx_l, ny_l = ext[1] - ext[0], ext[3] - ext[2]
        xc = [float(v) for v in vtr.find(".//DataArray[@Name='XCoordinates']").text.split()]
        yc = [float(v) for v in vtr.find(".//DataArray[@Name='YCoordinates']").text.split()]
        C = [float(v) for v in vtr.find(".//DataArray[@Name='PhaseIndicator']").text.split()]
        P = [float(v) for v in vtr.find(".//DataArray[@Name='Pressure']").text.split()]
        R = [float(v) for v in vtr.find(".//DataArray[@Name='Density']").text.split()]
        for j in range(ny_l):
            dy = yc[j + 1] - yc[j]
            for i in range(nx_l):
                dx = xc[i + 1] - xc[i]
                dvol = dx * dy
                k = j * nx_l + i
                sum_mass_air += (1.0 - C[k]) * R[k] * dvol
                if C[k] < 0.01:
                    sum_rho += R[k] * dvol
                    sum_p += P[k] * dvol
                    sum_vol += dvol
    return sum_rho / sum_vol, sum_p / sum_vol, sum_mass_air

if __name__ == "__main__":
    rho0, p0, v0, h0 = 1.176829268, 101325.0, 0.1, 0.075
    dt = 0.001
    for pvtr in sorted(Path("output/cpp/n128").glob("*.pvtr")):
        step = int(pvtr.stem.split("_")[-1])
        t = step * dt
        rho_m, p_m, m_air = eval_gas_state(str(pvtr))
        factor = 1.0 / (1.0 - v0 * t / h0)
        rho_th, p_th = rho0 * factor, p0 * factor
        print(f"t={t:.3f} s | rho={rho_m:.6f} (th={rho_th:.6f}, err={abs(rho_m-rho_th)/rho_th:.3e}) "
              f"| p={p_m:.1f} (th={p_th:.1f}) | M_air={m_air:.6e}")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Afficher `PhaseIndicator` avec un contour à $0.5$ pour vérifier la planéité parfaite de la surface libre ascendante.
3. Tracer un profil vertical via **Plot Over Line** ($x = 0.05\,\mathrm{m}$) sur `Pressure` et `Density` pour observer la pression uniforme dans la poche d'air et le saut de masse volumique à l'interface.
# Tutoriel 21 : Compression Isotrope d'une Bulle d'Air en Cavité Carrée

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_9_compression_bulle/` *(campagne avancée `Run_CASES`)*.

Ce tutoriel évalue la capacité de l'algorithme bas-Mach compressible (`C-AAL` / `MONOLITHIC` compressible) à traiter la compression multidimensionnelle d'une bulle de gaz immergée sous injection symétrique sur les quatre parois. L'évolution temporelle du rayon équivalent de la bulle, de sa masse volumique et de la pression est confrontée à la solution analytique exacte déduite du bilan de masse global.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre la dynamique de compression d'une bulle gazeuse dans un liquide faiblement compressible à bas nombre de Mach.
- Configurer correctement les lois d'état (`"eos": "LINEAR"`, `"chi_T"`), la géométrie circulaire explicite (`axis_a`, `axis_b`, `pert_amp: 0.0`) et les injections pariétales.
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Quantifier dynamiquement $R(t)$, $\overline{\rho}_{\mathrm{air}}(t)$, $\overline{p}_{\mathrm{air}}(t)$ et la circularité de la bulle depuis les fichiers `.vtr`.
- Visualiser le rétrécissement isotrope de la bulle sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `EOS`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_eos`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide la compression centripète d'une bulle d'air circulaire immergée dans de l'eau au moyen du modèle 1-fluide faiblement compressible (`"compressible": true`).

Dans une cavité carrée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ sans gravité ($\mathbf{g}=\mathbf{0}$) ni tension superficielle ($\sigma = 0$), une bulle d'air de rayon initial $R_0 = 0.03\,\mathrm{m}$ ($C=0$, $\rho_{2,0} = 1.176829268\,\mathrm{kg/m}^3$, $\mu_2 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$, $\chi_{T,2} = 9.8692327\times 10^{-6}\,\mathrm{Pa}^{-1}$) est centrée en $(0.05, 0.05)\,\mathrm{m}$ dans de l'eau ($C=1$, $\rho_{1,0} = 1000\,\mathrm{kg/m}^3$, $\mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$, $\chi_{T,1} = 4.44\times 10^{-10}\,\mathrm{Pa}^{-1}$), à la pression initiale $p_0 = 101\,325\,\mathrm{Pa}$. Sur chacune des 4 parois, de l'eau est injectée normalement vers l'intérieur à la vitesse $V_0 = 0.0025\,\mathrm{m/s}$, soit un débit surfacique entrant total :
$$q = 4 L V_0 = 4 \times 0.1 \times 0.0025 = 1.0\times 10^{-3}\,\mathrm{m}^2/\mathrm{s}.$$

### 2.2 Loi Analytique d'Évolution du Rayon, de la Densité et de la Pression
L'eau étant quasi-incompressible ($\chi_{T,1} \ll \chi_{T,2}$), le débit entrant $q$ comprime quasi intégralement le volume gazeux : $\pi R^2(t) = \pi R_0^2 - q\,t$, d'où un temps d'effondrement théorique $t_{\mathrm{collapse}} = \pi R_0^2 / q \approx 2.827433\,\mathrm{s}$ et les lois exactes pour $t < t_{\mathrm{collapse}}$ :
$$R_{\mathrm{th}}(t) = R_0\sqrt{1 - \frac{q\,t}{\pi R_0^2}}, \qquad \rho_{\mathrm{th}}(t) = \frac{\rho_{2,0}}{1 - \dfrac{q\,t}{\pi R_0^2}}, \qquad p_{\mathrm{th}}(t) = \frac{p_0}{1 - \dfrac{q\,t}{\pi R_0^2}}.$$
À $t = 1.5\,\mathrm{s}$, la solution analytique donne $R_{\mathrm{th}}(1.5\,\mathrm{s}) \approx 0.020554\,\mathrm{m}$, $\rho_{\mathrm{th}}(1.5\,\mathrm{s}) \approx 2.506647\,\mathrm{kg/m}^3$ et $p_{\mathrm{th}}(1.5\,\mathrm{s}) \approx 215\,821.8\,\mathrm{Pa}$.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et bulle initiale | $[0, L]^2$, $(x_c, y_c)$, $R_0$ | $L = 0.1\,\mathrm{m}$, centre $(0.05, 0.05)\,\mathrm{m}$, $R_0 = 0.03\,\mathrm{m}$ (`axis_a = axis_b = 0.03`, `pert_amp = 0.0`), $p_0 = 101\,325\,\mathrm{Pa}$ |
| Fluides 1 (Eau) et 2 (Air) | $\rho_{i,0}, \mu_i, \chi_{T,i}$, `eos` | Eau : $\rho_{1,0}=1000$, $\chi_{T,1}=4.44\times 10^{-10}$ ; Air : $\rho_{2,0}=1.176829268$, $\chi_{T,2}=9.8692327\times 10^{-6}\,\mathrm{Pa}^{-1}$ (`"LINEAR"`) |
| Conditions d'injection (4 bords) | Gauche, Droite, Bas, Haut | `DIRICHLET` entrant : $u_{\mathrm{left}}=+V_0$, $u_{\mathrm{right}}=-V_0$, $v_{\mathrm{bottom}}=+V_0$, $v_{\mathrm{top}}=-V_0$ ($V_0=0.0025\,\mathrm{m/s}$) |
| Séquence de grilles et $\Delta t$ | $N \in \{32, 64, 128, 256, 512\}$ | $\Delta t = 0.00125\,\mathrm{s}$, $t_{\mathrm{final}} = 1.5\,\mathrm{s}$ (1200 pas ; extension possible jusqu'à $t = 2.5\,\mathrm{s}$) |
| Solveurs et schémas VOF | `solver` & `vof` | `MONOLITHIC` ($r_{\mathrm{al}}=0$) ou `C-AAL` avec `BICGSTAB2` / `GMRES` (`BLOCK_TRIANGULAR`, `ILU0`, `SSOR`) ; `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC`, `VOF_TVD` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_11_compression_bulle_n128.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/n128"`.)*

```json
{
  "case_name": "cas_11_compression_bulle_n128",
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 0.1,
    "ly": 0.1
  },
  "physics": {
    "enable_two_phase": true,
    "compressible": true,
    "flow_regime": "COMPRESSIBLE",
    "enable_convection": true,
    "enable_surface_tension": false,
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": { "rho": 1000.0,      "mu": 1.0e-3, "eos": "LINEAR", "chi_T": 4.44e-10 },
    "fluid2": { "rho": 1.176829268, "mu": 1.85e-5, "eos": "LINEAR", "chi_T": 9.8692327e-6 }
  },
  "initial_condition": {
    "type": "BUBBLE",
    "center_x": 0.05,
    "center_y": 0.05,
    "radius": 0.03,
    "axis_a": 0.03,
    "axis_b": 0.03,
    "pert_amp": 0.0,
    "p_init": 101325.0
  },
  "time": {
    "dt": 0.00125,
    "t_final": 1.5,
    "time_scheme": "EULER",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u":  0.0025, "v":  0.0 },
    "right":  { "type": "DIRICHLET", "u": -0.0025, "v":  0.0 },
    "bottom": { "type": "DIRICHLET", "u":  0.0,    "v":  0.0025 },
    "top":    { "type": "DIRICHLET", "u":  0.0,    "v": -0.0025 }
  },
  "solver": {
    "type": "MONOLITHIC",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-10,
    "max_iter": 1000,
    "r_al": 0.0,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.5
  },
  "output": {
    "output_dir": "output/cpp/n128",
    "case_name": "cas_11_compression_bulle_n128",
    "vtk_interval": 100,
    "log_interval": 20,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur** (`mpirun -np 4`).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/n128 output/fortran/n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_11_compression_bulle_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_11_compression_bulle_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Convergence vers $\rho_{\mathrm{th}}(t)$, $p_{\mathrm{th}}(t)$ et $R_{\mathrm{th}}(t)$ :**
   À partir des fichiers `.pvtr`/`.vtr`, évaluer le rayon équivalent $R_{\mathrm{num}}(t) = \sqrt{\frac{1}{\pi}\int_\Omega (1-C)\,\mathrm{d}\Omega}$, la densité moyenne $\overline{\rho}_{\mathrm{air}}(t)$ et la pression moyenne $\overline{p}_{\mathrm{air}}(t)$ dans $\{C < 0.01\}$. Vérifier la décroissance des erreurs relatives lors du raffinement du maillage ($32^2 \to 512^2$).
2. **Conservation de la masse de gaz et isotropie :**
   Vérifier que la masse totale d'air $M_{\mathrm{air}}(t) = \int_\Omega (1-C)\rho\,\mathrm{d}\Omega$ reste conservée ($|M_{\mathrm{air}}(t_{\mathrm{final}}) - M_{\mathrm{air}}(0)|/M_{\mathrm{air}}(0) < 5\times 10^{-3}$) et que la bulle demeure symétrique par rapport aux axes $x = 0.05\,\mathrm{m}$ et $y = 0.05\,\mathrm{m}$.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier la concordance de $R_{\mathrm{num}}(t)$, $\overline{\rho}_{\mathrm{air}}(t)$ et $\overline{p}_{\mathrm{air}}(t)$ entre les deux moteurs.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` (`9` colonnes en C++ dont `bubble_volume` $= \pi R_{\mathrm{num}}^2(t)$ ; `11` colonnes en Fortran incluant `avg_rho_gas,avg_p_gas`).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 21 (Compression isotrope de bulle)."""
import math
import xml.etree.ElementTree as ET
from pathlib import Path

def eval_bubble(pvtr_path: str):
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    vol_b, sum_rho, sum_p, sum_v_pure = 0.0, 0.0, 0.0, 0.0
    for piece in root.findall(".//Piece"):
        vtr = ET.parse(pvtr.parent / piece.attrib["Source"]).getroot()
        ext = [int(v) for v in piece.attrib["Extent"].split()]
        nx_l, ny_l = ext[1] - ext[0], ext[3] - ext[2]
        xc = [float(v) for v in vtr.find(".//DataArray[@Name='XCoordinates']").text.split()]
        yc = [float(v) for v in vtr.find(".//DataArray[@Name='YCoordinates']").text.split()]
        C = [float(v) for v in vtr.find(".//DataArray[@Name='PhaseIndicator']").text.split()]
        P = [float(v) for v in vtr.find(".//DataArray[@Name='Pressure']").text.split()]
        R = [float(v) for v in vtr.find(".//DataArray[@Name='Density']").text.split()]
        for j in range(ny_l):
            dy = yc[j + 1] - yc[j]
            for i in range(nx_l):
                dx = xc[i + 1] - xc[i]
                dvol = dx * dy
                k = j * nx_l + i
                vol_b += (1.0 - C[k]) * dvol
                if C[k] < 0.01:
                    sum_rho += R[k] * dvol
                    sum_p += P[k] * dvol
                    sum_v_pure += dvol
    r_num = math.sqrt(vol_b / math.pi)
    return r_num, sum_rho / sum_v_pure, sum_p / sum_v_pure

if __name__ == "__main__":
    r0, rho0, p0, q, dt = 0.03, 1.176829268, 101325.0, 1.0e-3, 0.00125
    for pvtr in sorted(Path("output/cpp/n128").glob("*.pvtr")):
        t = int(pvtr.stem.split("_")[-1]) * dt
        rn, rhon, pn = eval_bubble(str(pvtr))
        f = 1.0 - q * t / (math.pi * r0**2)
        print(f"t={t:.3f} s | R={rn:.6f} (th={r0*math.sqrt(f):.6f}) "
              f"| rho={rhon:.5f} (th={rho0/f:.5f}) | p={pn:.1f} (th={p0/f:.1f})")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Appliquer le filtre **Contour** sur `PhaseIndicator = 0.5` et **Glyph** sur `Velocity` pour visualiser l'écoulement radial convergent qui comprime la bulle de manière isotrope.
3. Afficher `Pressure` et vérifier l'uniformité spatiale du champ de pression à chaque instant.
# Tutoriel 22 : Campagnes Hybrides FC-WY-CVI-Marker (`RUN_FC-WY-CVI-Marker`)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/RUN_FC-WY-CVI-Marker/` et `Guides_Tutos/Run_FC-WY-CVI-Marker/` *(8 études ciblées + `verify_cvi_markers.py`)*.

Ce tutoriel transverse présente l'utilisation des marqueurs lagrangiens non connectés couplés à l'interpolation conservative CVI (Variante B exacte) et au suivi VOF-PLIC de Weymouth-Yue (`use_markers: true`), ainsi que le protocole de non-régression automatisé `verify_cvi_markers.py`.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

- Activer et paramétrer le nuage de marqueurs lagrangiens (`vof.use_markers`, `markers_per_cell`, `marker_advection_scheme: "CVI_RK2"`, `marker_curvature`, `marker_normals_plic`, `marker_max_per_cell`, `write_markers_vtk`).
- Valider mathématiquement la divergence discrète nulle de CVI Variante B face aux schémas bilinéaires classiques via `verify_cvi_markers.py`.
- Comparer les simulations `VOF seul` (`use_markers: false`) aux simulations hybrides `FC-WY-CVI-Marker` (`use_markers: true`) sur l'ensemble des 8 études dédiées.
- Exploiter simultanément les champs VTK eulériens (`.pvtr`/`.vtr`) et les nuages lagrangiens (`_markers_<step:06d>.vtk`) sous ParaView.

---

## 2. Inventaire des Huit Campagnes CVI-Marqueurs

L'ensemble des études testant les marqueurs CVI est regroupé dans `VALIDATION_TOTALE/Run_CASES/RUN_FC-WY-CVI-Marker/` (avec miroir documentaire dans `Guides_Tutos/Run_FC-WY-CVI-Marker/`) :

| N° | Fichier d'Instructions | Phénomène Physique Étudié | Gain Ciblé par `FC-WY-CVI-Marker` |
| :-: | :--- | :--- | :--- |
| **01** | `01_cas_goutte_statique_young_laplace_cvi_markers.md` | Équilibre capillaire de Young-Laplace | Courbure parabolique $\kappa^{\text{marker}}$ et réduction des courants parasites |
| **02** | `02_cas_tourbillon_cisaille_rider_kothe_cvi_markers.md` | Tourbillon cisaillé réversible de Rider-Kothe | Advection CVI sans divergence et préservation de la queue spiralée |
| **03** | `03_cas_chute_cylindre_dense_cvi_markers.md` | Chute libre d'un cylindre dense ($\sigma=0$) | Préservation de la circularité sous-maille et conservation d'aire |
| **04** | `04_cas_oscillation_bulle_lamb_cvi_markers.md` | Oscillations capillaires amorties de Lamb | Fidélité de la fréquence propre $f_2$ et de l'amortissement visqueux |
| **05** | `05_cas_instabilite_rayleigh_taylor_cvi_markers.md` | Instabilité gravitaire de Rayleigh-Taylor | Résolution fine des volutes secondaires de Kelvin-Helmholtz |
| **06** | `06_cas_nappe_liquide_kelvin_helmholtz_cvi_markers.md` | Cisaillement d'une nappe liquide mince | Suivi continu des deux interfaces minces sans rupture artificielle |
| **07** | `07_cas_ascension_bulle_2d_cvi_markers.md` | Ascension de bulle 2D (benchmarks de Hysing) | Intégrité des jupes latérales à grand nombre d'Eötvös/Bond |
| **08** | `08_cas_flambement_jet_visqueux_cvi_markers.md` | Flambement gravitationnel d'un jet visqueux (`JET_INLET`) | Capture des lobes repliés et ré-ensemencement sous injection continue |

---

## 3. Protocole de Vérification Automatique (`verify_cvi_markers.py`)

Avant de lancer l'une des 8 études complètes, exécuter le banc de test à 3 paliers depuis la racine du projet :

```bash
python3 Guides_Tutos/Run_FC-WY-CVI-Marker/verify_cvi_markers.py
```

Ce script vérifie automatiquement :
1. **Palier 1 :** Annulation à la précision machine de la divergence analytique de l'interpolation CVI Variante B sur 5 000 blocs MAC aléatoires à divergence nulle.
2. **Palier 2 :** Réduction des courants parasites et précision de $\Delta P = \sigma / R_0$ sur la goutte statique de Young-Laplace ($32\times 32$ et $64\times 64$).
3. **Palier 3 :** Invariance stricte MPI `np=1` vs `mpirun -np 4` (`migrate_markers_mpi`, bande fantôme $2.5 h$, routage RK2 4 passes) et parité C++17 / Fortran 2018.

---

## 4. Paramétrage JSON Type du Bloc `vof` / `markers`

```json
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.5,
    "use_markers": true,
    "markers_per_cell": 6,
    "marker_advection_scheme": "CVI_RK2",
    "marker_curvature": true,
    "marker_normals_plic": true,
    "marker_max_per_cell": 12,
    "write_markers_vtk": true
  }
```

### Visualisation Duale Euler-Lagrange sous ParaView
1. Ouvrir simultanément la série eulérienne `<case_name>_..pvtr` (isocontour `PhaseIndicator = 0.5`) et la série lagrangienne `<case_name>_markers_..vtk`.
2. Représenter les marqueurs par des glyphes sphériques colorés par l'attribut `Curvature` pour inspecter la régularité géométrique le long de l'interface.
