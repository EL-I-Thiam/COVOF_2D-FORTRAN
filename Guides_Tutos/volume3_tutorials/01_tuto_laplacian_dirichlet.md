# Tutoriel 1 : Opérateur Laplacien avec Conditions de Dirichlet Homogènes

Ce premier tutoriel constitue le point de départ fondamental de la vérification de COVOF 2D. Il permet de valider la discrétisation spatiale à l'ordre 2 de l'opérateur Laplacien scalaire 2D (équation de Poisson) avec conditions aux limites de Dirichlet homogènes sur maillage cartésien décalé MAC.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements physiques et les équations aux dérivées partielles gouvernant l'écoulement.
- Configurer les paramètres numériques et conditions aux limites sur grille cartésienne décalée MAC.
- Exécuter le cas test avec les exécutables C++17 et Fortran 2018 (séquentiel et OpenMP).
- Quantifier avec rigueur la précision par rapport à la solution analytique exacte ($L_2$, $L_\infty$).
- Produire et exploiter les visualisations scientifiques volumiques sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 :** `Mesh2D`, `FlowFields`, `Config`, `PoissonSolver`, `BiCGSTAB2Solver`, `SSORPreconditioner`, `VTKWriter`.
- **Moteur Fortran 2018 :** `module_mesh`, `module_flow_fields`, `module_config`, `module_poisson`, `module_linear_solvers`, `module_preconditioners`, `module_vtk_writer`.

---

## 2. Formulation Physique et Solution Analytique Exacte

Équations Physiques et Configuration

Ce premier cas vérifie la consistance, la précision et l'ordre de convergence spatial de la discrétisation aux volumes finis de l'opérateur Laplacien bidimensionnel stationnaire (équation de Poisson) avec conditions aux limites de Dirichlet homogènes.

Sur le domaine rectangulaire $\Omega = [0, a] \times [0, b]$ de dimensions $a = L_x = 1.0\,\mathrm{m}$ et $b = L_y = 0.5\,\mathrm{m}$, on résout :
$$\Delta \Psi = \frac{\partial^2 \Psi}{\partial x^2} + \frac{\partial^2 \Psi}{\partial y^2} = S \quad \text{dans } \Omega = ]0, a[ \times ]0, b[, \qquad \text{avec } S = 1.0\,\mathrm{m}^{-2},$$
soumis à la condition de Dirichlet homogène sur l'ensemble du bord $\partial\Omega$ :
$$\Psi(0, y) = \Psi(a, y) = \Psi(x, 0) = \Psi(x, b) = 0, \qquad \forall (x,y) \in \partial\Omega.$$

### Solution analytique exacte
Par décomposition sur la base propre orthogonale de l'opérateur Laplacien avec conditions de Dirichlet sur $[0,a]\times[0,b]$, seuls les modes impairs $(2i+1, 2j+1)$ contribuent au développement en série de Fourier double :
$$\Psi_{\mathrm{ex}}(x,y) = -\frac{16 S}{\pi^4} \sum_{i=0}^{\infty} \sum_{j=0}^{\infty} \frac{\sin\left(\frac{(2i+1)\pi x}{a}\right)\sin\left(\frac{(2j+1)\pi y}{b}\right)}{(2i+1)(2j+1)\left[\left(\frac{2i+1}{a}\right)^2 + \left(\frac{2j+1}{b}\right)^2\right]}.$$
Dans les codes C++17 et Fortran 2018, cette série est évaluée aux centres des mailles $(x_I, y_J) = ((I-0.5)\Delta x, (J-0.5)\Delta y)$ avec $60 \times 60$ modes impairs ($i, j \in \{0, \dots, 59\}$). La valeur minimale au centre géométrique $(x,y) = (0.5, 0.25)\,\mathrm{m}$ vaut $\Psi_{\min} \approx -0.0366$.

Sur un maillage cartésien uniforme à mailles carrées $\Delta x = a/N_x = \Delta y = b/N_y$, le Laplacien est approché par le schéma centré à 5 points d'ordre 2 :
$$\frac{\Psi_{I+1,J} - 2\Psi_{I,J} + \Psi_{I-1,J}}{\Delta x^2} + \frac{\Psi_{I,J+1} - 2\Psi_{I,J} + \Psi_{I,J-1}}{\Delta y^2} = S,$$
et la condition $\Psi = 0$ sur les parois situées à $\Delta x/2$ du premier centre de maille intérieur est imposée par antisymétrie linéaire dans la maille fantôme ($\Psi_{0,J} = -\Psi_{1,J}$), modifiant le stencil de bord en $(\Psi_{2,J} - 3\Psi_{1,J})/\Delta x^2$.

---

---

## 3. Paramètres Physiques, Numériques et Discrétisation

Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Géométrie | Domaine $[0, L_x] \times [0, L_y]$ | $L_x = 1.0\,\mathrm{m}$, $L_y = 0.5\,\mathrm{m}$ (rapport d'aspect $2:1$) |
| Séquence de grilles | $(N_x \times N_y)$ à mailles carrées | $16\times 8$, $32\times 16$, $64\times 32$, $128\times 64$, $256\times 128$ |
| Pas d'espace | $\Delta x = \Delta y$ | $6.25\times 10^{-2}$, $3.125\times 10^{-2}$, $1.5625\times 10^{-2}$, $7.8125\times 10^{-3}$, $3.90625\times 10^{-3}\,\mathrm{m}$ |
| Terme source | $S$ (`source_term`) | $1.0\,\mathrm{m}^{-2}$ |
| Conditions aux limites | Gauche, Droite, Bas, Haut | `DIRICHLET` homogène ($\Psi = 0$) |
| Solveurs linéaires à tester | `solver.linear_solver` | `DIRECT` (Cholesky bande SPD), `CG`, `BICGSTAB2`, `BICGSTAB`, `GMRES`, `MINRES` ($\mathrm{tol}=10^{-12}$, $\mathrm{max\_iter}=5000$) |
| Préconditionneurs à tester | `solver.precond_type` | `NONE`, `JACOBI`, `SSOR`, `ILU0`, `BLOCK_TRIANGULAR` |

---

---

## 4. Fichier de Configuration JSON Pas à Pas

de Fichiers `.json` (C++17 et Fortran 2018)

### 5.1 Prototype Solveur Direct (`case/cpp/cas_1_grid64_direct.json`)
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
    "output_dir": "VALIDATION_TOTALE/validation_cas_1_laplacian/output/cpp/direct/grid_64x32",
    "case_name": "cas_1_grid64_direct",
    "vtk_interval": 1,
    "log_interval": 1,
    "write_csv_history": true
  }
}

```

### 5.2 Variante Solveur de Krylov Préconditionné (`case/cpp/cas_1_grid64_cg_ssor.json`)
Pour les variantes itératives et l'étude des préconditionneurs, remplacer le bloc `"solver"` par :
```json
  "solver": {
    "type": "MONOLITHIC",
    "linear_solver": "CG",
    "precond_type": "SSOR",
    "use_precond": true,
    "tol": 1e-12,
    "max_iter": 5000
  }

```
*(Pour Fortran 2018, générer les fichiers analogues dans `case/fortran/` avec `"output_dir": "VALIDATION_TOTALE/validation_cas_1_laplacian/output/fortran/..."`).*

---

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

L'exécution s'effectue indifféremment à l'aide du binaire C++17 (`covof2d_cpp`) ou du binaire Fortran 2018 (`covof2d_fortran`), garantissant une stricte équivalence des résultats numériques.

```bash
# 1. Se placer dans le répertoire du cas
cd VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_1_laplacian

# 2. Exécution avec le moteur C++17 (Séquentiel)
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_1_config.json

# 3. Exécution avec le moteur C++17 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_1_config.json

# 4. Exécution avec le moteur Fortran 2018 (Séquentiel)
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_1_config.json

# 5. Exécution avec le moteur Fortran 2018 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_1_config.json

```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification et Tolérances
à Mener (Précision, Intégrité et Performance)

1. **Précision et ordre de convergence spatial :**
   - Extraire les erreurs discrètes $E_{L_1}$, $E_{L_2}$ et $E_{L_\infty}$ depuis le fichier `<output_dir>/<case_name>_summary.csv` ou en comparant le champ `psi` du fichier `.vtk` à la série de Fourier $\Psi_{\mathrm{ex}}(x_I, y_J)$.
   - Vérifier que $E_{L_1}$ décroît de $9.88\times 10^{-5}$ ($16\times 8$) à $3.69\times 10^{-7}$ ($256\times 128$) et que l'ordre de convergence de Richardson vérifie :
     $$p_\alpha = \frac{\ln(E_{L_\alpha}(h_1)/E_{L_\alpha}(h_2))}{\ln(h_1/h_2)} \in [1.99, 2.03], \qquad \alpha \in \{1, 2, \infty\}.$$

2. **Intégrité et équivalence des solveurs linéaires et préconditionneurs :**
   - Vérifier que l'opérateur discret $-(\Delta_h \Psi)$ étant symétrique défini positif (SPD), les solveurs `DIRECT`, `CG`, `MINRES`, `BICGSTAB`, `BICGSTAB2` et `GMRES` convergent tous sans stagnation et fournissent exactement les mêmes normes d'erreur : $|E_{L_2}^{\mathrm{Krylov}} - E_{L_2}^{\mathrm{DIRECT}}| < 10^{-12}$.
   - Vérifier l'invariance de la solution finale vis-à-vis du préconditionneur (`NONE`, `JACOBI`, `SSOR`, `ILU0`, `BLOCK_TRIANGULAR`) : $\|\Psi_{P_1} - \Psi_{P_2}\|_\infty < 10^{-11}$.
   - Vérifier la symétrie géométrique quatre-quadrants de la solution : $|\Psi(x, y) - \Psi(a-x, y)| < 10^{-14}$ et $|\Psi(x, y) - \Psi(x, b-y)| < 10^{-14}$.

3. **Performance algorithmique et parité C++17 / Fortran 2018 :**
   - Dresser la table du nombre d'itérations $N_{\mathrm{iter}}$ et du temps CPU en fonction du conditionnement $\kappa(A) = \mathcal{O}(N_x^2)$ pour chaque couple `(linear_solver, precond_type)`.
   - Vérifier que le solveur `DIRECT` (Cholesky bande de demi-largeur $N_x$) affiche un nombre d'itérations égal à $1$.
   - Vérifier la parité stricte entre les binaires C++17 et Fortran 2018 : $\|\Psi^{\mathrm{C++}} - \Psi^{\mathrm{Fortran}}\|_\infty < 10^{-11}$.

---

### 6.2 Résultats Attendus et Visualisation ParaView

