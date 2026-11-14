# Tutoriel 2 : Opérateur Laplacien avec Conditions de Neumann Non Homogènes et Jauge de Pression

Ce tutoriel aborde la résolution d'une équation de Poisson sous conditions aux limites de Neumann pur. Il met en lumière deux aspects cruciaux du calcul CFD : le traitement rigoureux des flux pariétaux sur les faces décalées MAC et la levée du mode singulier à constante (jauge de pression) pour garantir l'unicité de la solution et la convergence du solveur.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements physiques et les équations aux dérivées partielles gouvernant l'écoulement.
- Configurer les paramètres numériques et conditions aux limites sur grille cartésienne décalée MAC.
- Exécuter le cas test avec les exécutables C++17 et Fortran 2018 (séquentiel et OpenMP).
- Quantifier avec rigueur la précision par rapport à la solution analytique exacte ($L_2$, $L_\infty$).
- Produire et exploiter les visualisations scientifiques volumiques sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 :** `Mesh2D`, `FlowFields`, `Config`, `PoissonSolver`, `PressureGauge`, `PCGSolver`, `ILU0Preconditioner`.
- **Moteur Fortran 2018 :** `module_mesh`, `module_flow_fields`, `module_config`, `module_poisson`, `module_pressure_gauge`, `module_linear_solvers`, `module_preconditioners`.

---

## 2. Formulation Physique et Solution Analytique Exacte

Équations Physiques et Configuration

Ce cas vérifie le traitement des conditions aux limites de Neumann non homogènes (flux normaux imposés sur $\partial\Omega$) et la régularisation de jauge d'une matrice Laplacienne singulière dont le noyau est engendré par les fonctions constantes.

Sur le domaine carré unitaire $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$, on résout l'équation de Laplace homogène ($S = 0$) :
$$\Delta \Psi = \frac{\partial^2 \Psi}{\partial x^2} + \frac{\partial^2 \Psi}{\partial y^2} = 0 \quad \text{dans } \Omega,$$
avec les conditions aux limites de Neumann suivantes sur les quatre frontières :
$$\left.\frac{\partial \Psi}{\partial n}\right|_{x=0} = -\left.\frac{\partial \Psi}{\partial x}\right|_{x=0} = 0, \qquad \left.\frac{\partial \Psi}{\partial n}\right|_{x=1} = +\left.\frac{\partial \Psi}{\partial x}\right|_{x=1} = +2,$$
$$\left.\frac{\partial \Psi}{\partial n}\right|_{y=0} = -\left.\frac{\partial \Psi}{\partial y}\right|_{y=0} = 0, \qquad \left.\frac{\partial \Psi}{\partial n}\right|_{y=1} = +\left.\frac{\partial \Psi}{\partial y}\right|_{y=1} = -2.$$
La condition de compatibilité $\oint_{\partial\Omega} \frac{\partial\Psi}{\partial n}\,\mathrm{d}\gamma = 0 + 2(1) + 0 - 2(1) = 0$ est exactement vérifiée. En fixant la constante de jauge par $\Psi(0,0) = 0$, la solution analytique exacte est le polynôme harmonique quadratique :
$$\Psi_{\mathrm{ex}}(x,y) = x^2 - y^2.$$

### Propriété théorique fondamentale (exactitude machine)
Le développement de Taylor du schéma aux volumes finis centré à 5 points montre que l'erreur de troncature locale s'écrit :
$$\tau_{I,J} = \frac{\Delta x^2}{12}\frac{\partial^4 \Psi_{\mathrm{ex}}}{\partial x^4} + \frac{\Delta y^2}{12}\frac{\partial^4 \Psi_{\mathrm{ex}}}{\partial y^4} + \mathcal{O}(\Delta x^4) \equiv 0$$
car toutes les dérivées d'ordre $\ge 3$ de $x^2 - y^2$ sont identiquement nulles. Après recalage de la moyenne de jauge :
$$\widetilde{\Psi}_{I,J} = \Psi_{I,J} - \frac{1}{N_x N_y}\sum_{p,q}\left(\Psi_{p,q} - \Psi_{\mathrm{ex}}(x_p, y_q)\right),$$
l'erreur discrète se situe au niveau du zéro machine ($\sim 10^{-15} - 10^{-12}$) quelle que soit la finesse du maillage.

---

---

## 3. Paramètres Physiques, Numériques et Discrétisation

Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Géométrie | Domaine $[0, L_x] \times [0, L_y]$ | $L_x = 1.0\,\mathrm{m}$, $L_y = 1.0\,\mathrm{m}$ |
| Séquence de grilles | $(N_x \times N_y)$ | $16\times 16$, $32\times 32$, $64\times 64$, $128\times 128$, $256\times 256$ |
| Terme source | $S$ (`source_term`) | $0.0\,\mathrm{m}^{-2}$ |
| Flux de Neumann | `dpsi_dn` (Gauche, Droite, Bas, Haut) | $0.0$, $+2.0$, $0.0$, $-2.0$ |
| Solveurs linéaires à tester | `solver.linear_solver` | `DIRECT`, `CG`, `BICGSTAB2`, `MINRES`, `GMRES` ($\mathrm{tol}=10^{-14}$, $\mathrm{max\_iter}=10000$) |
| Préconditionneurs à tester | `solver.precond_type` | `NONE`, `JACOBI`, `SSOR`, `ILU0`, `BLOCK_TRIANGULAR` |

---

---

## 4. Fichier de Configuration JSON Pas à Pas

de Fichiers `.json` (C++17 et Fortran 2018)

### Prototype `case/cpp/cas_2_neumann_grid64_direct.json`
```json
{
  "case_name": "cas_2_neumann_grid64_direct",
  "mode": "LAPLACIAN",
  "test_case": "CAS_2_NEUMANN",
  "source_term": 0.0,
  "mesh": {
    "nx": 64,
    "ny": 64,
    "lx": 1.0,
    "ly": 1.0
  },
  "boundary_conditions": {
    "left":   { "type": "NEUMANN", "dpsi_dn":  0.0 },
    "right":  { "type": "NEUMANN", "dpsi_dn":  2.0 },
    "bottom": { "type": "NEUMANN", "dpsi_dn":  0.0 },
    "top":    { "type": "NEUMANN", "dpsi_dn": -2.0 }
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
    "output_dir": "VALIDATION_TOTALE/validation_cas_2_neumann/output/cpp/direct/grid_64x64",
    "case_name": "cas_2_neumann_grid64_direct",
    "vtk_interval": 1,
    "log_interval": 1,
    "write_csv_history": true
  }
}

```

---

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

L'exécution s'effectue indifféremment à l'aide du binaire C++17 (`covof2d_cpp`) ou du binaire Fortran 2018 (`covof2d_fortran`), garantissant une stricte équivalence des résultats numériques.

```bash
# 1. Se placer dans le répertoire du cas
cd VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_2_neumann

# 2. Exécution avec le moteur C++17 (Séquentiel)
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_2_config.json

# 3. Exécution avec le moteur C++17 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_2_config.json

# 4. Exécution avec le moteur Fortran 2018 (Séquentiel)
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_2_config.json

# 5. Exécution avec le moteur Fortran 2018 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_2_config.json

```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification et Tolérances
à Mener (Précision, Intégrité et Performance)

1. **Exactitude machine sur polynôme quadratique :**
   - Vérifier qu'avec le solveur `DIRECT`, les erreurs $E_{L_1}, E_{L_2}, E_{L_\infty}$ sont inférieures à $3\times 10^{-14}$ sur toutes les grilles de $16\times 16$ à $256\times 256$ (par ex. $E_{L_1}(16\times 16) \approx 2.34\times 10^{-16}$).
2. **Intégrité de la régularisation du noyau singulier ($\ker(\Delta_N) = \mathrm{Vect}\{\mathbf{1}\}$) :**
   - Vérifier que l'ancrage de jauge dans `DirectSolver` élimine le mode constant sans perturber les flux de Neumann aux bords ($\sum_{I,J} (\Delta_h \Psi)_{I,J}\Delta x\Delta y = 0$ à $10^{-16}$ près).
   - Vérifier que les solveurs de Krylov (`CG`, `MINRES`, `BICGSTAB2`, `GMRES`) couplés aux préconditionneurs `BLOCK_TRIANGULAR` (qui projette explicitement le résidu sur l'hyperplan à moyenne nulle $\sum_c r_p[c] = 0$), `JACOBI`, `SSOR` et `ILU0` atteignent une erreur $E_{L_2} \le 10^{-12}$ sans dérive dans le noyau.
   - Vérifier l'antisymétrie diagonale $\Psi(y,x) = -\Psi(x,y)$ à $10^{-15}$ près.
3. **Performance et parité C++17 / Fortran 2018 :**
   - Comparer le nombre d'itérations et le temps CPU de chaque préconditionneur, et vérifier $\|\widetilde{\Psi}^{\mathrm{C++}} - \widetilde{\Psi}^{\mathrm{Fortran}}\|_\infty < 10^{-14}$.

---

### 6.2 Résultats Attendus et Visualisation ParaView

