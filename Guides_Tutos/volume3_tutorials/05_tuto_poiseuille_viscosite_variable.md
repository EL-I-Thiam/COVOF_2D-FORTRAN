# Tutoriel 5 : Écoulement de Poiseuille Plan à Viscosité Variable Exponentielle

Ce tutoriel explore l'effet d'une hétérogénéité spatiale continue de la viscosité dynamique mu(y) = mu_0 exp(beta y) sur un écoulement de Poiseuille. Il valide la moyenne harmonique de la viscosité aux faces des volumes de contrôle et vérifie le profil asymétrique de vitesse ainsi que le déplacement du pic de vitesse vers la zone de moindre viscosité.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements physiques et les équations aux dérivées partielles gouvernant l'écoulement.
- Configurer les paramètres numériques et conditions aux limites sur grille cartésienne décalée MAC.
- Exécuter le cas test avec les exécutables C++17 et Fortran 2018 (séquentiel et OpenMP).
- Quantifier avec rigueur la précision par rapport à la solution analytique exacte ($L_2$, $L_\infty$).
- Produire et exploiter les visualisations scientifiques volumiques sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 :** `Mesh2D`, `FlowFields`, `VariableViscosity`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `SSORPreconditioner`.
- **Moteur Fortran 2018 :** `module_mesh`, `module_flow_fields`, `module_viscosity`, `module_navier_stokes`, `module_linear_solvers`, `module_preconditioners`.

---

## 2. Formulation Physique et Solution Analytique Exacte

Équations Physiques et Configuration

Ce cas vérifie la discrétisation conservative du terme de diffusion visqueuse $\nabla \cdot [\mu(\mathbf{x})(\nabla \mathbf{u} + \nabla \mathbf{u}^T)]$ lorsque la viscosité dynamique varie fortement dans l'espace.

Dans un canal plan $\Omega = [0, L_x] \times [0, L_y]$ avec $L_x = 0.05\,\mathrm{m}$ et $L_y = 2e = 0.02\,\mathrm{m}$ (demi-largeur $e = 0.01\,\mathrm{m}$, coordonnée transversale centrée $z = y - e \in [-e, +e]$), la viscosité dynamique suit une loi exponentielle symétrique :
$$\mu(z) = \mu_0 \exp\left(a \frac{z^2}{e^2}\right), \qquad \text{avec } \mu_0 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}, \quad a = 2.0, \quad e = 0.01\,\mathrm{m}.$$
Le rapport de viscosité entre les parois ($|z|=e$) et l'axe central ($z=0$) atteint $\mu(\pm e)/\mu_0 = e^2 \approx 7.389$. Sous un forçage axial constant $G_x = -\partial p/\partial x = \rho g_x = 0.294207317\,\mathrm{Pa/m}$, l'équation différentielle transversale s'écrit :
$$\frac{\mathrm{d}}{\mathrm{d}z}\left(\mu(z)\frac{\mathrm{d}u}{\mathrm{d}z}\right) = -G_x, \qquad u(-e) = u(+e) = 0.$$

### Solution analytique exacte
Une première intégration avec $\frac{\mathrm{d}u}{\mathrm{d}z}(0) = 0$ donne la contrainte de cisaillement linéaire $\tau_{xz}(z) = \mu(z)\frac{\mathrm{d}u}{\mathrm{d}z} = -G_x z$. En remarquant que $\frac{\mathrm{d}}{\mathrm{d}z}\left[\exp(-a z^2/e^2)\right] = -\frac{2az}{e^2}\exp(-a z^2/e^2)$, une seconde intégration fournit le profil analytique gaussien tronqué :
$$u_{\mathrm{ex}}(z) = \frac{G_x e^2}{2 a \mu_0}\left[\exp\left(-a\frac{z^2}{e^2}\right) - \exp(-a)\right],$$
de vitesse maximale au centre $u_{\max} = u_{\mathrm{ex}}(0) = \frac{G_x e^2}{2 a \mu_0}(1 - e^{-2}) \approx 0.343771\,\mathrm{m/s}$.
Sur le maillage décalé MAC, les flux de cisaillement aux faces Nord et Sud des volumes de contrôle de $u_{i,j}$ sont évalués avec la viscosité aux nœuds de vorticité $y_j = j\Delta y$ :
$$\left.\frac{\partial}{\partial y}\left(\mu \frac{\partial u}{\partial y}\right)\right|_{i,j} \approx \frac{1}{\Delta y}\left[\mu_{\mathrm{N}}\frac{u_{i,j+1}-u_{i,j}}{\Delta y} - \mu_{\mathrm{S}}\frac{u_{i,j}-u_{i,j-1}}{\Delta y}\right],$$
ce qui assure une précision d'ordre 2 ($\mathcal{O}(\Delta y^2)$).

---

---

## 3. Paramètres Physiques, Numériques et Discrétisation

Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------------- | :----------------------- | :------------------------------------------------------------------ |
| Domaine | $[0, L_x] \times [0, L_y]$ | $L_x = 0.05\,\mathrm{m}$, $L_y = 0.02\,\mathrm{m}$ ($e = 0.01\,\mathrm{m}$) |
| Loi de viscosité | `viscosity_model: "EXPONENTIAL"` | $\mu_0 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$, $a = 2.0$, $e = 0.01\,\mathrm{m}$, $\rho = 1.0\,\mathrm{kg/m}^3$ |
| Forçage axial | $g_x, g_y$ | $g_x = 0.294207317\,\mathrm{m/s}^2$, $g_y = 0.0$ |
| Séquence de grilles | $N_x \times N_y$ | $N_x = 16$, $N_y \in \{16, 32, 64, 128, 256, 512\}$ |
| Intégration temporelle | $\Delta t$, $t_{\mathrm{final}}$ | $\Delta t = 0.5\,\mathrm{s}$, $t_{\mathrm{final}} = 50.0\,\mathrm{s}$, $\mathrm{steady\_state\_tol} = 10^{-12}$ |
| Solveurs testés | Formulations | `AL`, `AAL`, `MONO` avec `DIRECT`, `BICGSTAB2`, `BLOCK_TRIANGULAR`, `SSOR`, `ILU0` |

---

---

## 4. Fichier de Configuration JSON Pas à Pas

de Fichiers `.json` (C++17 et Fortran 2018)

### Prototype `case/cpp/cas_5_var_visc_al_direct_ny128.json`
```json
{
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
    "linear_solver": "DIRECT",
    "tol": 1e-12,
    "max_iter": 5000,
    "r_al": 0.0,
    "uzawa_iter": 1,
    "use_precond": false,
    "precond_type": "NONE",
    "steady_state_tol": 1e-12
  },
  "output": {
    "output_dir": "VALIDATION_TOTALE/validation_cas_5_poiseuille_visc_variable/output/cpp/al_direct/ny128",
    "case_name": "cas_5_var_visc_al_direct_ny128",
    "vtk_interval": 200,
    "log_interval": 50,
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
cd VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_5_poiseuille_visc_variable

# 2. Exécution avec le moteur C++17 (Séquentiel)
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_5_config.json

# 3. Exécution avec le moteur C++17 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_5_config.json

# 4. Exécution avec le moteur Fortran 2018 (Séquentiel)
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_5_config.json

# 5. Exécution avec le moteur Fortran 2018 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_5_config.json

```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification et Tolérances
à Mener (Précision, Intégrité et Performance)

1. **Convergence spatiale d'ordre 2 et linéarité du flux de cisaillement :**
   - Vérifier que $E_{L_1}$ décroît de $6.365\times 10^{-4}\,\mathrm{m/s}$ ($N_y=16$) à $5.487\times 10^{-7}\,\mathrm{m/s}$ ($N_y=512$), soit un ordre de Richardson $p_{L_1}, p_{L_2}, p_{L_\infty} \in [2.00, 2.07]$.
   - Vérifier que la contrainte de cisaillement reconstruite aux faces $\tau_{xz}(y_j) = \mu(y_j)\frac{u_{i,j+1}-u_{i,j}}{\Delta y}$ coïncide avec la droite exacte $-G_x (y_j - e)$ à $10^{-6}\,\mathrm{Pa}$ près.
2. **Sensibilité paramétrique à l'exposant d'hétérogénéité $a$ (`visc_a`) :**
   - Tester `"visc_a"` $\in \{0.0, 1.0, 2.0, 4.0\}$ (correspondant à des rapports $\mu_{\max}/\mu_{\min} \in \{1, 2.72, 7.39, 54.6\}$) et vérifier que pour $a \to 0$ on retrouve la solution parabolique exacte à la précision machine, tandis que l'ordre 2 est maintenu pour tout $a > 0$.
3. **Performance des préconditionneurs à coefficients variables et parité C++/Fortran :**
   - Vérifier que la prise en compte des viscosités variables $\mu_e, \mu_w, \mu_n, \mu_s$ dans les blocs $F_u, F_v$ de `BLOCK_TRIANGULAR` et dans la diagonale de `JACOBI`/`ILU0` maintient une convergence rapide de `BICGSTAB2`.
   - Vérifier la parité machine $|E_{L_1}^{\mathrm{C++}} - E_{L_1}^{\mathrm{Fortran}}| < 10^{-15}$.

---

### 6.2 Résultats Attendus et Visualisation ParaView

