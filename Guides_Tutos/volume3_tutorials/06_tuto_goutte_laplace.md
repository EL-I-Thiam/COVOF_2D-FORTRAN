# Tutoriel 6 : Goutte Statique en Équilibre Capillaire et Loi de Young-Laplace

Ce tutoriel constitue le banc d'essai incontournable de l'hydrostatique diphasique capillaire. Il valide le modèle de force superficielle continue CSF (Continuum Surface Force), la reconstruction géométrique de courbure kappa et démontre l'absence de courants parasites excessifs autour d'une goutte circulaire à l'équilibre.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements physiques et les équations aux dérivées partielles gouvernant l'écoulement.
- Configurer les paramètres numériques et conditions aux limites sur grille cartésienne décalée MAC.
- Exécuter le cas test avec les exécutables C++17 et Fortran 2018 (séquentiel et OpenMP).
- Quantifier avec rigueur la précision par rapport à la solution analytique exacte ($L_2$, $L_\infty$).
- Produire et exploiter les visualisations scientifiques volumiques sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 :** `Mesh2D`, `FlowFields`, `VOFTracker`, `CSFModel`, `NavierStokesSolver`, `UzawaSolver`.
- **Moteur Fortran 2018 :** `module_mesh`, `module_flow_fields`, `module_vof`, `module_csf`, `module_navier_stokes`, `module_uzawa`.

---

## 2. Formulation Physique et Solution Analytique Exacte

Équations Physiques et Configuration

Ce cas valide le couplage hydro-capillaire entre les équations de Navier-Stokes diphasiques et la force de tension superficielle *Continuum Surface Force* (CSF) discrétisée selon le schéma à équilibre exact de forces (*Balanced-Force*) sur maillage décalé MAC.

Dans un domaine carré fermé $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$ en apesanteur ($\mathbf{g} = \mathbf{0}$) avec conditions de glissement libre (`FREE_SLIP`) sur les 4 parois, une goutte circulaire de rayon $R_0 = 0.25\,\mathrm{m}$ est centrée en $(x_c, y_c) = (0.5, 0.5)\,\mathrm{m}$ au repos initial $\mathbf{u}(\mathbf{x}, 0) = \mathbf{0}$. Les deux fluides ont des propriétés identiques $\rho_1 = \rho_2 = 1.0\,\mathrm{kg/m}^3$, $\mu_1 = \mu_2 = 1.0\,\mathrm{Pa\cdot s}$, et un coefficient de tension superficielle $\sigma = 1.0\,\mathrm{N/m}$.

### Loi analytique de Young-Laplace et régularisation CSF
À l'équilibre exact, la vitesse est identiquement nulle ($\mathbf{u}_{\mathrm{ex}} = \mathbf{0}$) et le saut de pression à la traversée de l'interface circulaire de courbure $\kappa = 1/R_0 = 4.0\,\mathrm{m}^{-1}$ vaut :
$$\Delta P_{\mathrm{ex}} = P_{\mathrm{int}} - P_{\mathrm{ext}} = \frac{\sigma}{R_0} = \frac{1.0}{0.25} = 4.0000\,\mathrm{Pa}.$$
La force capillaire volumique $\mathbf{f}_{\mathrm{st}} = \sigma \kappa \nabla C$, où $\kappa = -\nabla \cdot (\nabla \widetilde{C}/\|\nabla \widetilde{C}\|)$, est évaluée aux centres des faces de vitesse exactement aux mêmes points que $-\nabla p$. Lorsque `enable_diffusion_smoothing: true`, l'indicatrice $C$ est préalablement régularisée par diffusion laplacienne afin d'éliminer les erreurs haute fréquence sur les dérivées secondes de la courbure et de supprimer les courants parasites (*spurious currents*).

---

---

## 3. Paramètres Physiques, Numériques et Discrétisation

Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et Goutte | $\Omega$, $(x_c, y_c)$, $R_0$ | $[0,1]\times[0,1]\,\mathrm{m}^2$, $(0.5, 0.5)\,\mathrm{m}$, $R_0 = 0.25\,\mathrm{m}$ |
| Propriétés physiques | $\rho_1=\rho_2, \mu_1=\mu_2, \sigma, \mathbf{g}$ | $\rho = 1.0\,\mathrm{kg/m}^3$, $\mu = 1.0\,\mathrm{Pa\cdot s}$, $\sigma = 1.0\,\mathrm{N/m}$, $\mathbf{g} = \mathbf{0}$ |
| Séquence de grilles et $\Delta t$ | $N \in \{32, 64, 128, 256\}$, $t_{\mathrm{final}}=0.05\,\mathrm{s}$ | $\Delta t \in \{0.005, 0.005, 0.0025, 0.00125\}\,\mathrm{s}$ |
| Axes d'étude | Axe A (3 VOF + lissage), Axe B (sans lissage), Axe C (solveurs) | `WY+PLIC`, `CIAM+PLIC`, `VOF_TVD` ; `AL` ($r_{\mathrm{al}}=100$), `AAL` ($K_{\mathrm{aal}}=100$), `MONOLITHIC` ($r_{\mathrm{al}}=0$), `STANDARD_PROJECTION` |
| Solveurs et Préconditionneurs | `linear_solver` & `precond_type` | `BICGSTAB2`, `BICGSTAB`, `GMRES`, `MINRES` avec `BLOCK_TRIANGULAR`, `JACOBI`, `SSOR`, `ILU0` |

---

---

## 4. Fichier de Configuration JSON Pas à Pas

de Fichiers `.json` (C++17 et Fortran 2018)

### Prototype `case/cpp/cas_6_axeA_wy_plic_n128.json`
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
    "enable_two_phase": true,
    "enable_surface_tension": true,
    "enable_convection": false,
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
    "radius": 0.25
  },
  "time": {
    "dt": 0.0025,
    "t_final": 0.05,
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
    "tol": 1e-07,
    "max_iter": 300,
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
    "output_dir": "VALIDATION_TOTALE/validation_cas_6_laplace_droplet/output/cpp/axeA_wy_plic/n128",
    "case_name": "cas_6_axeA_wy_plic_n128",
    "vtk_interval": 2,
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
cd VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_6_laplace_droplet

# 2. Exécution avec le moteur C++17 (Séquentiel)
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_6_config.json

# 3. Exécution avec le moteur C++17 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_6_config.json

# 4. Exécution avec le moteur Fortran 2018 (Séquentiel)
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_6_config.json

# 5. Exécution avec le moteur Fortran 2018 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_6_config.json

```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification et Tolérances
à Mener (Précision, Intégrité et Performance)

1. **Saut de pression de Young-Laplace et convergence d'ordre 2 :**
   - Calculer $P_{\mathrm{int}}$ sur $\{r < 0.18\,\mathrm{m}, C > 0.90\}$ et $P_{\mathrm{ext}}$ sur $\{r > 0.35\,\mathrm{m}, C < 0.05\}$.
   - Vérifier que $\Delta P = P_{\mathrm{int}} - P_{\mathrm{ext}}$ converge de $4.0405\,\mathrm{Pa}$ ($N=32$, erreur $1.01\,\%$) à $3.9985\,\mathrm{Pa}$ ($N=128$, erreur $0.04\,\%$) et $3.9966\,\mathrm{Pa}$ ($N=256$, erreur $0.09\,\%$), avec un ordre de convergence d'extrapolation $|\Delta P_h - \Delta P_{2h}|$ égal à $2.00$.
2. **Maîtrise des courants parasites et étude paramétrique de la régularisation CSF :**
   - Vérifier que $\|\mathbf{u}\|_\infty(N=256) \approx 3.34\times 10^{-3}\,\mathrm{m/s}$ ($\mathrm{Ca}_{\mathrm{parasite}} = \mu\|\mathbf{u}\|_\infty/\sigma \approx 3.34\times 10^{-3}$) avec `enable_diffusion_smoothing: true`, contre $0.190\,\mathrm{m/s}$ sans lissage (réduction d'un facteur $60\times$).
   - Étudier la sensibilité au paramètre `diffusion_smoothing_steps` $\in \{0, 1, 2, 3, 5\}$ et au filtre `enable_svof_smoothing: true` sur l'énergie cinétique parasite $E_k(t)$ et la largeur de la zone de transition de pression.
   - Vérifier que dans ce domaine fermé avec conditions `FREE_SLIP`, la projection à moyenne nulle $\sum_{c} r_p[c] = 0$ dans `BlockTriangularPreconditioner` garantit l'absence totale de dérive de la pression moyenne $\int_\Omega p\,\mathrm{d}\Omega = 0$.
3. **Invariance des solveurs et parité dual-code :**
   - Vérifier que `AL`, `AAL` et `MONOLITHIC` donnent le même $\Delta P$ à $6\times 10^{-4}\,\mathrm{Pa}$ près pour $\mathrm{tol} \in \{10^{-7}, 10^{-10}, 10^{-12}\}$, et que $|p_{\mathrm{C++}} - p_{\mathrm{Fortran}}| \le 1.2\times 10^{-14}\,\mathrm{Pa}$.

---

### 6.2 Résultats Attendus et Visualisation ParaView

