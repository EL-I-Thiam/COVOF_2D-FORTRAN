# Tutoriel 4 : Écoulement de Poiseuille Plan Instationnaire (Démarrage Impulsionnel)

Ce tutoriel aborde le régime transitoire d'établissement visqueux d'un écoulement de Poiseuille plan à partir du repos sous l'effet soudain d'un gradient de pression. La solution analytique en série de Fourier permet de valider avec une haute précision les schémas d'intégration en temps (Euler implicite, Crank-Nicolson) et la stabilité temporelle.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements physiques et les équations aux dérivées partielles gouvernant l'écoulement.
- Configurer les paramètres numériques et conditions aux limites sur grille cartésienne décalée MAC.
- Exécuter le cas test avec les exécutables C++17 et Fortran 2018 (séquentiel et OpenMP).
- Quantifier avec rigueur la précision par rapport à la solution analytique exacte ($L_2$, $L_\infty$).
- Produire et exploiter les visualisations scientifiques volumiques sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 :** `Mesh2D`, `FlowFields`, `TimeIntegrator`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `SSORPreconditioner`.
- **Moteur Fortran 2018 :** `module_mesh`, `module_flow_fields`, `module_time`, `module_navier_stokes`, `module_linear_solvers`, `module_preconditioners`.

---

## 2. Formulation Physique et Solution Analytique Exacte

Équations Physiques et Configuration

Ce cas étudie l'établissement transitoire d'un écoulement de Poiseuille plan partant du repos sous l'action brusque d'un terme source de pression axial constant. Il valide le couplage entre le schéma d'intégration temporelle (`EULER` d'ordre 1 ou `GEAR`/BDF2 d'ordre 2) et la discrétisation spatiale visqueuse d'ordre 2.

Dans un canal périodique en $x$ de dimensions $L_x = 0.16\,\mathrm{m}$ et $H = L_y = 2e = 0.16\,\mathrm{m}$ ($e = 0.08\,\mathrm{m}$), un fluide au repos initial $u(x,y,0) = 0, v(x,y,0) = 0$ de propriétés $\rho = 1.0\,\mathrm{kg/m}^3$ et $\mu = 0.01\,\mathrm{Pa\cdot s}$ ($\nu = 0.01\,\mathrm{m^2/s}$) est soumis à l'instant $t = 0^+$ à un forçage axial constant :
$$k = -\frac{\partial p}{\partial x} = \rho g_x = 0.0555\,\mathrm{N/m}^3.$$
L'équation d'évolution de la vitesse longitudinale $u(y,t)$ s'écrit :
$$\rho \frac{\partial u}{\partial t} = k + \mu \frac{\partial^2 u}{\partial y^2}, \qquad u(0,t) = u(H,t) = 0, \qquad u(y,0) = 0.$$

### Solution analytique exacte
Par séparation des variables, la solution exacte est la superposition du profil parabolique asymptotique ($u_{\max,\infty} = \frac{k H^2}{8\mu} = 0.01776\,\mathrm{m/s}$) et d'une série de Fourier amortie exponentiellement :
$$u_{\mathrm{ex}}(y, t) = \frac{k}{2\mu} y(H - y) - \frac{4 k H^2}{\mu \pi^3} \sum_{m=0}^{\infty} \frac{1}{(2m+1)^3} \sin\left(\frac{(2m+1)\pi y}{H}\right) \exp\left(-\nu \left(\frac{(2m+1)\pi}{H}\right)^2 t\right).$$
Avec le schéma d'Euler implicite d'ordre 1 en temps et le schéma centré d'ordre 2 en espace, l'erreur globale obéit à la loi de compétition asymptotique :
$$\epsilon_{\mathrm{EULER}}(\Delta t, \Delta y) = \alpha_1 \Delta t + \beta_1 \Delta t \Delta y + \gamma \Delta y^2,$$
tandis qu'avec le schéma de Gear/BDF2 d'ordre 2 en temps, le terme temporel dominant devient quadratique :
$$\epsilon_{\mathrm{GEAR}}(\Delta t, \Delta y) = \alpha_2 \Delta t^2 + \gamma \Delta y^2.$$

---

---

## 3. Paramètres Physiques, Numériques et Discrétisation

Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine | $[0, L_x] \times [0, L_y]$ | $L_x = 0.16\,\mathrm{m}$, $L_y = H = 0.16\,\mathrm{m}$ ($e = 0.08\,\mathrm{m}$) |
| Propriétés et forçage | $\rho, \mu, g_x$ | $\rho = 1.0\,\mathrm{kg/m}^3$, $\mu = 0.01\,\mathrm{Pa\cdot s}$, $g_x = 0.0555\,\mathrm{m/s}^2$, $g_y = 0.0$ |
| Conditions aux limites | Gauche/Droite ; Bas/Haut | `PERIODIC` en $x$ ; `DIRICHLET` ($u=v=0$) en $y=0$ et $y=H$ |
| Campagne 1 (Matrice $N_y \times \Delta t$) | $t_{\mathrm{final}} = 1.0\,\mathrm{s}$ | $N_y \in \{16, 32, 64, 128, 256\} \times \Delta t \in \{10^{-1}, 10^{-2}, 10^{-3}, 10^{-4}\}\,\mathrm{s}$ ($N_x=16$) |
| Campagne 2 (Profils transitoires) | $t_{\mathrm{final}} = 2.0\,\mathrm{s}$ | $N_y = 64$, $\Delta t \in \{0.2, 0.1, 0.05, 0.01\}\,\mathrm{s}$, instants $t \in \{0.5, 1.0, 1.5, 2.0\}\,\mathrm{s}$ |
| Schémas et Solveurs | Temps, Couplage, Krylov | `EULER` vs `GEAR` ; `AL+BICGSTAB2`, `AL+DIRECT`, `AAL+BICGSTAB2`, `MONO+BICGSTAB2`, `MONO+DIRECT` |

---

---

## 4. Fichier de Configuration JSON Pas à Pas

de Fichiers `.json` (C++17 et Fortran 2018)

### Prototype `case/cpp/cas_4_unsteady_al_bicgstab2_ny64_dt0_01.json`
```json
{
  "mode": "NAVIER_STOKES",
  "test_case": "POISEUILLE_unsteady",
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
    "r_al": 0.0,
    "uzawa_iter": 1,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "output": {
    "output_dir": "VALIDATION_TOTALE/validation_cas_4_poiseuille_instat/output/cpp/al_bicgstab2/ny64_dt0_01",
    "case_name": "cas_4_unsteady_al_bicgstab2_ny64_dt0_01",
    "vtk_interval": 25,
    "log_interval": 10,
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
cd VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_4_poiseuille_instat

# 2. Exécution avec le moteur C++17 (Séquentiel)
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_4_config.json

# 3. Exécution avec le moteur C++17 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_4_config.json

# 4. Exécution avec le moteur Fortran 2018 (Séquentiel)
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_4_config.json

# 5. Exécution avec le moteur Fortran 2018 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_4_config.json

```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification et Tolérances
à Mener (Précision, Intégrité et Performance)

1. **Matrice de convergence spatio-temporelle $5 \times 4$ à $t = 1.0\,\mathrm{s}$ :**
   - Pour $\Delta t = 10^{-4}\,\mathrm{s}$, vérifier que l'erreur spatiale décroît en $\mathcal{O}(\Delta y^2)$ de $6.68\times 10^{-5}$ ($N_y=16$) à $4.21\times 10^{-6}$ ($N_y=64$) et $2.84\times 10^{-7}$ ($N_y=256$).
   - Vérifier que l'ordre temporel apparent de `EULER` entre $\Delta t = 0.1\,\mathrm{s}$ et $\Delta t = 10^{-4}\,\mathrm{s}$ croît de $0.21$ ($N_y=16$, saturé par $\Delta y^2$) à $0.98$ ($N_y=256$, régime asymptotique d'ordre 1).
2. **Vérification comparative des schémas temporels `EULER` vs `GEAR` (BDF2) :**
   - Sur la grille fine $N_y = 256$, exécuter $\Delta t \in \{0.1, 0.05, 0.025, 0.01\}\,\mathrm{s}$ avec `"time_scheme": "GEAR"` et vérifier que l'erreur temporelle décroît avec un ordre proche de $2$ et divise l'erreur d'Euler par un facteur $> 10$ dès $\Delta t = 0.01\,\mathrm{s}$.
3. **Performance du cache de factorisation de `DirectSolver` et des préconditionneurs :**
   - Vérifier que `DirectSolver` réutilise la factorisation LU mise en cache lorsque $\Delta t$ et $\rho, \mu$ sont constants (`is_factored == true`), accélérant les $10\,000$ pas de temps à $\Delta t = 10^{-4}\,\mathrm{s}$.
   - Vérifier l'invariance des 4 combinaisons (`AL+BICGSTAB2`, `AL+DIRECT`, `MONO+BICGSTAB2`, `MONO+DIRECT`) et la parité C++/Fortran ($|E_{L_1}^{\mathrm{C++}} - E_{L_1}^{\mathrm{Fortran}}| < 10^{-15}$).

---

### 6.2 Résultats Attendus et Visualisation ParaView

