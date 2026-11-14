# Tutoriel 10 : Compression d'Air par Injection d'Eau (Piston Compressible Bas Mach)

Ce tutoriel valide la méthode du Lagrangien Augmenté Compressible (C-AAL) sur un cas d'injection d'un liquide incompressible comprimant une poche de gaz à bas nombre de Mach. La montée en pression est comparée à la loi analytique isotherme/adiabatique de compression de gaz parfait.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements physiques et les équations aux dérivées partielles gouvernant l'écoulement.
- Configurer les paramètres numériques et conditions aux limites sur grille cartésienne décalée MAC.
- Exécuter le cas test avec les exécutables C++17 et Fortran 2018 (séquentiel et OpenMP).
- Quantifier avec rigueur la précision par rapport à la solution analytique exacte ($L_2$, $L_\infty$).
- Produire et exploiter les visualisations scientifiques volumiques sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 :** `Mesh2D`, `FlowFields`, `CAALSolver`, `CompressibleGasEOS`, `UzawaPressureUpdate`.
- **Moteur Fortran 2018 :** `module_mesh`, `module_flow_fields`, `module_caal`, `module_gas_eos`, `module_uzawa`.

---

## 2. Formulation Physique et Solution Analytique Exacte

Équations Physiques et Configuration

Ce onzième cas valide l'extension du modèle diphasique 1-fluide aux écoulements faiblement compressibles à bas nombre de Mach (`compressible: true`), dans lesquels un liquide quasi-incompressible comprime un gaz compressible au sein d'une enceinte étanche.

Le système gouvernant couple l'équation de continuité compressible isotherme $\chi_T(C) = C\,\chi_{T,1} + (1-C)\,\chi_{T,2}$, l'équation de quantité de mouvement et l'advection de la fonction couleur avec terme de dilatation $C\,\nabla\cdot\mathbf{u}$ :
$$\chi_T(C)\frac{\partial p}{\partial t} + \nabla \cdot \mathbf{u} = 0, \qquad \rho(C)\left(\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u}\right) = -\nabla p + \nabla\cdot\left[\mu(C)\left(\nabla\mathbf{u} + \nabla\mathbf{u}^T - \frac{2}{3}(\nabla\cdot\mathbf{u})\mathbf{I}\right)\right],$$
$$\frac{\partial C}{\partial t} + \nabla\cdot(\mathbf{u}C) = C\,\nabla\cdot\mathbf{u}, \qquad \rho_i(p) = \rho_{i,0}\left[1 + \chi_{T,i}\,(p - p_0)\right].$$

#### Configuration géométrique et loi analytique exacte
Dans une cavité carrée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ sans gravité ($\mathbf{g}=\mathbf{0}$) ni capillarité ($\sigma=0$), l'eau ($C=1$) occupe initialement $y \in [0, 0.025]\,\mathrm{m}$, surmontée d'une hauteur $H_0 = 0.075\,\mathrm{m}$ d'air ($C=0$) à $p_0 = 101\,325\,\mathrm{Pa}$. De l'eau est injectée par la face inférieure $y=0$ à la vitesse verticale constante $V_0 = 0.1\,\mathrm{m/s}$, les trois autres parois étant à glissement libre (`FREE_SLIP`). L'interface plane monte selon $y_{\mathrm{int}}(t) = 0.025 + V_0 t$, et la masse volumique ainsi que la pression de l'air obéissent à la loi exacte :
$$\rho_{\mathrm{air}}(t) = \frac{\rho_{\mathrm{air},0}}{1 - \dfrac{V_0\,t}{H_0}}, \qquad p(t) = \frac{p_0}{1 - \dfrac{V_0\,t}{H_0}}.$$
Pour $\rho_{\mathrm{air},0} = 1.176829268\,\mathrm{kg/m}^3$, on obtient à $t = 0.3\,\mathrm{s}$ : $\rho_{\mathrm{th}} = 1.9613821138\,\mathrm{kg/m}^3$ ($p_{\mathrm{th}} = 168\,875\,\mathrm{Pa}$) ; et à $t = 0.5\,\mathrm{s}$ : $\rho_{\mathrm{th}} = 3.530487804\,\mathrm{kg/m}^3$ ($p_{\mathrm{th}} = 303\,975\,\mathrm{Pa}$).

---

---

## 3. Paramètres Physiques, Numériques et Discrétisation

Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et interface initiale | $[0, L_x]\times[0, L_y]$, $h_0$, $H_0$ | $L_x = L_y = 0.1\,\mathrm{m}$, eau sur $y\in[0, 0.025]\,\mathrm{m}$, air sur $H_0 = 0.075\,\mathrm{m}$ |
| Fluide 1 (Eau liquide) | $\rho_{1,0}, \mu_1, \chi_{T,1}$ | $\rho_{1,0} = 1000.0\,\mathrm{kg/m}^3$, $\mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$, $\chi_{T,1} = 0.444\times 10^{-9}\,\mathrm{Pa}^{-1}$ |
| Fluide 2 (Air compressible) | $\rho_{2,0}, \mu_2, \chi_{T,2}$ | $\rho_{2,0} = 1.176829268\,\mathrm{kg/m}^3$, $\mu_2 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$, $\chi_{T,2} = 9.8692327\times 10^{-6}\,\mathrm{Pa}^{-1}$ |
| Conditions aux limites | Bas ($y=0$), Haut/Gauche/Droite | Injection `DIRICHLET` $(u=0, v=V_0=0.1\,\mathrm{m/s})$ en bas ; `FREE_SLIP` sur les 3 autres bords |
| Séquence de grilles et $\Delta t$ | $N \in \{32, 64, 128, 256, 512\}$ | $\Delta t = 10^{-3}\,\mathrm{s}$ ($t_{\mathrm{final}} = 0.5\,\mathrm{s}$, soit $500$ pas de temps ; contrôle intermédiaire à $t = 0.3\,\mathrm{s}$) |
| Solveurs et schémas VOF | `solver` \& `vof` | `MONOLITHIC` ($r_{\mathrm{al}}=0$ requis en compressible) avec `DIRECT`, `BICGSTAB2`, `GMRES` ; `WEYMOUTH_YUE`, `CIAM`, `VOF_TVD` |

---

---

## 4. Fichier de Configuration JSON Pas à Pas

de Fichiers `.json` (C++17 et Fortran 2018)

### Prototype `case/cpp/cas_10_piston_compressible_n128.json`
*(Pour le moteur Fortran, remplacer `output/cpp/` par `output/fortran/` dans `output_dir`)*
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
    "enable_two_phase": true,
    "compressible": true,
    "enable_convection": true,
    "enable_surface_tension": false,
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": { "rho": 1000.0, "mu": 1.0e-3, "chi_t": 4.44e-10 },
    "fluid2": { "rho": 1.176829268, "mu": 1.85e-5, "chi_t": 9.8692327e-6 }
  },
  "initial_condition": {
    "type": "DAM_BREAK",
    "box_x0": 0.0,
    "box_y0": 0.0,
    "box_x1": 0.1,
    "box_y1": 0.025,
    "p_init": 101325.0,
    "u0": 0.0,
    "v0": 0.0
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
    "tvd_limiter": "SUPERBEE",
    "cfl_vof": 0.5
  },
  "output": {
    "output_dir": "VALIDATION_TOTALE/validation_cas_10_piston_compressible/output/cpp/n128",
    "case_name": "cas_10_piston_compressible_n128",
    "vtk_interval": 50,
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
cd VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_10_piston_compressible

# 2. Exécution avec le moteur C++17 (Séquentiel)
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_10_config.json

# 3. Exécution avec le moteur C++17 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_10_config.json

# 4. Exécution avec le moteur Fortran 2018 (Séquentiel)
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_10_config.json

# 5. Exécution avec le moteur Fortran 2018 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_10_config.json

```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification et Tolérances
à Mener (Précision, Intégrité et Performance)

1. **Précision sur la loi de compression analytique $\rho_{\mathrm{air}}(t)$ et $p(t)$ :** Extraire la masse volumique moyenne du gaz $\overline{\rho}_{\mathrm{air}}(t)$ et la pression moyenne $\overline{p}_{\mathrm{air}}(t)$ dans la zone $C < 0.01$ aux instants $t \in \{0.1, 0.2, 0.3, 0.4, 0.5\}\,\mathrm{s}$. Vérifier qu'à $t = 0.3\,\mathrm{s}$, $\overline{\rho}_{\mathrm{air}}$ converge vers $\rho_{\mathrm{th}} = 1.9613821138\,\mathrm{kg/m}^3$ avec une erreur relative décroissant de $0.50\,\%$ ($32^2$) à $0.031\,\%$ ($512^2$), soit un ordre de convergence spatial $p = 1.00$. Vérifier qu'à $t = 0.5\,\mathrm{s}$, $\overline{\rho}_{\mathrm{air}} \approx 3.5305\,\mathrm{kg/m}^3$ et $\overline{p}_{\mathrm{air}} \approx 303\,967\,\mathrm{Pa}$ (écart $< 0.01\,\%$).
2. **Intégrité du bloc diagonal compressible $C_p = -\frac{\chi_T}{\Delta t}I$ et des solveurs :** En mode `compressible: true`, le bloc $(2,2)$ du système selle vaut $C_{p,i,j} = -(C_{i,j}\chi_{T,1} + (1-C_{i,j})\chi_{T,2})/\Delta t$, ce qui lève la singularité de jauge de pression sans ancrage artificiel. Vérifier que `MONOLITHIC` ($r_{\mathrm{al}}=0$) résolu par `DIRECT`, `BICGSTAB2`, `BICGSTAB` et `GMRES` avec `BLOCK_TRIANGULAR`, `ILU0` et `SSOR` restitue exactement la même courbe $p(t)$ à $10^{-8}\,\mathrm{Pa}$ près.
3. **Comparaison des schémas VOF avec dilatation $C\nabla\cdot\mathbf{u}$ et parité dual-code :** Comparer `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC` et `VOF_TVD` sur la planéité de l'interface $y_{\mathrm{int}}(t) = 0.025 + 0.1 t$ et la conservation de la masse totale d'air $M_{\mathrm{air}}(t) = M_{\mathrm{air}}(0)$. Vérifier la parité entre `covof2d_cpp` et `covof2d_fortran` ($\Delta_{\mathrm{dual}}(p) < 10^{-8}\,\mathrm{Pa}$).

---

### 6.2 Résultats Attendus et Visualisation ParaView

