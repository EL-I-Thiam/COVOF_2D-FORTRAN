# Dictionnaire Encyclopédique des Paramètres d'Entrée

## Structure Générale du Fichier de Configuration JSON

COVOF 2D utilise un format de fichier d'entrée textuel standardisé en notation **JSON** (*JavaScript Object Notation*), tolérant les commentaires sur une seule ligne préfixés par `//`. Ce format offre une lisibilité parfaite, une hiérarchisation arborescente naturelle et une compatibilité immédiate avec les outils d'automatisation et d'exploration paramétrique en Python.

Le fichier d'entrée est structuré en **10 blocs fonctionnels principaux** :

1. `mode` et métadonnées globales : Choix du type de problème physique et nom du cas.
2. `mesh` : Définition de la grille cartésienne MAC décalée.
3. `time` : Paramètres d'intégration temporelle et de pas de temps adaptatif.
4. `physics` : Propriétés hydrodynamiques et thermodynamiques des fluides (phases 1 et 2).
5. `thermal` : Paramètres de l'équation d'énergie et couplage de Boussinesq.
6. `solver` : Méthode de couplage vitesse-pression, solveurs de Krylov et préconditionneurs.
7. `convection` : Schémas spatiaux et linéarisation de l'advection.
8. `obstacles` : Définition des corps solides ou poreux immergés par pénalisation de Brinkman.
9. `boundary_conditions` : Conditions aux limites sur les quatre frontières cardinales.
10. `output` : Fréquence et format des sorties volumiques VTK et journaux de convergence.

---

## Dictionnaire Exhaustif des Paramètres par Section

### 1. Paramètres Généraux et Régime Physique

| Clé JSON | Type | Défaut | Options / Plage | Description Physique et Numérique |
| :------------ | :------ | :-------- | :--------------------- | :--------------------------------------------------- |
| `mode` | `string` | `"NAVIER_STOKES"` | `"NAVIER_STOKES"`, `"TWO_PHASE"`, `"RAYLEIGH_BENARD"`, `"LAPLACIAN"`, `"VORTEX_ADVECTION"` | Définit le système d'équations résolu. `"LAPLACIAN"` résout l'équation elliptique pure $\Delta \phi = f$, `"VORTEX_ADVECTION"` active le transport cinématique pur de l'interface VOF. |
| `test_case` | `string` | `"SIMULATION"` | Texte libre | Chaîne descriptive du cas-test pour les en-têtes de journaux et métadonnées de sortie. |
| `steady_state_tol` | `double` | `0.0` | $[0.0, 1.0]$ | Seuil d'arrêt sur régime stationnaire : arrêt si $\frac{\|\mathbf{u}^{n+1} - \mathbf{u}^n\|_\infty}{\Delta t} \le \text{tol}$. Mettre `0.0` pour désactiver. |
| `source_term` | `double` | `0.0` | Réel quelconque | Terme source volumique uniforme utilisé en mode `"LAPLACIAN"`. |

### 2. Section `mesh` (Maillage Cartésien MAC)

| Clé JSON | Type | Défaut | Options / Plage | Description Physique et Numérique |
| :------------ | :------ | :-------- | :--------------------- | :--------------------------------------------------- |
| `mesh.nx` | `integer` | `100` | $\ge 4$ | Nombre total de cellules de calcul selon l'axe horizontal $X$. |
| `mesh.ny` | `integer` | `100` | $\ge 4$ | Nombre total de cellules de calcul selon l'axe vertical $Y$. |
| `mesh.lx` | `double` | `1.0` | $> 0.0\,\text{m}$ | Longueur physique totale du domaine de calcul selon $X$. |
| `mesh.ly` | `double` | `1.0` | $> 0.0\,\text{m}$ | Hauteur physique totale du domaine de calcul selon $Y$. |
| `mesh.x_min` | `double` | `0.0` | Réel quelconque | Coordonnée minimale d'origine selon $X$ $[\text{m}]$. |
| `mesh.y_min` | `double` | `0.0` | Réel quelconque | Coordonnée minimale d'origine selon $Y$ $[\text{m}]$. |
| `mesh.grid_type` | `string` | `"CARTESIAN"` | `"CARTESIAN"` | Type de maillage. Actuellement seul le maillage cartésien décalé est supporté. |

### 3. Section `time` (Intégration Temporelle)

| Clé JSON | Type | Défaut | Options / Plage | Description Physique et Numérique |
| :------------ | :------ | :-------- | :--------------------- | :--------------------------------------------------- |
| `time.dt` | `double` | `0.001` | $> 0.0\,\text{s}$ | Pas de temps d'intégration initial ou fixe $[\text{s}]$. |
| `time.t_end` | `double` | `1.0` | $> 0.0\,\text{s}$ | Temps physique final d'arrêt de la simulation $[\text{s}]$. |
| `time.cfl` | `double` | `0.5` | $]0.0, 1.0]$ | Nombre de Courant-Friedrichs-Lewy maximal autorisé pour le pas adaptatif. |
| `time.time_scheme` | `string` | `"EULER_EXPLICIT"` | `"EULER_EXPLICIT"`, `"GEAR"`, `"BDF2"`, `"SSP_RK3"` | Schéma d'intégration temporelle. `"GEAR"` (ordre 2 implicite stable) est recommandé pour Navier-Stokes. |
| `time.adaptative_dt` | `boolean` | `false` | `true`, `false` | Active le calcul dynamique du pas de temps sous contraintes CFL convective, visqueuse et capillaire. |

### 4. Section `physics` (Propriétés des Fluides et Interfaces)

| Clé JSON | Type | Défaut | Options / Plage | Description Physique et Numérique |
| :------------ | :------ | :-------- | :--------------------- | :--------------------------------------------------- |
| `physics.rho1` | `double` | `1.0` | $> 0.0\,\text{kg}\cdot\text{m}^{-3}$ | Masse volumique de la Phase 1 (fluide porteur ou liquide). |
| `physics.rho2` | `double` | `1.0` | $> 0.0\,\text{kg}\cdot\text{m}^{-3}$ | Masse volumique de la Phase 2 (goutte, bulle ou gaz). |
| `physics.mu1` | `double` | `0.01` | $> 0.0\,\text{Pa}\cdot\text{s}$ | Viscosité dynamique de la Phase 1. |
| `physics.mu2` | `double` | `0.01` | $> 0.0\,\text{Pa}\cdot\text{s}$ | Viscosité dynamique de la Phase 2. |
| `physics.sigma` | `double` | `0.0` | $\ge 0.0\,\text{N}\cdot\text{m}^{-1}$ | Coefficient de tension superficielle interfaciale (modèle CSF). |
| `physics.gx` | `double` | `0.0` | Réel $[\text{m}\cdot\text{s}^{-2}]$ | Composante horizontale du vecteur gravité. |
| `physics.gy` | `double` | `0.0` | Réel $[\text{m}\cdot\text{s}^{-2}]$ | Composante verticale du vecteur gravité (ex. `-9.81`). |
| `physics.enable_two_phase` | `boolean` | `false` | `true`, `false` | Active le modèle diphasique et le transport VOF de l'interface. |
| `physics.viscosity_model` | `string` | `"CONSTANT"` | `"CONSTANT"`, `"VARIABLE"` | Modèle de viscosité. `"VARIABLE"` active une viscosité non uniforme $\mu(y)$. |
| `physics.vof_model` | `string` | `"PLIC"` | `"PLIC"` | Modèle de suivi d'interface VOF. Reconstruction linéaire par morceaux. |
| `physics.compressible` | `boolean` | `false` | `true`, `false` | Active la formulation barotrope bas-Mach All-Mach ($\chi_T > 0$). |
| `physics.gamma` | `double` | `1.4` | $> 1.0$ | Rapport des chaleurs spécifiques $\gamma = C_p / C_v$ pour l'air ou le gaz. |
| `physics.p_inf` | `double` | `0.0` | $\ge 0.0\,\text{Pa}$ | Pression de cohésion pour l'équation d'état *Stiffened Gas*. |
| `physics.ic_type` | `string` | `"NONE"` | `"NONE"`, `"CIRCLE"`, `"RECTANGLE"`, `"DAM_BREAK"`, `"TAYLOR_GREEN"`, `"VORTEX"` | Géométrie de la condition initiale diphasique de fraction volumique $C(x, y)$. |
| `physics.ic_x0`, `ic_y0` | `double` | `0.0` | Coordonnées $[\text{m}]$ | Centre du disque ou coin inférieur gauche du rectangle initial. |
| `physics.ic_radius` | `double` | `0.1` | $> 0.0\,\text{m}$ | Rayon de la goutte ou bulle initiale (si `ic_type="CIRCLE"`). |
| `physics.ic_width`, `ic_height` | `double` | `0.1` | $> 0.0\,\text{m}$ | Largeur et hauteur du bloc fluide initial (si `ic_type="RECTANGLE"`). |

### 5. Section `thermal` (Énergie Thermique et Convection Naturelle)

| Clé JSON | Type | Défaut | Options / Plage | Description Physique et Numérique |
| :------------ | :------ | :-------- | :--------------------- | :--------------------------------------------------- |
| `thermal.enable_energy` | `boolean` | `false` | `true`, `false` | Active l'équation de transport-diffusion de la température. |
| `thermal.th_cond1` | `double` | `0.026` | $> 0.0\,\text{W}\cdot\text{m}^{-1}\cdot\text{K}^{-1}$ | Conductivité thermique de la Phase 1. |
| `thermal.th_cond2` | `double` | `0.026` | $> 0.0\,\text{W}\cdot\text{m}^{-1}\cdot\text{K}^{-1}$ | Conductivité thermique de la Phase 2. |
| `thermal.cp1` | `double` | `1005.0` | $> 0.0\,\text{J}\cdot\text{kg}^{-1}\cdot\text{K}^{-1}$ | Capacité thermique massique isobare de la Phase 1. |
| `thermal.cp2` | `double` | `1005.0` | $> 0.0\,\text{J}\cdot\text{kg}^{-1}\cdot\text{K}^{-1}$ | Capacité thermique massique isobare de la Phase 2. |
| `thermal.beta_th` | `double` | `0.0034` | $> 0.0\,\text{K}^{-1}$ | Coefficient de dilatation volumique thermique (approximation de Boussinesq). |
| `thermal.t_ref` | `double` | `293.15` | $> 0.0\,\text{K}$ | Température de référence pour la force de flottabilité. |
| `thermal.t_hot`, `t_cold` | `double` | `300.0`, `280.0` | $> 0.0\,\text{K}$ | Températures imposées sur les parois chaude et froide. |
| `thermal.t_init` | `double` | `293.15` | $> 0.0\,\text{K}$ | Température initiale uniforme dans le domaine de calcul. |

### 6. Section `solver` (Couplage Navier-Stokes et Solveurs Linéaires)

| Clé JSON | Type | Défaut | Options / Plage | Description Physique et Numérique |
| :------------ | :------ | :-------- | :--------------------- | :--------------------------------------------------- |
| `solver.type` | `string` | `"AUGMENTED_LAGRANGIAN"` | `"AUGMENTED_LAGRANGIAN"`, `"AAL"`, `"C-AAL"`, `"MONOLITHIC"`, `"STANDARD_PROJECTION"` | Formulation de couplage vitesse-pression. `"AAL"` est le mode adaptatif incompressible pur ; `"C-AAL"` active la thermodynamique bas-Mach ; `"MONOLITHIC"` résout le point-selle couplé. |
| `solver.linear_solver` | `string` | `"BICGSTAB2"` | `"BICGSTAB2"`, `"BICGSTAB"`, `"CG"`, `"GMRES"`, `"MINRES"`, `"DIRECT"` | Algorithme de résolution du système linéaire. `"MINRES"` est optimal pour opérateurs symétriques indéfinis. |
| `solver.preconditioner` | `string` | `"JACOBI"` | `"JACOBI"`, `"SSOR"`, `"ILU0"`, `"BLOCK_TRIANGULAR"`, `"CUSTOM"`, `"NONE"` | Préconditionneur linéaire. `"BLOCK_TRIANGULAR"` est recommandé pour `"MONOLITHIC"`. |
| `solver.tol` | `double` | `1.0e-5` | $]0.0, 1.0e-2]$ | Tolérance d'arrêt sur le résidu relatif du solveur linéaire de Krylov. |
| `solver.max_iter` | `integer` | `300` | $\ge 10$ | Nombre maximal d'itérations par résolution de système linéaire. |
| `solver.r_al` | `double` | `100.0` | $> 0.0$ | Paramètre scalaire de pénalité de Lagrangien augmenté standard $[\text{Pa}\cdot\text{s}\cdot\text{m}^{-1}]$. |
| `solver.adaptative_al` | `boolean` | `false` | `true`, `false` | Active le calcul adaptatif local de la pénalité $r(t, \mathbf{x})$ selon Vincent et al. (2004). |
| `solver.k_aal` | `double` | `1.0` | $]0.0, 10.0]$ | Facteur multiplicatif d'échelle pour la pénalité AAL / C-AAL ($K_{\text{aal}}$). |
| `solver.uzawa_iter` | `integer` | `1` | $\ge 1$ | Nombre maximal d'itérations dans la boucle d'Uzawa par pas de temps. |
| `solver.uzawa_tol` | `double` | `1.0e-4` | $> 0.0$ | Seuil d'arrêt d'Uzawa sur la divergence maximale $\|\nabla \cdot \mathbf{u}\|_\infty$. |
| `solver.use_precond` | `boolean` | `true` | `true`, `false` | Active ou désactive l'application du préconditionneur linéaire. |

### 7. Section `convection` (Discrétisation Spatiale des Flux)

| Clé JSON | Type | Défaut | Options / Plage | Description Physique et Numérique |
| :------------ | :------ | :-------- | :--------------------- | :--------------------------------------------------- |
| `convection.convective_scheme` | `string` | `"UPWIND"` | `"UPWIND"`, `"CENTERED"` | Schéma de flux d'advection spatiale. `"CENTERED"` est d'ordre 2 mais nécessite un maillage fin. |
| `convection.convective_linearization`| `string` | `"UPWIND"` | `"UPWIND"`, `"CENTERED"`, `"NEWTON"`, `"ADAMS_BASHFORTH"` | Linéarisation du terme convectif. `"NEWTON"` assure une convergence quadratique en monolithique. |

### 8. Section `obstacles` (Pénalisation de Brinkman pour Corps Immergés)

| Clé JSON | Type | Défaut | Options / Plage | Description Physique et Numérique |
| :------------ | :------ | :-------- | :--------------------- | :--------------------------------------------------- |
| `obstacles.enable_obstacles`| `boolean` | `false` | `true`, `false` | Active la présence d'obstacles solides ou poreux immergés. |
| `obstacles.type` | `string` | `"CYLINDER"` | `"CYLINDER"`, `"RECTANGLE"` | Géométrie de l'obstacle. |
| `obstacles.xc`, `yc` | `double` | `0.5`, `0.5` | Coordonnées $[\text{m}]$ | Centre du cylindre immergé (si `type="CYLINDER"`). |
| `obstacles.radius` | `double` | `0.1` | $> 0.0\,\text{m}$ | Rayon du cylindre immergé. |
| `obstacles.x0`, `x1` | `double` | `0.2`, `0.4` | Coordonnées $[\text{m}]$ | Limites en $X$ du rectangle immergé (si `type="RECTANGLE"`). |
| `obstacles.y0`, `y1` | `double` | `0.0`, `0.3` | Coordonnées $[\text{m}]$ | Limites en $Y$ du rectangle immergé. |
| `obstacles.permeability` | `double` | `1.0e-40` | $\ge 0.0\,\text{m}^2$ | Perméabilité de Darcy $K$. Valeur $\le 10^{-12}$ correspond à un solide imperméable adhérent. |

### 9. Section `boundary_conditions` (Frontières du Domaine)

Chaque frontière (`left`, `right`, `bottom`, `top`) est configurée indépendamment :

| Clé JSON | Type | Défaut | Options / Plage | Description Physique et Numérique |
| :------------ | :------ | :-------- | :--------------------- | :--------------------------------------------------- |
| `boundary_conditions.<face>.type` | `string` | `"DIRICHLET"` | `"DIRICHLET"`, `"NEUMANN"`, `"PERIODIC"`, `"FREESLIP"` | Type de condition aux limites hydrodynamique. |
| `boundary_conditions.<face>.u` | `double` | `0.0` | Réel $[\text{m}\cdot\text{s}^{-1}]$ | Valeur de vitesse horizontale prescrite si Dirichlet. |
| `boundary_conditions.<face>.v` | `double` | `0.0` | Réel $[\text{m}\cdot\text{s}^{-1}]$ | Valeur de vitesse verticale prescrite si Dirichlet. |
| `boundary_conditions.<face>.t` | `double` | `293.15` | Réel $[\text{K}]$ | Valeur de température prescrite si Dirichlet thermique. |

### 10. Section `output` (Visualisation et Fichiers de Données)

| Clé JSON | Type | Défaut | Options / Plage | Description Physique et Numérique |
| :------------ | :------ | :-------- | :--------------------- | :--------------------------------------------------- |
| `output.format` | `string` | `"VTK"` | `"VTK"`, `"CSV"`, `"NONE"` | Format des résultats de calcul. `"VTK"` génère `.vtr` et `.pvtr`. |
| `output.frequency` | `integer` | `10` | $\ge 1$ | Fréquence d'écriture des fichiers de résultats (toutes les $N$ itérations). |
| `output.directory` | `string` | `"output"` | Chemin valide | Répertoire de destination des fichiers de résultats. |
| `output.filename` | `string` | `"result"` | Nom valide | Préfixe des fichiers numériques de sortie (ex. `result_00010.vtr`). |
| `output.save_pvd` | `boolean` | `true` | `true`, `false` | Génère le fichier maître `.pvd` liant la série temporelle pour ParaView. |

---

## Fichier de Configuration Maître Intégral et Commenté

Voici la reproduction complète du fichier de référence maître (`template_parametres_covof2d.json`) utilisable comme canevas pour tout nouveau cas de simulation :

```json
{
  // =========================================================================================
  // COVOF 2D — FICHIER DE CONFIGURATION MAÎTRE ET CATALOGUE EXHAUSTIF DES PARAMÈTRES
  // =========================================================================================
  "mode": "NAVIER_STOKES",
  "test_case": "SIMULATION_MAITRE_COVOF2D",
  "steady_state_tol": 0.0,
  "source_term": 0.0,

  // -----------------------------------------------------------------------------------------
  // 1. MAILLAGE CARTÉSIEN DÉCALÉ (MAC)
  // -----------------------------------------------------------------------------------------
  "mesh": {
    "nx": 64,
    "ny": 128,
    "lx": 0.1,
    "ly": 0.2,
    "x_min": 0.0,
    "y_min": 0.0,
    "grid_type": "CARTESIAN"
  },

  // -----------------------------------------------------------------------------------------
  // 2. DISCRÉTISATION TEMPORELLE ET CFL ADAPTATIF
  // -----------------------------------------------------------------------------------------
  "time": {
    "dt": 0.0005,
    "t_end": 0.5,
    "cfl": 0.5,
    "time_scheme": "GEAR",
    "adaptative_dt": true
  },

  // -----------------------------------------------------------------------------------------
  // 3. PROPRIÉTÉS PHYSIQUES DES PHASES ET TENSION SUPERFICIELLE
  // -----------------------------------------------------------------------------------------
  "physics": {
    "rho1": 1000.0,
    "rho2": 1.2,
    "mu1": 1.0e-3,
    "mu2": 1.8e-5,
    "sigma": 0.0728,
    "gx": 0.0,
    "gy": -9.81,
    "enable_two_phase": true,
    "viscosity_model": "CONSTANT",
    "vof_model": "PLIC",
    "compressible": false,
    "ic_type": "DAM_BREAK",
    "ic_x0": 0.0,
    "ic_y0": 0.0,
    "ic_width": 0.035,
    "ic_height": 0.07
  },

  // -----------------------------------------------------------------------------------------
  // 4. TRANSFERT THERMIQUE (ÉNERGIE)
  // -----------------------------------------------------------------------------------------
  "thermal": {
    "enable_energy": false,
    "th_cond1": 0.6,
    "th_cond2": 0.026,
    "cp1": 4182.0,
    "cp2": 1005.0,
    "beta_th": 2.1e-4,
    "t_ref": 293.15,
    "t_hot": 303.15,
    "t_cold": 283.15,
    "t_init": 293.15
  },

  // -----------------------------------------------------------------------------------------
  // 5. COUPLAGE NAVIER-STOKES ET SOLVEURS LINÉAIRES
  // -----------------------------------------------------------------------------------------
  "solver": {
    "type": "AAL",
    "linear_solver": "BICGSTAB2",
    "preconditioner": "JACOBI",
    "tol": 1.0e-6,
    "max_iter": 500,
    "r_al": 100.0,
    "adaptative_al": true,
    "k_aal": 1.0,
    "uzawa_iter": 3,
    "uzawa_tol": 1.0e-5,
    "use_precond": true
  },

  // -----------------------------------------------------------------------------------------
  // 6. SCHÉMAS CONVECTIFS
  // -----------------------------------------------------------------------------------------
  "convection": {
    "convective_scheme": "UPWIND",
    "convective_linearization": "NEWTON"
  },

  // -----------------------------------------------------------------------------------------
  // 7. CORPS SOLIDES IMMERGÉS (PÉNALISATION DE BRINKMAN)
  // -----------------------------------------------------------------------------------------
  "obstacles": {
    "enable_obstacles": true,
    "type": "CYLINDER",
    "xc": 0.05,
    "yc": 0.04,
    "radius": 0.01,
    "permeability": 1.0e-40
  },

  // -----------------------------------------------------------------------------------------
  // 8. CONDITIONS AUX LIMITES HYDRODYNAMIQUES ET THERMIQUES
  // -----------------------------------------------------------------------------------------
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "right":  { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "FREESLIP",  "u": 0.0, "v": 0.0 }
  },

  // -----------------------------------------------------------------------------------------
  // 9. GESTION DES SORTIES NUMÉRIQUES ET LOGS
  // -----------------------------------------------------------------------------------------
  "output": {
    "format": "VTK",
    "frequency": 20,
    "directory": "output",
    "filename": "dambreak_obs",
    "save_pvd": true
  }
}
```
