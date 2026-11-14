# Architecture Logicielle Hybride C++ / Fortran

## Vision d'Ensemble et Paradigmes de Conception

Le projet **COVOF 2D** repose sur une stratégie de co-développement gémellaire originale au sein de la communauté de mécanique des fluides numérique (CFD). Plutôt que de contraindre un langage unique à assurer à la fois la haute abstraction système et le calcul tensoriel intensif, COVOF 2D met en œuvre deux implémentations parallèles indépendantes partageant exactement les mêmes formulations mathématiques, les mêmes disques de données et le même protocole de communication distribuée MPI.

### Les Paradigmes en Présence

1. **La Branche C++17 (`COVOF_2D-C++`) :**
   - **Paradigme orienté objet moderne :** Modélisation par classes distinctes à responsabilité unique (`Mesh`, `FlowFields`, `NavierStokesSolver`, `InterfaceTracker`, `LinearSolver`, `Simulation`).
   - **Gestion RAII stricte :** L'allocation et la libération des ressources mémoire dynamiques (`std::vector<double>`), des poignées de fichiers (`std::ofstream`) et des communicateurs MPI (`MPI_Comm`) sont strictement encapsulées dans les constructeurs et destructeurs des classes, éliminant tout risque de fuite mémoire (*memory leak*).
   - **Polymorphisme dynamique et généricité :** Les solveurs linéaires de Krylov et les préconditionneurs dérivent de classes de base abstraites (`LinearSolverBase`, `PreconditionerBase`), permettant le basculement dynamique à l'exécution via des pointeurs intelligents `std::unique_ptr` et l'injection aisée de nouvelles heuristiques.

2. **La Branche Fortran 2018 (`COVOF_2D-FORTRAN`) :**
   - **Paradigme modulaire fortement typé :** Découpage en 14 modules autonomes avec encapsulation `private` par défaut et directive stricte `implicit none` sur chaque unité de compilation.
   - **Types dérivés extensibles :** Définition de structures de données explicites (`mesh_2d_t`, `flow_fields_t`, `config_t`, `navier_stokes_t`, `simulation_2d_t`) avec procédures liées aux types (*type-bound procedures*).
   - **Vectorisation native des boucles imbriquées :** Exploitation optimale de l'ordonnancement mémoire en colonne majeure (*column-major order*), permettant au compilateur (`gfortran`, `ifx`) de dérouler les boucles les plus internes ($i$ contigu en mémoire) et de saturer les unités vectorielles SIMD (AVX2 / AVX-512) sans coût d'indirection.

---

## Cartographie Intégrale de l'Arborescence du Dépôt

Le tableau ci-après dresse l'inventaire exhaustif de chaque fichier source des deux dépôts jumeaux, en explicitant son rôle algorithmique précis :

### 1. Fichiers de la Branche C++ (`COVOF_2D-C++`)

| Chemin Relatif | Composant Logiciel | Responsabilité et Contenu Algorithmique |
| :--- | :--- | :--- |
| `include/covof/Types.hpp` | Types & Enums | Définition des énumérations globales : types de conditions aux limites (`BCType`), solveurs (`SolverType`), algorithmes de Krylov (`LinearSolverType`), préconditionneurs (`PrecondType`), schémas convectifs et constantes mathématiques ($\pi, \epsilon$). |
| `include/covof/Config.hpp` | Spécification Configuration | Structure de stockage complète des paramètres utilisateurs désérialisés. |
| `src/Config.cpp` | Analyseur JSON | Analyseur syntaxique natif tolérant les commentaires `//`, extraction sécurisée avec conversion de types et valeurs de repli par défaut. |
| `include/covof/Mesh.hpp` | En-tête Topologie MAC | Définition de la classe `Mesh2D` : coordonnées 1D des faces et centres, pas d'espace $\Delta x, \Delta y$, dimensions locales et globales. |
| `src/Mesh.cpp` | Implémentation Maillage | Calcul des métriques géométriques et initialisation des tableaux d'indices de mailles décalées. |
| `include/covof/Field.hpp` | Champs d'Écoulement | Structure `FlowFields` : stockage SoA (*Structure of Arrays*) contigu des composantes scalaires et vectorielles ($u, v, p, \rho, \mu, C, T, \chi_T, r_{\text{aal}}$). |
| `include/covof/MpiComm.hpp` | En-tête Parallélisme | Définition de la classe `MpiComm` : topologie cartésienne 2D, gestion des halos et réductions. |
| `src/MpiComm.cpp` | Implémentation MPI | Initialisation de `MPI_Cart_create`, échanges asynchrones non bloquants `MPI_Isend`/`MPI_Irecv` des mailles fantômes et réductions globales. |
| `include/covof/CSRMatrix.hpp` | Algèbre Creuse CSR | Définition de la structure `CSRMatrix` : stockage creux compressé par ligne (tableaux `values`, `col_indices`, `row_offsets`) pour opérations spécifiques. |
| `src/CSRMatrix.cpp` | Implémentation CSR | Produit matrice-vecteur creux optimisé et calcul de normes. |
| `include/covof/NavierStokesSolver.hpp` | En-tête Opérateur NS | Définition de la classe `NavierStokesSolver` : assemblage du second membre, application sans matrice de l'opérateur couplé, calcul de divergence et d'Uzawa. |
| `src/NavierStokesSolver.cpp` | Implémentation NS | Cœur du solveur hydrodynamique : pochoirs MAC pour la diffusion visqueuse $\mu_{xy}$, convection (Upwind, Centré, Newton, AB2), pénalisation de Brinkman et Lagrangien augmenté. |
| `include/covof/InterfaceReconstruction.hpp` | En-tête PLIC | Définition de la classe `InterfaceReconstruction` : reconstruction géométrique d'interface. |
| `src/InterfaceReconstruction.cpp` | Implémentation PLIC | Algorithme de Youngs / ELVIRA : calcul des normales unitaires $\mathbf{n}$, constante de plan $\alpha$ et calcul de courbure $\kappa$. |
| `include/covof/InterfaceTracker.hpp` | En-tête VOF | Définition de la classe `InterfaceTracker` : advection eulerienne conservatrice. |
| `src/InterfaceTracker.cpp` | Implémentation VOF | Algorithme de Weymouth-Yue (2010) : advection 1D directionnelle alternée sans perte de masse. |
| `include/covof/Preconditioner.hpp` | En-tête Préconditionneurs | Hiérarchie polymorphe : `JacobiPreconditioner`, `SSORPreconditioner`, `ILU0Preconditioner`, `BlockTriangularPreconditioner`, `CustomPreconditioner`. |
| `src/Preconditioner.cpp` | Implémentation Préconds | Calcul des inverses diagonaux enrichis, balayages de Gauss-Seidel et résolution du complément de Schur. |
| `include/covof/LinearSolver.hpp` | En-tête Solveurs Krylov | Hiérarchie : `BiCGSTAB2Solver`, `BiCGSTABSolver`, `CGSolver`, `GMRESSolver`, `MINRESSolver`, `DirectSolver`. |
| `src/LinearSolver.cpp` | Implémentation Krylov | Algorithmes de Krylov sans matrice avec réductions MPI et factorisation LU LAPACK (`dgetrf`/`dgetrs`). |
| `include/covof/VTKWriter.hpp` | En-tête Sorties | Classe `VTKWriter` pour l'exportation des résultats. |
| `src/VTKWriter.cpp` | Implémentation Sorties | Écriture XML Rectilinear Grid (`.vtr`, `.pvtr`, `.pvd`) et journaux CSV de convergence. |
| `include/covof/Simulation.hpp` | En-tête Orchestration | Définition de la classe `Simulation` : gestion de la boucle temporelle globale. |
| `src/Simulation.cpp` | Implémentation Simulation | Pilotage du pas de temps adaptatif, étapes d'Uzawa, couplage thermique Boussinesq et convergence. |
| `src/main.cpp` | Point d'Entrée Binaire | Programme principal : analyse CLI, initialisation MPI, exécution et mesure de temps CPU. |

### 2. Fichiers de la Branche Fortran (`COVOF_2D-FORTRAN`)

| Chemin Relatif | Module Fortran | Responsabilité et Contenu Algorithmique |
| :--- | :--- | :--- |
| `src/module_precision.f90` | `module_precision` | Définition du paramètre de précision double standardisé `dp = real64` issu de `iso_fortran_env`. |
| `src/module_constants.f90` | `module_constants` | Énumérations globales et constantes mathématiques ($\pi, \epsilon$) strictement typées. |
| `src/module_eos.f90` | `module_eos` | Équations d'état thermodynamiques : gaz parfaits, Stiffened Gas, Tait, calcul de compressibilité $\chi_T$. |
| `src/module_json.f90` | `module_json` | Analyseur syntaxique natif JSON en Fortran pur pour la lecture sans dépendance externe. |
| `src/module_config.f90` | `module_config` | Type dérivé `config_t` et extraction des paramètres du dictionnaire d'entrée. |
| `src/module_mesh.f90` | `module_mesh` | Type dérivé `mesh_2d_t` : discrétisation MAC décalée, coordonnées de mailles et métriques spatiales. |
| `src/module_fields.f90` | `module_fields` | Type dérivé `flow_fields_t` : allocations dynamiques 2D contiguës de tous les champs physiques. |
| `src/module_mpi_comm.f90` | `module_mpi_comm` | Gestionnaire MPI : topologie cartésienne 2D, échanges de cellules fantômes et réductions parallèles. |
| `src/module_interface_reconstruction.f90` | `module_interface_reconstruction` | Reconstruction VOF-PLIC : calcul des vecteurs normaux, constante de plan $\alpha$ et courbure $\kappa$. |
| `src/module_interface_vof.f90` | `module_interface_vof` | Advection géométrique VOF de Weymouth-Yue (2010) avec conservation de volume. |
| `src/module_navier_stokes.f90` | `module_navier_stokes` | Type dérivé `navier_stokes_t` : opérateur de diffusion, convection, AL, Brinkman et mise à jour d'Uzawa. |
| `src/module_linear_solvers.f90` | `module_linear_solvers` | Solveurs de Krylov sans matrice : BiCGSTAB(2), BiCGSTAB, PCG, GMRES, MINRES, Direct LU et préconditionneurs. |
| `src/module_vtk_writer.f90` | `module_vtk_writer` | Écriture parallèle des fichiers XML VTK (`.vtr`, `.pvtr`, `.pvd`) et fichiers CSV de suivi. |
| `src/module_simulation.f90` | `module_simulation` | Type dérivé `simulation_2d_t` : chef d'orchestre de la simulation instationnaire couplée. |
| `src/covof2d_main.f90` | `program covof2d_main` | Programme principal exécutable Fortran 2018. |

---

## Graphe de Dépendances et Cycle d'Exécution Complet

### Graphe de Dépendances Modulaire

Le diagramme suivant illustre la hiérarchie de dépendances entre les composants logiciels, commune aux deux codes :

```
                  +----------------------+
                  |   Types / Constants  |
                  +----------------------+
                             |
                             v
                  +----------------------+
                  |     Mesh / EOS       |
                  +----------------------+
                             |
                             v
         +-------------------+-------------------+
         |                                       |
         v                                       v
+------------------+                   +-------------------+
|  MpiComm / Fields|                   |   Config / JSON   |
+------------------+                   +-------------------+
         |                                       |
         +-------------------+-------------------+
                             |
                             v
         +-------------------+-------------------+
         |                                       |
         v                                       v
+------------------+                   +-------------------+
|   NavierStokes   |                   |  InterfaceTracker |
+------------------+                   +-------------------+
         |                                       |
         v                                       v
+------------------+                   +-------------------+
|   LinearSolver   |                   |     VTKWriter     |
+------------------+                   +-------------------+
         |                                       |
         +-------------------+-------------------+
                             |
                             v
                  +----------------------+
                  |  Simulation (Master) |
                  +----------------------+
                             |
                             v
                  +----------------------+
                  |   Programme Principal|
                  +----------------------+
```

### Cycle d'Exécution d'un Pas de Temps

Le déroulement chronologique d'un pas de temps complet $t^n \to t^{n+1}$ au sein de la classe `Simulation` (ou `simulation_2d_t`) respecte la séquence rigoureuse suivante :

```
+--------------------------------------------------------------------+
| Étape 1 : Calcul du Pas de Temps Adaptatif (CFL)                   |
| - Évaluation des vitesses maximales locales et capillaires         |
| - Réduction globale MPI : dt = min(dt_max, dt_cfl, dt_visc, dt_cap)|
+--------------------------------------------------------------------+
                                   |
                                   v
+--------------------------------------------------------------------+
| Étape 2 : Avancement Temporel de l'Interface VOF                   |
| - Advection conservative de C(x, y) via Weymouth-Yue               |
| - Reconstruction VOF-PLIC : calcul des normales et courbure kappa  |
| - Propriétés volumiques rho(C), mu(C) et force CSF f_sigma         |
| - Synchronisation MPI des mailles fantômes de rho, mu et f_sigma   |
+--------------------------------------------------------------------+
                                   |
                                   v
+--------------------------------------------------------------------+
| Étape 3 : Calcul des Paramètres de Pénalité et Préconditionneur    |
| - Si AAL / C-AAL : champ adaptatif r(t, x) selon temps physiques   |
| - Assemblage des blocs enrichis et complément de Schur             |
+--------------------------------------------------------------------+
                                   |
                                   v
+--------------------------------------------------------------------+
| Étape 4 : Boucle d'Uzawa / Résolution Navier-Stokes Couplée        |
| Pour chaque itération d'Uzawa m = 1 à uzawa_iter :                 |
|   a) Assemblage second membre (convection, diffusion, gravite, CSF)|
|   b) Résolution linéaire par Krylov (BiCGSTAB2 / MINRES / PCG) :   |
|        P^{-1} A x = P^{-1} b                                       |
|   c) Déballage des vitesses résolues u^{n+1}, v^{n+1}              |
|   d) Mise à jour Uzawa : p^{m+1} = p^m - r_eff * (div(u) + ...)    |
|   e) Test d'arrêt sur divergence : max|div(u)| <= tol -> SORTIE    |
+--------------------------------------------------------------------+
                                   |
                                   v
+--------------------------------------------------------------------+
| Étape 5 : Résolution de l'Équation Thermique (Si activée)          |
| - Transport advectif et diffusion de la température T              |
| - Calcul du terme de flottabilité de Boussinesq                    |
+--------------------------------------------------------------------+
                                   |
                                   v
+--------------------------------------------------------------------+
| Étape 6 : Sorties Numériques et Monitoring                         |
| - Écriture des fichiers XML VTK (.vtr/.pvtr)                       |
| - Enregistrement des résidus dans convergence.csv                  |
| - Avancement du temps physique : t = t + dt                        |
+--------------------------------------------------------------------+
```
