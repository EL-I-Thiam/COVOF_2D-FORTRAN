# Introduction Générale et Capacités du Solveur

## Présentation et Philosophie de COVOF 2D

Le code **COVOF 2D** (*Coupled Volume-Of-Fluid 2D*) est un environnement logiciel de recherche et d'ingénierie avancée dédié à la simulation numérique directe (DNS) et à la modélisation à haute fidélité des écoulements de fluides complexes, incompressibles ou faiblement compressibles, monophasiques ou diphasiques à surface libre. Issu des développements méthodologiques conduits dans le lignage des codes académiques de référence du laboratoire **I2M (UMR CNRS 5295, Université de Bordeaux)** — en particulier les plateformes **Thétis** et **Fugu** — et enrichi par les avancées récentes sur les méthodes monolithiques totalement couplées et les formulations adaptatives de Lagrangien augmenté, COVOF 2D constitue un banc d'essai et un outil de production scientifique de premier plan.

### Paradigme de Conception : L'Architecture Duale Jumelle C++17 / Fortran 2018

Contrairement aux logiciels CFD traditionnels reposant sur un seul langage ou intégrant du code hérité (*legacy code*) au sein d'un emballage hétérogène, COVOF 2D est structuré selon un principe novateur de **double implémentation miroir intégrale** :

1. **La branche `COVOF_2D-C++`** : Écrite en C++ moderne conforme au standard **ISO C++17**, elle privilégie une conception modulaire orientée objet avec gestion stricte des ressources par RAII (*Resource Acquisition Is Initialization*), polymorphisme dynamique pour les solveurs de Krylov et préconditionneurs, et encapsulation des structures de données topologiques et physiques.
2. **La branche `COVOF_2D-FORTRAN`** : Écrite en **Fortran 2018**, elle implémente rigoureusement les mêmes algorithmes au sein de modules fortement typés (`implicit none` systématique, précision double standardisée `real(c_double)` issue de `iso_c_binding`), exploitant la disposition mémoire native par colonnes (*column-major*) et la vectorisation optimale des boucles imbriquées.

Cette architecture duale n'a pas pour but de faire cohabiter des fragments de code non vérifiés, mais de garantir une **parité numérique stricte (bit-à-bit)** entre deux représentations logicielles indépendantes. Les deux exécutables :
- Partagent la même spécification du format d'entrée JSON (`template_parametres_covof2d.json`).
- Génèrent des résultats spatio-temporels identiques exportés au format standardisé XML VTK Rectilinear Grid (`.vtr` pour les sous-domaines locaux et `.pvtr` pour la vue globale assemblée sous ParaView).
- Utilisent une décomposition cartésienne 2D équivalente via la bibliothèque standard **MPI** (*Message Passing Interface*).
- Permettent aux chercheurs et aux ingénieurs de confronter les performances, l'empreinte mémoire et la portabilité des deux écosystèmes logiciels majeurs du calcul scientifique haute performance (HPC).

```
              +-------------------------------------+
              |    Fichier de Configuration JSON    |
              | (template_parametres_covof2d.json)  |
              +-------------------------------------+
                                 |
                 +---------------+---------------+
                 |                               |
                 v                               v
  +-------------------------------+ +-------------------------------+
  | Exécutable C++17              | | Exécutable Fortran 2018       |
  | bin/covof2d_cpp               | | bin/covof2d_fortran           |
  +-------------------------------+ +-------------------------------+
  | - Orchestration OO            | | - Orchestration Modulaire     |
  | - Grille MAC (Mesh.hpp)       | | - Grille MAC (module_mesh)    |
  | - Opérateur NS (NavierStokes) | | - Opérateur NS (module_ns)    |
  | - VOF-PLIC (InterfaceTracker) | | - VOF-PLIC (module_vof)       |
  | - Krylov & Schur (LinSolver)  | | - Krylov & Schur (module_solv)|
  | - MPI Cartésien 2D (MpiComm)  | | - MPI Cartésien (module_mpi)  |
  +-------------------------------+ +-------------------------------+
                 |                               |
                 +---------------+---------------+
                                 |
                                 v
              +-------------------------------------+
              | Sorties Numériques Unifiées         |
              | - Fichiers XML VTK (.vtr, .pvtr)    |
              | - Journaux de Convergence (.csv)    |
              | - Diagnostics Visuels ParaView      |
              +-------------------------------------+
```

### Domaines d'Application Visés

COVOF 2D est conçu pour traiter un spectre étendu de configurations de mécanique des fluides et de thermique appliquée :
- **Hydrodynamique à surface libre et écoulements diphasiques :** Rupture de barrage (*dam break*), ballottement de réservoirs (*sloshing*), impact et étalement de gouttes, dynamique capillaire et coalescence de bulles, instabilités interfaciales de Rayleigh-Taylor, Kelvin-Helmholtz et Plateau-Rayleigh.
- **Microfluidique et dynamique capillaire :** Écoulements dominés par les forces de tension superficielle, formation de jets, instabilité de rupture de filets fluides, oscillation de gouttes selon le mode propre de Rayleigh-Lamb.
- **Convection thermique et transferts couplés :** Convection naturelle dans des cavités fermées différentiellement chauffées, convection de Rayleigh-Bénard, panaches thermiques instationnaires sous approximation de Boussinesq.
- **Écoulements en milieux poreux et corps immergés :** Écoulements autour d'obstacles solides rigides ou poreux modélisés sans remaillage par pénalisation de Brinkman-Navier-Stokes (cylindre circulaire, plaque mince, marche descendante).
- **Acoustique numérique et régimes All-Mach :** Ondes de pression acoustiques à bas nombre de Mach, compression de bulles sous saut thermodynamique, régulation par méthode C-AAL barotrope.

---

## Matrice Complète des Capacités Physiques et Numériques

Le tableau suivant récapitule l'ensemble des modules physiques, algorithmes spatio-temporels, solveurs linéaires et modèles constitutifs implémentés et validés dans COVOF 2D.

| Catégorie | Fonctionnalité / Modèle | Options Disponibles | Description et Références |
| :--- | :--- | :--- | :--- |
| **Topologie & Maillage** | Type de grille | Grille Cartésienne Structurée Uniforme 2D | Maillage décalé MAC (*Marker-and-Cell*, Harlow & Welch 1965) |
| | Gestion des bords | Cellules fantômes (*Ghost cells*) | Couche de 1 à 2 mailles fantômes pour les CL d'ordre 2 |
| | Parallélisme spatial | Décomposition de domaine MPI 2D | Topologie cartésienne virtuelle `MPI_Cart_create` avec voisins E, W, N, S |
| **Régimes Physiques** | Incompressible | Monophasique & Diphasique | $\nabla \cdot \mathbf{u} = 0$, conservation exacte de la masse |
| | Faiblement compressible | Bas-Mach / All-Mach | Modèle barotrope $\chi_T \frac{\partial p}{\partial t} + \nabla \cdot \mathbf{u} = 0$ (El Ouafa et al. 2024) |
| | Thermique | Boussinesq / Rayleigh-Bénard | Terme de flottabilité $\rho \mathbf{g} \beta (T - T_{\text{ref}})$ |
| | Milieux poreux | Modèle de Darcy-Brinkman | Terme de traînée poreuse $-\frac{\mu}{K}\mathbf{u}$ pénalisé à $K \to 0$ pour corps solide |
| **Suivi d'Interface** | Modèle cinématique | VOF-PLIC Géométrique | Reconstruction linéaire par morceaux (*Piecewise Linear Interface Calculation*) |
| | Calcul des normales | Youngs (1982) / Elvira | Dérivées spatiales 2D de la fraction volumique $C$ |
| | Advection VOF | Weymouth & Yue (2010) | Schéma eulerien conservatif à séparation d'opérateurs 1D sans accumulation |
| | Tension superficielle | Modèle CSF (Brackbill et al. 1992) | Force volumique $\mathbf{f}_\sigma = \sigma \kappa \nabla C$, courbure $\kappa = -\nabla \cdot \mathbf{n}$ |
| **Couplage Vitesse-Pression** | Lagrangien Augmenté Standard | `AUGMENTED_LAGRANGIAN` | Pénalité $r \nabla(\nabla \cdot \mathbf{u})$ et boucle d'Uzawa $p^{k+1} = p^k - r \nabla \cdot \mathbf{u}$ |
| | Lagrangien Adaptatif Incompressible | `AAL` | Choix spatial adaptatif $r(t, M)$ basé sur les temps caractéristiques (Vincent 2004) |
| | Lagrangien Compressible Bas-Mach | `C-AAL` | Extension acoustique $r_{\text{eff}} = \frac{r}{1 + r \chi_T / \Delta t}$ (El Ouafa 2024) |
| | Résolution Monolithique | `MONOLITHIC` / `FULLY_COUPLED` | Système couplé point-selle global $\begin{pmatrix} A & B^T \\ B & -C_p \end{pmatrix} \begin{pmatrix} \mathbf{u} \\ p \end{pmatrix} = \begin{pmatrix} \mathbf{f} \\ g \end{pmatrix}$ |
| | Projection Standard | `STANDARD_PROJECTION` | Schéma à pas fractionnaire de Chorin-Temam (équation de Poisson en pression) |
| **Schémas Convectifs** | Advection de quantité de mouvement | `UPWIND` (Ordre 1) | Monotone, inconditionnellement stable mais diffusif |
| | | `CENTERED` (Ordre 2) | Précision spatiale $\mathcal{O}(\Delta x^2)$, non diffusif mais sujet aux oscillations de Peclet |
| | | `NEWTON` | Linéarisation complète de Newton pour convergence robuste en monolithique |
| | | `ADAMS_BASHFORTH` (AB2) | Extrapolation temporelle d'ordre 2 des vitesses convectives |
| **Intégration Temporelle** | Schémas temporels | `EULER_EXPLICIT` (Ordre 1) | Schéma basique pour l'analyse des petits pas de temps |
| | | `GEAR` / `BDF2` (Ordre 2) | Schéma à trois niveaux temporels implicite stable pour Navier-Stokes |
| | | `SSP_RK3` (Ordre 3) | Runge-Kutta à préservation de stabilité forte pour les fronts raides |
| **Solveurs Linéaires** | Algorithmes de Krylov Matrix-Free | `BICGSTAB2` | BiCGSTAB d'ordre 2 (Sleijpen & van der Vorst 1993) |
| | | `BICGSTAB` | BiCGSTAB standard (van der Vorst 1992) |
| | | `CG` / `CONJUGATE_GRADIENT` | Gradient Conjugué symétrique (Hestenes & Stiefel 1952) |
| | | `GMRES` | Generalized Minimal Residual redémarré (Saad & Schultz 1986) |
| | | `MINRES` | Minimum Residual pour opérateurs symétriques indéfinis (Paige & Saunders 1975) |
| | Solveur Direct Exact | `DIRECT` / `MUMPS` / `LAPACK` | Factorisation LU dense `dgetrf` / `dgetrs` avec mise en cache optimisée |
| **Préconditionneurs** | Méthodes parallèles | `JACOBI` / `DIAGONAL` | Inversion locale de la diagonale enrichie (visqueuse + AL + Brinkman + Newton) |
| | | `SSOR` / `GAUSS_SEIDEL` | Balayage avant/arrière symétrique avec facteur de relaxation $\omega = 1.0$ |
| | | `ILU0` | Factorisation incomplète Incomplete LU sans remplissage |
| | | `BLOCK_TRIANGULAR` / `BT_SBP` | Préconditionneur par blocs avec complément de Schur approché (El Ouafa 2023) |
| | | `CUSTOM` | Hook développeur permettant l'injection dynamique d'un opérateur externe |

---

## Architecture Fonctionnelle Vue par l'Utilisateur

L'utilisation de COVOF 2D s'articule autour d'un flux de travail unifié et standardisé en quatre étapes séquentielles :

```
+--------------------------------------------------------------------+
| 1. PRÉ-TRAITEMENT ET CONFIGURATION                                 |
| - Définition de la physique (viscosités, masses volumiques, sigma) |
| - Spécification de la grille MAC (nx, ny, dimensions lx, ly)       |
| - Définition des conditions initiales (interfaces VOF, vitesses)   |
| - Paramétrage des conditions aux limites (parois, entrées, sorties)|
| - Édition du fichier maître JSON : config.json                     |
+--------------------------------------------------------------------+
                                   |
                                   v
+--------------------------------------------------------------------+
| 2. EXÉCUTION DU CALCUL (Séquentiel ou Parallèle MPI)               |
| - Commande séquentielle : ./bin/covof2d_cpp config.json            |
| - Commande parallèle   : mpirun -np 4 ./bin/covof2d_fortran config |
| - Initialisation des structures et topologie cartésienne MPI       |
| - Boucle temporelle Navier-Stokes + VOF + Énergie + Uzawa          |
+--------------------------------------------------------------------+
                                   |
                                   v
+--------------------------------------------------------------------+
| 3. CONTRÔLE EN TEMPS RÉEL (Monitoring)                             |
| - Affichage console (pas de temps, temps physique, résidus L2)     |
| - Enregistrement de l'historique de convergence (convergence.csv)  |
| - Surveillance des critères de stabilité CFL convectif et diffusif |
+--------------------------------------------------------------------+
                                   |
                                   v
+--------------------------------------------------------------------+
| 4. POST-TRAITEMENT ET VALORISATION                                 |
| - Fichiers de champs volumiques XML VTK Rectilinear (.vtr, .pvtr)  |
| - Visualisation sous ParaView (vitesse, pression, fraction VOF)    |
| - Scripts d'analyse Python (profils de vitesse, bilans de masse)   |
+--------------------------------------------------------------------+
```

### Organisation des Fichiers et Conventions d'Exécution

À la racine de l'espace de travail, les répertoires et fichiers principaux sont structurés comme suit :

- `COVOF_2D-C++/` : Répertoire source du code C++17 (fichiers d'en-tête dans `include/covof/`, sources dans `src/`, exécutable dans `bin/covof2d_cpp`).
- `COVOF_2D-FORTRAN/` : Répertoire source du code Fortran 2018 (sources dans `src/`, exécutable dans `bin/covof2d_fortran`).
- `VALIDATION_TOTALE/` : Répertoire exhaustif regroupant les 24 cas-tests académiques de validation, avec leurs spécifications détaillées (`.md`), leurs fichiers JSON de configuration prêts à l'emploi et leurs scripts Python d'analyse comparative.
- `output/` : Répertoire généré automatiquement lors de l'exécution, recevant les fichiers `.vtr`, `.pvtr`, les fichiers `.pvd` de série temporelle et les fichiers `.csv` de convergence.
