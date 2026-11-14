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
# Installation, Compilation et Déploiement HPC

## Prérequis Logiciels et Matrice de Compatibilité

Le déploiement de **COVOF 2D** sur stations de travail, serveurs de calcul ou supercalculateurs massivement parallèles requiert un ensemble d'outils et de bibliothèques standards de l'écosystème scientifique. Les deux branches jumelles (C++ et Fortran) ont été conçues pour minimiser les dépendances externes tout en garantissant des performances optimales.

### Matrice de Compatibilité des Compilateurs

Les deux versions du solveur s'appuient sur les standards internationaux **ISO C++17** et **ISO Fortran 2018**. Le tableau suivant dresse la liste des suites de compilateurs testées et validées en production :

| Famille de Compilateurs | Version Minimale | Compilateur C++ | Compilateur Fortran | Drapeaux Recommandés |
| :--- | :--- | :--- | :--- | :--- |
| **GNU Compiler Collection (GCC)** | $\ge 9.3$ (recommandé $\ge 11.2$) | `g++` | `gfortran` | `-O3 -march=native -Wall -Wextra` |
| **LLVM / Clang + Flang** | $\ge 14.0$ | `clang++` | `flang` / `gfortran` | `-O3 -march=native -Wall` |
| **Intel oneAPI** | $\ge 2022.1$ | `icpx` | `ifx` / `ifort` | `-O3 -xHost -fp-model fast=2` |
| **Apple Clang (macOS)** | $\ge 14.0$ (Xcode CLT) | `clang++` | `gfortran` (Homebrew) | `-O3 -Wall -framework Accelerate` |

### Bibliothèques Mathématiques et de Communication Requises

1. **MPI (*Message Passing Interface*) :**
   - Implémentations compatibles : **OpenMPI** ($\ge 4.0$), **MPICH** ($\ge 3.4$), **Intel MPI** ($\ge 2021$), **MVAPICH2** ($\ge 2.3$).
   - Requis pour l'échange non bloquant des mailles fantômes (`MPI_Isend`, `MPI_Irecv`) et les réductions globales de convergence (`MPI_Allreduce`).
2. **Algèbre Linéaire Dense (BLAS / LAPACK) :**
   - Sous macOS : le framework natif vectorisé **Accelerate** (`-framework Accelerate`) est détecté et lié automatiquement.
   - Sous Linux / Clusters : **OpenBLAS**, **Intel MKL** (*Math Kernel Library*), ou les bibliothèques système de référence `liblapack-dev` et `libblas-dev`.
   - Requis pour le solveur direct exact par factorisation LU (`dgetrf`, `dgetrs`).
3. **Moteur de Compilation CMake :**
   - Version minimale requise : **CMake $\ge 3.16$** (pour C++) et **$\ge 3.18$** (pour Fortran 2018).

---

## Compilation Pas à Pas avec CMake

Les deux projets disposent chacun de leur fichier `CMakeLists.txt` autonome situé à la racine de leur répertoire respectif (`COVOF_2D-C++/CMakeLists.txt` et `COVOF_2D-FORTRAN/CMakeLists.txt`).

### Compilation de la Branche C++ (`COVOF_2D-C++`)

L'arborescence utilise le mode de compilation hors-source (*out-of-source build*) standard :

```bash
# 1. Se positionner dans le répertoire C++
cd /chemin/vers/PROJET_ANTIGRAVITY/COVOF_2D-C++

# 2. Créer et entrer dans le dossier de construction
mkdir -p build && cd build

# 3. Configurer avec CMake en mode Release (optimisations maximales)
cmake -DCMAKE_BUILD_TYPE=Release ..

# 4. Compiler en parallèle sur les cœurs disponibles
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)
```

À l'issue de la compilation, l'exécutable binaire est automatiquement placé dans `COVOF_2D-C++/bin/covof2d_cpp`, et un lien symbolique de commodité `COVOF_2D-C++/bin/covof2d` est créé.

#### Options CMake Avancées pour C++

| Option CMake | Valeur par Défaut | Description |
| :--- | :--- | :--- |
| `CMAKE_BUILD_TYPE` | `Release` | `Release` (active `-O3`), `Debug` (active `-g -O0 -fsanitize=address`), `RelWithDebInfo` |
| `CMAKE_CXX_STANDARD` | `17` | Standard de langage imposé (17 requis au minimum) |
| `CMAKE_EXPORT_COMPILE_COMMANDS`| `ON` | Génère `compile_commands.json` pour les analyseurs de code LSP / clangd |

### Compilation de la Branche Fortran (`COVOF_2D-FORTRAN`)

La procédure est strictement symétrique :

```bash
# 1. Se positionner dans le répertoire Fortran
cd /chemin/vers/PROJET_ANTIGRAVITY/COVOF_2D-FORTRAN

# 2. Créer et entrer dans le dossier de construction
mkdir -p build && cd build

# 3. Configurer avec CMake en mode Release
cmake -DCMAKE_BUILD_TYPE=Release ..

# 4. Compiler en parallèle
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)
```

Le fichier exécutable binaire produit est situé dans `COVOF_2D-FORTRAN/bin/covof2d_fortran`.

### Drapeaux d'Optimisation et Vectorisation SIMD

Pour exploiter les registres vectoriels étendus (AVX2, AVX-512) des processeurs modernes sur supercalculateurs, il est vivement recommandé d'ajouter les options d'architecture matérielle ciblée :

```bash
# Exemple pour nœuds de calcul Intel Xeon Skylake/Icelake :
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_FLAGS="-O3 -mavx512f -mavx512dq -mavx512cd -mavx512bw -mavx512vl" \
      -DCMAKE_Fortran_FLAGS="-O3 -mavx512f -mavx512dq -mavx512cd -mavx512bw -mavx512vl" ..

# Exemple générique auto-adaptatif pour le nœud courant :
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_FLAGS="-O3 -march=native -mtune=native" \
      -DCMAKE_Fortran_FLAGS="-O3 -march=native -mtune=native" ..
```

---

## Déploiement sur Supercalculateurs (Clusters HPC)

Dans un environnement de calcul partagé piloté par un gestionnaire de travaux (*workload manager*), COVOF 2D doit être exécuté en allouant précisément les cœurs de calcul, la mémoire et le réseau d'interconnexion rapide (InfiniBand / Omni-Path).

### Script de Soumission pour Gestionnaire Slurm

Voici un script de soumission complet et commenté pour **Slurm** (`run_covof_slurm.sh`), configuré pour une simulation parallèle distribuée sur 2 nœuds de calcul (64 cœurs par nœud, soit 128 cœurs au total) :

```bash
#!/bin/bash
#SBATCH --job-name=COVOF2D_PRODUCTION
#SBATCH --output=covof2d_%j.log
#SBATCH --error=covof2d_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --exclusive

# Chargement de l'environnement logiciel modulé (Modules Environment)
module purge
module load gcc/12.2.0
module load openmpi/4.1.4
module load openblas/0.3.21

# Affichage des informations sur les ressources allouées
echo "=========================================================="
echo "Démarrage du calcul COVOF 2D sur Slurm"
echo "Job ID              : $SLURM_JOB_ID"
echo "Nombre de nœuds     : $SLURM_NNODES"
echo "Nombre total de tâches MPI : $SLURM_NTASKS"
echo "Liste des nœuds     : $SLURM_JOB_NODELIST"
echo "=========================================================="

# Définition des variables d'environnement MPI et OpenMP
export OMP_NUM_THREADS=1
export UCX_TLS=rc,mm,sm
export OMPI_MCA_pml=ucx

# Répertoire de travail et fichier d'entrée
WORKDIR=$SLURM_SUBMIT_DIR
CASE_CONFIG="$WORKDIR/case_config.json"

cd $WORKDIR

# Exécution de la version C++ (ou remplacer par covof2d_fortran)
srun --mpi=pmix ./bin/covof2d_cpp $CASE_CONFIG

echo "=========================================================="
echo "Fin normale du calcul COVOF 2D - Statut: $?"
echo "=========================================================="
```

### Script de Soumission pour Gestionnaire PBS / Torque

Pour les centres de calcul s'appuyant sur **PBS Professional** ou **Torque**, le script équivalent `run_covof_pbs.sh` est structuré comme suit :

```bash
#!/bin/bash
#PBS -N COVOF2D_PBS
#PBS -l select=2:ncpus=64:mpiprocs=64
#PBS -l walltime=08:00:00
#PBS -q standard
#PBS -j oe
#PBS -o covof2d_pbs.log

cd $PBS_O_WORKDIR

# Chargement des modules Intel oneAPI
module purge
module load intel-oneapi-compilers/2023.1.0
module load intel-oneapi-mpi/2021.9.0
module load intel-oneapi-mkl/2023.1.0

# Contrôle de l'affinité des processus sur les sockets NUMA
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core
export I_MPI_FABRICS=shm:ofi

# Lancement MPI
mpirun -np 128 ./bin/covof2d_fortran case_config.json
```

### Gestion de l'Affinité Mémoire et Topologie NUMA

Pour les maillages fins ($N_x \times N_y \ge 1024^2$), la vitesse d'accès mémoire (*bandwidth-bound*) dicte la performance du produit matrice-vecteur creux et des balayages VOF. Il est fondamental d'éviter que les processus MPI ne migrent d'un cœur ou d'un socket à l'autre en cours de calcul :

```bash
# OpenMPI : Verrouillage strict processus par cœur CPU
mpirun -np 32 --bind-to core --map-by core ./bin/covof2d_cpp input.json

# Intel MPI : Placement NUMA optimisé
mpirun -np 32 -genv I_MPI_PIN_PROCESSOR_LIST=allcores ./bin/covof2d_fortran input.json
```

---

## Conteneurisation et Vérification de l'Installation

Pour garantir la reproductibilité absolue des résultats numériques sur n'importe quel système d'exploitation sans conflit de bibliothèques, COVOF 2D peut être encapsulé dans une image de conteneur standardisée.

### Définition Singularity / Apptainer (`Apptainer.def`)

Singularity est l'outil privilégié sur les supercalculateurs car il s'exécute sans privilèges administrateur (*rootless*) :

```dockerfile
Bootstrap: docker
From: ubuntu:22.04

%post
    apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        cmake \
        libopenmpi-dev \
        openmpi-bin \
        libopenblas-dev \
        liblapack-dev \
        git \
        ca-certificates
    apt-get clean && rm -rf /var/lib/apt/lists/*

%environment
    export LC_ALL=C
    export OMPI_MCA_btl_vader_single_copy_mechanism=none

%runscript
    exec "$@"
```

Pour compiler et exécuter l'image Singularity :

```bash
# 1. Construire l'image Singularity
apptainer build covof2d.sif Apptainer.def

# 2. Exécuter un calcul MPI à travers le conteneur
mpirun -np 8 apptainer exec covof2d.sif ./COVOF_2D-C++/bin/covof2d_cpp case.json
```

### Protocole de Vérification Immédiate Post-Installation

Afin de valider l'intégrité de la chaîne de compilation et du runtime MPI, il convient de lancer le cas de test de référence analytique (Cas 1 - Laplacian 2D) :

```bash
# Vérification rapide sur 1 processus
./COVOF_2D-C++/bin/covof2d_cpp VALIDATION_TOTALE/validation_cas_1_laplacian/cas1_laplacian.json

# Vérification sur 4 processus MPI
mpirun -np 4 ./COVOF_2D-FORTRAN/bin/covof2d_fortran VALIDATION_TOTALE/validation_cas_1_laplacian/cas1_laplacian.json
```

Si le message final `>>> SIMULATION COVOF 2D TERMINEE AVEC SUCCES <<<` s'affiche sans avertissement de fuite mémoire ou d'erreur MPI, le solveur est opérationnel pour la production scientifique.
# Formulation Mathématique et Modèles Physiques

## Système d'Équations Gouvernantes

Le solveur **COVOF 2D** résout les équations fondamentales de la mécanique des milieux continus régissant les écoulements de fluides visqueux, newtoniens, monophasiques ou diphasiques, sous diverses approximations physiques (incompressible, faiblement compressible à bas nombre de Mach, avec ou sans transfert thermique, et en présence d'obstacles solides ou poreux).

### Équations de Navier-Stokes sous Forme Continue

Dans un domaine spatial ouvert et borné $\Omega \subset \mathbb{R}^2$ de frontière $\partial \Omega$ et pour un intervalle de temps $t \in [0, t_{\text{end}}]$, les équations de conservation de la masse et de la quantité de mouvement s'écrivent :

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0$$

$$\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) = -\nabla p + \nabla \cdot \boldsymbol{\tau} + \rho \mathbf{g} + \mathbf{f}_\sigma + \mathbf{f}_{\text{brinkman}} + \mathbf{f}_{\text{boussinesq}}$$

où :
- $\rho = \rho(\mathbf{x}, t)$ est la masse volumique locale du milieu $[\text{kg}\cdot\text{m}^{-3}]$.
- $\mathbf{u} = (u, v)^T$ est le vecteur vitesse eulérien $[\text{m}\cdot\text{s}^{-1}]$.
- $p = p(\mathbf{x}, t)$ est le champ de pression hydrodynamique $[\text{Pa}]$.
- $\boldsymbol{\tau}$ est le tenseur des contraintes visqueuses de déviation $[\text{Pa}]$.
- $\mathbf{g} = (g_x, g_y)^T$ est le vecteur accélération de la pesanteur $[\text{m}\cdot\text{s}^{-2}]$.
- $\mathbf{f}_\sigma$ est la force volumique de tension superficielle concentrée à l'interface $[\text{N}\cdot\text{m}^{-3}]$.
- $\mathbf{f}_{\text{brinkman}}$ est la force de traînée volumique de pénalisation poreuse ou solide de Brinkman $[\text{N}\cdot\text{m}^{-3}]$.
- $\mathbf{f}_{\text{boussinesq}}$ est la force volumique de flottabilité thermique sous approximation de Boussinesq $[\text{N}\cdot\text{m}^{-3}]$.

### Formulation Advective Non Conservative

En exploitant l'équation de continuité, la dérivée matérielle de la vitesse s'exprime sous forme advective (forme utilisée pour la discrétisation par volumes finis décalés dans COVOF 2D) :

$$\rho \left( \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} \right) = -\nabla p + \nabla \cdot \boldsymbol{\tau} + \rho \mathbf{g} + \mathbf{f}_\sigma + \mathbf{f}_{\text{brinkman}} + \mathbf{f}_{\text{boussinesq}}$$

Soit, projetée selon les composantes cartésiennes horizontales $x$ et verticales $y$ :

$$\rho \left( \frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} \right) = -\frac{\partial p}{\partial x} + \frac{\partial \tau_{xx}}{\partial x} + \frac{\partial \tau_{xy}}{\partial y} + \rho g_x + f_{\sigma, x} + f_{\text{brinkman}, x}$$

$$\rho \left( \frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} \right) = -\frac{\partial p}{\partial y} + \frac{\partial \tau_{yx}}{\partial x} + \frac{\partial \tau_{yy}}{\partial y} + \rho g_y + f_{\sigma, y} + f_{\text{brinkman}, y} + f_{\text{boussinesq}, y}$$

---

## Modèle 1-Fluide Diphasique et Suivi d'Interface VOF

Dans la formulation 1-fluide de Kataoka (1986), le milieu diphasique est modélisé comme un fluide unique continu équivalent dont les propriétés thermodynamiques et de transport varient de manière discontinue à la traversée de l'interface séparant la phase 1 (fluide porteur ou liquide) de la phase 2 (goutte, bulle ou gaz).

### Fonction Indicatrice de Couleur (Fraction Volumique VOF)

La distribution spatiale des deux phases est représentée par la fonction indicatrice de couleur $C(\mathbf{x}, t) \in [0, 1]$ :

$$C(\mathbf{x}, t) = \begin{cases} 1 & \text{si } \mathbf{x} \in \text{Phase 1} \\ 0 & \text{si } \mathbf{x} \in \text{Phase 2} \\ ]0, 1[ & \text{dans les mailles traversées par l'interface} \end{cases}$$

L'équation cinématique de transport de l'interface sans changement de phase s'énonce :

$$\frac{\partial C}{\partial t} + \mathbf{u} \cdot \nabla C = 0$$

Sous forme conservative pour un écoulement à divergence nulle :

$$\frac{\partial C}{\partial t} + \nabla \cdot (C \mathbf{u}) = 0$$

### Lois de Mélange des Propriétés Matérielles

Dans chaque cellule de calcul, la masse volumique effective $\rho$ et la viscosité dynamique $\mu$ sont reconstruites à partir de la fraction volumique locale $C$ :

1. **Masse volumique cellulaire (Moyenne arithmétique) :**
   $$\rho(C) = C \rho_1 + (1 - C) \rho_2$$
   où $\rho_1$ est la masse volumique de la phase 1 (`cfg%rho1` / `config.rho1`) et $\rho_2$ celle de la phase 2 (`cfg%rho2` / `config.rho2`).

2. **Viscosité dynamique cellulaire et nodale :**
   Dans le centre de la maille :
   $$\mu(C) = C \mu_1 + (1 - C) \mu_2$$
   Pour le calcul des termes de cisaillement aux sommets du maillage décalé MAC ($\mu_{xy}$ aux nœuds $(i+\frac{1}{2}, j+\frac{1}{2})$), COVOF 2D utilise la **moyenne harmonique** des quatre mailles adjacentes, qui est la seule formulation rigoureusement conservative et stable pour les rapports de viscosité élevés ($\mu_1/\mu_2 \gg 1$) :
   $$\mu_{i+1/2, j+1/2} = \left[ \frac{1}{4} \left( \frac{1}{\mu_{i, j}} + \frac{1}{\mu_{i+1, j}} + \frac{1}{\mu_{i, j+1}} + \frac{1}{\mu_{i+1, j+1}} \right) \right]^{-1}$$

### Modèle de Tension Superficielle CSF (Continuum Surface Force)

La tension superficielle agissant sur l'interface est convertie en un terme source volumique équivalent selon la méthode CSF de Brackbill, Kothe et Zemach (1992) :

$$\mathbf{f}_\sigma = \sigma \kappa \nabla C$$

où :
- $\sigma$ est le coefficient de tension superficielle interfaciale $[\text{N}\cdot\text{m}^{-1}]$ (`cfg%sigma` / `config.sigma`).
- $\kappa$ est la courbure moyenne locale de l'interface $[\text{m}^{-1}]$.
- $\nabla C$ est le gradient spatial de la fraction volumique, orienté selon la normale unitaire à l'interface.

La normale locale unitaire $\mathbf{n}$ et la courbure $\kappa$ sont calculées par dérivation numérique du champ de fraction volumique régularisé :

$$\mathbf{n} = \frac{\nabla C}{\|\nabla C\| + \epsilon}, \quad \kappa = -\nabla \cdot \mathbf{n} = -\left( \frac{\partial n_x}{\partial x} + \frac{\partial n_y}{\partial y} \right)$$

où $\epsilon = 10^{-30}$ est un seuil de régularisation infinitésimal évitant les divisions par zéro dans les phases pures uniformes.

---

## Modèle Faiblement Compressible et Régime All-Mach

Pour les écoulements diphasiques présentant une compressibilité acoustique ou thermodynamique (notamment la compression d'inclusions gazeuses ou la propagation d'ondes de pression à bas nombre de Mach), COVOF 2D intègre la formulation **All-Mach C-AAL** (El Ouafa et al., *Fluids* 2024).

### Équation de Continuité Généralisée

L'équation de conservation de la masse est réécrite en introduisant le coefficient de compressibilité isotherme $\chi_T$ :

$$\chi_T = \frac{1}{\rho}\left( \frac{\partial \rho}{\partial p} \right)_T = \frac{1}{\rho c_s^2}$$

où $c_s$ désigne la célérité du son dans le fluide considéré $[\text{m}\cdot\text{s}^{-1}]$. L'équation de continuité prend alors la forme :

$$\chi_T \left( \frac{\partial p}{\partial t} + \mathbf{u} \cdot \nabla p \right) + \nabla \cdot \mathbf{u} = 0$$

À bas nombre de Mach ($Ma = \|\mathbf{u}\|/c_s \ll 1$), le terme convectif $\mathbf{u} \cdot \nabla p$ est d'ordre $\mathcal{O}(Ma^2)$ devant la dérivée instationnaire et est négligé, ce qui conduit à la relation de contrainte barotrope discrète :

$$\chi_T \frac{\partial p}{\partial t} + \nabla \cdot \mathbf{u} = 0 \quad \iff \quad \nabla \cdot \mathbf{u}^{n+1} + \frac{\chi_T}{\Delta t} (p^{n+1} - p^n) = 0$$

Cette relation modifie directement l'étape de mise à jour d'Uzawa dans la méthode du Lagrangien augmenté adaptatif compressible (C-AAL), régularisant le paramètre de pénalité local :

$$r_{\text{eff}}(\mathbf{x}) = \frac{r_{\text{AAL}}(\mathbf{x})}{1 + r_{\text{AAL}}(\mathbf{x}) \frac{\chi_T(\mathbf{x})}{\Delta t}}$$

---

## Modèle de Milieu Poreux et Pénalisation de Brinkman

Pour représenter des obstacles rigides immergés (pylônes, cylindres, marches, vannes) sans faire appel à un remaillage non structuré conforme, COVOF 2D implémente la méthode de pénalisation volumique de Brinkman (Angot, Bruneau & Fabrie, 1999).

### Équation de Brinkman-Navier-Stokes

Un terme de traînée de Darcy-Brinkman est ajouté au membre de droite de l'équation de quantité de mouvement :

$$\mathbf{f}_{\text{brinkman}} = -\frac{\mu}{\mathcal{K}(\mathbf{x})} \mathbf{u}$$

où $\mathcal{K}(\mathbf{x})$ est le champ spatial de perméabilité intrinsèque de Darcy $[\text{m}^2]$ :
- **Dans le fluide libre ($\mathbf{x} \in \Omega_f$) :** $\mathcal{K} \to +\infty \implies \mathbf{f}_{\text{brinkman}} = \mathbf{0}$, les équations de Navier-Stokes standard sont recouvrées exactement.
- **Dans un milieu poreux perméable ($\mathbf{x} \in \Omega_p$) :** $\mathcal{K} = K_{\text{poreux}} \in [10^{-6}, 10^{-2}]\,\text{m}^2$, la loi de Darcy avec cisaillement de Brinkman régit la filtration.
- **Dans un obstacle solide rigide imperméable ($\mathbf{x} \in \Omega_s$) :** $\mathcal{K} = \eta \to 0$ ($\eta \le 10^{-12}\,\text{m}^2$, paramètre `obs_perm` / `PERMEABILITY_SOLID_DEFAULT`). Par pénalisation asymptotique, l'équation impose mathématiquement :
  $$\mathbf{u} = \mathbf{0} + \mathcal{O}(\sqrt{\eta}) \quad \text{dans } \Omega_s$$
  garantissant le respect de la condition d'adhérence sans coût géométrique supplémentaire.

---

## Transfert Thermique et Approximation de Boussinesq

Pour les écoulements thermoconvectifs (cavité différentiellement chauffée, instabilités de Rayleigh-Bénard), le solveur couple l'équation de quantité de mouvement à une équation de transport-diffusion de l'énergie thermique.

### Équation de l'Énergie en Température

En négligeant la dissipation visqueuse et les variations de pression :

$$\rho C_p \left( \frac{\partial T}{\partial t} + \mathbf{u} \cdot \nabla T \right) = \nabla \cdot (k_{\text{th}} \nabla T) + S_T$$

où :
- $T(\mathbf{x}, t)$ est le champ de température absolue $[\text{K}]$ (`fields.T` / `fields%T`).
- $C_p$ est la capacité thermique isobare $[\text{J}\cdot\text{kg}^{-1}\cdot\text{K}^{-1}]$ (`cfg%cp1` / `config.cp1`).
- $k_{\text{th}}$ est la conductivité thermique du fluide $[\text{W}\cdot\text{m}^{-1}\cdot\text{K}^{-1}]$ (`cfg%th_cond1` / `config.th_cond1`).
- $S_T$ est un terme source volumique de chaleur $[\text{W}\cdot\text{m}^{-3}]$.

### Couplage de Flottabilité de Boussinesq

La variation de la masse volumique avec la température est linéarisée autour d'un état de référence $(T_{\text{ref}}, \rho_0)$ :

$$\rho(T) \approx \rho_0 \left( 1 - \beta (T - T_{\text{ref}}) \right)$$

où $\beta$ est le coefficient de dilatation volumique thermique $[\text{K}^{-1}]$ (`cfg%beta_th` / `config.beta_th`). Le terme source de gravité prend alors la forme :

$$\mathbf{f}_{\text{boussinesq}} = -\rho_0 \beta (T - T_{\text{ref}}) \mathbf{g}$$

Ce terme est orienté verticalement selon l'axe de la pesanteur, induisant la mise en mouvement naturelle du fluide chaud plus léger vers le haut et du fluide froid plus dense vers le bas.

---

## Nombres Adimensionnels et Régimes d'Écoulement

La physique des écoulements simulés dans COVOF 2D est caractérisée par les nombres sans dimension classiques, définis à partir d'une longueur de référence $L$, d'une vitesse de référence $U$, d'un écart de température $\Delta T = T_{\text{hot}} - T_{\text{cold}}$ et de la viscosité cinématique $\nu = \mu/\rho$ :

| Nombre Sans Dimension | Expression Mathématique | Signification Physique | Modules COVOF 2D Concernés |
| :--- | :--- | :--- | :--- |
| **Reynolds ($Re$)** | $Re = \frac{\rho U L}{\mu} = \frac{U L}{\nu}$ | Rapport forces d'inertie / forces visqueuses | Tous (Cavité, Poiseuille, Cylindre) |
| **Froude ($Fr$)** | $Fr = \frac{U}{\sqrt{g L}}$ | Rapport forces d'inertie / forces gravitationnelles | Surface libre, Rupture de barrage |
| **Weber ($We$)** | $We = \frac{\rho U^2 L}{\sigma}$ | Rapport forces d'inertie / forces capillaires | Rupture de nappe, Goutte, Jet |
| **Capillaire ($Ca$)** | $Ca = \frac{\mu U}{\sigma} = \frac{We}{Re}$ | Rapport forces visqueuses / forces capillaires | Microfluidique, Instabilité de jet |
| **Bond / Eötvös ($Bo$)** | $Bo = \frac{\Delta \rho g L^2}{\sigma}$ | Rapport forces gravitationnelles / capillaires | Ascension de bulle, Sédimentation |
| **Ohnesorge ($Oh$)** | $Oh = \frac{\mu}{\sqrt{\rho \sigma L}} = \frac{\sqrt{We}}{Re}$ | Dissipation visqueuse relative aux ondes capillaires | Oscillation de bulle, Atomisation |
| **Prandtl ($Pr$)** | $Pr = \frac{\nu}{\alpha} = \frac{\mu C_p}{k_{\text{th}}}$ | Rapport diffusivité visqueuse / diffusivité thermique | Rayleigh-Bénard, Convection naturelle |
| **Rayleigh ($Ra$)** | $Ra = \frac{g \beta \Delta T L^3}{\nu \alpha} = Gr \cdot Pr$ | Intensité motrice de la convection thermique | Cavité différentiellement chauffée |
| **Darcy ($Da$)** | $Da = \frac{K}{L^2}$ | Perméabilité adimensionnelle du milieu poreux | Milieux poreux, Pénalisation Brinkman |

---

## Table Générale de Nomenclature

Le tableau ci-dessous recense l'intégralité des grandeurs scalaires, vectorielles et tensorielles utilisées dans la formulation continue et leur dénomination exacte dans les codes sources C++ et Fortran :

| Symbole | Signification Physique | Unités SI | Variable C++ (`COVOF_2D-C++`) | Variable Fortran (`COVOF_2D-FORTRAN`) |
| :--- | :--- | :--- | :--- | :--- |
| $\mathbf{u}$ | Vecteur vitesse eulérien | $\text{m}\cdot\text{s}^{-1}$ | `fields.u`, `fields.v` | `fields%u`, `fields%v` |
| $p$ | Champ de pression scalaire | $\text{Pa} = \text{N}\cdot\text{m}^{-2}$ | `fields.p` | `fields%p` |
| $p^n$ | Pression au pas de temps antérieur | $\text{Pa}$ | `fields.p_old` | `fields%p_old` |
| $\rho$ | Masse volumique du fluide | $\text{kg}\cdot\text{m}^{-3}$ | `fields.rho` | `fields%rho` |
| $\mu$ | Viscosité dynamique | $\text{Pa}\cdot\text{s} = \text{kg}\cdot\text{m}^{-1}\cdot\text{s}^{-1}$ | `fields.mu` | `fields%mu` |
| $C$ | Fraction volumique de couleur (VOF) | Sans dimension | `fields.c` | `fields%c` |
| $\mathbf{n}$ | Vecteur normal unitaire à l'interface | Sans dimension | `fields.nx`, `fields.ny` | `fields%nx`, `fields%ny` |
| $\kappa$ | Courbure locale de l'interface | $\text{m}^{-1}$ | `fields.kappa` | `fields%kappa` |
| $\sigma$ | Tension superficielle | $\text{N}\cdot\text{m}^{-1}$ | `config.sigma` | `cfg%sigma` |
| $\mathbf{f}_\sigma$ | Force de tension superficielle CSF | $\text{N}\cdot\text{m}^{-3}$ | `fields.fx_surf`, `fields.fy_surf` | `fields%fx_surf`, `fields%fy_surf` |
| $\chi_T$ | Compressibilité isotherme | $\text{Pa}^{-1} = \text{m}\cdot\text{s}^2\cdot\text{kg}^{-1}$ | `fields.chi_T` | `fields%chi_T` |
| $r_{\text{al}}$ | Paramètre de pénalité Lagrangien | $\text{Pa}\cdot\text{s}\cdot\text{m}^{-1}$ | `config.r_al`, `fields.r_aal` | `cfg%r_al`, `fields%r_aal` |
| $K$ | Perméabilité de Darcy | $\text{m}^2$ | `config.obs_perm` | `cfg%obs_perm` |
| $T$ | Température absolue | $\text{K}$ | `fields.T` | `fields%T` |
| $k_{\text{th}}$ | Conductivité thermique | $\text{W}\cdot\text{m}^{-1}\cdot\text{K}^{-1}$ | `fields.k_therm` | `fields%k_therm` |
| $C_p$ | Capacité calorifique isobare | $\text{J}\cdot\text{kg}^{-1}\cdot\text{K}^{-1}$ | `fields.cp_therm` | `fields%cp_therm` |
| $\beta$ | Dilatation thermique volumique | $\text{K}^{-1}$ | `config.beta_th` | `cfg%beta_th` |
| $\mathbf{g}$ | Accélération de la pesanteur | $\text{m}\cdot\text{s}^{-2}$ | `config.gx`, `config.gy` | `cfg%gx`, `cfg%gy` |
| $\Delta t$ | Pas de temps d'intégration | $\text{s}$ | `dt` | `dt` |
| $\Delta x, \Delta y$ | Pas d'espace de la grille MAC | $\text{m}$ | `mesh.dx`, `mesh.dy` | `mesh%dx`, `mesh%dy` |
| $N_x, N_y$ | Nombre de mailles de calcul locales | Entier | `mesh.nx_local`, `mesh.ny_local` | `mesh%nx_local`, `mesh%ny_local` |
# Méthodes Numériques

## Discrétisation Spatiale sur Grille Cartésienne Décalée (MAC)

Le solveur **COVOF 2D** utilise la discrétisation spatiale sur grille décalée dite *Marker-and-Cell* (MAC), introduite par Harlow et Welch (1965). Ce choix topologique fondamental est universellement reconnu en mécanique des fluides incompressibles pour ses propriétés mathématiques intrinsèques d'inf-sup stabilité (condition de Ladyzhenskaya-Babuška-Brezzi ou LBB discrète) et l'élimination naturelle des oscillations parasites en damier (*checkerboard oscillations*) de pression, sans requérir de dissipation artificielle ni de correction de Rhie-Chow.

### Positionnement Staggered des Variables Primordiales

Le domaine de calcul bidimensionnel est pavé par un maillage régulier de pas $\Delta x$ et $\Delta y$. Les différentes grandeurs physiques sont localisées comme suit :

1. **Centres des cellules $(i, j)$ :**
   - Pression hydrodynamique $p_{i, j}$
   - Fraction volumique VOF $C_{i, j}$
   - Masse volumique $\rho_{i, j}$ et viscosité dynamique cellulaire $\mu_{i, j}$
   - Température scalaire $T_{i, j}$
   - Paramètre adaptatif de pénalité de Lagrangien augmenté $r_{i, j}$
   - Compressibilité isotherme $\chi_{T, i, j}$

2. **Faces verticales (Normales à $x$) $(i + \frac{1}{2}, j)$ :**
   - Composante horizontale de la vitesse $u_{i+1/2, j}$
   - Masse volumique transportée aux faces $\rho_{u, i+1/2, j}$
   - Force volumique superficielle horizontale $f_{\sigma, x, i+1/2, j}$

3. **Faces horizontales (Normales à $y$) $(i, j + \frac{1}{2})$ :**
   - Composante verticale de la vitesse $v_{i, j+1/2}$
   - Masse volumique transportée aux faces $\rho_{v, i, j+1/2}$
   - Force volumique superficielle verticale $f_{\sigma, y, i, j+1/2}$
   - Terme de flottabilité thermique de Boussinesq $f_{\text{boussinesq}, y, i, j+1/2}$

4. **Sommets / Nœuds de coin $(i + \frac{1}{2}, j + \frac{1}{2})$ :**
   - Viscosité dynamique harmonique nodale $\mu_{i+1/2, j+1/2}$
   - Taux de déformation croisé $\left( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right)$
   - Contrainte de cisaillement visqueuse $\tau_{xy, i+1/2, j+1/2}$

```
               (i, j+1)
                  |
        +---------v---------+
        |                   |
        |      p(i, j)      |
   u    |      C(i, j)      |    u
(i-1/2,j)----> T(i, j)     ----> (i+1/2, j)
        |     rho(i, j)     |
        |                   |
        +---------^---------+
                  |
                (i, j)
              v(i, j-1/2)
```

### Gradients Discrets et Conservation de la Divergence

Sur la grille MAC, l'opérateur de divergence discrète de la vitesse $(\nabla \cdot \mathbf{u})_{i, j}$ au centre de la cellule $(i, j)$ s'écrit de manière compacte et exacte au second ordre :

$$(\nabla \cdot \mathbf{u})_{i, j} = \frac{u_{i+1/2, j} - u_{i-1/2, j}}{\Delta x} + \frac{v_{i, j+1/2} - v_{i, j-1/2}}{\Delta y}$$

Réciproquement, les composantes du gradient de pression agissant respectivement sur les équations de quantité de mouvement en $u$ et en $v$ s'évaluent sans aucune interpolation :

$$\left( \frac{\partial p}{\partial x} \right)_{i+1/2, j} = \frac{p_{i+1, j} - p_{i, j}}{\Delta x}, \quad \left( \frac{\partial p}{\partial y} \right)_{i, j+1/2} = \frac{p_{i, j+1} - p_{i, j}}{\Delta y}$$

Cette dualité adjointe parfaite assure que le gradient discret $B^T$ est l'opposé exact de l'adjoint de l'opérateur de divergence discrète $B$, conférant au système de Stokes une structure point-selle symétrique parfaite.

---

## Schémas d'Advection et Linéarisations Convectives

Le terme d'advection non linéaire $(\mathbf{u} \cdot \nabla) \mathbf{u}$ dans les équations de quantité de mouvement représente le couplage cinématique d'inertie. COVOF 2D met à disposition quatre formulations sélectionnables par la clé `convective_linearization` dans le fichier JSON.

### Schéma Décentré Amont d'Ordre 1 (`UPWIND`)

Le schéma décentré amont utilise le sens de propagation de l'écoulement pour évaluer la vitesse sur les faces de contrôle :

$$\left( u \frac{\partial u}{\partial x} \right)_{i+1/2, j} \approx \begin{cases} u_{i+1/2, j} \frac{u_{i+1/2, j} - u_{i-1/2, j}}{\Delta x} & \text{si } u_{i+1/2, j} \ge 0 \\ u_{i+1/2, j} \frac{u_{i+3/2, j} - u_{i+1/2, j}}{\Delta x} & \text{si } u_{i+1/2, j} < 0 \end{cases}$$

- **Avantages :** Stabilité inconditionnelle, positivité stricte, absence d'oscillations numériques parasites même à très haut nombre de Reynolds de maille $Re_{\Delta} = \frac{\rho U \Delta x}{\mu} > 2$.
- **Inconvénients :** Diffusion numérique artificielle d'ordre 1 proportionnelle à $\frac{u \Delta x}{2}$, tronquant les tourbillons fins.

### Schéma Centré d'Ordre 2 (`CENTERED`)

Le schéma centré calcule la moyenne arithmétique de la vitesse aux interfaces de chaque volume de contrôle :

$$\left( u \frac{\partial u}{\partial x} \right)_{i+1/2, j} \approx \frac{1}{\Delta x} \left[ \left(\frac{u_{i+1/2, j} + u_{i+3/2, j}}{2}\right)^2 - \left(\frac{u_{i-1/2, j} + u_{i+1/2, j}}{2}\right)^2 \right]$$

- **Avantages :** Précision spatiale rigoureusement d'ordre 2 $\mathcal{O}(\Delta x^2, \Delta y^2)$, aucune diffusion numérique artificielle, conservation stricte de l'énergie cinétique globale.
- **Précautions :** Peut générer des oscillations de dispersion lorsque $Re_{\Delta} > 2$, nécessitant un maillage suffisamment résolu.

### Linéarisation Complète de Newton (`NEWTON`)

Pour les simulations monolithiques stationnaires ou à grands pas de temps, COVOF 2D dispose d'une linéarisation complète de Newton-Raphson de l'opérateur convectif :

$$\mathbf{u}^{n+1} = \mathbf{u}^* + \delta \mathbf{u}$$

$$(\mathbf{u}^{n+1} \cdot \nabla) \mathbf{u}^{n+1} \approx (\mathbf{u}^* \cdot \nabla) \delta \mathbf{u} + (\delta \mathbf{u} \cdot \nabla) \mathbf{u}^* + (\mathbf{u}^* \cdot \nabla) \mathbf{u}^*$$

Cette formulation introduit les dérivées croisées $\frac{\partial u}{\partial x}$, $\frac{\partial u}{\partial y}$, $\frac{\partial v}{\partial x}$, $\frac{\partial v}{\partial y}$ dans la matrice jacobienne du système linéaire, garantissant un taux de convergence asymptotique **quadratique** indispensable pour la convergence rapide des solveurs de Krylov couplés.

### Schéma d'Adams-Bashforth d'Ordre 2 (`ADAMS_BASHFORTH` / `AB2`)

Dans ce schéma, la vitesse advectante est extrapolée explicitement à partir des deux pas de temps antérieurs $(n)$ et $(n-1)$ :

$$\mathbf{u}_{\text{adv}} = \frac{3}{2}\mathbf{u}^n - \frac{1}{2}\mathbf{u}^{n-1}$$

L'équation en vitesse devient alors un problème linéaire à coefficients constants durant le pas de temps, découplant l'advection du problème de diffusion-continuité.

---

## Discrétisation Temporelle et Stabilité

Le choix de l'algorithme d'avancement en temps est piloté par la clé `time.scheme` dans le fichier de configuration.

### Schémas Temporels Implémentés

1. **Euler Explicite d'Ordre 1 (`EULER_EXPLICIT`) :**
   $$\frac{\mathbf{u}^{n+1} - \mathbf{u}^n}{\Delta t} = \mathbf{F}(\mathbf{u}^n)$$
   Utilisé principalement pour la vérification numérique ou les régimes purement dominés par la convection avec CFL strict.

2. **Schéma de Gear / BDF2 d'Ordre 2 (`GEAR` / `BDF2`) :**
   Schéma à trois pas de temps implicite, inconditionnellement stable pour les termes linéaires :
   $$\frac{3 \mathbf{u}^{n+1} - 4 \mathbf{u}^n + \mathbf{u}^{n-1}}{2 \Delta t} = \mathbf{F}(\mathbf{u}^{n+1})$$
   Dans COVOF 2D, le coefficient $\alpha_{\text{time}}$ est automatiquement basculé à $1.5$ dès que l'historique temporel $\mathbf{u}^{n-1}$ (`fields.has_nm1`) est disponible, avec repli sur Euler implicite ($\alpha_{\text{time}} = 1.0$) lors du tout premier pas de temps.

3. **Runge-Kutta à Préservation de Stabilité Forte (`SSP_RK3`) :**
   Schéma d'ordre 3 en trois sous-étapes pour l'advection conservative conjointe de la masse et de la quantité de mouvement (modèle Consistent Mass-Momentum d'El Ouafa et al., *Fluids* 2021) :
   $$\mathbf{u}^{(1)} = \mathbf{u}^n + \Delta t \mathcal{L}(\mathbf{u}^n)$$
   $$\mathbf{u}^{(2)} = \frac{3}{4}\mathbf{u}^n + \frac{1}{4}\mathbf{u}^{(1)} + \frac{1}{4}\Delta t \mathcal{L}(\mathbf{u}^{(1)})$$
   $$\mathbf{u}^{n+1} = \frac{1}{3}\mathbf{u}^n + \frac{2}{3}\mathbf{u}^{(2)} + \frac{2}{3}\Delta t \mathcal{L}(\mathbf{u}^{(2)})$$

### Critères de Stabilité et Pas de Temps Adaptatif (CFL)

Lorsque `time.adaptative_dt` est activé, COVOF 2D recalcule à chaque itération le pas de temps maximal admissible $\Delta t$ sous la triple contrainte de Courant-Friedrichs-Lewy (CFL) :

1. **Condition CFL Convective :**
   $$\Delta t_{\text{conv}} = \text{CFL} \cdot \left[ \max_{i, j} \left( \frac{|u_{i+1/2, j}|}{\Delta x} + \frac{|v_{i, j+1/2}|}{\Delta y} \right) \right]^{-1}$$

2. **Condition CFL Diffusive Visqueuse (pour intégration explicite) :**
   $$\Delta t_{\text{diff}} = \frac{1}{2} \min_{i, j} \left( \frac{\rho_{i, j}}{\mu_{i, j}} \right) \left[ \frac{1}{\Delta x^2} + \frac{1}{\Delta y^2} \right]^{-1}$$

3. **Condition Capillaire d'Ondes de Surface (Brackbill et al. 1992) :**
   $$\Delta t_{\text{cap}} = \sqrt{\frac{(\rho_1 + \rho_2) \min(\Delta x, \Delta y)^3}{4 \pi \sigma}}$$

Le pas de temps effectif retenu est :
$$\Delta t^{n+1} = \min \left( \Delta t_{\text{max}}, \Delta t_{\text{conv}}, \Delta t_{\text{diff}}, \Delta t_{\text{cap}} \right)$$

---

## Formulations du Couplage Vitesse-Pression

Le traitement de la contrainte d'incompressibilité $\nabla \cdot \mathbf{u} = 0$ constitue le cœur algorithmique de tout solveur Navier-Stokes. COVOF 2D intègre cinq formulations majeures.

### 1. Lagrangien Augmenté Standard (`AUGMENTED_LAGRANGIAN`)

La méthode du Lagrangien Augmenté (Fortin & Glowinski, 1982 ; Caltagirone & Breil, 1999) pénalise l'équation de quantité de mouvement par la divergence de la vitesse tout en conservant le gradient de pression :

$$\rho \left( \frac{\alpha_{\text{time}} \mathbf{u}^{k+1}}{\Delta t} + (\mathbf{u}^n \cdot \nabla)\mathbf{u}^{k+1} \right) - \nabla \cdot \boldsymbol{\tau}^{k+1} - r \nabla(\nabla \cdot \mathbf{u}^{k+1}) = -\nabla p^k + \mathbf{f}$$

où $r$ est un paramètre scalaire uniforme de pénalité de Lagrangien augmenté ($r \sim 10^2 - 10^4\,\text{Pa}\cdot\text{s}\cdot\text{m}^{-1}$, configuré via `solver.r_al`). La pression est ensuite réactualisée par une étape d'Uzawa :

$$p^{k+1} = p^k - r \nabla \cdot \mathbf{u}^{k+1}$$

Si `solver.uzawa_iter = 1`, la méthode se comporte comme un pseudo-pas fractionnaire très rapide. Si `solver.uzawa_iter > 1`, une boucle itérative d'Uzawa est exécutée jusqu'à ce que $\|\nabla \cdot \mathbf{u}^{k+1}\|_\infty \le \epsilon_{\text{uzawa}}$ (`solver.uzawa_tol`), assurant une divergence rigoureusement nulle au pas de temps $n+1$.

### 2. Lagrangien Augmenté Adaptatif Incompressible (`AAL`)

Dans les écoulements diphasiques présentant de forts contrastes de masse volumique ($\rho_1/\rho_2 \ge 10^3$) et de viscosité ($\mu_1/\mu_2 \ge 10^2$), une valeur scalaire uniforme de $r$ dégrade considérablement le conditionnement du système linéaire. La méthode AAL (*Adaptative Augmented Lagrangian*, Vincent et al., *Computers & Fluids* 2004) détermine automatiquement un champ local de pénalité $r(t, \mathbf{x})$ cellule par cellule :

$$r_{i, j} = K_{\text{aal}} \cdot \max \left( \tau_{\text{instat}}, \tau_{\text{inert}}, \tau_{\text{grav}}, \tau_{\text{pres}}, \tau_{\text{visc}}, \tau_{\text{cap}}, \tau_{\text{brink}} \right)$$

où les temps caractéristiques physiques sont construits à partir des échelles globales du problème :
- Terme instationnaire : $\tau_{\text{instat}} = \alpha_{\text{time}} \rho_{i, j} \frac{L_0^2}{t_0}$
- Terme inertiel : $\tau_{\text{inert}} = \rho_{i, j} u_0 L_0$
- Terme gravitationnel : $\tau_{\text{grav}} = \rho_{i, j} \frac{L_0^2}{u_0} \|\mathbf{g}\|$
- Terme de pression : $\tau_{\text{pres}} = \frac{p_0 L_0}{u_0}$
- Terme visqueux : $\tau_{\text{visc}} = \mu_{i, j}$
- Terme capillaire : $\tau_{\text{cap}} = \frac{\sigma}{u_0}$
- Terme de Brinkman : $\tau_{\text{brink}} = \mu_{i, j} \frac{\Delta x \Delta y}{\mathcal{K}_{i, j}}$

Ce préconditionnement physique optimal assure une convergence quasi insensible aux sauts de propriétés interfaciales.

### 3. Lagrangien Augmenté Compressible Bas-Mach (`C-AAL`)

Pour les régimes compressibles bas-Mach (`solver.type = "C-AAL"`), le terme de compressibilité isotherme $\chi_T = \frac{1}{\rho c_s^2}$ modifie la contrainte de fermeture. Le paramètre adaptatif de pénalité est régularisé selon :

$$r_{\text{eff}, i, j} = \frac{r_{\text{AAL}, i, j}}{1 + r_{\text{AAL}, i, j} \frac{\chi_{T, i, j}}{\Delta t}}$$

et la mise à jour d'Uzawa devient :

$$p^{k+1}_{i, j} = p^k_{i, j} - r_{\text{eff}, i, j} \left[ (\nabla \cdot \mathbf{u}^{k+1})_{i, j} + \frac{\chi_{T, i, j}}{\Delta t} (p^k_{i, j} - p^n_{i, j}) \right]$$

### 4. Résolution Monolithique Tout-Couplé (`MONOLITHIC` / `FULLY_COUPLED`)

Dans la formulation monolithique, la vitesse et la pression sont assemblées simultanément au sein d'un système point-selle global sans paramètre de pénalité ($r = 0$) :

$$\begin{pmatrix} A_{\mathbf{u}} & B^T \\ B & -C_p \end{pmatrix} \begin{pmatrix} \mathbf{u}^{n+1} \\ p^{n+1} \end{pmatrix} = \begin{pmatrix} \mathbf{f}_{\mathbf{u}} \\ g_p \end{pmatrix}$$

où :
- $A_{\mathbf{u}}$ est le bloc de quantité de mouvement vectoriel discret.
- $B^T$ est la matrice du gradient de pression discret.
- $B$ est la matrice de divergence discrète ($B = - (B^T)^T$).
- $C_p = \text{diag}\left( \frac{\chi_T \Delta x \Delta y}{\Delta t} \right)$ est la matrice de compressibilité barotrope ($C_p = 0$ en incompressible strict).

Ce système est résolu de manière unifiée à chaque pas de temps via les solveurs de Krylov couplés au préconditionneur par blocs triangulaires.

### 5. Méthode de Projection Standard (`STANDARD_PROJECTION`)

La méthode de projection de Chorin-Temam sépare la résolution en deux sous-étapes découplées :
1. **Étape prédictive :** Calcul d'un champ de vitesse intermédiaire $\mathbf{u}^*$ sans gradient de pression.
2. **Étape de projection :** Résolution d'une équation de Poisson pour le pseudo-potentiel de pression $\phi$ :
   $$\nabla \cdot \left( \frac{1}{\rho} \nabla \phi \right) = \frac{\nabla \cdot \mathbf{u}^*}{\Delta t}$$
3. **Étape de correction solénoïdale :**
   $$\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla \phi, \quad p^{n+1} = p^n + \phi$$

---

## Solveurs Linéaires de Krylov Matrix-Free

COVOF 2D met en œuvre une infrastructure de solveurs de Krylov **sans matrice explicite (*Matrix-Free*)**. L'opérateur complet $A \mathbf{x}$ n'est jamais assemblé sous forme de matrice creuse globale en mémoire, mais calculé à la volée par passage direct des pochoirs volumes finis décalés sur les tableaux locaux avec synchronisation MPI des mailles fantômes.

### Algorithmes de Résolution Linéaire

1. **BiCGSTAB(2) (`LINEAR_BICGSTAB2`) :**
   Algorithme de Sleijpen et van der Vorst (1993). Il combine deux étapes de BiCG avec une étape de minimisation polynomiale locale d'ordre 2, offrant une robustesse exceptionnelle contre les stagnations sur les systèmes asymétriques à fort gradient convectif.
2. **BiCGSTAB (`LINEAR_BICGSTAB`) :**
   Algorithme standard de van der Vorst (1992), rapide et efficace pour les matrices asymétriques modérément raides.
3. **Gradient Conjugué Préconditionné (`LINEAR_CG`) :**
   Algorithme classique de Hestenes et Stiefel (1952). Strictement dédié aux opérateurs symétriques définis positifs (système de Poisson en pression ou étape de diffusion pure).
4. **GMRES(m) (`LINEAR_GMRES`) :**
   *Generalized Minimal Residual* redémarré de Saad et Schultz (1986). Minimise la norme euclidienne du résidu dans le sous-espace de Krylov orthogonalisé par Modified Gram-Schmidt.
5. **MINRES (`LINEAR_MINRES`) :**
   Algorithme de Paige et Saunders (1975). Dédié aux systèmes linéaires **symétriques indéfinis** tels que les systèmes point-selle Navier-Stokes de Stokes ou monolithiques. Repose sur la tridiagonalisation de Lanczos à trois termes et les rotations de Givens, garantissant une minimisation monotone de la norme du résidu avec un coût mémoire minimal.
6. **Solveur Direct Exact (`LINEAR_DIRECT` / `MUMPS`) :**
   Factorisation LU complète par blocs via les routines LAPACK standards `dgetrf` (décomposition $P A = L U$) et `dgetrs` (descente et remontée triangulaires). Pour les régimes stationnaires linéaires, la factorisation factorisée est mise en cache pour éliminer tout surcoût de factorisation aux pas de temps ultérieurs.

---

## Préconditionneurs Linéaires Parallèles

Le choix du préconditionneur est essentiel pour accélérer la convergence des méthodes de Krylov. COVOF 2D propose cinq préconditionneurs natifs :

1. **Jacobi Diagonal (`PRECOND_JACOBI` / `DIAGONAL`) :**
   Inversion de la diagonale enrichie de l'opérateur local. Cette diagonale intègre la contribution instationnaire, la diffusion visqueuse $\mathcal{O}(\mu/\Delta x^2)$, la pénalité de Lagrangien augmenté adaptatif $\mathcal{O}(r/\Delta x^2)$, la résistance de Brinkman $\mu/\mathcal{K}$ et les termes diagonaux de Newton.
2. **SSOR / Gauss-Seidel Symétrique (`PRECOND_SSOR`) :**
   Balayage avant (*forward sweep*) puis balayage arrière (*backward sweep*) sur la grille locale avec coefficient de relaxation optimal $\omega = 1.0$, offrant un lissage haute fréquence efficace.
3. **ILU(0) (`PRECOND_ILU0`) :**
   Factorisation incomplète sans remplissage sur les blocs de vitesse et pression.
4. **Blocs Triangulaires avec Complément de Schur (`PRECOND_BLOCK_TRIANGULAR` / `BT_SBP`) :**
   Préconditionneur spécifique dédié au système monolithique tout-couplé (El Ouafa et al., 2021/2023). Il s'appuie sur une factorisation triangulaire supérieure par blocs :
   $$P = \begin{pmatrix} D_{\mathbf{u}} & B^T \\ 0 & S_p \end{pmatrix}$$
   où $D_{\mathbf{u}}$ est l'approximation diagonale du bloc de quantité de mouvement et $S_p$ est le complément de Schur approché $S_p \approx B D_{\mathbf{u}}^{-1} B^T - C_p$, résolu par balayages de Gauss-Seidel symétriques.
5. **Préconditionneur Personnalisé (`PRECOND_CUSTOM`) :**
   Point d'entrée permettant à l'utilisateur d'injecter une fonction lambda C++ ou un pointeur de procédure Fortran (`custom_precond_hook`) pour tester de nouvelles heuristiques spectrales ou des couplages multigrilles externes.
# Géométrie, Maillages et Conditions aux Limites

## Génération et Topologie du Maillage Cartésien MAC

Le solveur **COVOF 2D** opère sur des domaines géométriques rectangulaires discrétisés par un maillage cartésien structuré uniforme selon la topologie décalée MAC (*Marker-and-Cell*).

### Définition Géométrique du Domaine

Le domaine physique de calcul $\Omega = [x_{\text{min}}, x_{\text{max}}] \times [y_{\text{min}}, y_{\text{max}}]$ est caractérisé dans le bloc `mesh` du fichier JSON par :
- Les dimensions physiques : $L_x = x_{\text{max}} - x_{\text{min}}$ (`mesh.lx`) et $L_y = y_{\text{max}} - y_{\text{min}}$ (`mesh.ly`) $[\text{m}]$.
- L'origine des coordonnées cartésiennes : $(x_{\text{min}}, y_{\text{min}})$ (`mesh.x_min`, `mesh.y_min`) $[\text{m}]$ (défaut : $0.0, 0.0$).
- Le nombre total de subdivisions spatiales globales : $N_x$ mailles selon l'axe horizontal $X$ (`mesh.nx`) et $N_y$ mailles selon l'axe vertical $Y$ (`mesh.ny`).

Les pas d'espace constants du maillage sont calculés immédiatement selon :

$$\Delta x = \frac{L_x}{N_x}, \quad \Delta y = \frac{L_y}{N_y}$$

### Coordonnées Spatiales des Nœuds, Faces et Centres

Pour une maille repérée par les indices discrets globaux $(i, j)$ où $i \in [1, N_x]$ et $j \in [1, N_y]$ :
- **Coordonnées du centre de maille (Pression, VOF, Température) :**
  $$x_{c, i} = x_{\text{min}} + \left( i - \frac{1}{2} \right) \Delta x, \quad y_{c, j} = y_{\text{min}} + \left( j - \frac{1}{2} \right) \Delta y$$
- **Coordonnées des faces verticales (Vitesse $u$) :**
  $$x_{f, i} = x_{\text{min}} + (i - 1) \Delta x, \quad y_{c, j} = y_{\text{min}} + \left( j - \frac{1}{2} \right) \Delta y$$
  Il y a au total $(N_x + 1) \times N_y$ faces de vitesse horizontale dans le domaine complet.
- **Coordonnées des faces horizontales (Vitesse $v$) :**
  $$x_{c, i} = x_{\text{min}} + \left( i - \frac{1}{2} \right) \Delta x, \quad y_{f, j} = y_{\text{min}} + (j - 1) \Delta y$$
  Il y a au total $N_x \times (N_y + 1)$ faces de vitesse verticale dans le domaine complet.

---

## Décomposition de Domaine et Parallélisme Cartésien 2D (MPI)

Pour permettre l'exécution sur architectures distribuées (grappes de calcul et supercalculateurs), COVOF 2D partitionne le maillage global selon une grille 2D de processus MPI communicants.

```
       +-----------------------+-----------------------+
       |                       |                       |
       |     Rang MPI (0, 1)   |     Rang MPI (1, 1)   |
       |     Sous-domaine      |     Sous-domaine      |
       |                       |                       |
Ny / 2 +-----------------------+-----------------------+
       |                       |                       |
       |     Rang MPI (0, 0)   |     Rang MPI (1, 0)   |
       |     Sous-domaine      |     Sous-domaine      |
       |                       |                       |
       +-----------------------+-----------------------+
                   0                    Nx / 2                  Nx
```

### Topologie Virtuelle Cartésienne (`MPI_Cart_create`)

1. **Calcul automatique des dimensions de partitionnement :**
   À partir du nombre total de processus alloués $N_{\text{procs}}$, COVOF 2D détermine la décomposition optimale en $P_x$ tranches horizontales et $P_y$ tranches verticales via `MPI_Dims_create(n_procs, 2, dims)` :
   $$N_{\text{procs}} = P_x \times P_y \quad \text{avec} \quad \frac{P_x}{P_y} \approx \frac{N_x}{N_y}$$
2. **Identification des rangs voisins :**
   Grâce à `MPI_Cart_shift`, chaque processus détermine les rangs de ses quatre voisins directs dans la grille :
   - Voisin Ouest (`rank_west`) à la coordonnée $(p_x - 1, p_y)$
   - Voisin Est (`rank_east`) à la coordonnée $(p_x + 1, p_y)$
   - Voisin Sud (`rank_south`) à la coordonnée $(p_x, p_y - 1)$
   - Voisin Nord (`rank_north`) à la coordonnée $(p_x, p_y + 1)$
   Si un processus est situé sur une frontière physique sans condition périodique, le rang voisin correspondant est positionné à `MPI_PROC_NULL`, désactivant automatiquement les transferts réseau sans conditionnement lourd.

### Gestion des Cellules Fantômes (*Halos / Ghost Cells*)

Chaque sous-domaine local alloue une couche périphérique de cellules de garde :
- `nx_local` $\times$ `ny_local` : mailles intérieures possédées en propre par le processus.
- `nx_ghost = nx_local + 2` et `ny_ghost = ny_local + 2` : taille mémoire totale incluant les cellules fantômes d'indice $0$ et `nx_local + 1`.

À chaque étape de calcul nécessitant des dérivées spatiales ou des flux aux interfaces de sous-domaines :
1. Les données des mailles intérieures de bordure sont emballées dans des tampons contigus (*packing*).
2. Les communications asynchrones non bloquantes `MPI_Isend` et `MPI_Irecv` sont initiées vers les voisins E, W, N, S.
3. Après synchronisation par `MPI_Waitall`, les valeurs reçues sont déballées dans les cellules fantômes (*unpacking*).

---

## Catalogue Exhaustif des Conditions aux Limites (CL)

COVOF 2D dispose d'un catalogue complet de conditions aux limites physiques et numériques applicables indépendamment sur les quatre frontières cardinales du domaine : `left` (Ouest, $x = x_{\text{min}}$), `right` (Est, $x = x_{\text{max}}$), `bottom` (Sud, $y = y_{\text{min}}$) et `top` (Nord, $y = y_{\text{max}}$).

### 1. Vitesse Imposée (`DIRICHLET`)

La condition de Dirichlet prescrit explicitement la valeur vectorielle de la vitesse à la frontière :

$$\mathbf{u} = \mathbf{u}_{\text{wall}} = (u_{\text{bc}}, v_{\text{bc}})^T \quad \text{sur } \partial \Omega$$

#### Implémentation Numérique Décalée
- **Pour la composante normale :**
  La vitesse normale coïncide géométriquement avec la face frontière. Elle est imposée directement sans interpolation :
  $$u_{1/2, j} = u_{\text{bc}} \quad \text{(sur frontière gauche)}$$
- **Pour la composante tangentielle :**
  La vitesse tangentielle $v$ est localisée à une distance $\frac{\Delta x}{2}$ de la paroi. Pour imposer la valeur $v_{\text{bc}}$ sur la frontière physique avec une précision d'ordre 2, la cellule fantôme d'indice $i=0$ est extrapolée linéairement :
  $$\frac{v_{0, j+1/2} + v_{1, j+1/2}}{2} = v_{\text{bc}} \implies v_{0, j+1/2} = 2 v_{\text{bc}} - v_{1, j+1/2}$$

#### Cas d'Usage Typiques
- **Paroi solide immobile sans glissement (*No-Slip Wall*) :** $u_{\text{bc}} = 0, v_{\text{bc}} = 0$.
- **Paroi mobile de cavité entraînée (*Lid-Driven Cavity*) :** sur la frontière supérieure, $u_{\text{bc}} = U_{\text{lid}} = 1.0\,\text{m}\cdot\text{s}^{-1}, v_{\text{bc}} = 0$.
- **Profil d'entrée parabolique de Poiseuille :** $u(y) = 4 U_{\max} \frac{y (L_y - y)}{L_y^2}$.

### 2. Gradient Nul / Sortie Libre (`NEUMANN`)

La condition de Neumann prescrit un flux nul (dérivée normale nulle) des grandeurs de vitesse à travers la frontière :

$$\frac{\partial \mathbf{u}}{\partial \mathbf{n}} = \mathbf{0} \quad \iff \quad \frac{\partial u}{\partial n} = 0, \quad \frac{\partial v}{\partial n} = 0$$

#### Implémentation Numérique Décalée
Sur la frontière droite ($x = x_{\text{max}}$) :
$$u_{N_x+1/2, j} = u_{N_x-1/2, j}, \quad v_{N_x+1, j+1/2} = v_{N_x, j+1/2}$$
Pour la pression, la condition de Neumann homogène $\frac{\partial p}{\partial n} = 0$ est naturellement imposée pour assurer la compatibilité avec l'équation de quantité de mouvement.

#### Cas d'Usage Typiques
- **Sortie libre convective (*Outflow / Open Boundary*) :** Permet à l'écoulement et aux tourbillons de s'évacuer du domaine de calcul sans réflexion d'ondes numériques.

### 3. Périodicité Spatiale (`PERIODIC`)

La condition de périodicité relie topologiquement les frontières opposées du domaine pour simuler un milieu infini ou un écoulement en boucle fermée :

$$\mathbf{u}(x_{\text{min}}, y) = \mathbf{u}(x_{\text{max}}, y), \quad p(x_{\text{min}}, y) = p(x_{\text{max}}, y)$$

#### Implémentation Numérique Décalée
Les cellules fantômes de la frontière Ouest reçoivent les données des dernières cellules intérieures de la frontière Est, et réciproquement. Dans le module MPI, le communicateur cartésien active le drapeau de périodicité `periods = [true, false]`, de sorte que les processus de rang $0$ communiquent directement avec les processus de rang $P_x - 1$.

#### Cas d'Usage Typiques
- **Vortex de Taylor-Green et instabilités de cisaillement.**
- **Canaux infinis sous gradient de pression moteur.**

### 4. Glissement Parfait sans Frottement (`FREESLIP`)

La condition de glissement parfait impose une vitesse normale nulle (imperméabilité) et un cisaillement tangentiel nul (aucun frottement visqueux) :

$$\mathbf{u} \cdot \mathbf{n} = 0, \quad \frac{\partial (\mathbf{u} \times \mathbf{n})}{\partial \mathbf{n}} = \mathbf{0}$$

#### Implémentation Numérique Décalée
Sur la frontière inférieure ($y = y_{\text{min}}$, normale selon $y$) :
- Vitesse normale nulle sur la face frontière : $v_{i, 1/2} = 0$.
- Extrapolation symétrique de la vitesse tangentielle dans la maille fantôme :
  $$u_{i+1/2, 0} = u_{i+1/2, 1} \implies \left( \frac{\partial u}{\partial y} \right)_{y=y_{\text{min}}} = 0$$

#### Cas d'Usage Typiques
- **Plans de symétrie géométrique et hydrodynamique.**
- **Surfaces libres idéalisées et fonds marins sans couche limite visqueuse.**

---

## Conditions aux Limites Thermiques et VOF

En complément des conditions hydrodynamiques, les frontières physiques accueillent les conditions sur les scalaires passifs et actifs :

| Condition aux Limites | Mot-Clé JSON | Équation Continue | Implémentation Fantôme |
| :--- | :--- | :--- | :--- |
| **Température Imposée** | `THERMAL_DIRICHLET` | $T = T_{\text{wall}}$ | $T_{\text{ghost}} = 2 T_{\text{wall}} - T_{\text{interior}}$ |
| **Paroi Adiabatique** | `THERMAL_NEUMANN` | $\frac{\partial T}{\partial n} = 0$ | $T_{\text{ghost}} = T_{\text{interior}}$ |
| **Interface Imperméable** | `VOF_WALL` | $\mathbf{u} \cdot \mathbf{n}_{\text{wall}} = 0$ | Flux VOF nul à travers la face |
| **Sortie de Phase Libre** | `VOF_NEUMANN` | $\frac{\partial C}{\partial n} = 0$ | $C_{\text{ghost}} = C_{\text{interior}}$ |

---

## Modélisation Géométrique des Obstacles Immergés

COVOF 2D intègre un générateur géométrique interne d'obstacles solides immergés pénalisés par la méthode de Darcy-Brinkman (`obstacles.enable_obstacles = true`). Deux formes canoniques sont disponibles :

### 1. Obstacle Cylindrique (`CYLINDER`)

L'obstacle est défini par le centre de son disque $(x_c, y_c)$ (`obs_xc`, `obs_yc`) et son rayon $R$ (`obs_radius`) :

$$\Omega_s = \left\{ (x, y) \in \Omega \;\middle|\; (x - x_c)^2 + (y - y_c)^2 \le R^2 \right\}$$

- Pour chaque face verticale $u_{i+1/2, j}$, le test d'appartenance `is_u_face_in_obstacle` vérifie si le centre de face $(x_{f, i+1/2}, y_{c, j})$ est à l'intérieur du cylindre.
- Si le point est immergé et que la perméabilité `obs_perm` $\le 10^{-12}$, la composante de vitesse est directement annulée ($u = 0$, Dirichlet solide) ou pénalisée par le coefficient de Brinkman $\frac{\mu}{K}$.

### 2. Obstacle Rectangulaire / Plaque (`RECTANGLE`)

L'obstacle est défini par sa boîte englobante $[x_0, x_1] \times [y_0, y_1]$ (`obs_x0`, `obs_x1`, `obs_y0`, `obs_y1`) :

$$\Omega_s = \left\{ (x, y) \in \Omega \;\middle|\; x_0 \le x \le x_1 \quad \text{et} \quad y_0 \le y \le y_1 \right\}$$

Cette géométrie permet de représenter facilement des marches descendantes (*backward-facing step*), des chicanes, des parois internes minces ou des piliers carrés.
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
# Exécution, Post-Traitement et Dépannage

## Lancement d'une Simulation et Gestion du Calcul

Le lancement d'une simulation avec **COVOF 2D** est conçu pour être direct, rapide et unifié entre les deux implémentations C++17 et Fortran 2018.

### Lancement en Ligne de Commande

Les deux exécutables acceptent le chemin vers le fichier de configuration JSON comme argument unique en ligne de commande :

```bash
# Exécution séquentielle C++ :
./COVOF_2D-C++/bin/covof2d_cpp chemin/vers/mon_cas.json

# Exécution séquentielle Fortran :
./COVOF_2D-FORTRAN/bin/covof2d_fortran chemin/vers/mon_cas.json

# Exécution parallèle distribuée MPI (exemple sur 4 processus) :
mpirun -np 4 ./COVOF_2D-C++/bin/covof2d_cpp chemin/vers/mon_cas.json
mpirun -np 4 ./COVOF_2D-FORTRAN/bin/covof2d_fortran chemin/vers/mon_cas.json
```

Si aucun argument n'est fourni, les deux solveurs recherchent par défaut un fichier nommé `config.json` dans le répertoire d'exécution courant. En cas d'absence de ce fichier, ils affichent une notice d'aide claire récapitulant les options d'usage et les cas d'exemples.

---

## Interprétation des Journaux de Convergence (Logs)

Dès son initialisation, COVOF 2D affiche une bannière technique récapitulant les paramètres de la simulation : maillage local et global, topologie des processus MPI, propriétés thermodynamiques, choix des solveurs et pas de temps initial.

### Exemple de Journal Console au Démarrage

```text
================================================================================
   COVOF 2D — Solveur CFD Navier-Stokes Incompressible Parallèle (C++17 / MPI)   
================================================================================
[INFO] Lecture de la configuration : case_dambreak.json
[INFO] Topologie MPI 2D active     : 4 processus (Px = 2, Py = 2)
[INFO] Maillage global             : 128 x 256 mailles (Lx = 0.1 m, Ly = 0.2 m)
[INFO] Maillage local par rang     : 64 x 128 mailles (dx = 7.8125e-04 m, dy = 7.8125e-04 m)
[INFO] Modèle d'écoulement         : Diphasique VOF-PLIC (rho1=1000, rho2=1.2, mu1=1e-3, mu2=1.8e-5)
[INFO] Couplage Navier-Stokes      : AAL (Adaptative Augmented Lagrangian, Vincent et al. 2004)
[INFO] Solveur Linéaire            : BiCGSTAB(2) avec préconditionneur JACOBI (tol = 1.0e-06)
================================================================================
  Step        Time [s]        dt [s]    Iter_NS         Res_U        Res_Div         CFL
--------------------------------------------------------------------------------------
     1     5.00000e-04   5.00000e-04         12   1.45821e-07   8.91420e-06     0.0124
     2     1.00000e-03   5.00000e-04         14   2.11045e-07   7.23410e-06     0.0241
    10     5.00000e-03   5.00000e-04         15   3.45120e-07   4.12050e-06     0.0895
    50     2.50000e-02   5.00000e-04         18   5.12040e-07   2.89120e-06     0.2450
```

### Signification des Colonnes du Journal

- `Step` : Numéro de l'itération temporelle courante $n$.
- `Time [s]` : Temps physique total simulé écoulé depuis le début du calcul $t = \sum \Delta t$.
- `dt [s]` : Valeur du pas de temps d'avancement $\Delta t$ utilisé pour le pas courant.
- `Iter_NS` : Nombre cumulé d'itérations effectuées par le solveur linéaire de Krylov (ex. BiCGSTAB2 ou MINRES) sur l'ensemble des sous-itérations d'Uzawa du pas.
- `Res_U` : Norme relative euclidienne du résidu du système linéaire de quantité de mouvement à la convergence.
- `Res_Div` : Norme infinie maximale du résidu de continuité $\|\nabla \cdot \mathbf{u}^{n+1}\|_\infty$ sur l'ensemble du domaine de calcul $[\text{s}^{-1}]$.
- `CFL` : Nombre de Courant-Friedrichs-Lewy convectif maximal observé sur la grille : $\max_{i, j} \left( \frac{|u|\Delta t}{\Delta x} + \frac{|v|\Delta t}{\Delta y} \right)$.

### Fichier d'Historique de Convergence CSV

En parallèle de la sortie console, le solveur enregistre à chaque pas de temps une ligne dans le fichier `output/convergence.csv`, permettant le tracé dynamique en cours de calcul des courbes de résidus et du suivi du bilan de masse global.

---

## Visualisation et Extraction de Données

COVOF 2D produit des sorties volumiques au standard ouvert **XML VTK Rectilinear Grid**, optimisé pour les maillages cartésiens et directement exploitable sous **ParaView**, **VisIt** ou par les scripts Python de traitement d'images et de données scientifiques.

### Structure des Fichiers XML VTK Générés

Lorsqu'un pas d'écriture est franchi (selon la fréquence `output.frequency`), le répertoire `output/` reçoit :

1. **Les fichiers de sous-domaines `.vtr` (*VTK Rectilinear Grid*) :**
   Chaque processus MPI de rang $p$ écrit de manière autonome son sous-domaine local sans goulot d'étranglement centralisé :
   `output/nom_00100_piece_0.vtr`, `output/nom_00100_piece_1.vtr`, etc.
   Ces fichiers contiennent les tableaux de coordonnées 1D $x(i)$, $y(j)$ et les champs volumiques au format binaire ou base64 :
   - `Velocity` : vecteur vitesse $(u, v, 0)$ aux centres ou interpolé aux nœuds.
   - `Pressure` : champ scalaire de pression $p$.
   - `VOF_Color` : fraction volumique de liquide $C \in [0, 1]$.
   - `Density` et `Viscosity` : champs de masse volumique $\rho$ et viscosité $\mu$.
   - `Temperature` : champ thermique $T$ (si l'énergie est activée).

2. **Le fichier maître `.pvtr` (*Parallel VTK Rectilinear Grid*) :**
   Le processus de rang maître (rang 0) génère un fichier XML léger pointant vers les différentes pièces `.vtr` associées à chaque rang :
   ```xml
   <?xml version="1.0"?>
   <VTKFile type="PRectilinearGrid" version="0.1" byte_order="LittleEndian">
     <PRectilinearGrid WholeExtent="0 128 0 256 0 0" GhostLevel="0">
       <PCoordinates>
         <PDataArray type="Float64" Name="XCoordinates"/>
         <PDataArray type="Float64" Name="YCoordinates"/>
         <PDataArray type="Float64" Name="ZCoordinates"/>
       </PCoordinates>
       <PCellData Scalars="Pressure,VOF_Color" Vectors="Velocity">
         <PDataArray type="Float64" Name="Velocity" NumberOfComponents="3"/>
         <PDataArray type="Float64" Name="Pressure"/>
         <PDataArray type="Float64" Name="VOF_Color"/>
       </PCellData>
       <Piece Extent="0 64 0 128 0 0" Source="dambreak_00100_piece_0.vtr"/>
       <Piece Extent="64 128 0 128 0 0" Source="dambreak_00100_piece_1.vtr"/>
       <Piece Extent="0 64 128 256 0 0" Source="dambreak_00100_piece_2.vtr"/>
       <Piece Extent="64 128 128 256 0 0" Source="dambreak_00100_piece_3.vtr"/>
     </PRectilinearGrid>
   </VTKFile>
   ```

3. **Le fichier de collection temporelle `.pvd` (*ParaView Data*) :**
   Ce fichier liste la séquence temporelle complète, permettant à ParaView de charger l'animation en un seul clic :
   ```xml
   <?xml version="1.0"?>
   <VTKFile type="Collection" version="0.1">
     <Collection>
       <DataSet timestep="0.000" file="dambreak_00000.pvtr"/>
       <DataSet timestep="0.010" file="dambreak_00020.pvtr"/>
       <DataSet timestep="0.020" file="dambreak_00040.pvtr"/>
     </Collection>
   </VTKFile>
   ```

### Exploitation sous ParaView

1. **Chargement :** Ouvrir le fichier `.pvd` dans ParaView et cliquer sur **Apply**.
2. **Visualisation de l'interface liquide/gaz :**
   Appliquer le filtre **Contour** sur le champ scalaire `VOF_Color` avec la valeur seuil $C = 0.5$. La surface ainsi extraite représente fidèlement la frontière matérielle entre les deux phases.
3. **Vecteurs vitesse et lignes de courant :**
   Appliquer le filtre **Glyph** (flèches 2D orientées selon `Velocity`) ou **Stream Tracer** pour visualiser les tourbillons de recirculation dans la cavité ou le sillage d'obstacle.

### Script Python d'Extraction Automatisée et de Tracé de Profils

Voici un script Python complet (`extract_profiles.py`) s'appuyant sur les bibliothèques scientifiques `numpy`, `matplotlib` et `pyvista` pour extraire et tracer les profils de vitesse au centre du domaine :

```python
#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

def extract_mid_profile(pvtr_file):
    # Chargement du fichier parallèle VTK
    reader = pv.get_reader(pvtr_file)
    data = reader.read()

    # Définition de la ligne d'échantillonnage verticale au centre X = Lx / 2
    bounds = data.bounds  # [xmin, xmax, ymin, ymax, zmin, zmax]
    x_mid = 0.5 * (bounds[0] + bounds[1])
    p_start = [x_mid, bounds[2], 0.0]
    p_end   = [x_mid, bounds[3], 0.0]

    # Extraction des données le long du segment
    sampled_line = data.sample_over_line(p_start, p_end, resolution=500)
    
    y_coords = sampled_line.points[:, 1]
    u_velocity = sampled_line['Velocity'][:, 0]
    pressure = sampled_line['Pressure']

    # Tracé graphique haute définition
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), dpi=300)

    ax1.plot(u_velocity, y_coords, 'b-', lw=2, label='COVOF 2D')
    ax1.set_xlabel('Vitesse horizontale u [m/s]')
    ax1.set_ylabel('Hauteur y [m]')
    ax1.set_title(f'Profil de Vitesse à x = {x_mid:.3f} m')
    ax1.grid(True, ls='--', alpha=0.7)
    ax1.legend()

    ax2.plot(pressure, y_coords, 'r-', lw=2, label='Pression hydrodynamique')
    ax2.set_xlabel('Pression p [Pa]')
    ax2.set_ylabel('Hauteur y [m]')
    ax2.set_title(f'Profil de Pression à x = {x_mid:.3f} m')
    ax2.grid(True, ls='--', alpha=0.7)
    ax2.legend()

    plt.tight_layout()
    plt.savefig('profil_vertical_covof2d.png')
    print("Graphique sauvegardé avec succès dans profil_vertical_covof2d.png")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python3 extract_profiles.py <fichier.pvtr>")
        sys.exit(1)
    extract_mid_profile(sys.argv[1])
```

---

## Guide Complet de Dépannage (Troubleshooting)

Le tableau suivant présente une matrice exhaustive des anomalies numériques potentielles, leurs causes sous-jacentes et les actions correctives recommandées :

| Symptôme Observé | Cause Physique / Numérique Probable | Action Corrective Recommandée |
| :--- | :--- | :--- |
| **Divergence du solveur de Krylov (`Iter_NS` atteint `max_iter`, résidu NaN)** | 1. Pas de temps $\Delta t$ trop grand violant la condition CFL visqueuse ou convective.<br>2. Paramètre de pénalité de Lagrangien $r_{\text{al}}$ excessif détériorant le conditionnement spectral.<br>3. Absence de préconditionneur sur maillage fin. | 1. Réduire $\Delta t$ ou activer `time.adaptative_dt = true` avec `time.cfl = 0.3`.<br>2. Activer le mode adaptatif `solver.type = "AAL"` au lieu d'une valeur scalaire fixe `r_al`.<br>3. Activer `solver.preconditioner = "JACOBI"` ou `"BLOCK_TRIANGULAR"`. |
| **Résidu de continuité élevé ($\|\nabla \cdot \mathbf{u}\|_\infty > 10^{-3}$)** | Nombre insuffisant d'itérations d'Uzawa pour assurer la fermeture de l'incompressibilité. | Augmenter `solver.uzawa_iter` à $3$ ou $5$, et abaisser `solver.uzawa_tol` à $1.0\times 10^{-5}$. |
| **Courants parasites / Déformation artificielle de l'interface diphasique** | Mauvaise résolution de la courbure $\kappa$ dans le modèle CSF lorsque le maillage est trop grossier par rapport au rayon de courbure. | S'assurer que le diamètre de goutte ou bulle comporte au moins $20$ à $30$ cellules de calcul ($2R/\Delta x \ge 25$). Vérifier que le pas de temps respecte la limite d'ondes capillaires $\Delta t \le \sqrt{\frac{\rho \Delta x^3}{4\pi\sigma}}$. |
| **Oscillations de vitesse en damier (*checkerboard*)** | Présence d'un schéma convectif d'ordre 2 centré sans maillage suffisant ($Re_\Delta > 2$). | Basculer `convection.convective_scheme = "UPWIND"` ou raffiner le maillage localement pour ramener le nombre de Reynolds de maille sous la limite de Peclet $Pe \le 2$. |
| **Blocage ou dépassement de temps (*Timeout*) lors de l'exécution MPI** | 1. Décomposition de domaine incompatible avec les dimensions du maillage (mailles locales nulles : `nx_local = 0`).<br>2. Blocage sur barrière collective lors d'un calcul thermique ou projection. | 1. Vérifier que $N_x \ge 4 \times P_x$ et $N_y \ge 4 \times P_y$.<br>2. S'assurer que le communicateur MPI cartésien a été initialisé sans rangs orphelins. |
| **Perte de conservation de la masse globale de fluide** | Schéma d'advection VOF non conservatif ou pas de temps violant la condition de conservation volumique. | Vérifier que `physics.vof_model` est configuré sur `"PLIC"` et que le pas de temps respecte $u \Delta t / \Delta x \le 0.5$. |
