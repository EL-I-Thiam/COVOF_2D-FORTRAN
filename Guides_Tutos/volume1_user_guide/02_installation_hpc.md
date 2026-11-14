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
