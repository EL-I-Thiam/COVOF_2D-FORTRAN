# Installation, Compilation et Déploiement HPC

## Prérequis Logiciels et Matrice de Compatibilité

Le déploiement de **COVOF 2D** sur stations de travail, serveurs de calcul ou supercalculateurs massivement parallèles requiert uniquement des outils standards de l'écosystème scientifique. Les deux branches jumelles (`COVOF_2D-C++` et `COVOF_2D-FORTRAN`) ont été conçues sans dépendance JSON externe (analyseurs JSON natifs intégrés dans `Config.cpp` et `module_json.f90`) et reposent exclusivement sur **MPI pur** (*Message Passing Interface*) et une implémentation LAPACK/BLAS pour le solveur direct LU optionnel.

### Matrice de Compatibilité des Compilateurs

Les deux versions adhèrent strictement aux standards internationaux **ISO C++17** et **ISO Fortran 2018** :

| Famille de Compilateurs | Version Minimale | Compilateur C++ | Compilateur Fortran | Drapeaux CMake / Compilation |
| :--- | :--- | :--- | :--- | :--- |
| **GNU Compiler Collection (GCC)** | $\ge 9.3$ (recommandé $\ge 11.2$) | `g++` | `gfortran` | C++ : `-O3 -Wall -Wextra -Wno-unused-parameter` <br> Fortran : `-O3 -march=native` |
| **LLVM / Clang + Flang/GFortran** | $\ge 14.0$ | `clang++` | `gfortran` / `flang` | `-O3 -march=native` |
| **Intel oneAPI** | $\ge 2022.1$ | `icpx` | `ifx` | `-O3 -xHost` |
| **Apple Clang + Homebrew GCC (macOS)** | Xcode CLT $\ge 14$ | `clang++` | `gfortran` (Homebrew) | Liaison automatique `-framework Accelerate` |

### Bibliothèques Requises et Modèle de Parallélisme

1. **MPI pur (OpenMPI / MPICH / Intel MPI) :**
   - COVOF 2D utilise un **parallélisme distribué MPI pur** (`MPI::MPI_CXX` en C++, `MPI::MPI_Fortran` en Fortran), basé sur une topologie cartésienne 2D (`MPI_Cart_create`).
   - **Aucun thread OpenMP n'est utilisé par le solveur** : les fichiers `CMakeLists.txt` n'activent pas `-fopenmp`, et il n'existe aucune région parallèle `#pragma omp` ou `!$omp parallel` dans les noyaux de calcul. Chaque processus MPI gère un sous-domaine mémoire indépendant.
   - Implémentations validées : **OpenMPI** ($\ge 4.0$, référence principale du projet), **MPICH** ($\ge 3.4$), **Intel MPI** ($\ge 2021$).
2. **Algèbre Linéaire Dense (LAPACK) :**
   - Requise uniquement lorsque le solveur direct exact (`solver.linear_solver = "DIRECT"` ou `"MUMPS"` en repli autonome) invoque la factorisation LU dense `dgetrf_` / `dgetrs_`.
   - Sous **macOS (`APPLE`)** : `CMakeLists.txt` lie automatiquement le framework natif `-framework Accelerate`.
   - Sous **Linux / Clusters HPC** : lier `OpenBLAS` / `LAPACK` (`liblapack-dev`, `libopenblas-dev`) ou Intel MKL selon l'environnement système.
3. **CMake et Python (post-traitement / vérification) :**
   - **CMake $\ge 3.16$** pour `COVOF_2D-C++` et **$\ge 3.18$** pour `COVOF_2D-FORTRAN`.
   - **Python 3** avec `numpy`, `matplotlib`, `pyvista` (optionnel, pour le post-traitement et `verify_cvi_markers.py`).

---

## Compilation Pas à Pas avec CMake

Les deux dépôts disposent chacun de leur fichier `CMakeLists.txt` autonome à leur racine respective.

### Compilation de la Branche C++17 (`COVOF_2D-C++`)

Le fichier `COVOF_2D-C++/CMakeLists.txt` impose `CMAKE_CXX_STANDARD 17` (`REQUIRED ON`), active `CMAKE_EXPORT_COMPILE_COMMANDS ON`, recherche `MPI REQUIRED` et compile les 13 unités `.cpp` vers l'exécutable `COVOF_2D-C++/bin/covof2d_cpp` tout en créant automatiquement un lien symbolique de commodité `COVOF_2D-C++/bin/covof2d` :

```bash
# 1. Se positionner dans le répertoire C++
cd /Users/elhadjiibrahimathiam/Documents/PROJET_COVOF-2D/COVOF_2D-C++

# 2. Créer et entrer dans le dossier de construction hors-source
mkdir -p build && cd build

# 3. Configurer avec CMake en mode Release
cmake -DCMAKE_BUILD_TYPE=Release ..

# 4. Compiler en parallèle sur les cœurs disponibles
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)
```

À l'issue de l'édition de liens, vérifier la présence des binaires :

```bash
ls -la ../bin/
# Doit afficher :
#   covof2d_cpp
#   covof2d -> .../bin/covof2d_cpp
```

#### Options CMake et Drapeaux Utiles pour C++

| Option CMake | Valeur / Comportement | Rôle |
| :--- | :--- | :--- |
| `CMAKE_BUILD_TYPE` | `Release`, `Debug`, `RelWithDebInfo` | `Release` active `-O3` ; `Debug` permet l'ajout de `-g -O0 -fsanitize=address,undefined` |
| `CMAKE_CXX_STANDARD` | `17` (`REQUIRED ON`) | Standard ISO C++17 requis pour `std::filesystem`, `std::clamp`, etc. |
| `CMAKE_CXX_COMPILER` | ex. `mpicxx` ou `g++` | Surcharge explicite du compilateur C++ si plusieurs suites coexistent |
| `CMAKE_EXPORT_COMPILE_COMMANDS` | `ON` | Génère `build/compile_commands.json` pour les analyseurs LSP (`clangd`) |

### Compilation de la Branche Fortran 2018 (`COVOF_2D-FORTRAN`)

Le fichier `COVOF_2D-FORTRAN/CMakeLists.txt` impose `CMAKE_Fortran_STANDARD 2018` (`REQUIRED ON`), ajoute les drapeaux `-O3 -march=native`, recherche `MPI REQUIRED` et compile les 16 fichiers `.f90` dans l'ordre exact des dépendances modulaires vers l'exécutable `COVOF_2D-FORTRAN/bin/covof2d_fortran` :

```bash
# 1. Se positionner dans le répertoire Fortran
cd /Users/elhadjiibrahimathiam/Documents/PROJET_COVOF-2D/COVOF_2D-FORTRAN

# 2. Créer et entrer dans le dossier de construction
mkdir -p build && cd build

# 3. Configurer avec CMake en mode Release
cmake -DCMAKE_BUILD_TYPE=Release ..

# 4. Compiler en parallèle
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)
```

:::warning
Sur cluster hétérogène (nœuds de compilation frontaux d'architecture différente des nœuds de calcul), le drapeau `-march=native` de `COVOF_2D-FORTRAN/CMakeLists.txt` doit être remplacé par l'architecture cible des nœuds de calcul (par exemple `-march=skylake-avx512` ou `-xCORE-AVX512`) afin d'éviter une erreur `Illegal Instruction (SIGILL)` à l'exécution.
:::

### Drapeaux d'Optimisation SIMD et Résolution de LAPACK sous Linux

Sur architecture Linux x86_64 pour nœuds Xeon ou EPYC, si `dgetrf_`/`dgetrs_` nécessite une liaison explicite avec OpenBLAS ou LAPACK :

```bash
# Exemple C++ et Fortran avec drapeaux AVX-512 ciblés et LAPACK explicite :
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_FLAGS="-O3 -march=skylake-avx512" \
      -DCMAKE_EXE_LINKER_FLAGS="-llapack -lblas" ..
```

---

## Déploiement sur Supercalculateurs (Clusters HPC)

Dans un environnement distribué piloté par un ordonnanceur (*workload manager*), **COVOF 2D étant un code MPI pur**, chaque tâche Slurm ou PBS correspond à **1 processus MPI sur 1 cœur physique** (`--cpus-per-task=1`). L'utilisation de ressources multi-cœurs par rang n'apporte aucun gain, car le code ne crée pas de threads applicatifs OpenMP (seule l'éventuelle bibliothèque BLAS liée lors de l'appel à `DIRECT` doit être contrainte à un thread si elle est multi-threadée, par exemple `OPENBLAS_NUM_THREADS=1` ou `VECLIB_MAXIMUM_THREADS=1` sous macOS).

### Script de Soumission pour Gestionnaire Slurm

Voici un script de soumission complet pour **Slurm** (`run_covof_slurm.sh`), allouant 2 nœuds de 64 cœurs (soit 128 processus MPI purs) :

```bash
#!/bin/bash
#SBATCH --job-name=COVOF2D_PROD
#SBATCH --output=covof2d_%j.log
#SBATCH --error=covof2d_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --exclusive

# Chargement de l'environnement logiciel
module purge
module load gcc/12.2.0
module load openmpi/4.1.5
module load openblas/0.3.21

echo "=========================================================="
echo "Démarrage COVOF 2D (MPI pur) sur Slurm"
echo "Job ID                     : $SLURM_JOB_ID"
echo "Nombre de nœuds            : $SLURM_NNODES"
echo "Nombre total de rangs MPI  : $SLURM_NTASKS"
echo "Liste des nœuds            : $SLURM_JOB_NODELIST"
echo "=========================================================="

# Neutralisation d'un éventuel multi-threading implicite des bibliothèques BLAS
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Transport réseau haute vitesse UCX (InfiniBand)
export UCX_TLS=rc,mm,sm
export OMPI_MCA_pml=ucx

WORKDIR=$SLURM_SUBMIT_DIR
CASE_CONFIG="$WORKDIR/case_config.json"
cd "$WORKDIR"

# Lancement parallèle MPI pur (C++17 ou Fortran 2018)
srun --mpi=pmix --cpu-bind=cores ./COVOF_2D-C++/bin/covof2d_cpp "$CASE_CONFIG"

echo "=========================================================="
echo "Fin du calcul COVOF 2D - Code retour : $?"
echo "=========================================================="
```

### Script de Soumission pour Gestionnaire PBS / Torque

Pour les infrastructures s'appuyant sur **PBS Professional** ou **Torque** :

```bash
#!/bin/bash
#PBS -N COVOF2D_PBS
#PBS -l select=2:ncpus=64:mpiprocs=64
#PBS -l walltime=08:00:00
#PBS -q standard
#PBS -j oe
#PBS -o covof2d_pbs.log

cd "$PBS_O_WORKDIR"

module purge
module load intel-oneapi-compilers/2023.1.0
module load intel-oneapi-mpi/2021.9.0
module load intel-oneapi-mkl/2023.1.0

export MKL_NUM_THREADS=1
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core
export I_MPI_FABRICS=shm:ofi

mpirun -np 128 ./COVOF_2D-FORTRAN/bin/covof2d_fortran case_config.json
```

### Gestion de l'Affinité Mémoire et Équilibre de Charge Cartésien

Lors du partitionnement par `MPI_Cart_create` / `MPI_Dims_create`, il est crucial de choisir un nombre de processus $P = P_x \times P_y$ compatible avec la résolution $(N_x, N_y)$ du maillage :

1. **Verrouillage cœur par cœur (*CPU pinning*) :**
   ```bash
   # OpenMPI : 1 rang MPI épinglé sur 1 cœur physique, pas de migration cache
   mpirun -np 8 --bind-to core --map-by core ./COVOF_2D-C++/bin/covof2d_cpp input.json
   ```
2. **Contraintes sur $P_x$ et $P_y$ :**
   - Chaque sous-domaine local doit posséder au moins $N_{x,\text{local}} \ge 2$ et $N_{y,\text{local}} \ge 2$ mailles (idéalement $\ge 8$ à $16$ mailles par côté pour limiter le ratio communication/calcul des halos et des marqueurs CVI).
   - Éviter un nombre premier $P$ élevé (ex. $P = 17$) sur un domaine carré : `MPI_Dims_create` produirait des tranches 1D très allongées $(17 \times 1)$ augmentant la surface d'échange et le nombre de marqueurs migrant entre rangs. Privilégier des facteurs équilibrés ($P \in \{2, 4, 6, 8, 12, 16, 32, 64\}$).

---

## Conteneurisation et Vérification Post-Installation

### Définition Singularity / Apptainer (`Apptainer.def`)

Sur supercalculateurs où la compilation native est délicate, une image **Apptainer/Singularity** (*rootless*) garantit la portabilité complète :

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
        python3 \
        python3-numpy \
        ca-certificates
    apt-get clean && rm -rf /var/lib/apt/lists/*

%environment
    export LC_ALL=C
    export OPENBLAS_NUM_THREADS=1

%runscript
    exec "$@"
```

Construction et exécution :

```bash
apptainer build covof2d.sif Apptainer.def
mpirun -np 4 apptainer exec covof2d.sif ./COVOF_2D-C++/bin/covof2d_cpp cas.json
```

### Protocole de Vérification Immédiate Post-Installation

Après compilation des deux branches dans `/Users/elhadjiibrahimathiam/Documents/PROJET_COVOF-2D`, exécuter les vérifications suivantes depuis la racine du projet :

1. **Test d'aide en ligne de commande (`--help`) :**
   ```bash
   ./COVOF_2D-C++/bin/covof2d_cpp --help
   ./COVOF_2D-FORTRAN/bin/covof2d_fortran --help
   ```

2. **Test de vérification analytique rapide (Laplacien 2D en séquentiel et MPI) :**
   ```bash
   ./COVOF_2D-C++/bin/covof2d_cpp VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_1_laplacian/cas_1_laplacian_32x32.json
   mpirun -np 4 ./COVOF_2D-FORTRAN/bin/covof2d_fortran VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_1_laplacian/cas_1_laplacian_32x32.json
   ```
   *(Adapter le nom du fichier `.json` présent dans le dossier du cas de validation selon la grille souhaitée).*

3. **Suite de vérification automatique FC-WY-CVI-Marker (C++ & Fortran, `np=1` et `np=4`) :**
   ```bash
   python3 Guides_Tutos/Run_FC-WY-CVI-Marker/verify_cvi_markers.py
   ```
   Ce script exécute automatiquement les 3 paliers de vérification mathématique et physique (divergence analytique nulle de l'interpolation CVI Variante B, équilibre capillaire de Young-Laplace à $32\times 32$ et $64\times 64$, et parité stricte `np=1` vs `np=4` sur le tourbillon de Rider-Kothe et la goutte statique) pour les deux binaires `covof2d_cpp` et `covof2d_fortran`.

Lorsque la simulation Navier-Stokes se termine normalement, la console affiche :
`>>> SIMULATION COVOF 2D TERMINEE AVEC SUCCES A t = ... s <<<`
suivi du chemin du répertoire de sortie contenant les fichiers `.pvtr`, `.vtr` et `_history.csv`.
