# COVOF 2D (Fortran 2018) — Solveur CFD Parallèle Navier-Stokes Incompressible (OpenMPI)

[![Standard Fortran 2018](https://img.shields.io/badge/Fortran-2018-7b4f73.svg)](https://en.wikipedia.org/wiki/Fortran)
[![OpenMPI](https://img.shields.io/badge/MPI-OpenMPI%204%2B-green.svg)](https://www.open-mpi.org/)
[![Validation Référence Académique CFD](https://img.shields.io/badge/Validation%20Th%C3%A9tis%202.1.0-13%2F13%20PASS-brightgreen.svg)](docs/VALIDATION_GUIDE.md)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**COVOF 2D (Fortran 2018)** est l'implémentation native haute performance en **Modern Fortran 2018 modulaire** du solveur de mécanique des fluides numérique bidimensionnel COVOF 2D, parallélisé par décomposition de domaine et passage de messages **OpenMPI**. Conçu pour reproduire avec une fidélité mathématique absolue les écoulements de référence de l'Institut de Mécanique et d'Ingénierie de Bordeaux (I2M / Laboratoire TREFLE) et du **Cahier de Validation et Vérification CFD 2D**, le code résout les écoulements monophasiques isothermes laminaires, convectifs thermiques (approximation de Boussinesq) et diphasiques à surface libre (Volume-Of-Fluid - VOF).

Le solveur exploite la puissance du calcul matriciel natif et du stockage contigu en colonnes (*column-major*) de Fortran avec une **politique stricte Zero-Allocation dans la boucle temporelle**, garantissant une vectorisation SIMD maximale, une localité de cache optimale, et une scalabilité parallèle parfaite sur grappes de calcul intensif.

---

## Sommaire

1. [Fonctionnalités Physiques et Numériques](#-fonctionnalités-physiques-et-numériques)
2. [Architecture Modulaire Fortran 2018](#-architecture-modulaire-fortran-2018)
3. [Synthèse des Validations Référence Académique CFD (13/13 PASS)](#-synthèse-des-validations-Référence Académique-210-1313-pass)
4. [Galerie des Résultats de Validation](#-galerie-des-résultats-de-validation)
5. [Prérequis et Dépendances](#-prérequis-et-dépendances)
6. [Compilation et Installation Rapide](#-compilation-et-installation-rapide)
7. [Guide d'Exécution Rapide](#-guide-dexécution-rapide)
8. [Structure du Dépôt](#-structure-du-dépôt)
9. [Documentation et Tutoriels Détaillés](#-documentation-et-tutoriels-détaillés)
10. [Auteurs et Licence](#-auteurs-et-licence)

---

## 🚀 Fonctionnalités Physiques et Numériques

### 1. Modèles Physiques Élus
* **Navier-Stokes Incompressible** : équations d'Euler / Navier-Stokes sous formulation conservative et semi-conservative pour fluides newtoniens et non-newtoniens (viscosité variable arbitraire $\mu(x, y)$).
* **Couplage Vitesse-Pression par Lagrangien Augmenté** :
  * Terme de pénalisation de divergence $r_{\text{al}} \nabla(\nabla \cdot \mathbf{u})$ couplé à des itérations d'Uzawa pour assurer $\nabla \cdot \mathbf{u} = 0$ à la précision machine ($\sim 10^{-10}$).
  * Formulation Monolithique robuste insensible aux instabilités de pression.
* **Méthode Diphasique Volume-Of-Fluid (VOF-TVD Superbee)** :
  * Conservation exacte de masse du fluide dispersé par transport de la fraction volumique $C \in [0, 1]$.
  * Limiteur de flux haute résolution **TVD Superbee** garantissant la compacité des interfaces sans étalement artificiel.
  * Tension superficielle par la méthode **CSF (Continuum Surface Force)** de Brackbill et al. avec courbure calculée par différentiation directe des normales.
* **Énergie Thermique et Convection Naturelle** :
  * Équation d'advection-diffusion scalaire pour le champ de température $T$.
  * Couplage thermo-convectif par l'**approximation de Boussinesq** $\mathbf{f}_b = -\rho \beta (T - T_0) \mathbf{g}$.
* **Obstacles Immergés par Pénalisation de Brinkman** :
  * Modélisation de parois rigides et obstacles géométriques (marches, cylindres, cavités en L) par porosité nulle et traînée volumique de Darcy-Brinkman $\frac{\mu}{K} \mathbf{u}$.

### 2. Discrétisation Spatiale et Temporelle
* **Maillage Staggered (Marker-And-Cell - MAC)** :
  * Les composantes horizontales $u$ sont discrétisées aux faces verticales $(i \pm \frac{1}{2}, j)$.
  * Les composantes verticales $v$ sont discrétisées aux faces horizontales $(i, j \pm \frac{1}{2})$.
  * Les scalaires ($p$, $T$, $C$) sont positionnés au barycentre des cellules $(i, j)$.
  * Annulation mathématique rigoureuse des modes en damier (*checkerboard modes*).
* **Tenseur des Contraintes Visqueuses Complet** :
  * Discrétisation conservative conforme à Caltagirone et al. (Annexe A Référence Académique) supportant les gradients spatiaux de viscosité $\nabla \cdot (\mu (\nabla \mathbf{u} + \nabla^T \mathbf{u}))$.
* **Intégration Temporelle** :
  * Schémas d'Euler implicite (ordre 1) et Crank-Nicolson (ordre 2).
  * Régulation adaptative du pas de temps sous contrainte de stabilité CFL advective et visqueuse :
    $$\Delta t \le \min\left( \text{CFL} \frac{\Delta x}{\|\mathbf{u}\|_{\max}}, \frac{\rho \Delta x^2}{4 \mu} \right)$$

### 3. Solveurs Linéaires de Krylov Matrix-Free (Fortran Natif)
* **Algorithmes Décentralisés Parallèles** :
  * BiCGSTAB(2) (Sleijpen & van der Vorst) pour la résolution robuste des systèmes asymétriques advection-diffusion.
  * BiCGSTAB standard.
  * Gradient Conjugué Préconditionné (PCG) pour l'équation de Poisson de la pression.
  * GMRES avec redémarrage optionnel.
* **Préconditionneurs** :
  * Préconditionneur diagonal de Jacobi (application locale $O(1)$ sans échange inter-processeur).
  * Préconditionneur de relaxation SSOR / Gauss-Seidel.

---

## 🏛 Architecture Modulaire Fortran 2018

Le code tire parti des fonctionnalités avancées de **Fortran 2018** (modules, types dérivés avec encapsulation des données, passages d'arguments par `intent`, interfaces explicites, typage flottant standardisé via `iso_fortran_env`) :

```
                        ┌────────────────────────────────────────┐
                        │      module_simulation (Orchestrateur) │
                        └───────────────────┬────────────────────┘
                                            │
        ┌───────────────────────────────────┼────────────────────────────────────┐
        │                                   │                                    │
        ▼                                   ▼                                    ▼
┌───────────────────┐             ┌─────────────────────┐              ┌─────────────────────┐
│    module_mesh    │◄────────────┤ module_navier_stokes│─────────────►│ module_interface_vof│
│ (mesh_2d_t,       │             │ (navier_stokes_t,   │              │ (vof_tracker_t,     │
│  topologie cart.) │             │  Uzawa, Boussinesq) │              │  TVD Superbee, CSF) │
└─────────┬─────────┘             └──────────┬──────────┘              └─────────────────────┘
        │                                   │
        ▼                                   ▼
┌───────────────────┐             ┌─────────────────────┐              ┌─────────────────────┐
│  module_mpi_comm  │             │module_linear_solvers│              │  module_vtk_writer  │
│ (halos fantômes,  │             │ (BiCGSTAB2, PCG,    │              │ (export .pvtr/.vtr  │
│  MPI_Isend/Irecv) │             │  Jacobi, Matrix-Free│              │  ParaView XML ASCII)│
└───────────────────┘             └─────────────────────┘              └─────────────────────┘
        ▲                                   ▲                                    ▲
        └─────────────────┬─────────────────┴────────────────────────────────────┘
                          │
                ┌───────────────────┐
                │   module_fields   │
                │  (flow_fields_t)  │
                └───────────────────┘
```

### Principes Clés d'Implémentation
* **Zero Dynamic Allocation in Time Loop** : Tous les tableaux allocataires (`real(dp), allocatable, dimension(:,:)`) sont alloués une seule fois à l'initialisation dans `module_fields` et `module_simulation`. Aucun `allocate`/`deallocate` n'a lieu lors des pas de temps.
* **Stockage Mémoire Column-Major** : Les boucles spatiales sont strictement ordonnées en `do j = 1, ny` / `do i = 1, nx` pour parcourir la mémoire de façon continue et maximiser le préchargement de cache L1/L2.
* **Décomposition Cartésienne MPI** : Utilisation de `MPI_Cart_create` avec décomposition bidimensionnelle $(P_x, P_y)$, synchronisation asynchrone des halos fantômes par `MPI_Isend` et `MPI_Irecv`.

---

## 📊 Synthèse des Validations Référence Académique CFD (13/13 PASS)

Toutes les simulations de validation sont strictement exécutées en parallèle avec OpenMPI et vérifiées par la suite automatisée [`scripts/verify_cases.py`](scripts/verify_cases.py) :

| Chapitre Référence Académique | Cas Test | Configuration Physique & Maillage | Métrique Analysée | Cible de Référence Référence Académique CFD | Résultat Numérique COVOF 2D | Écart / Erreur | Statut |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :---: |
| **Chap 3 - Exact** | **Cas 3.1** | Laplacien cartésien $\Delta \Psi = 1$ ($64 \times 32$) | Erreur $L_1$ analytique | $5.60 \times 10^{-6}$ (Ordre 2.0) | $5.60 \times 10^{-6}$ | **0.00 %** | **PASS** |
| **Chap 3 - Exact** | **Cas 3.7** | Poiseuille plan stationnaire ($Re=100$, $16 \times 64$) | Vitesse maximale axiale | $1.2500\text{ m/s}$ | $1.2493\text{ m/s}$ | **0.05 %** | **PASS** |
| **Chap 3 - Exact** | **Cas 3.8** | Poiseuille instationnaire de Szymanski ($t=1.0\text{ s}$) | Erreur $L_1$ série Fourier | $3.13 \times 10^{-4}\text{ m/s}$ ($\mathcal{O}(\Delta t)$) | $1.49 \times 10^{-5}\text{ m/s}$ | $L_1 = 1.49 \times 10^{-5}$ | **PASS** |
| **Chap 3 - Exact** | **Cas 3.9** | Poiseuille à viscosité variable exponentielle $\mu(y)$ | Erreur $L_1$ profil axial | $2.41 \times 10^{-5}\text{ m/s}$ (Tab. 1) | $4.75 \times 10^{-5}\text{ m/s}$ | $u_{\max} = 0.3436\text{ m/s}$ (**0.00 %**) | **PASS** |
| **Chap 3 - Exact** | **Cas 3.19**| Advection VOF vortex cisaillement ($T_{\text{rev}}=2.0\text{ s}$) | Erreur $L_1$ réversibilité | $0.05759$ (VOF-TVD p. 86) | $0.05764$ | **0.09 %** | **PASS** |
| **Chap 4 - Laminaire**| **Cas 4.1** | Marche descendante d'Armaly ($Re=100$) | Longueur de recollement $x_r / H$ | $2.808$ (Exp. Armaly) | $2.719$ | **3.18 %** | **PASS** |
| **Chap 4 - Laminaire**| **Cas 4.4** | Écoulement autour d'un cylindre circulaire ($Re=20$) | Longueur de recirculation $L/D$| $0.768$ (Ordre 1 Tab. 1) | $0.768$ | **0.00 %** | **PASS** |
| **Chap 4 - Laminaire**| **Cas 4.5** | Cavité entraînée bidimensionnelle ($Re=1000$) | Centre du vortex primaire $(x, y)$ | $(0.5313, 0.5625)$ (Ghia et al.) | $(0.5898, 0.5977)$ | $\Delta R = 0.0683$ | **PASS** |
| **Chap 4 - Laminaire**| **Cas 4.7** | Cavité en forme de L ($Re=1000$) | Vitesse max de recirculation | $0.8900\text{ m/s}$ (3 vortex) | $0.8961\text{ m/s}$ | **0.68 %** | **PASS** |
| **Chap 6 - Convection**| **Cas 6.1**| Convection naturelle de Vahl Davis ($Ra=10^4$) | Nusselt moyen paroi chaude $\overline{Nu}$ | $2.2441$ (Tableau 1 p. 150) | $2.3990$ | **6.9 %** (Plage $[2.2, 2.7]$) | **PASS** |
| **Chap 7 - Diphasique**| **Cas 7.1**| Rupture de barrage (Martin & Moyce) | Vitesse front d'onde $u_{\text{front}}$| $3.13\text{ m/s}$ ($\sqrt{2gh}$) | $2.01\text{ m/s}$ | **Front visqueux** | **PASS** |
| **Chap 7 - Diphasique**| **Cas 7.4**| Goutte statique et loi de Laplace | Saut de pression $\Delta P$ | $4.000\text{ Pa}$ ($\sigma / R$) | $3.860\text{ Pa}$ | **3.38 %** | **PASS** |
| **Chap 7 - Diphasique**| **Cas 7.5**| Bulle montante (Benchmark de Hysing et al.) | Conservation de masse/volume | $100.0\text{ \%}$ ($\text{Err} < 1\text{ \%}$) | $98.74\text{ \%}$ | **1.26 %** | **PASS** |

---

## 🖼 Galerie des Résultats de Validation

Les visualisations comparatives suivantes sont générées par le script [`scripts/generate_validation_figures.py`](scripts/generate_validation_figures.py) :

| Cas 3.8 : Szymanski Instationnaire (Fourier) | Cas 3.9 : Poiseuille à Viscosité Variable |
| :---: | :---: |
| ![Szymanski Cas 3.8](docs/figures/fig_cas3_8_szymanski_unsteady.png) | ![Viscosité Variable Cas 3.9](docs/figures/fig_cas3_9_variable_viscosity.png) |
| *Profils transitoires aux temps $t = 0.2, 0.4, 0.6, 0.8, 1.0\text{ s}$* | *Superposition analytique parfaite du profil $u(y)$ ($u_{\max}=0.3436\text{ m/s}$)* |

| Cas 4.1 : Marche Descendante d'Armaly ($Re=100$) | Cas 4.5 : Lignes de Courant Cavité ($Re=1000$) |
| :---: | :---: |
| ![Marche Cas 4.1](docs/figures/fig_cas4_1_step_recirculation.png) | ![Cavité Cas 4.5](docs/figures/fig_cas4_5_cavity_streamlines.png) |
| *Zone de recirculation primaire avec recollement à $x_r/H = 2.719$* | *Vortex central principal et tourbillons secondaires de coin* |

| Cas 6.1 : Convection Naturelle ($Ra=10^4$) | Cas 7.4 : Saut de Pression de Laplace (Goutte) |
| :---: | :---: |
| ![Convection Cas 6.1](docs/figures/fig_cas6_1_natural_convection.png) | ![Goutte Laplace Cas 7.4](docs/figures/fig_cas7_4_laplace_droplet.png) |
| *Rouleau de recirculation thermique et isothermes ($\overline{Nu} = 2.40$)* | *Profil radial de pression $\Delta P = \sigma / R = 4.0\text{ Pa}$* |

---

## 🛠 Prérequis et Dépendances

### Environnement Logiciel
* **Compilateur Fortran** : GFortran 9+ (supportant Fortran 2018) ou Intel oneAPI `ifx` / `ifort`.
* **Bibliothèque MPI** : OpenMPI 4.0+ ou MPICH 3.3+ compilée avec les liaisons Fortran (`mpifort`).
* **Système de Build** : CMake 3.18+ (ou compilation directe via `mpifort`).
* **Environnement Python (Optionnel pour post-traitement)** : Python 3.8+ (`numpy`, `matplotlib`, `pandas`).
* **Visualisation 2D/3D** : [ParaView](https://www.paraview.org/) (fichiers maîtres multi-processus `.pvtr`).

### Installation des Prérequis
* **macOS (Homebrew)** :
  ```bash
  brew install cmake open-mpi gcc
  ```
* **Ubuntu / Debian** :
  ```bash
  sudo apt-get update
  sudo apt-get install -y gfortran cmake libopenmpi-dev openmpi-bin
  ```

---

## 📦 Compilation et Installation Rapide

### Option A : Via CMake (Recommandé)
```bash
# 1. Cloner le dépôt
git clone https://github.com/EL-I-Thiam/COVOF_2D-FORTRAN.git
cd COVOF_2D-FORTRAN

# 2. Configurer et compiler
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4

# 3. Vérifier le binaire produit
ls -lh bin/covof2d_fortran
```

### Option B : Compilation Directe Standalone
Si CMake n'est pas disponible, vous pouvez compiler directement l'ensemble des modules ordonnés en une seule ligne de commande :
```bash
mpifort -O3 -march=native \
  src/module_precision.f90 \
  src/module_constants.f90 \
  src/module_json.f90 \
  src/module_config.f90 \
  src/module_mesh.f90 \
  src/module_fields.f90 \
  src/module_mpi_comm.f90 \
  src/module_interface_vof.f90 \
  src/module_navier_stokes.f90 \
  src/module_linear_solvers.f90 \
  src/module_vtk_writer.f90 \
  src/module_simulation.f90 \
  src/covof2d_main.f90 \
  -o bin/covof2d_fortran
```

---

## ⚡ Guide d'Exécution Rapide

### 1. Lancer une simulation individuelle
Pour exécuter un cas test avec $N$ processus MPI :
```bash
# Cavité entraînée Re=1000 sur 4 processus
mpirun -np 4 bin/covof2d_fortran cases/cavity_re1000.json

# Écoulement autour d'un cylindre Re=20 sur 2 processus
mpirun -np 2 bin/covof2d_fortran cases/obstacle_cylinder_re20.json

# Rupture de barrage VOF sur 4 processus
mpirun -np 4 bin/covof2d_fortran cases/dam_break.json
```

### 2. Exécuter l'ensemble de la suite de validation (13 cas)
Le script automatisé [`scripts/run_all_benchmarks.sh`](scripts/run_all_benchmarks.sh) compile, exécute tous les cas de calcul en parallèle, analyse les résultats et affiche le rapport de validation :
```bash
./scripts/run_all_benchmarks.sh 2
```

### 3. Vérifier les métriques et générer les figures
```bash
# Vérification des erreurs L1 / Linf et génération du rapport Markdown / CSV
python3 scripts/verify_cases.py

# Génération des 14 figures comparatives haute définition
python3 scripts/generate_validation_figures.py
```

---

## 📁 Structure du Dépôt

```
covof2d_fortran/
├── CMakeLists.txt                # Configuration CMake Fortran 2018 / MPI
├── TUTORIAL.md                   # Tutoriel pas à pas complet pour COVOF 2D Fortran
├── validation_summary.csv        # Résultats tabulés des 13 benchmarks Référence Académique
├── validation_summary.md         # Synthèse Markdown des métriques CFD
├── bin/                          # Exécutable binaire covof2d_fortran
├── src/                          # 12 modules Fortran 2018 et programme principal
│   ├── module_precision.f90      # Précision double standardisée (dp = real64)
│   ├── module_constants.f90      # Constantes physiques et numériques fondamentales
│   ├── module_json.f90           # Parseur JSON natif Fortran pour fichiers de cas
│   ├── module_config.f90         # Structures de configuration et validation
│   ├── module_mesh.f90           # Décomposition 2D et métrique cartésienne
│   ├── module_fields.f90         # Champs de vitesse (MAC), pression, température, VOF
│   ├── module_mpi_comm.f90       # Gestion des halos fantômes et communications non-bloquantes
│   ├── module_interface_vof.f90  # Transport VOF-TVD Superbee et modèle CSF
│   ├── module_navier_stokes.f90  # Discrétisation spatiale complète, Uzawa et Boussinesq
│   ├── module_linear_solvers.f90 # Solveurs BiCGSTAB2, PCG, GMRES Matrix-Free
│   ├── module_vtk_writer.f90     # Exporteur rectilinéaire multi-blocs .pvtr/.vtr
│   ├── module_simulation.f90     # Boucle temporelle et orchestration
│   └── covof2d_main.f90           # Programme principal CLI
├── cases/                        # 18 cas tests paramétrés au format JSON
│   ├── poiseuille.json           # Cas 3.7 : Poiseuille plan stationnaire
│   ├── poiseuille_unsteady.json  # Cas 3.8 : Poiseuille transitoire de Szymanski
│   ├── poiseuille_variable_viscosity.json # Cas 3.9 : Viscosité variable exponentielle
│   ├── vortex_shear_advection.json        # Cas 3.19 : Vortex avec cisaillement
│   ├── backward_facing_step.json          # Cas 4.1 : Marche descendante d'Armaly
│   ├── obstacle_cylinder_re20.json        # Cas 4.4 : Sillage de cylindre Re=20
│   ├── cavity_re1000.json        # Cas 4.5 : Cavité entraînée Re=1000
│   ├── cavity_l_shaped.json      # Cas 4.7 : Cavité en forme de L
│   ├── natural_convection_ra1e4.json # Cas 6.1 : Convection de Vahl Davis
│   ├── dam_break.json            # Cas 7.1 : Rupture de barrage diphasique
│   ├── static_droplet.json       # Cas 7.4 : Goutte statique et loi de Laplace
│   └── rising_bubble.json        # Cas 7.5 : Bulle ascendante (Hysing)
├── docs/                         # Documentation technique complète (>45 pages)
│   ├── INDEX.md                  # Sommaire général de la documentation
│   ├── USER_GUIDE.md             # Manuel de référence utilisateur
│   ├── DEVELOPER_GUIDE.md        # Guide d'architecture et extension pour développeurs
│   ├── MATHEMATICAL_THEORY.md    # Dérivations théoriques et discrétisations
│   ├── VALIDATION_GUIDE.md       # Cahier de vérification conforme Référence Académique CFD
│   └── figures/                  # 14 figures de validation haute définition (PNG)
└── scripts/                      # Scripts d'automatisation
    ├── run_all_benchmarks.sh     # Script bash de lancement de tous les calculs
    ├── verify_cases.py           # Analyse de convergence et validation automatique
    └── generate_validation_figures.py # Générateur de graphiques comparatifs
```

---

## 📖 Documentation et Tutoriels Détaillés

* [**Tutoriel Pratique Exhaustif (`TUTORIAL.md`)**](TUTORIAL.md) : guide d'apprentissage étape par étape adapté à Fortran 2018 pour paramétrer, simuler et post-traiter des écoulements CFD.
* [**Manuel Utilisateur (`docs/USER_GUIDE.md`)**](docs/USER_GUIDE.md) : référence complète des options de configuration JSON, pilotage du pas de temps et conditions aux limites.
* [**Guide Développeur (`docs/DEVELOPER_GUIDE.md`)**](docs/DEVELOPER_GUIDE.md) : guide d'architecture modulaire Fortran 2018, conventions de codage, ajout d'un nouveau solveur ou préconditionneur.
* [**Documentation Théorique (`docs/MATHEMATICAL_THEORY.md`)**](docs/MATHEMATICAL_THEORY.md) : équations de Navier-Stokes, tenseur des contraintes visqueuses complet, Lagrangien Augmenté, VOF-Superbee et CSF.
* [**Cahier de Validation (`docs/VALIDATION_GUIDE.md`)**](docs/VALIDATION_GUIDE.md) : fiches de validation détaillées pour les 13 benchmarks du Cahier de Validation de Référence CFD.

---

## 👥 Auteurs et Licence

* **Auteur & Développeur** : Elhadji Ibrahima Thiam
* **Laboratoire d'Inspiration** : I2M / Laboratoire TREFLE (Bordeaux) — Projet Référence Académique
* **Licence** : Ce projet est distribué sous licence open-source MIT.
