# Standards de Codage, Profilage et Guide d'Extension

Pour pérenniser le co-développement gémellaire C++17 / Fortran 2018 et garantir une maintenabilité à long terme sur des décennies, COVOF 2D s'appuie sur des standards d'ingénierie logicielle rigoureux, une méthodologie d'analyse de performance fine et un protocole formel pour l'intégration de nouvelles fonctionnalités.

Ce chapitre formalise les règles de style de code, les outils de profilage HPC et guide le développeur pas à pas dans l'extension du solveur à travers un exemple pratique complet.

---

## 7.1. Standards de Codage et Conventions d'Ingénierie

### 7.1.1. Conventions C++17

Le code C++ respecte les directives du standard moderne C++17 et la norme *C++ Core Guidelines* de Bjarne Stroustrup et Herb Sutter :

1. **Formatage automatique :** Tout code doit être formaté via l'outil `clang-format` avec la configuration définie dans le fichier `.clang-format` racine (style basé sur Google C++ avec indentation à 4 espaces, limite de ligne à 100 colonnes).
2. **Gestion RAII stricte (*Resource Acquisition Is Initialization*) :**
   - Aucune gestion manuelle de mémoire via `new` ou `delete`.
   - Utilisation systématique de conteneurs standards (`std::vector`) ou de pointeurs intelligents à propriété exclusive (`std::unique_ptr`) et partagée (`std::shared_ptr`).
3. **Const-correctness et Noexcept :**
   - Toute méthode ne modifiant pas l'état interne de la classe doit être qualifiée `const`.
   - Les accesseurs triviaux et les destructeurs doivent être marqués `noexcept` pour assister le compilateur dans l'optimisation des sauts.
4. **Bannissement des macros pour les constantes :**
   - Remplacement de `#define PI 3.14159...` par `constexpr double PI = 3.14159265358979323846;`.
5. **Absence totale de variables globales mutables :**
   - Toutes les grandeurs d'état résident dans les instances explicites de `Simulation`, `Mesh2D` et `FlowFields`.

### 7.1.2. Conventions Fortran 2018

Le code Fortran adhère aux meilleures pratiques du calcul scientifique moderne :

1. **Directive `implicit none` obligatoire :**
   - Chaque module et chaque sous-programme doit débuter par `implicit none` sans exception.
2. **Encapsulation modulaire `private` par défaut :**
   - Tout module débute par la clause `private`, suivie de la liste explicite des symboles publics autorisés (`public :: ns_apply_system_matrix, ...`).
3. **Spécificateurs d'intention exhaustifs :**
   - Chaque argument formel de fonction ou sous-routine doit impérativement déclarer son intention :
     - `intent(in)` : Argument lu uniquement.
     - `intent(out)` : Argument initialisé et écrit par la routine.
     - `intent(inout)` : Argument lu et modifié.
4. **Interdiction des idiomes obsolètes :**
   - Bannissement formel des blocs `COMMON`, des instructions `EQUIVALENCE`, des `GOTO` calculés, et de l'attribut `SAVE` implicite.
5. **Types dérivés et procédures liées :**
   - Modélisation objet par types dérivés extensibles (`type, public :: mesh_2d_t`) avec interfaces génériques et surcharge d'opérateurs si pertinent.

---

## 7.2. Profilage HPC et Optimisation des Performances

Sur supercalculateur, l'optimisation ne procède jamais par intuition subjective mais par métrologie instrumentée rigoureuse.

### 7.2.1. Profilage Matériel sous Linux via `perf`

L'outil standard `perf` permet d'inspecter les compteurs de performance matérielle du processeur (*Hardware Performance Counters*) sans surcoût d'instrumentation :

```bash
# Mesure globale des compteurs matériels
perf stat -e cycles,instructions,cache-misses,L1-dcache-load-misses,LLC-load-misses \
     ./bin/covof_cpp params.json

# Analyse du nombre d'instructions par cycle (IPC)
# Un IPC < 1.0 indique un engorgement mémoire (Memory Stalls)
# Un IPC > 2.0 indique une excellente saturation des unités d'exécution vectorielles.
```

Pour localiser précisément les fonctions les plus consommatrices (*hot spots*) :
```bash
perf record -g --call-graph dwarf ./bin/covof_cpp params.json
perf report --stdio
```

### 7.2.2. Analyse Approfondie avec Intel VTune Profiler

Pour les architectures x86 Intel Xeon (Cascade Lake, Sapphire Rapids, Emerald Rapids), Intel VTune fournit une cartographie complète :
1. **Memory Access Analysis :** Mesure la bande passante consommée sur les canaux DDR5 / HBM et identifie les accès mémoire non séquentiels (*strided access*).
2. **Roofline Model :** Positionne les boucles de calcul critiques (diffusion visqueuse, projection de pression) sur le graphique intensité arithmétique (FLOPs/octet) versus performance crête (GFLOPs/s). Elle détermine immédiatement si l'algorithme est limité par la mémoire (*bandwidth bound*) ou par la puissance de calcul (*compute bound*).

### 7.2.3. Profilage des Communications MPI Distribuées avec Scalasca et Score-P

Pour traquer les goulets d'étranglement de scalabilité sur des centaines de cœurs :
```bash
# Compilation instrumentée Score-P
scorep mpicxx -O3 -std=c++17 src/*.cpp -o covof_scorep

# Exécution avec traçage OTF2
export SCOREP_EXPERIMENT_DIRECTORY=scorep_trace
mpirun -np 64 ./covof_scorep params.json

# Visualisation des déséquilibres avec Cube GUI
cube scorep_trace/profile.cubex
```
Ce profilage identifie immédiatement le temps d'attente parasite (*wait time*) induit par un déséquilibre de charge (*load imbalance*) dans les appels non bloquants `MPI_Waitall` des mailles fantômes.

---

## 7.3. Guide Développeur Pas-à-Pas : Ajout d'une Nouvelle Fonctionnalité

Pour illustrer de manière concrète et opérationnelle la procédure d'extension du solveur, nous détaillons ci-après l'implémentation complète de bout en bout de l'**équation de transport d'un scalaire passif $S(\mathbf{x}, t)$** (représentant par exemple la concentration d'un traceur chimique ou d'un polluant) :
$$
\frac{\partial S}{\partial t} + \nabla \cdot (\mathbf{u} S) = \nabla \cdot (D_s \nabla S)
$$

### Étape 1 : Ajout de la Clé de Configuration dans le Dictionnaire JSON

1. **Dans `include/covof/Config.hpp` :**
   ```cpp
   struct Config {
       // ... paramètres existants ...
       bool   enable_scalar_transport = false;
       double scalar_diffusivity      = 1.0e-5; // m^2/s
   };
   ```

2. **Dans `src/Config.cpp` :**
   ```cpp
   if (root.contains("enable_scalar_transport")) {
       cfg.enable_scalar_transport = root["enable_scalar_transport"].get<bool>();
   }
   if (root.contains("scalar_diffusivity")) {
       cfg.scalar_diffusivity = root["scalar_diffusivity"].get<double>();
   }
   ```

3. **Côté Fortran dans `COVOF_2D-FORTRAN/src/module_config.f90` :**
   Ajouter les champs dans le type dérivé `config_t` et leur lecture dans la routine d'analyse JSON.

### Étape 2 : Déclaration et Allocation du Champ $S$

1. **Dans `include/covof/Field.hpp` :**
   ```cpp
   struct FlowFields {
       // ... champs existants ...
       std::vector<double> S;     // Scalaire passif aux centres de mailles
       std::vector<double> S_old; // Étage temporel n
       
       void allocate(int nxt, int nyt) {
           // ...
           S.assign(nxt * nyt, 0.0);
           S_old.assign(nxt * nyt, 0.0);
       }
   };
   ```

2. **Dans `COVOF_2D-FORTRAN/src/module_flow_fields.f90` :**
   Ajouter `real(dp), allocatable :: S(:, :)` et allouer le tableau dans `fields_allocate_2d`.

### Étape 3 : Implémentation du Noyau Numérique d'Advection-Diffusion

1. **Dans `include/covof/NavierStokesSolver.hpp` et `src/NavierStokesSolver.cpp` :**
   Ajouter la méthode `advancePassiveScalar` :
   ```cpp
   void NavierStokesSolver::advancePassiveScalar(FlowFields& f, double dt) {
       if (!config_.enable_scalar_transport) return;
       
       const int nxt = mesh_.nx_total();
       const double idx = 1.0 / mesh_.dx();
       const double idy = 1.0 / mesh_.dy();
       const double Ds  = config_.scalar_diffusivity;
       
       // Boucle de mise à jour explicite du scalaire
       #pragma omp parallel for collapse(2)
       for (int j = 1; j <= mesh_.ny_local(); ++j) {
           for (int i = 1; i <= mesh_.nx_local(); ++i) {
               // Flux diffusifs de Fick aux faces
               double diff_x = Ds * ((f.S[IDX(i+1, j)] - f.S[IDX(i, j)]) * idx -
                                     (f.S[IDX(i, j)] - f.S[IDX(i-1, j)]) * idx) * idx;
               double diff_y = Ds * ((f.S[IDX(i, j+1)] - f.S[IDX(i, j)]) * idy -
                                     (f.S[IDX(i, j)] - f.S[IDX(i, j-1)]) * idy) * idy;
               
               // Flux convectifs Upwind
               double u_e = f.u[IDX_U(i, j)];
               double u_w = f.u[IDX_U(i-1, j)];
               double v_n = f.v[IDX_V(i, j)];
               double v_s = f.v[IDX_V(i, j-1)];
               
               double flux_e = (u_e > 0.0) ? u_e * f.S[IDX(i, j)] : u_e * f.S[IDX(i+1, j)];
               double flux_w = (u_w > 0.0) ? u_w * f.S[IDX(i-1, j)] : u_w * f.S[IDX(i, j)];
               double flux_n = (v_n > 0.0) ? v_n * f.S[IDX(i, j)] : v_n * f.S[IDX(i, j+1)];
               double flux_s = (v_s > 0.0) ? v_s * f.S[IDX(i, j-1)] : v_s * f.S[IDX(i, j)];
               
               double conv = (flux_e - flux_w) * idx + (flux_n - flux_s) * idy;
               
               // Avancement Euler
               f.S[IDX(i, j)] = f.S_old[IDX(i, j)] + dt * (diff_x + diff_y - conv);
           }
       }
   }
   ```

2. **Équivalent Fortran 2018 dans `module_navier_stokes.f90` :**
   Créer la sous-routine `subroutine ns_advance_passive_scalar(fields, mesh, cfg, dt)`.

### Étape 4 : Mise à Jour des Mailles Fantômes MPI

Dans la méthode `Simulation::timeStep()`, intercaler la synchronisation MPI du nouveau champ :
```cpp
mpi_comm_.exchangeHalos(fields_.S, mesh_.nx_total(), mesh_.ny_total(), mesh_.nghost());
```

### Étape 5 : Exportation dans les Fichiers VTK

Dans `src/VTKWriter.cpp`, ajouter la section de données scalaires `S` dans l'arbre XML `<CellData>` :
```cpp
if (config.enable_scalar_transport) {
    writeScalarField(vtr_file, "passive_scalar", fields.S);
}
```

### Étape 6 : Validation par Solution Analytique MMS

Créer un cas de test unitaire vérifiant la décroissance exponentielle d'un mode de Fourier sinusoïdal en diffusion pure ($u = v = 0$) :
$$S_{\text{anal}}(x, y, t) = S_0 \sin\left(\frac{2\pi x}{L}\right) \sin\left(\frac{2\pi y}{H}\right) e^{-D_s \left( \frac{4\pi^2}{L^2} + \frac{4\pi^2}{H^2} \right) t}$$
Vérifier que l'erreur relative entre la solution numérique et la solution analytique décroît au taux théorique attendu de $O(\Delta x^2 + \Delta t)$.
