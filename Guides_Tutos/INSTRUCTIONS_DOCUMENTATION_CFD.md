# Cahier des Charges et Instructions de Génération de la Documentation Complète d'un Code CFD Dual C++ / Fortran

## 0. Objectif et Périmètre de la Mission

Ce document constitue le **protocole directeur** pour l'analyse intégrale du dépôt source et la rédaction d'une suite documentaire de niveau industriel et académique pour un solveur de mécanique des fluides numérique (CFD) écrit en **architecture hybride C++ et Fortran**.

L'objectif est de produire **trois volumes distincts, autonomes et exhaustifs**, disponibles simultanément sous forme de **fichiers Markdown modulaires (`.md`)** et compilables en **manuels PDF haute typographie (`.pdf`)** d'au moins **40 pages chacun** (soit un minimum de 120 pages cumulées hors annexes générées) :

1. **Volume I — Guide Utilisateur (`user_guide.md` / `Guide_Utilisateur_CFD.pdf`)** : $\ge 40$ pages.
2. **Volume II — Guide de Développement (`dev_guide.md` / `Guide_Developpement_CFD.pdf`)** : $\ge 40$ pages.
3. **Volume III — Recueil de Tutoriels et Cas de Validation (`tutorials_guide.md` / `Recueil_Tutoriels_CFD.pdf`)** : $\ge 40$ pages.

---

## 1. Règles d'Or et Contraintes d'Exhaustivité

Toute rédaction effectuée selon ces instructions doit impérativement respecter les six règles fondamentales suivantes :

### Règle 1 : Principe du « Zéro Impasse » (Exhaustivité absolue)
* **Aucun élément du code ne doit être omis.** Chaque répertoire, chaque fichier source (`.cpp`, `.hpp`, `.h`, `.cxx`, `.f90`, `.F90`, `.f03`, `.f08`, `.cuh`), chaque classe C++, chaque module Fortran, chaque modèle physique, chaque schéma numérique et chaque paramètre d'entrée doit être recensé, expliqué et documenté.
* **Interdiction formelle des ellipses :** Les formulations telles que `// ... reste du code`, `etc.`, `autres paramètres similaires` ou les descriptions génériques d'une seule ligne sont strictement proscrites. Chaque mécanisme doit être décrit dans son intégralité.

### Règle 2 : Calibrage Volumétrique ($\ge 40$ pages par volume)
* Pour atteindre et dépasser le seuil de **40 pages au format A4** (police 10 ou 11 pt, interligne 1.15, marges standard 2.5 cm) par volume :
  * Chaque volume doit contenir entre **18 000 et 25 000 mots**, accompagnés d'équations mathématiques complètes, de tableaux de référence exhaustifs, d'extraits de code source commentés et de diagrammes d'architecture (Mermaid / TikZ).
  * Ne jamais compresser l'information : développer systématiquement le **pourquoi physique**, la **formulation mathématique continue**, la **discrétisation numérique**, l'**implémentation C++/Fortran** et les **précautions d'usage**.

### Règle 3 : Rigorisme Mathématique et Typographique
* Toutes les équations, variables scalaires, tenseurs, opérateurs différentiels et nombres adimensionnels doivent être rédigés exclusivement en syntaxe LaTeX :
  * Formules en ligne : `$ ... $` (ex. : nombre de Reynolds $Re = \frac{\rho U L}{\mu}$, pas de temps $\Delta t$).
  * Formules hors-texte numérotées ou centrées : `$$ ... $$`.
* Chaque symbole introduit dans une équation doit être défini immédiatement après (signification physique, unités SI, et nom de la variable correspondante dans le code C++ et/ou Fortran).

### Règle 4 : Traitement Spécifique de l'Architecture Duale C++ / Fortran
Le code reposant sur une cohabitation entre C++ et Fortran, la documentation doit systématiquement expliciter la frontière entre les deux langages :
* **Rôle de chaque langage :** par exemple, C++ pour l'orchestration orientée objet, la gestion des entrées/sorties (JSON/YAML/HDF5), les structures de maillage non structuré ou le pilotage MPI/GPU ; Fortran pour les noyaux de calcul intensif sur tableaux multidimensionnels, les opérateurs volumes finis/différences finies et les solveurs locaux.
* **Interopérabilité `iso_c_binding` :** documentation systématique des interfaces `extern "C"` côté C++ et `bind(C)` côté Fortran, de la correspondance des types (`double` $\leftrightarrow$ `real(c_double)`, `int` $\leftrightarrow$ `integer(c_int)`), de la gestion de l'ordre mémoire (**row-major** en C++ vs **column-major** en Fortran) et du passage des communicateurs MPI (`MPI_Comm_c2f` / `MPI_Comm_f2c`).

### Règle 5 : Double Livrable `.md` et `.pdf`
* Les fichiers `.md` doivent utiliser la syntaxe **MyST Markdown / Pandoc Markdown** enrichie (blocs d'admonition `:::note`, `:::warning`, tableaux Markdown étendus, blocs de code typés `cpp`, `fortran`, `cmake`, `bash`, `python`, `yaml`, diagrammes `mermaid`).
* Les fichiers `.md` doivent compiler sans aucune erreur vers un fichier `.pdf` professionnel via **Pandoc + XeLaTeX** ou **Sphinx + LaTeX** (avec table des matières automatique, liste des tableaux, liste des figures, en-têtes/pieds de page et numérotation des sections).

### Règle 6 : Protocole d'Inspection Préalable du Dépôt
Avant de rédiger chaque chapitre, l'auteur (ou l'agent IA) doit scanner systématiquement l'arborescence réelle du dépôt :
1. Lister tous les fichiers `CMakeLists.txt` et `Makefile` pour extraire les options de compilation, les macros préprocesseur (`#ifdef`, `#if defined`), et les bibliothèques liées.
2. Extraire tous les fichiers d'analyse d'entrées (parsers JSON/YAML/XML/Namelists) pour dresser la liste exhaustive à 100 % des mots-clés de configuration.
3. Cartographier toutes les déclarations `extern "C"` et `bind(C)` pour reconstituer la table complète des appels inter-langages.

---

## 2. Cahier des Charges Détaillé : Volume I — Guide Utilisateur ($\ge 40$ pages)

Le Guide Utilisateur s'adresse aux ingénieurs d'études, chercheurs et doctorants qui utilisent le code pour configurer, exécuter et analyser des simulations CFD sans modifier le noyau source.

### Structure obligatoire et pagination cible (Total : 42 à 50 pages)

#### Chapitre 1 : Introduction Générale et Capacités du Solveur (3 à 4 pages)
* **1.1. Présentation et philosophie du code :** historique, domaines d'application (aérodynamique, hydrodynamique, thermique, écoulements multiphasiques, milieux poreux, etc.).
* **1.2. Matrice complète des capacités physiques et numériques :** tableau récapitulatif exhaustif croisant dimensions (1D, 2D, 2D axisymétrique, 3D), régimes (stationnaire, instationnaire, incompressible, faiblement compressible, compressible), topologies de maillage et modèles physiques disponibles.
* **1.3. Architecture fonctionnelle vue par l'utilisateur :** flux de travail complet (Pré-traitement $\rightarrow$ Configuration $\rightarrow$ Exécution parallèle $\rightarrow$ Post-traitement).

#### Chapitre 2 : Installation, Compilation et Déploiement HPC (5 à 6 pages)
* **2.1. Prérequis logiciels et matrice de compatibilité :** compilateurs C++ (GCC, Clang, Intel `icpx`) et Fortran (`gfortran`, Intel `ifx`), standards requis (C++17/20, Fortran 2008/2018), MPI (OpenMPI, MPICH, Intel MPI), bibliothèques mathématiques (PETSc, HYPRE, Trilinos, LAPACK/BLAS), E/S (HDF5, CGNS, NetCDF, YAML-cpp).
* **2.2. Compilation pas à pas avec CMake :** explication détaillée de chaque option CMake (`-DENABLE_MPI=ON`, `-DENABLE_OPENMP=ON`, `-DENABLE_CUDA=OFF`, `-DCMAKE_BUILD_TYPE=Release`, drapeaux d'optimisation `-O3 -march=native`).
* **2.3. Déploiement sur supercalculateurs (Clusters HPC) :** scripts de soumission complets pour **Slurm** et **PBS/Torque**, gestion de l'affinité des processus MPI et des threads OpenMP (`OMP_NUM_THREADS`, `mpirun --bind-to core`).
* **2.4. Conteneurisation et vérification de l'installation :** utilisation via Docker / Singularity (Apptainer) et exécution de la suite de tests de validation rapide (`ctest`).

#### Chapitre 3 : Formulation Mathématique et Modèles Physiques (10 à 12 pages)
* **3.1. Système d'équations gouvernantes :** écriture complète sous forme conservative et non conservative (conservation de la masse, quantité de mouvement, énergie totale/interne, transport d'espèces). Exemple :
  $$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0$$
  $$\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) = -\nabla p + \nabla \cdot \boldsymbol{\tau} + \rho \mathbf{g} + \mathbf{f}_{\sigma}$$
* **3.2. Lois constitutives et thermodynamiques :** tenseur des contraintes visqueuses (fluides newtoniens et non newtoniens : loi puissance, Carreau, Bingham), lois d'état (gaz parfaits, gaz réels, fluides incompressibles à masse volumique variable), lois de viscosité et conductivité (Sutherland, polynômes).
* **3.3. Modèles de turbulence (RANS, URANS, LES, DES) :** formulation complète de chaque modèle implémenté ($k$-$\epsilon$, $k$-$\omega$ SST, Spalart-Allmaras, Smagorinsky, WALE), constantes de fermeture, lois de paroi (*wall functions*) et critères sur $y^+$.
* **3.4. Modèles multiphasiques et interfaces (si présents dans le code) :** méthodes VOF (*Volume of Fluid*), Level-Set, tension superficielle (modèle CSF de Brackbill), modèles diphasiques Euler-Euler ou Euler-Lagrange.
* **3.5. Transferts thermiques et couplages :** convection naturelle (approximation de Boussinesq), rayonnement, transfert thermique conjugué fluide-solide.
* **3.6. Table générale de nomenclature :** tableau recensant l'intégralité des symboles mathématiques, leurs unités SI, leur signification physique et leur identifiant dans le code.

#### Chapitre 4 : Méthodes Numériques (6 à 8 pages)
* **4.1. Discrétisation spatiale :** formulation Volumes Finis (centrée aux mailles / colocalisée / maillage décalé MAC) ou Éléments Finis ; calcul des gradients (Green-Gauss, moindres carrés pondérés).
* **4.2. Schémas d'advection et limiteurs de flux :** propriétés de conservation, monotonie (TVD), diffusion numérique et dispersion pour chaque schéma disponible (Upwind 1er ordre, Centré 2e ordre, QUICK, MUSCL, Van Leer, Superbee, Minmod, WENO).
* **4.3. Discrétisation temporelle et critère de stabilité :** schémas explicites (Euler, Runge-Kutta TVD d'ordre 2 et 3), implicites (Euler arrière, Crank-Nicolson, BDF2), calcul automatique du pas de temps adaptatif via les nombres $CFL$ convectif et diffusif.
* **4.4. Algorithmes de couplage vitesse-pression :** méthodes de projection (Chorin-Temam, projection incrémentale), algorithmes itératifs (SIMPLE, SIMPLEC, PISO, PIMPLE), interpolation de Rhie-Chow sur grilles colocalisées.
* **4.5. Solveurs linéaires et préconditionneurs :** méthodes de Krylov (CG, BiCGSTAB, GMRES), multigrille algébrique et géométrique (AMG, GAMG), critères d'arrêt absolus et relatifs.

#### Chapitre 5 : Géométrie, Maillages et Conditions aux Limites (5 à 6 pages)
* **5.1. Génération et formats de maillages supportés :** maillages cartésiens internes, import depuis Gmsh (`.msh`), CGNS (`.cgns`), OpenFOAM, MED ; métriques de qualité requises (non-orthogonalité, asymétrie *skewness*, rapport d'aspect).
* **5.2. Catalogue exhaustif des conditions aux limites (CL) :** pour chaque CL du code, fournir :
  * Le mot-clé exact dans le fichier d'entrée.
  * La formulation mathématique sur la frontière $\partial \Omega$ (Dirichlet, Neumann, Robin).
  * L'implémentation discrète via les cellules fantômes (*ghost cells*) ou flux de bord.
  * Les paramètres associés et cas d'usage typiques.

#### Chapitre 6 : Dictionnaire Encyclopédique des Paramètres d'Entrée (8 à 10 pages)
* Présentation exhaustive de **100 % des paramètres configurables** du code, classés par bloc/section logique du fichier d'entrée.
* Pour chaque bloc, fournir un tableau détaillé contenant :
  * **Nom exact de la clé** (sensible à la casse)
  * **Type informatique** (`integer`, `double`, `boolean`, `string`, `array<double>`)
  * **Valeur par défaut** (ou mention `OBLIGATOIRE`)
  * **Plage de valeurs autorisées / Options valides**
  * **Description physique et numérique détaillée** (impact sur la précision, la stabilité et le coût CPU).
* Fournir en fin de chapitre un **fichier de configuration maître commenté** regroupant toutes les sections.

#### Chapitre 7 : Exécution, Post-Traitement et Dépannage (5 à 6 pages)
* **7.1. Lancement d'une simulation et gestion des reprises (*Restart*) :** arguments en ligne de commande, signaux d'arrêt propre, fréquence de sauvegarde des points de contrôle (*checkpoints*).
* **7.2. Interprétation des journaux de convergence (*Logs*) :** explication colonne par colonne des sorties console et des fichiers d'historique de résidus ($L_1$, $L_2$, $L_\infty$), bilan de masse global et monitoring de sondes ponctuelles.
* **7.3. Visualisation et extraction de données :** structure des fichiers VTK/PVD, HDF5/XDMF, exploitation sous ParaView et scripts Python (`pyvista`, `matplotlib`, `h5py`).
* **7.4. Guide complet de dépannage (*Troubleshooting*) :** matrice symptômes $\rightarrow$ causes physiques/numériques $\rightarrow$ actions correctives (divergence de la pression, oscillation en damier, perte de conservation de la masse, blocage MPI).

---

## 3. Cahier des Charges Détaillé : Volume II — Guide de Développement ($\ge 40$ pages)

Le Guide de Développement s'adresse aux développeurs, chercheurs contributeurs et mainteneurs souhaitant comprendre l'architecture interne C++/Fortran, optimiser les performances ou implémenter de nouveaux modèles.

### Structure obligatoire et pagination cible (Total : 42 à 50 pages)

#### Chapitre 1 : Architecture Logicielle Hybride C++ / Fortran (5 à 6 pages)
* **1.1. Vision d'ensemble et paradigmes de conception :** séparation des responsabilités entre la couche haut niveau C++ (abstraction, polymorphisme, gestion du cycle de vie) et la couche bas niveau Fortran (noyaux numériques vectorisés).
* **1.2. Cartographie intégrale de l'arborescence :** description détaillée de chaque répertoire et de chaque fichier source du dépôt.
* **1.3. Graphe de dépendances et cycle d'exécution complet :** diagrammes de séquence détaillant l'initialisation, la boucle temporelle principale, un sous-pas de Runge-Kutta / une itération PIMPLE et la finalisation.

#### Chapitre 2 : Interopérabilité C++ $\leftrightarrow$ Fortran (6 à 8 pages)
* **2.1. Standards et conventions d'interface :** utilisation du module intrinsèque `iso_c_binding` en Fortran et des déclarations `extern "C"` en C++.
* **2.2. Table complète de correspondance des types et structures :** types scalaires, pointeurs (`type(c_ptr)` $\leftrightarrow$ `void*`), tableaux de taille connue et pointeurs de tableaux (`c_f_pointer`).
* **2.3. Gestion de la disposition mémoire (Row-Major vs Column-Major) :** règles d'indexation des champs 2D/3D (ex. : `field(i, j, k)` en Fortran correspondant à `field[k][j][i]` en C++), macros et accesseurs sécurisés.
* **2.4. Gestion des chaînes de caractères, des entrées/sorties et des erreurs :** propagation des codes d'erreur et exceptions entre Fortran et C++ sans fuite mémoire.
* **2.5. Catalogue des points d'entrée inter-langages :** tableau exhaustif de toutes les fonctions passerelles du projet avec leurs signatures C++ et Fortran exactes.

#### Chapitre 3 : Structures de Données, Mémoire et Parallélisme HPC (8 à 10 pages)
* **3.1. Représentation interne du maillage et des champs :** connectivité cellule-face-nœud, stockage AoS (*Array of Structures*) vs SoA (*Structure of Arrays*), alignement mémoire pour la vectorisation SIMD (AVX2/AVX-512).
* **3.2. Gestion des cellules fantômes (*Halo / Ghost Cells*) :** épaisseur des couches fantômes selon l'ordre du schéma spatial, ordre de mise à jour aux frontières physiques et aux interfaces entre sous-domaines.
* **3.3. Parallélisme distribué MPI :** stratégie de partitionnement de domaine (cartésien ou par graphe METIS/ParMETIS), communications point-à-point non bloquantes (`MPI_Isend` / `MPI_Irecv`), recouvrement calcul/communication, partage des communicateurs entre C++ et Fortran.
* **3.4. Parallélisme partagé et accélération GPU (OpenMP / CUDA / HIP) :** coloration de mailles ou opérations atomiques pour éviter les conditions de course (*race conditions*), gestion des transferts mémoire hôte-périphérique.

#### Chapitre 4 : Correspondance Théorie $\leftrightarrow$ Discrétisation $\leftrightarrow$ Code Source (10 à 12 pages)
* **Section pivot du guide :** pour chaque opérateur numérique et chaque modèle physique majeur du code, établir le lien explicite en trois volets :
  1. **Équation discrète exacte** (avec indices de mailles $i, j, k$ ou $P, N$, indices de faces $f$, et indices temporels $n, n+1$).
  2. **Localisation exacte dans le code** (`chemin/fichier.cpp` $\rightarrow$ `Classe::methode()`, puis appel vers `chemin/module.f90` $\rightarrow$ `subroutine nom_routine`).
  3. **Analyse ligne par ligne de l'implémentation** (boucles critiques, traitement des termes sources, régularisations numériques, astuces d'optimisation).

#### Chapitre 5 : Référence Exhaustive de l'API Interne (8 à 10 pages)
* Documentation technique complète générée et structurée de toutes les interfaces internes :
  * **Côté C++ :** hiérarchie des classes, membres de données, constructeurs/destructeurs, méthodes publiques et protégées, signatures complètes et invariants.
  * **Côté Fortran :** hiérarchie des modules (`use`), types dérivés (`type`), interfaces génériques, fonctions et sous-routines avec spécification systématique de `intent(in)`, `intent(out)`, `intent(inout)`.

#### Chapitre 6 : Stratégie de Vérification et Validation (V&V) et CI/CD (5 à 6 pages)
* **6.1. Architecture de la suite de tests :** tests unitaires C++ (GoogleTest / Catch2), tests unitaires Fortran (pFUnit / FUnit), tests d'intégration MPI.
* **6.2. Vérification formelle par la Méthode des Solutions Manufacturées (MMS) :** construction des termes sources analytiques, protocole de raffinement de maillage ($h, h/2, h/4, h/8$) et vérification automatique de l'ordre de convergence observé :
  $$p_{\text{obs}} = \frac{\ln\left(\|e_{2h}\|_{L_2} / \|e_h\|_{L_2}\right)}{\ln(2)}$$
* **6.3. Benchmarks de validation physique automatisés :** critères d'acceptation quantitatifs et tolérances numériques.
* **6.4. Pipeline d'Intégration Continue (CI/CD) :** étapes de compilation multi-compilateurs, analyse statique (`clang-tidy`, `-fcheck=all`), détection des fuites mémoire (Valgrind / AddressSanitizer) et tests de non-régression.

#### Chapitre 7 : Standards de Codage, Profilage et Guide d'Extension (4 à 5 pages)
* **7.1. Charte de codage C++ et Fortran :** règles de nommage, formatage (`clang-format`, `fprettify`), gestion RAII, interdiction des variables globales mutables (`SAVE` / `COMMON`).
* **7.2. Profilage et optimisation de performance :** utilisation de `perf`, Intel VTune, Score-P, et analyse des métriques roofline.
* **7.3. Tutoriel développeur pas à pas :** exemple complet d'ajout d'une nouvelle fonctionnalité de bout en bout (ex. : ajout d'un nouveau schéma convectif ou d'une équation de transport scalaire passif, de la lecture du paramètre en C++ jusqu'au noyau de flux en Fortran et son test MMS).

---

## 4. Cahier des Charges Détaillé : Volume III — Recueil de Tutoriels ($\ge 40$ pages)

Le Recueil de Tutoriels permet de parcourir **100 % des fonctionnalités du solveur** à travers une progression pédagogique allant des écoulements canoniques simples aux simulations multiphysiques et massivement parallèles.

### Gabarit obligatoire pour chaque tutoriel (4 à 5 pages par tutoriel)
Chaque tutoriel doit impérativement suivre cette structure en 7 points :
1. **Objectifs pédagogiques et fonctionnalités du code sollicitées** (liste des modules C++/Fortran activés).
2. **Description physique et formulation mathématique du cas** (géométrie, domaine $\Omega$, conditions aux limites exactes, nombres adimensionnels).
3. **Génération et analyse du maillage** (script de génération Gmsh ou paramètres du mailleur interne, nombre de cellules, raffinement près des parois).
4. **Fichier de configuration complet et commenté** (reproduction intégrale du fichier d'entrée avec explication de chaque choix numérique).
5. **Exécution de la simulation** (commande exacte séquentielle et MPI, temps de calcul indicatif, suivi des résidus).
6. **Analyse des résultats et post-traitement** (comparaison quantitative avec la solution analytique ou les données de la littérature, script Python complet de tracé et calcul d'erreurs).
7. **Exercices d'approfondissement** (étude de convergence en maillage, changement de schéma convectif ou variation du nombre de Reynolds/Rayleigh).

### Catalogue des 9 Tutoriels Couvrant l'Ensemble du Code (Total : 40 à 48 pages)

* **Tutoriel 1 : Prise en main et vérification analytique — Écoulement de Poiseuille et Couette plan**
  * *Fonctionnalités couvertes :* Compilation, structure d'un cas-test, maillage cartésien, conditions de Dirichlet/Neumann/périodicité, comparaison avec le profil parabolique exact $u(y)$.
* **Tutoriel 2 : Écoulement laminaire stationnaire — Cavité entraînée 2D et 3D (*Lid-Driven Cavity*)**
  * *Fonctionnalités couvertes :* Couplage vitesse-pression, schémas convectifs d'ordre 2, singularités de coin, validation quantitative à $Re = 100$, $Re = 1000$ et $Re = 3200$ par rapport aux profils de Ghia et al.
* **Tutoriel 3 : Écoulement instationnaire et détachement tourbillonnaire — Cylindre circulaire**
  * *Fonctionnalités couvertes :* Maillage non structuré / curviligne (Gmsh), pas de temps adaptatif ($CFL$), intégration des efforts surfaciques (coefficients de traînée $C_d$ et de portance $C_l$), analyse spectrale FFT et calcul du nombre de Strouhal $St$ à $Re = 100$.
* **Tutoriel 4 : Transfert thermique et convection naturelle — Cavité différentiellement chauffée**
  * *Fonctionnalités couvertes :* Couplage énergie / quantité de mouvement, terme source de flottabilité (Boussinesq), calcul du nombre de Nusselt local et moyen $Nu$ pour $Ra = 10^3$ à $Ra = 10^6$.
* **Tutoriel 5 : Modélisation de la turbulence — Écoulement sur une marche descendante (*Backward-Facing Step*)**
  * *Fonctionnalités couvertes :* Modèles RANS ($k$-$\omega$ SST), traitement de paroi ($y^+$), longueur de recollement et profils de contraintes de Reynolds.
* **Tutoriel 6 : Écoulements diphasiques à surface libre — Rupture de barrage (*Dam Break*) et ascension d'une bulle**
  * *Fonctionnalités couvertes :* Suivi d'interface (VOF / Level-Set), schémas compressifs d'interface, saut de densité/viscosité, tension superficielle (benchmark de Hysing).
* **Tutoriel 7 : Écoulements compressibles / Ondes de choc ou transport multi-espèces** *(adapter selon les modules exacts du dépôt)*
  * *Fonctionnalités couvertes :* Solveurs de Riemann, limiteurs TVD/WENO, lois d'état thermodynamiques, capture de discontinuités (tube à choc de Sod ou écoulement supersonique sur rampe).
* **Tutoriel 8 : Simulation Haute Performance (HPC) et passage à l'échelle (*Scaling*)**
  * *Fonctionnalités couvertes :* Décomposition de domaine MPI sur plusieurs dizaines/centaines de cœurs, écriture parallèle HDF5/XDMF, procédure d'arrêt et reprise (*Checkpoint/Restart*), mesure du *speedup* et de l'efficacité parallèle.
* **Tutoriel 9 : Couplage avancé et personnalisation utilisateur**
  * *Fonctionnalités couvertes :* Définition de conditions aux limites spatio-temporelles complexes, termes sources personnalisés, extraction de diagnostics avancés à la volée.

---

## 5. Architecture de Fichiers et Chaîne de Compilation `.md` $\rightarrow$ `.pdf`

### 5.1. Arborescence recommandée du répertoire `docs/`

```text
docs/
├── INSTRUCTIONS_DOCUMENTATION_CFD.md   # Le présent cahier des charges
├── Makefile                            # Automatisation de la génération MD + PDF
├── Doxyfile                            # Extraction API C++ et Fortran
├── templates/
│   └── eisvogel.latex                  # Gabarit LaTeX professionnel pour Pandoc
├── assets/
│   ├── figures/                        # Schémas vectoriels (SVG/PDF), courbes de validation
│   └── logos/                          # Logo du solveur
├── volume1_user_guide/
│   ├── 00_metadata.yaml                # Métadonnées Pandoc/LaTeX (Titre, Auteurs, Version)
│   ├── 01_introduction.md
│   ├── 02_installation_hpc.md
│   ├── 03_physics_models.md
│   ├── 04_numerical_methods.md
│   ├── 05_mesh_and_bc.md
│   ├── 06_input_dictionary.md
│   └── 07_execution_post_troubleshooting.md
├── volume2_dev_guide/
│   ├── 00_metadata.yaml
│   ├── 01_architecture.md
│   ├── 02_cpp_fortran_interop.md
│   ├── 03_data_structures_hpc.md
│   ├── 04_theory_to_code_map.md
│   ├── 05_api_reference.md
│   ├── 06_vandv_cicd.md
│   └── 07_coding_standards_extension.md
├── volume3_tutorials/
│   ├── 00_metadata.yaml
│   ├── 01_tuto_poiseuille_couette.md
│   ├── 02_tuto_lid_driven_cavity.md
│   ├── 03_tuto_cylinder_vortex_shedding.md
│   ├── 04_tuto_natural_convection.md
│   ├── 05_tuto_turbulent_step.md
│   ├── 06_tuto_multiphase_dambreak_bubble.md
│   ├── 07_tuto_compressible_or_species.md
│   ├── 08_tuto_hpc_scaling_restart.md
│   └── 09_tuto_custom_coupling.md
└── build/
    ├── Guide_Utilisateur_CFD.pdf       # Livrable PDF Volume I (>= 40 pages)
    ├── Guide_Developpement_CFD.pdf     # Livrable PDF Volume II (>= 40 pages)
    └── Recueil_Tutoriels_CFD.pdf       # Livrable PDF Volume III (>= 40 pages)
```

### 5.2. Configuration des métadonnées YAML pour Pandoc + XeLaTeX (`00_metadata.yaml`)

Placer en tête de chaque volume un fichier `00_metadata.yaml` garantissant une mise en page de livre technique de $\ge 40$ pages :

```yaml
---
title: "Solveur CFD — Guide Utilisateur Complet"
subtitle: "Formulation Physique, Méthodes Numériques et Référence des Paramètres"
author: "Équipe de Développement CFD"
date: "\\today"
lang: fr-FR
documentclass: report
classoption: [a4paper, 11pt, oneside]
geometry: "top=2.5cm, bottom=2.5cm, left=2.5cm, right=2.5cm"
toc: true
toc-depth: 3
lof: true
lot: true
numbersections: true
colorlinks: true
linkcolor: NavyBlue
urlcolor: RoyalBlue
toccolor: Black
highlight-style: tango
header-includes:
  - \usepackage{amsmath, amssymb, amsfonts, bm}
  - \usepackage{booktabs, longtable, array, multirow}
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhead[L]{\leftmark}
  - \fancyhead[R]{Documentation Officielle CFD}
  - \fancyfoot[C]{\thepage}
  - \usepackage{fvextra}
  - \DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,commandchars=\\\{\}}
---
```

### 5.3. Script de compilation et contrôle du volume de pages (`Makefile`)

```makefile
PANDOC = pandoc
PDF_ENGINE = xelatex
PANDOC_FLAGS = --from=markdown+tex_math_dollars+pipe_tables+grid_tables+fenced_code_attributes \
               --pdf-engine=$(PDF_ENGINE) \
               --listings \
               --number-sections \
               --toc

BUILD_DIR = build

all: user_guide dev_guide tutorials check_pages

user_guide:
	mkdir -p $(BUILD_DIR)
	$(PANDOC) volume1_user_guide/00_metadata.yaml volume1_user_guide/[0-9]*.md \
		$(PANDOC_FLAGS) -o $(BUILD_DIR)/Guide_Utilisateur_CFD.pdf

dev_guide:
	mkdir -p $(BUILD_DIR)
	$(PANDOC) volume2_dev_guide/00_metadata.yaml volume2_dev_guide/[0-9]*.md \
		$(PANDOC_FLAGS) -o $(BUILD_DIR)/Guide_Developpement_CFD.pdf

tutorials:
	mkdir -p $(BUILD_DIR)
	$(PANDOC) volume3_tutorials/00_metadata.yaml volume3_tutorials/[0-9]*.md \
		$(PANDOC_FLAGS) -o $(BUILD_DIR)/Recueil_Tutoriels_CFD.pdf

check_pages:
	@echo "=== Vérification du nombre de pages (Minimum requis : 40 pages/volume) ==="
	@for pdf in $(BUILD_DIR)/*.pdf; do \
		pages=$$(pdfinfo "$$pdf" | awk '/Pages:/ {print $$2}'); \
		echo "$$pdf : $$pages pages"; \
		if [ "$$pages" -lt 40 ]; then \
			echo "ALERTE : $$pdf contient moins de 40 pages !"; \
		fi; \
	done
```

---

## 6. Protocole d'Exécution Pas-à-Pas (Mode d'Emploi)

Pour produire l'ensemble des documents sans perte de qualité ni troncature, exécuter la génération selon le séquencement suivant :

1. **Phase 1 — Audit et cartographie du dépôt :**
   * Analyser `CMakeLists.txt`, l'arborescence `src/`, les parsers d'entrée et les interfaces `iso_c_binding` / `extern "C"`.
   * Produire une fiche d'inventaire exhaustive listant tous les fichiers C++, tous les modules Fortran, toutes les équations implémentées et tous les cas-tests existants.
2. **Phase 2 — Rédaction itérative du Volume I (Guide Utilisateur) :**
   * Rédiger chaque chapitre (chapitres 1 à 7) dans son fichier `.md` dédié en respectant strictement les quotas de pagination.
   * Vérifier que **100 % des paramètres du parser d'entrée** figurent dans le chapitre 6.
3. **Phase 3 — Rédaction itérative du Volume II (Guide de Développement) :**
   * Rédiger chaque chapitre (chapitres 1 à 7) dans son fichier `.md` dédié.
   * Vérifier que chaque schéma numérique et modèle physique du Volume I dispose de sa fiche de correspondance **Équation discrète $\leftrightarrow$ C++ $\leftrightarrow$ Fortran** dans le chapitre 4.
4. **Phase 4 — Rédaction itérative du Volume III (Recueil de Tutoriels) :**
   * Rédiger les 9 tutoriels complets avec les fichiers d'entrée réels issus du dépôt et les scripts Python de post-traitement fonctionnels.
5. **Phase 5 — Compilation PDF et contrôle qualité :**
   * Exécuter `make all` pour générer les trois PDF et vérifier via `pdfinfo` que chaque volume dépasse les **40 pages**, qu'aucune formule LaTeX n'est cassée et qu'aucun lien interne n'est rompu.