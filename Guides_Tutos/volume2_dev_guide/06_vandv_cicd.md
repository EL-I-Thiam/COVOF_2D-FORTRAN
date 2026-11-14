# Stratégie de Vérification et Validation (V&V) et CI/CD

Dans les codes de calcul scientifique de haute fidélité, la crédibilité des résultats repose sur une distinction méthodologique fondamentale établie par les normes ASME V&V 20 et AIAA :
- **La Vérification (*Verification*) :** Répond à la question *« Le code résout-il les équations correctement ? »* Elle quantifie les erreurs de programmation et vérifie que l'ordre de convergence formel théorique des discrétisations spatiales et temporelles est rigoureusement atteint.
- **La Validation (*Validation*) :** Répond à la question *« Le code résout-il les bonnes équations physiques ? »* Elle compare les prédictions numériques à des solutions analytiques exactes ou des mesures expérimentales de référence.

Ce chapitre décrit l'infrastructure complète de tests unitaires, la vérification formelle par la méthode des solutions manufacturées (MMS), la cartographie des 24 bancs de validation et le pipeline CI/CD automatisé de COVOF 2D.

---

## 6.1. Architecture de la Suite de Tests

Le harnais de test de COVOF 2D s'articule autour d'une hiérarchie à trois niveaux :

```
                  +-----------------------------------+
                  |     Tests Système & Validation    |  (24 cas complets,
                  |     (Python + MPI + PyVista)      |   comparaison analytique)
                  +-----------------+-----------------+
                                    |
                  +-----------------+-----------------+
                  |      Tests d'Intégration MPI      |  (Échanges halos,
                  |    (C++ ctest / Fortran ctest)    |   réductions cartésiennes)
                  +-----------------+-----------------+
                                    |
                  +-----------------+-----------------+
                  |    Tests Unitaires Bas Niveau     |  (Opérateurs MAC, CSR,
                  |   (Catch2 / GoogleTest / FUnit)   |   flux VOF, solveurs)
                  +-----------------------------------+
```

### 6.1.1. Tests Unitaires C++ (Catch2)

Les tests unitaires C++ valident indépendamment chaque composant isolé sans solliciter le sous-système MPI global :
- `test_mesh.cpp` : Vérification des métriques de pas d'espace $\Delta x, \Delta y$, des coordonnées de faces et de centres.
- `test_csrmatrix.cpp` : Test du produit matrice-vecteur creux SpMV contre une multiplication dense de référence.
- `test_linear_solvers.cpp` : Résolution d'un système linéaire test symétrique défini positif et asymétrique par CG, BiCGSTAB2, GMRES et MINRES avec validation du résidu $\|A \mathbf{x} - \mathbf{b}\|_2 / \|\mathbf{b}\|_2 < 10^{-12}$.
- `test_vof_geometry.cpp` : Calcul de la normale unitaire et constante de plan PLIC sur des fractions volumiques synthétiques connues avec conservation du volume à $10^{-15}$.

### 6.1.2. Tests Unitaires Fortran (pFUnit / Harness Natif)

Côté Fortran, chaque module est accompagné d'un programme de test dédié compilé avec l'exécutable `ctest` :
- `test_module_navier_stokes.f90` : Vérification de la conservation stricte de la divergence et de la symétrie de l'opérateur de diffusion visqueuse.
- `test_module_preconditioners.f90` : Évaluation du taux de réduction du conditionnement $\kappa(P^{-1} A)$ pour SSOR, ILU(0) et BT-SBP.

---

## 6.2. Vérification Formelle par la Méthode des Solutions Manufacturées (MMS)

La méthode des solutions manufacturées (MMS) constitue la référence absolue pour vérifier l'exactitude de l'implémentation des équations de Navier-Stokes.

### 6.2.1. Construction des Solutions Analytiques Manufacturées

On choisit délibérément un champ de vitesse et de pression régulier et infiniment différentiable $\mathcal{C}^\infty$ ne satisfaisant pas *a priori* les équations homogènes :
$$
\begin{aligned}
u_{\text{mms}}(x, y, t) &= U_0 \sin\left(\frac{2\pi x}{L}\right) \cos\left(\frac{2\pi y}{H}\right) e^{-\omega t} \\
v_{\text{mms}}(x, y, t) &= - U_0 \left(\frac{H}{L}\right) \cos\left(\frac{2\pi x}{L}\right) \sin\left(\frac{2\pi y}{H}\right) e^{-\omega t} \\
p_{\text{mms}}(x, y, t) &= P_0 \cos\left(\frac{2\pi x}{L}\right) \cos\left(\frac{2\pi y}{H}\right) e^{-\omega t}
\end{aligned}
$$
Par construction, ce champ satisfait rigoureusement la contrainte d'incompressibilité $\nabla \cdot \mathbf{u}_{\text{mms}} = 0$.

En injectant ces expressions dans l'équation de quantité de mouvement continue, on en déduit les termes sources analytiques exacts $\mathbf{S}_{\text{mms}} = [S_u, S_v]^T$ :
$$
\mathbf{S}_{\text{mms}} = \rho \left( \frac{\partial \mathbf{u}_{\text{mms}}}{\partial t} + \mathbf{u}_{\text{mms}} \cdot \nabla \mathbf{u}_{\text{mms}} \right) + \nabla p_{\text{mms}} - \mu \nabla^2 \mathbf{u}_{\text{mms}}
$$
Ces termes sources sont injectés comme forçage volumique dans le second membre du solveur discret.

### 6.2.2. Protocole de Raffinement et Calcul des Normes d'Erreur

La simulation est exécutée sur une séquence de quatre maillages cartésiens raffinés successivement par un facteur 2 :
$$h_1 = \frac{L}{32}, \quad h_2 = \frac{L}{64}, \quad h_3 = \frac{L}{128}, \quad h_4 = \frac{L}{256}$$

À l'instant final $t = T$, l'erreur discrète locale est calculée sur chaque maille :
$$e_{i,j} = \phi_{i,j}^{\text{num}} - \phi_{\text{mms}}(x_i, y_j, T)$$

On évalue les trois normes canoniques d'erreur :
$$
\begin{aligned}
\|e\|_{L_1} &= \frac{1}{N} \sum_{i,j} |e_{i,j}| \\
\|e\|_{L_2} &= \sqrt{\frac{1}{N} \sum_{i,j} |e_{i,j}|^2} \\
\|e\|_{L_\infty} &= \max_{i,j} |e_{i,j}|
\end{aligned}
$$

### 6.2.3. Formule du Taux de Convergence Asymptotique Observé

L'ordre de convergence effectif $p_{\text{obs}}$ entre deux niveaux de grille consécutifs de pas $2h$ et $h$ est calculé par :
$$
p_{\text{obs}} = \frac{\ln\left( \frac{\|e_{2h}\|}{\|e_h\|} \right)}{\ln(2)}
$$

Pour les schémas spatiaux d'ordre 2 (Centré, Upwind second ordre, BDF2 temporel), le critère d'acceptation de vérification impose :
$$1.95 \le p_{\text{obs}} \le 2.05$$
Toute valeur inférieure à 1.90 déclenche automatiquement un échec du test de vérification formelle.

---

## 6.3. Bancs d'Essai de Validation et Tolérances d'Acceptation

Le dossier `Validation_Totale` regroupe 24 cas tests couvrant 100 % des régimes d'écoulement du solveur. Chaque cas est muni de tolérances automatiques :

| N° Cas | Nom du Banc d'Essai | Régime / Modèles Activés | Grandeur de Validation | Tolérance Critique $\epsilon_{\text{tol}}$ |
| :---: | :--- | :--- | :--- | :---: |
| **01** | `test_taylor_green_vortex` | Incompressible, décroissance visqueuse | Énergie cinétique $E_k(t)$ | $\Delta E_k / E_k < 10^{-4}$ |
| **02** | `test_manufactured_solution` | Vérification MMS couplée | Ordre spatial $p_{\text{obs}}$ | $p_{\text{obs}} \ge 1.95$ |
| **03** | `test_poiseuille_channel` | Écoulement de canal laminaire | Profil parabolique $u(y)$ | $\|u - u_{\text{anal}}\|_{L_2} < 10^{-5}$ |
| **04** | `test_couette_flow` | Paroi supérieure mobile | Profil linéaire $u(y)$ | $\|u - u_{\text{anal}}\|_{L_\infty} < 10^{-6}$ |
| **05** | `test_poiseuille_periodic` | Périodicité $x$ avec gradient imposé | Débit massique $Q$ | $\Delta Q / Q < 10^{-5}$ |
| **06** | `test_laplace_static_bubble` | VOF-PLIC + Tension CSF | Saut de pression $\Delta p = \sigma / R$ | Erreur relative $< 1.5\%$ |
| **07** | `test_rider_kothe_vortex` | Advection VOF pure dans vortex | Conservation masse et forme | Perte volume $< 10^{-12}$ |
| **08** | `test_mpi_scaling_sync` | Parallélisme fort 1 à 64 cœurs | Débit & Vitesse | Écart séq/MPI $< 10^{-14}$ |
| **09** | `test_checkpoint_restart` | Sauvegarde binaire et reprise | Continuité solution à $t_{\text{restart}}$ | $\|u_{\text{repr}} - u_{\text{ref}}\|_\infty = 0$ |
| **10** | `test_multiphase_couette` | Deux fluides stratifiés $(\mu_1 \ne \mu_2)$ | Profil brisé analytique | $\|u - u_{\text{anal}}\|_{L_2} < 5 \times 10^{-4}$ |
| **11** | `test_compressible_piston` | Acoustique All-Mach (C-AAL) | Période oscillation piston | Écart théorie $< 0.8\%$ |
| **12** | `test_oscillating_bubble` | Oscillation diphasique de Lamb | Fréquence propre $\omega_2$ | Écart analytique $< 2.0\%$ |
| **13** | `test_bubble_compression` | Bulle compressible dans liquide | Rayon dynamique $R(t)$ | Écart Rayleigh-Plesset $< 2.5\%$ |
| **14** | `test_bicgstab2_krylov` | Algèbre linéaire de Krylov | Résidu linéaire $\|Ax - b\|$ | $\|r\|_2 / \|b\|_2 < 10^{-10}$ |
| **15** | `test_porous_drag_darcy` | Milieu poreux de Brinkman | Loi de Darcy $\Delta p = \frac{\mu L}{K} U$ | Écart théorique $< 0.5\%$ |
| **16** | `test_cylinder_flow_re100` | Écoulement autour d'un cylindre | Strouhal $St$ et traînée $C_d$ | $St = 0.165 \pm 0.005$ |
| **17** | `test_lid_driven_cavity` | Cavité entraînée $Re=1000$ | Profils $u(y)$ et $v(x)$ | Écart Ghia et al. $< 1.2\%$ |
| **18** | `test_natural_convection_ra1e5` | Convection de Rayleigh-Bénard | Nombre de Nusselt moyen $\overline{Nu}$ | Écart de Vahl Davis $< 1.0\%$ |
| **19** | `test_brinkman_cylinder_penalty` | Pénalisation solide immergé | Recirculation arrière $L_r$ | Écart corps maillé $< 2.0\%$ |
| **20** | `test_dam_break_water` | Rupture de barrage diphasique | Position front d'onde $X(t)$ | Écart expérimental Martin & Moyce $< 3\%$ |
| **21** | `test_rising_bubble_hysing` | Bulle ascendante (Benchmark Hysing) | Vitesse ascensionnelle terminale | Conforme benchmark à $1.5\%$ |
| **22** | `test_multiphase_natural_convec`| Convection thermique diphasique | Flux thermique interfacial | Bilan d'énergie stationnaire $< 10^{-4}$ |
| **23** | `test_krylov_solvers_preconds` | Comparaison croisée solveurs/préconds | Nombre d'itérations relatives | Convergence identique C++/Fortran |
| **24** | `test_viscous_jet_buckling` | Flambement d'un jet visqueux | Fréquence d'enroulement | Régime périodique stable |

---

## 6.4. Pipeline d'Intégration Continue (CI/CD)

Le pipeline d'intégration continue est orchestré par GitHub Actions (`.github/workflows/ci.yml`).

### 6.4.1. Matrice des Compilateurs et Environnements

À chaque *pull request* ou commit sur la branche principale, le code est compilé et testé sur une matrice croisée :
- **Compilateurs C++ :** GCC 11, GCC 12, GCC 13, Clang 14, Clang 16, Intel oneAPI ICX.
- **Compilateurs Fortran :** GNU Fortran (`gfortran` 11, 12, 13), Intel oneAPI IFX.
- **Bibliothèques MPI :** OpenMPI 4.1+, MPICH 4.0+.

### 6.4.2. Analyseurs Dynamiques de Mémoire (*Sanitizers* et Valgrind)

Pour interdire toute régression silencieuse :
1. **Compilation Instrumentée ASan / UBSan :**
   ```bash
   cmake -B build_san -DCMAKE_BUILD_TYPE=Debug \
         -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
         -DCMAKE_Fortran_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
   ```
   Valide l'absence totale de dépassements de tampons (*buffer overflows*) sur les mailles fantômes.
2. **Audit Valgrind Memcheck :**
   ```bash
   valgrind --leak-check=full --show-leak-kinds=all --error-exitcode=1 ./covof_tests
   ```
   Garantit une allocation mémoire nette de zéro octet lors de la destruction des objets (`Simulation::~Simulation()`).

### 6.4.3. Test de Parité Stricte C++ $\leftrightarrow$ Fortran

Le script d'intégration lance simultanément le solveur C++ et le solveur Fortran sur la même configuration JSON et vérifie que la différence maximale relative sur chaque champ scalaire reste inférieure à la tolérance binaire :
$$
\max_{i,j} \frac{|\phi_{\text{C++}}(i, j) - \phi_{\text{Fortran}}(i, j)|}{\max(|\phi_{\text{C++}}|, 10^{-12})} < 10^{-12}
$$
Ce test garantit que les deux codes demeurent des jumeaux numériques parfaits au cours de leur cycle de vie.
