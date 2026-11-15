
---

# 1. Principes Architecturaux et Garantie de Non-Régression

L'intégration de la méthode hybride **VOF-PLIC (Weymouth & Yue) + Interpolation CVI + Nuage de Marqueurs Lagrangiens Non Connectés** dans les deux moteurs jumeaux **`COVOF_2D-FORTRAN`** et **`COVOF_2D-C++`** obéit à quatre règles d'ingénierie strictes :

1. **Préservation à 100 % des méthodes existantes (Rétrocompatibilité bit-à-bit)** :
   - Tous les paramètres du nouveau module sont désactivés par défaut (`use_markers = .false.` / `false`, `marker_curvature = .false.` / `false`, `marker_normals_plic = .false.` / `false`).
   - Lorsque `use_markers = false`, aucun tableau de marqueurs n'est parcouru dans la boucle temporelle, et les sorties numériques (`.vtk` et `.csv`) des 14 cas de validation analytiques (`Run_Valide-CASES`) et des cas multiphysiques (`Run_CASES`) demeurent **strictement identiques bit-à-bit** à l'existant.
   - Toutes les méthodes déjà codées restent intactes et combinables :
     - **Solveurs vitesse-pression** : `MONOLITHIC` (FC), `AUGMENTED_LAGRANGIAN` (AL), `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (AAL), `CAAL`, `STANDARD_PROJECTION`.
     - **Solveurs linéaires et préconditionneurs** : `DIRECT_LU`, `BICGSTAB`, `BICGSTAB2`, `GMRES`, `MINRES`, `CG` avec `NONE`, `JACOBI`, `SSOR`, `ILU0`, `BLOCK_TRIANGULAR`, `BT_SBP`, `CUSTOM`.
     - **Schémas d'advection et de reconstruction VOF** : `WEYMOUTH_YUE`, `CIAM`, `VOF_TVD` ; `PLIC`, `ELVIRA`, `LVIRA`, `SLIC` ; limiteurs `SUPERBEE`, `VAN_LEER`, `MUSCL`, `HRIC`.
2. **Isolation modulaire** :
   - Toute la logique lagrangienne et géométrique des marqueurs est encapsulée dans un module dédié (`src/module_markers.f90` en Fortran 2018 ; `include/covof/InterfaceMarkers.hpp` et `src/InterfaceMarkers.cpp` en C++17).
3. **Conformité au halo de largeur 1 `(0:nx+1, 0:ny+1)`** :
   - L'interpolation CVI utilise exclusivement les indices `0` à `nx+1` (en $x$) et `0` à `ny+1` (en $y$) de la grille MAC existante, sans jamais élargir la couche de cellules fantômes MPI.
4. **Migration MPI en 2 passes dimensionnelles (4 voisins cardinaux)** :
   - La migration MPI des marqueurs réutilise les 4 rangs voisins cardinaux (`neighbor_left`, `neighbor_right`, `neighbor_bottom`, `neighbor_top`) de `mesh_2d_t` en deux passes successives (Passe X puis Passe Y), gérant automatiquement les traversées diagonales de coins.

---

# 2. Extension du Dictionnaire d'Entrée JSON (`vof` / `markers`)

Dans les fichiers de configuration `.json`, une sous-section optionnelle `"markers"` (ou des clés additionnelles dans `"vof"`) est ajoutée avec les valeurs par défaut suivantes :

```json
"vof": {
  "advection_method": "WEYMOUTH_YUE",
  "reconstruction_method": "PLIC",
  "tvd_limiter": "NONE",
  "cfl_vof": 0.3,
  "use_markers": true,
  "markers_per_cell": 4,
  "marker_advection_scheme": "CVI_RK2",
  "marker_curvature": true,
  "marker_normals_plic": false,
  "marker_max_per_cell": 8,
  "write_markers_vtk": true
}
```

| Clé JSON | Type | Défaut | Rôle physique et numérique |
| :--- | :---: | :---: | :--- |
| `use_markers` | `bool` | `false` | Active l'initialisation, l'advection CVI, la migration MPI et la maintenance du nuage de marqueurs. |
| `markers_per_cell` | `int` | `4` | Nombre cible $N_{\text{target}}$ de marqueurs semés sur le segment PLIC d'une cellule interfaciale nouvellement peuplée ($N_p = 0$). |
| `marker_advection_scheme` | `string` | `"CVI_RK2"` | Schéma d'advection des marqueurs : `"CVI_RK2"` (Heun d'ordre 2 avec $\mathbf{u}^n, \mathbf{u}^{n+1}$ et interpolation CVI à divergence nulle) ou `"BILINEAR_RK2"` (pour comparaison académique). |
| `marker_curvature` | `bool` | `false` | Active le calcul de la courbure capillaire $\kappa_{i,j}$ par ajustement parabolique local pondéré dans le repère $(\tau, \eta)$ sur les marqueurs du voisinage $3\times 3$ (avec repli automatique sur CSF si $M < 4$). |
| `marker_normals_plic` | `bool` | `false` | Permet d'orienter la normale de reconstruction PLIC $(n_x, n_y)$ par la normale moyenne des marqueurs présents dans la cellule. |
| `marker_max_per_cell` | `int` | `8` | Seuil maximal $N_{\max}$ au-delà duquel les marqueurs surnuméraires les plus proches dans une même cellule sont fusionnés/élagués. |
| `write_markers_vtk` | `bool` | `false` | Exporte le nuage de marqueurs au format VTK `POLYDATA` (`*_markers_XXXXXX.vtk`) à chaque pas d'écriture VTK. |

---

# 3. Feuille de Route Détaillée — Moteur Fortran 2018 (`COVOF_2D-FORTRAN`)

## Vue d'ensemble des fichiers impactés

| Étape | Fichier | Action | Impact sur les méthodes existantes |
| :---: | :--- | :--- | :--- |
| **F1** | `src/module_config.f90` | Ajout des 7 attributs dans `config_t` et lecture JSON dans `load_config` | Aucun (valeurs par défaut `.false.`) |
| **F2** | `src/module_fields.f90` | Ajout des tableaux optionnels `marker_nx(:, :)`, `marker_ny(:, :)`, `marker_kappa(:, :)`, `has_marker_kappa(:, :)` dans `flow_fields_t` | Aucun (tableaux passifs si `.not. use_markers`) |
| **F3** | `src/module_interface_reconstruction.f90` | Ajout de la routine utilitaire `clip_segment_in_cell(xmin, xmax, ymin, ymax, nx, ny, alpha, x1, y1, x2, y2, valid)` | Aucun (`cut_box_area`, `compute_alpha`, `reconstruct_plic` inchangés) |
| **F4** | `src/module_markers.f90` *(Nouveau)* | Création du module complet : `marker_t`, `marker_cloud_t`, `interpolate_velocity_cvi`, `advect_markers_rk2`, `migrate_markers_mpi`, `prune_and_reseed_markers`, `compute_curvature_from_markers` | Module indépendant |
| **F5** | `src/module_interface_vof.f90` | Ajout de l'instance `type(marker_cloud_t) :: cloud` dans `interface_vof_t` ; appel dans `advect_vof` avant `compute_curvature_and_csf` ; branchement dans `compute_curvature_and_csf` | Exécuté uniquement si `cfg%use_markers` |
| **F6** | `src/module_vtk_writer.f90` | Ajout de `write_markers_vtk(this, step, time, cloud)` (format ASCII VTK `POLYDATA`) | Exécuté uniquement si `cfg%write_markers_vtk` |
| **F7** | `src/module_simulation.f90` & `CMakeLists.txt` | Appel d'initialisation des marqueurs dans `sim_init`, écriture VTK dans `sim_run`, ajout de `src/module_markers.f90` dans `CMakeLists.txt` | Aucun |

---

## Spécification Technique Étape par Étape (Fortran 2018)

### Étape F1 : `src/module_config.f90`
Dans `type :: config_t`, ajouter :
```fortran
  ! --- Paramètres du module Marqueurs Lagrangiens & CVI ---
  logical            :: use_markers             = .false.
  integer            :: markers_per_cell        = 4
  character(len=32)  :: marker_advection_scheme = "CVI_RK2"
  logical            :: marker_curvature        = .false.
  logical            :: marker_normals_plic     = .false.
  integer            :: marker_max_per_cell     = 8
  logical            :: write_markers_vtk       = .false.
```
Dans `load_config`, lire ces clés depuis la section `"vof"` (et accepter également une section racine `"markers"`).

### Étape F3 : `src/module_interface_reconstruction.f90`
Ajouter la sous-routine publique `clip_segment_in_cell` qui calcule les deux extrémités exactes $(x_1, y_1)$ et $(x_2, y_2)$ de l'intersection de la droite $n_x x + n_y y + \alpha = 0$ avec le rectangle $[x_{\min}, x_{\max}] \times [y_{\min}, y_{\max}]$ :
```fortran
  pure subroutine clip_segment_in_cell(xmin, xmax, ymin, ymax, nx, ny, alpha, x1, y1, x2, y2, valid)
    real(dp), intent(in)  :: xmin, xmax, ymin, ymax, nx, ny, alpha
    real(dp), intent(out) :: x1, y1, x2, y2
    logical,  intent(out) :: valid
    ! Teste les 4 arêtes du rectangle, collecte les 2 points d'intersection où d(A)*d(B) <= 0
    ! et retourne valid = (n_pts == 2) avec la longueur du segment > 1.0e-14_dp.
  end subroutine clip_segment_in_cell
```

### Étape F4 : Nouveau fichier `src/module_markers.f90`

1. **Structures de données (`marker_t` et `marker_cloud_t`)** :
```fortran
  type :: marker_t
     real(dp) :: x      = 0.0_dp  ! Abscisse globale [m]
     real(dp) :: y      = 0.0_dp  ! Ordonnée globale [m]
     real(dp) :: nx     = 0.0_dp  ! Normale unitaire orientée vers C=1 (composante X)
     real(dp) :: ny     = 0.0_dp  ! Normale unitaire orientée vers C=1 (composante Y)
     real(dp) :: kappa  = 0.0_dp  ! Courbure locale estimée [1/m]
     integer  :: cell_i = 0       ! Indice local i in [1, nx_local]
     integer  :: cell_j = 0       ! Indice local j in [1, ny_local]
     logical  :: active = .false. ! Indicateur de validité
  end type marker_t

  type :: marker_cloud_t
     integer :: count    = 0
     integer :: capacity = 0
     type(marker_t), allocatable :: p(:)
  end type marker_cloud_t
```

2. **Interpolation Conservative de Vitesse (`interpolate_velocity_cvi`)** :
   Pour une position $(x, y)$ dans le sous-domaine local, identifier la cellule de pression $(i, j)$ :
   $$i = \max\left(1, \min\left(N_x, \left\lfloor \frac{x - x_{\min,\text{loc}}}{\Delta x} \right\rfloor + 1\right)\right), \quad j = \max\left(1, \min\left(N_y, \left\lfloor \frac{y - y_{\min,\text{loc}}}{\Delta y} \right\rfloor + 1\right)\right),$$
   et les coordonnées locales réduites $\xi = \frac{x - x_{\text{face}}(i)}{\Delta x} \in [0, 1]$, $\eta = \frac{y - y_{\text{face}}(j)}{\Delta y} \in [0, 1]$.
   Pour interpoler $u(x,y)$ et $v(x,y)$ avec continuité aux faces et divergence analytiquement nulle dans la cellule $(i,j)$ :
   - Sur la grille décalée, la divergence discrète dans $(i,j)$ est $D_{i,j} = \frac{u_{i+1,j}-u_{i,j}}{\Delta x} + \frac{v_{i,j+1}-v_{i,j}}{\Delta y} = 0$.
   - En utilisant les vitesses aux 4 coins de la cellule $(i,j)$ (Variante A, Section 2.5.C) ou les boîtes décalées (Variante B de Pusok et al., 2017) :
     ```fortran
     u_sw = 0.5_dp * (u_field(i, j-1)   + u_field(i, j))
     u_se = 0.5_dp * (u_field(i+1, j-1) + u_field(i+1, j))
     u_nw = 0.5_dp * (u_field(i, j)     + u_field(i, j+1))
     u_ne = 0.5_dp * (u_field(i+1, j)   + u_field(i+1, j+1))

     v_sw = 0.5_dp * (v_field(i-1, j)   + v_field(i, j))
     v_se = 0.5_dp * (v_field(i, j)     + v_field(i+1, j))
     v_nw = 0.5_dp * (v_field(i-1, j+1) + v_field(i, j+1))
     v_ne = 0.5_dp * (v_field(i, j+1)   + v_field(i+1, j+1))

     du_s = u_se - u_sw
     du_n = u_ne - u_nw
     dv_w = v_nw - v_sw
     dv_e = v_ne - v_se

     u_bilin = (1.0_dp - eta) * ((1.0_dp - xi) * u_sw + xi * u_se) + &
               eta            * ((1.0_dp - xi) * u_nw + xi * u_ne)
     v_bilin = (1.0_dp - xi)  * ((1.0_dp - eta) * v_sw + eta * v_nw) + &
               xi             * ((1.0_dp - eta) * v_se + eta * v_ne)

     if (use_cvi) then
        u_p = u_bilin + 0.5_dp * (dx / dy) * xi  * (1.0_dp - xi)  * (dv_e - dv_w)
        v_p = v_bilin + 0.5_dp * (dy / dx) * eta * (1.0_dp - eta) * (du_n - du_s)
     else
        u_p = u_bilin
        v_p = v_bilin
     end if
     ```
   *Vérification immédiate :* $\frac{\partial u_p}{\partial x} + \frac{\partial v_p}{\partial y} = \frac{\Delta u_S + \Delta u_N}{2\Delta x} + \frac{\Delta v_W + \Delta v_E}{2\Delta y} = \frac{1}{4}(D_{i,j-1} + 2D_{i,j} + D_{i,j+1}) + \dots = 0$, et tous les accès $(i\pm 1, j\pm 1)$ pour $i \in [1, N_x], j \in [1, N_y]$ restent dans `(0:nx+1, 0:ny+1)` !

3. **Advection RK2 (Heun) temporellement d'ordre 2 (`advect_markers_rk2`)** :
   Utilise `fields%u_old, fields%v_old` pour la vitesse $\mathbf{u}^n$ au point initial $\mathbf{x}_p^n$, et `fields%u, fields%v` pour la vitesse $\mathbf{u}^{n+1}$ au point prédicteur $\mathbf{x}_p^*$.

4. **Migration MPI en 2 passes (`migrate_markers_mpi`)** :
   Chaque marqueur migré est empaqueté en `5` réels double précision `[x, y, nx, ny, kappa]`.
   - **Passe X** : échange du nombre d'éléments `n_send_left, n_send_right` avec `mesh%neighbor_left` et `mesh%neighbor_right` via `MPI_Sendrecv`, puis transfert des tableaux de réels `buf_send_left(5*n_send_left)` et `buf_send_right(5*n_send_right)`.
   - **Passe Y** : idem avec `mesh%neighbor_bottom` et `mesh%neighbor_top` sur la liste incluant les marqueurs reçus à la Passe X.

5. **Calcul de courbure locale dans le repère tourné $(\tau, \eta)$ (`compute_curvature_from_markers`)** :
   Implémente l'algorithme exact de la Section 3.4 de la Synthèse (collecte sur voisinage $3\times 3$, filtrage angulaire $\mathbf{n}_p \cdot \mathbf{n}_{i,j}^{\text{ref}} \ge 0.5$, changement de repère $(\tau_p, \eta_p)$, résolution $3\times 3$ pondérée $\eta = a\tau^2 + b\tau + c$, extraction $\kappa_{i,j} = \frac{2a}{(1+b^2)^{3/2}}$ bornée par $\frac{2}{\min(\Delta x, \Delta y)}$).

### Étape F5 : Branchement dans `src/module_interface_vof.f90`
Dans `compute_curvature_and_csf` (ligne 620 de `module_interface_vof.f90`), après le calcul CSF standard de `kappa(i,j)` :
```fortran
    ! Si l'estimation de courbure par marqueurs est activée, substituer kappa(i,j)
    ! dans les cellules interfaciales disposant d'un ajustement parabolique valide
    if (cfg%use_markers .and. cfg%marker_curvature .and. allocated(fields%marker_kappa)) then
       do j = 1, ny
          do i = 1, nx
             if (fields%has_marker_kappa(i, j)) then
                kappa(i, j) = fields%marker_kappa(i, j)
             end if
          end do
       end do
       call exchange_halos_scalar(kappa, mesh)
    end if
```
Ainsi, le calcul CSF standard sert de **filet de sécurité automatique** (*fallback*) partout où l'interface n'a pas assez de marqueurs cohérents, et les formules *Balanced-Force* aux faces (lignes 643–661) restent **strictement inchangées** !

---

# 4. Feuille de Route Détaillée — Moteur C++17 (`COVOF_2D-C++`)

La parité architecturale entre `COVOF_2D-FORTRAN` et `COVOF_2D-C++` est maintenue classe par classe :

| Étape | Fichiers C++17 (`include/covof/` & `src/`) | Spécifications d'implémentation |
| :---: | :--- | :--- |
| **C1** | `Config.hpp`, `Config.cpp` | Ajout dans `struct VofConfig` (ou `Config`) de : `bool use_markers = false; int markers_per_cell = 4; std::string marker_advection_scheme = "CVI_RK2"; bool marker_curvature = false; bool marker_normals_plic = false; int marker_max_per_cell = 8; bool write_markers_vtk = false;` et parsing JSON correspondant. |
| **C2** | `InterfaceReconstruction.hpp`, `InterfaceReconstruction.cpp` | Ajout de la méthode statique `clipSegmentInCell(xmin, xmax, ymin, ymax, seg, x1, y1, x2, y2)`. |
| **C3** | `InterfaceMarkers.hpp`, `InterfaceMarkers.cpp` *(Nouveaux)* | Définition de `struct Marker { double x, y, nx, ny, kappa; int cell_i, cell_j; bool active; };` et de la classe `InterfaceMarkers` contenant `std::vector<Marker> markers_`, `interpolateVelocityCVI(...)`, `advectRK2(...)`, `migrateMPI(...)`, `pruneAndReseed(...)`, `computeCurvatureField(...)`. |
| **C4** | `InterfaceTracker.hpp`, `InterfaceTracker.cpp` | Ajout d'un membre `InterfaceMarkers marker_cloud_;` dans `InterfaceTracker`. Dans `InterfaceTracker::advect(...)`, après `applyVofBC` et `smoothDiffusion`, appeler `marker_cloud_.advectRK2(...)`, `marker_cloud_.migrateMPI(...)`, `marker_cloud_.pruneAndReseed(...)`, et `marker_cloud_.computeCurvatureField(...)` avant `computeCurvatureAndCSF(...)`. |
| **C5** | `VTKWriter.hpp`, `VTKWriter.cpp` | Ajout de `VTKWriter::writeMarkers(int step, double time, const InterfaceMarkers& cloud)` écrivant un fichier ASCII Legacy VTK `DATASET POLYDATA` avec `POINTS`, `VERTICES` et les attributs `Normals` (vecteur) et `Curvature` (scalaire). |
| **C6** | `CMakeLists.txt` | Ajout de `src/InterfaceMarkers.cpp` à la liste des sources de la cible `covof2d_cpp`. |

---

# 5. Protocole de Vérification et de Validation de l'Implémentation

Pour garantir une intégration sans défaut avant le lancement des grandes campagnes multiphysiques, l'implémentation devra franchir **4 paliers de vérification** :

1. **Palier 0 — Non-régression stricte (`use_markers: false`)** :
   Exécuter un cas analytique rapide (ex. Cas 3 Poiseuille et Cas 6 Goutte de Laplace avec `use_markers: false`) sur `covof2d_fortran` et `covof2d_cpp` et vérifier que les fichiers `*_history.csv` sont identiques à $10^{-16}$ près aux résultats antérieurs.

2. **Palier 1 — Vérification mathématique de l'incompressibilité CVI ($\nabla \cdot \mathbf{u}_{\text{CVI}} = 0$)** :
   Sur le champ tourbillonnaire de Rider-Kothe (Cas 27) et sur une solution Navier-Stokes monolithique FC, évaluer par différences finies à $\epsilon = 10^{-7}\,\text{m}$ la divergence numérique de l'interpolateur :
   $$\mathcal{D}_{\text{num}}(x,y) = \frac{u_{\text{CVI}}(x+\epsilon, y) - u_{\text{CVI}}(x-\epsilon, y)}{2\epsilon} + \frac{v_{\text{CVI}}(x, y+\epsilon) - v_{\text{CVI}}(x, y-\epsilon)}{2\epsilon}.$$
   Vérifier que $\max_{(x,y)} |\mathcal{D}_{\text{num}}(x,y)| \le 10^{-11}\,\text{s}^{-1}$ (contre $\sim 10^{-1}\,\text{s}^{-1}$ pour l'interpolation bilinéaire standard).

3. **Palier 2 — Précision géométrique de la courbure sur cercle exact** :
   Pour une goutte circulaire de rayon $R_0 = 0.25\,\text{m}$ ($\kappa_{\text{ex}} = 4.0\,\text{m}^{-1}$), vérifier que l'erreur relative $\|\kappa^{\text{marker}} - \kappa_{\text{ex}}\|_\infty / \kappa_{\text{ex}}$ décroît à l'ordre 2 sur la séquence de maillages $32^2, 64^2, 128^2$ ($< 0.5\,\%$ dès $64^2$).

4. **Palier 3 — Invariance MPI ($NP = 1$ vs $NP = 4$)** :
   Exécuter 100 pas de temps d'advection d'un cylindre traversant le point quadruple central $(x = L_x/2, y = L_y/2)$ sur 1 processus et sur $2\times 2 = 4$ processus MPI, et vérifier que le nombre total de marqueurs et la courbure moyenne sont invariants.
