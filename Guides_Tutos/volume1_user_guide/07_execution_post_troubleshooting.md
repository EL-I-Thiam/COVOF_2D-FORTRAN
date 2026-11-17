# Exécution, Post-Traitement et Dépannage

## Lancement d'une Simulation en Ligne de Commande

Les exécutables **C++17** (`covof2d_cpp` et son lien `covof2d`) et **Fortran 2018** (`covof2d_fortran`) acceptent exactement la même interface en ligne de commande et fonctionnent en **MPI pur** (aucun thread OpenMP).

### Syntaxe d'Exécution

```bash
# Exécution séquentielle (1 processus MPI) :
./COVOF_2D-C++/bin/covof2d_cpp chemin/vers/mon_cas.json
./COVOF_2D-FORTRAN/bin/covof2d_fortran chemin/vers/mon_cas.json

# Exécution parallèle distribuée MPI (exemple sur 4 rangs MPI) :
mpirun -np 4 ./COVOF_2D-C++/bin/covof2d_cpp chemin/vers/mon_cas.json
mpirun -np 4 ./COVOF_2D-FORTRAN/bin/covof2d_fortran chemin/vers/mon_cas.json

# Affichage de l'aide intégrée :
./COVOF_2D-C++/bin/covof2d_cpp --help
./COVOF_2D-FORTRAN/bin/covof2d_fortran --help
```

Si aucun argument n'est fourni sur la ligne de commande, les deux programmes tentent de charger le fichier par défaut `cases/cavity.json`. En cas d'échec d'ouverture du fichier JSON, le rang 0 affiche un message d'erreur et déclenche `MPI_Abort(MPI_COMM_WORLD, 1)`.

---

## Interprétation des Journaux Console (`main.cpp` / `Config::print()` / `Simulation::run()`)

Dès le lancement, seul le processus maître (`rank == 0`) affiche l'en-tête d'accueil, le résumé complet de la configuration (`Config::print()` en C++), puis l'en-tête de la boucle temporelle.

### Exemple Réel de Sortie Console (C++17)

```text
================================================================================
   COVOF 2D — Solveur CFD Navier-Stokes Incompressible Parallèle (C++17/MPI)   
================================================================================
  Nombre de processus MPI  : 4
  Fichier de configuration : rising_bubble.json
================================================================================
========================================================
       COVOF 2D PARALLELE - CONFIGURATION DU CAS         
========================================================
 Mode calcul : NAVIER_STOKES
 Maillage : 64 x 128 cellules (Lx = 1 m, Ly = 2 m)
 dx = 0.015625 m, dy = 0.015625 m
 Physique : Diphasique | Convection : Oui | Tension de surface : Oui
 Fluide 1 : rho=1000 kg/m3, mu=10 Pa.s
 Fluide 2 : rho=100 kg/m3, mu=1 Pa.s
 Tension superficielle sigma = 24.5 N/m
 VOF : Advection=WEYMOUTH_YUE | Reconstruction=ELVIRA | Marqueurs CVI=CVI_RK2 (N_init=4, Courbure=Oui, NormalesPLIC=Oui)
 Gravite : g = (0, -0.98) m/s^2
   Schema temporel : GEAR | Linearisation convective : ADAMS_BASHFORTH | Lissage SVOF : Oui (k=0.500000)
 Temps : dt = 0.001 s, t_final = 3 s
 Solveur : MONOLITHIC | BiCGSTAB(2), tol=1e-10, max_iter=1000, r_al=500, uzawa_iter=1
 Sorties : Dossier 'output/bulle_cvi', intervalle VTK=50
========================================================

>>> DEMARRAGE DU CALCUL COVOF 2D PARALLELE <<<
Pas     Temps [s]     dt [s]      Iter Krylov Vitesse Max     Div Max         
-------------------------------------------------------------------------------------
10      0.01000       1.000e-03   14          1.2410e-02      3.1842e-11      
20      0.02000       1.000e-03   15          2.4635e-02      2.9410e-11      
...
------------------------------------------------------------------------------
>>> SIMULATION COVOF 2D TERMINEE AVEC SUCCES A t = 3.00000 s <<<
Fichiers de visualisation VTK et CSV disponibles dans : output/bulle_cvi
```

*(Note : en Fortran, la bannière de démarrage de boucle affiche `>>> DEMARRAGE DU CALCUL COVOF 2D FORTRAN 2018 PARALLELE <<<` et inclut systématiquement la colonne `Nu Moyen`, qui vaut `0.0000` lorsque `enable_energy = false` ; en C++, la colonne `Nu Moyen` n'est affichée que si `enable_energy = true`.)*

### Signification des Colonnes de la Boucle Temporelle

| Colonne | Signification Physique et Numérique |
| :--- | :--- |
| `Pas` | Indice de l'itération temporelle courante (`step`). |
| `Temps [s]` | Temps physique cumulé $t = \sum \Delta t$ $[\text{s}]$. |
| `dt [s]` | Pas de temps effectif $\Delta t$ utilisé à l'itération courante $[\text{s}]$. |
| `Iter Krylov` | Nombre cumulé d'itérations du solveur linéaire (`total_iters`) sur toutes les itérations d'Uzawa du pas. |
| `Vitesse Max` | Norme maximale de vitesse observée sur la grille globale (`global_max_u`, $[\text{m}\cdot\text{s}^{-1}]$). |
| `Div Max` | Divergence discrète maximale $\max_{i,j} \left|\frac{u_{i+1,j}-u_{i,j}}{\Delta x} + \frac{v_{i,j+1}-v_{i,j}}{\Delta y}\right|$ $[\text{s}^{-1}]$. |
| `Nu Moyen` | Nombre de Nusselt moyen à la paroi chaude $x=0$ (si `enable_energy = true`, ou toujours affiché en Fortran). |

---

## Fichiers de Sortie, Visualisation ParaView et Post-Traitement

Tous les fichiers générés sont placés dans le répertoire `output.output_dir` en utilisant le préfixe `output.case_name` (noté `<case>` ci-dessous).

### 1. Historique Temporel CSV : `<case>_history.csv`

Écrit à chaque pas de temps lorsque `output.write_csv_history = true` (`VTKWriter.cpp` / `module_vtk_writer.f90`).

- **Format C++ (9 colonnes) :**
  `step,time,dt,center_of_mass_y,rise_velocity,bubble_volume,kinetic_energy,max_u,max_div`
- **Format Fortran (11 colonnes) :**
  `step,time,dt,center_of_mass_y,rise_velocity,bubble_volume,kinetic_energy,max_u,max_div,avg_rho_gas,avg_p_gas`

Définition exacte des métriques calculées sur l'indicateur de la phase 2 (`1.0 - C`) :
- `center_of_mass_y` : ordonnée du barycentre de la phase 2, $y_c = \frac{\int y (1-C)\,d\Omega}{\int (1-C)\,d\Omega}$.
- `rise_velocity` : vitesse verticale moyenne de la phase 2, $v_c = \frac{\int v (1-C)\,d\Omega}{\int (1-C)\,d\Omega}$.
- `bubble_volume` : aire totale (volume 2D) de la phase 2, $\int (1-C)\,d\Omega$.
- `kinetic_energy` : énergie cinétique discrète aux nœuds de pression, $\frac{1}{2}\sum_{i,j}\rho_{i,j}\left(\frac{u_{i,j}^2 + u_{i+1,j}^2}{2} + \frac{v_{i,j}^2 + v_{i,j+1}^2}{2}\right)\Delta x \Delta y$.
- `max_u`, `max_div` : vitesse maximale au centre des mailles et divergence discrète maximale.
- `avg_rho_gas`, `avg_p_gas` (Fortran) : masse volumique gaz et pression moyenne dans la phase $1-C$ (pour suivi compressible).

### 2. Champs Volumiques Parallèles XML VTK : `.vtr` et `.pvtr`

Écrits au pas `0`, tous les `output.vtk_interval` pas, et au pas final :

1. **Pièces locales par rang MPI :** `<case>_p<rank>_<step:06d>.vtr` (ex. `rising_bubble_p0_000050.vtr`). Chaque processus écrit son propre sous-domaine `RectilinearGrid` avec ses coordonnées `XCoordinates`, `YCoordinates`, `ZCoordinates`.
2. **Fichier maître parallèle :** `<case>_<step:06d>.pvtr` (ex. `rising_bubble_000050.pvtr`), écrit par le rang 0 et référençant l'ensemble des `<Piece Extent="..." Source="<case>_p<rank>_<step:06d>.vtr"/>`.

Les tableaux `<CellData>` exportés dans chaque `.vtr` et `.pvtr` sont exactement :
- `PhaseIndicator` : fraction volumique brute $C_{i,j}$
- `PhaseIndicator_Smoothed` : fraction volumique lissée SVOF $C_{s,i,j}$ (présente si `C_s` est alloué/non vide)
- `Pressure` : pression $p_{i,j}$
- `Density` : masse volumique $\rho_{i,j}$
- `Viscosity` : viscosité dynamique $\mu_{i,j}$
- `Divergence` : divergence discrète $(\nabla_h \cdot \mathbf{u})_{i,j}$
- `Temperature` : température $T_{i,j}$
- `Curvature` : courbure locale $\kappa_{i,j}$
- `Velocity` : vecteur vitesse $(u_c, v_c, 0.0)$ à 3 composantes interpolé aux centres de mailles

### 3. Nuage de Marqueurs Lagrangiens CVI : `<case>_markers_<step:06d>.vtk`

Lorsque `vof.use_markers = true` et `vof.write_markers_vtk = true`, les rangs MPI rassemblent par `MPI_Gatherv` vers le rang 0 l'ensemble des marqueurs actifs et écrivent un fichier VTK Legacy ASCII `DATASET POLYDATA` (`InterfaceMarkers::write_vtk` / `write_markers_vtk`) :
- `POINTS <N> double` : coordonnées $(x_p, y_p, 0.0)$
- `VERTICES <N> <2N>` : sommets ponctuels pour affichage direct dans ParaView
- `POINT_DATA <N>` contenant :
  - `VECTORS MarkerNormal double` : normale unitaire $(n_{x,p}, n_{y,p}, 0.0)$ orientée vers $C=1$
  - `SCALARS MarkerCurvature double 1` : courbure locale $\kappa_p$ portée par le marqueur

### 4. Sorties Analytiques Spécifiques aux Modes de Vérification

Selon la configuration (`mode` ou `test_case`), les routines dédiées de `Simulation.cpp` et `module_simulation.f90` exportent automatiquement des fichiers CSV de métrologie exacte dans `output_dir` :

| Cas / Mode | Fichiers Générés | Colonnes Principales |
| :--- | :--- | :--- |
| `mode = "LAPLACIAN"` (`CAS_3_1`, `CAS_3_2_NEUMANN`) | `laplacian_error.csv`, `solver_history.csv`, `solution_field_2d.csv`, `profile_x_mid.csv`, `profile_y_mid.csv`, `solution_field.vtk` | `nx,ny,dx,dy,iters,err_l1,err_l2,err_linf` et `iter,residual` |
| `test_case = "CAS_3_7"` / `"POISEUILLE"` | `poiseuille_error.csv`, `profile_u_y.csv`, `profile_p_x.csv`, `solution_field_2d.csv`, `solution_field.vtk` | `nx,ny,dx,dy,iters,err_l1,err_l2,err_linf,err_grad_p,cpu_time` |
| `test_case = "CAS_3_8"` / `"CAS_3_8_POISEUILLE_UNSTEADY"` | `poiseuille_unsteady_error.csv`, `profile_u_y.csv`, `solution_field_2d.csv`, `solution_field.vtk` | `nx,ny,dx,dy,total_step,final_time,err_l1,err_l2,err_linf,cpu_time` |
| `mode = "VORTEX_ADVECTION"` | `vortex_error.csv` (+ `.vtr`/`.pvtr` et `_markers_*.vtk`) | `nx,ny,dt,T,advection,reconstruction,limiter,err_l1,err_l2,err_linf,mass_init,mass_final,mass_conservation` |

### Exploitation sous ParaView

1. **Champs volumiques :** Dans ParaView, ouvrir le groupe temporel `<case>_*.pvtr` et cliquer sur **Apply**.
2. **Iso-contour d'interface VOF :** Appliquer le filtre **Cell Data to Point Data**, puis **Contour** sur `PhaseIndicator` avec la valeur `0.5`.
3. **Visualisation des marqueurs CVI :** Ouvrir simultanément le groupe `<case>_markers_*.vtk`, appliquer une taille de point visible (*Point Size = 4*) ou le filtre **Glyph** orienté selon `MarkerNormal` et coloré par `MarkerCurvature`.
4. **Vecteurs et lignes de courant :** Utiliser **Glyph** ou **Stream Tracer** sur le champ vectoriel `Velocity`.

### Scripts Python de Post-Traitement (`pyvista`, `numpy`, `matplotlib`)

#### Script 1 : Extraction de Profils et Superposition des Marqueurs depuis `.pvtr` et `_markers_*.vtk`

```python
#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

def plot_fields_and_markers(pvtr_file, markers_vtk=None):
    data = pv.get_reader(pvtr_file).read()
    bounds = data.bounds  # [xmin, xmax, ymin, ymax, zmin, zmax]
    x_mid = 0.5 * (bounds[0] + bounds[1])

    # 1. Extraction d'un profil vertical au centre X = Lx / 2
    line = data.sample_over_line([x_mid, bounds[2], 0.0],
                                 [x_mid, bounds[3], 0.0], resolution=400)
    y = line.points[:, 1]
    u = line["Velocity"][:, 0]
    p = line["Pressure"]

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), dpi=200)
    axes[0].plot(u, y, "b-", lw=2, label="u(x_mid, y)")
    axes[0].set_xlabel("Vitesse horizontale u [m/s]")
    axes[0].set_ylabel("Ordonnée y [m]")
    axes[0].grid(True, ls="--", alpha=0.6)
    axes[0].legend()

    axes[1].plot(p, y, "r-", lw=2, label="p(x_mid, y)")
    axes[1].set_xlabel("Pression p [Pa]")
    axes[1].set_ylabel("Ordonnée y [m]")
    axes[1].grid(True, ls="--", alpha=0.6)
    axes[1].legend()

    if markers_vtk is not None:
        mk = pv.read(markers_vtk)
        print(f"Nombre de marqueurs CVI actifs : {mk.n_points}")
        print(f"Courbure moyenne |kappa|       : {np.mean(np.abs(mk['MarkerCurvature'])):.4f} m^-1")

    plt.tight_layout()
    plt.savefig("covof_postproc_profiles.png")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 postproc.py <case_000050.pvtr> [case_markers_000050.vtk]")
        sys.exit(1)
    mk_arg = sys.argv[2] if len(sys.argv) > 2 else None
    plot_fields_and_markers(sys.argv[1], mk_arg)
```

#### Script 2 : Suivi de l'Historique CSV (`<case>_history.csv`)

```python
#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

csv_path = sys.argv[1]
data = np.genfromtxt(csv_path, delimiter=",", names=True)

fig, ax = plt.subplots(1, 2, figsize=(10, 4), dpi=200)
ax[0].plot(data["time"], data["rise_velocity"], "k-", lw=1.8)
ax[0].set_xlabel("Temps [s]")
ax[0].set_ylabel("Vitesse d'ascension [m/s]")
ax[0].grid(True, ls="--", alpha=0.6)

ax[1].semilogy(data["time"], np.maximum(data["max_div"], 1e-16), "m-", lw=1.8)
ax[1].set_xlabel("Temps [s]")
ax[1].set_ylabel("Divergence max [s^-1]")
ax[1].grid(True, ls="--", alpha=0.6)

plt.tight_layout()
plt.savefig("history_diagnostics.png")
```

---

## Guide Complet de Dépannage (*Troubleshooting*)

| Symptôme Observé | Cause Physique / Numérique Probable | Action Corrective Recommandée |
| :--- | :--- | :--- |
| **Divergence de Krylov (`Iter Krylov` = `max_iter`, `NaN` dans les logs)** | 1. Pas de temps $\Delta t$ violant la CFL convective ou capillaire ($\Delta t_\sigma$).<br>2. Saut $\rho_1/\rho_2$ ou $\mu_1/\mu_2$ brutal sans lissage SVOF.<br>3. Préconditionneur inadapté au mode `MONOLITHIC`. | 1. Activer `time.adaptive_dt = true` avec `time.cfl_max = 0.35`.<br>2. Activer `physics.enable_svof_smoothing = true` (`svof_filter_width = 0.5`).<br>3. Avec `solver.type = "MONOLITHIC"`, utiliser impérativement `solver.precond_type = "BLOCK_TRIANGULAR"`. |
| **Divergence maximale `Div Max` élevée ($> 10^{-4}$)** | Trop peu d'itérations d'Uzawa en mode AL/AAL ou tolérance linéaire trop relâchée. | Passer `solver.type = "MONOLITHIC"` (ou `"AAL"` avec `uzawa_iter = 5`, `uzawa_tol = 1.0e-8`) et abaisser `solver.tol` à `1.0e-9` ou `1.0e-10`. |
| **Courants parasites intenses autour d'une goutte statique** | Courbure CSF bruitée par différentiation finie de $C$. | Activer la méthode hybride CVI-Marqueurs (`vof.use_markers = true`, `vof.marker_curvature = true`, `vof.marker_advection_scheme = "CVI_RK2"`) et la reconstruction `"ELVIRA"`. |
| **Fragmentation d'un filament mince ou d'une jupe de bulle** | Sous-résolution sous-maille ($\delta < \Delta x$) ou perte de cohérence de normale. | Activer `vof.use_markers = true`, `vof.marker_normals_plic = true` et `vof.advection_method = "WEYMOUTH_YUE"` avec `cfl_vof <= 0.35`. |
| **Oscillations spatiales en damier dans les zones à fort gradient** | Schéma convectif `CENTERED` utilisé alors que $Re_\Delta > 2$. | Basculer `physics.convective_scheme = "UPWIND"` ou `"TVD_SUPERBEE"`, ou raffiner `(mesh.nx, mesh.ny)`. |
| **Dépassement ou blocage MPI lors du lancement** | Nombre de processus MPI incompatible avec `nx`, `ny`, ou rangs bloqués sur barrière suite à une erreur JSON. | Vérifier que $N_x / P_x \ge 2$ et $N_y / P_y \ge 2$, et contrôler le JSON (clés exactes, pas de virgule traînante en fin de bloc). |
| **Clé JSON ignorée silencieusement** | Utilisation d'un ancien nom de clé (ex. `physics.rho1`, `time.t_end`, `thermal.*`, `convection.*`, `NO_SLIP`). | Consulter le Chapitre 6 : utiliser `physics.fluid1.rho`, `time.t_final`, `physics.enable_energy`, `physics.convective_*`, `DIRICHLET`, `FREE_SLIP`, `BOX`. |
