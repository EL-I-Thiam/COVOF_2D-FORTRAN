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
