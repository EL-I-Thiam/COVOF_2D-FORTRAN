# Tutoriel 22 : Campagnes Hybrides FC-WY-CVI-Marker (`RUN_FC-WY-CVI-Marker`)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/RUN_FC-WY-CVI-Marker/` et `Guides_Tutos/Run_FC-WY-CVI-Marker/` *(8 études ciblées + `verify_cvi_markers.py`)*.

Ce tutoriel transverse présente l'utilisation des marqueurs lagrangiens non connectés couplés à l'interpolation conservative CVI (Variante B exacte) et au suivi VOF-PLIC de Weymouth-Yue (`use_markers: true`), ainsi que le protocole de non-régression automatisé `verify_cvi_markers.py`.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

- Activer et paramétrer le nuage de marqueurs lagrangiens (`vof.use_markers`, `markers_per_cell`, `marker_advection_scheme: "CVI_RK2"`, `marker_curvature`, `marker_normals_plic`, `marker_max_per_cell`, `write_markers_vtk`).
- Valider mathématiquement la divergence discrète nulle de CVI Variante B face aux schémas bilinéaires classiques via `verify_cvi_markers.py`.
- Comparer les simulations `VOF seul` (`use_markers: false`) aux simulations hybrides `FC-WY-CVI-Marker` (`use_markers: true`) sur l'ensemble des 8 études dédiées.
- Exploiter simultanément les champs VTK eulériens (`.pvtr`/`.vtr`) et les nuages lagrangiens (`_markers_<step:06d>.vtk`) sous ParaView.

---

## 2. Inventaire des Huit Campagnes CVI-Marqueurs

L'ensemble des études testant les marqueurs CVI est regroupé dans `VALIDATION_TOTALE/Run_CASES/RUN_FC-WY-CVI-Marker/` (avec miroir documentaire dans `Guides_Tutos/Run_FC-WY-CVI-Marker/`) :

| N° | Fichier d'Instructions | Phénomène Physique Étudié | Gain Ciblé par `FC-WY-CVI-Marker` |
| :-: | :--- | :--- | :--- |
| **01** | `01_cas_goutte_statique_young_laplace_cvi_markers.md` | Équilibre capillaire de Young-Laplace | Courbure parabolique $\kappa^{\text{marker}}$ et réduction des courants parasites |
| **02** | `02_cas_tourbillon_cisaille_rider_kothe_cvi_markers.md` | Tourbillon cisaillé réversible de Rider-Kothe | Advection CVI sans divergence et préservation de la queue spiralée |
| **03** | `03_cas_chute_cylindre_dense_cvi_markers.md` | Chute libre d'un cylindre dense ($\sigma=0$) | Préservation de la circularité sous-maille et conservation d'aire |
| **04** | `04_cas_oscillation_bulle_lamb_cvi_markers.md` | Oscillations capillaires amorties de Lamb | Fidélité de la fréquence propre $f_2$ et de l'amortissement visqueux |
| **05** | `05_cas_instabilite_rayleigh_taylor_cvi_markers.md` | Instabilité gravitaire de Rayleigh-Taylor | Résolution fine des volutes secondaires de Kelvin-Helmholtz |
| **06** | `06_cas_nappe_liquide_kelvin_helmholtz_cvi_markers.md` | Cisaillement d'une nappe liquide mince | Suivi continu des deux interfaces minces sans rupture artificielle |
| **07** | `07_cas_ascension_bulle_2d_cvi_markers.md` | Ascension de bulle 2D (benchmarks de Hysing) | Intégrité des jupes latérales à grand nombre d'Eötvös/Bond |
| **08** | `08_cas_flambement_jet_visqueux_cvi_markers.md` | Flambement gravitationnel d'un jet visqueux (`JET_INLET`) | Capture des lobes repliés et ré-ensemencement sous injection continue |

---

## 3. Protocole de Vérification Automatique (`verify_cvi_markers.py`)

Avant de lancer l'une des 8 études complètes, exécuter le banc de test à 3 paliers depuis la racine du projet :

```bash
python3 Guides_Tutos/Run_FC-WY-CVI-Marker/verify_cvi_markers.py
```

Ce script vérifie automatiquement :
1. **Palier 1 :** Annulation à la précision machine de la divergence analytique de l'interpolation CVI Variante B sur 5 000 blocs MAC aléatoires à divergence nulle.
2. **Palier 2 :** Réduction des courants parasites et précision de $\Delta P = \sigma / R_0$ sur la goutte statique de Young-Laplace ($32\times 32$ et $64\times 64$).
3. **Palier 3 :** Invariance stricte MPI `np=1` vs `mpirun -np 4` (`migrate_markers_mpi`, bande fantôme $2.5 h$, routage RK2 4 passes) et parité C++17 / Fortran 2018.

---

## 4. Paramétrage JSON Type du Bloc `vof` / `markers`

```json
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.5,
    "use_markers": true,
    "markers_per_cell": 6,
    "marker_advection_scheme": "CVI_RK2",
    "marker_curvature": true,
    "marker_normals_plic": true,
    "marker_max_per_cell": 12,
    "write_markers_vtk": true
  }
```

### Visualisation Duale Euler-Lagrange sous ParaView
1. Ouvrir simultanément la série eulérienne `<case_name>_..pvtr` (isocontour `PhaseIndicator = 0.5`) et la série lagrangienne `<case_name>_markers_..vtk`.
2. Représenter les marqueurs par des glyphes sphériques colorés par l'attribut `Curvature` pour inspecter la régularité géométrique le long de l'interface.
