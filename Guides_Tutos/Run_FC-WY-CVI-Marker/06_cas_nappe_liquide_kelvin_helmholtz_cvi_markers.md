# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 6 (FC-WY-CVI-MARKER) : CISAILLEMENT PÉRIODIQUE DE KELVIN-HELMHOLTZ D'UNE NAPPE LIQUIDE (CAS 3.2 *FLUIDS* 2021) AVEC MARQUEURS CVI (`06_cas_nappe_liquide_kelvin_helmholtz_cvi_markers`)

## 1. Mission et Rôle de l'Agent

Vous êtes un agent autonome chargé d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 6 de la méthode hybride FC-WY-CVI-Marker (`06_cas_nappe_liquide_kelvin_helmholtz_cvi_markers`)** sur les moteurs **C++17** (`covof2d_cpp`) et **Fortran 2018** (`covof2d_fortran`) de la plateforme **COVOF 2D**.

Ce cas correspond au **Cas 3.2 de [El Ouafa, Vincent & Le Chenadec (*Fluids* 2021)](https://doi.org/10.3390/fluids6010023)** couplé à l'architecture **Solveur Monolithique (FC) + VOF-PLIC (Weymouth & Yue) + Interpolation CVI + Marqueurs Non Connectés**. Il valide simultanément :
1. La conservation de l'énergie cinétique $E_k(t)/E_k(0)$ en écoulement d'Euler diphasique à grand rapport de densité ($r_\rho = \rho_2/\rho_1 \in \{100, 1000\}$) grâce au solveur **Monolithique FC** et au schéma convectif centré (`CENTERED`) ;
2. Le maintien de l'intégrité géométrique des deux interfaces opposées de la nappe liquide mince ($\delta_0 = 3\times 10^{-4}\,\mathrm{m}$, soit 6 mailles à $t=0$ sur grille $64\times 64$) lorsqu'elle s'étire sous la taille d'une maille ($\delta_{\mathrm{loc}}(t) < \Delta x$) grâce au **filtrage directionnel des marqueurs CVI** ($\mathbf{n}_p \cdot \mathbf{n}_{i,j}^{\mathrm{ref}} \ge 0.5$).

---

## 2. Théorie : Nappe Liquide Cisaillée et Protection des Films Minces par Filtrage Angulaire

Dans une boîte carrée périodique en $x$ de côté $L = 3\times 10^{-3}\,\mathrm{m}$ ($\Omega = [0, L] \times [0, L]$), une nappe liquide horizontale d'épaisseur $\delta = 3\times 10^{-4}\,\mathrm{m}$ est centrée en $y_c = L/2 = 1.5\times 10^{-3}\,\mathrm{m}$.
Le champ de vitesse initial est analytiquement solénoïdal ($\nabla \cdot \mathbf{u}_0 \equiv 0$) :
$$u(x, y, 0) = A - 0.04 \cos\left(\frac{2\pi x}{L}\right) \left(\frac{L}{2\pi}\right) \left(-\frac{2}{\delta}\right) \exp\left(-\frac{2|y - y_c|}{\delta}\right),$$
$$v(x, y, 0) = 0.04 \sin\left(\frac{2\pi x}{L}\right) \exp\left(-\frac{2|y - y_c|}{\delta}\right),$$
avec $A = 2.0\,\mathrm{m/s}$ dans le liquide et $A = 30.0\,\mathrm{m/s}$ dans le gaz.

### 2.1 Pourquoi le couplage CVI + Marqueurs à filtrage angulaire est décisif sur la nappe liquide
Lorsque l'instabilité de Kelvin-Helmholtz s'amplifie ($t \in [2\times 10^{-4}, 5\times 10^{-4}]\,\mathrm{s}$), la nappe liquide s'affine localement jusqu'à $\delta_{\mathrm{loc}} \le \Delta x$.
- L'interpolateur **CVI** (`CVI_RK2`) sur le champ de vitesse monolithique **FC** ($\|B_u \mathbf{u}^{n+1}\|_\infty \le 10^{-10}$) advecte les deux rideaux de marqueurs (interface supérieure de normale $\mathbf{n}^+ \approx (0, +1)$ et interface inférieure de normale $\mathbf{n}^- \approx (0, -1)$) le long des lignes de courant exactes, y compris à travers la frontière périodique $x=0 \leftrightarrow x=L$.
- Le **filtrage angulaire** ($\mathbf{n}_p \cdot \mathbf{n}_{i,j}^{\mathrm{ref}} \ge 0.5$) garantit que seules les particules appartenant au bon côté de la nappe sont utilisées pour reconstruire la normale PLIC (`marker_normals_plic: true`), empêchant l'annulation mutuelle des deux interfaces opposées.

---

## 3. Matrice des Paramètres et Prototype `.json`

| Paramètre | Spécification |
| :--- | :--- |
| Domaine et Nappe | $L_x = L_y = 0.003\,\mathrm{m}$, $y_c = 0.0015\,\mathrm{m}$, $\delta = 0.0003\,\mathrm{m}$ |
| Rapports de densité | $r_\rho = \rho_2/\rho_1 \in \{100, 1000\}$, $\mu_1 = \mu_2 = 10^{-12}\,\mathrm{Pa\cdot s}$, $\sigma = 0$ |
| Maillage et Temps | $64\times 64$ et $128\times 128$, $\Delta t = 10^{-6}\,\mathrm{s}$, $t_{\mathrm{final}} = 3\times 10^{-4}\,\mathrm{s}$ (et $5\times 10^{-4}\,\mathrm{s}$) |
| Configurations comparées | `MONOLITHIC + WY + PLIC` seul vs `MONOLITHIC + WY + PLIC + CVI_RK2 (marker_normals_plic: true)` |

### Prototype `case/cpp/cas_06_liquid_sheet_fc_cvi_markers_r1000.json`
```json
{
  "mode": "NAVIER_STOKES",
  "test_case": "CAS_10_LIQUID_SHEET",
  "mesh": {
    "nx": 64,
    "ny": 64,
    "lx": 0.003,
    "ly": 0.003
  },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": false,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_svof_smoothing": false,
    "sigma": 0.0,
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": { "rho": 1.0,    "mu": 1e-12 },
    "fluid2": { "rho": 1000.0, "mu": 1e-12 }
  },
  "initial_condition": {
    "type": "LIQUID_SHEET",
    "center_y": 0.0015,
    "thickness": 0.0003
  },
  "boundary_conditions": {
    "left":   { "type": "PERIODIC" },
    "right":  { "type": "PERIODIC" },
    "bottom": { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "top":    { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 }
  },
  "time": {
    "time_scheme": "GEAR",
    "dt": 1e-06,
    "t_final": 0.0003,
    "cfl_max": 0.45,
    "adaptive_dt": false
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.45,
    "use_markers": true,
    "markers_per_cell": 4,
    "marker_advection_scheme": "CVI_RK2",
    "marker_curvature": false,
    "marker_normals_plic": true,
    "marker_max_per_cell": 8,
    "write_markers_vtk": true
  },
  "solver": {
    "type": "MONOLITHIC",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-10,
    "max_iter": 800,
    "r_al": 0.0,
    "uzawa_iter": 1,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "output": {
    "output_dir": "Validation_Totale/Run_CASES/Run_FC-WY-CVI-Marker/06_cas_nappe_liquide_kelvin_helmholtz_cvi_markers/output/cpp/fc_cvi_r1000",
    "case_name": "cas_06_liquid_sheet_cvi_r1000",
    "vtk_interval": 50,
    "log_interval": 25,
    "write_csv_history": true
  }
}
```

---

## 4. Livrables et Diagnostics Attendus

1. **Courbes d'énergie cinétique normalisée $E_k(t)/E_k(0)$** (`figures/conservation_energie_cinetique_nappe.pdf`) pour $r_\rho \in \{100, 1000\}$.
2. **Superposition du contour VOF $C=0.5$ et des deux rideaux de marqueurs CVI** (`figures/nappe_liquide_rideaux_marqueurs_cvi.png`) montrant la préservation du ligament mince et la traversée périodique $x=0 \leftrightarrow x=L$.
3. **Rapport PDF complet (`rapport/rapport_validation_06_nappe_liquide_kelvin_helmholtz_cvi_markers.pdf`)**.
