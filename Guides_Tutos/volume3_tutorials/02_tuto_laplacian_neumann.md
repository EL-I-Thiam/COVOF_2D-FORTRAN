# Tutoriel 2 : Opérateur Laplacien avec Conditions de Neumann Non Homogènes et Jauge de Pression

Ce tutoriel aborde la résolution d'une équation de Laplace sous conditions aux limites de Neumann pur. Il met en lumière deux aspects cruciaux du calcul CFD : le traitement rigoureux des flux pariétaux sur grille structurée et la levée du mode singulier à constante (jauge de pression) pour garantir l'unicité de la solution et la convergence du solveur.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre la compatibilité de flux d'un problème de Neumann pur et l'ancrage de la constante de jauge $\Psi(0,0) = 0$.
- Configurer un fichier `.json` conforme en mode `"LAPLACIAN"` avec `"test_case": "CAS_3_2_NEUMANN"`.
- Exécuter le cas test avec les binaires précompilés `covof2d_cpp` et `covof2d_fortran` via `mpirun`.
- Vérifier dynamiquement la propriété d'exactitude polynomiale machine sur $\Psi_{\mathrm{ex}}(x,y) = x^2 - y^2$.
- Exploiter les fichiers `.csv` et `.vtk` sous Python et ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `CSRMatrix`, `Simulation` (`Simulation::run_laplacian_cas3_2`).
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_simulation` (`sim_run_laplacian_cas3_2`).

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas vérifie le traitement des conditions aux limites de Neumann non homogènes (flux normaux imposés sur $\partial\Omega$) et la régularisation de jauge d'une matrice Laplacienne singulière dont le noyau est engendré par les fonctions constantes.

Sur le domaine carré unitaire $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$, on résout l'équation de Laplace homogène ($S = 0$) :
$$\Delta \Psi = \frac{\partial^2 \Psi}{\partial x^2} + \frac{\partial^2 \Psi}{\partial y^2} = 0 \quad \text{dans } \Omega,$$
avec les conditions aux limites de Neumann suivantes sur les quatre frontières :
$$\left.\frac{\partial \Psi}{\partial n}\right|_{x=0} = -\left.\frac{\partial \Psi}{\partial x}\right|_{x=0} = 0, \qquad \left.\frac{\partial \Psi}{\partial n}\right|_{x=1} = +\left.\frac{\partial \Psi}{\partial x}\right|_{x=1} = +2,$$
$$\left.\frac{\partial \Psi}{\partial n}\right|_{y=0} = -\left.\frac{\partial \Psi}{\partial y}\right|_{y=0} = 0, \qquad \left.\frac{\partial \Psi}{\partial n}\right|_{y=1} = +\left.\frac{\partial \Psi}{\partial y}\right|_{y=1} = -2.$$
La condition de compatibilité $\oint_{\partial\Omega} \frac{\partial\Psi}{\partial n}\,\mathrm{d}\gamma = 0 + 2(1) + 0 - 2(1) = 0$ est exactement satisfaite. En fixant la constante de jauge par la condition d'ancrage $\Psi(0,0) = 0$ (ligne $0$ de la matrice fixée à l'identité dans `run_laplacian_cas3_2` / `sim_run_laplacian_cas3_2`), la solution analytique exacte est le polynôme harmonique quadratique :
$$\Psi_{\mathrm{ex}}(x,y) = x^2 - y^2.$$

### 2.2 Propriété Théorique d'Exactitude Machine
Le développement de Taylor du schéma centré à 5 points montre que l'erreur de troncature locale s'écrit :
$$\tau_{i,j} = \frac{\Delta x^2}{12}\frac{\partial^4 \Psi_{\mathrm{ex}}}{\partial x^4} + \frac{\Delta y^2}{12}\frac{\partial^4 \Psi_{\mathrm{ex}}}{\partial y^4} + \mathcal{O}(\Delta x^4) \equiv 0,$$
car toutes les dérivées d'ordre supérieur ou égal à $3$ de $\Psi_{\mathrm{ex}}(x,y) = x^2 - y^2$ sont identiquement nulles, et le traitement par miroir d'ordre 2 aux nœuds de Neumann ($2/\Delta x^2$ et second membre $\pm 4/\Delta x$) est également exact pour les polynômes quadratiques. L'erreur discrète $E_{L_1}, E_{L_2}, E_{L_\infty}$ doit donc atteindre le zéro machine ($\sim 10^{-14} - 10^{-12}$) quelle que soit la finesse du maillage.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Géométrie | Domaine $[0, L_x] \times [0, L_y]$ | $L_x = 1.0\,\mathrm{m}$, $L_y = 1.0\,\mathrm{m}$ |
| Séquence de grilles | $(N_x \times N_y)$ | $16\times 16$, $32\times 32$, $64\times 64$, $128\times 128$, $256\times 256$ |
| Terme source | $S$ (`source_term`) | $0.0\,\mathrm{m}^{-2}$ |
| Dispatching cas analytique | `test_case` | `"CAS_3_2_NEUMANN"` (active les flux de Neumann exacts et l'ancrage $\Psi(0,0)=0$) |
| Conditions aux limites | `boundary_conditions` | `NEUMANN` sur les 4 bords |
| Solveurs disponibles | `solver.linear_solver` | `DIRECT` (factorisation LU bande de largeur $2N_x + 3$) et `BICGSTAB2` (préconditionneur ligne tridiagonal Thomas) |
| Tolérance et itérations | `solver.tol`, `solver.max_iter` | $\mathrm{tol} = 10^{-14}$, $\mathrm{max\_iter} = 10000$ |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype Solveur Direct (`case/cpp/cas_2_neumann_grid64_direct.json`)
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/grid64_direct"`.)*

```json
{
  "case_name": "cas_2_neumann_grid64_direct",
  "mode": "LAPLACIAN",
  "test_case": "CAS_3_2_NEUMANN",
  "source_term": 0.0,
  "mesh": {
    "nx": 64,
    "ny": 64,
    "lx": 1.0,
    "ly": 1.0
  },
  "boundary_conditions": {
    "left":   { "type": "NEUMANN", "u": 0.0, "v": 0.0 },
    "right":  { "type": "NEUMANN", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "NEUMANN", "u": 0.0, "v": 0.0 },
    "top":    { "type": "NEUMANN", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "DIRECT",
    "linear_solver": "DIRECT",
    "precond_type": "NONE",
    "use_precond": false,
    "tol": 1e-14,
    "max_iter": 10000
  },
  "output": {
    "output_dir": "output/cpp/grid64_direct",
    "case_name": "cas_2_neumann_grid64_direct",
    "vtk_interval": 1,
    "log_interval": 1,
    "write_csv_history": true
  }
}
```

### 4.2 Variante Solveur Itératif `BICGSTAB2` (`case/cpp/cas_2_neumann_grid64_bicgstab2.json`)
Pour activer le solveur itératif `BiCGSTAB(2)` couplé au préconditionneur ligne tridiagonal (algorithme de Thomas ligne par ligne), remplacer le bloc `"solver"` par :

```json
  "solver": {
    "type": "MONOLITHIC",
    "linear_solver": "BICGSTAB2",
    "precond_type": "BLOCK_TRIANGULAR",
    "use_precond": true,
    "tol": 1e-14,
    "max_iter": 10000
  }
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur**. Comme pour le Tutoriel 1, le mode `"LAPLACIAN"` assemble la matrice globale complète sans décomposition par sous-domaines MPI et s'exécute avec `mpirun -np 1`.

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/grid64_direct output/fortran/grid64_direct

# 2. Exécution avec le moteur C++17 (MPI séquentiel)
mpirun -np 1 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_2_neumann_grid64_direct.json

# 3. Exécution avec le moteur Fortran 2018 (MPI séquentiel)
mpirun -np 1 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_2_neumann_grid64_direct.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Exactitude machine sur polynôme quadratique $\Psi_{\mathrm{ex}}(x,y) = x^2 - y^2$ :**
   À partir de `laplacian_error.csv`, vérifier que pour `DIRECT` les erreurs discrètes $E_{L_1}, E_{L_2}, E_{L_\infty}$ restent au niveau de la précision machine ($< 10^{-12}$) sur toutes les grilles testées ($16^2$ à $256^2$).
2. **Convergence du solveur itératif `BICGSTAB2` et intégrité de la jauge :**
   Vérifier dans `solver_history.csv` que le résidu préconditionné décroît de manière monotone sous `tol = 1e-14`, que l'erreur finale satisfait $E_{L_2} \le 10^{-11}$ sans dérive de constante, et que $\Psi(0,0) = 0$.
3. **Antisymétrie diagonale et parité C++17 / Fortran 2018 :**
   Vérifier à partir de `solution_field_2d.csv` que l'antisymétrie $\Psi(y,x) = -\Psi(x,y)$ est respectée à $< 10^{-12}$ près et que l'écart maximal entre C++17 et Fortran 2018 vérifie $\|\Psi^{\mathrm{C++}} - \Psi^{\mathrm{Fortran}}\|_\infty < 10^{-12}$.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `laplacian_error.csv` : `nx,ny,dx,dy,iters,err_l1,err_l2,err_linf,cpu_time`.
- `solver_history.csv` : historique `iter,residual` des itérations de résolution.
- `solution_field_2d.csv` : grille nodale $(N_x+1)\times(N_y+1)$ avec colonnes `i,j,x,y,psi_num,psi_ex,error_abs`.
- `profile_x_mid.csv` et `profile_y_mid.csv` : profils 1D le long de $y = 0.5$ et $x = 0.5$.
- `solution_field.vtk` : fichier VTK `STRUCTURED_GRID` contenant les champs scalaires `psi_numeric`, `psi_exact`, `error_abs`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 2 (Laplacien Neumann) sur fichiers réels."""
import csv
from pathlib import Path

def verify_neumann_case(dir_cpp: str, dir_ftn: str | None = None) -> None:
    with (Path(dir_cpp) / "laplacian_error.csv").open(newline="") as f:
        err = next(csv.DictReader(f))
    print(f"[C++] nx={err['nx']}, iters={err['iters']}, L1={float(err['err_l1']):.3e}, "
          f"L2={float(err['err_l2']):.3e}, Linf={float(err['err_linf']):.3e}")

    pts = {}
    with (Path(dir_cpp) / "solution_field_2d.csv").open(newline="") as f:
        for r in csv.DictReader(f):
            pts[(int(r["i"]), int(r["j"]))] = float(r["psi_num"])
    antisym = max(abs(pts[(i, j)] + pts[(j, i)]) for i, j in pts if (j, i) in pts)
    print(f"[C++] Antisymétrie max |Psi(x,y)+Psi(y,x)| = {antisym:.3e} | Psi(0,0) = {pts[(0,0)]:.3e}")

    if dir_ftn and (Path(dir_ftn) / "solution_field_2d.csv").exists():
        pts_f = {}
        with (Path(dir_ftn) / "solution_field_2d.csv").open(newline="") as f:
            for r in csv.DictReader(f):
                pts_f[(int(r["i"]), int(r["j"]))] = float(r["psi_num"])
        parity = max(abs(pts[k] - pts_f[k]) for k in pts)
        print(f"[Parité C++/Fortran] Linf = {parity:.3e}")

if __name__ == "__main__":
    verify_neumann_case("output/cpp/grid64_direct", "output/fortran/grid64_direct")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir `output/cpp/grid64_direct/solution_field.vtk` avec **File -> Open** puis **Apply**.
2. Afficher `psi_numeric` en palette divergente (`Cool to Warm`) pour mettre en évidence la selle harmonique $x^2 - y^2$ (positif vers $(1,0)$, négatif vers $(0,1)$, nul sur la diagonale $x=y$).
3. Appliquer le filtre **Contour** sur `psi_numeric` avec des isovaleurs régulières et vérifier l'orthogonalité des isolignes vis-à-vis des frontières à flux nul ($x=0$ et $y=0$).
