#!/usr/bin/env python3
import os
import sys
import json
import math
import random
import shutil
import subprocess
import xml.etree.ElementTree as ET
import numpy as np

ROOT = "/Users/elhadjiibrahimathiam/Documents/PROJET_COVOF-2D"
BIN_FORTRAN = os.path.join(ROOT, "COVOF_2D-FORTRAN/bin/covof2d_fortran")
BIN_CPP = os.path.join(ROOT, "COVOF_2D-C++/bin/covof2d_cpp")
WORK_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "build_verify_runs")

os.makedirs(WORK_DIR, exist_ok=True)


def run_cmd(cmd, cwd=WORK_DIR):
    res = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if res.returncode != 0:
        print(f"ERROR running {' '.join(cmd)}:\n{res.stdout}")
        sys.exit(1)
    return res.stdout


def parse_vtr_global_fields(out_dir, case_name, step, nx_global, ny_global, np_procs):
    """Assemble global 2D cell scalar fields (ny_global, nx_global) from MPI .vtr pieces."""
    fields = {}
    for r in range(np_procs):
        vtr_path = os.path.join(out_dir, f"{case_name}_p{r}_{step:06d}.vtr")
        tree = ET.parse(vtr_path)
        root = tree.getroot()
        piece = root.find(".//Piece")
        ext = [int(x) for x in piece.attrib["Extent"].split()]
        x0, x1, y0, y1 = ext[0], ext[1], ext[2], ext[3]
        nx_loc = x1 - x0
        ny_loc = y1 - y0

        cell_data = piece.find("CellData")
        for da in cell_data.findall("DataArray"):
            name = da.attrib["Name"]
            n_comp = int(da.attrib.get("NumberOfComponents", "1"))
            if n_comp != 1:
                continue
            if name not in fields:
                fields[name] = np.zeros((ny_global, nx_global), dtype=float)
            vals = np.fromstring(da.text, sep=" ", dtype=float)
            fields[name][y0:y1, x0:x1] = vals.reshape((ny_loc, nx_loc))
    return fields


def parse_markers_vtk(vtk_path):
    """Parse points and curvatures from *_markers_*.vtk POLYDATA."""
    with open(vtk_path, "r") as f:
        lines = f.readlines()
    pts = []
    kappas = []
    i = 0
    n_pts = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("POINTS"):
            n_pts = int(line.split()[1])
            i += 1
            while len(pts) < n_pts and i < len(lines):
                row = [float(x) for x in lines[i].split()]
                if len(row) >= 2:
                    pts.append((row[0], row[1]))
                i += 1
        elif line.startswith("SCALARS MarkerCurvature"):
            i += 2  # skip LOOKUP_TABLE
            while len(kappas) < n_pts and i < len(lines):
                kappas.extend([float(x) for x in lines[i].split()])
                i += 1
        else:
            i += 1
    return np.array(pts), np.array(kappas[:n_pts])


def test_palier_1_cvi_divergence():
    print("=" * 72)
    print("PALIER 1 : Verification mathematique CVI Variante B (Wang 2015 / Pusok 2017)")
    print("=" * 72)
    random.seed(42)
    dx, dy = 0.03125, 0.025
    max_div_err_cvi_analytic = 0.0
    max_div_err_cvi_fd = 0.0
    max_div_err_bilin = 0.0
    max_div_err_var_a = 0.0

    for _ in range(5000):
        # Construire un bloc 3x3 de cellules MAC (i-1..i+1, j-1..j+1) a divergence discrete exacte D_{k,l} = 0
        # u de taille (4, 3) pour i in 0..3 (faces), j in 0..2 (cellules)
        # v de taille (3, 4) pour i in 0..2 (cellules), j in 0..3 (faces)
        u = np.random.uniform(-2.0, 2.0, size=(4, 3))
        v = np.random.uniform(-2.0, 2.0, size=(3, 4))
        for kl_i in range(3):
            for kl_j in range(3):
                # Imposer D_{kl_i, kl_j} = (u[kl_i+1, kl_j] - u[kl_i, kl_j])/dx + (v[kl_i, kl_j+1] - v[kl_i, kl_j])/dy == 0
                v[kl_i, kl_j + 1] = v[kl_i, kl_j] - (dy / dx) * (u[kl_i + 1, kl_j] - u[kl_i, kl_j])

        # Cellule centrale (i=1, j=1) : tester un point aleatoire dans l'un des 4 sous-quadrants
        # (en evitant de chevaucher les lignes medianes x=x_c, y=y_c pour le stencil DF eps=1e-7)
        quad_x = random.choice([0, 1])
        quad_y = random.choice([0, 1])
        xi = random.uniform(0.05, 0.45) if quad_x == 0 else random.uniform(0.55, 0.95)
        eta = random.uniform(0.05, 0.45) if quad_y == 0 else random.uniform(0.55, 0.95)

        def eval_var_b(xi_loc, eta_loc, use_cvi=True):
            i, j = 1, 1
            j_u = j - 1 if eta_loc < 0.5 else j
            i_v = i - 1 if xi_loc < 0.5 else i
            xi_u = xi_loc
            eta_u = eta_loc + 0.5 if j_u == j - 1 else eta_loc - 0.5
            xi_v = xi_loc + 0.5 if i_v == i - 1 else xi_loc - 0.5
            eta_v = eta_loc

            u_sw, u_se = u[i, j_u], u[i + 1, j_u]
            u_nw, u_ne = u[i, j_u + 1], u[i + 1, j_u + 1]
            v_sw, v_se = v[i_v, j], v[i_v + 1, j]
            v_nw, v_ne = v[i_v, j + 1], v[i_v + 1, j + 1]

            du_s = u_se - u_sw
            du_n = u_ne - u_nw
            dv_w = v_nw - v_sw
            dv_e = v_ne - v_se

            u_b = (1.0 - eta_u) * ((1.0 - xi_u) * u_sw + xi_u * u_se) + eta_u * ((1.0 - xi_u) * u_nw + xi_u * u_ne)
            v_b = (1.0 - xi_v) * ((1.0 - eta_v) * v_sw + eta_v * v_nw) + xi_v * ((1.0 - eta_v) * v_se + eta_v * v_ne)
            if not use_cvi:
                return u_b, v_b
            u_c = u_b + 0.5 * (dx / dy) * xi_u * (1.0 - xi_u) * (dv_e - dv_w)
            v_c = v_b + 0.5 * (dy / dx) * eta_v * (1.0 - eta_v) * (du_n - du_s)
            return u_c, v_c

        # Divergence analytique de la Variante B
        i, j = 1, 1
        j_u = j - 1 if eta < 0.5 else j
        i_v = i - 1 if xi < 0.5 else i
        xi_u = xi
        eta_u = eta + 0.5 if j_u == j - 1 else eta - 0.5
        xi_v = xi + 0.5 if i_v == i - 1 else xi - 0.5
        eta_v = eta
        du_s = u[i + 1, j_u] - u[i, j_u]
        du_n = u[i + 1, j_u + 1] - u[i, j_u + 1]
        dv_w = v[i_v, j + 1] - v[i_v, j]
        dv_e = v[i_v + 1, j + 1] - v[i_v + 1, j]

        du_bilin_dx = ((1.0 - eta_u) * du_s + eta_u * du_n) / dx
        dv_bilin_dy = ((1.0 - xi_v) * dv_w + xi_v * dv_e) / dy
        div_bilin = du_bilin_dx + dv_bilin_dy

        dcorr_u_dx = (0.5 / dy) * (1.0 - 2.0 * xi_u) * (dv_e - dv_w)
        dcorr_v_dy = (0.5 / dx) * (1.0 - 2.0 * eta_v) * (du_n - du_s)
        div_cvi_analytic = div_bilin + dcorr_u_dx + dcorr_v_dy

        # Verification par differences finies centrees (sans erreur de troncature sur polynome quadratique)
        eps_x = 1.0e-3 * dx
        eps_y = 1.0e-3 * dy
        u_xp, _ = eval_var_b(xi + eps_x / dx, eta, True)
        u_xm, _ = eval_var_b(xi - eps_x / dx, eta, True)
        _, v_yp = eval_var_b(xi, eta + eps_y / dy, True)
        _, v_ym = eval_var_b(xi, eta - eps_y / dy, True)
        div_cvi_fd = (u_xp - u_xm) / (2.0 * eps_x) + (v_yp - v_ym) / (2.0 * eps_y)

        # Comparaison avec Variante A (moyenne de 2 faces adjacentes aux coins)
        u_sw_a = 0.5 * (u[i, j - 1] + u[i, j])
        u_se_a = 0.5 * (u[i + 1, j - 1] + u[i + 1, j])
        u_nw_a = 0.5 * (u[i, j] + u[i, j + 1])
        u_ne_a = 0.5 * (u[i + 1, j] + u[i + 1, j + 1])
        v_sw_a = 0.5 * (v[i - 1, j] + v[i, j])
        v_se_a = 0.5 * (v[i, j] + v[i + 1, j])
        v_nw_a = 0.5 * (v[i - 1, j + 1] + v[i, j + 1])
        v_ne_a = 0.5 * (v[i, j + 1] + v[i + 1, j + 1])
        div_var_a = ((u_se_a - u_sw_a) + (u_ne_a - u_nw_a)) / (2.0 * dx) + ((v_nw_a - v_sw_a) + (v_ne_a - v_se_a)) / (2.0 * dy)

        max_div_err_cvi_analytic = max(max_div_err_cvi_analytic, abs(div_cvi_analytic))
        max_div_err_cvi_fd = max(max_div_err_cvi_fd, abs(div_cvi_fd))
        max_div_err_bilin = max(max_div_err_bilin, abs(div_bilin))
        max_div_err_var_a = max(max_div_err_var_a, abs(div_var_a))

    print(f"  Erreur max |div(u_Bilineaire)|          : {max_div_err_bilin:.6e} s^-1 (Violation O(h))")
    print(f"  Erreur max |div(u_CVI Variante A)|      : {max_div_err_var_a:.6e} s^-1 (Non nulle : moyenne 2 faces)")
    print(f"  Erreur max |div(u_CVI Variante B)| (DF) : {max_div_err_cvi_fd:.6e} s^-1 (Exact a < 1e-11)")
    print(f"  Erreur max |div(u_CVI Variante B)| (An) : {max_div_err_cvi_analytic:.6e} s^-1 (Zero machine exact)")
    assert max_div_err_cvi_analytic < 1.0e-13, "CVI Variante B analytical divergence is not machine zero!"
    assert max_div_err_cvi_fd < 1.0e-10, "CVI Variante B finite-difference divergence exceeds tolerance!"


def run_young_laplace_case(bin_path, engine_name, nx, use_markers, enable_svof=False, np_procs=1):
    svof_tag = "_SVOF" if enable_svof else ""
    tag = f"{engine_name}_YL_N{nx}_{'MARKERS' if use_markers else 'CSF_STD'}{svof_tag}_np{np_procs}"
    out_dir = os.path.join(WORK_DIR, tag)
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    cfg = {
        "mode": "NAVIER_STOKES",
        "mesh": {"nx": nx, "ny": nx, "lx": 1.0, "ly": 1.0},
        "physics": {
            "enable_two_phase": True,
            "enable_surface_tension": True,
            "enable_convection": False,
            "enable_energy": False,
            "enable_svof_smoothing": enable_svof,
            "svof_filter_width": 0.5,
            "enable_diffusion_smoothing": False,
            "sigma": 1.0,
            "gx": 0.0,
            "gy": 0.0,
            "fluid1": {"rho": 1.0, "mu": 0.1},
            "fluid2": {"rho": 1.0, "mu": 0.1}
        },
        "initial_condition": {
            "type": "DROPLET",
            "center_x": 0.5,
            "center_y": 0.5,
            "radius": 0.25
        },
        "boundary_conditions": {
            "left":   {"type": "DIRICHLET", "u": 0.0, "v": 0.0},
            "right":  {"type": "DIRICHLET", "u": 0.0, "v": 0.0},
            "bottom": {"type": "DIRICHLET", "u": 0.0, "v": 0.0},
            "top":    {"type": "DIRICHLET", "u": 0.0, "v": 0.0}
        },
        "vof": {
            "advection_method": "WEYMOUTH_YUE",
            "reconstruction_method": "ELVIRA",
            "cfl_vof": 0.25,
            "use_markers": use_markers,
            "markers_per_cell": 4,
            "marker_advection_scheme": "CVI_RK2",
            "marker_curvature": True,
            "marker_normals_plic": True,
            "write_markers_vtk": True
        },
        "time": {"dt": 1.0e-4, "t_final": 5.0e-3, "adaptive_dt": False},
        "solver": {
            "type": "AUGMENTED_LAGRANGIAN",
            "linear_solver": "BICGSTAB2",
            "precond_type": "JACOBI",
            "tol": 1.0e-10,
            "max_iter": 1000,
            "r_al": 100.0,
            "uzawa_iter": 2
        },
        "output": {
            "output_dir": tag,
            "case_name": "yl",
            "vtk_interval": 50,
            "log_interval": 25,
            "write_csv_history": True
        }
    }

    cfg_path = os.path.join(out_dir, "config.json")
    with open(cfg_path, "w") as f:
        json.dump(cfg, f, indent=2)

    cfg_rel = os.path.join(tag, "config.json")
    cmd = ["mpirun", "-np", str(np_procs), bin_path, cfg_rel] if np_procs > 1 else [bin_path, cfg_rel]
    run_cmd(cmd)

    # Parse final step VTR (step 50)
    fields = parse_vtr_global_fields(out_dir, "yl", 50, nx, nx, np_procs)
    C = fields["PhaseIndicator"]
    P = fields["Pressure"]
    K = fields["Curvature"]

    # Pressure jump: inside (r < 0.15) vs outside (r > 0.38)
    dx = 1.0 / nx
    xs = (np.arange(nx) + 0.5) * dx
    ys = (np.arange(nx) + 0.5) * dx
    XX, YY = np.meshgrid(xs, ys)
    RR = np.hypot(XX - 0.5, YY - 0.5)

    p_in = np.mean(P[RR < 0.15])
    p_out = np.mean(P[RR > 0.38])
    dp_num = p_in - p_out
    dp_err_rel = abs(dp_num - 4.0) / 4.0

    # Interface curvature in mixed cells (0.05 < C < 0.95)
    mixed = (C > 0.05) & (C < 0.95)
    k_mixed = np.abs(K[mixed])
    k_mean = np.mean(k_mixed)
    k_linf_err = np.max(np.abs(k_mixed - 4.0)) / 4.0
    k_l2_err = np.sqrt(np.mean((k_mixed - 4.0) ** 2)) / 4.0

    # Parse max velocity from yl_history.csv
    diag_file = os.path.join(out_dir, "yl_history.csv")
    with open(diag_file, "r") as f:
        dlines = [l.strip() for l in f.readlines() if l.strip()]
    last_row = [float(x) for x in dlines[-1].split(",")]
    u_max = last_row[7]

    n_markers = 0
    if use_markers:
        m_vtk = os.path.join(out_dir, "yl_markers_000050.vtk")
        pts, mk_kap = parse_markers_vtk(m_vtk)
        n_markers = len(pts)

    return {
        "dp_num": dp_num,
        "dp_err_rel": dp_err_rel,
        "k_mean": k_mean,
        "k_l2_err": k_l2_err,
        "k_linf_err": k_linf_err,
        "u_max": u_max,
        "n_markers": n_markers
    }


def run_vortex_case(bin_path, engine_name, nx, use_markers, np_procs=1):
    tag = f"{engine_name}_VORTEX_N{nx}_{'MARKERS' if use_markers else 'WY_STD'}_np{np_procs}"
    out_dir = os.path.join(WORK_DIR, tag)
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    cfg = {
        "mode": "VORTEX_ADVECTION",
        "mesh": {"nx": nx, "ny": nx, "lx": 1.0, "ly": 1.0},
        "physics": {"enable_two_phase": True},
        "initial_condition": {
            "type": "VORTEX_SHEAR",
            "center_x": 0.5,
            "center_y": 0.75,
            "radius": 0.15
        },
        "vof": {
            "advection_method": "WEYMOUTH_YUE",
            "reconstruction_method": "ELVIRA",
            "use_markers": use_markers,
            "markers_per_cell": 4,
            "marker_advection_scheme": "CVI_RK2",
            "marker_curvature": True,
            "marker_normals_plic": True,
            "write_markers_vtk": True
        },
        "time": {"dt": 0.005, "t_final": 1.0},
        "output": {
            "output_dir": tag,
            "case_name": "vortex",
            "vtk_interval": 100,
            "log_interval": 100
        }
    }
    cfg_path = os.path.join(out_dir, "config.json")
    with open(cfg_path, "w") as f:
        json.dump(cfg, f, indent=2)

    cfg_rel = os.path.join(tag, "config.json")
    cmd = ["mpirun", "-np", str(np_procs), bin_path, cfg_rel] if np_procs > 1 else [bin_path, cfg_rel]
    run_cmd(cmd)

    err_csv = os.path.join(out_dir, "vortex_error.csv")
    with open(err_csv, "r") as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]
    vals = lines[-1].split(",")
    err_l1 = float(vals[7])
    err_l2 = float(vals[8])
    mass_init = float(vals[10])
    mass_final = float(vals[11])
    mass_err = abs(mass_final - mass_init) / mass_init

    n_markers = 0
    r_err_mean = 0.0
    if use_markers:
        m_vtk = os.path.join(out_dir, "vortex_markers_000200.vtk")
        pts, _ = parse_markers_vtk(m_vtk)
        n_markers = len(pts)
        # At t = T = 1.0 s, the vortex has reversed and markers should lie on the circle (0.5, 0.75), R=0.15
        radii = np.hypot(pts[:, 0] - 0.5, pts[:, 1] - 0.75)
        r_err_mean = np.mean(np.abs(radii - 0.15))

    return {
        "err_l1": err_l1,
        "err_l2": err_l2,
        "mass_err": mass_err,
        "n_markers": n_markers,
        "r_err_mean": r_err_mean
    }


def main():
    test_palier_1_cvi_divergence()

    print("\n" + "=" * 72)
    print("PALIER 2 : Goutte statique de Young-Laplace (R=0.25, kappa_ex=4.0, dP_ex=4.0)")
    print("=" * 72)
    for eng, bpath in [("FORTRAN", BIN_FORTRAN), ("CPP", BIN_CPP)]:
        for nx in [32, 64]:
            r_std = run_young_laplace_case(bpath, eng, nx, use_markers=False, enable_svof=False, np_procs=1)
            r_svof = run_young_laplace_case(bpath, eng, nx, use_markers=False, enable_svof=True, np_procs=1)
            r_hyb = run_young_laplace_case(bpath, eng, nx, use_markers=True, enable_svof=False, np_procs=1)
            r_hyb_svof = run_young_laplace_case(bpath, eng, nx, use_markers=True, enable_svof=True, np_procs=1)
            print(f"[{eng} {nx}x{nx}]")
            print(f"  CSF Brut (sans lissage)          : <kappa>={r_std['k_mean']:.4f} (Err L2={r_std['k_l2_err']*100:.2f}%, Linf={r_std['k_linf_err']*100:.2f}%) | dP={r_std['dp_num']:.4f} (Err={r_std['dp_err_rel']*100:.2f}%) | ||u||_inf={r_std['u_max']:.3e}")
            print(f"  CSF + Lissage SVOF (k=0.5)       : <kappa>={r_svof['k_mean']:.4f} (Err L2={r_svof['k_l2_err']*100:.2f}%, Linf={r_svof['k_linf_err']*100:.2f}%) | dP={r_svof['dp_num']:.4f} (Err={r_svof['dp_err_rel']*100:.2f}%) | ||u||_inf={r_svof['u_max']:.3e}")
            print(f"  Hybride CVI+Marqueurs            : <kappa>={r_hyb['k_mean']:.4f} (Err L2={r_hyb['k_l2_err']*100:.2f}%, Linf={r_hyb['k_linf_err']*100:.2f}%) | dP={r_hyb['dp_num']:.4f} (Err={r_hyb['dp_err_rel']*100:.2f}%) | ||u||_inf={r_hyb['u_max']:.3e} | Marqueurs={r_hyb['n_markers']}")
            print(f"  Hybride CVI+Marqueurs + SVOF     : <kappa>={r_hyb_svof['k_mean']:.4f} (Err L2={r_hyb_svof['k_l2_err']*100:.2f}%, Linf={r_hyb_svof['k_linf_err']*100:.2f}%) | dP={r_hyb_svof['dp_num']:.4f} (Err={r_hyb_svof['dp_err_rel']*100:.2f}%) | ||u||_inf={r_hyb_svof['u_max']:.3e} | Marqueurs={r_hyb_svof['n_markers']}")
            assert r_hyb["k_l2_err"] < r_std["k_l2_err"], f"{eng}: Marker curvature L2 error should be smaller than unsmoothed CSF!"
            assert r_hyb["u_max"] < r_std["u_max"], f"{eng}: Spurious currents should be reduced with marker curvature!"

    print("\n" + "=" * 72)
    print("PALIER 3 : Tourbillon cisaille de Rider-Kothe & Coherence MPI (np=1 vs np=4)")
    print("=" * 72)
    for eng, bpath in [("FORTRAN", BIN_FORTRAN), ("CPP", BIN_CPP)]:
        v_np1 = run_vortex_case(bpath, eng, 64, use_markers=True, np_procs=1)
        v_np4 = run_vortex_case(bpath, eng, 64, use_markers=True, np_procs=4)
        print(f"[{eng} Vortex 64x64 T=1.0s]")
        print(f"  np=1 : L1(C)={v_np1['err_l1']:.4e} | |dM|/M0={v_np1['mass_err']:.3e} | N_mk={v_np1['n_markers']} | Err_rayon_marqueurs={v_np1['r_err_mean']:.4e}")
        print(f"  np=4 : L1(C)={v_np4['err_l1']:.4e} | |dM|/M0={v_np4['mass_err']:.3e} | N_mk={v_np4['n_markers']} | Err_rayon_marqueurs={v_np4['r_err_mean']:.4e}")

        yl_np1 = run_young_laplace_case(bpath, eng, 32, use_markers=True, np_procs=1)
        yl_np4 = run_young_laplace_case(bpath, eng, 32, use_markers=True, np_procs=4)
        print(f"[{eng} Young-Laplace 32x32 MPI np=1 vs np=4 (Interface traversant les 4 sous-domaines en (0.5,0.5))]")
        print(f"  np=1 : <kappa>={yl_np1['k_mean']:.6f} | dP={yl_np1['dp_num']:.6f} | ||u||_inf={yl_np1['u_max']:.4e} | N_mk={yl_np1['n_markers']}")
        print(f"  np=4 : <kappa>={yl_np4['k_mean']:.6f} | dP={yl_np4['dp_num']:.6f} | ||u||_inf={yl_np4['u_max']:.4e} | N_mk={yl_np4['n_markers']}")
        assert abs(yl_np1["k_mean"] - yl_np4["k_mean"]) < 1.0e-5, f"{eng}: MPI np=1 vs np=4 curvature mismatch!"
        assert yl_np1["n_markers"] == yl_np4["n_markers"], f"{eng}: MPI np=1 vs np=4 marker count mismatch!"

    print("\n>>> TOUS LES PALIERS DE VERIFICATION (FORTRAN & C++, SEQ & MPI) SONT VALIDES AVEC SUCCES ! <<<")


if __name__ == "__main__":
    main()
