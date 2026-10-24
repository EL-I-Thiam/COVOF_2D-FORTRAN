#!/bin/bash
# ==============================================================================
# COVOF 2D (C++17 + OpenMPI) - Lancement de la Suite Complète de Validation
# Couvre TOUS les cas des Chapitres 3, 4, 6, 7 et les 3 Benchmarks Fluids 2021
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${ROOT_DIR}"

EXEC="${ROOT_DIR}/bin/covof2d"
if [ ! -f "${EXEC}" ]; then
    echo "Compilation de COVOF 2D C++17..."
    cmake -B build && cmake --build build -j4
fi

NP=2
if [ -n "$1" ]; then
    NP=$1
fi

echo "================================================================================"
echo "          LANCEMENT DE LA SUITE DE BENCHMARKS COVOF 2D C++17                     "
echo "          Processus MPI : ${NP}                                                 "
echo "================================================================================"

echo ""
echo ">>> [1/16] Cas 3.1 : Laplacien Cartesien Delta(Psi) = 1 <<<"
mpirun -np ${NP} "${EXEC}" cases/laplacian.json

echo ""
echo ">>> [2/16] Cas 3.7 : Poiseuille Plan Stationnaire <<<"
mpirun -np ${NP} "${EXEC}" cases/poiseuille.json

echo ""
echo ">>> [3/16] Cas 3.8 : Poiseuille Instationnaire (Szymanski) <<<"
mpirun -np ${NP} "${EXEC}" cases/poiseuille_unsteady.json

echo ""
echo ">>> [4/16] Cas 3.9 : Poiseuille a Viscosite Variable Exponentielle <<<"
mpirun -np ${NP} "${EXEC}" cases/poiseuille_variable_viscosity.json

echo ""
echo ">>> [5/16] Cas 3.19 : Advection VOF Vortex Cisaillement Reversible <<<"
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_advection.json
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_ciam_plic.json
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_ciam_slic.json
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_ciam_elvira.json
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_ciam_lvira.json
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_weymouth_yue_plic.json
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_weymouth_yue_slic.json
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_weymouth_yue_elvira.json
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_weymouth_yue_lvira.json
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_vof_tvd_hric.json
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_vof_tvd_muscl.json
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_vof_tvd_vanleer.json

echo ""
echo ">>> [6/16] Cas 4.1 : Marche Descendante (Armaly Re=100) <<<"
mpirun -np ${NP} "${EXEC}" cases/backward_facing_step.json

echo ""
echo ">>> [7/16] Cas 4.4 : Ecoulement Autour d'un Cylindre (Re=20) <<<"
mpirun -np ${NP} "${EXEC}" cases/obstacle_cylinder_re20.json

echo ""
echo ">>> [8/16] Cas 4.5 : Cavite Entrainee (Ghia Re=1000) <<<"
mpirun -np ${NP} "${EXEC}" cases/cavity_re1000.json

echo ""
echo ">>> [9/16] Cas 4.7 : Cavite en Forme de L (Re=1000) <<<"
mpirun -np ${NP} "${EXEC}" cases/cavity_l_shaped.json

echo ""
echo ">>> [10/16] Cas 6.1 : Convection Naturelle Differentiellement Chauffee (Ra=1e4) <<<"
mpirun -np ${NP} "${EXEC}" cases/natural_convection_ra1e4.json

echo ""
echo ">>> [11/16] Cas 7.1 : Rupture de Barrage (Dam Break Martin & Moyce) <<<"
mpirun -np ${NP} "${EXEC}" cases/dam_break.json

echo ""
echo ">>> [12/16] Cas 7.4 : Goutte Statique (Loi de Laplace) <<<"
mpirun -np ${NP} "${EXEC}" cases/static_droplet.json

echo ""
echo ">>> [13/16] Cas 7.5 : Bulle Montante (Benchmark Hysing) <<<"
mpirun -np ${NP} "${EXEC}" cases/rising_bubble.json

echo ""
echo ">>> [14/16] Cas 7.10 : Instabilite de Rayleigh-Taylor Diphasique <<<"
mpirun -np ${NP} "${EXEC}" cases/rayleigh_taylor.json

echo ""
echo ">>> [15/16] Fluids 2021 Cas 3.1 & 3.2 : Cylindre Dense & Nappe Liquide <<<"
mpirun -np ${NP} "${EXEC}" cases/free_fall_dense_cylinder.json
mpirun -np ${NP} "${EXEC}" cases/liquid_sheet_ratio100.json

echo ""
echo ">>> [16/16] Fluids 2021 Cas 3.3 : Flambement de Jet Visqueux (3 Regimes Re) <<<"
mpirun -np ${NP} "${EXEC}" cases/viscous_jet_case1.json
mpirun -np ${NP} "${EXEC}" cases/viscous_jet_full.json
mpirun -np ${NP} "${EXEC}" cases/viscous_jet_case2.json
mpirun -np ${NP} "${EXEC}" cases/viscous_jet_case3.json

echo ""
echo "================================================================================"
echo "          EXECUTION DES SCRIPTS DE VERIFICATION ANALYTIQUE PYTHON               "
echo "================================================================================"
python3 scripts/verify_cases.py
python3 scripts/verify_dense_cylinder.py
python3 scripts/verify_liquid_sheet.py
python3 scripts/verify_viscous_jet.py
