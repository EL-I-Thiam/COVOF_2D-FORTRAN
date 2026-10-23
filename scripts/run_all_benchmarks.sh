#!/bin/bash
# ==============================================================================
# COVOF 2D (Fortran 2018 + OpenMPI) - Lancement de la Suite Complète de Validation
# Couvre TOUS les cas des Chapitres 3, 4, 6 et 7 du Cahier de Validation et Vérification CFD 2D
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${ROOT_DIR}"

EXEC="${ROOT_DIR}/bin/covof2d_fortran"
if [ ! -f "${EXEC}" ]; then
    echo "Compilation de COVOF 2D Fortran 2018..."
    mkdir -p bin
    mpifort -O3 src/module_precision.f90 src/module_constants.f90 src/module_json.f90 \
            src/module_config.f90 src/module_mesh.f90 src/module_fields.f90 \
            src/module_mpi_comm.f90 src/module_interface_vof.f90 \
            src/module_navier_stokes.f90 src/module_linear_solvers.f90 \
            src/module_vtk_writer.f90 src/module_simulation.f90 src/covof2d_main.f90 \
            -o "${EXEC}"
fi

NP=2
if [ -n "$1" ]; then
    NP=$1
fi

echo "================================================================================"
echo "          LANCEMENT DE LA SUITE DE BENCHMARKS COVOF 2D FORTRAN 2018              "
echo "          Processus MPI : ${NP}                                                 "
echo "================================================================================"

echo ""
echo ">>> [1/13] Cas 3.1 : Laplacien Cartesien Delta(Psi) = 1 <<<"
mpirun -np ${NP} "${EXEC}" cases/laplacian.json

echo ""
echo ">>> [2/13] Cas 3.7 : Poiseuille Plan Stationnaire <<<"
mpirun -np ${NP} "${EXEC}" cases/poiseuille.json

echo ""
echo ">>> [3/13] Cas 3.8 : Poiseuille Instationnaire (Szymanski) <<<"
mpirun -np ${NP} "${EXEC}" cases/poiseuille_unsteady.json

echo ""
echo ">>> [4/13] Cas 3.9 : Poiseuille a Viscosite Variable Exponentielle <<<"
mpirun -np ${NP} "${EXEC}" cases/poiseuille_variable_viscosity.json

echo ""
echo ">>> [5/13] Cas 3.19 : Advection VOF Vortex Cisaillement Reversible <<<"
echo "    -> [5.1/13] VOF-TVD Superbee (Reference)..."
rm -rf output_vortex_shear/*
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_advection.json

echo "    -> [5.2/13] Geometrique CIAM + PLIC..."
rm -rf output_vortex_shear_ciam_plic/*
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_ciam_plic.json

echo "    -> [5.3/13] Multidimensionnel Weymouth & Yue + ELVIRA..."
rm -rf output_vortex_shear_weymouth_yue_elvira/*
mpirun -np ${NP} "${EXEC}" cases/vortex_shear_weymouth_yue_elvira.json

echo ""
echo ">>> [6/13] Cas 4.1 : Marche Descendante (Armaly Re=100) <<<"
mpirun -np ${NP} "${EXEC}" cases/backward_facing_step.json

echo ""
echo ">>> [7/13] Cas 4.4 : Ecoulement Autour d'un Cylindre (Re=20) <<<"
mpirun -np ${NP} "${EXEC}" cases/obstacle_cylinder_re20.json

echo ""
echo ">>> [8/13] Cas 4.5 : Cavite Entrainee (Ghia Re=1000) <<<"
mpirun -np ${NP} "${EXEC}" cases/cavity_re1000.json

echo ""
echo ">>> [9/13] Cas 4.7 : Cavite en Forme de L (Re=1000) <<<"
mpirun -np ${NP} "${EXEC}" cases/cavity_l_shaped.json

echo ""
echo ">>> [10/13] Cas 6.1 : Convection Naturelle Differentiellement Chauffee (Ra=1e4) <<<"
mpirun -np ${NP} "${EXEC}" cases/natural_convection_ra1e4.json

echo ""
echo ">>> [11/13] Cas 7.1 : Rupture de Barrage (Dam Break Martin & Moyce) <<<"
mpirun -np ${NP} "${EXEC}" cases/dam_break.json

echo ""
echo ">>> [12/13] Cas 7.4 : Goutte Statique (Loi de Laplace) <<<"
mpirun -np ${NP} "${EXEC}" cases/static_droplet.json

echo ""
echo ">>> [13/13] Cas 7.5 : Bulle Montante (Benchmark Hysing) <<<"
mpirun -np ${NP} "${EXEC}" cases/rising_bubble.json

echo ""
echo "================================================================================"
echo "          EXECUTION DU SCRIPT DE VERIFICATION ANALYTIQUE PYTHON                 "
echo "================================================================================"
python3 scripts/verify_cases.py

echo ""
echo "================================================================================"
echo "          GENERATION DES FIGURES COMPARATIVES DU Cahier de Validation de Référence                 "
echo "================================================================================"
python3 scripts/generate_validation_figures.py

echo ""
echo ">>> SUITE DE BENCHMARKS ET VALIDATIONS TERMINEE AVEC SUCCES <<<"
