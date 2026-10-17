cmake_minimum_required(VERSION 3.18)
project(fugu2d_fortran LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native")

set(SOURCES
    src/module_precision.f90
    src/module_constants.f90
    src/module_json.f90
    src/module_config.f90
    src/module_mesh.f90
    src/module_fields.f90
    src/module_mpi_comm.f90
    src/module_interface_vof.f90
    src/module_navier_stokes.f90
    src/module_linear_solvers.f90
    src/module_vtk_writer.f90
    src/module_simulation.f90
    src/fugu2d_main.f90
)

add_executable(fugu2d_fortran ${SOURCES})

target_link_libraries(fugu2d_fortran PRIVATE MPI::MPI_Fortran)

set_target_properties(fugu2d_fortran PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
