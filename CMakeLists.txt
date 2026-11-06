cmake_minimum_required(VERSION 3.18)
project(COVOF_2D_FORTRAN LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native")

set(SOURCES
    src/module_precision.f90
    src/module_constants.f90
    src/module_eos.f90
    src/module_json.f90
    src/module_config.f90
    src/module_mesh.f90
    src/module_fields.f90
    src/module_mpi_comm.f90
    src/module_interface_reconstruction.f90
    src/module_interface_vof.f90
    src/module_navier_stokes.f90
    src/module_linear_solvers.f90
    src/module_vtk_writer.f90
    src/module_simulation.f90
    src/covof2d_main.f90
)

# Executable target: covof2d_fortran
add_executable(covof2d_fortran ${SOURCES})
if(APPLE)
    target_link_libraries(covof2d_fortran PRIVATE MPI::MPI_Fortran "-framework Accelerate")
else()
    target_link_libraries(covof2d_fortran PRIVATE MPI::MPI_Fortran)
endif()

set_target_properties(covof2d_fortran PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)
