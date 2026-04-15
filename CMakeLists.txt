cmake_minimum_required(VERSION 3.10)
project(packmol LANGUAGES Fortran CXX)

include(CTest)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    add_compile_options(
        -Wall "$<$<CONFIG:Debug>:-Werror>"
    )
endif()

set(PACKMOL_CORE_SOURCES
       src/exit_codes.f90
       src/cenmass.f90
       src/gencan/eval_callbacks.f
       src/gencan/easygencan.f
       src/gencan/gencan_driver.f
       src/gencan/gencan_bridge.f90
       src/gencan/spgls.f
       src/gencan/spgls_bridge.f90
       src/gencan/cg.f
       src/gencan/cg_bridge.f90
       src/gencan/tnls.f
       src/gencan/tnls_bridge.f90
       src/gencan/reduced_wrappers.f
       src/gencan/vector_ops.f
       src/gencan/numeric_utils.f
       src/gencan_ieee_signal_routines.f90
       src/pgencan.f90
       src/initial.f90
       src/title.f90
       src/setsizes.f90
       src/getinp.f90
       src/strlength.f90
       src/output.f90
       src/cli_parser.f90
       src/checkpoint.f90
       src/writesuccess.f90
       src/fparc.f90
       src/gparc.f90
       src/gwalls.f90
       src/comprest.f90
       src/comparegrad.f90
       src/polartocart.f90
       src/resetcells.f90
       src/tobar.f90
       src/cell_indexing.f90
       src/restmol.f90
       src/swaptype.f90
       src/swaptypemod.f90
       src/ahestetic.f90
       src/heuristics.f90
       src/flashsort.f90
       src/jacobi.f90
       src/random.f90
       src/sizes.f90
       src/pbc.f90
       src/usegencan.f90
       src/gencan_cpp_bridge.f90
       src/compute_data.f90
       src/flashmod.f90
       src/computef.f90
       src/computeg.f90
       src/eval_state_bridge.f90
       src/input.f90
)

add_library(packmol_gencan_cpp_obj OBJECT
    src/cpp/gencan/gencan_interface.cpp
    src/cpp/gencan/gencan_runtime.cpp
    src/cpp/gencan/gencan_eval_interface.cpp
    src/cpp/gencan/gencan_post.cpp
    src/cpp/gencan/gencan_eval.cpp
    src/cpp/gencan/gencan_eval_state.cpp
    src/cpp/gencan/gencan_eval_objective.cpp
    src/cpp/gencan/gencan_spgls.cpp
    src/cpp/gencan/gencan_tnls.cpp
    src/cpp/gencan/gencan_cg.cpp
    src/cpp/gencan/gencan_numeric.cpp
    src/cpp/gencan/gencan_probe.cpp
    src/cpp/gencan/gencan_easy.cpp
    src/cpp/gencan/gencan_solver.cpp
)
target_compile_features(packmol_gencan_cpp_obj PUBLIC cxx_std_17)
target_include_directories(packmol_gencan_cpp_obj PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp
)

add_library(packmol_core STATIC
    ${PACKMOL_CORE_SOURCES}
    $<TARGET_OBJECTS:packmol_gencan_cpp_obj>
)

target_include_directories(packmol_core PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Build the executable
add_executable(packmol
       app/packmol.f90
)

target_link_libraries(packmol PRIVATE packmol_core)

# Installation directive
install(TARGETS packmol DESTINATION bin)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
