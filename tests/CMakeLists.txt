find_package(Python3 COMPONENTS Interpreter REQUIRED)

add_library(packmol_test_support STATIC
    support/test_assertions.f90
    support/test_runtime_setup.f90
)

target_link_libraries(packmol_test_support PUBLIC packmol_core)

function(packmol_add_test name source label)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE packmol_test_support packmol_core)
    add_test(NAME ${name} COMMAND $<TARGET_FILE:${name}>)
    set_tests_properties(${name} PROPERTIES LABELS "${label}")
endfunction()

packmol_add_test(test_strings unit/test_strings.f90 unit)
packmol_add_test(test_geometry unit/test_geometry.f90 unit)
packmol_add_test(test_cell_indexing unit/test_cell_indexing.f90 unit)
packmol_add_test(test_random unit/test_random.f90 unit)
packmol_add_test(test_flashsort unit/test_flashsort.f90 unit)
packmol_add_test(test_output_helpers unit/test_output_helpers.f90 unit)
packmol_add_test(test_pbc unit/test_pbc.f90 unit)
packmol_add_test(test_ieee_signal unit/test_ieee_signal.f90 unit)
packmol_add_test(test_jacobi unit/test_jacobi.f90 unit)
packmol_add_test(test_gencan_vector_ops gencan/test_gencan_vector_ops.f90 gencan)
packmol_add_test(test_gencan_wrappers gencan/test_gencan_wrappers.f90 gencan)
packmol_add_test(test_molecule_centers component/test_molecule_centers.f90 component)
packmol_add_test(test_input_helpers component/test_input_helpers.f90 component)
packmol_add_test(test_getkeywords component/test_getkeywords.f90 component)
add_executable(test_getinp component/test_getinp.f90)
target_link_libraries(test_getinp PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_getinp
    COMMAND $<TARGET_FILE:test_getinp> ${CMAKE_SOURCE_DIR}/tests/fixtures/minimal_packmol.inp
)
set_tests_properties(test_getinp PROPERTIES LABELS "component")
set_tests_properties(test_getinp PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/fixtures")
packmol_add_test(test_discret component/test_discret.f90 component)
packmol_add_test(test_resetcells component/test_resetcells.f90 component)
packmol_add_test(test_packmolprecision component/test_packmolprecision.f90 component)
packmol_add_test(test_core_evaluators component/test_core_evaluators.f90 component)
packmol_add_test(test_gencan_cpp_bridge component/test_gencan_cpp_bridge.f90 component)
packmol_add_test(test_constraint_kernels component/test_constraint_kernels.f90 component)
packmol_add_test(test_restmol component/test_restmol.f90 component)
packmol_add_test(test_swaptype component/test_swaptype.f90 component)
packmol_add_test(test_pgencan component/test_pgencan.f90 component)
add_executable(test_initial component/test_initial.f90)
target_link_libraries(test_initial PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_initial
    COMMAND $<TARGET_FILE:test_initial> ${CMAKE_SOURCE_DIR}/tests/fixtures/tiny_packing.inp
)
set_tests_properties(test_initial PROPERTIES LABELS "component")
set_tests_properties(test_initial PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/fixtures")

add_executable(test_movebad component/test_movebad.f90)
target_link_libraries(test_movebad PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_movebad
    COMMAND $<TARGET_FILE:test_movebad> ${CMAKE_SOURCE_DIR}/tests/fixtures/tiny_packing.inp
)
set_tests_properties(test_movebad PROPERTIES LABELS "component")
set_tests_properties(test_movebad PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/fixtures")
packmol_add_test(test_write_connect component/test_write_connect.f90 component)

add_executable(test_cli_parser_probe component/test_cli_parser_probe.f90)
target_link_libraries(test_cli_parser_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_cli_parser_with_input
    COMMAND test_cli_parser_probe -i fixtures/input.inp
)
set_tests_properties(test_cli_parser_with_input PROPERTIES LABELS "component")

add_test(
    NAME test_cli_parser_with_input_output
    COMMAND test_cli_parser_probe -i fixtures/input.inp -o fixtures/output.pdb
)
set_tests_properties(test_cli_parser_with_input_output PROPERTIES LABELS "component")

add_executable(test_cli_parser_functions_probe component/test_cli_parser_functions_probe.f90)
target_link_libraries(test_cli_parser_functions_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_cli_parser_functions
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/component/check_cli_parser_functions.py $<TARGET_FILE:test_cli_parser_functions_probe>
)
set_tests_properties(test_cli_parser_functions PROPERTIES LABELS "component")

add_test(
    NAME inventory_check
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/generate_inventory.py --root ${CMAKE_SOURCE_DIR} --check
)
set_tests_properties(inventory_check PROPERTIES LABELS "component")

add_executable(test_failopen_probe component/test_failopen_probe.f90)
target_link_libraries(test_failopen_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_failopen
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/component/check_failopen.py $<TARGET_FILE:test_failopen_probe>
)
set_tests_properties(test_failopen PROPERTIES LABELS "component")

add_executable(test_output_messages_probe component/test_output_messages_probe.f90)
target_link_libraries(test_output_messages_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_output_messages
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/component/check_output_messages.py $<TARGET_FILE:test_output_messages_probe>
)
set_tests_properties(test_output_messages PROPERTIES LABELS "component")

add_executable(test_comparegrad_probe component/test_comparegrad_probe.f90)
target_link_libraries(test_comparegrad_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_comparegrad
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/component/check_comparegrad.py $<TARGET_FILE:test_comparegrad_probe>
)
set_tests_properties(test_comparegrad PROPERTIES LABELS "component")

add_executable(test_output_probe component/test_output_probe.f90)
target_link_libraries(test_output_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_output
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/component/check_output_files.py output $<TARGET_FILE:test_output_probe> ${CMAKE_SOURCE_DIR}/tests/fixtures/minimal_packmol.inp
)
set_tests_properties(test_output PROPERTIES LABELS "component")

add_executable(test_checkpoint_probe component/test_checkpoint_probe.f90)
target_link_libraries(test_checkpoint_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_checkpoint
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/component/check_output_files.py checkpoint $<TARGET_FILE:test_checkpoint_probe> ${CMAKE_SOURCE_DIR}/tests/fixtures/minimal_packmol.inp
)
set_tests_properties(test_checkpoint PROPERTIES LABELS "component")

add_executable(test_setsizes component/test_setsizes.f90)
target_link_libraries(test_setsizes PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_setsizes
    COMMAND $<TARGET_FILE:test_setsizes> ${CMAKE_SOURCE_DIR}/tests/fixtures/minimal_packmol.inp
)
set_tests_properties(test_setsizes PROPERTIES LABELS "component")
set_tests_properties(test_setsizes PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/fixtures")

packmol_add_test(test_gencan_hessian_diff gencan/test_gencan_hessian_diff.f90 gencan)
packmol_add_test(test_gencan_eval_wrappers gencan/test_gencan_eval_wrappers.f90 gencan)
packmol_add_test(test_gencan_calchd gencan/test_gencan_calchd.f90 gencan)
set_tests_properties(
    test_gencan_wrappers
    test_gencan_hessian_diff
    test_gencan_eval_wrappers
    test_gencan_calchd
    PROPERTIES LABELS "gencan;gencan_wrappers"
)
packmol_add_test(test_spgls gencan/test_spgls.f90 gencan)
packmol_add_test(test_spgls_limits gencan/test_spgls_limits.f90 gencan)
add_executable(test_spgls_ab_probe gencan/test_spgls_ab_probe.f90)
target_link_libraries(test_spgls_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_spgls_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_spgls_ab.py
        $<TARGET_FILE:test_spgls_ab_probe>
)
set_tests_properties(test_spgls_ab PROPERTIES LABELS "gencan;gencan_ab")
packmol_add_test(test_cg gencan/test_cg.f90 gencan)
packmol_add_test(test_cg_limits gencan/test_cg_limits.f90 gencan)
add_executable(test_cg_cpp_evalhd_callback gencan/test_cg_cpp_evalhd_callback.f90)
target_link_libraries(test_cg_cpp_evalhd_callback PRIVATE packmol_test_support packmol_core)
add_test(NAME test_cg_cpp_evalhd_callback COMMAND $<TARGET_FILE:test_cg_cpp_evalhd_callback>)
set_tests_properties(
    test_cg_cpp_evalhd_callback
    PROPERTIES
        LABELS "gencan"
        ENVIRONMENT "PACKMOL_GENCAN_IMPL=cpp"
)
add_executable(test_cg_ab_probe gencan/test_cg_ab_probe.f90)
target_link_libraries(test_cg_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_cg_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_cg_ab.py
        $<TARGET_FILE:test_cg_ab_probe>
)
set_tests_properties(test_cg_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_cg_nqmp_ab_probe gencan/test_cg_nqmp_ab_probe.f90)
target_link_libraries(test_cg_nqmp_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_cg_nqmp_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_cg_ab.py
        $<TARGET_FILE:test_cg_nqmp_ab_probe>
)
set_tests_properties(test_cg_nqmp_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_cg_general_ab_probe gencan/test_cg_general_ab_probe.f90)
target_link_libraries(test_cg_general_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_cg_general_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_cg_ab.py
        $<TARGET_FILE:test_cg_general_ab_probe>
)
set_tests_properties(test_cg_general_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_cg_samep_ab_probe gencan/test_cg_samep_ab_probe.f90)
target_link_libraries(test_cg_samep_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_cg_samep_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_cg_ab.py
        $<TARGET_FILE:test_cg_samep_ab_probe>
)
set_tests_properties(test_cg_samep_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_cg_nqmpstop_ab_probe gencan/test_cg_nqmpstop_ab_probe.f90)
target_link_libraries(test_cg_nqmpstop_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_cg_nqmpstop_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_cg_ab.py
        $<TARGET_FILE:test_cg_nqmpstop_ab_probe>
)
set_tests_properties(test_cg_nqmpstop_ab PROPERTIES LABELS "gencan;gencan_ab")
packmol_add_test(test_tnls gencan/test_tnls.f90 gencan)
packmol_add_test(test_tnls_limits gencan/test_tnls_limits.f90 gencan)
add_executable(test_tnls_ab_probe gencan/test_tnls_ab_probe.f90)
target_link_libraries(test_tnls_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_ab.py
        $<TARGET_FILE:test_tnls_ab_probe>
)
set_tests_properties(test_tnls_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_tnls_nobound_ab_probe gencan/test_tnls_nobound_ab_probe.f90)
target_link_libraries(test_tnls_nobound_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_nobound_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_ab.py
        $<TARGET_FILE:test_tnls_nobound_ab_probe>
)
set_tests_properties(test_tnls_nobound_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_tnls_nobound_interior_ab_probe gencan/test_tnls_nobound_interior_ab_probe.f90)
target_link_libraries(test_tnls_nobound_interior_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_nobound_interior_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_ab.py
        $<TARGET_FILE:test_tnls_nobound_interior_ab_probe>
)
set_tests_properties(test_tnls_nobound_interior_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_tnls_nobound_interior_noextrap_ab_probe gencan/test_tnls_nobound_interior_noextrap_ab_probe.f90)
target_link_libraries(test_tnls_nobound_interior_noextrap_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_nobound_interior_noextrap_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_ab.py
        $<TARGET_FILE:test_tnls_nobound_interior_noextrap_ab_probe>
)
set_tests_properties(test_tnls_nobound_interior_noextrap_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_tnls_nobound_interior_ex1_ab_probe gencan/test_tnls_nobound_interior_ex1_ab_probe.f90)
target_link_libraries(test_tnls_nobound_interior_ex1_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_nobound_interior_ex1_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_ab.py
        $<TARGET_FILE:test_tnls_nobound_interior_ex1_ab_probe>
)
set_tests_properties(test_tnls_nobound_interior_ex1_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_tnls_nobound_interior_ex2_ab_probe gencan/test_tnls_nobound_interior_ex2_ab_probe.f90)
target_link_libraries(test_tnls_nobound_interior_ex2_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_nobound_interior_ex2_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_ab.py
        $<TARGET_FILE:test_tnls_nobound_interior_ex2_ab_probe>
)
set_tests_properties(test_tnls_nobound_interior_ex2_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_tnls_nobound_interior_ex3_ab_probe gencan/test_tnls_nobound_interior_ex3_ab_probe.f90)
target_link_libraries(test_tnls_nobound_interior_ex3_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_nobound_interior_ex3_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_ab.py
        $<TARGET_FILE:test_tnls_nobound_interior_ex3_ab_probe>
)
set_tests_properties(test_tnls_nobound_interior_ex3_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_tnls_nobound_interior_ex4_ab_probe gencan/test_tnls_nobound_interior_ex4_ab_probe.f90)
target_link_libraries(test_tnls_nobound_interior_ex4_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_nobound_interior_ex4_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_ab.py
        $<TARGET_FILE:test_tnls_nobound_interior_ex4_ab_probe>
)
set_tests_properties(test_tnls_nobound_interior_ex4_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_tnls_nobound_extrap_ab_probe gencan/test_tnls_nobound_extrap_ab_probe.f90)
target_link_libraries(test_tnls_nobound_extrap_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_nobound_extrap_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_ab.py
        $<TARGET_FILE:test_tnls_nobound_extrap_ab_probe>
)
set_tests_properties(test_tnls_nobound_extrap_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_tnls_nobound_mininterp2_ab_probe gencan/test_tnls_nobound_mininterp2_ab_probe.f90)
target_link_libraries(test_tnls_nobound_mininterp2_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_nobound_mininterp2_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_ab.py
        $<TARGET_FILE:test_tnls_nobound_mininterp2_ab_probe>
)
set_tests_properties(test_tnls_nobound_mininterp2_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_tnls_extrap_ab_probe gencan/test_tnls_extrap_ab_probe.f90)
target_link_libraries(test_tnls_extrap_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_extrap_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_extrap_ab.py
        $<TARGET_FILE:test_tnls_extrap_ab_probe>
)
set_tests_properties(test_tnls_extrap_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_tnls_interior_ab_probe gencan/test_tnls_interior_ab_probe.f90)
target_link_libraries(test_tnls_interior_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_interior_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_ab.py
        $<TARGET_FILE:test_tnls_interior_ab_probe>
)
set_tests_properties(test_tnls_interior_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_tnls_interior_interp_ab_probe gencan/test_tnls_interior_interp_ab_probe.f90)
target_link_libraries(test_tnls_interior_interp_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_interior_interp_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_ab.py
        $<TARGET_FILE:test_tnls_interior_interp_ab_probe>
)
set_tests_properties(test_tnls_interior_interp_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_tnls_interior_beta_extrap_ab_probe gencan/test_tnls_interior_beta_extrap_ab_probe.f90)
target_link_libraries(test_tnls_interior_beta_extrap_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_interior_beta_extrap_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_ab.py
        $<TARGET_FILE:test_tnls_interior_beta_extrap_ab_probe>
)
set_tests_properties(test_tnls_interior_beta_extrap_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_tnls_gtype2_ab_probe gencan/test_tnls_gtype2_ab_probe.f90)
target_link_libraries(test_tnls_gtype2_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_tnls_gtype2_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tnls_ab.py
        $<TARGET_FILE:test_tnls_gtype2_ab_probe>
)
set_tests_properties(test_tnls_gtype2_ab PROPERTIES LABELS "gencan;gencan_ab")
add_executable(test_easygencan_ab_probe gencan/test_easygencan_ab_probe.f90)
target_link_libraries(test_easygencan_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_easygencan_ab
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_easygencan_ab.py
        $<TARGET_FILE:test_easygencan_ab_probe>
)
set_tests_properties(test_easygencan_ab PROPERTIES LABELS "gencan;gencan_ab")
add_test(
    NAME test_easygencan_ab_draft
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_easygencan_ab.py
        $<TARGET_FILE:test_easygencan_ab_probe> ${CMAKE_SOURCE_DIR}/tests/fixtures/tiny_packing.inp
)
set_tests_properties(
    test_easygencan_ab_draft
    PROPERTIES
        LABELS "gencan;gencan_ab;gencan_draft"
)
add_executable(test_gencan_runtime_driver_probe gencan/test_gencan_runtime_driver_probe.f90)
target_link_libraries(test_gencan_runtime_driver_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_gencan_runtime_driver_probe
    COMMAND $<TARGET_FILE:test_gencan_runtime_driver_probe> ${CMAKE_SOURCE_DIR}/tests/fixtures/tiny_packing.inp
)
set_tests_properties(test_gencan_runtime_driver_probe PROPERTIES LABELS "gencan")
set_tests_properties(test_gencan_runtime_driver_probe PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/fixtures")
add_executable(test_gencan_ab_probe gencan/test_gencan_ab_probe.f90)
target_link_libraries(test_gencan_ab_probe PRIVATE packmol_test_support packmol_core)
add_executable(test_gencan_spg_post_ab_probe gencan/test_gencan_spg_post_ab_probe.f90)
target_link_libraries(test_gencan_spg_post_ab_probe PRIVATE packmol_test_support packmol_core)
add_executable(test_gencan_entry_stop_ab_probe gencan/test_gencan_entry_stop_ab_probe.f90)
target_link_libraries(test_gencan_entry_stop_ab_probe PRIVATE packmol_test_support packmol_core)
add_executable(test_gencan_entry_stop_gtype2_ab_probe gencan/test_gencan_entry_stop_gtype2_ab_probe.f90)
target_link_libraries(test_gencan_entry_stop_gtype2_ab_probe PRIVATE packmol_test_support packmol_core)
add_executable(test_gencan_entry_stop_inform2_ab_probe gencan/test_gencan_entry_stop_inform2_ab_probe.f90)
target_link_libraries(test_gencan_entry_stop_inform2_ab_probe PRIVATE packmol_test_support packmol_core)
add_executable(test_gencan_entry_stop_inform3_ab_probe gencan/test_gencan_entry_stop_inform3_ab_probe.f90)
target_link_libraries(test_gencan_entry_stop_inform3_ab_probe PRIVATE packmol_test_support packmol_core)
add_test(
    NAME test_gencan_tn_draft_smoke
    COMMAND $<TARGET_FILE:test_gencan_runtime_driver_probe> ${CMAKE_SOURCE_DIR}/tests/fixtures/tiny_packing.inp
)
set_tests_properties(
    test_gencan_tn_draft_smoke
    PROPERTIES
        LABELS "gencan;gencan_draft"
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/fixtures"
        ENVIRONMENT "PACKMOL_GENCAN_IMPL=cpp"
)
add_test(
    NAME test_gencan_spg_post_no_fortran_tail_cpp_closure
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_spg_post_no_fortran_tail_cpp_closure.py
)
set_tests_properties(
    test_gencan_spg_post_no_fortran_tail_cpp_closure
    PROPERTIES
        LABELS "gencan;gencan_closure"
)
add_test(
    NAME test_gencan_tn_post_no_fortran_tail_cpp_closure
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_tn_post_no_fortran_tail_cpp_closure.py
)
set_tests_properties(
    test_gencan_tn_post_no_fortran_tail_cpp_closure
    PROPERTIES
        LABELS "gencan;gencan_closure"
)
add_test(
    NAME test_gencan_no_fallback_matrix_closure
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_gencan_no_fallback_matrix.py
)
set_tests_properties(
    test_gencan_no_fallback_matrix_closure
    PROPERTIES
        LABELS "gencan;gencan_closure"
)
add_test(
    NAME test_gencan_no_fallback_matrix_ab_drift
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/gencan/check_gencan_no_fallback_matrix.py
)
set_tests_properties(
    test_gencan_no_fallback_matrix_ab_drift
    PROPERTIES
        LABELS "gencan;gencan_ab"
        ENVIRONMENT "PACKMOL_GENCAN_AB_COMPARE=1"
)
packmol_add_test(test_easygencan gencan/test_easygencan.f90 gencan)
packmol_add_test(test_gencan_driver gencan/test_gencan_driver.f90 gencan)

add_test(
    NAME regression_cli
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR}/testing ./test_cli.sh $<TARGET_FILE:packmol>
)
set_tests_properties(regression_cli PROPERTIES LABELS "regression")

add_test(
    NAME regression_failed_run
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR}/testing ./test_failed.sh ./input_files/water_box_failed.inp packmol.log FORCED
)
set_tests_properties(regression_failed_run PROPERTIES LABELS "regression")
set_tests_properties(regression_failed_run PROPERTIES ENVIRONMENT "PACKMOL_BIN=$<TARGET_FILE:packmol>")

add_test(
    NAME regression_connectivity
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR}/testing ./test_connectivity.sh
)
set_tests_properties(regression_connectivity PROPERTIES LABELS "regression")
set_tests_properties(regression_connectivity PROPERTIES ENVIRONMENT "PACKMOL_BIN=$<TARGET_FILE:packmol>")
