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

function(packmol_add_gencan_compare_test name source)
    cmake_parse_arguments(ARG "" "FIXTURE" "INT_KEYS;FLOAT_KEYS;LABELS" ${ARGN})

    add_executable(${name}_probe ${source})
    target_link_libraries(${name}_probe PRIVATE packmol_test_support packmol_core)

    set(command
        ${Python3_EXECUTABLE}
        ${CMAKE_SOURCE_DIR}/tests/gencan/check_gencan_compare.py
        $<TARGET_FILE:${name}_probe>
    )

    if(ARG_FIXTURE)
        list(APPEND command --fixture ${ARG_FIXTURE})
    endif()

    foreach(key IN LISTS ARG_INT_KEYS)
        list(APPEND command --int ${key})
    endforeach()

    foreach(spec IN LISTS ARG_FLOAT_KEYS)
        list(APPEND command --float ${spec})
    endforeach()

    add_test(NAME ${name} COMMAND ${command})

    if(ARG_LABELS)
        set(labels "${ARG_LABELS}")
    else()
        set(labels "gencan;gencan_compare")
    endif()
    set_tests_properties(${name} PROPERTIES LABELS "${labels}")
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

add_executable(test_eval_objective_compare_probe component/test_eval_objective_compare_probe.f90)
target_link_libraries(test_eval_objective_compare_probe PRIVATE packmol_test_support packmol_core)

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
packmol_add_gencan_compare_test(
    test_spgls_ab
    gencan/test_spgls_ab_probe.f90
    INT_KEYS "inform;fcnt"
    FLOAT_KEYS "f;xsum;dnorm2"
    LABELS "gencan;gencan_compare"
)
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
packmol_add_gencan_compare_test(
    test_cg_ab
    gencan/test_cg_ab_probe.f90
    INT_KEYS "inform;iter"
    FLOAT_KEYS "snorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_cg_nqmp_ab
    gencan/test_cg_nqmp_ab_probe.f90
    INT_KEYS "inform;iter"
    FLOAT_KEYS "snorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_cg_general_ab
    gencan/test_cg_general_ab_probe.f90
    INT_KEYS "inform;iter"
    FLOAT_KEYS "snorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_cg_samep_ab
    gencan/test_cg_samep_ab_probe.f90
    INT_KEYS "inform;iter"
    FLOAT_KEYS "snorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_cg_nqmpstop_ab
    gencan/test_cg_nqmpstop_ab_probe.f90
    INT_KEYS "inform;iter"
    FLOAT_KEYS "snorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_test(test_tnls gencan/test_tnls.f90 gencan)
packmol_add_test(test_tnls_limits gencan/test_tnls_limits.f90 gencan)
packmol_add_gencan_compare_test(
    test_tnls_ab
    gencan/test_tnls_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_tnls_nobound_ab
    gencan/test_tnls_nobound_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_tnls_nobound_interior_ab
    gencan/test_tnls_nobound_interior_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_tnls_nobound_interior_noextrap_ab
    gencan/test_tnls_nobound_interior_noextrap_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_tnls_nobound_interior_ex1_ab
    gencan/test_tnls_nobound_interior_ex1_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_tnls_nobound_interior_ex2_ab
    gencan/test_tnls_nobound_interior_ex2_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_tnls_nobound_interior_ex3_ab
    gencan/test_tnls_nobound_interior_ex3_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_tnls_nobound_interior_ex4_ab
    gencan/test_tnls_nobound_interior_ex4_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_tnls_nobound_extrap_ab
    gencan/test_tnls_nobound_extrap_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_tnls_nobound_mininterp2_ab
    gencan/test_tnls_nobound_mininterp2_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_tnls_extrap_ab
    gencan/test_tnls_extrap_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_tnls_interior_ab
    gencan/test_tnls_interior_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_tnls_interior_interp_ab
    gencan/test_tnls_interior_interp_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_tnls_interior_beta_extrap_ab
    gencan/test_tnls_interior_beta_extrap_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_tnls_gtype2_ab
    gencan/test_tnls_gtype2_ab_probe.f90
    INT_KEYS "inform;fcnt;gcnt;intcnt;exgcnt;exbcnt"
    FLOAT_KEYS "f;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_easygencan_ab
    gencan/test_easygencan_ab_probe.f90
    INT_KEYS "inform;iter;fcnt;gcnt"
    FLOAT_KEYS "f;gpsupn;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
)
packmol_add_gencan_compare_test(
    test_gencan_ab
    gencan/test_gencan_compare_probe.f90
    INT_KEYS "inform;iter;fcnt;gcnt;cgcnt"
    FLOAT_KEYS "f;gpsupn;xsum;gnorm2"
    LABELS "gencan;gencan_compare"
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
