set(erf_exe_name erf_elmfire_constant)

add_executable(${erf_exe_name} "")

target_sources(${erf_exe_name}
  PRIVATE
    ERF_Prob.cpp
    ${CMAKE_SOURCE_DIR}/Source/main.cpp
)

target_include_directories(${erf_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(${erf_exe_name} PRIVATE ${erf_lib_name})

install(TARGETS ${erf_exe_name}
        RUNTIME DESTINATION bin)

if(ERF_ENABLE_REGRESSION_TESTS OR ERF_ENABLE_TESTS)
  add_test(NAME erf_elmfire_constant COMMAND ${erf_exe_name} inputs_elmfire_constant max_step=10)
  set_tests_properties(erf_elmfire_constant PROPERTIES
                       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                       LABELS "unit"
                       TIMEOUT 5400)
endif()
