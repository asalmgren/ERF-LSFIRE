set(erf_exe_name erf_elmfire_constant)

add_executable(${erf_exe_name} "")

target_sources(${erf_exe_name}
  PRIVATE
    ERF_Prob.cpp
)

target_include_directories(${erf_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

include(${CMAKE_SOURCE_DIR}/CMake/BuildERFExe.cmake)
build_erf_exe(${erf_exe_name})

if(ERF_ENABLE_REGRESSION_TESTS OR ERF_ENABLE_TESTS)
  add_test(NAME erf_elmfire_constant COMMAND ${erf_exe_name} inputs_elmfire_constant max_step=10)
  set_tests_properties(erf_elmfire_constant PROPERTIES
                       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                       LABELS "unit"
                       TIMEOUT 5400)
endif()
