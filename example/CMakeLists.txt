cmake_minimum_required(VERSION 3.14)

project(elecstructExamples CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)

if(PROJECT_IS_TOP_LEVEL)
  find_package(elecstruct REQUIRED)
endif()

add_custom_target(run-examples)

function(add_example NAME)
  add_executable("${NAME}" "${NAME}.cpp")
  target_link_libraries("${NAME}" PRIVATE elecstruct::elecstruct)
  target_compile_features("${NAME}" PRIVATE cxx_std_20)
  add_custom_target("run_${NAME}" COMMAND "${NAME}" VERBATIM)
  add_dependencies("run_${NAME}" "${NAME}")
  add_dependencies(run-examples "run_${NAME}")
endfunction()

add_example(water_sto3g_example)
target_include_directories(
    water_sto3g_example
    ${warning_guard}
    SYSTEM PRIVATE "${EIGEN_PATH}"
)

add_example(main)
target_include_directories(
    main
    ${warning_guard}
    SYSTEM PRIVATE "${EIGEN_PATH}"
)

add_folders(Example)
