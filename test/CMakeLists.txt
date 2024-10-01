cmake_minimum_required(VERSION 3.14)

project(elecstructTests LANGUAGES CXX)

# ---- Get access to Catch2 ----
 
include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY git@github.com:catchorg/Catch2.git
    GIT_TAG        v3.4.0
)
FetchContent_MakeAvailable(Catch2)
include(Catch)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)
include(../cmake/windows-set-path.cmake)

# ---- Dependencies ----

if(PROJECT_IS_TOP_LEVEL)
    find_package(elecstruct REQUIRED)
    enable_testing()
endif()

# ---- Get access to multiple argument parser ----
include(CMakeParseArguments)

# ---- Tests ----

cmake_path(GET PROJECT_SOURCE_DIR PARENT_PATH PIMC_SIM_SOURCE_DIR)
set(SOURCE_FILES_DIR "${PIMC_SIM_SOURCE_DIR}/source")

function(add_test_target)
    set(ONE_VALUE_ARGS TARGET)
    set(MULTI_VALUE_ARGS SOURCES)
    cmake_parse_arguments(add_test_target "" "${ONE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}" "${ARGN}")

    add_executable(
        ${add_test_target_TARGET}
        ${add_test_target_SOURCES}
    )

    target_compile_features(
        ${add_test_target_TARGET}
        PRIVATE cxx_std_20
    )

    target_include_directories(
        ${add_test_target_TARGET}
        ${warning_guard}
        PRIVATE "${SOURCE_FILES_DIR}"
    )

    target_link_libraries(
        ${add_test_target_TARGET}
        PRIVATE Catch2::Catch2WithMain
        PRIVATE elecstruct::elecstruct
    )

    catch_discover_tests(${add_test_target_TARGET})
endfunction()

# ---- Tests ----

add_test_target(TARGET cartesian3d_test SOURCES "source/cartesian3d_test.cpp")
add_test_target(TARGET grid2d_test SOURCES "source/grid2d_test.cpp")
add_test_target(TARGET grid3d_test SOURCES "source/grid3d_test.cpp")
add_test_target(TARGET grid4d_test SOURCES "source/grid4d_test.cpp")
add_test_target(TARGET geometry_test SOURCES "source/geometry_test.cpp")

# ---- End-of-file commands ----

add_folders(Test)
