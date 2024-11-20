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
    set(EIGEN_OPTION_ARG ENABLE_EIGEN)
    set(ONE_VALUE_ARGS TARGET)
    set(MULTI_VALUE_ARGS SOURCES)
    cmake_parse_arguments(add_test_target "${EIGEN_OPTION_ARG}" "${ONE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}" "${ARGN}")

    add_executable(
        ${add_test_target_TARGET}
        ${add_test_target_SOURCES}
    )

    target_compile_features(
        ${add_test_target_TARGET}
        PRIVATE cxx_std_20
    )

    if(add_test_target_ENABLE_EIGEN)
        target_include_directories(
            ${add_test_target_TARGET}
            ${warning_guard}
            PRIVATE "${SOURCE_FILES_DIR}"
        )

        # separating the statements makes the warning guards work; otherwise we get lots of compilation warnings
        # from these external libraries
        target_include_directories(
            ${add_test_target_TARGET}
            ${warning_guard}
            SYSTEM PRIVATE "${PIMC_SIM_SOURCE_DIR}"
        )

        target_include_directories(
            ${add_test_target_TARGET}
            ${warning_guard}
            SYSTEM PRIVATE "${EIGEN_PATH}"
        )
    else()
        target_include_directories(
            ${add_test_target_TARGET}
            ${warning_guard}
            PRIVATE "${SOURCE_FILES_DIR}"
        )
    endif()

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
add_test_target(TARGET gaussian_test SOURCES "source/gaussian_test.cpp")
add_test_target(TARGET mathtools_test SOURCES "source/mathtools_test.cpp")
add_test_target(TARGET nuclear_electron_index_iterator_test SOURCES "source/nuclear_electron_index_iterator_test.cpp")
add_test_target(TARGET electron_electron_index_iterator_test SOURCES "source/electron_electron_index_iterator_test.cpp")
add_test_target(TARGET kinetic_integrals_test SOURCES "source/kinetic_integrals_test.cpp")
add_test_target(ENABLE_EIGEN TARGET overlap_matrix_test SOURCES "source/overlap_matrix_test.cpp")

# ---- End-of-file commands ----

add_folders(Test)
