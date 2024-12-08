cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# ---- Declare library ----

add_library(
    elecstruct_elecstruct
    atoms.cpp
    basis/basis_sets/sto3g.cpp
    geometry.cpp
    input_file_parser/input_file_parser.cpp
    input_file_parser/parsed_information.cpp
    integrals/boys.cpp
    integrals/electron_electron_index_iterator.cpp
    integrals/electron_electron_integrals.cpp
    integrals/f_coefficient.cpp
    integrals/kinetic_integrals.cpp
    integrals/nuclear_electron_index_iterator.cpp
    integrals/nuclear_electron_integrals.cpp
    integrals/overlap_integrals.cpp
    integrals/two_electron_integral_grid.cpp
    mathtools/gaussian.cpp
    mathtools/misc.cpp
    restricted_hartree_fock/initial_density_matrix.cpp
    restricted_hartree_fock/restricted_hartree_fock.cpp
    restricted_hartree_fock/step.cpp
    matrices.cpp
    orbitals.cpp
)
add_library(
    elecstruct::elecstruct
    ALIAS elecstruct_elecstruct
)

# ---- Modify Library Properties ----

include(GenerateExportHeader)
generate_export_header(
    elecstruct_elecstruct
    BASE_NAME elecstruct
    EXPORT_FILE_NAME export/elecstruct/elecstruct_export.hpp
    CUSTOM_CONTENT_FROM_VARIABLE pragma_suppress_c4251
)

if(NOT BUILD_SHARED_LIBS)
    target_compile_definitions(
        elecstruct_elecstruct
        PUBLIC ELECSTRUCT_STATIC_DEFINE
    )
endif()

set_target_properties(
    elecstruct_elecstruct PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
    VERSION "${PROJECT_VERSION}"
    SOVERSION "${PROJECT_VERSION_MAJOR}"
    EXPORT_NAME elecstruct
    OUTPUT_NAME elecstruct
)

target_include_directories(
    elecstruct_elecstruct
    ${warning_guard}
    PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern>"
)

target_include_directories(
    elecstruct_elecstruct
    SYSTEM
    PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/source/export>"
)

target_include_directories(
    elecstruct_elecstruct
    ${warning_guard}
    SYSTEM PRIVATE "${EIGEN_PATH}"
)

target_compile_features(elecstruct_elecstruct PUBLIC cxx_std_20)
