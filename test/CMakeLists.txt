cmake_minimum_required(VERSION 3.14)

project(elecstructTests LANGUAGES CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)
include(../cmake/windows-set-path.cmake)

# ---- Dependencies ----

if(PROJECT_IS_TOP_LEVEL)
  find_package(elecstruct REQUIRED)
  enable_testing()
endif()

# ---- Tests ----

add_executable(elecstruct_test source/elecstruct_test.cpp)
target_link_libraries(elecstruct_test PRIVATE elecstruct::elecstruct)
target_compile_features(elecstruct_test PRIVATE cxx_std_20)

add_test(NAME elecstruct_test COMMAND elecstruct_test)
windows_set_path(elecstruct_test elecstruct::elecstruct)

# ---- End-of-file commands ----

add_folders(Test)
