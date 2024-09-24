if(PROJECT_IS_TOP_LEVEL)
  set(
      CMAKE_INSTALL_INCLUDEDIR "include/elecstruct-${PROJECT_VERSION}"
      CACHE PATH ""
  )
endif()

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# find_package(<package>) call for consumers to find this project
set(package elecstruct)

install(
    DIRECTORY
    include/
    "${PROJECT_BINARY_DIR}/export/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    COMPONENT elecstruct_Development
)

install(
    TARGETS elecstruct_elecstruct
    EXPORT elecstructTargets
    RUNTIME #
    COMPONENT elecstruct_Runtime
    LIBRARY #
    COMPONENT elecstruct_Runtime
    NAMELINK_COMPONENT elecstruct_Development
    ARCHIVE #
    COMPONENT elecstruct_Development
    INCLUDES #
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

write_basic_package_version_file(
    "${package}ConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion
)

# Allow package maintainers to freely override the path for the configs
set(
    elecstruct_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/${package}"
    CACHE PATH "CMake package config location relative to the install prefix"
)
mark_as_advanced(elecstruct_INSTALL_CMAKEDIR)

install(
    FILES cmake/install-config.cmake
    DESTINATION "${elecstruct_INSTALL_CMAKEDIR}"
    RENAME "${package}Config.cmake"
    COMPONENT elecstruct_Development
)

install(
    FILES "${PROJECT_BINARY_DIR}/${package}ConfigVersion.cmake"
    DESTINATION "${elecstruct_INSTALL_CMAKEDIR}"
    COMPONENT elecstruct_Development
)

install(
    EXPORT elecstructTargets
    NAMESPACE elecstruct::
    DESTINATION "${elecstruct_INSTALL_CMAKEDIR}"
    COMPONENT elecstruct_Development
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()
