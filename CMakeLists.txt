cmake_minimum_required(VERSION 3.21)

# ---------------------------------------------------------------------------
# Version
# Priority: cmake/version.cmake > ENV{LEMON_VERSION} > FATAL_ERROR
# ---------------------------------------------------------------------------
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.cmake")
  include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.cmake")
elseif(DEFINED ENV{LEMON_VERSION})
  set(LEMON_VERSION "$ENV{LEMON_VERSION}" CACHE STRING "LEMON version string.")
else()
  message(FATAL_ERROR
    "LEMON version could not be determined.\n"
    "Please set the LEMON_VERSION environment variable before configuring.\n"
    "Example: LEMON_VERSION=1.3.1 cmake ..."
  )
endif()

# ---------------------------------------------------------------------------
# Project
# ---------------------------------------------------------------------------
project(LEMON
  VERSION   "${LEMON_VERSION}"
  LANGUAGES C CXX
)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# ---------------------------------------------------------------------------
# C++ standard
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD          20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

# ---------------------------------------------------------------------------
# Optional dependencies (documentation)
# ---------------------------------------------------------------------------
find_package(Doxygen OPTIONAL_COMPONENTS dot)
find_package(Ghostscript QUIET)

# ---------------------------------------------------------------------------
# Optional solver back-ends
# ---------------------------------------------------------------------------
option(LEMON_ENABLE_GLPK   "Enable GLPK solver backend."         ON)
option(LEMON_ENABLE_ILOG   "Enable ILOG (CPLEX) solver backend." ON)
option(LEMON_ENABLE_COIN   "Enable COIN solver backend."         ON)
option(LEMON_ENABLE_SOPLEX "Enable SoPlex solver backend."       ON)

if(LEMON_ENABLE_GLPK)
  find_package(GLPK 4.33 QUIET)
endif()
if(LEMON_ENABLE_ILOG)
  find_package(ILOG QUIET)
endif()
if(LEMON_ENABLE_COIN)
  find_package(COIN QUIET)
endif()
if(LEMON_ENABLE_SOPLEX)
  find_package(SOPLEX QUIET)
endif()

# Capability flags derived from found solvers
if(GLPK_FOUND)
  set(LEMON_HAVE_LP   TRUE)
  set(LEMON_HAVE_MIP  TRUE)
  set(LEMON_HAVE_GLPK TRUE)
endif()
if(ILOG_FOUND)
  set(LEMON_HAVE_LP    TRUE)
  set(LEMON_HAVE_MIP   TRUE)
  set(LEMON_HAVE_CPLEX TRUE)
endif()
if(COIN_FOUND)
  set(LEMON_HAVE_LP  TRUE)
  set(LEMON_HAVE_MIP TRUE)
  set(LEMON_HAVE_CLP TRUE)
  set(LEMON_HAVE_CBC TRUE)
endif()
if(SOPLEX_FOUND)
  set(LEMON_HAVE_LP     TRUE)
  set(LEMON_HAVE_SOPLEX TRUE)
endif()

# Default LP / MIP solver selection (priority: CPLEX > CLP/CBC > GLPK > SOPLEX)
if(ILOG_FOUND)
  set(_default_lp  "CPLEX")
  set(_default_mip "CPLEX")
elseif(COIN_FOUND)
  set(_default_lp  "CLP")
  set(_default_mip "CBC")
elseif(GLPK_FOUND)
  set(_default_lp  "GLPK")
  set(_default_mip "GLPK")
elseif(SOPLEX_FOUND)
  set(_default_lp  "SOPLEX")
  set(_default_mip "")
else()
  set(_default_lp  "")
  set(_default_mip "")
endif()

# Validate / reset cached LP choice if the corresponding library is not found
set(_lp_needs_reset FALSE)
if(NOT LEMON_DEFAULT_LP)
  set(_lp_needs_reset TRUE)
elseif((LEMON_DEFAULT_LP STREQUAL "CPLEX"  AND NOT ILOG_FOUND)   OR
       (LEMON_DEFAULT_LP STREQUAL "CLP"    AND NOT COIN_FOUND)   OR
       (LEMON_DEFAULT_LP STREQUAL "GLPK"   AND NOT GLPK_FOUND)   OR
       (LEMON_DEFAULT_LP STREQUAL "SOPLEX" AND NOT SOPLEX_FOUND))
  set(_lp_needs_reset TRUE)
endif()
set(LEMON_DEFAULT_LP "${_default_lp}" CACHE STRING
  "Default LP solver backend (GLPK, CPLEX, CLP or SOPLEX)"
  FORCE)

set(_mip_needs_reset FALSE)
if(NOT LEMON_DEFAULT_MIP)
  set(_mip_needs_reset TRUE)
elseif((LEMON_DEFAULT_MIP STREQUAL "CPLEX" AND NOT ILOG_FOUND) OR
       (LEMON_DEFAULT_MIP STREQUAL "CBC"   AND NOT COIN_FOUND) OR
       (LEMON_DEFAULT_MIP STREQUAL "GLPK"  AND NOT GLPK_FOUND))
  set(_mip_needs_reset TRUE)
endif()
set(LEMON_DEFAULT_MIP "${_default_mip}" CACHE STRING
  "Default MIP solver backend (GLPK, CPLEX or CBC)"
  FORCE)

unset(_default_lp)
unset(_default_mip)
unset(_lp_needs_reset)
unset(_mip_needs_reset)

# ---------------------------------------------------------------------------
# Compiler warning flags
# ---------------------------------------------------------------------------
if(DEFINED ENV{LEMON_CXX_WARNING})
  # La variable d'env est une string — separate_arguments la découpe en liste CMake
  separate_arguments(_cxx_warning_list UNIX_COMMAND "$ENV{LEMON_CXX_WARNING}")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Chaque flag est un élément de liste séparé : pas de string(JOIN), pas de quotes
  set(_cxx_warning_list
    -Wall
    -Wextra
    -Wunused
    -Wformat=2
    -Wctor-dtor-privacy
    -Wnon-virtual-dtor
    -Wno-char-subscripts
    -Wwrite-strings
    -Wreturn-type
    -Wcast-qual
    -Wcast-align
    -Wsign-promo
    -Woverloaded-virtual
    -fno-strict-aliasing
    -Wold-style-cast
    -Wno-unknown-pragmas
  )
elseif(MSVC)
  # C4250: inherits via dominance  C4355: 'this' in initializer list
  # C4503: name truncated           C4800: bool coercion perf warning
  # C4996: deprecated function
  set(_cxx_warning_list /wd4250 /wd4355 /wd4503 /wd4800 /wd4996)
else()
  set(_cxx_warning_list -Wall)
endif()

# Cache : string lisible dans cmake-gui / ccmake
string(JOIN " " LEMON_CXX_WARNING_FLAGS ${_cxx_warning_list})
set(LEMON_CXX_WARNING_FLAGS "${LEMON_CXX_WARNING_FLAGS}"
  CACHE STRING "LEMON warning flags.")

# Interface target — transporte la liste de flags (éléments séparés, pas une string)
add_library(lemon_warnings INTERFACE)
target_compile_options(lemon_warnings INTERFACE
  $<BUILD_INTERFACE:${_cxx_warning_list}>
)
unset(_cxx_warning_list)

# ---------------------------------------------------------------------------
# Build type & Maintainer configuration
# ---------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Build type: None Debug Release RelWithDebInfo MinSizeRel Maintainer." FORCE)
endif()

if(MSVC)
  set(_maintainer_cxx  "/WX ${CMAKE_CXX_FLAGS_DEBUG}")
  set(_maintainer_c    "/WX ${CMAKE_C_FLAGS_DEBUG}")
else()
  set(_maintainer_cxx  "-Werror -ggdb -O0")
  set(_maintainer_c    "-Werror -O0")
endif()

set(CMAKE_CXX_FLAGS_MAINTAINER    "${_maintainer_cxx}"
  CACHE STRING "C++ flags for maintainer builds." FORCE)
set(CMAKE_C_FLAGS_MAINTAINER      "${_maintainer_c}"
  CACHE STRING "C flags for maintainer builds."   FORCE)
set(CMAKE_EXE_LINKER_FLAGS_MAINTAINER
  "${CMAKE_EXE_LINKER_FLAGS_DEBUG}"
  CACHE STRING "Linker flags for maintainer builds." FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_MAINTAINER
  "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}"
  CACHE STRING "Shared linker flags for maintainer builds." FORCE)

mark_as_advanced(
  CMAKE_CXX_FLAGS_MAINTAINER
  CMAKE_C_FLAGS_MAINTAINER
  CMAKE_EXE_LINKER_FLAGS_MAINTAINER
  CMAKE_SHARED_LINKER_FLAGS_MAINTAINER
)

if(CMAKE_CONFIGURATION_TYPES)
  list(APPEND CMAKE_CONFIGURATION_TYPES Maintainer)
  list(REMOVE_DUPLICATES CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}"
    CACHE STRING "Available build configurations." FORCE)
endif()

unset(_maintainer_cxx)
unset(_maintainer_c)

# ---------------------------------------------------------------------------
# Platform / type introspection
# ---------------------------------------------------------------------------
include(CheckTypeSize)
check_type_size("long long" LONG_LONG)
set(LEMON_HAVE_LONG_LONG ${HAVE_LONG_LONG})

find_package(Threads QUIET)

if(NOT LEMON_THREADING)
  if(CMAKE_USE_PTHREADS_INIT)
    set(LEMON_THREADING "Pthread")
  elseif(CMAKE_USE_WIN32_THREADS_INIT)
    set(LEMON_THREADING "Win32")
  else()
    set(LEMON_THREADING "None")
  endif()
endif()

set(LEMON_THREADING "${LEMON_THREADING}" CACHE STRING
  "Threading library: Pthread | Win32 | None." FORCE)

if(LEMON_THREADING STREQUAL "Pthread")
  set(LEMON_USE_PTHREAD TRUE)
elseif(LEMON_THREADING STREQUAL "Win32")
  set(LEMON_USE_WIN32_THREADS TRUE)
endif()

# ---------------------------------------------------------------------------
# Testing
# ---------------------------------------------------------------------------
enable_testing()

if(CMAKE_BUILD_TYPE STREQUAL "Maintainer")
  add_custom_target(check ALL COMMAND ${CMAKE_CTEST_COMMAND})
else()
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
endif()

# ---------------------------------------------------------------------------
# Sub-directories
# ---------------------------------------------------------------------------
add_subdirectory(lemon)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  add_subdirectory(contrib)
  add_subdirectory(demo)
  add_subdirectory(tools)
  add_subdirectory(doc)
  add_subdirectory(test)
endif()

# ---------------------------------------------------------------------------
# CMake config file installation
# ---------------------------------------------------------------------------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/LEMONConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/LEMONConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
  PATH_VARS
    CMAKE_INSTALL_LIBDIR
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/LEMONConfigVersion.cmake"
  VERSION       "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/LEMONConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/LEMONConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# ---------------------------------------------------------------------------
# Source archive target (git-based, no hg)
# ---------------------------------------------------------------------------
string(TOLOWER "${PROJECT_NAME}" _archive_base)
set(_archive_name "${_archive_base}-${PROJECT_VERSION}")

add_custom_target(dist
  COMMAND ${CMAKE_COMMAND} -E remove_directory "${_archive_name}"
  COMMAND git archive --prefix="${_archive_name}/" --format=tar HEAD
          | tar -x
  COMMAND ${CMAKE_COMMAND} -E copy
          cmake/version.cmake "${_archive_name}/cmake/version.cmake"
  COMMAND tar -czf "${_archive_base}-nodoc-${PROJECT_VERSION}.tar.gz" "${_archive_name}"
  COMMAND zip  -r  "${_archive_base}-nodoc-${PROJECT_VERSION}.zip"    "${_archive_name}"
  COMMAND ${CMAKE_COMMAND} -E copy_directory doc/html "${_archive_name}/doc/html"
  COMMAND tar -czf "${_archive_name}.tar.gz" "${_archive_name}"
  COMMAND zip  -r  "${_archive_name}.zip"    "${_archive_name}"
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          doc/html "${_archive_base}-doc-${PROJECT_VERSION}"
  COMMAND tar -czf "${_archive_base}-doc-${PROJECT_VERSION}.tar.gz"
          "${_archive_base}-doc-${PROJECT_VERSION}"
  COMMAND zip  -r  "${_archive_base}-doc-${PROJECT_VERSION}.zip"
          "${_archive_base}-doc-${PROJECT_VERSION}"
  COMMAND ${CMAKE_COMMAND} -E remove_directory "${_archive_name}"
  COMMAND ${CMAKE_COMMAND} -E remove_directory "${_archive_base}-doc-${PROJECT_VERSION}"
  DEPENDS html
  WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
)

unset(_archive_base)
unset(_archive_name)

# ---------------------------------------------------------------------------
# CPack / NSIS (top-level build only)
# ---------------------------------------------------------------------------
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(CPACK_PACKAGE_NAME                 "${PROJECT_NAME}")
  set(CPACK_PACKAGE_VENDOR               "EGRES")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "LEMON - Library for Efficient Modeling and Optimization in Networks")
  set(CPACK_RESOURCE_FILE_LICENSE        "${PROJECT_SOURCE_DIR}/LICENSE")
  set(CPACK_PACKAGE_VERSION              "${PROJECT_VERSION}")
  set(CPACK_PACKAGE_INSTALL_DIRECTORY    "${PROJECT_NAME} ${PROJECT_VERSION}")
  set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${PROJECT_NAME} ${PROJECT_VERSION}")

  set(CPACK_COMPONENTS_ALL headers library html_documentation bin)

  set(CPACK_COMPONENT_HEADERS_DISPLAY_NAME          "C++ headers")
  set(CPACK_COMPONENT_LIBRARY_DISPLAY_NAME           "Dynamic-link library")
  set(CPACK_COMPONENT_BIN_DISPLAY_NAME               "Command line utilities")
  set(CPACK_COMPONENT_HTML_DOCUMENTATION_DISPLAY_NAME "HTML documentation")

  set(CPACK_COMPONENT_HEADERS_DESCRIPTION            "C++ header files")
  set(CPACK_COMPONENT_LIBRARY_DESCRIPTION            "DLL and import library")
  set(CPACK_COMPONENT_BIN_DESCRIPTION                "Command line utilities")
  set(CPACK_COMPONENT_HTML_DOCUMENTATION_DESCRIPTION "Doxygen generated documentation")

  set(CPACK_COMPONENT_HEADERS_DEPENDS  library)
  set(CPACK_COMPONENT_HEADERS_GROUP    "Development")
  set(CPACK_COMPONENT_LIBRARY_GROUP    "Development")
  set(CPACK_COMPONENT_HTML_DOCUMENTATION_GROUP "Documentation")

  set(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION
    "Components needed to develop software using LEMON")
  set(CPACK_COMPONENT_GROUP_DOCUMENTATION_DESCRIPTION
    "Documentation of LEMON")

  set(CPACK_ALL_INSTALL_TYPES Full Developer)
  set(CPACK_COMPONENT_HEADERS_INSTALL_TYPES          Developer Full)
  set(CPACK_COMPONENT_LIBRARY_INSTALL_TYPES          Developer Full)
  set(CPACK_COMPONENT_HTML_DOCUMENTATION_INSTALL_TYPES Full)

  set(CPACK_GENERATOR             "NSIS")
  set(CPACK_NSIS_MUI_ICON         "${PROJECT_SOURCE_DIR}/cmake/nsis/lemon.ico")
  set(CPACK_NSIS_MUI_UNIICON      "${PROJECT_SOURCE_DIR}/cmake/nsis/uninstall.ico")
  set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\lemon.ico")
  set(CPACK_NSIS_DISPLAY_NAME     "${CPACK_PACKAGE_INSTALL_DIRECTORY} ${PROJECT_NAME}")
  set(CPACK_NSIS_HELP_LINK        "http:\\\\\\\\lemon.cs.elte.hu")
  set(CPACK_NSIS_URL_INFO_ABOUT   "http:\\\\\\\\lemon.cs.elte.hu")
  set(CPACK_NSIS_CONTACT          "lemon-user@lemon.cs.elte.hu")
  set(CPACK_NSIS_CREATE_ICONS_EXTRA
    "CreateShortCut \\\"$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\Documentation.lnk\\\" \
\\\"$INSTDIR\\\\share\\\\doc\\\\index.html\\\"")
  set(CPACK_NSIS_DELETE_ICONS_EXTRA
    "!insertmacro MUI_STARTMENU_GETFOLDER Application $MUI_TEMP\n\
    Delete \\\"$SMPROGRAMS\\\\$MUI_TEMP\\\\Documentation.lnk\\\"")

  include(CPack)
endif()
