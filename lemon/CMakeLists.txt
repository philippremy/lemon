# ---------------------------------------------------------------------------
# lemon/CMakeLists.txt
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Configure generated headers
# ---------------------------------------------------------------------------
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/config.h"
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/lemon.pc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/lemon.pc"
  @ONLY
)

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------
set(LEMON_SOURCES
  arg_parser.cc
  base.cc
  color.cc
  lp_base.cc
  lp_skeleton.cc
  random.cc
  bits/windows.cc
)

if(LEMON_HAVE_GLPK)
  list(APPEND LEMON_SOURCES glpk.cc)
endif()
if(LEMON_HAVE_CPLEX)
  list(APPEND LEMON_SOURCES cplex.cc)
endif()
if(LEMON_HAVE_CLP)
  list(APPEND LEMON_SOURCES clp.cc)
endif()
if(LEMON_HAVE_CBC)
  list(APPEND LEMON_SOURCES cbc.cc)
endif()
if(LEMON_HAVE_SOPLEX)
  list(APPEND LEMON_SOURCES soplex.cc)
endif()

# ---------------------------------------------------------------------------
# Library type  (contrôlé par BUILD_SHARED_LIBS, standard CMake)
#   cmake -DBUILD_SHARED_LIBS=ON  ...   → shared (.so / .dll)
#   cmake -DBUILD_SHARED_LIBS=OFF ...   → static (.a / .lib)  [défaut]
# ---------------------------------------------------------------------------
add_library(lemon ${LEMON_SOURCES})

# Alias pour que les consumers internes utilisent le même nom que l'extérieur
add_library(LEMON::lemon ALIAS lemon)

# ---------------------------------------------------------------------------
# Propriétés de la cible
# ---------------------------------------------------------------------------
# Note: CXX_VISIBILITY_PRESET hidden n'est pas utilisé ici car les headers
# publics de lemon n'ont pas de macros d'export (LEMON_API) sur leurs
# déclarations. Activer hidden visibility sans annoter les headers cache tous
# les symboles dans la .so → undefined reference à la link.
# Pour une shared lib, on laisse la visibilité par défaut (public).
set_target_properties(lemon PROPERTIES
  OUTPUT_NAME   "emon"
  VERSION       "${PROJECT_VERSION}"
  SOVERSION     "${PROJECT_VERSION_MAJOR}"        # ABI : majeur seulement
  # Position-independent code : obligatoire pour les .so, inoffensif pour .a
  POSITION_INDEPENDENT_CODE ON
)

# ---------------------------------------------------------------------------
# Warning flags — copiés directement depuis lemon_warnings (PRIVATE/BUILD only)
# On n'utilise pas target_link_libraries(lemon PRIVATE lemon_warnings) car
# CMake inclut quand même les INTERFACE targets PRIVATE dans l'export set,
# ce qui provoque : "requires target lemon_warnings that is not in any export set"
# ---------------------------------------------------------------------------
target_compile_options(lemon
  PRIVATE
    $<BUILD_INTERFACE:$<TARGET_PROPERTY:lemon_warnings,INTERFACE_COMPILE_OPTIONS>>
)

# ---------------------------------------------------------------------------
# Includes
# ---------------------------------------------------------------------------
target_include_directories(lemon
  PUBLIC
    # consumers (build tree)
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>   # config.h + lemon_export.h
    # consumers (install tree)
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    # solvers — uniquement nécessaires lors de la compilation de lemon lui-même
    $<$<BOOL:${LEMON_HAVE_GLPK}>:${GLPK_INCLUDE_DIRS}>
    $<$<BOOL:${LEMON_HAVE_CPLEX}>:${ILOG_INCLUDE_DIRS}>
    $<$<BOOL:${LEMON_HAVE_CLP}>:${COIN_INCLUDE_DIRS}>
    $<$<BOOL:${LEMON_HAVE_CBC}>:${COIN_INCLUDE_DIRS}>
    $<$<BOOL:${LEMON_HAVE_SOPLEX}>:${SOPLEX_INCLUDE_DIRS}>
)

# ---------------------------------------------------------------------------
# Link
# ---------------------------------------------------------------------------
target_link_libraries(lemon
  PRIVATE
    $<$<BOOL:${GLPK_LIBRARIES}>:${GLPK_LIBRARIES}>
    $<$<BOOL:${COIN_LIBRARIES}>:${COIN_LIBRARIES}>
    $<$<BOOL:${ILOG_LIBRARIES}>:${ILOG_LIBRARIES}>
    $<$<BOOL:${SOPLEX_LIBRARIES}>:${SOPLEX_LIBRARIES}>
    $<$<BOOL:${LEMON_USE_PTHREAD}>:Threads::Threads>
)

# ---------------------------------------------------------------------------
# Install — cible
# ---------------------------------------------------------------------------

include(GNUInstallDirs)

install(
  TARGETS lemon
  EXPORT  LEMONTargets
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"   COMPONENT library
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"   COMPONENT library
  RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}"   COMPONENT library  # .dll Windows
)

# ---------------------------------------------------------------------------
# Install — headers publics
# ---------------------------------------------------------------------------
install(
  DIRECTORY   "${CMAKE_CURRENT_SOURCE_DIR}/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/lemon"
  COMPONENT   headers
  FILES_MATCHING PATTERN "*.h"
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/config.h"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/lemon"
  COMPONENT   headers
)

# ---------------------------------------------------------------------------
# Install — pkg-config
# ---------------------------------------------------------------------------
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/lemon.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
)

# ---------------------------------------------------------------------------
# Install — export set (consommé par LEMONConfig.cmake)
# ---------------------------------------------------------------------------
install(
  EXPORT      LEMONTargets
  NAMESPACE   LEMON::
  DESTINATION "${CMAKE_INSTALL_DATADIR}/lemon/cmake"
)
