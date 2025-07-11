INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
)

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/lemon.pc.in
  ${CMAKE_CURRENT_BINARY_DIR}/lemon.pc
  @ONLY
)

SET(LEMON_SOURCES
  arg_parser.cc
  base.cc
  color.cc
  lp_base.cc
  lp_skeleton.cc
  random.cc
  bits/windows.cc
)

IF(LEMON_HAVE_GLPK)
  SET(LEMON_SOURCES ${LEMON_SOURCES} glpk.cc)
  INCLUDE_DIRECTORIES(${GLPK_INCLUDE_DIRS})
  IF(WIN32)
    INSTALL(FILES ${GLPK_BIN_DIR}/glpk.dll DESTINATION bin)
    INSTALL(FILES ${GLPK_BIN_DIR}/libltdl3.dll DESTINATION bin)
    INSTALL(FILES ${GLPK_BIN_DIR}/zlib1.dll DESTINATION bin)
  ENDIF()
ENDIF()

IF(LEMON_HAVE_CPLEX)
  SET(LEMON_SOURCES ${LEMON_SOURCES} cplex.cc)
  INCLUDE_DIRECTORIES(${ILOG_INCLUDE_DIRS})
ENDIF()

IF(LEMON_HAVE_CLP)
  SET(LEMON_SOURCES ${LEMON_SOURCES} clp.cc)
  INCLUDE_DIRECTORIES(${COIN_INCLUDE_DIRS})
ENDIF()

IF(LEMON_HAVE_CBC)
  SET(LEMON_SOURCES ${LEMON_SOURCES} cbc.cc)
  INCLUDE_DIRECTORIES(${COIN_INCLUDE_DIRS})
ENDIF()

IF(LEMON_HAVE_SOPLEX)
  SET(LEMON_SOURCES ${LEMON_SOURCES} soplex.cc)
  INCLUDE_DIRECTORIES(${SOPLEX_INCLUDE_DIRS})
ENDIF()

if(LEMON_BUILD_SHARED_LIBRARY)
    ADD_LIBRARY(lemon SHARED ${LEMON_SOURCES})
else()
    ADD_LIBRARY(lemon STATIC ${LEMON_SOURCES})
endif()

TARGET_LINK_LIBRARIES(lemon
  ${GLPK_LIBRARIES} ${COIN_LIBRARIES} ${ILOG_LIBRARIES} ${SOPLEX_LIBRARIES} ${ZLIB_LIBRARIES}
)

IF(UNIX)
  SET_TARGET_PROPERTIES(lemon PROPERTIES OUTPUT_NAME lemon VERSION ${LEMON_VERSION} SOVERSION ${LEMON_VERSION} MACOSX_RPATH ON INSTALL_NAME_DIR "@rpath")
ENDIF()

INSTALL(
  TARGETS lemon
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  COMPONENT library
)

INSTALL(
  DIRECTORY . bits concepts
  DESTINATION include/lemon
  COMPONENT headers
  FILES_MATCHING PATTERN "*.h"
)

INSTALL(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/config.h
  DESTINATION include/lemon
  COMPONENT headers
)

INSTALL(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/lemon.pc
  DESTINATION lib/pkgconfig
)
