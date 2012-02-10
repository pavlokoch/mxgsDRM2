cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(mxgsDRM)

#SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
#SET(LINK_SEARCH_END_STATIC ON)
#SET(LINK_SEARCH_START_STATIC ON)
#SET(CMAKE_EXE_LINKER_FLAGS -static)

# messing about with rpath.
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build mxgsDRM with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})
#set(LIBS ${LIBS} ${GSL_LIBRARIES})

set(GDML_FILES
  aces.gdml
  asim.gdml
  bgo.gdml
  bgoSingle.gdml
  codedMask.gdml
  columbus.gdml
  czt.gdml
  instrument.gdml
  mmia.gdml
  mxgs.gdml
  )

foreach(_gdml ${GDML_FILES})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_gdml}
    ${PROJECT_BINARY_DIR}/${_gdml}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(mxgsDRM mxgsDRM.cc ${sources} ${headers})
target_link_libraries(mxgsDRM ${Geant4_LIBRARIES})
target_link_libraries(mxgsDRM ${GSL_LIBRARIES})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS mxgsDRM DESTINATION bin)


