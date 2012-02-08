cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(mxgsDRM)

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
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

set(GDML_FILES
  aces.gdml
  asim.gdml
  bgo.gdml
  bgoSingle.gdml
  codedMask.gdml
  columbus.gdml
  columbus_bk.gdml
  czt.gdml
  instrument.gdml
  iss_C2I1.gdml
  iss_C3I2.gdml
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


