cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(PROJECT_OT)
find_package(ROOT)
message("FOUND ROOT >>>>>>" ${ROOT_FOUND})
#----------------------------------------------------------------------------

set(useROOT true)
set(OT_SOURCE_DIRECTORIES src)
set(OT_EXECUTABLES lampslow)
set(OT_SCRIPTS vis.mac)

#----------------------------------------------------------------------------

option(WITH_GEANT4_UIVIS "Build lampslow with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

if(useROOT)
  list(APPEND OT_SOURCE_DIRECTORIES)
  list(APPEND OT_EXECUTABLES)
endif(useROOT)

set(OT_SOURCES_DIRECTORIES_FULL)
foreach(_directory ${OT_SOURCE_DIRECTORIES})
  list(APPEND OT_SOURCE_DIRECTORIES_FULL ${PROJECT_SOURCE_DIR}/${_directory})
endforeach(_directory)

set(OT_SOURCES)
foreach(_directory ${OT_SOURCE_DIRECTORIES_FULL})
  list(APPEND OT_SOURCES ${_directory}/*.cc)
endforeach(_directory)

set(OT_HEADERS)
foreach(_directory ${OT_SOURCE_DIRECTORIES_FULL})
  list(APPEND OT_HEADERS ${_directory}/../include/*.hh)
endforeach(_directory)

foreach(_a ${OT_SOURCES})
  message(${_a})
endforeach(_a)

include(${Geant4_USE_FILE})
include_directories(${OT_SOURCE_DIRECTORIES_FULL}/../include)

#if(useROOT)
#find_package(ROOT REQUIRED)

#include_directories(${ROOT_INCLUDE_DIR})
#endif(useROOT)

if(useROOT)
  EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
endif(useROOT)

file(GLOB sources ${OT_SOURCES})
file(GLOB headers ${OT_HEADERS})

foreach(_executable ${OT_EXECUTABLES})
  add_executable(${_executable} ${_executable}.cc ${sources} ${headers})
if(useROOT)
  EXECUTE_PROCESS(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CMAKE_EXE_LINKER_FLAGS ${ROOT_LD_FLAGS})
endif(useROOT)
  target_link_libraries(${_executable} ${Geant4_LIBRARIES})
endforeach(_executable)

foreach(_script ${OT_SCRIPTS})
  configure_file( ${PROJECT_SOURCE_DIR}/${_script} ${PROJECT_BINARY_DIR}/${_script} COPYONLY)
endforeach(_script)

add_custom_target(PROJECT_OT DEPENDS ${OT_EXECUTABLES})

foreach(_executable ${OT_EXECUTABLES})
  install(TARGETS ${_executable} DESTINATION bin)
endforeach(_executable)
