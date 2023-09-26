cmake_minimum_required(VERSION 3.20)

project(2S-ESS 
  VERSION 2.4.2
  LANGUAGES Fortran)

option (ENABLE_SHARED            "Enable shared libraries" OFF)

enable_testing()


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
include_directories(${CMAKE_BINARY_DIR}/mod_files)
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
message(STATUS "Bin Dir " ${CMAKE_BINARY_DIR})
#=====================================================
# Compiler and Options
#=====================================================
# ... write stuff here
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-fpp -D IFORT -O3 -heap-arrays")
  set(CMAKE_Fortran_FLAGS_DEBUG "-fpp -D IFORT -g -warn all -check all -zero")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-c -std=legacy -mcmodel=large -O3 -ffree-form -ffree-line-length-none -fall-intrinsics -fimplicit-none")
  set(CMAKE_Fortran_FLAGS_DEBUG "-c -std=legacy -mcmodel=large -g -O0 -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow -Wall -fcheck=all -fbacktrace -fall-intrinsics -Werror=underflow")
endif()


#=====================================================
# twostr
#=====================================================
file(GLOB REG_2S_SRC sourcecode/2stream_2p4_feb15/*.f90)

#=====================================================
# ESS
#=====================================================
set(FO_PATH_1P4 sourcecode/FO_1p4_sourcecode)

# set(REG_FO_SRC 
#    ${FO_PATH_1P4}/FO_geometry_Generic.f90 
#    ${FO_PATH_1P4}/FO_geometry_Routines.f90 
#    ${FO_PATH_1P4}/FO_SSgeometry_master.f90
#    ${FO_PATH_1P4}/FO_ScalarSS_Spherfuncs.f90
#    ${FO_PATH_1P4}/FO_ScalarSS_RTCalcs_I.f90
# )

file(GLOB REG_FO_SRC ${FO_PATH_1P4}/*.f90)

#=====================================================
# twostr library (static)
#=====================================================
add_library(TWOSTR STATIC ${REG_2S_SRC})

set(MODLIST)
foreach (_fsource ${REG_2S_SRC})
  # extract program name from source
  execute_process(COMMAND grep -i -P "^\\s*module\\s+(\\w*)" ${_fsource}
                  COMMAND awk "{ print tolower(\$2) }"
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                  OUTPUT_VARIABLE processOutput
                  )

  # strip terminating newline
  string(REGEX REPLACE "\n$" "" modName "${processOutput}")

  if (NOT ${modName} STREQUAL "")
      #message( "Found ${modName} in ${_fsource}" )
      list(APPEND MODLIST "${CMAKE_BINARY_DIR}/${modName}.mod")
  endif ()
endforeach(_fsource ${REG_2S_SRC})
list(REMOVE_DUPLICATES MODLIST)
#message(STATUS "arg list: ${MODLIST}")

set_target_properties (TWOSTR
      PROPERTIES
      PUBLIC_HEADER "${MODLIST}"
      POSITION_INDEPENDENT_CODE ON)

install(
  TARGETS TWOSTR
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )

#=====================================================
# twostr library (Shared, optional)
#=====================================================
if(ENABLE_SHARED)
  add_library(TWOSTR_shared SHARED ${REG_2S_SRC})
  #target_link_libraries(TWOSTR_shared PRIVATE ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})
  install(
    TARGETS TWOSTR_shared
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    )
endif()

#=====================================================
# ESS library (static)
#=====================================================
add_library(ESS STATIC ${REG_FO_SRC})

set(MODLIST_FO)
foreach (_fsource ${REG_FO_SRC})
  # extract program name from source
  execute_process(COMMAND grep -i -P "^\\s*module\\s+(\\w*)" ${_fsource}
                  COMMAND awk "{ print tolower(\$2) }"
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                  OUTPUT_VARIABLE processOutput
                  )

  # strip terminating newline
  string(REGEX REPLACE "\n$" "" modName "${processOutput}")

  if (NOT ${modName} STREQUAL "")
      #message( "Found ${modName} in ${_fsource}" )
      list(APPEND MODLIST_FO "${CMAKE_BINARY_DIR}/${modName}.mod")
  endif ()
endforeach(_fsource ${REG_FO_SRC})
list(REMOVE_DUPLICATES MODLIST)
#message(STATUS "arg list: ${MODLIST}")

set_target_properties (ESS
      PROPERTIES
      PUBLIC_HEADER "${MODLIST_FO}"
      POSITION_INDEPENDENT_CODE ON)

install(
  TARGETS ESS
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )

#=====================================================
# ESS library (Shared, optional)
#=====================================================
if(ENABLE_SHARED)
  add_library(ESS_shared SHARED ${REG_FO_SRC})
  #target_link_libraries(TWOSTR_shared PRIVATE ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})
  install(
    TARGETS ESS_shared
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    )
endif()

#====================================================
# Add test
#====================================================
set(EXE regtest)
add_executable(${EXE} TEST/test_regular.f90)
target_link_libraries(${EXE} TWOSTR ESS)
install(TARGETS ${EXE} RUNTIME DESTINATION bin)

add_test(
  NAME "regular run"
  COMMAND $<TARGET_FILE:${EXE}>
)
set_tests_properties("regular run" PROPERTIES LABELS "run")

add_test(
  NAME "regular output"
  COMMAND numdiff -S -r 1e-14 ${CMAKE_SOURCE_DIR}/TEST/results_saved_output.dat test_regular_output.dat 
)
set_tests_properties("regular output" PROPERTIES LABELS "test")