cmake_minimum_required(VERSION 3.12.0)
include("cmake/warning_flags.cmake")
option(HAS_CLANG_TIDY "use clang-tidy" OFF)
include("cmake/clang_tidy.cmake")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CC_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -pedantic -Wextra -g -ggdb ")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -Ofast -march=native -Wall -std=c++17")

project(GENOMESEQ)
include("cmake/GetEigen.cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(TARGET_NAME parseqal)
option(DEBUG "debug mode" OFF)
set(CMAKE_BUILD_TYPE "Release")
include("cmake/GetEigen.cmake")

if (DEBUG)
  set(CMAKE_BUILD_TYPE "Debug")
  set(TARGET_NAME parseqal_debug)
  add_definitions(-DDEBUG)
endif ()

aux_source_directory(${PROJECT_SOURCE_DIR}/src/aligner ALIGNER_SRCS)
add_library(libaligner ${ALIGNER_SRCS})
target_link_libraries(libaligner PUBLIC Eigen)
include_directories(${PROJECT_SOURCE_DIR}/src/aligner)

# googletest
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

add_executable(${TARGET_NAME} src/main.cpp)
target_link_libraries(${TARGET_NAME} PUBLIC Eigen)
target_link_libraries(${TARGET_NAME} PRIVATE libaligner)
target_link_libraries(${TARGET_NAME} PRIVATE warning_flags)
target_link_libraries(${TARGET_NAME} PRIVATE gtest_main)
#target_compile_definitions(${TARGET_NAME} PUBLIC -D_USE_MATH_DEFINES)

add_executable(bwt_fm src/bwt_fm.cpp)
target_link_libraries(bwt_fm PUBLIC Eigen)
target_link_libraries(bwt_fm PRIVATE libaligner)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

if(CLANG_TIDY_EXE)
  set_target_properties(${TARGET_NAME} PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif()
