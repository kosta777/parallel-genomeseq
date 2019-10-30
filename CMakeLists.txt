cmake_minimum_required(VERSION 3.12.0)
include("cmake/warning_flags.cmake")
option(HAS_CLANG_TIDY "use clang-tidy" OFF)
include("cmake/clang_tidy.cmake")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CC_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb ")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall -std=c++14")

project(GENOMESEQ)
set(TARGET_NAME parseqal)
option(DEBUG "debug mode" OFF)
set(CMAKE_BUILD_TYPE "Release")

if (DEBUG)
  set(CMAKE_BUILD_TYPE "Debug")
  set(TARGET_NAME parseqal_debug)
endif ()


aux_source_directory(${PROJECT_SOURCE_DIR}/src/aligner ALIGNER_SRCS)
add_library(libaligner ${ALIGNER_SRCS})
include_directories(${PROJECT_SOURCE_DIR}/src/aligner)

add_executable(${TARGET_NAME} src/main.cpp)
target_link_libraries(${TARGET_NAME} PRIVATE libaligner)
target_link_libraries(${TARGET_NAME} PRIVATE warning_flags)
#target_compile_definitions(${TARGET_NAME} PUBLIC -D_USE_MATH_DEFINES)

add_executable(bwt_fm src/bwt_fm.cpp)
target_link_libraries(bwt_fm PRIVATE libaligner)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

if(CLANG_TIDY_EXE)
  set_target_properties(${TARGET_NAME} PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif()