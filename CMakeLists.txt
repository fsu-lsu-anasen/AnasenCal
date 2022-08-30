cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE On)

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE "Release")
    message("Building for release")
else()
    message("Building for debug")
endif()

project(AnasenCal)

find_package(ROOT REQUIRED COMPONENTS Spectrum GenVector)

set(ANASENCAL_BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(ANASENCAL_LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)

add_subdirectory(src)
