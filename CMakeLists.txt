# -- CMAKE Config for TOMOMPI -- #

cmake_minimum_required (VERSION 2.8.8 FATAL_ERROR)
project (TomoMPI)

# The version number.
SET (TomoMPI_VERSION_MAJOR "2")
SET (TomoMPI_VERSION_MINOR "0")
SET (TomoMPI_VERSION_PATCH "0")

# Auto check
include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckIncludeFiles)
include(ExternalProject)

# build type release by defualt
set(CMAKE_BUILD_TYPE "Release")

# External Dependencies 
# -- MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
# -- FFTW
include(thirdparty/fftw3.cmake)
# -- NAPI
include(thirdparty/napi.cmake)
# -- MXML
include(thirdparty/mxml.cmake)
# -- HDF5
include(thirdparty/hdf5.cmake)
# -- HDF4
include(thirdparty/hdf4.cmake)

# Add the files in the include directory to be included
# while compiling all the source code header files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# add the sub-directories that we want CMake to scan
add_subdirectory(include)
add_subdirectory(thirdparty)
add_subdirectory(src)
