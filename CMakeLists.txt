# -- CMAKE Config for TOMOMPI -- #

cmake_minimum_required (VERSION 2.8.8 FATAL_ERROR)

# force to use Homebrew GNU compiler
if(APPLE)
    set(CMAKE_C_COMPILER   "gcc-7")
    set(CMAKE_CXX_COMPILER "g++-7")
    set(CMAKE_LINKER       "gcc-7")
endif(APPLE)

# project def
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

# add the sub-directories that we want CMake to scan
add_subdirectory(include)
add_subdirectory(thirdparty)
add_subdirectory(src)

# build type release by defualt
set(CMAKE_BUILD_TYPE "Release")

# External Dependencies 
# -- MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
# -- Parallel built CPU number 
set(NCPU 1)
# -- FFTW
include(thirdparty/fftw3.cmake)
# -- MXML
include(thirdparty/mxml.cmake)
# -- HDF4
include(thirdparty/hdf4.cmake)
# -- HDF5
include(thirdparty/hdf5.cmake)
# -- NAPI
include(thirdparty/napi.cmake)


# Add the files in the include directory to be included
# while compiling all the source code header files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

