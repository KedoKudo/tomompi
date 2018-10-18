# -- CMAKE Config for TOMOMPI -- #

cmake_minimum_required (VERSION 2.8.8 FATAL_ERROR)

# force to use Homebrew GNU compiler on MacOS
if(APPLE)
    SET (MACOSX_RPATH        FALSE )
    SET (CMAKE_MACOSX_RPATH  OFF   )
endif(APPLE)

# project def
project (TomoMPI)

# The version number.
SET (TomoMPI_VERSION_MAJOR "2")
SET (TomoMPI_VERSION_MINOR "0")
SET (TomoMPI_VERSION_PATCH "0")

# stop build shared libraries
SET (BUILD_SHARED_LIBS OFF)

# Auto check and Exteral project building support
include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckIncludeFiles)
include(ExternalProject)

# add the sub-directories that we want CMake to scan
add_subdirectory(include)
add_subdirectory(thirdparty)
add_subdirectory(src)

# build type release by defualt
SET (CMAKE_BUILD_TYPE "Release")

# -- External Dependencies -- # 
# -- MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
# -- Parallel built CPU number 
set(NCPU 4)
# -- FFTW
include(thirdparty/fftw3.cmake)
# -- MXML
include(thirdparty/mxml.cmake)
# -- ZLIB
# include(thirdparty/zlib.cmake)
# -- JPEG
include(thirdparty/jpeg.cmake)
# -- szlib
include(thirdparty/szlib.cmake)
# -- HDF4
include(thirdparty/hdf4.cmake)
# -- HDF5
include(thirdparty/hdf5.cmake)
# -- NAPI
include(thirdparty/napi.cmake)

# -- BUILD TOMOMPI -- #
# Add the files in the include directory to be included
# while compiling all the source code header files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/nexus)

# -- build Nexuslibrary
add_library(nexuslibrary STATIC 
    src/nexusattribute.cpp
    src/nexusbox.cpp 
    src/nexusapi.cpp
    src/nexusdata.cpp
    src/nexusexceptionclass.cpp
    src/nexusfield.cpp
    src/nexusgroup.cpp
)
set_property(
    TARGET nexuslibrary 
    PROPERTY COMPILE_FLAGS "-D__unix__"
)
