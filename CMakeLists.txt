# -- CMAKE Config for TOMOMPI -- #

cmake_minimum_required (VERSION 2.8.8 FATAL_ERROR)

# force to use Homebrew GNU compiler on MacOS
if(APPLE)
    SET (MACOSX_RPATH        FALSE )
    SET (CMAKE_MACOSX_RPATH  OFF   )
endif(APPLE)


# -- project def --
project (TomoMPI)

# The version number.
SET (TomoMPI_VERSION_MAJOR "2")
SET (TomoMPI_VERSION_MINOR "0")
SET (TomoMPI_VERSION_PATCH "0")

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
SET (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# -- External Dependencies -- # 
# -- MPI
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_DIRS})
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
link_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib
    /usr/lib64
    /usr/lib
)

# -- compiler flags
set (CMAKE_CXX_FLAGS 
    "-fPIC -D__unix__ -DHDF4 -DHDF5"
)
# -- linker flags
set (CMAKE_EXE_LINKER_FLAGS
    "-lNeXus -lpthread -lm  -lfftw3f -lmfhdf -ldf -lhdf5 -lmxml"
)

message(${CMAKE_EXE_LINKER_FLAGS})

add_library(nexuslibrary STATIC
    src/nexusapi.cpp    
    src/nexusattribute.cpp
    src/nexusbox.cpp 
    src/nexusdata.cpp
    src/nexusexceptionclass.cpp
    src/nexusfield.cpp
    src/nexusgroup.cpp
)
add_dependencies(nexuslibrary ${NAPI_PREFIX})
add_custom_command(
    TARGET nexuslibrary POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:nexuslibrary> "${CMAKE_CURRENT_SOURCE_DIR}/lib/${targetfile}"
)

add_executable(tomompi
    src/tomompi.cpp
    src/tomompi_client.cpp
    src/tomompi_sinogramserver.cpp
    src/centeringclass.cpp
    src/errorlogclass.cpp
    src/filteredbackprojection.cpp
    src/gridrec.cpp
    src/linkedlistclass.cpp
    src/logfileclass.cpp
    src/recon_algorithm.cpp
    src/fft_fftw.c
)
add_dependencies(tomompi nexuslibrary)
target_link_libraries(tomompi 
    nexuslibrary
    ${MPI_CXX_LIBRARIES}
)

install(TARGETS tomompi nexuslibrary
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin
    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
    ARCHIVE DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
)