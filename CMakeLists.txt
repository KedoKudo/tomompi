# -- CMAKE Config for TOMOMPI -- #

cmake_minimum_required (VERSION 2.8.8 FATAL_ERROR)

# force to use Homebrew GNU compiler on MacOS
if(APPLE)
    set (MACOSX_RPATH        FALSE )
    set (CMAKE_MACOSX_RPATH  OFF   )
endif(APPLE)


# -- project def --
project (TomoMPI)

# The version number.
set (TomoMPI_VERSION_MAJOR "2")
set (TomoMPI_VERSION_MINOR "0")
set (TomoMPI_VERSION_PATCH "0")

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
set (CMAKE_BUILD_TYPE                  "Release")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE     )

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
include_directories(${HOME}/include)
link_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib
    /usr/lib64
    /usr/lib32
    /usr/lib
    ${HOME}/lib
    ${HOME}/lib64
)

# -- compiler flags
# NOTE:
#  Not the best way to set include directory, but somehow Orthros just keep
#  breaks at the linking for fftw...
set (CMAKE_CXX_FLAGS 
    "-fPIC -std=c++11 -pthread \
    -D__unix__ -DHDF4 -DHDF5 \
    -I${CMAKE_CURRENT_SOURCE_DIR}/include"
)
# -- linker flags
# package specific
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/lib64)
    set (NAPI_LINKER_FLAGS "-L${CMAKE_CURRENT_SOURCE_DIR}/lib64 -lNeXus")
    set (NAPI_LIB "${CMAKE_CURRENT_SOURCE_DIR}/lib64/libNeXus.a")
else()
    set (NAPI_LINKER_FLAGS "-L${CMAKE_CURRENT_SOURCE_DIR}/lib   -lNeXus")
    set (NAPI_LIB "${CMAKE_CURRENT_SOURCE_DIR}/lib/libNeXus.a")
endif()

set (HDF4_LINKER_FLAGS "-L${CMAKE_CURRENT_SOURCE_DIR}/lib -lmfhdf -ldf")
set (MFHDF4_LIB "${CMAKE_CURRENT_SOURCE_DIR}/lib/libmfhdf.a")
set (HDF4_LIB   "${CMAKE_CURRENT_SOURCE_DIR}/lib/libdf.a")

set (JPEG_LINKER_FLAGS "-L${CMAKE_CURRENT_SOURCE_DIR}/lib -ljpeg")
set (JPEG_LIB "${CMAKE_CURRENT_SOURCE_DIR}/lib/libjpeg.a")

set (SZLB_LINKER_FLAGS "-L${CMAKE_CURRENT_SOURCE_DIR}/lib -lsz")
set (SZLB_LIB "${CMAKE_CURRENT_SOURCE_DIR}/lib/libsz.a")

set (HDF5_LINKER_FLAGS "-L${CMAKE_CURRENT_SOURCE_DIR}/lib -lhdf5")
set (HDF5_LIB "${CMAKE_CURRENT_SOURCE_DIR}/lib/libhdf5.a")

set (MXML_LINKER_FLAGS "-L${CMAKE_CURRENT_SOURCE_DIR}/lib -lmxml")
set (MXML_LIB "${CMAKE_CURRENT_SOURCE_DIR}/lib/libmxml.a")

set (FFTW_LINKER_FLAGS "-L${CMAKE_CURRENT_SOURCE_DIR}/lib -lfftw3f")
string (CONCAT THIRDPARTY_LIBS  "${NAPI_LINKER_FLAGS} "
                                "${HDF4_LINKER_FLAGS} "
                                "${JPEG_LINKER_FLAGS} "
                                "${SZLB_LINKER_FLAGS} "
                                "${HDF5_LINKER_FLAGS} "
                                "${MXML_LINKER_FLAGS} "
                                "${FFTW_LINKER_FLAGS} "
)
# set the actual linker flags
set (CMAKE_EXE_LINKER_FLAGS "-lpthread -lm -lz ${THIRDPARTY_LIBS}")

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
target_link_libraries(nexuslibrary 
    ${NAPI_LIB}
    ${MFHDF4_LIB}
    ${HDF4_LIB}
    ${JPEG_LIB}
    ${SZLB_LIB}
    ${HDF5_LIB}
    ${MXML_LIB}
    )
add_custom_command(
    TARGET nexuslibrary POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:nexuslibrary> "${CMAKE_CURRENT_SOURCE_DIR}/lib/${targetfile}"
)

add_library(fft_wrapper STATIC
    src/fft_fftw.c
)
add_dependencies(fft_wrapper ${FFTW_PREFIX})
target_link_libraries(fft_wrapper 
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/libfftw3f.a
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
)
add_dependencies(tomompi nexuslibrary fft_wrapper)
target_link_libraries(tomompi 
    nexuslibrary
    fft_wrapper
    ${MPI_CXX_LIBRARIES}
)

install(TARGETS tomompi nexuslibrary
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin
    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
    ARCHIVE DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
)