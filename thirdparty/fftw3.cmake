# FFTW Build Info
message("****************")
message("-- BUILD FFTW --")
message("****************")

# the folder where to compile FFTW3.3.8
set(FFTW_PREFIX fftw)

# set the source location
set(FFTW_URL ${CMAKE_CURRENT_LIST_DIR}/fftw-3.3.8.tar.gz)

# check MD5
set(FFTW_URL_MD5  8aac833c943d8e90d51b697b27d4384d)

# FFTW examples
# -- somehow the examples are not in the tarball
set(BUILD_FFTW_EXAMPLES FALSE) 
set(FFTW_EXAMPLES_STEP ${FFTW_PREFIX}_examples)

# make system
set(FFTW_MAKE make)
set(NCPU      4   )

set(FFTW_SRC ${CMAKE_SOURCE_DIR}/build/${FFTW_PREFIX}/src/${FFTW_PREFIX})
# add instructions to build the FFTW source
# -- build float precision (required by napi)
ExternalProject_Add(${FFTW_PREFIX}
    PREFIX              ${FFTW_PREFIX}
    URL                 ${FFTW_URL}
    URL_MD5             ${FFTW_URL_MD5}
    CONFIGURE_COMMAND   ${FFTW_SRC}/configure --enable-float --enable-mpi
    BUILD_COMMAND       ${FFTW_MAKE} -j${NCPU}
	BUILD_IN_SOURCE     1
	INSTALL_COMMAND     ""
	LOG_DOWNLOAD        1
	LOG_BUILD           1
)

# get the unpacked source directory path
ExternalProject_Get_Property(${FFTW_PREFIX} SOURCE_DIR)
message(STATUS "Source directory of ${FFTW_PREFIX} ${SOURCE_DIR}")

# build example
if (BUILD_FFTW_EXAMPLES)
    ExternalProject_Add_Step(${FFTW_PREFIX} 
        ${FFTW_PREFIX}_examples
	    COMMAND make -j${NCPU} examples fftw_build_prefix=${FFTW_PREFIX}
	    DEPENDEES build
	    WORKING_DIRECTORY ${SOURCE_DIR}
	    LOG 1
	)
endif (BUILD_FFTW_EXAMPLES)

# build both debug and release
set(FFTW_DEBUG_DIR      ${SOURCE_DIR}/build/${FFTW_PREFIX}_debug)
set(FFTW_RELEASE_DIR    ${SOURCE_DIR}/build/${FFTW_PREFIX}_release)
message(STATUS "FFTW Debug directory ${FFTW_DEBUG_DIR}")
message(STATUS "FFTW Release directory ${FFTW_RELEASE_DIR}")

# set the include directory variable and include it
set(FFTW_INCLUDE_DIRS ${SOURCE_DIR}/include)
include_directories(${FFTW_INCLUDE_DIRS})

# link the correct FFTW directory when the project is 
# in Debug or Release mode
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    link_directories(${FFTW_RELEASE_DIR})
	set(FFTW_LIBS fftw_debug fftwmalloc_debug)
	set(FFTW_LIBRARY_DIRS ${FFTW_DEBUG_DIR})
else (CMAKE_BUILD_TYPE STREQUAL "Debug")
	# in Release mode
	link_directories(${FFTW_RELEASE_DIR})
	set(FFTW_LIBS fftw fftwmalloc)
	set(FFTW_LIBRARY_DIRS ${FFTW_RELEASE_DIR})
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

# verify that the FFTW header files can be included
set(CMAKE_REQUIRED_INCLUDES_SAVE ${CMAKE_REQUIRED_INCLUDES})
set(CMAKE_REQUIRED_INCLUDES      ${CMAKE_REQUIRED_INCLUDES} ${FFTW_INCLUDE_DIRS})
check_include_file_cxx("fftw/fftw.h" HAVE_FFTW)
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES_SAVE})
if (NOT HAVE_FFTW)
	message(STATUS "Did not build FFTW correctly as cannot find fftw.h. Will build it.")
	set(HAVE_FFTW 1)
endif (NOT HAVE_FFTW)

message("")
