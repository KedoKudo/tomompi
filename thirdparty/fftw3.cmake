# FFTW Build Info
message("****************")
message("-- BUILD FFTW --")
message("****************")

# the folder where to compile FFTW3.3.8
# NOTE:
#  Only the float precision of FFTW is called from TomoMPI, therefore we only
#  build the float precision here.
set(FFTW_PREFIX fftw)

# set the source location
set(FFTW_URL ${CMAKE_CURRENT_LIST_DIR}/fftw-3.3.8.tar.gz)

# check MD5
set(FFTW_URL_MD5  8aac833c943d8e90d51b697b27d4384d)

# FFTW examples
# -- somehow the examples are not in the tarball
set(BUILD_FFTW_EXAMPLES FALSE) 
set(FFTW_EXAMPLES_STEP ${FFTW_PREFIX}_examples)

# build system
set(FFTW_MAKE       make)
set(FFTW_DIR        ${CMAKE_SOURCE_DIR})
set(FFTW_SRC        ${FFTW_DIR}/build/${FFTW_PREFIX}/src/${FFTW_PREFIX})
set(FFTW_CONFIG_OPT "--enable-float --enable-mpi --prefix=${FFTW_DIR}")
# add instructions to build the FFTW source
# -- build float precision (required by napi)
ExternalProject_Add(${FFTW_PREFIX}
    PREFIX              ${FFTW_PREFIX}
    URL                 ${FFTW_URL}
    URL_MD5             ${FFTW_URL_MD5}
    CONFIGURE_COMMAND   ${FFTW_SRC}/configure --enable-float --enable-mpi --prefix=${FFTW_DIR}
    BUILD_COMMAND       ${FFTW_MAKE} -j${NCPU}
	BUILD_IN_SOURCE     1
	INSTALL_COMMAND     ${FFTW_MAKE} install
	LOG_DOWNLOAD        1
	LOG_BUILD           1
)

# set the include directory variable and include it
set(FFTW_INCLUDE_DIRS ${FFTW_DIR}/include)
include_directories(${FFTW_INCLUDE_DIRS})

# set the library directory variable and link it
set(FFTW_LIBRARY_DIRS ${FFTW_DIR}/lib)
link_directories(${FFTW_LIBRARY_DIRS})
set(FFTW_LIBS fftw3f)
set(FFTW_LIBRARY_DIRS ${FFTW_LIBRARY_DIRS})

# display info
message("Build FFTW in ${FFTW_SRC} with CONFIGURE OPTS:")
message(">> ${FFTW_CONFIG_OPT}")
message("FFTW_INCLUDE_DIRS=${FFTW_INCLUDE_DIRS}")
message("FFTW_LIBRARY_DIRS=${FFTW_LIBRARY_DIRS}")

message("")
