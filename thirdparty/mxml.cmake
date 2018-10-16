message("****************")
message("-- BUILD MXML --")
message("****************")

# the folder where to compile mxml-2.12
set(MXML_PREFIX mxml)

# set the source location
set(MXML_URL ${CMAKE_CURRENT_LIST_DIR}/mxml-2.12.tar)

# check MD5
set(MXML_URL_MD5  83f7494458b67e4514c8a005b77ef90a)

# build system
set(MXML_MAKE       make)
set(MXML_DIR        ${CMAKE_SOURCE_DIR}/build)
set(MXML_SRC        ${MXML_DIR}/${MXML_PREFIX}/src/${MXML_PREFIX})
set(MXML_CONFIG_OPT "--prefix=${MXML_DIR}")
ExternalProject_Add(${MXML_PREFIX}
    PREFIX              ${MXML_PREFIX}
    URL                 ${MXML_URL}
    URL_MD5             ${MXML_URL_MD5}
    CONFIGURE_COMMAND   ${MXML_SRC}/configure --prefix=${MXML_DIR}
    BUILD_COMMAND       ${MXML_MAKE} -j${NCPU}
	BUILD_IN_SOURCE     1
    INSTALL_COMMAND     ${MXML_MAKE} install
    DEPENDS             ${FFTW_PREFIX}
	LOG_DOWNLOAD        1
	LOG_BUILD           1
)

# post-build setup
set(MXML_INCLUDE_DIRS ${MXML_DIR}/include)
include_directories(${MXML_INCLUDE_DIRS})

# set the library directory variable and link it
set(MXML_LIBRARY_DIRS ${MXML_DIR}/lib)
link_directories(${MXML_LIBRARY_DIRS})
set(MXML_LIBS mxml)
set(MXML_LIBRARY_DIRS ${MXML_LIBRARY_DIRS})

# display info
message("Build MXML in ${MXML_SRC} with CONFIGURE OPTS:")
message(">> ${MXML_CONFIG_OPT}")
message("MXML_INCLUDE_DIRS=${MXML_INCLUDE_DIRS}")
message("MXML_LIBRARY_DIRS=${MXML_LIBRARY_DIRS}")

message("")
