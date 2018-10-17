message("****************")
message("-- BUILD SZLB --")
message("****************")

# the folder where to compile mxml-2.12
set(SZLIB_PREFIX szlib)

# set the source location
set(SZLIB_URL ${CMAKE_CURRENT_LIST_DIR}/szip-2.1.1.tar)

# check MD5
set(SZLIB_URL_MD5  5addbf2a5b1bf928b92c47286e921f72)

# build system
set(SZLIB_MAKE       make)
set(SZLIB_DIR        ${CMAKE_SOURCE_DIR}/build)
set(SZLIB_SRC        ${SZLIB_DIR}/${SZLIB_PREFIX}/src/${SZLIB_PREFIX})
set(SZLIB_CONFIG_OPT "--prefix=${SZLIB_DIR}")
ExternalProject_Add(${SZLIB_PREFIX}
    PREFIX              ${SZLIB_PREFIX}
    URL                 ${SZLIB_URL}
    URL_MD5             ${SZLIB_URL_MD5}
    CONFIGURE_COMMAND   ${SZLIB_SRC}/configure --prefix=${SZLIB_DIR}
    BUILD_COMMAND       ${SZLIB_MAKE} -j${NCPU}
	BUILD_IN_SOURCE     1
    INSTALL_COMMAND     ${SZLIB_MAKE} install
    DEPENDS             ${MXML_PREFIX}
	LOG_DOWNLOAD        1
	LOG_BUILD           1
)

# post-build setup
set(SZLIB_INCLUDE_DIRS ${SZLIB_DIR}/include)
include_directories(${SZLIB_INCLUDE_DIRS})

# set the library directory variable and link it
set(SZLIB_LIBRARY_DIRS ${SZLIB_DIR}/lib)
link_directories(${SZLIB_LIBRARY_DIRS})
set(SZLIB_LIBS mxml)
set(SZLIB_LIBRARY_DIRS ${SZLIB_LIBRARY_DIRS})

# display info
message("Build SZLIB in ${SZLIB_SRC} with CONFIGURE OPTS:")
message(">> ${SZLIB_CONFIG_OPT}")
message("SZLIB_INCLUDE_DIRS=${SZLIB_INCLUDE_DIRS}")
message("SZLIB_LIBRARY_DIRS=${SZLIB_LIBRARY_DIRS}")

message("")