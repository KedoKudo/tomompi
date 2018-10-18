message("****************")
message("-- BUILD ZLIB --")
message("****************")

# the folder where to compile zlib
# NOTE:
#  it is recommended to use the system zlib, so this is just a backup plan
set(ZLIB_PREFIX zlib)

# set the source location
set(ZLIB_URL ${CMAKE_CURRENT_LIST_DIR}/zlib-1.2.11.tar)

# check MD5
set(ZLIB_URL_MD5  8efb9889dfcc5c7f083b90013f70d942)

# build system
set(ZLIB_MAKE       make)
set(ZLIB_DIR        ${CMAKE_SOURCE_DIR}/build)
set(ZLIB_SRC        ${ZLIB_DIR}/${ZLIB_PREFIX}/src/${ZLIB_PREFIX})
ExternalProject_Add(${ZLIB_PREFIX}
    PREFIX              ${ZLIB_PREFIX}
    URL                 ${ZLIB_URL}
    URL_MD5             ${ZLIB_URL_MD5}
    CONFIGURE_COMMAND   ${ZLIB_SRC}/configure --prefix=${ZLIB_DIR}
    BUILD_COMMAND       ${ZLIB_MAKE} -j${NCPU}
	BUILD_IN_SOURCE     1
    INSTALL_COMMAND     ${ZLIB_MAKE} install
    DEPENDS             ${MXML_PREFIX}
	LOG_DOWNLOAD        1
	LOG_BUILD           1
)

# post-build setup
set(ZLIB_INCLUDE_DIRS ${ZLIB_DIR}/include)
include_directories(${ZLIB_INCLUDE_DIRS})

# set the library directory variable and link it
set(ZLIB_LIBRARY_DIRS ${ZLIB_DIR}/lib)
link_directories(${ZLIB_LIBRARY_DIRS})
set(ZLIB_LIBS mxml)
set(ZLIB_LIBRARY_DIRS ${ZLIB_LIBRARY_DIRS})
