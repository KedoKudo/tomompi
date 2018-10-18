message("****************")
message("-- BUILD JPEG --")
message("****************")

# the folder where to compile jpeg
# NOTE:
#  Although it is recommended to use system JPEG, building HDF4 on Orthros 
#  using system JPEG somehow leads to strange linking error.  Therefore we 
#  use a customized JPEG build here.
set(JPEG_PREFIX jpeg)

# set the source location
set(JPEG_URL ${CMAKE_CURRENT_LIST_DIR}/jpegsrc.v9c.tar)

# check MD5
set(JPEG_URL_MD5  fae140877c90c3b58dfb7b4ed4814b4a)

# build system
set(JPEG_MAKE       make)
set(JPEG_DIR        ${CMAKE_SOURCE_DIR})
set(JPEG_SRC        ${JPEG_DIR}/build/${JPEG_PREFIX}/src/${JPEG_PREFIX})
ExternalProject_Add(${JPEG_PREFIX}
    PREFIX              ${JPEG_PREFIX}
    URL                 ${JPEG_URL}
    URL_MD5             ${JPEG_URL_MD5}
    CONFIGURE_COMMAND   ${JPEG_SRC}/configure --with-pic --prefix=${JPEG_DIR}
    BUILD_COMMAND       ${JPEG_MAKE} -j${NCPU}
	BUILD_IN_SOURCE     1
    INSTALL_COMMAND     ${JPEG_MAKE} install
    DEPENDS             ${MXML_PREFIX}
	LOG_DOWNLOAD        1
	LOG_BUILD           1
)

# post-build setup
set(JPEG_INCLUDE_DIRS ${JPEG_DIR}/include)
include_directories(${JPEG_INCLUDE_DIRS})

# set the library directory variable and link it
set(JPEG_LIBRARY_DIRS ${JPEG_DIR}/lib)
link_directories(${JPEG_LIBRARY_DIRS})
set(JPEG_LIBS mxml)
set(JPEG_LIBRARY_DIRS ${JPEG_LIBRARY_DIRS})
