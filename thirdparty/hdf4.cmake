message("****************")
message("-- BUILD HDF4 --")
message("****************")

# before building HDF4, szlib needs to be built
message("HDF4 requires szlib")
# -- szlib
include(${CMAKE_CURRENT_LIST_DIR}/szlib.cmake)

# the folder where to compile mxml-2.12
set(HDF4_PREFIX hdf4)

# set the source location
set(HDF4_URL ${CMAKE_CURRENT_LIST_DIR}/hdf-4.2.13.tar.bz2)

# check MD5
set(HDF4_URL_MD5  2c1b6c7fdf97738251154680b37bd86a)

# build system
set(HDF4_MAKE       make)
set(NCPU            4   )
set(HDF4_DIR        ${CMAKE_SOURCE_DIR}/build/)
set(HDF4_SRC        ${HDF4_DIR}/${HDF4_PREFIX}/src/${HDF4_PREFIX})
set(HDF4_CONFIG_OPT "--with-szlib=${HDF4_DIR}/lib --prefix=${HDF4_DIR}")
ExternalProject_Add(${HDF4_PREFIX}
    PREFIX              ${HDF4_PREFIX}
    URL                 ${HDF4_URL}
    URL_MD5             ${HDF4_URL_MD5}
    CONFIGURE_COMMAND   ${HDF4_SRC}/configure --with-szlib=${HDF4_DIR}/lib --enable-cxx --prefix=${HDF4_DIR}
    BUILD_COMMAND       ${HDF4_MAKE} -j${NCPU}
	BUILD_IN_SOURCE     1
	INSTALL_COMMAND     ${HDF4_MAKE} install
	LOG_DOWNLOAD        1
	LOG_BUILD           1
)

message("")
