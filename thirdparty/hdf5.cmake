message("****************")
message("-- BUILD HDF5 --")
message("****************")

# the folder where to compile HDF5-1.8.10
# NOTE:
#   TomoMPI is using some deprecated API when calling HDF5, therefore we need
#   to use a relatively older version of HDF5 here.
set(HDF5_PREFIX hdf5)

# set the source location
set(HDF5_URL ${CMAKE_CURRENT_LIST_DIR}/hdf5-1.8.10.tar.gz)

# check MD5
set(HDF5_URL_MD5  95da64de3113d940acbe11b9f2f7a67a)

# build system
set(HDF5_MAKE       make)
set(HDF5_DIR        ${CMAKE_SOURCE_DIR})
set(HDF5_SRC        ${HDF5_DIR}/build/${HDF5_PREFIX}/src/${HDF5_PREFIX})
set(HDF5_CONFIG_OPT "--with-szlib=${HDF5_DIR}/lib --enable-cxx --prefix=${HDF5_DIR}")
ExternalProject_Add(${HDF5_PREFIX}
    PREFIX              ${HDF5_PREFIX}
    URL                 ${HDF5_URL}
    URL_MD5             ${HDF5_URL_MD5}
    CONFIGURE_COMMAND   ${HDF5_SRC}/configure --with-szlib=${HDF5_DIR}/lib --enable-cxx --prefix=${HDF5_DIR}
    BUILD_COMMAND       ${HDF5_MAKE} -j${NCPU}
	BUILD_IN_SOURCE     1
    INSTALL_COMMAND     ${HDF5_MAKE} install
    DEPENDS             ${HDF4_PREFIX}
	LOG_DOWNLOAD        1
	LOG_BUILD           1
)

# post-build setup
set(HDF5_INCLUDE_DIRS ${HDF5_DIR}/include)
include_directories(${HDF5_INCLUDE_DIRS})

# set the library directory variable and link it
set(HDF5_LIBRARY_DIRS ${HDF5_DIR}/lib)
link_directories(${HDF5_LIBRARY_DIRS})
set(HDF5_LIBS mxml)
set(HDF5_LIBRARY_DIRS ${HDF5_LIBRARY_DIRS})

# display info
message("Build HDF5 in ${HDF5_SRC} with CONFIGURE OPTS:")
message(">> ${HDF5_CONFIG_OPT}")
message("HDF5_INCLUDE_DIRS=${HDF5_INCLUDE_DIRS}")
message("HDF5_LIBRARY_DIRS=${HDF5_LIBRARY_DIRS}")

message("")
