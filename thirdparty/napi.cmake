message("****************")
message("-- BUILD NAPI --")
message("****************")

# the folder where to compile napi
set(NAPI_PREFIX napi)

# set the source location
# NOTE: 
#  The original napi project is no longer supported, which leaves some 
#  broken cmake options. 
#  The tar ball used here has been modified to address this issue, and 
#  building process has been tested and passed on both MacOS (local) and 
#  CentOS7 (Orthros).
set(NAPI_URL ${CMAKE_CURRENT_LIST_DIR}/napi-v4.4.3.tar.gz) 


# check MD5
set(NAPI_URL_MD5 ddd8886196f37b3b067c1b80ec7e2e7d)

# build system
set(NAPI_MAKE       make)
set(NAPI_DIR        ${CMAKE_SOURCE_DIR}/build)
set(NAPI_SRC        ${NAPI_DIR}/${NAPI_PREFIX}/src/${NAPI_PREFIX})
ExternalProject_Add(${NAPI_PREFIX}
    PREFIX              ${NAPI_PREFIX}
    URL                 ${NAPI_URL}
    URL_MD5             ${NAPI_URL_MD5}
    CMAKE_ARGS
        -DBUILD_SHARED_LIBS=OFF    
        -DCMAKE_BUILD_TYPE="Release"          
        -DENABLE_HDF4:BOOL=ON
        -DENABLE_HDF5:BOOL=ON
        -DENABLE_MXML:BOOL=ON
        -DJPEG_INCLUDE_DIR=${NAPI_DIR}/include
        -DJPEG_LIBRARY_RELEASE=${NAPI_DIR}/lib
        -DCMAKE_MACOSX_RPATH:BOOL=OFF 
        -DHDF5_C_COMPILER_EXECUTABLE:FILEPATH=${NAPI_DIR}/bin/h5cc
        -DCMAKE_INSTALL_PREFIX=${NAPI_DIR}
    INSTALL_DIR         ${NAPI_DIR}
    INSTALL_COMMAND     ${NAPI_MAKE} install
    DEPENDS             ${HDF5_PREFIX}
	LOG_DOWNLOAD        1
	LOG_BUILD           1
)

# post-build setup
set(NAPI_INCLUDE_DIRS ${NAPI_DIR}/include)
include_directories(${NAPI_INCLUDE_DIRS})

# set the library directory variable and link it
set(NAPI_LIBRARY_DIRS ${NAPI_DIR}/lib)
link_directories(${NAPI_LIBRARY_DIRS})
set(NAPI_LIBS mxml)
set(NAPI_LIBRARY_DIRS ${NAPI_LIBRARY_DIRS})

# display info
message("Build NAPI in ${NAPI_SRC} with CMAKE")
message("NAPI_INCLUDE_DIRS=${NAPI_INCLUDE_DIRS}")
message("NAPI_LIBRARY_DIRS=${NAPI_LIBRARY_DIRS}")


message("")
