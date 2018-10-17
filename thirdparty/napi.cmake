message("****************")
message("-- BUILD NAPI --")
message("****************")

# the folder where to compile mxml-2.12
set(NAPI_PREFIX napi)

# set the source location
set(NAPI_URL ${CMAKE_CURRENT_LIST_DIR}/napi_v4.4.3.tar.gz)

# check MD5
set(NAPI_URL_MD5  51016616c0af18d14b34018da476c326)

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
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_LINKER=${CMAKE_LINKER}      
        -DCMAKE_BUILD_TYPE="Release"          
        -DENABLE_HDF4:BOOL=ON
        -DENABLE_HDF5:BOOL=ON
        -DENABLE_MXML:BOOL=ON
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
