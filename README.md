# tomompi

This is a fork of the original TomoMPI code developed back in 2010.
The aim of this project is to convert the previous system dependent TomoMPI into a more self-contained system.


## Installation

* Clone this repository

* Make the building directory and run cmake to configure the project using the commands below

```bash
make -p build
cmake ..
```

* Once the CMake config completes, run the following commands within the building directory

```bahs
make; make install
```
* Now all encessary libraries are under ${PROJECT_DIR}/lib, ${PROJECT_DIR}/lib64 (linux), and the header files are in ${PROJECT_DIR}/include.

> NOTE: you might need to install __MPICH__, __bison__ and __flex__ if these packages are not already installed on your system by your admin.

> NOTE: you might have to add ${PROJECT_DIR}/lib to your LD_LIBRARY_PATH if you cannot move these libraries to system default location safely.

## Implementation

A new branch "dev" is used to restructure the code, along with the necessary build scripts for thrid party librareis.

> The current third party software build order is 
>
> fftw  
> -> mxml  
> -> (zlib, skipped)  
> -> jpeg  
> -> szlib (only for encoding) 
> -> HDF4  
> -> HDF5  
> -> napi

This build order is realized through a pseudo-dependices check provided by CMAKE.

## NOTE

* A constant (131072) is hard coded in teh source code (filteredBackProject.cpp). What is the reason behind this number?

* A constant varialbe, GRID_PRECISION, is declared to be 1000000. What is the reason behind this particular number?

* The reconstruction algorithm code is doing way too much work. For the next iteration of refactoring, it worthwell to modulize this part.
