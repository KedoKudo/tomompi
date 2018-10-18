# tomompi

This is a fork of the original TomoMPI code developed back in 2010.
The aim of this project is to convert the previous system dependent TomoMPI into a more self-contained system.

## Proposed plan

* Moving all neceesary source code into a single src directory
* Use cmake build system for better modular control
 1. the original code use folder as modular control
 2. the new system will relies on cmake to make it happen
* The previous system requires the system admin to install serveral dynamic libraries in a shared folder for TomoMPI to link during runtim.
 1. The new system would instead using sandbox design where third party libraries are within TomoMPI itself (waste a bit of disk space, but should improve mobility)
 2. A set of build scripts will be provided for building the thrid party libraries

## Implementation

A new branch "dev" is used to restructure the code, along with the necessary build scripts for thrid party librareis.

> The current third party software build order is 
>
> fftw  
> -> mxml  
> -> (zlib, skipped)  
> -> jpeg  
> -> szlib  
> -> HDF4  
> -> HDF5  
> -> napi

This build order is realized through a pseudo-dependices check provided by CMAKE.

## NOTE

* A constant (131072) is hard coded in teh source code (filteredBackProject.cpp). What is the reason behind this number?

* A constant varialbe, GRID_PRECISION, is declared to be 1000000. What is the reason behind this particular number?

* The reconstruction algorithm code is doing way too much work. For the next iteration of refactoring, it worthwell to modulize this part.