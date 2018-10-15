### Generated automatically from Makefile.in by configure.
all: fft recon support main link install reconsolo
reconsolo : fft recon support tomosolo link_tomosolo

##### User configurable options #####

MPI_ROOT    = /clhome/aps_tools/mpi/mpi
USR_HOME		= /clhome/B224538
SYS_HOME		= /clhome/aps_tools
USR_LOCAL_SHARED	= /clhome/aps_tools/shared

# sys env
SHELL       = /bin/sh
ARCH        = linux
CC          = $(MPI_ROOT)/bin/mpicc 
CCC         = $(MPI_ROOT)/bin/mpiCC 
CLINKER     = $(CCC)
CCLINKER    = $(CCC)
AR          = ar crl
RANLIB      = ranlib
#OPTFLAGS	= -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX

#Special flags
HDF_FLAGS = -D__unix__ -DHDF4 -DHDF5 
#-DHDF5

#Root directory for various code modules

# TOMOMPI_ROOT 	= 	$(USR_HOME)/workspace/TomoMPI
TOMOMPI_ROOT 	= 	$(USR_HOME)/tomompi
# NEXUSLIB_ROOT	=	$(USR_HOME)/workspace/NexusLibrary
NEXUSLIB_ROOT	=	$(USR_HOME)/workspace/NexusLibrary/trunk

HDF_ROOT	=	$(USR_LOCAL_SHARED)
HDF5_ROOT	=	$(USR_LOCAL_SHARED)
XML_ROOT	=	$(USR_LOCAL_SHARED)

# UTILITY_ROOT	=	$(USR_HOME)/workspace/utility
UTILITY_ROOT	=	$(USR_HOME)/workspace/Utility/trunk

NAPI_ROOT	=	$(SYS_HOME)/nexus-4.1.0

#Include directories for various code modules
TOMOMPI_INC 	= 	-I$(TOMOMPI_ROOT)/include
MPI_INC 	= 	-I$(MPI_ROOT)/include
NEXUSLIB_INC	=	-I$(NEXUSLIB_ROOT)/include
NAPI_INC	=	-I$(NAPI_ROOT)/include
HDF_INC		=	-I$(HDF_ROOT)/include
HDF5_INC	=	-I$(HDF5_ROOT)/include
XML_INC		=	-I$(XML_ROOT)
UTILITY_INC	=	-I$(UTILITY_ROOT)
FFTW_INC	=	-I$(SYS_ROOT)/fftw-3.1.2/include

#Library directories of various code modules
HDF_LIBS	=	-L$(HDF_ROOT)/lib -lmfhdf -ldf -ljpeg -lz \
			-L/usr/local/szip/szip/lib -lsz
HDF5_LIBS	=	-L$(HDF5_ROOT)/lib -lhdf5
XML_LIBS	=	-L$(XML_ROOT) -lmxml
NEXUS_LIBS	=	-L$(NEXUSLIB_ROOT)/lib -lnexuslibrary
FFTW_LIBS	= 	-L$(USR_LOCAL_SHARED)/lib -lfftw3f
PROF_LIBS   	= 	-L$(MPE_ROOT)/lib -lmpe 
LOG_LIBS    	= 	-L$(MPE_ROOT)/lib -llmpe -lmpe 
TRACE_LIBS 	= 	-L$(MPE_ROOT)/lib -ltmpe 
OTHER_LIBS 	= 	-L/usr/lib64 -lpthread -lm
#OTHER_LIBS	=	-L/usr/lib -lpthread -lm

#Additional Object references
UTILITY_OBJS	=	$(UTILITY_ROOT)/logfileclass.o \
				$(UTILITY_ROOT)/linkedlistclass.o
NAPI_OBJS		=	$(NAPI_ROOT)/src/napi.o \
				$(NAPI_ROOT)/src/napi4.o \
				$(NAPI_ROOT)/src/napi5.o \
				$(NAPI_ROOT)/src/nxdataset.o \
				$(NAPI_ROOT)/src/nxio.o \
				$(NAPI_ROOT)/src/nxxml.o \
				$(NAPI_ROOT)/src/nxstack.o \
				$(NAPI_ROOT)/src/stptok.o \
				$(NAPI_ROOT)/src/napiu.o 

### End User configurable options ###

CCFLAGS	  	= $(CFLAGS)
EXECS	  	= 
MPE_CFLAGS  = -DMPI_LINUX -DUSE_STDARG -DHAVE_PROTOTYPES
#OPTFLAGS	= $(MPE_CFLAGS) -O3 

-include fft/fft.mk
-include recon/recon.mk

-include support/support.mk
-include main/main.mk



fft: $(FFT_SRC)

recon: $(RECON_SRC)

support: $(SUPPORT_SRC)

main: $(MAIN_SRC)

link: $(MAIN_PROGS)

install: $(MAIN_INSTALL)

clean: fft_clean recon_clean support_clean main_clean

