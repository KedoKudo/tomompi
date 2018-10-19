#ifndef MPIDEAMON_H
#define MPIDEAMON_H

//_____________________________________________________________________________________

#include <stdio.h>

#ifdef _WIN32_
#include <windows.h>
#include <sys\timeb.h>
#else
#include <sys/timeb.h>
#endif

#include <time.h>
#include <pthread.h>
#include <stdio.h>
#include <sys/stat.h>

#include "mpi_communication_structures.h"

#include "nexusbox.h"
#include "linkedlistclass.h"
#include "logfileclass.h"
#include "errorlogclass.h"
#include "mpi.h"
#include "tomocommands.h"

//_____________________________________________________________________________________
//Defines for version
#define VERSION	"V2.5"

//_____________________________________________________________________________________
//Defines for sinogram status

#define	NOT_DONE_YET	0
#define	BEING_PROCESSED	1
#define	COMPLETED		2

//_____________________________________________________________________________________
//File format

#define	HDF4	0
#define	HDF5	1
#define BIN		2

//_____________________________________________________________________________________
//Defines for reconstruction algorithm

#define         RECONSTRUCTION_FBP_DELAY_ONLY                  0
#define         RECONSTRUCTION_GRIDREC_DELAY_ONLY              1
#define         RECONSTRUCTION_FBP_NO_OPTIMIZATION             2
#define         RECONSTRUCTION_FBP_OPTIMIZED                   3
#define         RECONSTRUCTION_FBP_CYLINDER_ONLY               4
#define         RECONSTRUCTION_GRIDREC                         5

//_____________________________________________________________________________________
//Defines for zero padding options for GRIDREC

#define         GRIDREC_PADDING_NONE                           0
#define         GRIDREC_PADDING_HALF                           1
#define         GRIDREC_PADDING_ONE_AND_HALF                   2
#define         GRIDREC_PADDING_BOUNDARY                       3

//_____________________________________________________________________________________
//Defines for Debug Type

#define DEBUG_NONE                  0
#define DEBUG_WHITEFIELD            1
#define DEBUG_DARKFIELD             2
#define DEBUG_PRENORM               3
#define DEBUG_POSTNORM              4
#define DEBUG_POSTCENTERING         5
#define DEBUG_POSTRING              6
#define DEBUG_FULL                  7
#define DEBUG_NO_OUTPUT             8
#define DEBUG_NO_RECONSTRUCTION		9

//_____________________________________________________________________________________


typedef struct {
	int		sinogram_number,
			status,
			processing_client;
} SINOGRAMINFO;

typedef struct {
	float			*reconstruction;
	int				recon_number;
} RECONINFO;

typedef struct {
	int		start_sinogram_number,
			num_sinograms_to_create;
} CREATESINOINFO;

//_____________________________________________________________________________________

class FileListClass : public LinkedListClass
{
public:
	char	file_name[256];
    FileListClass (char *file_name){
		strcpy(this->file_name, file_name);
		next_in_list = NULL;
	};

	void Print (void){
		printf ("Found file %s\n", file_name);
	};

    ~FileListClass (void){};
};

//_____________________________________________________________________________________

int StartTomoMPIServer (int argc, char* argv[], char *log_file_path);
int StartTomoMPIClient (int argc, char* argv[], char *log_file_name, char *err_file_name);


void RingCorrectionSingle (float *data, float ring_coeff);
void LogProj(float *data, int xdim, int ydim);

//_____________________________________________________________________________________

extern char                         processor_name[MPI_MAX_PROCESSOR_NAME];
extern char                         msg[256];
extern char                         *sinogram_data_set;


extern int                          my_id;
extern int                          num_processes;
extern int                          processor_name_len;

extern float						data_range_min;
extern float						data_range_max;

extern RECON_INFO                   recon_info_record;

extern LogFileClass                 *log_file;
extern errorlogclass                *error_log;
//_____________________________________________________________________________________

#endif
