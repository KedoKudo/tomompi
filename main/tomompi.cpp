//_____________________________________________________________________________________
//_____________________________________________________________________________________
//_____________________________________________________________________________________

#include "tomompi.h"

//_____________________________________________________________________________________

char processor_name[MPI_MAX_PROCESSOR_NAME], msg[256], *sinogram_data_set;
    
int my_id = 0, num_processes = 0, processor_name_len = 0;
                
float data_range_min, data_range_max;
                
RECON_INFO      recon_info_record;

LogFileClass    *log_file;
errorlogclass   *error_log;

//_____________________________________________________________________________________

void MPIOpen (int argc, char* argv[])
{

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &num_processes);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_id);
  MPI_Get_processor_name (processor_name, &processor_name_len);
}

//_____________________________________________________________________________________

void MPIClose (void)
{
  MPI_Finalize ();
}

//_____________________________________________________________________________________

int main (int argc, char* argv[])
{
  MPI_Comm    comm_world;
  MPI_Group   group_world, group_worker;
  int         log_index;
  char        sample_path[256], log_dir[256], check_log_dir[256], log_file_name[256], err_file_name[256],  msg[256];
  struct stat stat_buffer;

  MPIOpen (argc, argv);
            
  if (my_id == 0){
    
    strcpy (sample_path, argv[1]);
    if (sample_path[strlen (sample_path)] != '/')
      strcat (sample_path, "/");

    strcpy (log_dir, sample_path);

    //create new log directory
    log_index = 1;
    sprintf (check_log_dir, "%s%s%02d", sample_path, "logs_", log_index);
    while (stat (check_log_dir, &stat_buffer) == 0){
      
      log_index++;
      sprintf (check_log_dir, "%s%s%02d", sample_path, "logs_", log_index);
    }        
    sprintf (log_dir, "%s%s%02d/", sample_path, "logs_", log_index);
    mkdir (log_dir, 0777);    
    
    sprintf (err_file_name, "%s_%d.err", processor_name, my_id);
    error_log = new errorlogclass ();
    error_log->setErrorFileLocation (log_dir, err_file_name);

    sprintf (log_file_name, "%s_%d.log", processor_name, my_id);
    log_file = new LogFileClass (log_dir, log_file_name);
    
    StartTomoMPIServer (argc, argv, log_dir);

    delete (error_log);
    delete (log_file);
  }
  else {
    sprintf (log_file_name, "%s_%d.log", processor_name, my_id);
    sprintf (err_file_name, "%s_%d.log", processor_name, my_id);
        
    StartTomoMPIClient (argc, argv, log_file_name, err_file_name);
  }
    

  MPIClose ();

  return 0;
}

//_____________________________________________________________________________________
