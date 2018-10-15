//_____________________________________________________________________________________
//_____________________________________________________________________________________
#include <fstream>
using namespace std;

#include <math.h>

#include "tomompi.h"

#include "mpi_communication_structures.h"
#include "recon_algorithm.h"
#include "filteredbackprojection.h"
#include "gridrec.h"

#include "centeringclass.h"

//_____________________________________________________________________________________

typedef struct {
  int
    sinogram_size,
    reconstruction_size,
    sinogram_adjusted_xdim,
    sinogram_adjusted_size,
    *sinogram_numbers,
    *sinogram_mpi_numbers,
    *sinogram_calc_numbers;

  unsigned short int  *short_sinogram,
                      *short_white_field_sino,
                      *int_scale_buffer;
  float		*sinograms,
          *reconstructions,
          *sino_mpi_buffer,
          *sino_calc_buffer,
          *recon_mpi_buffer,
          *recon_calc_buffer,
          *shifted_recon,
          *shifted_sinogram,    
          *dark_field_sino_ave,
          *float_white_field_sino,
          *debug_post_centering_sinogram,
          *debug_post_ring_sinogram,
          *debug_post_centering_sinogram_mpi_buffer,
          *debug_post_centering_sinogram_calc_buffer,
          *debug_post_ring_sinogram_mpi_buffer,
          *debug_post_ring_sinogram_calc_buffer;

  // buffers for boundary padding
  float *sinograms_boundary_padding;
  float *reconstructions_boundary_padding;

  float *sinogram_shifts;
  float *sinogram_mpi_shifts;
  float *sinogram_calc_shifts;

  // buffers for ring correction, borrowed from centeringclass.cpp
  float* mean_vect;
  float* low_pass_sino_lines_data;
  float*  mean_sino_line_data;

} INFO_RECORD;

//_____________________________________________________________________________________

int processing_slice_timer       = 0,
    total_processing_slice_timer = 0,
    MPI_requests_timer           = 0,
    MPI_sinogram_timer           = 0,
    single_MPI_sinogram_timer    = 0,
    write_reconstruction_timer   = 0,
    delay_time                   = 0,
    num_reconstructions_processed;

NexusBoxClass *recon_file;

INFO_RECORD information;

long int  image_size;

bool  recon_thread_running;

ReconAlgorithm      *recon_algorithm;
CenteringClass	    centering_processor;

char  log_file_name  [256],
      error_file_name[256];

//_____________________________________________________________________________________

void ClientAcknowledgements (void){
  LogFileClass *acknowledge_file;

  acknowledge_file = new LogFileClass ("./", "client_acknowledge.txt");

  acknowledge_file->Message ("__________________________________________________________________");
  acknowledge_file->Message ("TomoMPI_Client.cpp");
  acknowledge_file->Message ("");
  acknowledge_file->Message ("Controlling code for TomoMPI Client calculations.");
  acknowledge_file->Message ("Developed and maintained by:");
  acknowledge_file->Message ("       Brian Tieman");
  acknowledge_file->Message ("       Argonne National Laboratory");
  acknowledge_file->Message ("       tieman@aps.anl.gov");
  acknowledge_file->Message ("");
  acknowledge_file->Message ("8/20/2003   V1.0   BT  First version with acknowledgements");
  acknowledge_file->Message ("9/4/2003    V1.0   BT  Integer sinograms and now passed to the client.");
  acknowledge_file->Message ("        The client does the normalization step upon recieving the integer");
  acknowledge_file->Message ("        sinograms and a sinogram of the white fields and dark fields.");
  acknowledge_file->Message ("9/12/2003   V1.2   BT  Added code to shift back the reconstruction by the");
  acknowledge_file->Message ("        number of pixels it was shifted during the centering step.  This");
  acknowledge_file->Message ("        appears to fix the jumping around of the image that was occuring");
  acknowledge_file->Message ("        with some samples.");
  acknowledge_file->Message ("9/12/2003   V1.2   BT  Fixed a few minor de-allocation bugs that were causing");
  acknowledge_file->Message ("        the system to fail to shut down properly--sometimes loosing the last");
  acknowledge_file->Message ("        few reconstructions.");
  acknowledge_file->Message ("9/12/2003   V1.2   BT  Fixed a bug that would cause the system to hang if there");
  acknowledge_file->Message ("        were fewer reconstructions to do than processors to do them on.  This");
  acknowledge_file->Message ("        mostly makes debugging easier as I can now reconstruct fewer images to");
  acknowledge_file->Message ("        test the system.");
  acknowledge_file->Message ("1/30/2004   V1.3   BT  Upgraded cluster to run on RedHat 9.0.  This version fails");
  acknowledge_file->Message ("        to run properly when linked with the MKL.  A few memory issues were fixed");
  acknowledge_file->Message ("        as a result of this upgrade.");
  acknowledge_file->Message ("2/26/2004   V1.4   BT  Dark fields are now averaged before being sent to the client.");
  acknowledge_file->Message ("        This allows for more flexibility in the acquisition--for example, it's now");
  acknowledge_file->Message ("        possible to take a number of darks at the end or take darks at regular");
  acknowledge_file->Message ("        intervals throughout the dat set and the same reconstruction code can");
  acknowledge_file->Message ("        handle it transparently.");
  acknowledge_file->Message ("3/15/2004   V1.5   BT  The client now recieves the centering related fields--");
  acknowledge_file->Message ("        fixed_shift and fixed_shift_value.  If fixed_shift is non-zero, the auto");
  acknowledge_file->Message ("        findcenter routine will be ignored and the value in fixed_shift_value will");
  acknowledge_file->Message ("        be used for the shift value for all sinograms.");
  acknowledge_file->Message ("3/26/2004   V1.5   BT  The client will now write out the filter it is using to the");
  acknowledge_file->Message ("        log file.");
  acknowledge_file->Message ("5/20/2004   V1.5   BT  There is now support for intermediate calculation debug info.");
  acknowledge_file->Message ("        There are several points in the calculation chain where sinograms may be");
  acknowledge_file->Message ("        written out to disk.  The debug type is specified in the recdonstruction.cfg");
  acknowledge_file->Message ("        file by dsetting the <Debug> field to one of the following types:");
  acknowledge_file->Message ("        None -- no debuging");
  acknowledge_file->Message ("        SinoPreNorm -- sinograms before normalization");
  acknowledge_file->Message ("        SinoPostNorm -- sinograms post normlalization");
  acknowledge_file->Message ("        PostCentering -- sinograms after centering");
  acknowledge_file->Message ("        PostRingRemoval -- sinograms post ring removal");
  acknowledge_file->Message ("6/1/2004    V1.5  BT  Added a new debug mode");
  acknowledge_file->Message ("        MarkSuspicious -- look for suspicious reconstructions and note them in the");
  acknowledge_file->Message ("        log.  Also, try and output the whitefield, darkfield, raw sinogram, normalized");
  acknowledge_file->Message ("        sinogram, post centering sinogram, and post ring removal sinogram for that slice.");
  acknowledge_file->Message ("        This mode chews up a lot of memory so it might be better to run just a single");
  acknowledge_file->Message ("        client on each node.");
  acknowledge_file->Message ("6/7/2004    V1.6  BT  Added another new debug mode");
  acknowledge_file->Message ("        I broke down and did the accounting that needed to be done to write all intermediate");
  acknowledge_file->Message ("        results from outside the reconstruction thread.  This means that all intermediate");
  acknowledge_file->Message ("        results may now be written at once while avoiding thread issues with HDF not being");
  acknowledge_file->Message ("        thread safe.  The new debug mode is:");
  acknowledge_file->Message ("        Full -- write all intermediate results.  Note, this mode creates additional client");
  acknowledge_file->Message ("        memory and uses a ton of disk space!");
  acknowledge_file->Message ("6/22/2004   V1.7  BT  Finally!  I finally found the bug that was causing random slices to ");
  acknowledge_file->Message ("        be reconstructed incorrectly!  The problem turned out to be in the interplay of the");
  acknowledge_file->Message ("        ring removal routine and the normalize routine.  The normalize routine had an");
  acknowledge_file->Message ("        improperly structured while loop--while (frame_number <= recon_info_record.sinogram_ydim)--");
  acknowledge_file->Message ("        instead of--while (frame_number < recon_info_record.sinogram_ydim).  The = sign was causing");
  acknowledge_file->Message ("        the loop to execute 1 extra time and thus it overwrote part of the buffer that");
  acknowledge_file->Message ("        happened to be used by the ring removal routine.  As the two routines are in different");
  acknowledge_file->Message ("        threads, the result was a random, unpredictable, poorly reconstructed slice.");
  acknowledge_file->Message ("12/10/2004	V1.8  BT Built tomompi on development cluster (AMD opteron machines).  Turned on");
  acknowledge_file->Message ("		compiler optimizations.  ALso built with MPE profiling libraries so I could start");
  acknowledge_file->Message ("		trying to use jumpshot to profile our performance.");
  acknowledge_file->Message ("12/10/2004	V1.8  BT Created a recon info record structure that contains all the initialization");
  acknowledge_file->Message ("		parameters needed by the clients.  This allowed for a significant reduction of MPI");
  acknowledge_file->Message ("		sends to the clients upon initialization as nearly the entire structure can be sent");
  acknowledge_file->Message ("		at once--as opposed to each datum being sent seperately.  Also did this for the sinogram");
  acknowledge_file->Message ("		data that was being sent each time a sinogram was requested.  This appears to have about");
  acknowledge_file->Message ("		a 10% impact on performance.");
  acknowledge_file->Message ("8/9/2005	V1.8  BT Modified nromalize routine to take into account sinograms that are not a");
  acknowledge_file->Message ("		power of 2 in size (as required by gridrec).  This was done by introducing an");
  acknowledge_file->Message ("		sinogram_adjusted_xdim and sinogram_adjusted_size to the information structure");
  acknowledge_file->Message ("		and subsequently using that information to create appropriate space for the larger");
  acknowledge_file->Message ("		data structures.  The raw sinogram is then centered in the larger power of 2 sinogram");
  acknowledge_file->Message ("		with the dead space before the raw sinogram being filled with the results for the first");
  acknowledge_file->Message ("		pixel and the dead space after the raw sinogram being filled with the results from the");
  acknowledge_file->Message ("		last pixel in the raw sinogram.  This is handled in the normalization routine while the");
  acknowledge_file->Message ("		sinograms are being  normalized.  It's possible that in the future we will also crop");
  acknowledge_file->Message ("		the resulting reconstruction to remove the dead space and save disk space/write time.");
  acknowledge_file->Message ("1/12/2006	V1.9  BT A problem was found where the sinogram_ydim and theta_list_size values were one");
  acknowledge_file->Message ("		bigger than correct.  This was due to a miscount of the files.  To figure out the sinogram");
  acknowledge_file->Message ("		size the cluster reads the exp file and counts the projections.  But the projections were");
  acknowledge_file->Message ("		being counted in a while loop where the counter was being incremented at the end of the loop.");
  acknowledge_file->Message ("		It wasn't until count = projections+1 that the while condition failed and exitied then");
  acknowledge_file->Message ("		sinogram_ydim was set to counter.  This then was propogated throughout the system.");
  acknowledge_file->Message ("1/24/2006	V1.9  BT Added a reconstruction.cfg flag <Ring Removal> which can be 0 or 1.  1 means");
  acknowledge_file->Message ("		use ring removal.  0 means do not use ring removal.");
  acknowledge_file->Message ("1/31/2006	V1.10 BT Added a reconstruction.cfg parameter <Ring Removal Coeff> which can be 0= < coeff <= 5.");
  acknowledge_file->Message ("		This adjusts a coefficient used in the ring removal code.  See the CenteringCalss comments for");
  acknowledge_file->Message ("		more specific information about what this does.");
  acknowledge_file->Message ("2/13/2006	V1.10 BT Fixed a problem with the output of the Post Center and Post Ring buffers that");
  acknowledge_file->Message ("		crept in after the code was modified to padd the data to a multiple of 2 wide.  Basically,");
  acknowledge_file->Message ("		the Post Center and Post Ring buffers were being written to disk with an incorrect width dimension.");
  acknowledge_file->Message ("3/23/2006	V1.11 BT Added a methode to compress reconstructed files.  Just using the built in loss-less");
  acknowledge_file->Message ("		compression in HDF provides poor compression results on the floating point reconstructed values.");
  acknowledge_file->Message ("		However, as the reconstruction only occupies a circle in the image, it's possible to set the corners");
  acknowledge_file->Message ("		to 0.0 without loosing any data.  All these 0.0 values then compress rather well yielding an");
  acknowledge_file->Message ("		approximate 25% reduction in file size for the saved reconstruction files.");
  acknowledge_file->Message ("		The appropriate reconstruction.cfg flag is <Copression Type> and its values may be:");
  acknowledge_file->Message ("			NONE -- no compression (corners will not even be set to 0.0");
  acknowledge_file->Message ("			LZW -- gzip type compression (this appears to be the best for most data sets");
  acknowledge_file->Message ("			RLE -- Run Length Encoding");
  acknowledge_file->Message ("			NONE -- Skipping-Huffman");
  acknowledge_file->Message ("3/23/2006	V1.11 BT While debugging the problem noted in the centeringclass.cpp on 3/22/2006 it");
  acknowledge_file->Message ("		was discovered that averaging the white fields can sometimes give a far superior result.");
  acknowledge_file->Message ("		A flag was added that allows the user to switch from using a single white field per projection");
  acknowledge_file->Message ("		to averaging all the white fields and using the average white field for all projections.  The");
  acknowledge_file->Message ("		flag is <Average White Fields> and it should be set to 1 to average all the whites or to 0");
  acknowledge_file->Message ("		to use a single white field (last acquired white field before acquisition of projection)");
  acknowledge_file->Message ("		per projection");
  acknowledge_file->Message ("4/23/2006	V1.13 BT Added a new debug mode:");
  acknowledge_file->Message ("		No Reconstruction -- perform all the steps of a reconstruction, but, to save time, do not");
  acknowledge_file->Message ("		actually reconstruct the data.");
  acknowledge_file->Message ("4/23/2006	V1.13 BT Removed MarkSuspicious from the possible debug modes--it wasn't useful.");
  acknowledge_file->Message ("4/23/2006   V1.13 BT Modified handling of filter type to be done by int rather than by string.");
  acknowledge_file->Message ("8/1/2006	V1.14 BT Modified to use FFTW as the fft algorithm.  FFTW does an algorithm performance");
  acknowledge_file->Message ("		comparison the first time it is called to determine what the fastest algorithm to use is.");
  acknowledge_file->Message ("		It then uses that algorithm for subsequent fft calls--provided the input data is the same");
  acknowledge_file->Message ("		size.  FFTW appears to be about 3x faster than Numerical Recipies for 2048 sinograms.");
  acknowledge_file->Message ("10/20/2006  V1.15 BT Fixed a few array indexing errors in the normalization routine.  These errors");
  acknowledge_file->Message ("		were causing normalization to be done improperly by not indexing the white fields as they should.");
  acknowledge_file->Message ("		This resulted in reconstructions where the ring removal routine failed to work as well as it");
  acknowledge_file->Message ("		should.  There is still one last problem we know of with the normalization:  the first band");
  acknowledge_file->Message ("		of sinogram lines are still not being normalized properly--they appear brighter than the rest");
  acknowledge_file->Message ("		of the sinogam.  It is not known where this error comes from or if it has a serious impact on");
  acknowledge_file->Message ("		the reconstructions.");
  acknowledge_file->Message ("1/10/2007   V2.0 BT Change of pace--a new major version number!  This version is in subversion control");
  acknowledge_file->Message ("		at svn+ssh://m81-cluster.aps.anl.gov/home/svn/repository/cpp.  The other major change is the");
  acknowledge_file->Message ("		upgrade to Nexus V3.0.0.  This brings with it the ability to write HDF5 files.  HDF5 will be the");
  acknowledge_file->Message ("		default file format and will have the extension .h5--hdf4 will retain the .hdf extension.  It");
  acknowledge_file->Message ("		is possible to use the override file to force HDF4.  When saving in HDF5, the only comrpession");
  acknowledge_file->Message ("		format available is NX_COMP_LZW--all other compression formats are not supported.  As an aside");
  acknowledge_file->Message ("		this is the first version built completely under Eclipse.");
  acknowledge_file->Message ("4/12/2007   V2.1 BT Added calculations to the clients to compute the data min and max for the entire");
  acknowledge_file->Message ("		reconstruction.  The clients track the min/max for all the slices they calculate then report");
  acknowledge_file->Message ("		their local min/max to the server who then computes the full reconstruction min/max.  The");
  acknowledge_file->Message ("		results are then stored in the exp file at ;experiment;reconstruction;data_range_min and");
  acknowledge_file->Message ("		;experiment;reconstruction;data_range_max.");
  acknowledge_file->Message ("4/12/2007	V2.1 BT Added ability to scale data into short ints by using a user supplied data_range_min");
  acknowledge_file->Message ("		and data_range_max values.  These values must be supplied prior to reconstructing as the full");
  acknowledge_file->Message ("		data range of the reconstruction can not be calculated until completion.  The data_range_min/");
  acknowledge_file->Message ("		data_range_max values of a previous reconstruction are generally good values to use.  The actual");
  acknowledge_file->Message ("		data_range_min/max values are stored in the exp file as ;experiment;reconstruction;scaled_data_min");
  acknowledge_file->Message ("		and ;experiment;reconstruction;scaled_data_max for reference.");
  acknowledge_file->Message ("11/14/2007	V2.3 BT Added a BIN file type to the output file types.  If a user selects BIN as the");
  acknowledge_file->Message ("		<output file type> the output files will be written as binary without any sort of header.");
  acknowledge_file->Message ("		Eventually, we will need a bit of a header, but we haven't settled on one yet.  Bin files");
  acknowledge_file->Message ("		can not be compressed");
  acknowledge_file->Message ("1/28/2008   V2.3 BT Modified code to not read a reconstruction.cfg file but to get that information");
  acknowledge_file->Message ("        from the command line.  This will allow for the program to be queued in a system like N1");
  acknowledge_file->Message ("        Grid Engine.  As queued applications can be reprioritized by the system, having a single");
  acknowledge_file->Message ("        reconstruction.cfg file is problematic as it's impossible to keep the file in sync with");
  acknowledge_file->Message ("        the correct application in the queue--now each instance can get that info from the");
  acknowledge_file->Message ("        command prompt.");
  acknowledge_file->Message ("1/28/2008   V2.3 BT Modified the location of the log files to be in the smaple tree.  Also allow");
  acknowledge_file->Message ("        for versioning of the log and override files.  This is done by creating a new log directory");
  acknowledge_file->Message ("        \"logs_##\" each time the sample is reconstructed.  This way, we can keep a running list");
  acknowledge_file->Message ("        of logs and configurations in case the sample is reconstructed multiple times.");
  acknowledge_file->Message ("2/4/2008    V2.3 BT Added an error log file that will only be created if there is an error.  The");
  acknowledge_file->Message ("        error log will log all errros encountered by the system.  This provides a convenient location");
  acknowledge_file->Message ("        to look for problems if reconstructions fail or do not look quite right.");
  acknowledge_file->Message ("2/4/2008    V2.3 BT The system should no longer crash due to a raw data file missing.  Instead, the");
  acknowledge_file->Message ("        system will mark in the error log that the file is missing and will grab the previous file");
  acknowledge_file->Message ("        of the same type to fill in.  If a projection is missing, the previous projection will");
  acknowledge_file->Message ("        be used to complete the reconstruction.");
  acknowledge_file->Message ("2/11/2008  V2.3 BT Added a small header to the binary output files.  The header is xdim<32bits>,");
  acknowledge_file->Message ("       ydim<32bits>, byte_depth<32bits>.  This is the bare minimum information needed to understand");
  acknowledge_file->Message ("       how to open the file if someone forgets or can not look it up elsewhere.");
  acknowledge_file->Message ("7/29/2009  V2.4 BT Added two new parameters to the cluster override_exp_file.config.  <start_fixed_shift");
  acknowledge_file->Message ("	   is the start value of a range of shifts.  <end_shift_value> is the end value of the range.");
  acknowledge_file->Message ("	   If <use_slices_file> is true, the slices in the slice file will be reconstructed at each of the");
  acknowledge_file->Message ("	   shifts in the range start-end and save in the reconstructed directory with a unique file name");
  acknowledge_file->Message ("	   containing the shift value.");

  acknowledge_file->Message ("6/30/2011 V2.5 Yongsheng Pan and BT Added the zero padding functionality to the GridRec algorithm. An entry named");
  acknowledge_file->Message ("          <GRIDREC Zero Padding> is added to the overwrite_exp_file.config, with GRIDREC_PADDING_NONE");
  acknowledge_file->Message ("          GRIDREC_PADDING_HALF(default) and GRIDREC_PADDING_ONE_AND_HALF as options. This functionality");
  acknowledge_file->Message ("          removes artifacts from GridRec.");

  acknowledge_file->Message ("10/30/2011 V2.5 Yongsheng Pan updated the padding functionality using boundary values from each sinogram.");
  acknowledge_file->Message ("          This functionality help removes artifacts caused by zero padding.");

  acknowledge_file->Message ("");
  acknowledge_file->Message ("");
  acknowledge_file->Message ("__________________________________________________________________");

  NexusBoxClass::acknowledgements (acknowledge_file);
  ReconAlgorithm::acknowledgements (acknowledge_file);
  FBP::acknowledgements (acknowledge_file);
  OptimizedFBP::acknowledgements (acknowledge_file);
  CircleFBP::acknowledgements (acknowledge_file);
  GridRec::acknowledgements (acknowledge_file);
  CenteringClass::acknowledgements (acknowledge_file);

  delete (acknowledge_file);
}

//_____________________________________________________________________________________

void FirstContact (void){
  MPI_Status		status;

  char			msg[256];

  int       rank,
            type,
            dims[2],
            temp;

  int				power,
            size;

  bool			still_smaller;

  MPI_Recv (&recon_info_record, 
            sizeof (recon_info_record), 
            MPI_BYTE, 
            MPI_ANY_SOURCE, 
            0, 
            MPI_COMM_WORLD, 
            &status);

  log_file = new LogFileClass (recon_info_record.log_path, log_file_name);

  error_log = new errorlogclass ();
  error_log->setErrorFileLocation (recon_info_record.log_path, error_file_name);

  sprintf (msg, "I am processor %s with id %d", processor_name, my_id);
  log_file->TimeStamp (msg);
  log_file->Message ("I am but a humble servant...");

  information.sinogram_size = recon_info_record.sinogram_xdim*recon_info_record.sinogram_ydim;
  information.short_sinogram = (unsigned short int *) malloc (sizeof(unsigned short int)*information.sinogram_size);
  if (information.short_sinogram == NULL){
    sprintf (msg, "Could not alocate memory for information.short_sinogram.");
    error_log->addError (msg, "FirstContact ()");
  }

  information.short_white_field_sino = (unsigned short int *) malloc (sizeof(unsigned short int)*recon_info_record.white_size);
  if (information.short_white_field_sino == NULL){
    sprintf (msg, "Could not alocate memory for information.short_white_field_sino.");
    error_log->addError (msg, "FirstContact ()");
  }

  information.float_white_field_sino = (float *) malloc (sizeof (float)*recon_info_record.sinogram_xdim);
  if (information.float_white_field_sino == NULL){
    sprintf (msg, "Could not alocate memory for information.float_white_field_sino.");
    error_log->addError (msg, "FirstContact ()");
  }

  information.dark_field_sino_ave = (float *) malloc (sizeof(float)*recon_info_record.dark_size);
  if (information.dark_field_sino_ave == NULL){
    sprintf (msg, "Could not alocate memory for information.dark_field_sino_ave.");
    error_log->addError (msg, "FirstContact ()");
  }

  sinogram_data_set = (char *) malloc (recon_info_record.sinogram_set_size);
  if (sinogram_data_set == NULL){
    sprintf (msg, "Could not alocate memory for sinogram_data_set.");
    error_log->addError (msg, "FirstContact ()");
  }

  recon_info_record.theta_list = (float *) malloc (sizeof(float)*recon_info_record.theta_list_size);
  if (recon_info_record.theta_list == NULL){
    sprintf (msg, "Could not alocate memory for recon_info_record.theta_list.");
    error_log->addError (msg, "FirstContact ()");
  }

  MPI_Recv (recon_info_record.theta_list, 
            recon_info_record.theta_list_size, 
            MPI_FLOAT, 
            recon_info_record.mayor_id, 
            0, 
            MPI_COMM_WORLD, 
            &status);

  sprintf (msg, "First contact recieved!  The server is process %d", recon_info_record.mayor_id);
  log_file->Message (msg);

  switch (recon_info_record.debug){
    case 0  : sprintf (msg, "Debug Type: None");               break;
    case 1  : sprintf (msg, "Debug Type: White Field");        break;
    case 2  : sprintf (msg, "Debug Type: Dark Field");         break;
    case 3  : sprintf (msg, "Debug Type: Pre-Normalization");  break;
    case 4  : sprintf (msg, "Debug Type: Post-Normalization"); break;
    case 5  : sprintf (msg, "Debug Type: Post-Centering");     break;
    case 6  : sprintf (msg, "Debug Type: Post-Ring Removal");  break;
    case 7  : sprintf (msg, "Debug Type: Mark Suspicious");    break;
    case 8  : sprintf (msg, "Debug Type: Full");               break;
    case 9  : sprintf (msg, "Debug Type: No Output");          break;
    case 10 : sprintf (msg, "Debug Type: No Reconstruction");  break;
    default : sprintf (msg, "Debug Type: Unknown");            break;
  }
  log_file->Message (msg);

  switch (recon_info_record.filter){
    case FILTER_NONE        : log_file->Message ("Filter: None");
    case FILTER_SHEPP_LOGAN : log_file->Message ("Filter: Shepp/Logan");
    case FILTER_HANN        : log_file->Message ("Filter: Hann");
    case FILTER_HAMMING     : log_file->Message ("Filter: Hamming");
    case FILTER_RAMP        : log_file->Message ("Filter: Ramp");
    case FILTER_FBP         : log_file->Message ("Filter: FBP");
    default                 : log_file->Message ("Filter: unknown");
  }

  sprintf (msg, "%s%d", "Theta list size: ", recon_info_record.theta_list_size);
  log_file->Message(msg);

  sprintf (msg, "Sinogram dimensions: %d x %d", recon_info_record.sinogram_xdim, recon_info_record.sinogram_ydim);
  log_file->Message (msg);

  still_smaller = true;
  power = 0;
  while (still_smaller)
    if (recon_info_record.sinogram_xdim > pow (2, power)){
	    power++;
	    still_smaller = true;
    }
    else
      still_smaller = false;

  if (recon_info_record.sinogram_xdim == pow (2, power)){
      // it seems that reconstruction_ydim and reconstruction_xdim are assumed to be equal here. Y. Pan. 6/20/2011
    information.sinogram_adjusted_xdim = recon_info_record.sinogram_xdim;

    information.sinogram_adjusted_size = information.sinogram_adjusted_xdim * recon_info_record.sinogram_ydim;
    information.reconstruction_size = recon_info_record.reconstruction_xdim*recon_info_record.reconstruction_ydim;

    sprintf (msg, "Sinograms are a power of 2!");
    log_file->Message (msg);
  }
  else{
    size = (int) pow (2, power);
    information.sinogram_adjusted_xdim = size;

    information.sinogram_adjusted_size = information.sinogram_adjusted_xdim * recon_info_record.sinogram_ydim;
    recon_info_record.reconstruction_xdim = size;
    recon_info_record.reconstruction_ydim = size;
    information.reconstruction_size = recon_info_record.reconstruction_xdim*recon_info_record.reconstruction_ydim;

    sprintf (msg, "Sinograms are not a power of 2.  They will be increased to %d", information.sinogram_adjusted_xdim);
    log_file->Message (msg);
  }

  sprintf (msg, "White field size: %d dark_field_size: %d", recon_info_record.white_size, recon_info_record.dark_size);
  log_file->Message (msg);

  sprintf (msg, "White/Dark interval: %d", recon_info_record.whitedark_interval);
  log_file->Message (msg);

  sprintf (msg, "Reconstruction dimensions are %d x %d", recon_info_record.reconstruction_xdim, recon_info_record.reconstruction_ydim);
  log_file->Message (msg);

  if (recon_info_record.centering == 1)
    if (recon_info_record.use_fixed_shift == 1)
      if (recon_info_record.use_slices_file == 1){
        sprintf (msg, "I will be centering with fixed shifts from %f -> %f.", recon_info_record.start_fixed_shift, recon_info_record.end_fixed_shift);
    	  log_file->Message(msg);
      }
      else{
	      sprintf (msg, "I will be centering with a fixed shift of %f.", recon_info_record.fixed_shift_value);
        log_file->Message(msg);
	    }
    else
      log_file->Message("I will be attempting to auto-center.");
  else
    log_file->Message("I will not be auto-centering.");

  switch (recon_info_record.file_format){
    case HDF4 : {
      log_file->Message ("File format: HDF4");
      recon_file->SetFileMode (HDF4_MODE);
      break;
    }
    case HDF5 : {
      log_file->Message ("File format: HDF5");
      recon_file->SetFileMode (HDF5_MODE);
      break;
    }
    case BIN : {
      log_file->Message ("File format: BIN");
      break;
    }
    default : {
      log_file->Message ("File format: UNKNOWN--defaulting to HDF5");
      recon_file->SetFileMode (HDF5_MODE);
      break;
    }
  }

  if (recon_info_record.use_ring_removal){
    sprintf (msg, "Ring removal is: ON");
    log_file->Message(msg);
    sprintf (msg, "Ring removal coefficient is: %f", recon_info_record.ring_removal_coeff);
    log_file->Message(msg);
  }
  else{
    sprintf (msg, "Ring removal is: OFF");
    log_file->Message(msg);
  }

  if (recon_info_record.average_white_fields)
    log_file->Message("Averaging white fields: ON");
  else
    log_file->Message("Averaging white fields: OFF");

  switch (recon_info_record.compression_type){
    case NX_COMP_NONE : log_file->Message ("Compression type: NONE");                            break;
    case NX_COMP_LZW  : log_file->Message ("Compression type: LZW");                             break;
    case NX_COMP_RLE  : log_file->Message ("Compression type: RLE");                             break;
    case NX_COMP_HUF  : log_file->Message ("Compression type: HUF");                             break;
    default           : log_file->Message ("Compression type: unknown--this may be a problem!"); break;
    }
  recon_file->CompressionScheme (recon_info_record.compression_type);

  //if we're rescaling the data, we need an int buffer to scale to...
  if (recon_info_record.rescale_to_int){
    information.int_scale_buffer = (unsigned short int *) malloc (sizeof (unsigned short int) * recon_info_record.reconstruction_ydim * recon_info_record.reconstruction_ydim);
    if (information.int_scale_buffer == NULL){
      sprintf (msg, "Could not alocate memory for int_scale_buffer.");
      error_log->addError (msg, "FirstContact ()");
    }
  }

  if (recon_info_record.rescale_to_int)
    type = NX_UINT16;
  else
    type = NX_FLOAT32;
  
  rank = 2;
  dims[0] = recon_info_record.reconstruction_ydim;
  dims[1] = recon_info_record.reconstruction_xdim;
  recon_file->RegisterVar ("Image", rank, dims, type, information.reconstructions);
  recon_file->InitTemplate ("./", "output.template");

  data_range_min = 1000000.0;
  data_range_max = -1000000.0;

  log_file->TimeStamp ("Now I know everything");

}

//_____________________________________________________________________________________

void WriteBinaryData (void *data, int slice, char *prefix, int xdim, int ydim, int type){
  char    		file_name[256];
  ofstream		output_file;
  int				datum_size;

  sprintf (file_name, "%s%s_%s_%05d.bin", recon_info_record.reconstruction_path, prefix, recon_info_record.base_name, slice);

  output_file.open (file_name);

  switch (type){
    case NX_INT8  :;
    case NX_UINT8 : datum_size = 1; break;

    case NX_INT16  :;
    case NX_UINT16 : datum_size = 2; break;

    case NX_INT32  :;
    case NX_UINT32 : datum_size = 4; break;

    case NX_FLOAT32 : datum_size = 4; break;
    }

  output_file.write ((char *) &xdim, sizeof (xdim));
  output_file.write ((char *) &ydim, sizeof (ydim));
  output_file.write ((char *) &datum_size, sizeof (datum_size));

  output_file.write ((char *) data, datum_size*xdim*ydim);

  output_file.close ();

}

//_____________________________________________________________________________________

void WriteHDFData (void *data, int slice, char *prefix, int xdim, int ydim, int type){
  char    file_name[256];
  int     dims[2];
  float   *temp;

  dims[0] = ydim;
  dims[1] = xdim;
  recon_file->UpdateVarInfo ("Image", dims, type, data);

  recon_file->UpdateVars ();

  strcpy (file_name, prefix);
  strcat (file_name, "_");
  strcat (file_name, recon_info_record.base_name);
  switch (recon_info_record.file_format){
    case HDF4 : sprintf (file_name, "%s_%05d.hdf", file_name, slice); break;
    case HDF5 :	sprintf (file_name, "%s_%05d.h5", file_name, slice);  break;
  }

  recon_file->WriteAll (recon_info_record.reconstruction_path, file_name);
}

//_____________________________________________________________________________________

void WriteData (void *data, int slice, char *prefix, int xdim, int ydim, int type){
  if (recon_info_record.file_format == BIN)
    WriteBinaryData (data, slice, prefix, xdim, ydim, type);
  else
    WriteHDFData (data, slice, prefix, xdim, ydim, type);
}

//_____________________________________________________________________________________
void WriteBinaryReconstruction (void *reconstruction, int slice, float shift){
  char    		file_name[256];
  ofstream		output_file;
  int				datum_size;

  if (!recon_info_record.use_slices_file){
    sprintf (file_name, "%srec_%s_%05d.bin", recon_info_record.reconstruction_path, recon_info_record.base_name, slice);
  }
  else{
    sprintf (file_name, "%srec_%s_%05d_%f.bin", recon_info_record.reconstruction_path, recon_info_record.base_name, slice, shift);
  }

  output_file.open (file_name);

  output_file.write ((char *) reconstruction, sizeof(float)*recon_info_record.reconstruction_xdim*recon_info_record.reconstruction_ydim);
  output_file.flush ();

  output_file.close ();

}

//_____________________________________________________________________________________
void WriteHDFReconstruction (void *reconstruction, int slice, float shift){
  char    file_name[256];
  int     dims[2];
  float   *temp,
          source_range,
          scale_factor;
  int		r,
        destination_range;

  temp = (float *) reconstruction;

  if (recon_info_record.compression_type != NX_COMP_NONE){
    r = recon_info_record.reconstruction_xdim/2;
    sprintf (msg, "Radius: %d", r);
    log_file->Message (msg);

    int count = 0;
    for (int loop=-r;loop<r;loop++)
      for (int loop2=-r;loop2<r;loop2++)
        if ((loop*loop+loop2*loop2) > r*r){
          count++;
          temp[(loop+r)*recon_info_record.reconstruction_xdim+(loop2+r)] = 0.0;
        }
  }

  if (recon_info_record.rescale_to_int){
    destination_range = 65535;
    source_range = recon_info_record.scale_data_range_max - recon_info_record.scale_data_range_min;

    sprintf (msg, "dest range: %d, source range: %f", destination_range, source_range);
    log_file->Message (msg);

    for (int loop=0;loop<recon_info_record.reconstruction_ydim;loop++)
      for (int loop2=0;loop2<recon_info_record.reconstruction_xdim;loop2++){

        if (temp[loop*recon_info_record.reconstruction_xdim+loop2] > recon_info_record.scale_data_range_max)
          temp[loop*recon_info_record.reconstruction_xdim+loop2] = recon_info_record.scale_data_range_max;

        if (temp[loop*recon_info_record.reconstruction_xdim+loop2] < recon_info_record.scale_data_range_min)
          temp[loop*recon_info_record.reconstruction_xdim+loop2] = recon_info_record.scale_data_range_min;

        scale_factor = (temp[loop*recon_info_record.reconstruction_xdim+loop2] - data_range_min) / source_range;
        information.int_scale_buffer[loop*recon_info_record.reconstruction_xdim+loop2] = (unsigned short int) floor (scale_factor * destination_range);
      }

    dims[0] = recon_info_record.reconstruction_ydim;
    dims[1] = recon_info_record.reconstruction_xdim;
    recon_file->UpdateVarInfo ("Image", dims, information.int_scale_buffer);
  }
  else{
    dims[0] = recon_info_record.reconstruction_ydim;
    dims[1] = recon_info_record.reconstruction_xdim;
    recon_file->UpdateVarInfo ("Image", dims, reconstruction);
  }

  recon_file->UpdateVars ();

  if (!recon_info_record.use_slices_file)
    switch (recon_info_record.file_format){
      case HDF4 : sprintf (file_name, "rec_%s_%05d.hdf", recon_info_record.base_name, slice); break;
      case HDF5 :	sprintf (file_name, "rec_%s_%05d.h5", recon_info_record.base_name, slice);  break;
    }
  else
    switch (recon_info_record.file_format){
      case HDF4 : sprintf (file_name, "rec_%s_%05d_%f.hdf", recon_info_record.base_name, slice, shift); break;
      case HDF5 :	sprintf (file_name, "rec_%s_%05d_%f.h5", recon_info_record.base_name, slice, shift);  break;
    }

  recon_file->WriteAll (recon_info_record.reconstruction_path, file_name);

}

//_____________________________________________________________________________________

void WriteReconstruction (void *reconstruction, int slice, float shift){
  sprintf (msg, "WriteReconstruction--slice: %d, shift: %f", slice, shift);

  log_file->Message (msg);

  log_file->StartTimer (write_reconstruction_timer);

  if (recon_info_record.file_format == BIN)
    WriteBinaryReconstruction (reconstruction, slice, shift);
  else
    WriteHDFReconstruction (reconstruction, slice, shift);

  log_file->StopTimer (write_reconstruction_timer);
  log_file->AccumulateTimer (write_reconstruction_timer);

  num_reconstructions_processed++;

}

//_____________________________________________________________________________________

void Normalize (unsigned short int *short_sino, 
                unsigned short int *short_white_field_sino, 
                float *dark_field_sino_ave, 
                float *norm_sino){
  int frame_interval,
      frame_number,
      white_dark_number,
      loop,
      loop2,
      pad_size,
      front_pad_size,
      back_pad_size;

  float temp,
       *float_white_field_sino;

  float_white_field_sino = information.float_white_field_sino;

  pad_size = information.sinogram_adjusted_xdim - recon_info_record.sinogram_xdim;
  front_pad_size = pad_size / 2;
  back_pad_size = pad_size - front_pad_size;

  if (recon_info_record.average_white_fields){
    for (loop=0;loop<recon_info_record.sinogram_xdim;loop++)
      float_white_field_sino[loop] = 0;

    for (loop=0;loop<recon_info_record.num_white_fields;loop++)
      for (loop2=0;loop2<recon_info_record.sinogram_xdim;loop2++)
        float_white_field_sino[loop2] = float_white_field_sino[loop2] + (((float) short_white_field_sino[(loop*recon_info_record.sinogram_xdim)+loop2]) / recon_info_record.num_white_fields);
  }
  else
    for (loop=0;loop<recon_info_record.sinogram_xdim;loop++)
      float_white_field_sino[loop] = (float) short_white_field_sino[loop];

  frame_interval = 0;
  frame_number = 0;
  white_dark_number = 0;
  while (frame_number < recon_info_record.sinogram_ydim){
    if (!recon_info_record.average_white_fields)
	    if (frame_interval == recon_info_record.whitedark_interval){
        frame_interval = 0;
        white_dark_number++;
        for (loop=0;loop<recon_info_record.sinogram_xdim;loop++)
          float_white_field_sino[loop] = (float) short_white_field_sino[white_dark_number*recon_info_record.sinogram_xdim+loop];
	    }

    //if sinogram is smaller than a power of 2, pad front with results of first real sinogram pixel
    temp = float_white_field_sino[0] - dark_field_sino_ave[0];
    if (temp != 0)
	    temp = ((float) short_sino[frame_number*recon_info_record.sinogram_xdim] - dark_field_sino_ave[0]) / temp;
    else
	    temp = 1e-3;

    for (loop=0;loop<front_pad_size;loop++)
	    norm_sino[frame_number*information.sinogram_adjusted_xdim+loop] = temp;

    for (loop=front_pad_size;loop<front_pad_size+recon_info_record.sinogram_xdim;loop++){
	    temp = float_white_field_sino[(loop-front_pad_size)] - dark_field_sino_ave[(loop-front_pad_size)];

	  if (temp != 0)
	    norm_sino[frame_number*information.sinogram_adjusted_xdim+loop] = ((float) short_sino[frame_number*recon_info_record.sinogram_xdim+(loop-front_pad_size)] - dark_field_sino_ave[(loop-front_pad_size)]) / temp;
	  else
	    norm_sino[frame_number*information.sinogram_adjusted_xdim+loop] = 1e-3;
	  }

    //if sinogram is smaller than a power of 2, pad back with results of last real sinogram pixel
    temp = float_white_field_sino[(recon_info_record.sinogram_xdim-1)] - dark_field_sino_ave[recon_info_record.sinogram_xdim-1];
    if (temp != 0)
      temp = ((float) short_sino[frame_number*recon_info_record.sinogram_xdim+(recon_info_record.sinogram_xdim-1)] - dark_field_sino_ave[recon_info_record.sinogram_xdim-1]) / temp;
    else
      temp = 1e-3;
    for (loop=front_pad_size+recon_info_record.sinogram_xdim;loop<information.sinogram_adjusted_xdim;loop++)
      norm_sino[frame_number*information.sinogram_adjusted_xdim+loop] = temp;

    frame_interval++;
    frame_number++;
  }

}

//---------------------------------------------------------------------------
// from CenteringClass::LogProj
void LogProj(float *data, int xdim, int ydim) { 
  int     i, k; 
  float   mean, max; 
  
  for (i=0;i<ydim;i++) {
    max = data[i*xdim]; 
    for (k=0;k<xdim;k++) {
      if (data[i*xdim+k] > max) 
	      max = data[i*xdim+k]; 
    }

    for (k=0;k<xdim;k++) { 
      if (data[i*xdim+k] <= 0.0) 
	      data[i*xdim+k] = 1.0; // this is really only to check if is == 0 
        data[i*xdim+k] = log (max/data[i*xdim+k]); 
    } 
  } 
} 
 

void LogSinogram (float *data, int xdim, int ydim){
  int     i, k;

  for (i=0;i<ydim;i++){
    for (k=0;k<xdim;k++){
      if (data[i*xdim+k] > 0)
	      data[i*xdim+k] = -1 * log (data[i*xdim+k]);
      else
	      data[i*xdim+k] = 0;
    }
  }
}

//_____________________________________________________________________________________

int ReconstructSinograms (void){
  int			return_val;

  float   shift[2];

  float  *recon_buffer;

  return_val = 0;

  log_file->ResetTimer (processing_slice_timer);
  log_file->StartTimer (processing_slice_timer);

  log_file->StartTimer (total_processing_slice_timer);

  if (recon_info_record.debug == DEBUG_NO_RECONSTRUCTION)
    return (return_val);

  // These data assignment codes may be useful for the future use of centering_processor.FindCenter (which is not used
  // in the current framework - recon_info_record.centering is always true). Yongsheng Pan, 7/6/2011
  if ( recon_info_record.recon_algorithm == RECONSTRUCTION_GRIDREC_DELAY_ONLY || recon_info_record.recon_algorithm == RECONSTRUCTION_GRIDREC ) { // GridRec

    for (int loop=0;loop<=recon_algorithm->numberOfSinogramsNeeded()-1;loop++){
      if(recon_info_record.gridrec_padding == GRIDREC_PADDING_HALF || recon_info_record.gridrec_padding == GRIDREC_PADDING_BOUNDARY ){
          for( int j = 0; j < recon_info_record.sinogram_ydim; j++ ){
            memcpy( &information.sinograms_boundary_padding[ loop * information.sinogram_adjusted_size * 2 + j * information.sinogram_adjusted_xdim * 2 + information.sinogram_adjusted_xdim / 2 ],
                    &information.sino_calc_buffer[ loop * information.sinogram_adjusted_size + j * information.sinogram_adjusted_xdim ], 
                    sizeof(float) * information.sinogram_adjusted_xdim 
            );
          }

      log_file->Message("Sinogram processed for GridRec using half boundary padding!"); 
      }
      else if(recon_info_record.gridrec_padding == GRIDREC_PADDING_ONE_AND_HALF){
        for( int j = 0; j < recon_info_record.sinogram_ydim; j++ ){
          memcpy( &information.sinograms_boundary_padding[ loop * information.sinogram_adjusted_size * 4 + j * information.sinogram_adjusted_xdim * 4 + information.sinogram_adjusted_xdim * 3 / 2 ],
                  &information.sino_calc_buffer[ loop * information.sinogram_adjusted_size + j * information.sinogram_adjusted_xdim ], 
                  sizeof(float) * information.sinogram_adjusted_xdim 
          );
        }
        log_file->Message("Sinogram processed for GridRec using one and half boundary padding!"); 
      }
      else{
        log_file->Message("Sinogram processed for GridRec using no boundary padding!"); 
      }
    }
  }

  if (recon_info_record.centering){
    if (recon_info_record.use_fixed_shift){
      shift[0] = recon_info_record.fixed_shift_value;
      shift[1] = recon_info_record.fixed_shift_value;

      char msg[256]; 

      sprintf (msg, "recon_info_record.use_fixed_shift is %d, shift[0] is %f ",  
      	       recon_info_record.use_fixed_shift, shift[0]);

      log_file->Message(msg); // debug only      
    }
    else{
      // It seems these codes are not used because recon_info_record.use_fixed_shift is always true in the current framework
      if(recon_info_record.gridrec_padding == GRIDREC_PADDING_HALF ||recon_info_record.gridrec_padding == GRIDREC_PADDING_BOUNDARY ){
	      centering_processor.FindCenter (information.sinograms_boundary_padding, 
                                       &information.sinograms_boundary_padding[information.sinogram_adjusted_size * 2], 
                                       &shift[0], 
                                       &shift[1], 
                                       recon_info_record.ring_removal_coeff
        );
      } 
      else if(recon_info_record.gridrec_padding == GRIDREC_PADDING_ONE_AND_HALF){
	      centering_processor.FindCenter (information.sinograms_boundary_padding, 
                                       &information.sinograms_boundary_padding[information.sinogram_adjusted_size * 4], 
                                       &shift[0], 
                                       &shift[1], 
                                       recon_info_record.ring_removal_coeff
        );
      }
      else{
	      centering_processor.FindCenter (information.sino_calc_buffer, 
                                       &information.sino_calc_buffer[information.sinogram_adjusted_size], 
                                       &shift[0], 
                                       &shift[1], 
                                       recon_info_record.ring_removal_coeff
        );
      }

      // char msg[256]; 
      // sprintf (msg, "centering_processor.FindCenter called");
      // log_file->Message(msg); // debug only      

    }
  }
  else{
    shift[0] = 0;
    shift[1] = 0;

    if(recon_info_record.gridrec_padding == GRIDREC_PADDING_HALF || recon_info_record.gridrec_padding == GRIDREC_PADDING_BOUNDARY ){
      LogProj(information.sinograms_boundary_padding, 
              information.sinogram_adjusted_xdim * 2,  
	            recon_info_record.sinogram_ydim
      );
      LogProj (&information.sinograms_boundary_padding[information.sinogram_adjusted_size * 2],
	              information.sinogram_adjusted_xdim * 2, 
                recon_info_record.sinogram_ydim
      );
    }
    else if( recon_info_record.gridrec_padding == GRIDREC_PADDING_ONE_AND_HALF){
      LogProj(information.sinograms_boundary_padding, 
              information.sinogram_adjusted_xdim * 4,  
	            recon_info_record.sinogram_ydim
      );
      LogProj(&information.sinograms_boundary_padding[information.sinogram_adjusted_size * 4],
	             information.sinogram_adjusted_xdim * 4, 
               recon_info_record.sinogram_ydim
      );
    }
    else{
      LogProj(information.sino_calc_buffer, 
              information.sinogram_adjusted_xdim, 
	            recon_info_record.sinogram_ydim
      );
      LogProj(&information.sino_calc_buffer[information.sinogram_adjusted_size],
	             information.sinogram_adjusted_xdim, 
               recon_info_record.sinogram_ydim
      );
    }

    char msg[256]; 

    sprintf (msg, "recon_info_record.centering is %d, shift[0] is %f, LogProj called ",  
    	     recon_info_record.use_fixed_shift, shift[0]);

    log_file->Message(msg); // debug only      

  }

  for (int loop=0;loop<=recon_algorithm->numberOfSinogramsNeeded()-1;loop++){
    if (recon_info_record.centering){
      if (   (recon_info_record.recon_algorithm != RECONSTRUCTION_GRIDREC_DELAY_ONLY && recon_info_record.recon_algorithm != RECONSTRUCTION_GRIDREC)
      	  || (recon_info_record.recon_algorithm == RECONSTRUCTION_GRIDREC            && recon_info_record.gridrec_padding == GRIDREC_PADDING_NONE) 
      ){ 
        centering_processor.OffCenterCorrSingleManual (&information.sino_calc_buffer[loop*information.sinogram_adjusted_size], (int)shift[loop]);
        log_file->Message("centering_processor.OffCenterCorrSingleManual performed for sino_calc_buffer");      
      }
      else{
	      // Perform centering_processor.OffCenterCorrSingleManual() explicitly for boundary padding
	      LogProj(&information.sino_calc_buffer[loop * information.sinogram_adjusted_size],
		             information.sinogram_adjusted_xdim, 
                 recon_info_record.sinogram_ydim
        );

	      for( int j = 0; j < recon_info_record.sinogram_ydim; j++ ){
	        for( int k = 0; k < information.sinogram_adjusted_xdim; k++ ){
	          information.shifted_recon[j * information.sinogram_adjusted_xdim+ k] = 0.0f;
	        }
	      }

	      for( int j = 0; j < recon_info_record.sinogram_ydim; j++ ){
	        for( int k = 0; k < information.sinogram_adjusted_xdim; k++ ){
            float kk = k - shift[loop]; 
	          int nkk = (int)floor(kk);

	          float fInterpPixel = 0.0f;
	          float fInterpWeight = 0.0f;
	    
            if( recon_info_record.gridrec_padding == GRIDREC_PADDING_BOUNDARY ){
	            if( nkk >= 0 && nkk < information.sinogram_adjusted_xdim ){
		            fInterpPixel += information.sino_calc_buffer[loop * information.sinogram_adjusted_size + j * information.sinogram_adjusted_xdim + nkk ] * (nkk + 1 - kk);
	            }
	            else if( nkk < 0 ){
		            fInterpPixel += information.sino_calc_buffer[loop * information.sinogram_adjusted_size + j * information.sinogram_adjusted_xdim ] * (nkk + 1 - kk);
	            }
	            else if( nkk > information.sinogram_adjusted_xdim  ){
		            fInterpPixel += information.sino_calc_buffer[loop * information.sinogram_adjusted_size + j * information.sinogram_adjusted_xdim + information.sinogram_adjusted_xdim - 1 ] * (nkk + 1 - kk);
	            }

	            if( nkk + 1 >= 0 && nkk + 1 < information.sinogram_adjusted_xdim ){
		            fInterpPixel += information.sino_calc_buffer[loop * information.sinogram_adjusted_size + j * information.sinogram_adjusted_xdim + nkk + 1] * (kk - nkk);
              }
	            else if( nkk + 1 < 0 ){
		            fInterpPixel += information.sino_calc_buffer[loop * information.sinogram_adjusted_size + j * information.sinogram_adjusted_xdim ] * (kk - nkk);
              }
	            else if( nkk + 1 > information.sinogram_adjusted_xdim ){
		            fInterpPixel += information.sino_calc_buffer[loop * information.sinogram_adjusted_size + j * information.sinogram_adjusted_xdim + information.sinogram_adjusted_xdim - 1] * (kk - nkk);
              }
	          }
	          else { // no boundary padding
	            if( nkk >= 0 && nkk < information.sinogram_adjusted_xdim ){
		            fInterpPixel += information.sino_calc_buffer[loop * information.sinogram_adjusted_size + j * information.sinogram_adjusted_xdim + nkk ] * (nkk + 1 - kk);
		            fInterpWeight = nkk + 1 - kk;
	            }

	            // 
	            if( nkk + 1 >= 0 && nkk + 1 < information.sinogram_adjusted_xdim ){
		            fInterpPixel += information.sino_calc_buffer[loop * information.sinogram_adjusted_size + j * information.sinogram_adjusted_xdim + nkk + 1] * (kk - nkk);
		            fInterpWeight += kk - nkk;
              }

	            if( fInterpWeight < 1e-5 )
		            fInterpPixel = 0.0f;
	            else
		            fInterpPixel /= fInterpWeight;
            
	          }
	  
	          information.shifted_sinogram[ j * information.sinogram_adjusted_xdim + k ] = fInterpPixel;
	        }
	      }

	      memcpy(&information.sino_calc_buffer[loop * information.sinogram_adjusted_size], 
	              information.shifted_sinogram, 
	              sizeof(float) * information.sinogram_adjusted_size
        );

	      log_file->Message("Local OffCenterCorrSingleManual performed for sino_calc_buffer");      
      }
      
      if ((   (recon_info_record.debug == DEBUG_FULL) 
           || (recon_info_record.debug == DEBUG_POSTCENTERING) 
           || (recon_info_record.debug == DEBUG_POSTRING)
          ) && (recon_info_record.debug != DEBUG_NO_OUTPUT)){
	      memcpy( &information.debug_post_centering_sinogram_calc_buffer[loop*information.sinogram_adjusted_size],
		            &information.sino_calc_buffer[ loop * information.sinogram_adjusted_size],
		            sizeof(float) * information.sinogram_adjusted_size 
        );
      }

      if (recon_info_record.use_ring_removal){

	      log_file->Message("Removing rings!");

	      if (   (recon_info_record.recon_algorithm != RECONSTRUCTION_GRIDREC_DELAY_ONLY && recon_info_record.recon_algorithm != RECONSTRUCTION_GRIDREC)
	          || ( recon_info_record.recon_algorithm == RECONSTRUCTION_GRIDREC           && recon_info_record.gridrec_padding == GRIDREC_PADDING_NONE  ) 
           ){ 
          centering_processor.RingCorrectionSingle (&information.sino_calc_buffer[loop*information.sinogram_adjusted_size], 
                                                    recon_info_record.ring_removal_coeff
          );
	        log_file->Message("centering_processor.RingCorrectionSingle performed for sino_calc_buffer"); 
	      }
	      else{  // use local version
	        RingCorrectionSingle (&information.sino_calc_buffer[loop*information.sinogram_adjusted_size], 
                                recon_info_record.ring_removal_coeff
          );
	        log_file->Message("Local RingCorrectionSingle performed for sino_calc_buffer");      
	      }

	      if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_POSTCENTERING) || (recon_info_record.debug == DEBUG_POSTRING)) && (recon_info_record.debug != DEBUG_NO_OUTPUT)){
	        memcpy (&information.debug_post_ring_sinogram_calc_buffer[loop*information.sinogram_adjusted_size], 
                  &information.sino_calc_buffer[loop*information.sinogram_adjusted_size], 
                  sizeof(float)*information.sinogram_adjusted_size
          );
	      }
      }
    }

    sprintf (msg, "sinogram number: %d, shift %f", information.sinogram_calc_numbers[loop], shift[loop]);

    log_file->Message(msg);

    // data transfer (again) centering and ring removal
    if ( recon_info_record.recon_algorithm == RECONSTRUCTION_GRIDREC_DELAY_ONLY || recon_info_record.recon_algorithm == RECONSTRUCTION_GRIDREC ) { 
      if(recon_info_record.gridrec_padding == GRIDREC_PADDING_HALF || recon_info_record.gridrec_padding == GRIDREC_PADDING_BOUNDARY ){
	      
        for( int j = 0; j < recon_info_record.sinogram_ydim; j++ ){

	        memcpy( &information.sinograms_boundary_padding[ loop * information.sinogram_adjusted_size * 2 + j * information.sinogram_adjusted_xdim * 2 + information.sinogram_adjusted_xdim / 2 ],
		              &information.sino_calc_buffer[ loop * information.sinogram_adjusted_size + j * information.sinogram_adjusted_xdim ], 
		              sizeof(float) * information.sinogram_adjusted_xdim 
          );

	        for( int k = 0; k < information.sinogram_adjusted_xdim /2; k++ ){
	          information.sinograms_boundary_padding[ loop * information.sinogram_adjusted_size * 2 + j * information.sinogram_adjusted_xdim * 2 + k ] = information.sinograms_boundary_padding[ loop * information.sinogram_adjusted_size * 2 + j * information.sinogram_adjusted_xdim * 2 + information.sinogram_adjusted_xdim / 2 ];
          }

	        for( int k = 0; k < information.sinogram_adjusted_xdim /2; k++ ){
	          information.sinograms_boundary_padding[ loop * information.sinogram_adjusted_size * 2 + j * information.sinogram_adjusted_xdim * 2 + information.sinogram_adjusted_xdim / 2 + information.sinogram_adjusted_xdim + k ] = information.sinograms_boundary_padding[ loop * information.sinogram_adjusted_size * 2 + j * information.sinogram_adjusted_xdim * 2 + information.sinogram_adjusted_xdim / 2 + information.sinogram_adjusted_xdim - 1];
          }

	      }
  
        log_file->Message("Sinogram processed for GridRec using half boundary padding!"); 
      }
      else if(recon_info_record.gridrec_padding == GRIDREC_PADDING_ONE_AND_HALF){
	
        for( int j = 0; j < recon_info_record.sinogram_ydim; j++ ){
	        memcpy( &information.sinograms_boundary_padding[ loop * information.sinogram_adjusted_size * 4 + j * information.sinogram_adjusted_xdim * 4 + information.sinogram_adjusted_xdim * 3 / 2 ],
		              &information.sino_calc_buffer[ loop * information.sinogram_adjusted_size + j * information.sinogram_adjusted_xdim ], 
		              sizeof(float) * information.sinogram_adjusted_xdim 
          );

	        for( int k = 0; k < information.sinogram_adjusted_xdim * 3/2; k++ ){
	        information.sinograms_boundary_padding[ loop * information.sinogram_adjusted_size * 4 + j * information.sinogram_adjusted_xdim * 4 + k ] = information.sinograms_boundary_padding[ loop * information.sinogram_adjusted_size * 4 + j * information.sinogram_adjusted_xdim * 4 + information.sinogram_adjusted_xdim * 3 / 2 ];
          }

	        for( int k = 0; k < information.sinogram_adjusted_xdim *3/2; k++ ){
	          information.sinograms_boundary_padding[ loop * information.sinogram_adjusted_size * 4 + j * information.sinogram_adjusted_xdim * 4 + information.sinogram_adjusted_xdim * 3/ 2 + information.sinogram_adjusted_xdim + k ] = information.sinograms_boundary_padding[ loop * information.sinogram_adjusted_size * 4 + j * information.sinogram_adjusted_xdim * 4 + information.sinogram_adjusted_xdim * 3 / 2 + information.sinogram_adjusted_xdim - 1];
	        }
        }

	      log_file->Message("Sinogram processed for GridRec using one and half boundary padding!"); 
      }
      else{
	      log_file->Message("Sinogram processed for GridRec using no boundary padding!"); 
      }

    }

    if ( recon_info_record.recon_algorithm == RECONSTRUCTION_GRIDREC_DELAY_ONLY || recon_info_record.recon_algorithm == RECONSTRUCTION_GRIDREC ) { // GridRec
      	
      if(recon_info_record.gridrec_padding == GRIDREC_PADDING_HALF || recon_info_record.gridrec_padding == GRIDREC_PADDING_BOUNDARY ){
	      recon_algorithm->setSinoAndReconBuffers(loop+1, 
                                                &information.sinograms_boundary_padding[loop*information.sinogram_adjusted_size * 2], 
                                                &information.reconstructions_boundary_padding[loop*information.reconstruction_size * 4]
        );
      }
      else if(recon_info_record.gridrec_padding == GRIDREC_PADDING_ONE_AND_HALF){
	      recon_algorithm->setSinoAndReconBuffers(loop+1, 
                                                &information.sinograms_boundary_padding[loop*information.sinogram_adjusted_size * 4], 
                                                &information.reconstructions_boundary_padding[loop*information.reconstruction_size * 16]
        );
      }
      else{
	      recon_algorithm->setSinoAndReconBuffers(loop+1, 
                                                &information.sino_calc_buffer[loop*information.sinogram_adjusted_size], 
                                                &information.recon_calc_buffer[loop*information.reconstruction_size]
        );
      }
    }
    else{ // FBP
      recon_algorithm->setSinoAndReconBuffers(loop+1, 
                                              &information.sino_calc_buffer[loop*information.sinogram_adjusted_size], 
                                              &information.recon_calc_buffer[loop*information.reconstruction_size]
      );
    }
  }

  recon_algorithm->reconstruct ();

  if ( recon_info_record.recon_algorithm == RECONSTRUCTION_GRIDREC_DELAY_ONLY || recon_info_record.recon_algorithm == RECONSTRUCTION_GRIDREC ) {
    if(recon_info_record.gridrec_padding == GRIDREC_PADDING_HALF || recon_info_record.gridrec_padding == GRIDREC_PADDING_BOUNDARY ){
      for (int loop=0;loop<=recon_algorithm->numberOfSinogramsNeeded()-1;loop++){
        for (int j=0;j<recon_info_record.reconstruction_ydim;j++){
          // assume reconstruction_xdim = reconstruction_ydim here. Y. Pan 6/20/2011
	        if (shift[loop] >= 0){
	          memcpy(&information.recon_calc_buffer[information.reconstruction_size*loop + j * recon_info_record.reconstruction_xdim ],
		               &information.reconstructions_boundary_padding[ loop * information.reconstruction_size * 4 + ( j + recon_info_record.reconstruction_xdim / 2 ) * recon_info_record.reconstruction_xdim * 2 + recon_info_record.reconstruction_xdim / 2 ], 
		               sizeof(float) * (recon_info_record.reconstruction_xdim) 
            ); 
	        }
	        else{
	          memcpy(&information.recon_calc_buffer[information.reconstruction_size*loop + j * recon_info_record.reconstruction_xdim],
		               &information.reconstructions_boundary_padding[ loop * information.reconstruction_size * 4 + ( j + recon_info_record.reconstruction_xdim / 2 ) * recon_info_record.reconstruction_xdim * 2 + recon_info_record.reconstruction_xdim / 2 ], 
	    	           sizeof(float) * (recon_info_record.reconstruction_xdim) 
            ); 
	        }
	      }
      }
    }

    if(recon_info_record.gridrec_padding == GRIDREC_PADDING_ONE_AND_HALF){
      for (int loop=0;loop<=recon_algorithm->numberOfSinogramsNeeded()-1;loop++){
	      for (int j=0;j<recon_info_record.reconstruction_ydim;j++){

	        if (shift[loop] >= 0){
	          memcpy(&information.recon_calc_buffer[information.reconstruction_size*loop + j * recon_info_record.reconstruction_xdim ],
		               &information.reconstructions_boundary_padding[ loop * information.reconstruction_size * 16 + ( j + recon_info_record.reconstruction_xdim * 3 / 2 ) * recon_info_record.reconstruction_xdim * 4 + recon_info_record.reconstruction_xdim * 3 / 2 + (int)round(shift[loop]) ], 
		               sizeof(float) * (recon_info_record.reconstruction_xdim - (int)round(shift[loop]) ) 
            ); 
	        }
	        else{
	          memcpy(&information.recon_calc_buffer[information.reconstruction_size*loop + j * recon_info_record.reconstruction_xdim + abs( (int)round(shift[loop])) ],
		               &information.reconstructions_boundary_padding[ loop * information.reconstruction_size * 16 + ( j + recon_info_record.reconstruction_xdim * 3 / 2 ) * recon_info_record.reconstruction_xdim * 4 + recon_info_record.reconstruction_xdim * 3 / 2 ], 
		               sizeof(float) * (recon_info_record.reconstruction_xdim - abs( (int)round(shift[loop]) ) ) 
            ); 
	        }
	      }
      }
    }
  }

  // //if we centered, unshift the reconstruction
  for( int j = 0; j < recon_info_record.sinogram_ydim; j++ ){
    for( int k = 0; k < recon_info_record.reconstruction_xdim; k++ ){
      information.shifted_recon[j * recon_info_record.reconstruction_xdim + k] = 0.0f;
    }
  }

  if (recon_info_record.centering){
    for (int loop=0;loop<=recon_algorithm->numberOfSinogramsNeeded()-1;loop++){
      recon_buffer = &information.recon_calc_buffer[information.reconstruction_size*loop];
      if (shift[loop] >= 0){
  	    for (int j=0;j<recon_info_record.reconstruction_ydim;j++)
  	      memcpy (           &information.shifted_recon[j*recon_info_record.reconstruction_xdim], 
		                (void *) &recon_buffer[(j*recon_info_record.reconstruction_xdim)+ (int)round(shift[loop]) ], 
		                sizeof(float)*(recon_info_record.reconstruction_xdim- (int)round(shift[loop]) )
          );
      }
      else {
  	    for (int j=0;j<recon_info_record.reconstruction_ydim;j++)
  	      memcpy (           &information.shifted_recon[(j*recon_info_record.reconstruction_xdim)+abs ((int)round(shift[loop]))], 
                    (void *) &recon_buffer[j*recon_info_record.reconstruction_xdim], 
                    sizeof(float)*(recon_info_record.reconstruction_xdim-abs ((int)round(shift[loop]) ))
          );
      }

      memcpy ((void *) recon_buffer, 
              information.shifted_recon, 
              sizeof(float)*information.reconstruction_size
      );
    }
  }

  // thresholding
  for (int loop=0;loop<=recon_algorithm->numberOfSinogramsNeeded()-1;loop++){
    recon_buffer = &information.recon_calc_buffer[information.reconstruction_size*loop];
    for (int loopy=0;loopy<recon_info_record.reconstruction_ydim;loopy++){
      for (int loopx=0;loopx<recon_info_record.reconstruction_xdim;loopx++){

      	if (recon_buffer[loopy*recon_info_record.reconstruction_xdim+loopx] > data_range_max)
	        data_range_max = recon_buffer[loopy*recon_info_record.reconstruction_xdim+loopx];
	      if (recon_buffer[loopy*recon_info_record.reconstruction_xdim+loopx] < data_range_min)
	        data_range_min = recon_buffer[loopy*recon_info_record.reconstruction_xdim+loopx];
      }
    }
  }

  log_file->StopTimer (total_processing_slice_timer);
  log_file->AccumulateTimer (total_processing_slice_timer);

  log_file->StopTimer (processing_slice_timer);
  log_file->AccumulateTimer (processing_slice_timer);
  //	log_file->TimerMessage (processing_slice_timer);
  log_file->Message ("\n");

  return (return_val);
}

//_____________________________________________________________________________________

void *ReconstructionThread (void *){
  ReconstructSinograms ();
  recon_thread_running = false;

  return (0);
}

//_____________________________________________________________________________________

void MPISendSinogramRequest (int number_to_request, int process_id){
  int				send_command;

  log_file->StartTimer (MPI_requests_timer);

  //Send command
  MPI_Send (&my_id, 1, MPI_INT, process_id, 0, MPI_COMM_WORLD);
  send_command = CLIENT__Request_New_Sinogram;
  MPI_Send (&send_command, 1, MPI_INT, process_id, 0, MPI_COMM_WORLD);
  MPI_Send (&number_to_request, 1, MPI_INT, process_id, 0, MPI_COMM_WORLD);

  log_file->StopTimer (MPI_requests_timer);
  log_file->AccumulateTimer (MPI_requests_timer);

}

//_____________________________________________________________________________________
void MPIRecieveNewSinogram (int *sinogram_number, float *sinogram_shift, float *sinogram_buffer){
  MPI_Status		    status;
  char			    msg[256];
  unsigned short int  white_field_max;
  int					offset;

  log_file->StartTimer (MPI_sinogram_timer);

  log_file->ResetTimer (single_MPI_sinogram_timer);
  log_file->StartTimer (single_MPI_sinogram_timer);

  MPI_Recv (sinogram_shift, 1, MPI_FLOAT, recon_info_record.mayor_id, 0, MPI_COMM_WORLD, &status);

  MPI_Recv (sinogram_data_set, recon_info_record.sinogram_set_size, MPI_BYTE, recon_info_record.mayor_id, 0, MPI_COMM_WORLD, &status);
  recon_info_record.fixed_shift_value = *sinogram_shift;

  log_file->StopTimer (single_MPI_sinogram_timer);
  log_file->AccumulateTimer (single_MPI_sinogram_timer);
  //    log_file->TimerMessage (single_MPI_sinogram_timer);

  offset = 0;
  memcpy (sinogram_number, &sinogram_data_set[offset], sizeof (int));

  offset += sizeof (int);
  memcpy (information.short_sinogram, &sinogram_data_set[offset], sizeof (short) * information.sinogram_size);

  offset += sizeof (short) * information.sinogram_size;
  memcpy (information.short_white_field_sino, &sinogram_data_set[offset], sizeof (short) * recon_info_record.white_size);

  offset += sizeof (short) * recon_info_record.white_size;
  memcpy (information.dark_field_sino_ave, &sinogram_data_set[offset], sizeof (float) * recon_info_record.sinogram_xdim);

  if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_WHITEFIELD)) && (recon_info_record.debug != DEBUG_NO_OUTPUT))
    WriteData (information.short_white_field_sino, *sinogram_number, "whitesino", recon_info_record.sinogram_xdim, recon_info_record.white_size / recon_info_record.sinogram_xdim, NX_UINT16);
  if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_DARKFIELD)) && (recon_info_record.debug != DEBUG_NO_OUTPUT))
    WriteData (information.dark_field_sino_ave, *sinogram_number, "darksinoave", recon_info_record.sinogram_xdim, 1, NX_FLOAT32);
  if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_PRENORM)) && (recon_info_record.debug != DEBUG_NO_OUTPUT))
    WriteData (information.short_sinogram, *sinogram_number, "prenorm", recon_info_record.sinogram_xdim, recon_info_record.sinogram_ydim, NX_UINT16);

  Normalize (information.short_sinogram, information.short_white_field_sino, information.dark_field_sino_ave, sinogram_buffer);

  if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_POSTNORM)) && (recon_info_record.debug != DEBUG_NO_OUTPUT))
    WriteData (sinogram_buffer, *sinogram_number, "postnorm", information.sinogram_adjusted_xdim, recon_info_record.sinogram_ydim, NX_FLOAT32);

  sprintf (msg, "Recieved Sinogram %d with shift %f", *sinogram_number, *sinogram_shift);

  log_file->Message (msg);

  log_file->StopTimer (MPI_sinogram_timer);
  log_file->AccumulateTimer (MPI_sinogram_timer);

}

//_____________________________________________________________________________________

void MPISendClientExiting (int process_id){
  int				send_command;

  MPI_Send (&my_id, 1, MPI_INT, recon_info_record.mayor_id, 0, MPI_COMM_WORLD);
  send_command = CLIENT__Exiting;
  MPI_Send (&send_command, 1, MPI_INT, recon_info_record.mayor_id, 0, MPI_COMM_WORLD);

  //Some closing information for the server to track
  MPI_Send (&data_range_min, 1, MPI_FLOAT, recon_info_record.mayor_id, 0, MPI_COMM_WORLD);
  MPI_Send (&data_range_max, 1, MPI_FLOAT, recon_info_record.mayor_id, 0, MPI_COMM_WORLD);

}

//_____________________________________________________________________________________

int MPIRecieveCommand (int process_id){
  int				recieve_command;
  MPI_Status		status;

  MPI_Recv (&recieve_command, 1, MPI_INT, process_id, 0, MPI_COMM_WORLD, &status);

  return (recieve_command);
}

//_____________________________________________________________________________________

void InitFBP (void){
  //Twice the required space is allocated for buffering purposes
  switch (recon_info_record.recon_algorithm){
    case RECONSTRUCTION_FBP_DELAY_ONLY : 
      recon_algorithm = new FBP (); 
      break;
    case RECONSTRUCTION_FBP_NO_OPTIMIZATION : 
      recon_algorithm = new FBP (); 
      break;
    case RECONSTRUCTION_FBP_OPTIMIZED : 
      recon_algorithm = new OptimizedFBP (); 
      break;
    case RECONSTRUCTION_FBP_CYLINDER_ONLY : 
      recon_algorithm = new CircleFBP (); 
      break;
    default : 
      recon_algorithm = new OptimizedFBP (); 
      break;
  }

  information.sinogram_numbers = (int *) malloc (sizeof(int)*recon_algorithm->numberOfSinogramsNeeded()*2);
  if (information.sinogram_numbers == NULL){
    sprintf (msg, "Could not alocate memory for information.sinogram_numbers.");
    error_log->addError (msg, "InitFBP ()");
  }

  for (int loop=0;loop<recon_algorithm->numberOfSinogramsNeeded()*2;loop++)
    information.sinogram_numbers[loop] = -123;

  information.sinogram_shifts = (float *) malloc (sizeof(float)*recon_algorithm->numberOfSinogramsNeeded()*2);
  if (information.sinogram_shifts == NULL){
    sprintf (msg, "Could not alocate memory for information.sinogram_shifts.");
    error_log->addError (msg, "InitFBP ()");
  }

  for (int loop=0;loop<recon_algorithm->numberOfSinogramsNeeded()*2;loop++)
    information.sinogram_shifts[loop] = 0.0f;

  information.sinograms = (float *) malloc (sizeof(float)*information.sinogram_adjusted_size*(recon_algorithm->numberOfSinogramsNeeded()*2));
  if (information.sinograms == NULL){
    sprintf (msg, "Could not alocate memory for information.sinograms.");
    error_log->addError (msg, "InitFBP ()");
  }

  information.reconstructions = (float *) malloc (sizeof(float)*information.reconstruction_size*(recon_algorithm->numberOfSinogramsNeeded()*2));
  if (information.reconstructions == NULL){
    sprintf (msg, "Could not alocate memory for information.reconstructions.");
    error_log->addError (msg, "InitFBP ()");
  }

  information.sinogram_mpi_numbers = information.sinogram_numbers;
  information.sinogram_calc_numbers = &information.sinogram_numbers[recon_algorithm->numberOfSinogramsNeeded()];

  information.sinogram_mpi_shifts = information.sinogram_shifts;
  information.sinogram_calc_shifts = &information.sinogram_shifts[recon_algorithm->numberOfSinogramsNeeded()];

  information.sino_mpi_buffer = information.sinograms;
  information.sino_calc_buffer = &information.sinograms[information.sinogram_adjusted_size*recon_algorithm->numberOfSinogramsNeeded()];

  information.recon_mpi_buffer = information.reconstructions;
  information.recon_calc_buffer = &information.reconstructions[information.reconstruction_size*recon_algorithm->numberOfSinogramsNeeded()];

  if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_POSTCENTERING) || (recon_info_record.debug == DEBUG_POSTRING)) && (recon_info_record.debug != DEBUG_NO_OUTPUT)){
    
    information.debug_post_centering_sinogram = (float *) malloc (sizeof(float)*information.sinogram_adjusted_size*(recon_algorithm->numberOfSinogramsNeeded()*2));
    if (information.debug_post_centering_sinogram == NULL){
	    sprintf (msg, "Could not alocate memory for information.debug_post_centering_sinogram.");
	    error_log->addError (msg, "InitFBP ()");
    }

    information.debug_post_ring_sinogram = (float *) malloc (sizeof(float)*information.sinogram_adjusted_size*(recon_algorithm->numberOfSinogramsNeeded()*2));
    if (information.debug_post_ring_sinogram == NULL){
	    sprintf (msg, "Could not alocate memory for information.debug_post_ring_sinogram.");
	    error_log->addError (msg, "InitFBP ()");
    }

    information.debug_post_centering_sinogram_mpi_buffer = information.debug_post_centering_sinogram;
    information.debug_post_centering_sinogram_calc_buffer = &information.debug_post_centering_sinogram[information.sinogram_adjusted_size*recon_algorithm->numberOfSinogramsNeeded()];

    information.debug_post_ring_sinogram_mpi_buffer = information.debug_post_ring_sinogram;
    information.debug_post_ring_sinogram_calc_buffer = &information.debug_post_ring_sinogram[information.sinogram_adjusted_size*recon_algorithm->numberOfSinogramsNeeded()];
  }

  recon_algorithm->setSinogramDimensions(information.sinogram_adjusted_xdim, recon_info_record.sinogram_ydim);
  recon_algorithm->setThetaList (recon_info_record.theta_list, recon_info_record.theta_list_size);
  recon_algorithm->setFilter (recon_info_record.filter);

  recon_algorithm->init();

  centering_processor.init (recon_algorithm);

}

//_____________________________________________________________________________________

void InitGridrec (void){
  int				loop;

  recon_algorithm = new GridRec ();

  information.sinogram_numbers = (int *) malloc (sizeof(int)*recon_algorithm->numberOfSinogramsNeeded()*2);
  if (information.sinogram_numbers == NULL){
    sprintf (msg, "Could not alocate memory for information.sinogram_numbers.");
    error_log->addError (msg, "InitGridrec ()");
  }

  for (loop=0;loop<recon_algorithm->numberOfSinogramsNeeded()*2;loop++)
    information.sinogram_numbers[loop] = -123;

  information.sinogram_shifts = (float *) malloc (sizeof(float)*recon_algorithm->numberOfSinogramsNeeded()*2);
  if (information.sinogram_shifts == NULL){
    sprintf (msg, "Could not alocate memory for information.sinogram_shifts.");
    error_log->addError (msg, "InitGridrec ()");
  }

  for (int loop=0;loop<recon_algorithm->numberOfSinogramsNeeded()*2;loop++)
    information.sinogram_shifts[loop] = 0.0f;

  //Double the required space for buffering purposes
  information.sinograms = (float *) malloc (sizeof(float)*information.sinogram_adjusted_size*(recon_algorithm->numberOfSinogramsNeeded()*2));
  if (information.sinograms == NULL){
    sprintf (msg, "Could not alocate memory for information.sinograms.");
    error_log->addError (msg, "InitGridrec ()");
  }

  information.reconstructions = (float *) malloc (sizeof(float)*information.reconstruction_size*(recon_algorithm->numberOfSinogramsNeeded()*2));
  if (information.reconstructions == NULL){
    sprintf (msg, "Could not alocate memory for information.reconstructions.");
    error_log->addError (msg, "InitGridrec ()");
  }

  information.sinogram_mpi_numbers = information.sinogram_numbers;
  information.sinogram_calc_numbers = &information.sinogram_numbers[recon_algorithm->numberOfSinogramsNeeded()];

  information.sinogram_mpi_shifts = information.sinogram_shifts;
  information.sinogram_calc_shifts = &information.sinogram_shifts[recon_algorithm->numberOfSinogramsNeeded()];

  information.sino_mpi_buffer = information.sinograms;
  information.sino_calc_buffer = &information.sinograms[information.sinogram_adjusted_size*recon_algorithm->numberOfSinogramsNeeded()];

  information.recon_mpi_buffer = information.reconstructions;
  information.recon_calc_buffer = &information.reconstructions[information.reconstruction_size*recon_algorithm->numberOfSinogramsNeeded()];

  if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_POSTCENTERING) || (recon_info_record.debug == DEBUG_POSTRING)) && (recon_info_record.debug != DEBUG_NO_OUTPUT)){
    information.debug_post_centering_sinogram = (float *) malloc (sizeof(float)*information.sinogram_adjusted_size*(recon_algorithm->numberOfSinogramsNeeded()*2));
    if (information.debug_post_centering_sinogram == NULL){
      sprintf (msg, "Could not alocate memory for information.debug_post_centering_sinogram.");
      error_log->addError (msg, "InitGridrec ()");
    }

    information.debug_post_ring_sinogram = (float *) malloc (sizeof(float)*information.sinogram_adjusted_size*(recon_algorithm->numberOfSinogramsNeeded()*2));
    if (information.debug_post_ring_sinogram == NULL){
      sprintf (msg, "Could not alocate memory for information.debug_post_ring_sinogram.");
      error_log->addError (msg, "InitGridrec ()");
    }

    information.debug_post_centering_sinogram_mpi_buffer = information.debug_post_centering_sinogram;
    information.debug_post_centering_sinogram_calc_buffer = &information.debug_post_centering_sinogram[information.sinogram_adjusted_size*recon_algorithm->numberOfSinogramsNeeded()];

    information.debug_post_ring_sinogram_mpi_buffer = information.debug_post_ring_sinogram;
    information.debug_post_ring_sinogram_calc_buffer = &information.debug_post_ring_sinogram[information.sinogram_adjusted_size*recon_algorithm->numberOfSinogramsNeeded()];
  }

  //Double the required space for 0.5 boundary padding
  if(recon_info_record.gridrec_padding == GRIDREC_PADDING_HALF || recon_info_record.gridrec_padding == GRIDREC_PADDING_BOUNDARY ){

    information.sinograms_boundary_padding = (float *) malloc (sizeof(float)*information.sinogram_adjusted_size*(recon_algorithm->numberOfSinogramsNeeded()*2));
    if (information.sinograms_boundary_padding == NULL){
      sprintf (msg, "Could not alocate memory for information.sinograms_boundary_padding.");
      error_log->addError (msg, "InitGridrec ()");
    }

    information.reconstructions_boundary_padding = (float *) malloc (sizeof(float)*information.reconstruction_size*(recon_algorithm->numberOfSinogramsNeeded()*4));
    if (information.reconstructions_boundary_padding == NULL){
      sprintf (msg, "Could not alocate memory for information.reconstructions.");
      error_log->addError (msg, "InitGridrec ()");
    }
  }

  //quadruple the required space for 1.5 boundary padding
  char msg1[256]; 
  sprintf (msg1, "In InitGridrec, recon_info_record.gridrec_padding is %d", recon_info_record.gridrec_padding);	
  log_file->Message (msg1);

  if(recon_info_record.gridrec_padding == GRIDREC_PADDING_ONE_AND_HALF){

    information.sinograms_boundary_padding = (float *) malloc (sizeof(float)*information.sinogram_adjusted_size*(recon_algorithm->numberOfSinogramsNeeded()*4));
    if (information.sinograms_boundary_padding == NULL){
    	sprintf (msg, "Could not alocate memory for information.sinograms_boundary_padding.");
    	error_log->addError (msg, "InitGridrec ()");
    }

    information.reconstructions_boundary_padding = (float *) malloc (sizeof(float)*information.reconstruction_size*(recon_algorithm->numberOfSinogramsNeeded()*16));
    if (information.reconstructions_boundary_padding == NULL){
    	sprintf (msg, "Could not alocate memory for information.reconstructions_boundary_padding.");
    	error_log->addError (msg, "InitGridrec ()");
    }
  }

  // used for ring correction of size information.sinogram_adjusted_size
  information.mean_vect = (float  *) malloc (sizeof(float)*recon_info_record.sinogram_ydim); 
  information.low_pass_sino_lines_data = (float  *) malloc (sizeof(float)*information.sinogram_adjusted_xdim); 
  information.mean_sino_line_data = (float *) malloc (sizeof(float)*information.sinogram_adjusted_xdim); 

  if(recon_info_record.gridrec_padding == GRIDREC_PADDING_HALF || recon_info_record.gridrec_padding == GRIDREC_PADDING_BOUNDARY ){
    recon_algorithm->setSinogramDimensions(information.sinogram_adjusted_xdim * 2, 
                                           recon_info_record.sinogram_ydim);
  }
  else if(recon_info_record.gridrec_padding == GRIDREC_PADDING_ONE_AND_HALF){
    recon_algorithm->setSinogramDimensions(information.sinogram_adjusted_xdim * 4, 
                                           recon_info_record.sinogram_ydim);
  }
  else{
    recon_algorithm->setSinogramDimensions(information.sinogram_adjusted_xdim, 
                                           recon_info_record.sinogram_ydim);
  }

  recon_algorithm->setThetaList (recon_info_record.theta_list, recon_info_record.theta_list_size);
  recon_algorithm->setFilter (recon_info_record.filter);

  recon_algorithm->init();

  centering_processor.init (recon_algorithm);

}

//_____________________________________________________________________________________

void InitializeClient (void){
  information.shifted_recon = (float *) malloc (sizeof (float)*information.reconstruction_size);
  if (information.shifted_recon == NULL){
    sprintf (msg, "Could not alocate memory for information.shifted_recon.");
    error_log->addError (msg, "InitializeClient ()");
  }

  information.shifted_sinogram = (float *) malloc (sizeof (float)*information.sinogram_adjusted_size);
  if (information.shifted_sinogram == NULL){
    sprintf (msg, "Could not alocate memory for information.shifted_sinogram.");
    error_log->addError (msg, "InitializeClient ()");
  }

  switch (recon_info_record.recon_algorithm){

    case RECONSTRUCTION_FBP_DELAY_ONLY      : InitFBP();      break;
    case RECONSTRUCTION_GRIDREC_DELAY_ONLY  : InitGridrec();  break;
    case RECONSTRUCTION_FBP_NO_OPTIMIZATION : InitFBP();      break;
    case RECONSTRUCTION_FBP_OPTIMIZED       : InitFBP();      break;
    case RECONSTRUCTION_FBP_CYLINDER_ONLY   : InitFBP();      break;
    case RECONSTRUCTION_GRIDREC             : InitGridrec (); break;

    default : break;
    }

}

//_____________________________________________________________________________________

void ProcessingLoop (void){
  void  *ptr_temp_buffer;

  int   mpi_command,
        sinograms_recieved,
        reconstructions_sent;

  pthread_t   reconstruction_thread_handle;
  
  float   *saved_sinograms;

  //Alert server process we're ready for a sinogram
  //request num_buffers sinograms to get us going
  MPISendSinogramRequest (recon_algorithm->numberOfSinogramsNeeded(), recon_info_record.mayor_id);

  saved_sinograms = NULL;

  sinograms_recieved = 0;
  while(sinograms_recieved != recon_algorithm->numberOfSinogramsNeeded()){
    //This shoulde be the Sending_Sinogram command
    mpi_command = MPIRecieveCommand (recon_info_record.mayor_id);
    //If we recieve a stop here, we're not being sent anything...

    if (mpi_command == STOP){
      log_file->Message("Recieved a stop command in ProcessingLoop location 1.");
      MPISendClientExiting (recon_info_record.mayor_id);
      return;
	  }
      
    //We'll put the appropriate buffer
    MPIRecieveNewSinogram (&information.sinogram_mpi_numbers[sinograms_recieved], 
                           &information.sinogram_mpi_shifts[sinograms_recieved], 
                           &information.sino_mpi_buffer[sinograms_recieved*information.sinogram_adjusted_size]
    );

    sinograms_recieved++;
  }

  //this should be the start calculations commend...
  mpi_command = MPIRecieveCommand (recon_info_record.mayor_id);
  if (mpi_command == STOP){
    log_file->Message("Recieved a stop command in ProcessingLoop location 2.");
    MPISendClientExiting (recon_info_record.mayor_id);
    return;
  }

  //Swap buffers...
  ptr_temp_buffer = (void *) information.sinogram_mpi_numbers;
  information.sinogram_mpi_numbers = information.sinogram_calc_numbers;
  information.sinogram_calc_numbers = (int *) ptr_temp_buffer;

  ptr_temp_buffer = (void *) information.sinogram_mpi_shifts;
  information.sinogram_mpi_shifts = information.sinogram_calc_shifts;
  information.sinogram_calc_shifts = (float *) ptr_temp_buffer;

  ptr_temp_buffer = (void *) information.sino_mpi_buffer;
  information.sino_mpi_buffer = information.sino_calc_buffer;
  information.sino_calc_buffer = (float *) ptr_temp_buffer;

  ptr_temp_buffer = (void *) information.recon_mpi_buffer;
  information.recon_mpi_buffer = information.recon_calc_buffer;
  information.recon_calc_buffer = (float *) ptr_temp_buffer;

  if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_POSTCENTERING) || (recon_info_record.debug == DEBUG_POSTRING)) && (recon_info_record.debug != DEBUG_NO_OUTPUT)){
    ptr_temp_buffer = (void *) information.debug_post_centering_sinogram_mpi_buffer;
    information.debug_post_centering_sinogram_mpi_buffer = information.debug_post_centering_sinogram_calc_buffer;
    information.debug_post_centering_sinogram_calc_buffer = (float *) ptr_temp_buffer;

    ptr_temp_buffer = (void *) information.debug_post_ring_sinogram_mpi_buffer;
    information.debug_post_ring_sinogram_mpi_buffer = information.debug_post_centering_sinogram_calc_buffer;
    information.debug_post_ring_sinogram_calc_buffer = (float *) ptr_temp_buffer;
  }

  //Start processing
  log_file->Message("Starting reconstruction thread!");
  recon_thread_running = true;
  pthread_create (&reconstruction_thread_handle, NULL, ReconstructionThread, NULL);

  //request sinograms to buffer
  MPISendSinogramRequest (recon_algorithm->numberOfSinogramsNeeded(), recon_info_record.mayor_id);

  sinograms_recieved = 0;
  while(sinograms_recieved != recon_algorithm->numberOfSinogramsNeeded()){
    //This shoulde be the Sending_Sinogram command
    mpi_command = MPIRecieveCommand (recon_info_record.mayor_id);
    //If we recieve a stop here, we should have recieved 2 sinograms to reconstruct--but then should exit--we won't be sent any more
    if (mpi_command != STOP){
      //We'll put the appropriate buffer
      MPIRecieveNewSinogram (&information.sinogram_mpi_numbers[sinograms_recieved], &information.sinogram_mpi_shifts[sinograms_recieved], &information.sino_mpi_buffer[sinograms_recieved*information.sinogram_adjusted_size]);
      sinograms_recieved++;
	  }
    else{
      log_file->Message("Recieved a stop command in ProcessingLoop location 3.");
      sinograms_recieved = recon_algorithm->numberOfSinogramsNeeded();
	  }
  }

  //if mpi_command == STOP, there will be no new command to kick off the loop...otherwise, recieve a new command to kick off the loop...
  if (mpi_command != STOP)
    mpi_command = MPIRecieveCommand (recon_info_record.mayor_id);

  while (mpi_command != STOP){
      if (mpi_command == SERVER__Sending_Sinogram){
	      log_file->Message("Trying to recieve new sinograms");

	      sinograms_recieved = 0;
	      while(sinograms_recieved != recon_algorithm->numberOfSinogramsNeeded()){

          //We'll put the appropriate buffer
          MPIRecieveNewSinogram (&information.sinogram_mpi_numbers[sinograms_recieved], 
                                 &information.sinogram_mpi_shifts[sinograms_recieved], 
                                 &information.sino_mpi_buffer[sinograms_recieved*information.sinogram_adjusted_size]
          );

          //This shoulde be the Sending_Sinogram command
          mpi_command = MPIRecieveCommand (recon_info_record.mayor_id);

          sinograms_recieved++;
	      }

	      log_file->Message("Recieved all sinograms");
	    }

      if (mpi_command == SERVER__Start_Calculations){
	      if (recon_thread_running){
	        pthread_join (reconstruction_thread_handle, NULL);
	        recon_thread_running = false;
	      }

        //Swap buffers...
        ptr_temp_buffer = (void *) information.sinogram_mpi_numbers;
        information.sinogram_mpi_numbers = information.sinogram_calc_numbers;
        information.sinogram_calc_numbers = (int *) ptr_temp_buffer;

        ptr_temp_buffer = (void *) information.sinogram_mpi_shifts;
        information.sinogram_mpi_shifts = information.sinogram_calc_shifts;
        information.sinogram_calc_shifts = (float *) ptr_temp_buffer;

        ptr_temp_buffer = (void *) information.sino_mpi_buffer;
        information.sino_mpi_buffer = information.sino_calc_buffer;
        information.sino_calc_buffer = (float *) ptr_temp_buffer;

        ptr_temp_buffer = (void *) information.recon_mpi_buffer;
        information.recon_mpi_buffer = information.recon_calc_buffer;
        information.recon_calc_buffer = (float *) ptr_temp_buffer;

	      if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_POSTCENTERING) || (recon_info_record.debug == DEBUG_POSTRING)) && (recon_info_record.debug != DEBUG_NO_OUTPUT)){
          ptr_temp_buffer = (void *) information.debug_post_centering_sinogram_mpi_buffer;
          information.debug_post_centering_sinogram_mpi_buffer = information.debug_post_centering_sinogram_calc_buffer;
          information.debug_post_centering_sinogram_calc_buffer = (float *) ptr_temp_buffer;

          ptr_temp_buffer = (void *) information.debug_post_ring_sinogram_mpi_buffer;
          information.debug_post_ring_sinogram_mpi_buffer = information.debug_post_centering_sinogram_calc_buffer;
          information.debug_post_ring_sinogram_calc_buffer = (float *) ptr_temp_buffer;
	      }

        //Start processing
        log_file->Message("Starting reconstruction thread!");
        recon_thread_running = true;
        pthread_create (&reconstruction_thread_handle, NULL, ReconstructionThread, NULL);

	      reconstructions_sent = 0;
	      while(reconstructions_sent != recon_algorithm->numberOfSinogramsNeeded()){
	        WriteReconstruction (&information.recon_mpi_buffer[reconstructions_sent*information.reconstruction_size], 
                                information.sinogram_mpi_numbers[reconstructions_sent], 
                                information.sinogram_mpi_shifts[reconstructions_sent]
          );

	        if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_POSTCENTERING)) && (recon_info_record.debug != DEBUG_NO_OUTPUT))
		        WriteData (&information.debug_post_centering_sinogram_mpi_buffer[reconstructions_sent*information.sinogram_adjusted_size], 
                        information.sinogram_mpi_numbers[reconstructions_sent], "postcenter", 
                        information.sinogram_adjusted_xdim, 
                        recon_info_record.sinogram_ydim, 
                        NX_FLOAT32
            );

	        if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_POSTRING)) && (recon_info_record.debug != DEBUG_NO_OUTPUT))
		        WriteData (&information.debug_post_ring_sinogram_mpi_buffer[reconstructions_sent*information.sinogram_adjusted_size], 
                        information.sinogram_mpi_numbers[reconstructions_sent], 
                        "postring", 
                        information.sinogram_adjusted_xdim,
                        recon_info_record.sinogram_ydim, 
                        NX_FLOAT32
            );

	        reconstructions_sent++;
	      }

	      //Request new sinogram
	      MPISendSinogramRequest (recon_algorithm->numberOfSinogramsNeeded(), recon_info_record.mayor_id);
	    }

      //Recieve server command
      mpi_command = MPIRecieveCommand (recon_info_record.mayor_id);
    }

  if (recon_thread_running){

      log_file->TimeStamp ("Waiting for reconstruction thread to exit");
      if (recon_thread_running){
        pthread_join (reconstruction_thread_handle, NULL);
        recon_thread_running = false;
      }

      reconstructions_sent = 0;
      while (reconstructions_sent != recon_algorithm->numberOfSinogramsNeeded()){
	      
        WriteReconstruction (&information.recon_calc_buffer[reconstructions_sent*information.reconstruction_size], 
                              information.sinogram_calc_numbers[reconstructions_sent], 
                              information.sinogram_calc_shifts[reconstructions_sent]
        );

	      if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_POSTCENTERING)) && (recon_info_record.debug != DEBUG_NO_OUTPUT))
	        WriteData (&information.debug_post_centering_sinogram_mpi_buffer[reconstructions_sent*information.sinogram_adjusted_size], 
                      information.sinogram_mpi_numbers[reconstructions_sent], 
                      "postcenter", 
                      recon_info_record.sinogram_xdim, 
                      recon_info_record.sinogram_ydim, 
                      NX_FLOAT32
          );

	      if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_POSTRING)) && (recon_info_record.debug != DEBUG_NO_OUTPUT))
	        WriteData (&information.debug_post_ring_sinogram_mpi_buffer[reconstructions_sent*information.sinogram_adjusted_size], 
                      information.sinogram_mpi_numbers[reconstructions_sent], 
                      "postring", 
                      recon_info_record.sinogram_xdim, 
                      recon_info_record.sinogram_ydim, 
                      NX_FLOAT32
          );

	      reconstructions_sent++;
	    }
  }
  else{
    reconstructions_sent = 0;
    while (reconstructions_sent != recon_algorithm->numberOfSinogramsNeeded()){
	    WriteReconstruction (&information.recon_calc_buffer[reconstructions_sent*information.reconstruction_size], 
                            information.sinogram_calc_numbers[reconstructions_sent], 
                            information.sinogram_calc_shifts[reconstructions_sent]
      );

	  if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_POSTCENTERING)) && (recon_info_record.debug != DEBUG_NO_OUTPUT))
	    WriteData (&information.debug_post_centering_sinogram_mpi_buffer[reconstructions_sent*information.sinogram_adjusted_size], 
                  information.sinogram_mpi_numbers[reconstructions_sent], 
                  "postcenter", 
                  recon_info_record.sinogram_xdim, 
                  recon_info_record.sinogram_ydim, 
                  NX_FLOAT32
      );

	  if (((recon_info_record.debug == DEBUG_FULL) || (recon_info_record.debug == DEBUG_POSTRING)) && (recon_info_record.debug != DEBUG_NO_OUTPUT))
	    WriteData (&information.debug_post_ring_sinogram_mpi_buffer[reconstructions_sent*information.sinogram_adjusted_size], 
                  information.sinogram_mpi_numbers[reconstructions_sent], 
                  "postring", 
                  recon_info_record.sinogram_xdim, 
                  recon_info_record.sinogram_ydim, 
                  NX_FLOAT32
      );

	  reconstructions_sent++;
	  }
  }

  if (saved_sinograms != NULL)
    free (saved_sinograms);

  MPISendClientExiting (recon_info_record.mayor_id);
}

//_____________________________________________________________________________________

void DestroyClient (void){
  log_file->Message ("Destroying client process...");

  if (information.debug_post_centering_sinogram != NULL)
    free (information.debug_post_centering_sinogram);
  if (information.debug_post_ring_sinogram != NULL)
    free (information.debug_post_ring_sinogram);

  if (information.sinogram_numbers != NULL)
    free (information.sinogram_numbers);

  if (information.sinogram_shifts != NULL)
    free (information.sinogram_shifts);

  log_file->Message ("sinogram_shifts is destroyed");

  if (recon_info_record.theta_list != NULL)
    free (recon_info_record.theta_list);

  if (information.sinograms != NULL)
    free (information.sinograms);

  if (information.shifted_recon != NULL)
    free (information.shifted_recon);

  log_file->Message ("shifted_recon is destroyed");

  if (information.shifted_sinogram != NULL) 
    free (information.shifted_sinogram);

  log_file->Message ("shifted_sinogram is destroyed");

  recon_algorithm->destroy();

  log_file->Message ("recon_algorithm is destroyed");

  if (information.short_sinogram != NULL)
    free (information.short_sinogram);
  if (information.short_white_field_sino != NULL)
    free (information.short_white_field_sino);
  if (information.float_white_field_sino != NULL)
    free (information.float_white_field_sino);
  if (information.dark_field_sino_ave != NULL)
    free (information.dark_field_sino_ave);

  if (sinogram_data_set != NULL)
    free (sinogram_data_set);

  log_file->Message ("sinogram_data_set is destroying...");

  if (information.int_scale_buffer != NULL)
    free (information.int_scale_buffer);

  delete (recon_file);

  if(recon_info_record.gridrec_padding == GRIDREC_PADDING_HALF || 
     recon_info_record.gridrec_padding == GRIDREC_PADDING_BOUNDARY || 
     recon_info_record.gridrec_padding == GRIDREC_PADDING_ONE_AND_HALF){  

    if (information.sinograms_boundary_padding != NULL)
      free (information.sinograms_boundary_padding);
    if (information.reconstructions_boundary_padding != NULL)
      free (information.reconstructions_boundary_padding);
  }

  log_file->Message ("Before destroying mean_vect...");

  if(information.mean_vect != NULL)
    free(information.mean_vect);

  if(information.low_pass_sino_lines_data != NULL)
    free(information.low_pass_sino_lines_data);
  
  if(information.mean_sino_line_data != NULL)
    free(information.mean_sino_line_data);

  if (information.reconstructions != NULL)
    free (information.reconstructions);

  log_file->Message ("Client destroyed");
}

//_____________________________________________________________________________________

void CreateClientTimers (void){
  processing_slice_timer = log_file->CreateTimer ("Single_Slice");
  log_file->ResetTimer (processing_slice_timer);

  total_processing_slice_timer = log_file->CreateTimer ("All_Slices");
  log_file->ResetTimer (total_processing_slice_timer);

  MPI_requests_timer = log_file->CreateTimer ("MPI_Requests");
  log_file->ResetTimer (MPI_requests_timer);

  MPI_sinogram_timer = log_file->CreateTimer ("MPI_Recieve_Sinogram");
  log_file->ResetTimer (MPI_sinogram_timer);

  single_MPI_sinogram_timer = log_file->CreateTimer ("Single_MPI_Recieve_Sinogram");
  log_file->ResetTimer (single_MPI_sinogram_timer);

  write_reconstruction_timer = log_file->CreateTimer ("Write_Reconstruction");
  log_file->ResetTimer (write_reconstruction_timer);
}

//_____________________________________________________________________________________

void DestroyClientTimers (void){
  log_file->DestroyTimer (processing_slice_timer      );
  log_file->DestroyTimer (total_processing_slice_timer);
  log_file->DestroyTimer (MPI_requests_timer          );
  log_file->DestroyTimer (MPI_sinogram_timer          );
  log_file->DestroyTimer (single_MPI_sinogram_timer   );
  log_file->DestroyTimer (write_reconstruction_timer  );
}

//_____________________________________________________________________________________

void ClientProcess (void){
  char    method[256];

  recon_file = new NexusBoxClass ();
  recon_file->SetFileMode (HDF5_MODE);

  FirstContact ();

  CreateClientTimers ();

  switch (recon_info_record.recon_algorithm){
    case RECONSTRUCTION_FBP_DELAY_ONLY      : strcpy (method, "RECONSTRUCTION_FBP_DELAY_ONLY"    ); break;
    case RECONSTRUCTION_GRIDREC_DELAY_ONLY  : strcpy (method, "RECONSTRUCTION_GRIDREC_DELAY_ONLY"); break;
    case RECONSTRUCTION_FBP_NO_OPTIMIZATION : strcpy (method, "RECONSTRUCTION_NO_OPTIMIZATION"   ); break;
    case RECONSTRUCTION_FBP_OPTIMIZED       : strcpy (method, "RECONSTRUCTION_OPTIMIZED"         ); break;
    case RECONSTRUCTION_FBP_CYLINDER_ONLY   : strcpy (method, "RECONSTRUCTION_CYLINDER_ONLY"     ); break;
    case RECONSTRUCTION_GRIDREC             : strcpy (method, "RECONSTRUCTION_GRIDREC"           ); break;
  }
  sprintf (msg, "I will be using method %s.\n", method);
  log_file->Message (msg);

  InitializeClient ();

  log_file->TimeStamp ("ClientInitialized");  // debug

  ProcessingLoop ();

  log_file->TimeStamp ("Loop processed");  // debug

  DestroyClient ();

  log_file->TimeStamp ("Exiting client process");
}

//_____________________________________________________________________________________

int StartTomoMPIClient (int argc, char* argv[], char *this_log_file_name, char * this_err_file_name){
  char        msg[256];

  strcpy (log_file_name, this_log_file_name);
  strcpy (error_file_name, this_err_file_name);

  ClientAcknowledgements ();

  num_reconstructions_processed = 0;

  recon_info_record.theta_list = NULL;

  sinogram_data_set = NULL;
  information.sinogram_numbers = NULL;
  information.sinogram_shifts = NULL;
  information.sinograms = NULL;
  information.reconstructions = NULL;

  information.sinograms_boundary_padding = NULL;
  information.reconstructions_boundary_padding = NULL;

  information.mean_vect = NULL;  
  information.low_pass_sino_lines_data = NULL;  
  information.mean_sino_line_data = NULL; 

  information.short_sinogram = NULL;
  information.short_white_field_sino = NULL;
  information.dark_field_sino_ave = NULL;
  information.shifted_recon = NULL;
  information.shifted_sinogram = NULL; 
  information.debug_post_centering_sinogram = NULL;
  information.debug_post_ring_sinogram = NULL;
  information.int_scale_buffer = NULL;

  ClientProcess();

  sprintf (msg, "Process %s with id %d has had enough and quits.", processor_name, my_id);
  log_file->TimeStamp (msg);

  sprintf (msg, "Reconstructions processed: %d", num_reconstructions_processed);
  log_file->TimeStamp (msg);

  log_file->TimerMessage (write_reconstruction_timer);
  log_file->TimerMessage (MPI_requests_timer);
  log_file->TimerMessage (MPI_sinogram_timer);
  log_file->TimerMessage (total_processing_slice_timer);

  DestroyClientTimers ();

  delete (log_file);

  return 0;
}
//_____________________________________________________________________________________


void RingCorrectionSingle (float *data, float ring_coeff) { 
  int         i, j, m; 
  float       mean_total; 
  float       tmp; 
 
  for (m=0;m<20;m++) {

    // normalization of each projection: mean values estimation 
    for (i=0;i<recon_info_record.sinogram_ydim;i++) 
      information.mean_vect[i] = 0.0; 
    mean_total = 0.0; 
 
    for (i=0;i<recon_info_record.sinogram_ydim;i++) {
      for (j=0;j<information.sinogram_adjusted_xdim;j++) {
	      information.mean_vect[i] += data[i*information.sinogram_adjusted_xdim+j]; 
      }

      information.mean_vect[i] /= information.sinogram_adjusted_xdim; 
      mean_total += information.mean_vect[i]; 
    } 
    mean_total /= recon_info_record.sinogram_ydim; 
 
    // renormalization of each projection to the global mean 
    for (i=0;i<recon_info_record.sinogram_ydim;i++) {
      for (j=0;j<information.sinogram_adjusted_xdim;j++) {
  	    if (information.mean_vect[i] != 0.0) {
	        data[i*information.sinogram_adjusted_xdim+j] = data[i*information.sinogram_adjusted_xdim+j]*mean_total/information.mean_vect[i];        // ring filtering: sum of projection and low-pass filter of the result 
        }
      }
    }

    for (i=0;i<information.sinogram_adjusted_xdim;i++) 
      information.mean_sino_line_data[i] = 0.0; 
 
    for (i=0;i<recon_info_record.sinogram_ydim;i++) 
      for (j=0;j<information.sinogram_adjusted_xdim;j++) 
	      information.mean_sino_line_data[j] += data[i*information.sinogram_adjusted_xdim+j]; 
 
    for (i=0;i<information.sinogram_adjusted_xdim;i++) 
      information.mean_sino_line_data[i] /= recon_info_record.sinogram_ydim; 
 
    for (j=1;j<information.sinogram_adjusted_xdim-1;j++) {
      information.low_pass_sino_lines_data[j] = (information.mean_sino_line_data[j-1]+information.mean_sino_line_data[j]+information.mean_sino_line_data[j+1])/3.0; 
    }

    information.low_pass_sino_lines_data[0] = information.mean_sino_line_data[0]; 
    information.low_pass_sino_lines_data[information.sinogram_adjusted_xdim-1] = information.mean_sino_line_data[information.sinogram_adjusted_xdim-1]; 
 
    // ring corrections 
    for (i=0;i<recon_info_record.sinogram_ydim;i++) {
      for (j=0;j<information.sinogram_adjusted_xdim;j++) { 
	      tmp = information.mean_sino_line_data[j]-information.low_pass_sino_lines_data[j]; 
	    
        if ((data[i*information.sinogram_adjusted_xdim+j] - (tmp * ring_coeff) ) > 0.0) 
	        data[i*information.sinogram_adjusted_xdim+j] -= (tmp * ring_coeff); 
	      else 
	        data[i*information.sinogram_adjusted_xdim+j] = 0.0; 
      
      } 
    } 
  }
} 
