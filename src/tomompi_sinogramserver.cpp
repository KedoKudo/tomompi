//_____________________________________________________________________________________ 
//_____________________________________________________________________________________ 
//_____________________________________________________________________________________ 
#include <iostream> 
#include <fstream> 
using namespace std; 
 
#include "tomompi.h" 

#include "mpi_communication_structures.h" 
#include "recon_algorithm.h" 
#include "filteredbackprojection.h" 
#include "gridrec.h" 
 
#include "centeringclass.h" 

//_____________________________________________________________________________________ 
#define DEFAULT_FILES_PER_PASS		128 
 
#define	TIME_FORMAT			"%04d-%02d-%02d %02d:%02d:%02d" 
//_____________________________________________________________________________________ 
//_____________________________________________________________________________________ 

int				front_end_id = 0, 
				sum_of_processes = 0, 
				num_files_processed = 0, 
				num_files_handled = 0, 
				num_sinograms, 
				data_xdim, 
				data_ydim, 
				sinogram_size, 
				num_white_fields, 
				num_dark_fields, 
				reconstruction_size, 
				reconstruction_count = 0, 
				files_per_pass, 
				sinograms_to_create, 
				pass_number, 
				sinogram_buffer; 

bool			sinogram_thread_running; 

char			exp_file_path[256], 
				exp_file_name[256], 
				data_group_index[256], 
				file_version[10]; 

FileListClass 		*top_projection_file_list = NULL, 
					*top_white_file_list      = NULL, 
					*top_dark_file_list       = NULL; 

NexusBoxClass           data_file; 

float		       	*dark_field_sino_ave_buf1    = NULL, 
					*dark_field_sino_ave_buf2    = NULL, 
					*dark_field_buffer_sinograms = NULL, 
					*dark_field_server_sinograms = NULL; 
 
unsigned short int	*temp_buffer                  = NULL, 
					*server_sinograms             = NULL, 
					*buffer_sinograms             = NULL, 
					*sino_buf1                    = NULL, 
					*sino_buf2                    = NULL, 
					*white_field_server_sinograms = NULL, 
					*white_field_buffer_sinograms = NULL, 
					*white_field_sino_buf1        = NULL, 
					*white_field_sino_buf2        = NULL;

SINOGRAMINFO	*sinogram_list = NULL; 
 
int				total_reconstruction_timer, 
				sinogram_pass_timer, 
				waiting_on_clients_timer, 
				waiting_on_MPI_timer; 

//_____________________________________________________________________________________ 
 
void ServerAcknowledgements (void) { 
    LogFileClass *acknowledge_file = NULL; 
 
	acknowledge_file = new LogFileClass ("./", "server_acknowledge.txt"); 
 
	acknowledge_file->Message ("__________________________________________________________________"); 
	acknowledge_file->Message ("TomoMPI_SinogramServer.cpp"); 
	acknowledge_file->Message (""); 
	acknowledge_file->Message ("Controlling code for TomoMPI cluster."); 
	acknowledge_file->Message ("Developed and maintained by:"); 
	acknowledge_file->Message ("       Brian Tieman"); 
	acknowledge_file->Message ("       Argonne National Laboratory"); 
	acknowledge_file->Message ("       tieman@aps.anl.gov"); 
    acknowledge_file->Message (""); 
	acknowledge_file->Message ("8/20/2003   V1.0   BT  First version with acknowledgements"); 
    acknowledge_file->Message ("9/4/2003    V1.0   BT  Server now generates integer sinograms"); 
    acknowledge_file->Message ("        to send to the clients.  The clients then perform the"); 
	acknowledge_file->Message ("        normalization."); 
	acknowledge_file->Message ("9/10/2003   V1.2   BT  Server now determines whether or not to autocenter"); 
	acknowledge_file->Message ("        based on a flag set in the experiment file."); 
	acknowledge_file->Message ("1/30/2004   V1.3   BT  Upgraded cluster to run on RedHat 9.0.  This version fails"); 
	acknowledge_file->Message ("        to run properly when linked with the MKL.  A few memory issues were fixed"); 
	acknowledge_file->Message ("        as a result of this upgrade."); 
	acknowledge_file->Message ("2/26/2004   V1.4   BT  Dark fields are now averaged before being sent to the client."); 
    acknowledge_file->Message ("        This allows for more flexibility in the acquisition--for example, it's now"); 
	acknowledge_file->Message ("        possible to take a number of darks at the end or take darks at regular"); 
	acknowledge_file->Message ("        intervals throughout the dat set and the same reconstruction code can"); 
    acknowledge_file->Message ("        handle it transparently."); 
    acknowledge_file->Message ("2/26/2004   V1.4   BT  Started adding ERROR lines to the log files.  This will"); 
	acknowledge_file->Message ("        hopefully make it easier to track down problems in the future. Grep all the"); 
	acknowledge_file->Message ("        log files for \"ERROR\" to find the error lines.  All error lines should"); 
	acknowledge_file->Message ("        reference the routine where they occured."); 
	acknowledge_file->Message ("3/15/2004   V1.5   BT  Fixed a few problems with the darkfield average.  One was that"); 
	acknowledge_file->Message ("        there were still a few places where dark_field_sino_ave was being interpreted as"); 
	acknowledge_file->Message ("        unsigned short instead of float.  This may have been causing some problems with"); 
    acknowledge_file->Message ("        data overruns.  Also, the wrong dark_fiel_sino_ave was being forwarded to the clients."); 
	acknowledge_file->Message ("3/15/2004   V1.5   BT  The client now recieves the centering related fields--"); 
	acknowledge_file->Message ("        fixed_shift and fixed_shift_value.  If fixed_shift is non-zero, the auto"); 
	acknowledge_file->Message ("        findcenter routine will be ignored and the value in fixed_shift_value will"); 
    acknowledge_file->Message ("        be used for the shift value for all sinograms.  For the moment, these values"); 
    acknowledge_file->Message ("        are read from reconstruction.cfg as the fields <Fixed Shift> and <Fixed Shift Value>"); 
    acknowledge_file->Message ("3/23/2004   V1.5   BT  The array used to store the average dark field was not being"); 
	acknowledge_file->Message ("        initialized prior to calculating a new average dark field.  This caused problems"); 
    acknowledge_file->Message ("        as the dark field average is calculated progressively."); 
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
	acknowledge_file->Message ("6/7/2004    V1.6  BT  The new debug mode allowed me to very quickly find a fairly serious"); 
	acknowledge_file->Message ("        bug.  The dark fields were not being handled properly.  Only the first row of dark"); 
	acknowledge_file->Message ("        field was ever being handed out.  In addition, the buffer the dark field average"); 
	acknowledge_file->Message ("        was being performed into was not being adaquately protected from multi-thread access."); 
	acknowledge_file->Message ("        In some cases, teh dark fields handed out may have been very wrong due to averages"); 
	acknowledge_file->Message ("        either not being complete, or the averaging buffer being set back to 0 prior to reading"); 
    acknowledge_file->Message ("        new buffers.  The result is that most sinograms were getting the wrong dark field "); 
    acknowledge_file->Message ("        information with occasional sinograms recieving dark fields with very wrong values."); 
	acknowledge_file->Message ("        Thi issue has now been resolved."); 
	acknowledge_file->Message ("12/10/2004	V1.8  BT Built tomompi on development cluster (AMD opteron machines).  Turned on"); 
	acknowledge_file->Message ("		compiler optimizations.  Also built with MPE profiling libraries so I could start"); 
	acknowledge_file->Message ("		trying to use jumpshot to profile our performance."); 
	acknowledge_file->Message ("12/10/2004	V1.8  BT Created a recon info record structure that contains all the initialization"); 
	acknowledge_file->Message ("		parameters needed by the clients.  This allowed for a significant reduction of MPI"); 
    acknowledge_file->Message ("		sends to the clients upon initialization as nearly the entire structure can be sent"); 
	acknowledge_file->Message ("		at once--as opposed to each datum being sent seperately.  Also did this for the sinogram"); 
	acknowledge_file->Message ("		data that was being sent each time a sinogram was requested.  This appears to have about"); 
	acknowledge_file->Message ("		a 10% impact on performance."); 
	acknowledge_file->Message ("1/24/2006	V1.9  BT Added a reconstruction.cfg flag <Ring Removal> which can be 0 or 1.  1 means"); 
	acknowledge_file->Message ("		use ring removal.  0 means do not use ring removal."); 
	acknowledge_file->Message ("1/31/2006	V1.10 BT Added a reconstruction.cfg parameter <Ring Removal Coeff> which can be 0= < coeff <= 5."); 
	acknowledge_file->Message ("		This adjusts a coefficient used in the ring removal code.  See the CenteringCalss comments for"); 
    acknowledge_file->Message ("		more specific information about what this does."); 
    acknowledge_file->Message ("3/1/2006	V1.10 BT Added a feature where random slices can be entered into a file called slices.list"); 
    acknowledge_file->Message ("		If you set the parameter <Use Slices File> in reconstruction.cfg to 1, the cluster will only"); 
	acknowledge_file->Message ("		reconstruct the files listed in slices.list.  Setting <Use Slices File> to 0 will cause the"); 
    acknowledge_file->Message ("		cluster to perform a full reconstruction as normal."); 
    acknowledge_file->Message ("3/1/2006	V1.10 BT Fixed a bug where the first row of a sinogram was actually the wrong row--it was"); 
    acknowledge_file->Message ("		the last row of the previous sinogram.  It turns out that we should generally ignore this last"); 
	acknowledge_file->Message ("		row as it is redundant data.  It's the 180 degree angle from the first row.  The code was"); 
	acknowledge_file->Message ("		modified to ignore the last file (ie sinograms are 720 pixels high rather than 721 pixels."); 
    acknowledge_file->Message ("3/23/2006	V1.11 BT Added a methode to compress reconstructed files.  Just using the built in loss-less"); 
	acknowledge_file->Message ("		compression in HDF provides poor compression results on the floating point reconstructed values."); 
	acknowledge_file->Message ("		However, as the reconstruction only occupies a circle in the image, it's possible to set the corners"); 
	acknowledge_file->Message ("		to 0.0 without loosing any data.  All these 0.0 values then compress rather well yielding an"); 
	acknowledge_file->Message ("		approximate 25% reduction in file size for the saved reconstruction files."); 
	acknowledge_file->Message ("		The appropriate reconstruction.cfg flag is <Copression Type> and its values may be:"); 
	acknowledge_file->Message ("			NONE -- no compression (corners will not even be set to 0.0"); 
	acknowledge_file->Message ("			LZW -- gzip type compression (this appears to be the best for most data sets"); 
	acknowledge_file->Message ("			RLE -- Run Length Encoding"); 
	acknowledge_file->Message ("			HUF -- Skipping-Huffman"); 
	acknowledge_file->Message ("3/23/2006	V1.11 BT While debugging the problem noted in the centeringclass.cpp on 3/22/2006 it"); 
	acknowledge_file->Message ("		was discovered that averaging the white fields can sometimes give a far superior result."); 
	acknowledge_file->Message ("		A flag was added that allows the user to switch from using a single white field per projection"); 
	acknowledge_file->Message ("		to averaging all the white fields and using the average white field for all projections.  The"); 
	acknowledge_file->Message ("		flag is <Average White Fields> and it should be set to 1 to average all the whites or to 0"); 
	acknowledge_file->Message ("		to use a single white field (last acquired white field before acquisition of projection)"); 
	acknowledge_file->Message ("		per projection"); 
	acknowledge_file->Message ("5/12/2006   V1.12 BT Fixed a bug that cropped up with performing reconstructions with projections"); 
	acknowledge_file->Message ("		whose Y dimension is not equal to a multiple of files_per_pass.  If the projections were"); 
	acknowledge_file->Message ("		not an even multiple of files_per_pass, only the full passes were being done, the \"extra\""); 
	acknowledge_file->Message ("		files that did not fill up a full pass were not being reconstructed.  It should now be possible"); 
	acknowledge_file->Message ("		to reconstruct samples of _ANY_ size."); 
	acknowledge_file->Message ("4/23/2006	V1.13 BT Added a new debug mode:"); 
	acknowledge_file->Message ("		No Reconstruction -- perform all the steps of a reconstruction, but, to save time, do not"); 
	acknowledge_file->Message ("		actually reconstruct the data."); 
	acknowledge_file->Message ("4/23/2006	V1.13 BT Removed MarkSuspicious from the possible debug modes--it wasn't useful."); 
	acknowledge_file->Message ("4/23/2006	V1.13 BT  Modified handling of filter type to be done by int rather than by string."); 
	acknowledge_file->Message ("6/15/2006	V1.13 BT  Modified the way parameters are passed into the reconstruction code as follows:"); 
	acknowledge_file->Message ("		reconstruction.cfg -- ths file only contains the data top directory, reconstruction save directory,"); 
	acknowledge_file->Message ("		and the experiment file to use.  All other configuration parameters were moved into the exp_file"); 
	acknowledge_file->Message ("		where they belong."); 
	acknowledge_file->Message ("		override_exp_file.config -- if this file exists, it is read in _AFTER_ the exp_file is read."); 
	acknowledge_file->Message ("		Parameters listed in the override_exp_file.config will over ride those found in the experiment"); 
	acknowledge_file->Message ("		file."); 
	acknowledge_file->Message ("		Parameters used for the reconstruction--whether read in from the exp_file or from the"); 
	acknowledge_file->Message ("		override file--are now rewritten to the experiment file plus a <sample name>.config file that"); 
	acknowledge_file->Message ("		is saved in the same directory as the experiment file.  The intent is that the exp_file is"); 
	acknowledge_file->Message ("		always up to dat with what was done last and there is a file that can be copied to the"); 
	acknowledge_file->Message ("		override file as a simple way to batch a number of samples with the same processing."); 
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
	acknowledge_file->Message ("4/13/2007	V2.1 BT fixed a typo on the reading of ring removal coefficient.  The code was doing atoi"); 
	acknowledge_file->Message ("		when it should have been doing atof."); 
	acknowledge_file->Message ("4/12/2007	V2.1 BT Added ability to scale data into short ints by using a user supplied data_range_min"); 
	acknowledge_file->Message ("		and data_range_max values.  These values must be supplied prior to reconstructing as the full"); 
	acknowledge_file->Message ("		data range of the reconstruction can not be calculated until completion.  The data_range_min/"); 
	acknowledge_file->Message ("		data_range_max values of a previous reconstruction are generally good values to use.  The actual"); 
	acknowledge_file->Message ("		data_range_min/max values are stored in the exp file as ;experiment;reconstruction;scaled_data_min"); 
	acknowledge_file->Message ("		and ;experiment;reconstruction;scaled_data_max for reference."); 
    acknowledge_file->Message ("5/17/2007   V2.2 BT Fixed some types where the Sinogram Pre-Norm and Post-Norm debug types did not work."); 
    acknowledge_file->Message ("6/05/2007   V2.2 BT Modified the way sinograms were being read in CreateSinograms to be more efficient."); 
    acknowledge_file->Message ("        This is still nowhere near as efficient as reading a large group of sequential lines at once."); 
    acknowledge_file->Message ("6/05/2007   V2.2 BT Fixed a problem with memory allocation for the sinogram buffers.  Turns out the"); 
    acknowledge_file->Message ("        buffers were being allocated for sinograms with 1 line less than needed (eg: 719 instead of 720)"); 
    acknowledge_file->Message ("        In general, this did not seem to cause problems, but when using the slice file, the system"); 
    acknowledge_file->Message ("        would segment fault or cause some other kind of memory error when freeing the sinogram buffer"); 
    acknowledge_file->Message ("        memory."); 
    acknowledge_file->Message ("6/22/2007   V2.2 BT Fixed a problem with the reconstruction file list files not having the correct"); 
    acknowledge_file->Message ("		extension when HDF5 was the selected output file type."); 
    acknowledge_file->Message ("6/22/2007	V2.2 BT Fixed a problem where the total reconstruction time was not being put into the"); 
    acknowledge_file->Message ("		experiment file correctly."); 
    acknowledge_file->Message ("11/14/2007	V2.3 BT Added a BIN file type to the output file types.  If a user selects BIN as the"); 
    acknowledge_file->Message ("		<output file type> the output files will be written as binary without any sort of header."); 
    acknowledge_file->Message ("		Eventually, we will need a bit of a header, but we haven't settled on one yet.  Bin files"); 
	acknowledge_file->Message ("		can not be compressed"); 
	acknowledge_file->Message ("11/14/2007	V2.3 BT Made some hacks to test the performance of reading binary projection files"); 
	acknowledge_file->Message ("		vs reading HDF projection files.  The results are phenomonal!  The sinogram creation time"); 
	acknowledge_file->Message ("		plummets by a factor of 10!  I was able to reconstruct a full 2048x2048x1440 data set in"); 
	acknowledge_file->Message ("		< 5 minutes!  A lot still needs to be done to clean this up for general use if we decide to"); 
	acknowledge_file->Message ("		abandon HDF as an intermediary file format, but the performance advantage is huge!"); 
	acknowledge_file->Message ("1/28/2008	V2.3 BT	Modified code to not read a reconstruction.cfg file but to get that information"); 
	acknowledge_file->Message ("		from the command line.  This will allow for the program to be queued in a system like N1"); 
	acknowledge_file->Message ("		Grid Engine.  As queued applications can be reprioritized by the system, having a single"); 
	acknowledge_file->Message ("		reconstruction.cfg file is problematic as it's impossible to keep the file in sync with"); 
	acknowledge_file->Message ("		the correct application in the queue--now each instance can get that info from the"); 
	acknowledge_file->Message ("		command prompt."); 
	acknowledge_file->Message ("1/28/2008	V2.3 BT Modified the location of the log files to be in the smaple tree.  Also allow"); 
	acknowledge_file->Message ("		for versioning of the log and override files.  This is done by creating a new log directory"); 
	acknowledge_file->Message ("		\"logs_##\" each time the sample is reconstructed.  This way, we can keep a running list"); 
	acknowledge_file->Message ("		of logs and configurations in case the sample is reconstructed multiple times."); 
	acknowledge_file->Message ("1/30/2008  V2.3 BT Made it so the logs_ directory is prepended with 0 as needed so ls orders properly."); 
	acknowledge_file->Message ("2/4/2008   V2.3 BT Added an error log file that will only be created if there is an error.  The"); 
	acknowledge_file->Message ("		error log will log all errros encountered by the system.  This provides a convenient location"); 
	acknowledge_file->Message ("		to look for problems if reconstructions fail or do not look quite right."); 
	acknowledge_file->Message ("2/4/2008   V2.3 BT The system should no longer crash due to a raw data file missing.  Instead, the"); 
	acknowledge_file->Message ("		system will mark in the error log that the file is missing and will grab the previous file"); 
	acknowledge_file->Message ("		of the same type to fill in.  If a projection is missing, the previous projection will"); 
	acknowledge_file->Message ("		be used to complete the reconstruction."); 
	acknowledge_file->Message ("2/11/2008  V2.3 BT CreateSinograms will now look to the file extension to decide if the file format"); 
    acknowledge_file->Message ("        is hdf or binary.  It expects .bin for binary and anything else defaults to HDF."); 
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
 
	acknowledge_file->Message ("__________________________________________________________________"); 
 
	delete (acknowledge_file); 
} 
 
//_____________________________________________________________________________________ 
 
void InitSampleLocation (int argc, char *argv[]) { 
	strcpy (exp_file_path, argv[1]); 
	if (exp_file_path[strlen (exp_file_path)-1] != '/') 
		strcat (exp_file_path, "/"); 
 
	strcpy (recon_info_record.reconstruction_path, exp_file_path); 
	strcat (recon_info_record.reconstruction_path, "reconstructed/"); 
	mkdir (recon_info_record.reconstruction_path, 0777); 
 
	strcpy (exp_file_name, argv[2]); 
} 
 
//_____________________________________________________________________________________ 
 
int ReadOverrideConfigFile (void) { 

	int 	error_code, 
			index; 

	FILE 	*setup_file = NULL; 

	char    override_file_name[256], 
			parameter_name[256], 
			parameter_value[256]; 

	char 	ch; 

	error_code = 0; 
 
	recon_info_record.recon_algorithm = RECONSTRUCTION_GRIDREC; 
	recon_info_record.filter = FILTER_NONE; 
 
	strcpy (override_file_name, exp_file_path); 
	strcat (override_file_name, "override_exp_file.config"); 
	if ((setup_file = fopen (override_file_name, "rt")) == NULL){ 
		sprintf(msg, "override_exp_file.config could not be opened.  We will be using the information from the experiment file."); 
		log_file->Message (msg); 
		return (1); 
	} 
 
	sprintf(msg, "override_exp_file.config exists!  We will be OVERRIDING the information from the experiment file."); 
	log_file->Message (msg); 
 
	strcpy (parameter_name, "empty"); 
	while (strcmp (parameter_name, "End")) { 
		//Get parameter_name 
		ch = fgetc (setup_file); 
		if (ch == EOF) 
			return(-1); 
 
		index = 0; 
		while (ch != '<') { 
			//ignore comment lines 
			if (ch == '#') 
				while ((ch != '\n') && (ch != '\r')) { 
					ch = fgetc (setup_file); 
					if (ch == EOF) 
						return (-1); 
				} 
 
			ch = fgetc (setup_file); 
			if (ch == EOF) 
				return(-1); 
		}

		ch = fgetc (setup_file); 
		if (ch == EOF) 
			return(-1); 
		while (ch != '>') { 
			parameter_name[index] = ch; 
			index++; 
 
			ch = fgetc (setup_file); 
			if (ch == EOF) 
				return(-1); 
		} 
		parameter_name[index] = '\0'; 
 
		//Get parameter_value 
		if (strcmp (parameter_name, "End")) { 
			ch = fgetc (setup_file); 
			if (ch == EOF) 
				return(-1); 
 
			index = 0; 
			while ((ch != '\n') && (ch != '\r')) { 
				parameter_value[index] = ch; 
				index++; 
 
				ch = fgetc (setup_file); 
				if (ch == EOF) 
					return(-1); 
			} 
			parameter_value[index] = '\0'; 
		} 
 
		//Reconstruction algorithm 
		if (!strcmp (parameter_name, "Reconstruction Algorithm")) { 
			if (!strcmp (parameter_value, "FBP Delay Only")) { 
				sprintf (msg, "Reconstruction algorithm: FBP Delay Only"); 
				recon_info_record.recon_algorithm = RECONSTRUCTION_FBP_DELAY_ONLY; 
			} 
			if (!strcmp (parameter_value, "GRIDREC Delay Only")) { 
				sprintf (msg, "Reconstruction algorithm: Gridrec Delay Only"); 
				recon_info_record.recon_algorithm = RECONSTRUCTION_GRIDREC_DELAY_ONLY; 
			} 
			if (!strcmp (parameter_value, "FBP No Optimization")) { 
				sprintf (msg, "Reconstruction algorithm: FBP No Optimization"); 
				recon_info_record.recon_algorithm = RECONSTRUCTION_FBP_NO_OPTIMIZATION; 
			} 
			if (!strcmp (parameter_value, "FBP Optimized")) { 
				sprintf (msg, "Reconstruction algorithm: FBP Optimized"); 
				recon_info_record.recon_algorithm = RECONSTRUCTION_FBP_OPTIMIZED; 
			} 
			if (!strcmp (parameter_value, "FBP Cylinder Only")) { 
				sprintf (msg, "Reconstruction algorithm: FBP Cylinder Only"); 
				recon_info_record.recon_algorithm = RECONSTRUCTION_FBP_CYLINDER_ONLY; 
			} 
			if (!strcmp (parameter_value, "GRIDREC")) { 
				sprintf (msg, "Reconstruction algorithm: GRIDREC"); 
				recon_info_record.recon_algorithm = RECONSTRUCTION_GRIDREC; 
			} 

		} 

		// Zero padding for GRIDREC
		if (!strcmp (parameter_name, "GRIDREC Padding")) {

			if (!strcmp (parameter_value, "GRIDREC_PADDING_NONE")) { 
				sprintf (msg, "GRIDREC Padding: GRIDREC_PADDING_NONE"); 
				recon_info_record.gridrec_padding = GRIDREC_PADDING_NONE; 
			} 

			if (!strcmp (parameter_value, "GRIDREC_PADDING_HALF")) { 
				sprintf (msg, "GRIDREC Padding: GRIDREC_PADDING_HALF"); 
				recon_info_record.gridrec_padding = GRIDREC_PADDING_HALF; 
			} 

			if (!strcmp (parameter_value, "GRIDREC_PADDING_ONE_AND_HALF")) { 
				sprintf (msg, "GRIDREC Padding: GRIDREC_PADDING_ONE_AND_HALF"); 
				recon_info_record.gridrec_padding = GRIDREC_PADDING_ONE_AND_HALF; 
			} 

			if (!strcmp (parameter_value, "GRIDREC_PADDING_BOUNDARY")) { 
				sprintf (msg, "GRIDREC Padding: GRIDREC_PADDING_BOUNDARY"); 
				recon_info_record.gridrec_padding = GRIDREC_PADDING_BOUNDARY; 
			} 

		 } 

		//Filter algorithm 
		if (!strcmp (parameter_name, "Filter Algorithm")) { 

			if (!strcmp (parameter_value, "None")) { 
				sprintf (msg, "Filter Algorithm: None"); 
				recon_info_record.filter = FILTER_NONE; 
			} 

			if (!strcmp (parameter_value, "Shepp/Logan")) { 
				sprintf (msg, "Filter Algorithm: Shepp/Logan"); 
				recon_info_record.filter = FILTER_SHEPP_LOGAN; 
			}

			if (!strcmp (parameter_value, "Hann")) { 
				sprintf (msg, "Filter Algorithm: Hann"); 
				recon_info_record.filter = FILTER_HANN; 
			}

			if (!strcmp (parameter_value, "Hamming")) { 
				sprintf (msg, "Filter Algorithm: Hamming"); 
				recon_info_record.filter = FILTER_HAMMING; 
			} 

			if (!strcmp (parameter_value, "Ramp")) { 
				sprintf (msg, "Filter Algorithm: Ramp"); 
				recon_info_record.filter = FILTER_RAMP; 
			} 

			if (!strcmp (parameter_value, "FBP")) { 
				sprintf (msg, "Filter Algorithm: FBP"); 
				recon_info_record.filter = FILTER_FBP; 
			} 
 
		} 
 
		//Output File Format 
		if (!strcmp (parameter_name, "Output File Format")) {

			if (!strcmp (parameter_value, "HDF4")) { 
				sprintf (msg, "Output File Format: HDF4"); 
				recon_info_record.file_format = HDF4; 
			}

			if (!strcmp (parameter_value, "HDF5")) { 
				sprintf (msg, "Output File Format: HDF5"); 
				recon_info_record.file_format = HDF5; 
			}

			if (!strcmp (parameter_value, "BIN")) { 
				sprintf (msg, "Output File Format: BIN"); 
				recon_info_record.file_format = BIN; 
			}

		} 
 
		//Use Slices File 
		if (!strcmp (parameter_name, "Use Slices File")) { 
			
			if (!strcmp (parameter_value, "true")) { 
				sprintf (msg, "Use Slices File: true"); 
				recon_info_record.use_slices_file = true; 
			} 

			if (!strcmp (parameter_value, "false")) { 
				sprintf (msg, "Use Slices File: false"); 
				recon_info_record.use_slices_file = false; 
			}

		} 
 
		//Files per Pass 
		if (!strcmp (parameter_name, "Files Per Pass")) { 
			sprintf (msg, "Files Per Pass: %s", parameter_value); 
			files_per_pass = atoi (parameter_value); 
		} 
 
		//Fixed Shift 
		if (!strcmp (parameter_name, "Use Fixed Shift")) { 
			
			if (!strcmp (parameter_value, "true")) { 
				sprintf (msg, "Use Fixed Shift: true"); 
				recon_info_record.use_fixed_shift = true; 
			} 

			if (!strcmp (parameter_value, "false")) { 
				sprintf (msg, "Use Fixed Shift: false"); 
				recon_info_record.use_fixed_shift = false; 
			} 
		} 
		if (!strcmp (parameter_name, "Fixed Shift Value")) { 
			sprintf (msg, "Fixed Shift Value: %s", parameter_value); 
			recon_info_record.fixed_shift_value = atof (parameter_value); 
		} 
 
		//Start fixed shift 
		if (!strcmp (parameter_name, "Start Fixed Shift")) { 
			sprintf (msg, "Start Fixed Shift: %s", parameter_value); 
			recon_info_record.start_fixed_shift = atof (parameter_value);  
		} 
 
		//End fixed shift 
		if (!strcmp (parameter_name, "End Fixed Shift")) { 
			sprintf (msg, "End Fixed Shift: %s", parameter_value); 
			recon_info_record.end_fixed_shift = atof (parameter_value);   
		} 

		//Fixed shift interval
		if (!strcmp (parameter_name, "Fixed Shift Interval")) { 
			sprintf (msg, "Fixed Shift Interval: %s", parameter_value); 
			recon_info_record.fixed_shift_interval = atof (parameter_value); 
		} 

 
		//Auto-Centering 
		if (!strcmp (parameter_name, "Auto-Centering")) { 
			if (!strcmp (parameter_value, "true")) { 
				sprintf (msg, "Auto-Centering: true"); 
				recon_info_record.centering = true; 
			} 
			
			if (!strcmp (parameter_value, "false")) { 
				sprintf (msg, "Auto-Centering: false"); 
				recon_info_record.centering = false; 
			} 
		} 
 
		//Ring Removal 
		if (!strcmp (parameter_name, "Use Ring Removal")) { 

			if (!strcmp (parameter_value, "true")) { 
				sprintf (msg, "Use Ring Removal: true"); 
				recon_info_record.use_ring_removal = true; 
			} 

			if (!strcmp (parameter_value, "false")) { 
				sprintf (msg, "Use Ring Removal: false"); 
				recon_info_record.use_ring_removal = false; 
			} 
		} 
		if (!strcmp (parameter_name, "Ring Removal Coefficient")) { 
			sprintf (msg, "Ring Removal Coefficient: %s", parameter_value); 
			recon_info_record.ring_removal_coeff = atof (parameter_value); 
		} 
 
		//Average White Fields 
		if (!strcmp (parameter_name, "Average White Fields")) { 
			
			if (!strcmp (parameter_value, "true")) { 
				sprintf (msg, "Average White Fields: true"); 
				recon_info_record.average_white_fields = true; 
			} 

			if (!strcmp (parameter_value, "false")) { 
				sprintf (msg, "Average White Fields: false"); 
				recon_info_record.average_white_fields = false; 
			}

		} 
 
		//Data Index 
		if (!strcmp (parameter_name, "Data Index")) { 
			sprintf (msg, "Data Index: %s", parameter_value); 
			strcpy (data_group_index, parameter_value); 
		} 
 
		//Compression Type 
		if (!strcmp (parameter_name, "Compression Type")) { 
			
			if (!strcmp (parameter_value, "None")) { 
				sprintf (msg, "Compression Type: None"); 
				recon_info_record.compression_type = NX_COMP_NONE; 
			} 

			if (!strcmp (parameter_value, "LZW")) { 
				sprintf (msg, "Compression Type: LZW"); 
				recon_info_record.compression_type = NX_COMP_LZW; 
			}

			if (!strcmp (parameter_value, "RLE")) { 
				sprintf (msg, "Compression Type: RLE"); 
				recon_info_record.compression_type = NX_COMP_RLE; 
			}

			if (!strcmp (parameter_value, "HUF")) { 
				sprintf (msg, "Compression Type: HUF"); 
				recon_info_record.compression_type = NX_COMP_HUF; 
			}

		} 
 
		//Rescale To Integer 
		if (!strcmp (parameter_name, "Rescale To Int")) { 
			
			if (!strcmp (parameter_value, "true")) { 
				sprintf (msg, "Rescale To Int: true"); 
				recon_info_record.rescale_to_int = true; 
			}

			if (!strcmp (parameter_value, "false")) { 
				sprintf (msg, "Rescale To Int: false"); 
				recon_info_record.rescale_to_int = false; 
			}

		}

		if (!strcmp (parameter_name, "Data Range Min")) { 
			sprintf (msg, "Data Range Min: %s", parameter_value); 
			recon_info_record.scale_data_range_min = atof (parameter_value); 
		}

		if (!strcmp (parameter_name, "Data Range Max")) { 
			sprintf (msg, "Data Range Max: %s", parameter_value); 
			recon_info_record.scale_data_range_max = atof (parameter_value); 
		} 
 
		//Debug 
		if (!strcmp (parameter_name, "Debug")) { 

			if (!strcmp (parameter_value, "None")) { 
				sprintf (msg, "Debug: None"); 
				recon_info_record.debug = DEBUG_NONE; 
			}

			if (!strcmp (parameter_value, "White Field")) { 
				sprintf (msg, "White Field"); 
				recon_info_record.debug = DEBUG_WHITEFIELD; 
			}

			if (!strcmp (parameter_value, "Dark Field")) { 
				sprintf (msg, "Debug: Dark Field"); 
				recon_info_record.debug = DEBUG_DARKFIELD; 
			}

			if (!strcmp (parameter_value, "Sinogram Pre-Norm")) { 
				sprintf (msg, "Debug: Sinogram Pre-Norm"); 
				recon_info_record.debug = DEBUG_PRENORM; 
			}

			if (!strcmp (parameter_value, "Sinogram Post-Norm")) { 
				sprintf (msg, "Debug: Sinogram Post-Norm"); 
				recon_info_record.debug = DEBUG_POSTNORM; 
			} 

			if (!strcmp (parameter_value, "Post Centering")) { 
				sprintf (msg, "Debug: Post Centering"); 
				recon_info_record.debug = DEBUG_POSTCENTERING; 
			} 

			if (!strcmp (parameter_value, "Post Ring Removal")) { 
				sprintf (msg, "Debug: Post Ring Removal"); 
				recon_info_record.debug = DEBUG_POSTRING; 
			} 

			if (!strcmp (parameter_value, "Full")) { 
				sprintf (msg, "Debug: Full"); 
				recon_info_record.debug = DEBUG_FULL; 
			} 

			if (!strcmp (parameter_value, "No Output")) { 
				sprintf (msg, "Debug: No Output"); 
				recon_info_record.debug = DEBUG_NO_OUTPUT; 
			}

			if (!strcmp (parameter_value, "No Reconstruction")) { 
				sprintf (msg, "Debug: No Reconstruction"); 
				recon_info_record.debug = DEBUG_NO_RECONSTRUCTION; 
			} 

		} 
 
		log_file->Message(msg); 
	} 
 
	fclose (setup_file); 
 
	//Validate Ring Removal Coefficient 
	if (recon_info_record.ring_removal_coeff < 0) { 
		recon_info_record.ring_removal_coeff = 0; 
		sprintf (msg, "WARNING! Ring Removal Coefficient must be > 0--I will adjust it for you."); 
		log_file->Message (msg); 
	} 
 
	if (recon_info_record.ring_removal_coeff > 5) { 
		recon_info_record.ring_removal_coeff = 5; 
		sprintf (msg, "WARNING! Ring Removal Coefficient must be < 5--I will adjust it for you."); 
		log_file->Message (msg); 
	} 
 
	return (error_code); 
 
} 
 
//_____________________________________________________________________________________ 
 
void WriteConfigFile () { 
	FILE 	*setup_file = NULL; 
	char	*configuration_file = NULL; 
 
	configuration_file = (char *) malloc (strlen (recon_info_record.log_path) + strlen (recon_info_record.base_name) + strlen (".config") + 1); 
    if (configuration_file == NULL) { 
        sprintf (msg, "Could not allocat memory for configuration_file."); 
        error_log->addError (msg, "WriteConfigFile ()"); 
    } 
 
	sprintf (configuration_file, "%s%s.config", recon_info_record.log_path, recon_info_record.base_name); 
 
	if ((setup_file = fopen (configuration_file, "w")) == NULL) { 
		sprintf(msg, "%s%s", configuration_file, " could not be opened."); 
		log_file->ErrorMessage (msg, "WriteConfigFile"); 
        error_log->addError (msg, "WriteConfigFile ()"); 
        error_log->addAutoResolution ("Configuration file was not written!"); 
		return; 
	} 
 
	fprintf (setup_file, "##%s\n", exp_file_path); 
	fprintf (setup_file, "##%s\n", exp_file_name); 
	fprintf (setup_file, "##TomoMPI Version:%s\n", VERSION); 
 
	//Reconstruction algorithm 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##Reconstruction Algorithm\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		FBP Delay Only\n"); 
	fprintf (setup_file, "##		GRIDREC Delay Only\n"); 
	fprintf (setup_file, "##		FBP No Optimization\n"); 
	fprintf (setup_file, "##		FBP Optimized\n"); 
	fprintf (setup_file, "##		FBP Cylinder Only\n"); 
	fprintf (setup_file, "##		GRIDREC\n"); 
	switch (recon_info_record.recon_algorithm) { 

		case RECONSTRUCTION_FBP_DELAY_ONLY      : fprintf (setup_file, "<Reconstruction Algorithm>FBP Delay Only\n"     ); break; 
		case RECONSTRUCTION_GRIDREC_DELAY_ONLY  : fprintf (setup_file, "<Reconstruction Algorithm>GRIDREC Delay Only\n" ); break; 
		case RECONSTRUCTION_FBP_NO_OPTIMIZATION : fprintf (setup_file, "<Reconstruction Algorithm>FBP No Optimization\n"); break; 
		case RECONSTRUCTION_FBP_OPTIMIZED       : fprintf (setup_file, "<Reconstruction Algorithm>FBP Optimized\n"      ); break; 
		case RECONSTRUCTION_FBP_CYLINDER_ONLY   : fprintf (setup_file, "<Reconstruction Algorithm>FBP Cylinder Only\n"  ); break; 
		case RECONSTRUCTION_GRIDREC             : fprintf (setup_file, "<Reconstruction Algorithm>GRIDREC\n"            ); break; 
		
		default : fprintf (setup_file, "<Reconstruction Algorithm>Gridrec\n"); break; 
	} 

	//GRIDREC Padding

	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##GRIDREC Padding\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		GRIDREC_PADDING_NONE\n"); 
	fprintf (setup_file, "##		GRIDREC_PADDING_HALF\n"); 
	fprintf (setup_file, "##		GRIDREC_PADDING_ONE_AND_HALF\n"); 
	fprintf (setup_file, "##		GRIDREC_PADDING_BOUNDARY\n"); 
	switch (recon_info_record.gridrec_padding) { 

		case GRIDREC_PADDING_NONE         : fprintf (setup_file, "<GRIDREC Padding>GRIDREC_PADDING_NONE\n"        ); break; 
		case GRIDREC_PADDING_HALF         : fprintf (setup_file, "<GRIDREC Padding>GRIDREC_PADDING_HALF\n"        ); break; 
		case GRIDREC_PADDING_ONE_AND_HALF : fprintf (setup_file, "<GRIDREC Padding>GRIDREC_PADDING_ONE_AND_HALF\n"); break; 
		case GRIDREC_PADDING_BOUNDARY     : fprintf (setup_file, "<GRIDREC Padding>GRIDREC_PADDING_BOUNDARY\n"    ); break; 
		
		default : fprintf (setup_file, "<GRIDREC Padding>GRIDREC_PADDING_BOUNDARY\n"); break; 
	 } 

	//Filter algorithm 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##Filter Algorithm\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##	None\n"); 
	fprintf (setup_file, "##	Shepp/Logan\n"); 
	fprintf (setup_file, "##	Hann\n"); 
	fprintf (setup_file, "##	Hamming\n"); 
	fprintf (setup_file, "##	Ramp\n"); 
	fprintf (setup_file, "##	FBP\n"); 
	switch (recon_info_record.filter) {

		case FILTER_NONE        : fprintf (setup_file, "<Filter Algorithm>None\n"       ); break; 
		case FILTER_SHEPP_LOGAN : fprintf (setup_file, "<Filter Algorithm>Shepp/Logan\n"); break; 
		case FILTER_HANN        : fprintf (setup_file, "<Filter Algorithm>Hann\n"       ); break; 
		case FILTER_HAMMING     : fprintf (setup_file, "<Filter Algorithm>Hamming\n"    ); break; 
		case FILTER_RAMP        : fprintf (setup_file, "<Filter Algorithm>Ramp\n"       ); break; 
		case FILTER_FBP         : fprintf (setup_file, "<Filter Algorithm>FBP\n"        ); break; 

		default : fprintf (setup_file, "<Filter Algorithm>None\n"); break; 
	 } 
 
	//Output File Format 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##Output File Format\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		HDF4\n"); 
	fprintf (setup_file, "##		HDF5\n"); 
	fprintf (setup_file, "##		BIN\n"); 
	switch (recon_info_record.file_format) { 

		case HDF4 : fprintf (setup_file, "<Output File Format>HDF4\n"); break; 
		case HDF5 : fprintf (setup_file, "<Output File Format>HDF5\n"); break; 
		case BIN  : fprintf (setup_file, "<Output File Format>BIN\n" ); break; 

		default : fprintf (setup_file, "<Output File Format>HDF5\n"); break; 
	} 
 
	//Use Slices File 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##Use Slices File\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		true\n"); 
	fprintf (setup_file, "##		false\n"); 
	fprintf (setup_file, "##	Requires:\n"); 
	fprintf (setup_file, "##		<Start Fixed Shift>int\n"); 
	fprintf (setup_file, "##		<End Fixed Shift>int\n"); 
	if (recon_info_record.use_slices_file) 
		fprintf (setup_file, "<Use Slices File>true\n"); 
	else 
		fprintf (setup_file, "<Use Slices File>false\n"); 
 
	//Start fixed shift 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##Start Fixed Shift\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		start fixed shift value<int>\n"); 
	fprintf (setup_file, "##	Requires:\n"); 
	fprintf (setup_file, "##		<Use Slices File>true\n"); 
	fprintf (setup_file, "##		<Use Fixed Shift>true\n"); 
	fprintf (setup_file, "<Start Fixed Shift>%f\n", recon_info_record.start_fixed_shift); 
 
	//End fixed shift 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##End Fixed Shift\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		end fixed shift value<int>\n"); 
	fprintf (setup_file, "##	Requires:\n"); 
	fprintf (setup_file, "##		<Use Slices File>true\n"); 
	fprintf (setup_file, "##		<Use Fixed Shift>true\n"); 
	fprintf (setup_file, "<End Fixed Shift>%f\n", recon_info_record.end_fixed_shift); 

	//Fixed shift interval                     
	fprintf (setup_file, "\n");                             
	fprintf (setup_file, "##Fixed Shift Interval\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		fixed shift interval value<int>\n"); 
	fprintf (setup_file, "##	Requires:\n"); 
	fprintf (setup_file, "##		<Use Slices File>true\n"); 
	fprintf (setup_file, "##		<Use Fixed Shift>true\n"); 
	fprintf (setup_file, "<Fixed Shift Interval>%f\n", recon_info_record.fixed_shift_interval); 
  
	//Files per Pass 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##Files Per Pass\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		power of 2 prefered\n"); 
	fprintf (setup_file, "<Files Per Pass>%d\n", files_per_pass); 
 
	//Fixed Shift 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##Use Fixed Shift\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		true\n"); 
	fprintf (setup_file, "##		false\n"); 
	fprintf (setup_file, "##	Requires:\n"); 
	fprintf (setup_file, "##		<Slices File>false\n"); 
	fprintf (setup_file, "##		<Fixed Shift Value>int\n"); 
	if (recon_info_record.use_fixed_shift) 
		fprintf (setup_file, "<Use Fixed Shift>true\n"); 
	else 
		fprintf (setup_file, "<Use Fixed Shift>false\n"); 
	fprintf (setup_file, "##Fixed Shift Value\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		fixed shift value<int>\n"); 
	fprintf (setup_file, "##	Requires:\n"); 
	fprintf (setup_file, "##		<Use Fixed Shift>true\n"); 

	// fprintf (setup_file, "<Fixed Shift Value>%d\n", recon_info_record.fixed_shift_value); 
	fprintf (setup_file, "<Fixed Shift Value>%f\n", recon_info_record.fixed_shift_value); 
 
	//Auto-Centering 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##Auto-Centering\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		true\n"); 
	fprintf (setup_file, "##		false\n"); 
	if (recon_info_record.centering) 
		fprintf (setup_file, "<Auto-Centering>true\n"); 
	else 
		fprintf (setup_file, "<Auto-Centering>false\n"); 
 
	//Ring Removal 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##Use Ring Removal\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		true\n"); 
	fprintf (setup_file, "##		false\n"); 
	fprintf (setup_file, "##	Requires:\n"); 
	fprintf (setup_file, "##		<Ring Removal Coefficient>float\n"); 
	if (recon_info_record.use_ring_removal) 
		fprintf (setup_file, "<Use Ring Removal>true\n"); 
	else 
		fprintf (setup_file, "<Use Ring Removal>false\n"); 
	fprintf (setup_file, "##Ring Removal Coefficient\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		coefficient value<float>\n"); 
	fprintf (setup_file, "##	Requires:\n"); 
	fprintf (setup_file, "##		<Use Ring Removal>true\n"); 
	fprintf (setup_file, "<Ring Removal Coefficient>%f\n", recon_info_record.ring_removal_coeff); 
 
	//Average White Fields 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##Average White Fields\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		true\n"); 
	fprintf (setup_file, "##		false\n"); 
	if (recon_info_record.average_white_fields) 
		fprintf (setup_file, "<Average White Fields>true\n"); 
	else 
		fprintf (setup_file, "<Average White Fields>false\n"); 
 
	//Data Index 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##Data Index\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		example: ;entry1;data;data\n"); 
	fprintf (setup_file, "<Data Index>%s\n", data_group_index); 
 
	//Compression Type 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##Compressiong Type\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		NONE\n"); 
	fprintf (setup_file, "##		LZW\n"); 
	fprintf (setup_file, "##		RLE\n"); 
	fprintf (setup_file, "##		HUF\n"); 
	fprintf (setup_file, "##	Requires:\n"); 
	fprintf (setup_file, "##		<Output File Format>HDF4\n"); 
	fprintf (setup_file, "##		<Output File Format>HDF5\n"); 
	switch (recon_info_record.compression_type) { 

		case NX_COMP_NONE : fprintf (setup_file, "<Compression Type>None\n"); break; 
		case NX_COMP_LZW  : fprintf (setup_file, "<Compression Type>LZW\n" ); break; 
		case NX_COMP_RLE  : fprintf (setup_file, "<Compression Type>RLE\n" ); break; 
		case NX_COMP_HUF  : fprintf (setup_file, "<Compression Type>HUF\n" ); break; 

		default : fprintf (setup_file, "<Compression Type>None\n"); break; 
	} 
 
//Rescale To Int 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##Rescale To Int\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		true\n"); 
	fprintf (setup_file, "##		false\n"); 
	fprintf (setup_file, "##	Requires:\n"); 
	fprintf (setup_file, "##		<Data Range Min>float\n"); 
	fprintf (setup_file, "##		<Data Range Max>float\n"); 
	if (recon_info_record.rescale_to_int) 
		fprintf (setup_file, "<Rescale To Int>true\n"); 
	else 
		fprintf (setup_file, "<Rescale To Int>false\n"); 
	fprintf (setup_file, "##Data Range Min\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		minimum value to scale from<float>\n"); 
	fprintf (setup_file, "##	Requires:\n"); 
	fprintf (setup_file, "##		<Rescale To Int>true\n"); 
	fprintf (setup_file, "<Data Range Min>%f\n", recon_info_record.scale_data_range_min); 
	fprintf (setup_file, "##Data Range Max\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		maximum value to scale from<float>\n"); 
	fprintf (setup_file, "##	Requires:\n"); 
	fprintf (setup_file, "##		<Rescale To Int>true\n"); 
	fprintf (setup_file, "<Data Range Max>%f\n", recon_info_record.scale_data_range_max); 
 
//Debug 
	fprintf (setup_file, "\n"); 
	fprintf (setup_file, "##Debug\n"); 
	fprintf (setup_file, "##	Options:\n"); 
	fprintf (setup_file, "##		None\n"); 
	fprintf (setup_file, "##		White Field\n"); 
	fprintf (setup_file, "##		Dark Field\n"); 
	fprintf (setup_file, "##		Sinogram Pre-Norm\n"); 
	fprintf (setup_file, "##		Sinogram Post-Norm\n"); 
	fprintf (setup_file, "##		Post Centering\n"); 
	fprintf (setup_file, "##		Post Ring Removal\n"); 
	fprintf (setup_file, "##		Full\n"); 
	fprintf (setup_file, "##		No Output\n"); 
	fprintf (setup_file, "##		No Reconstruction\n"); 
	switch (recon_info_record.debug) { 

		case DEBUG_NONE              : fprintf (setup_file, "<Debug>None\n"              ); break; 
		case DEBUG_WHITEFIELD        : fprintf (setup_file, "<Debug>White Field\n"       ); break; 
		case DEBUG_DARKFIELD         : fprintf (setup_file, "<Debug>Dark Field\n"        ); break; 
		case DEBUG_PRENORM           : fprintf (setup_file, "<Debug>Sinogram Pre-Norm\n" ); break; 
		case DEBUG_POSTNORM          : fprintf (setup_file, "<Debug>Sinogram Post-Norm\n"); break; 
		case DEBUG_POSTCENTERING     : fprintf (setup_file, "<Debug>Post Centering\n"    ); break; 
		case DEBUG_POSTRING          : fprintf (setup_file, "<Debug>Post Ring Removal\n" ); break; 
		case DEBUG_FULL              : fprintf (setup_file, "<Debug>Full\n"              ); break; 
		case DEBUG_NO_OUTPUT         : fprintf (setup_file, "<Debug>No Output\n"         ); break; 
		case DEBUG_NO_RECONSTRUCTION : fprintf (setup_file, "<Debug>No Reconstruction\n" ); break; 

		default : fprintf (setup_file, "<Debug>None\n"); break; 
	} 
 
	fprintf (setup_file, "\n\n\n"); 
	fprintf (setup_file, "<End>\n"); 
 
	fclose (setup_file); 
 
	delete (configuration_file); 
 
	log_file->Message ("Done writing configuration file to sample directory."); 
 
} 
 
//_____________________________________________________________________________________ 
 
void ReadExpFile (void) { 
	NexusBoxClass	exp_file; 

	int	rank, 
		dims[2], 
		type, 
		loop; 

	char	index[256], 
			temp_str[256], 
			file_name[256], 
			*file_list = NULL, 
			*temp_file_list = NULL, 
			recon_algorithm[50], 
			error_msg[256]; 

	float	start_angle, 
			end_angle, 
			angle_interval; 
	
	file_list = NULL; 
 
	exp_file.SetReadScheme (ENTIRE_CONTENTS); 
	exp_file.ReadAll (exp_file_path, exp_file_name); 
 
	//File Version 
	strcpy (index, ";tomompi_version;tomompi_version"); 
	if (!exp_file.IndexExists(index)) { 
		strcpy (file_version, "V1.3"); 
	} 
	else { 
		exp_file.GetDatumInfo (index, &rank, dims, &type); 
		exp_file.GetDatum (index, file_version); 
	} 
 
	//Base Name 
	strcpy (index, ";experiment;setup;sample;base_name"); 
	if (!exp_file.IndexExists(index)) { 
		sprintf (error_msg, "Index %s does not exist.", index); 
		log_file->ErrorMessage (error_msg, "ReadExpFile"); 
        error_log->addError (error_msg, "ReadExpFile ()"); 
	} 
	exp_file.GetDatumInfo (index, &rank, dims, &type); 
	exp_file.GetDatum (index, recon_info_record.base_name); 
 
	//Reconstruction Algorithm 
	strcpy (index, ";experiment;reconstruction;algorithm"); 
	if (exp_file.IndexExists(index)) { 
		exp_file.GetDatumInfo (index, &rank, dims, &type); 
		exp_file.GetDatum (index, recon_algorithm); 
		if ((!strcmp (recon_algorithm, "FBP_DELAY_ONLY")) || (!strcmp (recon_algorithm, "FBP Delay Only"))) 
			recon_info_record.recon_algorithm = RECONSTRUCTION_FBP_DELAY_ONLY; 
		if ((!strcmp (recon_algorithm, "GRIDREC_DELAY_ONLY")) || (!strcmp (recon_algorithm, "GRIDREC Delay Only"))) 
			recon_info_record.recon_algorithm = RECONSTRUCTION_GRIDREC_DELAY_ONLY; 
		if ((!strcmp (recon_algorithm, "FBP1")) || (!strcmp (recon_algorithm, "FBP No Optimization"))) 
			recon_info_record.recon_algorithm = RECONSTRUCTION_FBP_NO_OPTIMIZATION; 
		if ((!strcmp (recon_algorithm, "FBP2")) || (!strcmp (recon_algorithm, "FBP Optimized"))) 
			recon_info_record.recon_algorithm = RECONSTRUCTION_FBP_OPTIMIZED; 
		if ((!strcmp (recon_algorithm, "FBP3")) || (!strcmp (recon_algorithm, "FBP Cylinder Only"))) 
			recon_info_record.recon_algorithm = RECONSTRUCTION_FBP_CYLINDER_ONLY; 
		if (!strcmp (recon_algorithm, "GRIDREC")) 
			recon_info_record.recon_algorithm = RECONSTRUCTION_GRIDREC; 
	} 
	else { 
		sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
		log_file->WarningMessage (error_msg, "ReadExpFile"); 
	} 
 
	//Filter Algorithm 
	strcpy (index, ";experiment;reconstruction;filter"); 
	if (!exp_file.IndexExists(index)) { 
		exp_file.GetDatumInfo (index, &rank, dims, &type); 
		exp_file.GetDatum (index, temp_str); 
		if (!strcmp (temp_str, "None")) 
			recon_info_record.filter = FILTER_NONE; 
		if ((!strcmp (temp_str, "Shepp_Logan")) || (!strcmp (temp_str, "Shepp/Logan"))) 
			recon_info_record.filter = FILTER_SHEPP_LOGAN; 
		if (!strcmp (temp_str, "Hann")) 
			recon_info_record.filter = FILTER_HANN; 
		if (!strcmp (temp_str, "Hamming")) 
			recon_info_record.filter = FILTER_HAMMING; 
		if (!strcmp (temp_str, "Ramp")) 
			recon_info_record.filter = FILTER_RAMP; 
		if (!strcmp (temp_str, "FBP")) 
			recon_info_record.filter = FILTER_FBP; 
	} 
	else { 
		sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
		log_file->WarningMessage (error_msg, "ReadExpFile"); 
	} 
 
	//Use Slices File 
	strcpy (index, ";experiment;reconstruction;cluster_config;use_slices_file"); 
	if (exp_file.IndexExists(index)) { 
		exp_file.GetDatumInfo (index, &rank, dims, &type); 
		exp_file.GetDatum (index, temp_str); 
		if (!strcmp (temp_str, "true")) 
			recon_info_record.use_slices_file = 1; 
		else 
			recon_info_record.use_slices_file = 0; 
	} 
	else { 
		sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
		log_file->WarningMessage (error_msg, "ReadExpFile"); 
	} 
 
	//Files per Pass 
	sprintf (error_msg, "Index %s: %d", index, files_per_pass); 
	log_file->WarningMessage (error_msg, "ReadExpFile"); 
 
	strcpy (index, ";experiment;reconstruction;cluster_config;files_per_pass"); 
	if (exp_file.IndexExists(index)) { 
		exp_file.GetDatumInfo (index, &rank, dims, &type); 
		exp_file.GetDatum (index, &files_per_pass); 
	} 
	else { 
		sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
		log_file->WarningMessage (error_msg, "ReadExpFile"); 
	} 
 
	sprintf (error_msg, "Index %s: %d", index, files_per_pass); 
	log_file->WarningMessage (error_msg, "ReadExpFile"); 
 
	//Fixed Shift 
	strcpy (index, ";experiment;reconstruction;use_fixed_shift"); 
	if (exp_file.IndexExists(index)) { 
		exp_file.GetDatumInfo (index, &rank, dims, &type); 
		exp_file.GetDatum (index, temp_str); 
		if (!strcmp (temp_str, "true")) { 
			recon_info_record.use_fixed_shift = 1; 
 
			strcpy (index, ";experiment;reconstruction;fixed_shift_value"); 
			if (exp_file.IndexExists(index)) { 
				exp_file.GetDatumInfo (index, &rank, dims, &type); 
				exp_file.GetDatum (index, &recon_info_record.fixed_shift_value); 
			} 
			else { 
				sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
				log_file->WarningMessage (error_msg, "ReadExpFile"); 
			} 
		} 
		else 
			recon_info_record.use_fixed_shift = 0; 
	} 
	else { 
		sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
		log_file->WarningMessage (error_msg, "ReadExpFile"); 
	} 
 
	//Auto-Centering 
	strcpy (index, ";experiment;reconstruction;auto_centering"); 
	if (exp_file.IndexExists(index)) { 
		exp_file.GetDatumInfo (index, &rank, dims, &type); 
		exp_file.GetDatum (index, temp_str); 
		if ((!strcmp (temp_str, "on")) || (!strcmp (temp_str, "true"))) 
			recon_info_record.centering = 1; 
		else 
			recon_info_record.centering = 0; 
	} 
	else { 
		sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
		log_file->WarningMessage (error_msg, "ReadExpFile"); 
	} 
 
	//Ring Removal 
	strcpy (index, ";experiment;reconstruction;use_ring_removal"); 
	if (exp_file.IndexExists(index)) { 
		exp_file.GetDatumInfo (index, &rank, dims, &type); 
		exp_file.GetDatum (index, temp_str); 
		if (!strcmp (temp_str, "true")) { 
			recon_info_record.use_ring_removal = 1; 
 
			strcpy (index, ";experiment;reconstruction;ring_removal_coefficient"); 
			if (exp_file.IndexExists(index)) { 
				exp_file.GetDatumInfo (index, &rank, dims, &type); 
				exp_file.GetDatum (index, &recon_info_record.ring_removal_coeff); 
			} 
			else { 
				sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
				log_file->WarningMessage (error_msg, "ReadExpFile"); 
			} 
		} 
		else 
			recon_info_record.use_ring_removal = 0; 
	} 
	else { 
		sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
		log_file->WarningMessage (error_msg, "ReadExpFile"); 
	} 
 
	//Average White Fields 
	strcpy (index, ";experiment;reconstruction;average_white_fields"); 
	if (exp_file.IndexExists(index)) { 
		exp_file.GetDatumInfo (index, &rank, dims, &type); 
		exp_file.GetDatum (index, temp_str); 
		if (!strcmp (temp_str, "true")) 
			recon_info_record.average_white_fields = 1; 
		else 
			recon_info_record.centering = 0; 
	} 
	else { 
		sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
		log_file->WarningMessage (error_msg, "ReadExpFile"); 
	} 
 
	//Data Group Index 
	strcpy (index, ";experiment;reconstruction;cluster_config;data_group_index"); 
	if (exp_file.IndexExists(index)) { 
		exp_file.GetDatumInfo (index, &rank, dims, &type); 
		exp_file.GetDatum (index, data_group_index); 
	} 
	else { 
		sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
		log_file->WarningMessage (error_msg, "ReadExpFile"); 
	} 
 
	//Generate Thetas 
	strcpy (index, ";experiment;acquisition;parameters;start_angle"); 
	if (!exp_file.IndexExists(index)) { 
		sprintf (error_msg, "Index %s does not exist.", index); 
		log_file->ErrorMessage (error_msg, "ReadExpFile"); 
		error_log->addError (error_msg, "ReadExpFile ()"); 
	} 
	exp_file.GetDatum (index, &start_angle); 
	sprintf (msg, "%s%f", "Start Angle: ", start_angle); 
	log_file->Message(msg); 
 
	strcpy (index, ";experiment;acquisition;parameters;end_angle"); 
	if (!exp_file.IndexExists(index)) { 
		sprintf (error_msg, "Index %s does not exist.", index); 
		log_file->ErrorMessage (error_msg, "ReadExpFile"); 
        error_log->addError (error_msg, "ReadExpFile ()"); 
	} 
	exp_file.GetDatum (index, &end_angle); 
	sprintf (msg, "%s%f", "End Angle: ", end_angle); 
	log_file->Message(msg); 
 
	strcpy (index, ";experiment;acquisition;parameters;angle_interval"); 
	if (!exp_file.IndexExists(index)) { 
		sprintf (error_msg, "Index %s does not exist.", index); 
		log_file->ErrorMessage (error_msg, "ReadExpFile"); 
        error_log->addError (error_msg, "ReadExpFile ()"); 
	} 
	exp_file.GetDatum (index, &angle_interval); 
	sprintf (msg, "%s%f", "Angle Interval: ", angle_interval); 
	log_file->Message(msg); 
 
	//In file version "V1.3" and lower (no file version listed in file), the correct index is 
	//      ";experiment;acquisition;parameters;whitedark_interval" 
	//In file version "V1.4" and higher (version saved in ";tomompi_version", the correct index is 
	//      ";experiment;acquisition;parameters;white_interval" 
	if (!strcmp (file_version, "V1.3")) 
		strcpy (index, ";experiment;acquisition;parameters;whitedark_interval"); 
	else 
		strcpy (index, ";experiment;acquisition;parameters;white_interval"); 
	if (!exp_file.IndexExists(index)) { 
		sprintf (error_msg, "Index %s does not exist.", index); 
		log_file->ErrorMessage (error_msg, "ReadExpFile"); 
        error_log->addError (error_msg, "ReadExpFile ()"); 
	} 
	exp_file.GetDatum (index, &recon_info_record.whitedark_interval); 
	sprintf (msg, "%s%d", "White/Dark Interval: ", recon_info_record.whitedark_interval); 
	log_file->Message(msg); 
 
	if (start_angle < end_angle) { 
		recon_info_record.theta_list_size = (int) ((end_angle - start_angle) / angle_interval); 
        if (recon_info_record.theta_list != NULL) { 
            free (recon_info_record.theta_list); 
            recon_info_record.theta_list = NULL; 
        } 
		recon_info_record.theta_list = (float *) malloc (sizeof(float)*recon_info_record.theta_list_size); 
        if (recon_info_record.theta_list == NULL) { 
            sprintf (msg, "Could not allocat memory for recon_info_record.theta_list."); 
            error_log->addError (msg, "ReadExpFile ()"); 
        } 
        memset (recon_info_record.theta_list, 0, sizeof(float)*recon_info_record.theta_list_size); 
 
		int theta_index = 0; 
		while (start_angle < end_angle) { 
	    	recon_info_record.theta_list[theta_index] = start_angle; 
	    	start_angle += angle_interval; 
	    	theta_index++; 
		} 
	} 
	else { 
		recon_info_record.theta_list_size = (int) ((start_angle - end_angle) / angle_interval); 
		if (recon_info_record.theta_list != NULL) { 
			free (recon_info_record.theta_list); 
			recon_info_record.theta_list = NULL;
			} 
	  	recon_info_record.theta_list = (float *) malloc (sizeof(float)*recon_info_record.theta_list_size); 
	  	if (recon_info_record.theta_list == NULL) {
			sprintf (msg, "Could not allocat memory for recon_info_record.theta_list"); 
	      	error_log->addError (msg, "ReadExpFile ()"); 
		} 
	  	memset (recon_info_record.theta_list, 0, sizeof(float)*recon_info_record.theta_list_size); 
 
	  	int theta_index = 0; 
	  	while (end_angle < start_angle) { 
	    	recon_info_record.theta_list[theta_index] = end_angle; 
	      	end_angle += angle_interval; 
	      	theta_index++; 
	    } 
 	} 
 
	sprintf (msg, "List of %d theta values generated...", recon_info_record.theta_list_size); 
	log_file->Message(msg); 
 
	//Projection File Names 
	strcpy (index, ";experiment;acquisition;projections;names"); 
	if (!exp_file.IndexExists(index)) { 
		sprintf (error_msg, "Index %s does not exist.", index); 
		log_file->ErrorMessage (error_msg, "ReadExpFile"); 
		error_log->addError (error_msg, "ReadExpFile ()"); 
	} 
	exp_file.GetDatumInfo (index, &rank, dims, &type); 
 
	sprintf (msg, "Dims[0]=%d, dims[1]=%d", dims[0], dims[1]); 
	log_file->Message (msg); 
 
    if (file_list != NULL) { 
        free (file_list); 
        file_list = NULL; 
    } 
    file_list = (char *) malloc (sizeof(char)*dims[0]*dims[1]); 
    if (file_list == NULL) { 
        sprintf (msg, "Could not allocat memory for file_list."); 
        error_log->addError (msg, "ReadExpFile ()"); 
    } 
	exp_file.GetDatum (index, file_list); 
 
	temp_file_list = file_list; 
 
	strncpy (file_name, temp_file_list, dims[1]); 
	file_name[dims[1]] = '\0'; 
	temp_file_list = &temp_file_list[dims[1]]; 
 
	top_projection_file_list = new FileListClass(file_name); 
 
	for (loop=1;loop<dims[0];loop++) { 
		strncpy (file_name, temp_file_list, dims[1]); 
		file_name[dims[1]] = '\0'; 
		temp_file_list = &temp_file_list[dims[1]]; 
 
		top_projection_file_list->AddToList (new FileListClass(file_name)); 
	} 
 
	sprintf (msg, "List of %d projection file names read...", dims[0]); 
	log_file->Message(msg); 
 
	//White Field File Names 
	strcpy (index, ";experiment;acquisition;white_field;names"); 
	if (!exp_file.IndexExists(index)) { 
		sprintf (error_msg, "Index %s does not exist.", index); 
		log_file->ErrorMessage (error_msg, "ReadExpFile"); 
        error_log->addError (error_msg, "ReadExpFile ()"); 
	} 
	exp_file.GetDatumInfo (index, &rank, dims, &type); 
    if (file_list != NULL) { 
        free (file_list); 
        file_list = NULL; 
    } 
	file_list = (char *) malloc (sizeof(char)*dims[0]*dims[1]); 
    if (file_list == NULL) { 
        sprintf (msg, "Could not allocat memory for file_list."); 
        error_log->addError (msg, "ReadExpFile ()"); 
    } 
	exp_file.GetDatum (index, file_list); 
 
	temp_file_list = file_list; 
 
	strncpy (file_name, temp_file_list, dims[1]); 
	file_name[dims[1]] = '\0'; 
	temp_file_list = &temp_file_list[dims[1]]; 
 
	top_white_file_list = new FileListClass(file_name); 
 
	for (loop=1;loop<dims[0];loop++) { 
		strncpy (file_name, temp_file_list, dims[1]); 
		file_name[dims[1]] = '\0'; 
		temp_file_list = &temp_file_list[dims[1]]; 
 
		top_white_file_list->AddToList (new FileListClass(file_name)); 
	} 
 
	num_white_fields = dims[0]; 
 
	sprintf (msg, "List of %d white field file names read...", num_white_fields); 
	log_file->Message(msg); 
 
	//Dark Field File Names 
	strcpy (index, ";experiment;acquisition;black_field;names"); 
	if (!exp_file.IndexExists(index)) { 
		sprintf (error_msg, "Index %s does not exist.", index); 
		log_file->ErrorMessage (error_msg, "ReadExpFile"); 
		error_log->addError (error_msg, "ReadExpFile ()"); 
	} 
	exp_file.GetDatumInfo (index, &rank, dims, &type); 
	if (file_list != NULL) { 
	    free (file_list); 
	    file_list = NULL; 
	} 
	file_list = (char *) malloc (sizeof(char)*dims[0]*dims[1]); 
	if (file_list == NULL) { 
	    sprintf (msg, "Could not allocat memory for file_list."); 
	    error_log->addError (msg, "ReadExpFile ()"); 
	} 
	exp_file.GetDatum (index, file_list); 
 
	temp_file_list = file_list; 
 
	strncpy (file_name, temp_file_list, dims[1]); 
	file_name[dims[1]] = '\0'; 
	temp_file_list = &temp_file_list[dims[1]]; 
 
	top_dark_file_list = new FileListClass(file_name); 
 
	num_dark_fields = dims[0]; 
 
	for (loop=1;loop<num_dark_fields;loop++) { 
		strncpy (file_name, temp_file_list, dims[1]); 
		file_name[dims[1]] = '\0'; 
		temp_file_list = &temp_file_list[dims[1]]; 
 
		top_dark_file_list->AddToList (new FileListClass(file_name)); 
	} 
 
	if (file_list != NULL) { 
	    free (file_list); 
	    file_list == NULL; 
	} 
 
	sprintf (msg, "List of %d dark field file names read...", num_dark_fields); 
	log_file->Message(msg); 
 
	//Compression Type 
	strcpy (index, ";experiment;reconstruction;cluster_config;compression_type"); 
	if (exp_file.IndexExists(index)) {

		exp_file.GetDatumInfo (index, &rank, dims, &type); 
		exp_file.GetDatum (index, temp_str); 

		if (!strcmp (temp_str, "None")) 
			recon_info_record.compression_type = NX_COMP_NONE; 
		if (!strcmp (temp_str, "LZW")) 
			recon_info_record.compression_type = NX_COMP_RLE; 
		if (!strcmp (temp_str, "RLE")) 
			recon_info_record.compression_type = NX_COMP_RLE; 
		if (!strcmp (temp_str, "HUF")) 
			recon_info_record.compression_type = NX_COMP_HUF; 
	} 
	else { 
		sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
		log_file->WarningMessage (error_msg, "ReadExpFile"); 
	} 
 
	//Debug Type 
	strcpy (index, ";experiment;reconstruction;cluster_config;debug_type"); 
	if (exp_file.IndexExists(index)) { 

		exp_file.GetDatumInfo (index, &rank, dims, &type); 
		exp_file.GetDatum (index, temp_str); 

		if (!strcmp (temp_str, "None")) 
			recon_info_record.debug = DEBUG_NONE; 
		if (!strcmp (temp_str, "White Field")) 
			recon_info_record.debug = DEBUG_WHITEFIELD; 
		if (!strcmp (temp_str, "Dark Field")) 
			recon_info_record.debug = DEBUG_DARKFIELD; 
		if (!strcmp (temp_str, "Sinogram Pre-Norm")) 
			recon_info_record.debug = DEBUG_PRENORM; 
		if (!strcmp (temp_str, "Sinogram Post-Norm")) 
			recon_info_record.debug = DEBUG_POSTNORM; 
		if (!strcmp (temp_str, "Post Centering")) 
			recon_info_record.debug = DEBUG_POSTCENTERING; 
		if (!strcmp (temp_str, "Post Ring Removal")) 
			recon_info_record.debug = DEBUG_POSTRING; 
		if (!strcmp (temp_str, "Full")) 
			recon_info_record.debug = DEBUG_FULL; 
		if (!strcmp (temp_str, "No Output")) 
			recon_info_record.debug = DEBUG_NO_OUTPUT; 
		if (!strcmp (temp_str, "No Reconstruction")) 
			recon_info_record.debug = DEBUG_NO_RECONSTRUCTION; 
	} 
	else { 
		sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
		log_file->WarningMessage (error_msg, "ReadExpFile"); 
	} 
 
	//Read Data Dims 
	strcpy (index, ";experiment;setup;detector;size_x"); 
	if (exp_file.IndexExists(index)) {

		exp_file.GetDatumInfo (index, &rank, dims, &type); 
		exp_file.GetDatum (index, &data_xdim); 
	} 
	else { 
		sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
		log_file->WarningMessage (error_msg, "ReadExpFile"); 
	}

	strcpy (index, ";experiment;setup;detector;size_y"); 
	if (exp_file.IndexExists(index)) { 
		exp_file.GetDatumInfo (index, &rank, dims, &type); 
		exp_file.GetDatum (index, &data_ydim); 
	} 
	else { 
		sprintf (error_msg, "Index %s does not exist.  Using default.", index); 
		log_file->WarningMessage (error_msg, "ReadExpFile"); 
	} 
 
} 
 
//_____________________________________________________________________________________ 
 
void UpdateExpFile (void) { 
	int			        file_name_size; 
	NexusBoxClass		exp_file; 
	int    				dims[2], 
						loop; 
	char    			file_name[256], 
						value_str[256], 
						*file_list = NULL; 
	float		    	time; 
 
	log_file->Message ("Updating experiment file..."); 
 
	exp_file.SetFileMode (HDF4_MODE); 
	exp_file.SetReadScheme (ENTIRE_CONTENTS); 
	exp_file.ReadAll (exp_file_path, exp_file_name); 
 
	loop = 0; 
 
	switch (recon_info_record.file_format) { 
		case HDF4 : { 
			sprintf (file_name, "rec_%s_%05d.hdf", recon_info_record.base_name, loop); 
			break; 
		} 
		case HDF5 : { 
			sprintf (file_name, "rec_%s_%05d.h5", recon_info_record.base_name, loop); 
			break; 
		} 
		default : { 
			sprintf (file_name, "rec_%s_%05d.hdf", recon_info_record.base_name, loop); 
			break; 
		} 
	} 
 
	file_name_size = strlen (file_name); 
    if (file_list != NULL) { 
        free (file_list); 
        file_list = NULL; 
    } 
	file_list = (char *) malloc (num_files_handled * file_name_size); 
    if (file_list == NULL) { 
        sprintf (msg, "Could not allocat memory for file_list."); 
        error_log->addError (msg, "UpdateExpFile ()"); 
    } 
 
	for (loop=0;loop<num_files_handled;loop++) { 
		switch (recon_info_record.file_format) { 
			case HDF4 : { 
				sprintf (file_name, "rec_%s_%05d.hdf", recon_info_record.base_name, sinogram_list[loop].sinogram_number); 
				break; 
			} 
			case HDF5 : { 
				sprintf (file_name, "rec_%s_%05d.h5", recon_info_record.base_name, sinogram_list[loop].sinogram_number); 
				break; 
			} 
			default : { 
				sprintf (file_name, "rec_%s_%05d.hdf", recon_info_record.base_name, sinogram_list[loop].sinogram_number); 
				break; 
			} 
		} 
		strncpy (&file_list[loop*file_name_size], file_name, file_name_size); 
	} 
 
	//If the "reconstruction" group does not exist--create it... 
	if (!exp_file.IndexExists (";experiment;reconstruction")) 
		exp_file.CreateGroup (";experiment", "reconstruction", "Reconstruction"); 
	
	//If the "cluster_config" group does not exist--create it... 
	if (!exp_file.IndexExists (";experiment;reconstruction;cluster_config")) 
		exp_file.CreateGroup (";experiment;reconstruction", "cluster_config", "Cluster"); 
 
	//Write reconstructed list 
	dims[1] = file_name_size; 
	dims[0] = num_files_handled; 
	if (!exp_file.IndexExists (";experiment;reconstruction;slices")) 
		exp_file.CreateGroup (";experiment;reconstruction", "slices", "FileList"); 
	if (!exp_file.IndexExists (";experiment;reconstruction;slices;names")) 
		exp_file.CreateField (";experiment;reconstruction;slices", "names", 2, dims, NX_CHAR, file_list); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;slices;names", file_list, 2, dims, NX_CHAR); 
 
	//Write reconstruction time 
	dims[1] = 0; 
	dims[0] = 1; 
	time = (float) log_file->GetTime(total_reconstruction_timer); 
	if (!exp_file.IndexExists (";experiment;reconstruction;reconstruction_time")) 
		exp_file.CreateField (";experiment;reconstruction", "reconstruction_time", 1, dims, NX_FLOAT32, (void *) &time); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;reconstruction_time", (void *) &time, 1, dims, NX_FLOAT32); 
 
	//Write min and max data range values 
	dims[1] = 0; 
	dims[0] = 1; 
	if (!exp_file.IndexExists (";experiment;reconstruction;data_range_min")) 
		exp_file.CreateField (";experiment;reconstruction", "data_range_min", 1, dims, NX_FLOAT32, (void *) &data_range_min); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;data_range_min", (void *) &data_range_min, 1, dims, NX_FLOAT32); 
 
	if (!exp_file.IndexExists (";experiment;reconstruction;data_range_max")) 
		exp_file.CreateField (";experiment;reconstruction", "data_range_max", 1, dims, NX_FLOAT32, (void *) &data_range_max); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;data_range_max", (void *) &data_range_max, 1, dims, NX_FLOAT32); 
 
	//Write scaled min and max data range values and if scaling was even used 
	if (recon_info_record.rescale_to_int) 
		strcpy (value_str, "true"); 
	else 
		strcpy (value_str, "false"); 
	dims[1] = 0; 
	dims[0] = strlen (value_str); 
	if (!exp_file.IndexExists (";experiment;reconstruction;rescaled_to_int")) 
		exp_file.CreateField (";experiment;reconstruction", "rescaled_to_int", 1, dims, NX_CHAR, value_str); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;rescaled_to_int", value_str, 1, dims, NX_CHAR); 
 
	dims[1] = 0; 
	dims[0] = 1; 
	if (!exp_file.IndexExists (";experiment;reconstruction;scale_data_min")) 
		exp_file.CreateField (";experiment;reconstruction", "scale_data_min", 1, dims, NX_FLOAT32, (void *) &recon_info_record.scale_data_range_min); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;scale_data_min", (void *) &recon_info_record.scale_data_range_min, 1, dims, NX_FLOAT32); 
 
	if (!exp_file.IndexExists (";experiment;reconstruction;scale_data_max")) 
		exp_file.CreateField (";experiment;reconstruction", "scale_data_max", 1, dims, NX_FLOAT32, (void *) &recon_info_record.scale_data_range_max); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;scale_data_max", (void *) &recon_info_record.scale_data_range_max, 1, dims, NX_FLOAT32); 
 
	//Write number of processes used 
	if (!exp_file.IndexExists (";experiment;reconstruction;cluster_config;num_processes")) 
		exp_file.CreateField (";experiment;reconstruction;cluster_config", "num_processes", 1, dims, NX_UINT32, &num_processes); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;cluster_config;num_processes", &num_processes, 1, dims, NX_UINT32); 
 
	//Write TomoMPI Version 
	strcpy (value_str, VERSION); 
	dims[1] = 0; 
	dims[0] = strlen (value_str); 
	if (!exp_file.IndexExists (";experiment;reconstruction;cluster_config;version")) 
		exp_file.CreateField (";experiment;reconstruction;cluster_config", "version", 1, dims, NX_CHAR, value_str); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;cluster_config;version", value_str, 1, dims, NX_CHAR); 
 
	//Write reconstruction algorithm 
	switch (recon_info_record.recon_algorithm) { 

		case RECONSTRUCTION_FBP_DELAY_ONLY      : sprintf (value_str, "FBP Delay Only"     ); break; 
		case RECONSTRUCTION_GRIDREC_DELAY_ONLY  : sprintf (value_str, "GRIDREC Delay Only" ); break; 
		case RECONSTRUCTION_FBP_NO_OPTIMIZATION : sprintf (value_str, "FBP No Optimization"); break; 
		case RECONSTRUCTION_FBP_OPTIMIZED       : sprintf (value_str, "FBP Optimized"      ); break; 
		case RECONSTRUCTION_FBP_CYLINDER_ONLY   : sprintf (value_str, "FBP Cylinder Only"  ); break; 
		case RECONSTRUCTION_GRIDREC             : sprintf (value_str, "GRIDREC"            ); break;

		default : sprintf (value_str, "Unknown"); break; 
	}

	dims[1] = 0; 
	dims[0] = strlen (value_str); 
	if (!exp_file.IndexExists (";experiment;reconstruction;algorithm")) 
		exp_file.CreateField (";experiment;reconstruction", "algorithm", 1, dims, NX_CHAR, value_str); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;algorithm", value_str, 1, dims, NX_CHAR); 
 
	//Write filter 
	switch (recon_info_record.filter) {

		case FILTER_NONE        : sprintf (value_str, "None"       ); break; 
		case FILTER_SHEPP_LOGAN : sprintf (value_str, "Shepp/Logan"); break; 
		case FILTER_HANN        : sprintf (value_str, "Hann"       ); break; 
		case FILTER_HAMMING     : sprintf (value_str, "Hamming"    ); break; 
		case FILTER_RAMP        : sprintf (value_str, "Ramp"       ); break; 
		case FILTER_FBP         : sprintf (value_str, "FBP"        ); break; 

		default : sprintf (value_str, "Unknown"); break; 
	} 

	dims[1] = 0; 
	dims[0] = strlen (value_str); 
	if (!exp_file.IndexExists (";experiment;reconstruction;filter")) 
		exp_file.CreateField (";experiment;reconstruction", "filter", 1, dims, NX_CHAR, value_str); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;filter", value_str, 1, dims, NX_CHAR); 
 
	//Write use_slices_file 
	if (recon_info_record.use_slices_file) 
		strcpy (value_str, "true"); 
	else 
		strcpy (value_str, "false"); 
	dims[1] = 0; 
	dims[0] = strlen (value_str); 
	if (!exp_file.IndexExists (";experiment;reconstruction;cluster_config;use_slices_file")) 
		exp_file.CreateField (";experiment;reconstruction;cluster_config", "use_slices_file", 1, dims, NX_CHAR, value_str); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;cluster_config;use_slices_file", value_str, 1, dims, NX_CHAR); 
 
	//Files per Pass 
	dims[1] = 0; 
	dims[0] = 1; 
	if (!exp_file.IndexExists (";experiment;reconstruction;cluster_config;files_per_pass")) 
		exp_file.CreateField (";experiment;reconstruction;cluster_config", "files_per_pass", 1, dims, NX_INT32, (void *) &files_per_pass); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;cluster_config;files_per_pass", (void *) &files_per_pass, 1, dims, NX_INT32); 
 
	//Fixed Shift 
	if (recon_info_record.use_fixed_shift) 
		strcpy (value_str, "true"); 
	else 
		strcpy (value_str, "false"); 

	dims[1] = 0; 
	dims[0] = strlen (value_str); 
	if (!exp_file.IndexExists (";experiment;reconstruction;use_fixed_shift")) 
		exp_file.CreateField (";experiment;reconstruction", "use_fixed_shift", 1, dims, NX_CHAR, value_str); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;use_fixed_shift", value_str, 1, dims, NX_CHAR); 

	dims[1] = 0; 
	dims[0] = 1; 
	if (!exp_file.IndexExists (";experiment;reconstruction;fixed_shift_value"))
	  exp_file.CreateField (";experiment;reconstruction", "fixed_shift_value", 1, dims, NX_FLOAT32, (void *) &recon_info_record.fixed_shift_value); 
	else
	  exp_file.PutDatum (";experiment;reconstruction;fixed_shift_value", (void *) &recon_info_record.fixed_shift_value, 1, dims, NX_FLOAT32); 
	
 
	//Auto-Centering 
	if (recon_info_record.centering) 
		strcpy (value_str, "true"); 
	else 
		strcpy (value_str, "false"); 
	
	dims[1] = 0; 
	dims[0] = strlen (value_str); 
	
	if (!exp_file.IndexExists (";experiment;reconstruction;auto_centering")) 
		exp_file.CreateField (";experiment;reconstruction", "auto_centering", 1, dims, NX_CHAR, value_str); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;auto_centering", value_str, 1, dims, NX_CHAR); 
 
	//Ring Removal 
	if (recon_info_record.use_ring_removal) 
		strcpy (value_str, "true"); 
	else 
		strcpy (value_str, "false"); 
	
	dims[1] = 0; 
	dims[0] = strlen (value_str); 
	
	if (!exp_file.IndexExists (";experiment;reconstruction;use_ring_removal")) 
		exp_file.CreateField (";experiment;reconstruction", "use_ring_removal", 1, dims, NX_CHAR, value_str); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;use_ring_removal", value_str, 1, dims, NX_CHAR); 
	
	dims[1] = 0; 
	dims[0] = 1; 
	
	if (!exp_file.IndexExists (";experiment;reconstruction;ring_removal_coefficient")) 
		exp_file.CreateField (";experiment;reconstruction", "ring_removal_coefficient", 1, dims, NX_FLOAT32, (void *) &recon_info_record.ring_removal_coeff); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;ring_removal_coefficient", (void *) &recon_info_record.ring_removal_coeff, 1, dims, NX_FLOAT32); 
 
	//Average White Fields 
	if (recon_info_record.average_white_fields) 
		strcpy (value_str, "true"); 
	else 
		strcpy (value_str, "false"); 
	
	dims[1] = 0; 
	dims[0] = strlen (value_str); 
	
	if (!exp_file.IndexExists (";experiment;reconstruction;average_white_fields")) 
		exp_file.CreateField (";experiment;reconstruction", "average_white_fields", 1, dims, NX_CHAR, value_str); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;average_white_fields", value_str, 1, dims, NX_CHAR); 
 
	//Data Index 
	strcpy (value_str, data_group_index); 
	dims[1] = 0; 
	dims[0] = strlen (value_str); 
	
	if (!exp_file.IndexExists (";experiment;reconstruction;cluster_config;data_group_index")) 
		exp_file.CreateField (";experiment;reconstruction;cluster_config", "data_group_index", 1, dims, NX_CHAR, value_str); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;cluster_config;data_group_index", value_str, 1, dims, NX_CHAR); 
 
	//Compression Type 
	switch (recon_info_record.compression_type) {

		case NX_COMP_NONE : sprintf (value_str, "None"); break; 
		case NX_COMP_LZW  : sprintf (value_str, "LZW" ); break; 
		case NX_COMP_RLE  : sprintf (value_str, "RLE" ); break; 
		case NX_COMP_HUF  : sprintf (value_str, "HUF" ); break;

		default : sprintf (value_str, "Unknown"); break; 
	}

	dims[1] = 0; 
	dims[0] = strlen (value_str); 
	if (!exp_file.IndexExists (";experiment;reconstruction;cluster_config;compression_type")) 
		exp_file.CreateField (";experiment;reconstruction;cluster_config", "compression_type", 1, dims, NX_CHAR, value_str); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;cluster_config;compression_type", value_str, 1, dims, NX_CHAR); 
 
	//Debug Type 
	switch (recon_info_record.debug) {

		case DEBUG_NONE              : sprintf (value_str, "None"              ); break; 
		case DEBUG_WHITEFIELD        : sprintf (value_str, "White Field"       ); break; 
		case DEBUG_DARKFIELD         : sprintf (value_str, "Dark Field"        ); break; 
		case DEBUG_PRENORM           : sprintf (value_str, "Sinogram Pre-Norm" ); break; 
		case DEBUG_POSTNORM          : sprintf (value_str, "Sinogram Post-Norm"); break; 
		case DEBUG_POSTCENTERING     : sprintf (value_str, "Post Centering"    ); break; 
		case DEBUG_POSTRING          : sprintf (value_str, "Post Ring Removal" ); break; 
		case DEBUG_FULL              : sprintf (value_str, "Full"              ); break; 
		case DEBUG_NO_OUTPUT         : sprintf (value_str, "No Output"         ); break; 
		case DEBUG_NO_RECONSTRUCTION : sprintf (value_str, "No Reconstruction" ); break; 
		
		default : sprintf (value_str, "Unknown"); break; 
	} 

	dims[1] = 0; 
	dims[0] = strlen (value_str);

	if (!exp_file.IndexExists (";experiment;reconstruction;cluster_config;debug_type")) 
		exp_file.CreateField (";experiment;reconstruction;cluster_config", "debug_type", 1, dims, NX_CHAR, value_str); 
	else 
		exp_file.PutDatum (";experiment;reconstruction;cluster_config;debug_type", value_str, 1, dims, NX_CHAR); 
 
	exp_file.WriteAll (exp_file_path, exp_file_name); 
 
    if (file_list != NULL) { 
	   free (file_list); 
       file_list = NULL; 
    } 
 
	log_file->Message ("Experiment file updated..."); 
} 
 
//_____________________________________________________________________________________ 
 
void CreateSinogramsFromBinary (void) { 

	FileListClass		*current_file = NULL; 
	int					current_line, 
						data_size, 
						data_offset, 
						loop, 
						header_size; 
	char    			temp_path[256], 
						file_name[256]; 
	ifstream			input_file; 
	struct stat 		stat_buffer; 
 
	sprintf (msg, "\nLooks like files are binary.  Starting pass %d", pass_number); 
	log_file->Message(msg); 
	sprintf (msg, "Starting to create %d sinograms", files_per_pass); 
	log_file->TimeStamp (msg); 
 
	log_file->StartTimer (sinogram_pass_timer); 
 
	strcpy (temp_path, exp_file_path); 
	strcat (temp_path, "raw/"); 
 
    log_file->Message ("Generating first pass of sinograms..."); 
 
    if (recon_info_record.use_slices_file) { 
        log_file->Message ("Generating requested sinograms only..."); 
         
        //buffer projections 
        header_size = 1024; 
         
        current_file = top_projection_file_list; 
        current_line = 0; 
        while (current_file->NextInList () != NULL) { 
			sprintf (file_name, "%s%s", temp_path, current_file->file_name);

			if (stat (file_name, &stat_buffer) == 0) { 
				input_file.open (file_name, ios::binary); 
	
				data_size = sizeof (unsigned short) * recon_info_record.sinogram_xdim; 
				
				for (loop=0;loop<files_per_pass;loop++) { 
					data_offset = header_size + sizeof (unsigned short) * sinogram_list[loop].sinogram_number * recon_info_record.sinogram_xdim; 
					input_file.seekg (data_offset, ios::beg); 
					input_file.read ((char *) &buffer_sinograms[(loop*sinogram_size)+(current_line*recon_info_record.sinogram_xdim)], data_size); 
				} 
	
				input_file.close (); 
			} 
			else { 
				for (loop=0;loop<sinograms_to_create;loop++) 
					memcpy (&buffer_sinograms[(loop*sinogram_size)+(current_line*recon_info_record.sinogram_xdim)], 
					        &buffer_sinograms[(loop*sinogram_size)+((current_line-1)*recon_info_record.sinogram_xdim)], 
							sizeof (unsigned short int) * recon_info_record.sinogram_xdim
					); 
	
				sprintf (msg, "Could not read file %s!", file_name); 
				error_log->addError (msg, "StartTomoMPIServer ()"); 
				error_log->addAutoResolution ("Using previous projection instead..."); 
			} 
				
			current_file = (FileListClass *) current_file->NextInList (); 
			
			current_line++; 
	    } 
 
        log_file->Message ("Generating requested white fields only..."); 
 
		//buffer white fields 
		memset (white_field_buffer_sinograms, 0, sizeof(float) * recon_info_record.white_size); 
		current_file = top_white_file_list; 
		current_line = 0; 
		while (current_file != NULL) { 
			sprintf (file_name, "%s%s", temp_path, current_file->file_name); 
			if (stat (file_name, &stat_buffer) == 0) { 
				input_file.open (file_name, ios::binary); 
 
				data_size = sizeof (unsigned short) * recon_info_record.sinogram_xdim; 
    	        for (loop=0;loop<files_per_pass;loop++) { 
           	    	data_offset = header_size + sizeof (unsigned short) * sinogram_list[loop].sinogram_number * recon_info_record.sinogram_xdim; 
					input_file.seekg (data_offset, ios::beg); 
					input_file.read ((char *) &white_field_buffer_sinograms[(loop*recon_info_record.white_size)+(current_line*recon_info_record.sinogram_xdim)], data_size); 
               	} 
					 
				input_file.close (); 
			} 
       		else { 
	   	        for (loop=0;loop<sinograms_to_create;loop++) 
					memcpy (&white_field_buffer_sinograms[(loop*recon_info_record.white_size)+(current_line*recon_info_record.sinogram_xdim)], 
					        &white_field_buffer_sinograms[(loop*recon_info_record.white_size)+((current_line-1)*recon_info_record.sinogram_xdim)], 
							sizeof (unsigned short int) * recon_info_record.sinogram_xdim
					); 
 
				sprintf (msg, "Could not read file %s!", file_name); 
				error_log->addError (msg, "StartTomoMPIServer ()"); 
				error_log->addAutoResolution ("Using previous white field instead..."); 
       		} 
 
			current_file = (FileListClass *) current_file->NextInList ();

			current_line++; 
		} 
 
        log_file->Message ("Generating requested dark fields only..."); 
 
		//now average the dark fields--we'll send only the average to the clients. 
		memset (dark_field_buffer_sinograms, 0, sizeof(float) * recon_info_record.dark_size); 
		current_file = top_dark_file_list; 
		while (current_file != NULL) { 
			sprintf (file_name, "%s%s", temp_path, current_file->file_name); 
			if (stat (file_name, &stat_buffer) == 0) { 
				input_file.open (file_name, ios::binary); 
 
				data_size = sizeof (unsigned short) * recon_info_record.sinogram_xdim; 
   	            for (loop=0;loop<files_per_pass;loop++) { 
           	    	data_offset = header_size + sizeof (unsigned short) * sinogram_list[loop].sinogram_number * recon_info_record.sinogram_xdim; 
					input_file.seekg (data_offset, ios::beg); 
					input_file.read ((char *) temp_buffer, data_size); 
					for (int loop2=0;loop2<recon_info_record.sinogram_xdim;loop2++) 
						dark_field_buffer_sinograms[loop*recon_info_record.sinogram_xdim] += (float) temp_buffer[loop2] / (float) num_dark_fields; 
               	} 
 
				input_file.close (); 
			} 
        	else { 
				sprintf (msg, "Could not read file %s!", file_name); 
	            error_log->addError (msg, "StartTomoMPIServer ()"); 
				error_log->addAutoResolution ("This dark field will not be averaged in..."); 
        	} 
 
			current_file = (FileListClass *) current_file->NextInList (); 
		} 
		 
		log_file->Message ("Done generating requested slices."); 
		 
	} 
	else { 
        log_file->Message ("Generating first pass of sinograms..."); 
 
		header_size = 1024; 
		data_size = sizeof (unsigned short)*sinograms_to_create*recon_info_record.sinogram_xdim; 
		data_offset = header_size + sizeof (unsigned short)*pass_number*sinograms_to_create*recon_info_record.sinogram_xdim; 
 
		//buffer projections 
		current_file = top_projection_file_list; 
		current_line = 0; 

		//this let's us stop 1 file before the end--which is a redundant file anyway
		while (current_file->NextInList () != NULL) { 
			sprintf (file_name, "%s%s", temp_path, current_file->file_name); 
			if (stat (file_name, &stat_buffer) == 0) { 
				input_file.open (file_name, ios::binary); 
				input_file.seekg (data_offset, ios::beg); 
 
				input_file.read ((char *) temp_buffer, data_size); 
				input_file.close (); 
 
				for (loop=0;loop<sinograms_to_create;loop++) 
					memcpy (&buffer_sinograms[(loop*sinogram_size)+(current_line*recon_info_record.sinogram_xdim)], 
					        &temp_buffer[loop*recon_info_record.sinogram_xdim], 
							sizeof(unsigned short int)*recon_info_record.sinogram_xdim
					); 
			} 
        	else { 
   	        	for (loop=0;loop<sinograms_to_create;loop++) 
					memcpy (&buffer_sinograms[(loop*sinogram_size)+(current_line*recon_info_record.sinogram_xdim)], 
					        &buffer_sinograms[(loop*sinogram_size)+((current_line-1)*recon_info_record.sinogram_xdim)], 
							sizeof (unsigned short int) * recon_info_record.sinogram_xdim
					); 
 
				sprintf (msg, "Could not read file %s!", file_name); 
				error_log->addError (msg, "StartTomoMPIServer ()"); 
				error_log->addAutoResolution ("Using previous projection instead..."); 
        	} 
 
			current_file = (FileListClass *) current_file->NextInList (); 
			current_line++; 
		} 
 
		//buffer white fields 
		memset (white_field_buffer_sinograms, 0, sizeof(float) * recon_info_record.white_size); 
		current_file = top_white_file_list; 
		current_line = 0; 
		while (current_file != NULL) { 
			sprintf (file_name, "%s%s", temp_path, current_file->file_name); 
			if (stat (file_name, &stat_buffer) == 0) { 
				input_file.open (file_name, ios::binary); 
				input_file.seekg (data_offset, ios::beg); 
 
				input_file.read ((char *) temp_buffer, data_size); 
				input_file.close (); 
 
				for (loop=0;loop<sinograms_to_create;loop++) 
					memcpy (&white_field_buffer_sinograms[(loop*recon_info_record.white_size)+(current_line*recon_info_record.sinogram_xdim)], 
					        &temp_buffer[loop*recon_info_record.sinogram_xdim], 
							sizeof(unsigned short int)*recon_info_record.sinogram_xdim
					); 
			} 
        	else { 
	   	        for (loop=0;loop<sinograms_to_create;loop++) 
					memcpy (&white_field_buffer_sinograms[(loop*recon_info_record.white_size)+(current_line*recon_info_record.sinogram_xdim)], 
					        &white_field_buffer_sinograms[(loop*recon_info_record.white_size)+((current_line-1)*recon_info_record.sinogram_xdim)], 
							sizeof (unsigned short int) * recon_info_record.sinogram_xdim
					); 
 
				sprintf (msg, "Could not read file %s!", file_name); 
				error_log->addError (msg, "StartTomoMPIServer ()"); 
				error_log->addAutoResolution ("Using previous white field instead..."); 
        	} 
 
			current_file = (FileListClass *) current_file->NextInList (); 
			current_line++; 
		} 
 
		//now average the dark fields--we'll send only the average to the clients. 
		memset (dark_field_buffer_sinograms, 0, sizeof(float) * recon_info_record.dark_size); 
		current_file = top_dark_file_list; 
		while (current_file != NULL) { 
			sprintf (file_name, "%s%s", temp_path, current_file->file_name); 
			if (stat (file_name, &stat_buffer) == 0) { 
				input_file.open (file_name, ios::binary); 
				input_file.seekg (data_offset, ios::beg); 
 
				input_file.read ((char *) temp_buffer, data_size); 
				input_file.close (); 
 
				for (int loop=0;loop<recon_info_record.dark_size;loop++) 
					dark_field_buffer_sinograms[loop] += (float) temp_buffer[loop] / (float) num_dark_fields; 
			} 
        	else { 
				sprintf (msg, "Could not read file %s!", file_name); 
    	        error_log->addError (msg, "StartTomoMPIServer ()"); 
				error_log->addAutoResolution ("This dark field will not be averaged in..."); 
        	} 
 
			current_file = (FileListClass *) current_file->NextInList (); 
		} 
	} 
	 
	log_file->StopTimer (sinogram_pass_timer); 
	log_file->AccumulateTimer (sinogram_pass_timer); 
 
	log_file->TimeStamp ("Finished creating sinograms"); 
} 
 
//_____________________________________________________________________________________ 
 
void CreateSinogramsFromHDF (void) {

	FileListClass       *current_file = NULL; 

	int                 start_dims[2], 
						length_dims[2], 
						current_line, 
						loop; 

	char                temp_path[256], 
						full_file_name[512]; 

	struct stat         stat_buffer; 
 
 
    sprintf (msg, "\nLooks like files are HDF.  Starting pass %d", pass_number); 
    log_file->Message(msg); 
    sprintf (msg, "Starting to create %d sinograms", files_per_pass); 
    log_file->TimeStamp (msg); 
 
    log_file->StartTimer (sinogram_pass_timer); 
 
    strcpy (temp_path, exp_file_path); 
    strcat (temp_path, "raw/"); 
 
    //all files should have same template--so create it here 
    data_file.CreateTemplate (temp_path, top_projection_file_list->file_name); 
 
    if (recon_info_record.use_slices_file) { 
        log_file->Message ("Generating requested sinograms only..."); 
 
        //buffer projections 
        current_file = top_projection_file_list; 
        current_line = 0; 
        while (current_file->NextInList () != NULL) { 
            sprintf (full_file_name, "%s%s", temp_path, current_file->file_name); 
            if (stat (full_file_name, &stat_buffer) == 0) { 
                data_file.ChangeFile (temp_path, current_file->file_name); 
 
                if (!data_file.IndexExists(data_group_index)) { 
                    sprintf (msg, "ERROR: Data Group Index %s not found in file %s%s!", data_group_index, temp_path, current_file->file_name); 
                    log_file->Message (msg); 
                    error_log->addError (msg, "CreateSinogramsFromHDF ()"); 
                } 
 
                for (loop=0;loop<files_per_pass;loop++) { 
                    start_dims[0] = sinogram_list[loop].sinogram_number;    //this is the y dim in a 2D array 
                    start_dims[1] = 0;                                      //this is the x dim in a 2D array 
                    length_dims[0] = 1; 
                    length_dims[1] = recon_info_record.sinogram_xdim; 
 
                    data_file.GetDatumSlab (data_group_index, &buffer_sinograms[(loop*sinogram_size)+(current_line*recon_info_record.sinogram_xdim)], start_dims, length_dims); 
                } 
            } 
            else { 
                for (loop=0;loop<files_per_pass;loop++) 
                    memcpy (&buffer_sinograms[(loop*sinogram_size)+(current_line*recon_info_record.sinogram_xdim)], 
					        &buffer_sinograms[(loop*sinogram_size)+((current_line-1)*recon_info_record.sinogram_xdim)], 
							sizeof (unsigned short) * recon_info_record.sinogram_xdim
					); 
 
                sprintf (msg, "Could not read file %s!", full_file_name); 
                error_log->addError (msg, "StartTomoMPIServer ()"); 
                error_log->addAutoResolution ("Using previous projection instead..."); 
            } 
 
            current_file = (FileListClass *) current_file->NextInList (); 
            current_line++; 
        } 
 
        //buffer white fields 
		memset (white_field_buffer_sinograms, 0, sizeof(float) * recon_info_record.white_size); 
        current_file = top_white_file_list; 
        current_line = 0; 
        while (current_file != NULL) { 
            sprintf (full_file_name, "%s%s", temp_path, current_file->file_name); 
            if (stat (full_file_name, &stat_buffer) == 0) { 
                data_file.ChangeFile (temp_path, current_file->file_name); 
 
                if (!data_file.IndexExists(data_group_index)) { 
                    sprintf (msg, "ERROR: Data Group Index %s not found in file %s%s!", data_group_index, temp_path, current_file->file_name); 
                    log_file->Message (msg); 
                    error_log->addError (msg, "CreateSinogramsFromHDF ()"); 
                } 
 
                for (int loop=0;loop<files_per_pass;loop++) { 
                    start_dims[0] = sinogram_list[loop].sinogram_number;    //this is the y dim in a 2D array 
                    start_dims[1] = 0;                                      //this is the x dim in a 2D array 
                    length_dims[0] = 1; 
                    length_dims[1] = recon_info_record.sinogram_xdim; 
 
                    data_file.GetDatumSlab (data_group_index, &white_field_buffer_sinograms[(loop*recon_info_record.white_size)+(current_line*recon_info_record.sinogram_xdim)], start_dims, length_dims); 
                } 
            } 
            else { 
                for (loop=0;loop<sinograms_to_create;loop++) 
                    memcpy (&white_field_buffer_sinograms[(loop*recon_info_record.white_size)+(current_line*recon_info_record.sinogram_xdim)], 
					        &white_field_buffer_sinograms[(loop*recon_info_record.white_size)+((current_line-1)*recon_info_record.sinogram_xdim)], 
							sizeof (unsigned short int) * recon_info_record.sinogram_xdim
					); 
 
                sprintf (msg, "Could not read file %s!", full_file_name); 
                error_log->addError (msg, "StartTomoMPIServer ()"); 
                error_log->addAutoResolution ("Using previous white field instead..."); 
            } 
 
            current_file = (FileListClass *) current_file->NextInList (); 
            current_line++; 
        } 
 
        //now average the dark fields--we'll send only the average to the clients. 
		memset (dark_field_buffer_sinograms, 0, sizeof(float) * recon_info_record.dark_size); 
        current_file = top_dark_file_list; 
        current_line = 0; 
        while (current_file != NULL) { 
            sprintf (full_file_name, "%s%s", temp_path, current_file->file_name); 
            if (stat (full_file_name, &stat_buffer) == 0) { 
                data_file.ChangeFile (temp_path, current_file->file_name); 
 
                if (!data_file.IndexExists(data_group_index)) { 
                    sprintf (msg, "ERROR: Data Group Index %s not found in file %s%s!", data_group_index, temp_path, current_file->file_name); 
                    log_file->Message (msg); 
                    error_log->addError (msg, "CreateSinogramsFromHDF ()"); 
                } 
 
                for (int loop=0;loop<files_per_pass;loop++) { 
                    memset (&dark_field_buffer_sinograms[loop*recon_info_record.sinogram_xdim], 0, sizeof(float) * recon_info_record.sinogram_xdim); 
 
                    start_dims[0] = sinogram_list[loop].sinogram_number;    //this is the y dim in a 2D array 
                    start_dims[1] = 0;                                      //this is the x dim in a 2D array 
                    length_dims[0] = 1; 
                    length_dims[1] = recon_info_record.sinogram_xdim; 
 
                    data_file.GetDatumSlab (data_group_index, temp_buffer, start_dims, length_dims); 
 
                    for (int loop2=0;loop2<recon_info_record.sinogram_xdim;loop2++) 
                       dark_field_buffer_sinograms[(loop*recon_info_record.sinogram_xdim)+loop2] += (float) temp_buffer[loop2] / (float) num_dark_fields; 
                } 
            } 
            else { 
                sprintf (msg, "Could not read file %s!", full_file_name); 
                error_log->addError (msg, "StartTomoMPIServer ()"); 
                error_log->addAutoResolution ("This dark field will not be averaged in..."); 
            } 
 
            current_file = (FileListClass *) current_file->NextInList (); 
            current_line++; 
        } 
 
    } 
    else { 
        log_file->Message ("Generating first pass of sinograms..."); 
 
        start_dims[0] = pass_number*sinograms_to_create;    //this is the y dim in a 2D array 
        start_dims[1] = 0;                                  //this is the x dim in a 2D array 
        length_dims[0] = sinograms_to_create; 
        length_dims[1] = recon_info_record.sinogram_xdim; 
 
        //buffer projections 
        current_file = top_projection_file_list; 
        current_line = 0; 
		//this let's us stop 1 file before the end--which is a redundant file anyway
        while (current_file->NextInList () != NULL) { 
            sprintf (full_file_name, "%s%s", temp_path, current_file->file_name); 
            if (stat (full_file_name, &stat_buffer) == 0) { 
                data_file.ChangeFile (temp_path, current_file->file_name); 
 
                if (!data_file.IndexExists(data_group_index)) { 
                    sprintf (msg, "ERROR: Data Group Index %s not found in file %s%s!", data_group_index, temp_path, current_file->file_name); 
                    log_file->Message (msg); 
                    error_log->addError (msg, "CreateSinogramsFromHDF ()"); 
                } 
 
                data_file.GetDatumSlab (data_group_index, temp_buffer, start_dims, length_dims); 
 
                for (loop=0;loop<sinograms_to_create;loop++) 
                    memcpy (&buffer_sinograms[(loop*sinogram_size)+(current_line*recon_info_record.sinogram_xdim)], 
					        &temp_buffer[loop*recon_info_record.sinogram_xdim], 
							sizeof(unsigned short int)*recon_info_record.sinogram_xdim
					); 
            } 
            else { 
                for (loop=0;loop<sinograms_to_create;loop++) 
                    memcpy (&buffer_sinograms[(loop*sinogram_size)+(current_line*recon_info_record.sinogram_xdim)], 
					        &buffer_sinograms[(loop*sinogram_size)+((current_line-1)*recon_info_record.sinogram_xdim)], 
							sizeof (unsigned short int) * recon_info_record.sinogram_xdim
					); 
 
                sprintf (msg, "Could not read file %s!", full_file_name); 
                error_log->addError (msg, "StartTomoMPIServer ()"); 
                error_log->addAutoResolution ("Using previous projection instead..."); 
            } 
 
            current_file = (FileListClass *) current_file->NextInList (); 
            current_line++; 
        } 
 
        //buffer white fields 
		memset (white_field_buffer_sinograms, 0, sizeof(float) * recon_info_record.white_size); 
        current_file = top_white_file_list; 
        current_line = 0; 
        while (current_file != NULL) { 
            sprintf (full_file_name, "%s%s", temp_path, current_file->file_name); 
            if (stat (full_file_name, &stat_buffer) == 0) { 
                data_file.ChangeFile (temp_path, current_file->file_name); 
 
                if (!data_file.IndexExists(data_group_index)) { 
                    sprintf (msg, "ERROR: Data Group Index %s not found in file %s%s!", data_group_index, temp_path, current_file->file_name); 
                    log_file->Message (msg); 
                    error_log->addError (msg, "CreateSinogramsFromHDF ()"); 
                } 
 
                data_file.GetDatumSlab (data_group_index, temp_buffer, start_dims, length_dims); 
 
                for (loop=0;loop<sinograms_to_create;loop++) 
                    memcpy (&white_field_buffer_sinograms[(loop*recon_info_record.white_size)+(current_line*recon_info_record.sinogram_xdim)], 
					        &temp_buffer[loop*recon_info_record.sinogram_xdim], 
							sizeof(unsigned short int)*recon_info_record.sinogram_xdim
					); 
            } 
            else { 
                for (loop=0;loop<sinograms_to_create;loop++) 
                    memcpy (&white_field_buffer_sinograms[(loop*recon_info_record.white_size)+(current_line*recon_info_record.sinogram_xdim)], 
					        &white_field_buffer_sinograms[(loop*recon_info_record.white_size)+((current_line-1)*recon_info_record.sinogram_xdim)], 
							sizeof (unsigned short int) * recon_info_record.sinogram_xdim
					); 
 
                sprintf (msg, "Could not read file %s!", full_file_name); 
                error_log->addError (msg, "StartTomoMPIServer ()"); 
                error_log->addAutoResolution ("Using previous white field instead..."); 
            } 
 
            current_file = (FileListClass *) current_file->NextInList (); 
            current_line++; 
        } 
 
 
        //now average the dark fields--we'll send only the average to the clients. 
        memset (dark_field_buffer_sinograms, 0, sizeof(float) * recon_info_record.dark_size); 
        current_file = top_dark_file_list; 
        current_line = 0; 
        while (current_file != NULL) { 
            sprintf (full_file_name, "%s%s", temp_path, current_file->file_name); 
            if (stat (full_file_name, &stat_buffer) == 0) { 
                data_file.ChangeFile (temp_path, current_file->file_name); 
 
                if (!data_file.IndexExists(data_group_index)) { 
                    sprintf (msg, "ERROR: Data Group Index %s not found in file %s%s!", data_group_index, temp_path, current_file->file_name); 
                    log_file->Message (msg); 
                    error_log->addError (msg, "CreateSinogramsFromHDF ()"); 
                } 
 
                data_file.GetDatumSlab (data_group_index, temp_buffer, start_dims, length_dims); 
 
                for (int loop=0;loop<recon_info_record.dark_size;loop++) 
                    dark_field_buffer_sinograms[loop] += (float) temp_buffer[loop] / (float) num_dark_fields; 
            } 
            else { 
                sprintf (msg, "Could not read file %s!", full_file_name); 
                error_log->addError (msg, "StartTomoMPIServer ()"); 
                error_log->addAutoResolution ("This dark field will not be averaged in..."); 
            } 
 
            current_file = (FileListClass *) current_file->NextInList (); 
            current_line++; 
        } 
    } 
 
    log_file->StopTimer (sinogram_pass_timer); 
    log_file->AccumulateTimer (sinogram_pass_timer); 
 
    log_file->TimeStamp ("Finished creating sinograms"); 
 
} 
 
//_____________________________________________________________________________________ 
 
void *CreateSinograms (void *) { 
	char    file_name[256]; 
 
    sprintf (file_name, "%s", top_projection_file_list->file_name); 
 
    if (!strcmp (&file_name[strlen (file_name)-3], "bin")) 
        CreateSinogramsFromBinary (); 
    else 
        CreateSinogramsFromHDF (); 
 
	return (NULL); 
} 
 
//_____________________________________________________________________________________ 
 
void InitializeServer (void) { 

	FileListClass	*current_file = NULL; 
	int				loop, 
					num_files, 
					slice_index; 
	char			temp_path[256], 
					ch, 
					slice[10], 
					slices_file[512]; 
	FILE			*slices_list_file = NULL; 
	bool			end_file; 
 
    if (sinogram_list != NULL) { 
        free (sinogram_list); 
        sinogram_list = NULL; 
    } 

	sinogram_list = (SINOGRAMINFO *) malloc (sizeof (SINOGRAMINFO) * data_ydim); 
    if (sinogram_list == NULL) { 
        sprintf (msg, "Could not allocat memory for sinogram_list."); 
        error_log->addError (msg, "InitializeServer ()"); 
    } 
 
	if (recon_info_record.use_slices_file) { 
		log_file->Message ("I will be reconstructing only the files listed in slices.list."); 
 
		sprintf (slices_file, "%s%s", exp_file_path, "slices.list"); 
		if ((slices_list_file = fopen (slices_file, "rt")) == NULL) { 
			sprintf(msg, "%s%s", slices_file, " could not be opened."); 
			log_file->ErrorMessage (msg, "ReadConfigFile"); 
			log_file->Message ("...Ooops--I guess I will be reconstructing the entire data set after all..."); 
			error_log->addError (msg, "ReadConfigFile ()"); 
 
			recon_info_record.use_slices_file = false; 
		} 
 
		num_sinograms = 0; 
		ch = ' '; 
		end_file = false; 
		while (!end_file) { 
			while (ch != '<') 
				ch = fgetc (slices_list_file); 
 
			ch = fgetc (slices_list_file); 
 
			slice_index = 0; 
			while (ch != '>') { 
				slice[slice_index] = ch; 
				slice_index++; 
 
				ch = fgetc (slices_list_file); 
			} 
			slice[slice_index] = '\0'; 
 
			if (strcmp (slice, "END") == 0) 
				end_file = true; 
			else { 
				sinogram_list[num_sinograms].sinogram_number = atoi (slice); 
				sinogram_list[num_sinograms].status = NOT_DONE_YET; 
				sinogram_list[num_sinograms].processing_client = my_id; 
 
				num_sinograms++; 
			} 
		} 
 
		fclose (slices_list_file); 
 
		sprintf(msg, "Number of slices found in slices.list: %d", num_sinograms); 
		log_file->Message (msg); 
 
		//If an odd number of slices--boost to an even number by doing the last slice twice 
		if ((num_sinograms % 2) == 1) { 
			sprintf(msg, "Number of slices found is odd--I will do the last slice twice..."); 
			log_file->Message (msg); 
 
			sinogram_list[num_sinograms].sinogram_number = sinogram_list[num_sinograms-1].sinogram_number; 
			sinogram_list[num_sinograms].status = NOT_DONE_YET; 
			sinogram_list[num_sinograms].processing_client = my_id; 
 
            num_sinograms++; 
		} 
 
		if (files_per_pass > num_sinograms) 
			files_per_pass = num_sinograms; 
 
	} 
 
	if (!recon_info_record.use_slices_file) { 
		log_file->Message ("I will be reconstructing the entire data set."); 
		num_sinograms = data_ydim; 
 
		for (loop=0;loop<num_sinograms;loop++) { 
			sinogram_list[loop].sinogram_number = loop; 
			sinogram_list[loop].status = NOT_DONE_YET; 
			sinogram_list[loop].processing_client = my_id; 
		} 
 
		recon_info_record.start_fixed_shift = recon_info_record.fixed_shift_value; 
		recon_info_record.end_fixed_shift = recon_info_record.fixed_shift_value; 
	} 
 
    if (temp_buffer != NULL) { 
        free (temp_buffer); 
        temp_buffer = NULL; 
    } 
	temp_buffer = (unsigned short int *) malloc (sizeof(unsigned short int)*files_per_pass*data_xdim); 
    if (temp_buffer == NULL) { 
        sprintf (msg, "Could not allocat memory for temp_buffer."); 
        error_log->addError (msg, "InitializeServer ()"); 
    } 
    memset (temp_buffer, 0, sizeof(unsigned short int)*files_per_pass*data_xdim); 
 
	num_files = 0; 
	current_file = top_projection_file_list; 
	while (current_file != NULL) { 
		num_files++; 
		current_file = (FileListClass *) current_file->NextInList (); 
	} 
	num_files--;	//skip the last file--it's 180 degrees off the first file and so redundant 
 
	recon_info_record.sinogram_xdim = data_xdim; 
	recon_info_record.sinogram_ydim = num_files; 
	sinogram_size = recon_info_record.sinogram_xdim * recon_info_record.sinogram_ydim; 
 
	recon_info_record.reconstruction_xdim = data_xdim; 
	recon_info_record.reconstruction_ydim = data_xdim;  //Reconstruction should be square--for now. 
	reconstruction_size = recon_info_record.reconstruction_xdim * recon_info_record.reconstruction_ydim; 
 
    if (sino_buf1 != NULL) { 
        free (sino_buf1); 
        sino_buf1 = NULL; 
    } 
	sino_buf1 = (unsigned short int *) malloc (sizeof(unsigned short int)*files_per_pass*(sinogram_size+1)); 
    if (sino_buf1 == NULL) { 
        sprintf (msg, "Could not allocat memory for sino_buf1"); 
        error_log->addError (msg, "InitializeServer ()"); 
    } 
    memset (sino_buf1, 0, sizeof(unsigned short int)*files_per_pass*(sinogram_size+1)); 
 
    if (sino_buf2 != NULL) { 
        free (sino_buf2); 
        sino_buf2 = NULL; 
    } 
	sino_buf2 = (unsigned short int *) malloc (sizeof(unsigned short int)*files_per_pass*(sinogram_size+1)); 
    if (sino_buf2 == NULL) { 
        sprintf (msg, "Could not allocat memory for sino_buf2"); 
        error_log->addError (msg, "InitializeServer ()"); 
    } 
    memset (sino_buf2, 0, sizeof(unsigned short int)*files_per_pass*(sinogram_size+1)); 
 
	recon_info_record.num_white_fields = num_white_fields; 
	recon_info_record.white_size = num_white_fields * recon_info_record.sinogram_xdim; 
    if (white_field_sino_buf1 != NULL) { 
        free (white_field_sino_buf1); 
        white_field_sino_buf1 = NULL; 
    } 
	white_field_sino_buf1 = (unsigned short int *) malloc (files_per_pass*sizeof(unsigned short int)*recon_info_record.white_size); 
    if (white_field_sino_buf1 == NULL) { 
        sprintf (msg, "Could not allocat memory for white_field_sino_buf1"); 
        error_log->addError (msg, "InitializeServer ()"); 
    } 
    memset (white_field_sino_buf1, 0, files_per_pass*sizeof(unsigned short int)*recon_info_record.white_size); 
 
    if (white_field_sino_buf2 != NULL) { 
        free (white_field_sino_buf2); 
        white_field_sino_buf2 = NULL; 
    } 
	white_field_sino_buf2 = (unsigned short int *) malloc (files_per_pass*sizeof(unsigned short int)*recon_info_record.white_size); 
    if (white_field_sino_buf2 == NULL) { 
        sprintf (msg, "Could not allocat memory for white_field_sino_buf2"); 
        error_log->addError (msg, "InitializeServer ()"); 
    } 
    memset (white_field_sino_buf2, 0, files_per_pass*sizeof(unsigned short int)*recon_info_record.white_size); 
 
	recon_info_record.dark_size = files_per_pass * recon_info_record.sinogram_xdim; 
    if (dark_field_sino_ave_buf1 != NULL) { 
        free (dark_field_sino_ave_buf1); 
        dark_field_sino_ave_buf1 = NULL; 
    } 
	dark_field_sino_ave_buf1 = (float *) malloc (recon_info_record.dark_size*sizeof(float)); 
    if (dark_field_sino_ave_buf1 == NULL) { 
        sprintf (msg, "Could not allocat memory for dark_field_sino_ave_buf1"); 
        error_log->addError (msg, "InitializeServer ()"); 
    } 
    memset (dark_field_sino_ave_buf1, 0, recon_info_record.dark_size*sizeof(float)); 
 
    if (dark_field_sino_ave_buf2 != NULL) { 
        free (dark_field_sino_ave_buf2); 
        dark_field_sino_ave_buf2 = NULL; 
    } 
	dark_field_sino_ave_buf2 = (float *) malloc (recon_info_record.dark_size*sizeof(float)); 
    if (dark_field_sino_ave_buf2 == NULL) { 
        sprintf (msg, "Could not allocat memory for dark_field_sino_ave_buf2"); 
        error_log->addError (msg, "InitializeServer ()"); 
    } 
    memset (dark_field_sino_ave_buf2, 0, recon_info_record.dark_size*sizeof(float)); 
 
} 
 
//_____________________________________________________________________________________ 
 
void MakeFirstContact (void) { 
	int	loop, 
		temp; 
 
	log_file->Message ("Attempting first contact."); 
 
	for (loop=0;loop<num_processes;loop++) { 
		if (loop != my_id) { 
			//Send each client the following info 
			sprintf (msg, "Trying to contact %d.", loop); 
			log_file->Message (msg); 
 
			recon_info_record.sinogram_set_size  = sizeof (int); //sinogram_number 
			recon_info_record.sinogram_set_size += sizeof (short) * sinogram_size; //sinogram size 
			recon_info_record.sinogram_set_size += sizeof (short) * sinogram_size; //white field size 
			recon_info_record.sinogram_set_size += sizeof (float) * recon_info_record.sinogram_xdim; //dark field size 
            if (sinogram_data_set != NULL) { 
                free (sinogram_data_set); 
                sinogram_data_set = NULL; 
            } 
			sinogram_data_set = (char *) malloc (recon_info_record.sinogram_set_size); 
            if (sinogram_data_set == NULL) { 
                sprintf (msg, "Could not allocat memory for sinogram_data_set"); 
                error_log->addError (msg, "MakeFirstcontact ()"); 
            } 
    		memset (sinogram_data_set, 0, recon_info_record.sinogram_set_size); 
 
			MPI_Send ((void *) &recon_info_record, sizeof (recon_info_record), MPI_BYTE, loop, 0, MPI_COMM_WORLD); 
 
			//theta list must be seperate so clients can allocate memory... 
			MPI_Send (recon_info_record.theta_list, recon_info_record.theta_list_size, MPI_FLOAT, loop, 0, MPI_COMM_WORLD); 
		} 
	} 
 
	log_file->TimeStamp ("Everyone is awake"); 
 
} 
 
//_____________________________________________________________________________________ 

void MPIRequestingSinograms (int process_id, int sino_pass_number, int *current_sinogram, float current_shift) { 

	int	sinogram_number, 
		num_sinograms_requested, 
		send_command, 
		sinogram_to_send, 
		slot, 
		num_sinos_loop, 
		offset;

	unsigned short int	*sinogram_to_send_ptr    = NULL, 
						*white_field_to_send_ptr = NULL; 

	float	*dark_field_to_send_ptr = NULL; 
			MPI_Status	status; 
 
//	log_file->ResetTimer (waiting_on_MPI_timer); 
	log_file->StartTimer (waiting_on_MPI_timer); 
 
	MPI_Recv (&num_sinograms_requested, 1, MPI_INT, process_id, 0, MPI_COMM_WORLD, &status); 
 
	for (num_sinos_loop=0;num_sinos_loop<num_sinograms_requested;num_sinos_loop++) { 
	  sinogram_number = sinogram_list[sino_pass_number*files_per_pass+(*current_sinogram)].sinogram_number; 
 
	  sinogram_to_send = sinogram_number; 
	  sinogram_to_send_ptr = &server_sinograms[(*current_sinogram)*sinogram_size]; 
	  white_field_to_send_ptr = &white_field_server_sinograms[(*current_sinogram)*recon_info_record.white_size]; 
	  dark_field_to_send_ptr = &dark_field_server_sinograms[(*current_sinogram)*recon_info_record.sinogram_xdim]; 
 
	  *current_sinogram = *current_sinogram + 1; 
 
	  send_command = SERVER__Sending_Sinogram; 
	  MPI_Send ((void *) &send_command, 1, MPI_INT, process_id, 0, MPI_COMM_WORLD); 
 
	  offset = 0; 
	  memcpy (&sinogram_data_set[offset], &sinogram_to_send, sizeof (int)); 
 
	  offset += sizeof (int); 
	  memcpy (&sinogram_data_set[offset], sinogram_to_send_ptr, sizeof (short) * sinogram_size); 
 
	  offset += sizeof (short) * sinogram_size; 
	  memcpy (&sinogram_data_set[offset], white_field_to_send_ptr, sizeof (short) * recon_info_record.white_size); 
 
	  offset += sizeof (short) * recon_info_record.white_size; 
	  memcpy (&sinogram_data_set[offset], dark_field_to_send_ptr, sizeof (float) * recon_info_record.sinogram_xdim); 
 
	  sprintf (msg, "Handing out sinogram %d with shift %f", sinogram_to_send, current_shift); 

	  log_file->TimeStamp (msg); 

	  MPI_Send (&current_shift, 1, MPI_FLOAT, process_id, 0, MPI_COMM_WORLD); 
	  MPI_Send (sinogram_data_set, recon_info_record.sinogram_set_size, MPI_BYTE, process_id, 0, MPI_COMM_WORLD); 
 
	  num_files_handled++; 
 
	  sinogram_list[sinogram_to_send].status = BEING_PROCESSED; 
	  sinogram_list[sinogram_to_send].processing_client = process_id; 
 
	  if (*current_sinogram == files_per_pass) 
	    break; 
	} 
 
	send_command = SERVER__Start_Calculations; 
	MPI_Send (&send_command, 1, MPI_INT, process_id, 0, MPI_COMM_WORLD); 
 
	log_file->StopTimer (waiting_on_MPI_timer); 
	log_file->AccumulateTimer (waiting_on_MPI_timer); 
//	log_file->TimerMessage (waiting_on_MPI_timer); 
} 
 
//_____________________________________________________________________________________ 
 
int MPIWaitForCommand (int *process_id) { 
	int			recieve_command; 
	MPI_Status	status; 
 
//	log_file->ResetTimer (waiting_on_clients_timer); 
	log_file->StartTimer (waiting_on_clients_timer); 
 
	MPI_Recv (process_id, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status); 
	MPI_Recv (&recieve_command, 1, MPI_INT, *process_id, 0, MPI_COMM_WORLD, &status); 
 
	log_file->StopTimer (waiting_on_clients_timer); 
	log_file->AccumulateTimer (waiting_on_clients_timer); 
//	log_file->TimerMessage (waiting_on_clients_timer); 
 
	return (recieve_command); 
} 
 
//_____________________________________________________________________________________ 
 
void Reconstruct (void) { 

	int	current_sinogram, 
    	process_id, 
    	send_command, 
    	process_count, 
    	mpi_command, 
    	requested_slice, 
    	slot, 
    	num_passes, 
    	temp;
 
	float current_shift;

	pthread_t	sinogram_thread_handle;

	MPI_Status	status; 

	sinogram_thread_running = false; 

	pass_number = 0; 

	buffer_sinograms = sino_buf1; 
	server_sinograms = sino_buf2; 
	white_field_buffer_sinograms = white_field_sino_buf1; 
	dark_field_buffer_sinograms = dark_field_sino_ave_buf1; 

	sinograms_to_create = files_per_pass; 
	CreateSinograms (NULL); 

	server_sinograms = sino_buf1; 
	buffer_sinograms = sino_buf2; 
	white_field_server_sinograms = white_field_sino_buf1; 
	dark_field_server_sinograms = dark_field_sino_ave_buf1; 
	sinogram_buffer = 1; 
 
	if ((num_sinograms % files_per_pass) == 0) {
		num_passes = (num_sinograms/files_per_pass); 
		sprintf (msg, "num_sinograms is a multiple of files_per_pass.  Num_passes: %d", num_passes); 
		log_file->Message (msg); 
	} 
	else { 
		num_passes = (num_sinograms/files_per_pass)+1;	//1 extra pass to clean up extra files 
		sprintf (msg, "num_sinograms is not a multiple of files_per_pass.  Num_passes: %d", num_passes); 
		log_file->Message (msg); 
	} 
 
	for (pass_number=1;pass_number<=num_passes;pass_number++) {
		if (sinogram_thread_running) {
			pthread_join (sinogram_thread_handle, NULL);
			sinogram_thread_running = false;
		}
	
		if (sinogram_buffer == 1) {
			server_sinograms = sino_buf1; 
			white_field_server_sinograms = white_field_sino_buf1; 
			dark_field_server_sinograms = dark_field_sino_ave_buf1; 
	
			buffer_sinograms = sino_buf2; 
			white_field_buffer_sinograms = white_field_sino_buf2; 
			dark_field_buffer_sinograms = dark_field_sino_ave_buf2; 
	
			sinogram_buffer = 2; 
		} 
		else { 
			server_sinograms = sino_buf2; 
			white_field_server_sinograms = white_field_sino_buf2; 
			dark_field_server_sinograms = dark_field_sino_ave_buf2; 
	
			buffer_sinograms = sino_buf1; 
			white_field_buffer_sinograms = white_field_sino_buf1; 
			dark_field_buffer_sinograms = dark_field_sino_ave_buf1; 
	
			sinogram_buffer = 1; 
		} 

		//these should be full passes guaranteed
		if (pass_number < num_passes) {  
		
			sinograms_to_create = files_per_pass; 
	
			sinogram_thread_running = true; 
			pthread_create (&sinogram_thread_handle, NULL, CreateSinograms, NULL); 
		} 
		else if (pass_number = num_passes) { //maybe a full pass/maybe not 
		
			sinograms_to_create = num_sinograms - ((pass_number-1) * files_per_pass); 
	
			sprintf (msg, "Handling extra files.  files_per_pass: %d", files_per_pass); 
			log_file->Message (msg); 
		} 
	
		current_shift = recon_info_record.start_fixed_shift; 

		// Note: be sure to set recon_info_record.fixed_shift_interval positive for full reconstruction
		while (current_shift <= recon_info_record.end_fixed_shift + recon_info_record.fixed_shift_interval / 2.0) { 
		
			current_sinogram = 0; 
			while (current_sinogram < sinograms_to_create) {
				mpi_command = MPIWaitForCommand (&process_id); 
				switch (mpi_command) {

					case CLIENT__Request_New_Sinogram : { 
						MPIRequestingSinograms (process_id, (pass_number-1), &current_sinogram, current_shift); 
						sprintf (msg, "MPIRequestingSinograms called! current_shift is %f", current_shift);  // test
						log_file->Message (msg);
					}; break; 
					
					default : {
						sprintf (msg, "Recieved invalid command %d", mpi_command); 
						log_file->TimeStamp (msg);
					}; break; 

				} 
			} 
	
			current_shift += recon_info_record.fixed_shift_interval;   
		} 
	} 
 
	log_file->TimeStamp ("Sent last sinogram."); 
} 
 
//_____________________________________________________________________________________ 
 
void DestroyServer (void) {

	int	process_id, 
		send_command, 
		process_count, 
		mpi_command, 
		requested_slice, 
		slot, 
		temp;

	float	temp_data_range_min, 
			temp_data_range_max;

	MPI_Status	status; 

	data_range_min =  1000000.0; 
	data_range_max = -1000000.0; 
 
	//Wait for everyone to be told to stop. 
	process_count = num_processes - 1; //subtract 1 for the sinogram server 

	while (process_count != 0) { 
		mpi_command = MPIWaitForCommand (&process_id); 

		switch (mpi_command) { 
			case CLIENT__Exiting : { 
				MPI_Recv (&temp_data_range_min, 1, MPI_FLOAT, process_id, 0, MPI_COMM_WORLD, &status); 
				MPI_Recv (&temp_data_range_max, 1, MPI_FLOAT, process_id, 0, MPI_COMM_WORLD, &status); 

				if (temp_data_range_min < data_range_min) 
					data_range_min = temp_data_range_min; 
		
				if (temp_data_range_max > data_range_max) 
					data_range_max = temp_data_range_max; 

				process_count--; 

				sprintf (msg, "Process %d is exiting--process count = %d", process_id, process_count); 
				log_file->TimeStamp (msg); 
			}; break; 

			case CLIENT__Request_New_Sinogram : { 
				MPI_Recv (&temp, 1, MPI_INT, process_id, 0, MPI_COMM_WORLD, &status); 

				send_command = STOP; 
				MPI_Send (&send_command, 1, MPI_INT, process_id, 0, MPI_COMM_WORLD); 

				sprintf (msg, "Requested process %d to exit", process_id); 
				log_file->TimeStamp (msg); 
			}; break; 

			default : { 
				sprintf (msg, "Recieved invalid command %d from %d", mpi_command, process_id); 
				log_file->TimeStamp (msg); 
			}; break; 
		} 

	} 

	sprintf (msg, "Min data range: %e", data_range_min); 
	log_file->TimeStamp (msg); 
	sprintf (msg, "Max data range: %e", data_range_max); 
	log_file->TimeStamp (msg); 

	UpdateExpFile (); 

	if (sino_buf1 != NULL) { 
		free (sino_buf1); 
		sino_buf1 = NULL; 
	} 

	if (sino_buf2 != NULL) { 
		free (sino_buf2); 
		sino_buf2 = NULL; 
	} 

	if (sinogram_list != NULL) { 
		free (sinogram_list); 
		sinogram_list = NULL; 
	} 

	if (temp_buffer != NULL) { 
		free (temp_buffer); 
		temp_buffer = NULL; 
	} 

	if (recon_info_record.theta_list != NULL) { 
		free (recon_info_record.theta_list); 
		recon_info_record.theta_list = NULL; 
	} 

	if (white_field_sino_buf2 != NULL) { 
		free (white_field_sino_buf1); 
		white_field_sino_buf1 = NULL; 
	} 
	if (white_field_sino_buf2 != NULL) { 
		free (white_field_sino_buf2); 
		white_field_sino_buf2 = NULL; 
	} 

	if (dark_field_sino_ave_buf1 != NULL) { 
		free (dark_field_sino_ave_buf1); 
		dark_field_sino_ave_buf1 = NULL; 
	} 
	if (dark_field_sino_ave_buf2 != NULL) { 
		free (dark_field_sino_ave_buf2); 
		dark_field_sino_ave_buf2 = NULL; 
	} 

	if (sinogram_data_set != NULL); { 
	free (sinogram_data_set); 
	sinogram_data_set = NULL; 
	} 

	if (top_projection_file_list != NULL) { 
		delete (top_projection_file_list); 
		top_projection_file_list = NULL; 
	} 
	if (top_white_file_list != NULL) { 
		delete (top_white_file_list); 
		top_white_file_list = NULL; 
	} 
	if (top_dark_file_list != NULL) { 
		delete (top_dark_file_list); 
		top_dark_file_list = NULL; 
	} 
 
} 
 
//_____________________________________________________________________________________ 
 
void ServerProcess (int argc, char *argv[]) {

	log_file->ResetTimer (total_reconstruction_timer); 
	log_file->StartTimer (total_reconstruction_timer); 

	sprintf (msg, "I am processor %s with id %d", processor_name, my_id); 
	log_file->Message (msg); 
	log_file->Message ("I am in control of your reconstruction..."); 

	InitSampleLocation (argc, argv); 

	sprintf (msg, "Today we will be trying to reconstruct %s%s.", exp_file_path, exp_file_name); 
	log_file->Message (msg); 

	num_files_handled = 0; 
	log_file->Message ("Starting to initialize..."); 
	ReadExpFile(); 
	ReadOverrideConfigFile(); 
	InitializeServer (); 

	WriteConfigFile (); 

	log_file->Message ("Initialization complete--starting to make first contact..."); 
	log_file->Message(" "); 

	MakeFirstContact (); 

	log_file->Message ("First contact made--starting to perform reconstruction..."); 
	Reconstruct (); 

	log_file->StopTimer (total_reconstruction_timer); 
	log_file->AccumulateTimer (total_reconstruction_timer); 

	DestroyServer (); 
 
} 
 
//_____________________________________________________________________________________ 
 
void CreateServerTimers (void) {

	total_reconstruction_timer = log_file->CreateTimer ("Full_Reconstruction"); 
	log_file->ResetTimer (total_reconstruction_timer); 
 
	sinogram_pass_timer = log_file->CreateTimer ("Sinogram_Pass"); 
	log_file->ResetTimer (sinogram_pass_timer); 
 
	waiting_on_clients_timer = log_file->CreateTimer ("Wait_On_Client"); 
	log_file->ResetTimer (waiting_on_clients_timer); 
 
	waiting_on_MPI_timer = log_file->CreateTimer ("Wait_On_MPI_Send_Sinogram"); 
	log_file->ResetTimer (waiting_on_MPI_timer); 
} 
 
//_____________________________________________________________________________________ 
 
void DestroyServerTimers (void) {

	log_file->DestroyTimer(total_reconstruction_timer); 
	log_file->DestroyTimer(sinogram_pass_timer); 
	log_file->DestroyTimer(waiting_on_clients_timer); 
	log_file->DestroyTimer(waiting_on_MPI_timer); 
} 
 
//_____________________________________________________________________________________ 
 
int StartTomoMPIServer (int argc, char* argv[], char *log_path) { 

	log_file->TimeStamp ("Server process starting"); 
	sprintf (msg, "<Number_of_Processes>%d", num_processes); 
	log_file->Message(msg); 

	log_file->Message (""); 
	log_file->Message ("Command line arguments:"); 
	for (int loop=0;loop<argc;loop++)  
	log_file->Message (argv[loop]); 

	log_file->Message (""); 

	ServerAcknowledgements (); 

	temp_buffer              = NULL; 
	sino_buf1                = NULL; 
	sino_buf2                = NULL; 
	sinogram_list            = NULL; 
	white_field_sino_buf1    = NULL; 
	white_field_sino_buf2    = NULL; 
	dark_field_sino_ave_buf1 = NULL; 
	dark_field_sino_ave_buf2 = NULL; 

	recon_info_record.mayor_id = my_id; 

	//set defaults 
	strcpy (recon_info_record.log_path, log_path); 
	recon_info_record.recon_algorithm = RECONSTRUCTION_GRIDREC; 
	recon_info_record.gridrec_padding = GRIDREC_PADDING_HALF;  
	recon_info_record.filter = FILTER_NONE; 
	recon_info_record.use_slices_file = false; 
	files_per_pass = DEFAULT_FILES_PER_PASS; 

	recon_info_record.start_fixed_shift = 0.0f; 
	recon_info_record.end_fixed_shift = 0.0f; 
	recon_info_record.fixed_shift_interval = 1.0f;   
	recon_info_record.fixed_shift_value = 0.0f; 

	recon_info_record.centering = 0; 
	recon_info_record.use_ring_removal = 1; 
	recon_info_record.ring_removal_coeff = 1; 
	recon_info_record.average_white_fields = false; 
	strcpy (data_group_index, ";entry1;data;data"); 
	recon_info_record.file_format = HDF5; 
	recon_info_record.compression_type = NX_COMP_NONE; 
	recon_info_record.debug = DEBUG_NONE; 

	recon_info_record.sinogram_xdim = 0; 
	recon_info_record.sinogram_ydim =0; 
	recon_info_record.white_size = 0; 
	recon_info_record.dark_size = 0; 
	recon_info_record.whitedark_interval = 0; 
	recon_info_record.reconstruction_xdim = 0; 
	recon_info_record.reconstruction_ydim = 0; 
	recon_info_record.theta_list_size = 0; 
	sprintf (recon_info_record.reconstruction_path, "empty"); 
	sprintf (recon_info_record.base_name, "empty"); 
	recon_info_record.theta_list = NULL; 

	CreateServerTimers (); 

	ServerProcess (argc, argv); 

	log_file->TimerMessage(total_reconstruction_timer); 
	log_file->TimerMessage(sinogram_pass_timer); 
	log_file->TimerMessage(waiting_on_clients_timer); 
	log_file->TimerMessage(waiting_on_MPI_timer); 

	DestroyServerTimers (); 

	sprintf (msg, "Process %d exiting.", my_id); 
	log_file->Message (msg); 

	log_file->TimeStamp ("Server process exiting"); 

	return 0; 
} 
//_____________________________________________________________________________________ 
 
