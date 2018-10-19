#ifndef MPICOMMUNICATIONSTRUCTS_H
#define MPICOMMUNICATIONSTRUCTS_H

//________________________________________________________________________________________

typedef struct RECON_INFO {
int	mayor_id,
	sinogram_set_size,
	debug,
	recon_algorithm,
        gridrec_padding, 
	sinogram_xdim,
	sinogram_ydim,
	num_white_fields,
	white_size,
	dark_size,
	whitedark_interval,
	reconstruction_xdim,
	reconstruction_ydim,
	theta_list_size,
	centering,
	use_slices_file, 
	use_fixed_shift, 
	use_ring_removal,
	file_format,
	compression_type,
	average_white_fields,
	filter,
	rescale_to_int;

float   fixed_shift_value,
	start_fixed_shift, 
        end_fixed_shift,
        fixed_shift_interval;   // NEW

char	reconstruction_path[256],
        log_path[256],
	base_name[256];

float	*theta_list, 
	ring_removal_coeff, 
	scale_data_range_min,
	scale_data_range_max;
};

//________________________________________________________________________________________

#endif
