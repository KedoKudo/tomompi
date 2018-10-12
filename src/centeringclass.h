#ifndef centeringclassH
#define centeringclassH
//---------------------------------------------------------------------------

#include "math.h"
#include "gridrec.h"

//---------------------------------------------------------------------------

#define MAX_SHIFT 25

//---------------------------------------------------------------------------

class CenteringClass
{
private:
ReconAlgorithm      *recon_algorithm;
int                 sinogram_x_dim, sinogram_y_dim;
float               *sinogram1,
                    *sinogram2,
                    *recon1,
                    *recon2,
                    *shifted_data,
                    *mean_vect,
                    *low_pass_sino_lines_data,
                    *mean_sino_line_data;
unsigned short      *wtemp;

	void LogProj (float *data);
	void Gridrec ();

public:

    CenteringClass (void):
		recon_algorithm(NULL),
		sinogram1(NULL),sinogram2(NULL), 
		recon1(NULL),recon2(NULL), 
		shifted_data(NULL), 
		wtemp(NULL), 
		mean_vect(NULL),
		low_pass_sino_lines_data(NULL), 
    mean_sino_line_data(NULL){}

	void init (ReconAlgorithm *recon_algorithm);

	void FindCenter (float *origional_sinogram1, 
	                 float *origional_sinogram2, 
					 float *shift_1, 
					 float *shift_2, 
					 float ring_coeff
					 );

	void FindCenter (float *origional_sinogram1, 
	                 float *origional_sinogram2, 
					 int *shift_1, 
					 int *shift_2, 
					 float ring_coeff
					 );

	void OffCenterCorrSingleManual (float *data, 
	                                int   shift
									);

	void RingCorrectionSingle (float *data, 
	                           float ring_coeff
							   );

	static void acknowledgements (LogFileClass *acknowledge_file);

	~CenteringClass ();
};

//---------------------------------------------------------------------------
#endif
