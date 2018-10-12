//---------------------------------------------------------------------------

#pragma hdrstop

#include "centeringclass.h"
//---------------------------------------------------------------------------

#pragma package(smart_init)

//---------------------------------------------------------------------------

void CenteringClass::acknowledgements (LogFileClass *acknowledge_file)
{
    acknowledge_file->Message ("__________________________________________________________________");
    acknowledge_file->Message ("Centering class");
    acknowledge_file->Message ("");
    acknowledge_file->Message ("Class for performing auto-centering of sinograms as well as ring removal.");
    acknowledge_file->Message ("Origional algorithm developed and maintened in C by:"); 
    acknowledge_file->Message ("       Antonio Brunetti"); 
    acknowledge_file->Message ("C++ class developed and maintained by:"); 
    acknowledge_file->Message ("       Brian Tieman & Francesco DeCarlo"); 
    acknowledge_file->Message ("       Argonne National Laboratory"); 
    acknowledge_file->Message ("       tieman@aps.anl.gov"); 
    acknowledge_file->Message (""); 
    acknowledge_file->Message ("8/20/2003  V1.0   BT  First version with acknowledgements"); 
    acknowledge_file->Message ("8/20/2003  V1.0   BT  Ported C code to a CPP object structure"); 
    acknowledge_file->Message ("9/10/2003  V1.2   BT  In the LogProj routine, added check that if data is < 0"); 
    acknowledge_file->Message ("        set it to 1.0.  This avoids taking the log of a negative number which"); 
    acknowledge_file->Message ("        was occasionally causing the reconstruction routine to blow up.  Should"); 
    acknowledge_file->Message ("        look into why data is negative in the first place..."); 
    acknowledge_file->Message ("9/11/2003  V1.2   BT  Class now centers two sinograms--taking advantage"); 
    acknowledge_file->Message ("        of gridrec's ability to reconstruct 2 sinograms.  If you only need one"); 
    acknowledge_file->Message ("        sinogram centered--pass the same sinogram twice."); 
    acknowledge_file->Message ("2/26/2004  V1.3   BT  The centering routine will now abort if the shift of"); 
    acknowledge_file->Message ("        either sinogram is > 25.  This is hard coded in the define MAX_SHIFT."); 
    acknowledge_file->Message ("3/15/2004  V1.5   BT  Fixed and error with the aborting of centering when MAX_SHIFT"); 
    acknowledge_file->Message ("        is exceeded.  Turns out I was only abborting on positive shifts.  The code"); 
    acknowledge_file->Message ("        now compares abs (shift) to MAX_SHIFT."); 
	acknowledge_file->Message ("1/31/2006  V1.10  BT  Added a ring coefficient parameter which can be 0= < coeff <= 5."); 
	acknowledge_file->Message ("		This was done in void CenteringClass::RingCorrectionSingle (float *data, float ring_coeff)"); 
    acknowledge_file->Message ("		Specifically the code"); 
    acknowledge_file->Message ("			if ((data[i*sinogram_x_dim+j]-tmp)>0.0)"); 
    acknowledge_file->Message ("			    data[i*sinogram_x_dim+j] -= tmp;"); 
    acknowledge_file->Message ("		    else"); 
	acknowledge_file->Message ("		        data[i*sinogram_x_dim+j] = 0.0;"); 
	acknowledge_file->Message ("		was changed to"); 
    acknowledge_file->Message ("			if ((data[i*sinogram_x_dim+j]-(tmp * ring_coeff))>0.0)"); 
    acknowledge_file->Message ("			    data[i*sinogram_x_dim+j] -= (tmp * ring_coeff);"); 
    acknowledge_file->Message ("		    else"); 
	acknowledge_file->Message ("		        data[i*sinogram_x_dim+j] = 0.0;"); 
    acknowledge_file->Message ("3/22/2006 V1.10  BT  A serious problem was found in the LogProj routine.  The LogProj routine"); 
    acknowledge_file->Message ("		attempted to normalize the data to the average pixel value of the first 10 pixels plus"); 
    acknowledge_file->Message ("		the last ten pixels over every row in the sinogram.  Due to background intensity variations"); 
    acknowledge_file->Message ("		from projection to projection, this mean value would cause data to be clipped in data"); 
    acknowledge_file->Message ("		sets with low contrast.  The resulting reconstruction looked like a reconstruction with"); 
    acknowledge_file->Message ("		slices missing--which it was due to the data being clipped.  The best fix found was to"); 
    acknowledge_file->Message ("		normalize each row to the maximum data value in that row.  This is how the new LogProj"); 
    acknowledge_file->Message ("		works."); 
    acknowledge_file->Message (""); 
    acknowledge_file->Message (""); 
    acknowledge_file->Message ("__________________________________________________________________"); 
}

//---------------------------------------------------------------------------

void CenteringClass::init (ReconAlgorithm *recon_algorithm)
{
    this->recon_algorithm = recon_algorithm;

    this->sinogram_x_dim = recon_algorithm->getSinogramXDimension (); 
    this->sinogram_y_dim = recon_algorithm->getSinogramYDimension (); 
 
    sinogram1= (float *) malloc(sizeof(float)*sinogram_x_dim*sinogram_y_dim); 
    sinogram2= (float *) malloc(sizeof(float)*sinogram_x_dim*sinogram_y_dim); 
 
    recon1 = (float *) malloc (sizeof(float)*sinogram_x_dim*sinogram_x_dim); 
    recon2 = (float *) malloc (sizeof(float)*sinogram_x_dim*sinogram_x_dim); 
 
    shifted_data = (float *) malloc (sizeof(float)*sinogram_x_dim*sinogram_y_dim+1); 
 
    wtemp = (unsigned short *) malloc(sizeof(unsigned short)*sinogram_x_dim); 
 
    mean_vect = (float  *) malloc (sizeof(float)*sinogram_y_dim); 
 
    low_pass_sino_lines_data = (float  *) malloc (sizeof(float)*sinogram_x_dim); 
    mean_sino_line_data = (float *) malloc (sizeof(float)*sinogram_x_dim); 
} 

//---------------------------------------------------------------------------

void CenteringClass::FindCenter (float *origional_sinogram1, float *origional_sinogram2, 
				 int *shift_1, int *shift_2, float ring_coeff)
{

// ********************************************************************* 
// important: sinogram cointains the logarithmic sinogram data 
int                 i, 
                    j; 
int                 shift_value1[3], 
                    shift_value2[3],
                    pixel_count1[3],
                    pixel_count2[3],
                    best_count1,
                    best_count2,
                    shift1,
                    shift2,
                    step1,
                    step2,
                    count;
bool                shift1_found,
                    shift2_found;
 
    shift1_found = false; 
    shift2_found = false; 
 
    for (i=-1;i<=1;i++) 
    { 
        memcpy (sinogram1, origional_sinogram1, sizeof(float)*sinogram_x_dim*sinogram_y_dim); 
        memcpy (sinogram2, origional_sinogram2, sizeof(float)*sinogram_x_dim*sinogram_y_dim); 
 
        OffCenterCorrSingleManual (sinogram1, i); 
        RingCorrectionSingle (sinogram1, ring_coeff); 
 
        OffCenterCorrSingleManual (sinogram2, i); 
        RingCorrectionSingle (sinogram2, ring_coeff); 
 
        //reconstruct both sinograms 
        Gridrec (); 
 
        count=0; 
        for (j=0;j<(sinogram_x_dim*sinogram_x_dim);j++) 
            if (recon1[j] != 0) 
                count += 1; 
 
        shift_value1[i+1] = i; 
        pixel_count1[i+1] = count; 
 
        count=0; 
        for (j=0;j<(sinogram_x_dim*sinogram_x_dim);j++) 
            if (recon2[j] != 0) 
                count+=1; 
 
        shift_value2[i+1] = i; 
        pixel_count2[i+1] = count; 
    } 
 
    //if 0 is the best shift, the shift has been found! 
    if ((pixel_count1[1] <= pixel_count1[0]) && (pixel_count1[1] <= pixel_count1[2])) 
    { 
        shift1 = shift_value1[1]; 
        shift1_found = true; 
    } 
 
    if ((pixel_count2[1] <= pixel_count2[0]) && (pixel_count2[1] <= pixel_count2[2])) 
    { 
        shift2 = shift_value2[1]; 
        shift2_found = true; 
    } 
 
    //if -1 is the best shift--continue to step negative--otherwise step positive 
    if (pixel_count1[0] <= pixel_count1[1]) 
    { 
        shift1 = shift_value1[0]; 
        best_count1 = pixel_count1[0]; 
        step1 = -1; 
    } 
    else 
    { 
        shift1 = shift_value1[2]; 
        best_count1 = pixel_count1[2]; 
        step1 = 1; 
    } 
 
    if (pixel_count2[0] <= pixel_count2[1]) 
    { 
        shift2 = shift_value2[0]; 
        best_count2 = pixel_count2[0]; 
        step2 = -1; 
    } 
    else 
    { 
        shift2 = shift_value2[2]; 
        best_count2 = pixel_count2[2]; 
        step2 = 1; 
    } 
 
    //continue to step until both are done 
    while ((!shift1_found || !shift2_found) && 
            (abs (shift1) < MAX_SHIFT) && 
            (abs (shift2) < MAX_SHIFT)) 
    { 
        if (!shift1_found) 
            shift1 = shift1 + step1; 
 
        if (!shift2_found) 
            shift2 = shift2 + step2; 
 
        memcpy (sinogram1, origional_sinogram1, sizeof(float)*sinogram_x_dim*sinogram_y_dim); 
        memcpy (sinogram2, origional_sinogram2, sizeof(float)*sinogram_x_dim*sinogram_y_dim); 
 
        OffCenterCorrSingleManual (sinogram1, shift1); 
        RingCorrectionSingle (sinogram1, ring_coeff); 
 
        OffCenterCorrSingleManual (sinogram2, shift2); 
        RingCorrectionSingle (sinogram2, ring_coeff); 
 
        //reconstruct both sinograms 
        Gridrec (); 
 
        count=0; 
        for (j=0;j<(sinogram_x_dim*sinogram_x_dim);j++) 
            if (recon1[j] != 0) 
                count+=1; 
 
        if (count >= best_count1) 
        { 
            shift1 = shift1 - step1; 
            shift1_found = true; 
        } 
        else 
            best_count1 = count; 
 
        count=0; 
        for (j=0;j<(sinogram_x_dim*sinogram_x_dim);j++) 
            if (recon2[j] != 0) 
                count+=1; 
 
        if (count >= best_count2) 
        { 
            shift2 = shift2 - step2; 
            shift2_found = true; 
        } 
        else 
            best_count2 = count; 
    } 
 
    *shift_1 = shift1; 
    *shift_2 = shift2; 
} 


void CenteringClass::FindCenter (float *origional_sinogram1, float *origional_sinogram2, 
				 float *shift_1, float *shift_2, float ring_coeff)
{

// ********************************************************************* 
// important: sinogram cointains the logarithmic sinogram data 
int                 i, 
                    j; 
int                 shift_value1[3], 
                    shift_value2[3],
                    pixel_count1[3],
                    pixel_count2[3],
                    best_count1,
                    best_count2,
                    shift1,
                    shift2,
                    step1,
                    step2,
                    count;
bool                shift1_found,
                    shift2_found;
 
    shift1_found = false; 
    shift2_found = false; 
 
    for (i=-1;i<=1;i++) 
    { 
        memcpy (sinogram1, origional_sinogram1, sizeof(float)*sinogram_x_dim*sinogram_y_dim); 
        memcpy (sinogram2, origional_sinogram2, sizeof(float)*sinogram_x_dim*sinogram_y_dim); 
 
        OffCenterCorrSingleManual (sinogram1, i); 
        RingCorrectionSingle (sinogram1, ring_coeff); 
 
        OffCenterCorrSingleManual (sinogram2, i); 
        RingCorrectionSingle (sinogram2, ring_coeff); 
 
        //reconstruct both sinograms 
        Gridrec (); 
 
        count=0; 
        for (j=0;j<(sinogram_x_dim*sinogram_x_dim);j++) 
            if (recon1[j] != 0) 
                count += 1; 
 
        shift_value1[i+1] = i; 
        pixel_count1[i+1] = count; 
 
        count=0; 
        for (j=0;j<(sinogram_x_dim*sinogram_x_dim);j++) 
            if (recon2[j] != 0) 
                count+=1; 
 
        shift_value2[i+1] = i; 
        pixel_count2[i+1] = count; 
    } 
 
    //if 0 is the best shift, the shift has been found! 
    if ((pixel_count1[1] <= pixel_count1[0]) && (pixel_count1[1] <= pixel_count1[2])) 
    { 
        shift1 = shift_value1[1]; 
        shift1_found = true; 
    } 
 
    if ((pixel_count2[1] <= pixel_count2[0]) && (pixel_count2[1] <= pixel_count2[2])) 
    { 
        shift2 = shift_value2[1]; 
        shift2_found = true; 
    } 
 
    //if -1 is the best shift--continue to step negative--otherwise step positive 
    if (pixel_count1[0] <= pixel_count1[1]) 
    { 
        shift1 = shift_value1[0]; 
        best_count1 = pixel_count1[0]; 
        step1 = -1; 
    } 
    else 
    { 
        shift1 = shift_value1[2]; 
        best_count1 = pixel_count1[2]; 
        step1 = 1; 
    } 
 
    if (pixel_count2[0] <= pixel_count2[1]) 
    { 
        shift2 = shift_value2[0]; 
        best_count2 = pixel_count2[0]; 
        step2 = -1; 
    } 
    else 
    { 
        shift2 = shift_value2[2]; 
        best_count2 = pixel_count2[2]; 
        step2 = 1; 
    } 
 
    //continue to step until both are done 
    while ((!shift1_found || !shift2_found) && 
            (abs (shift1) < MAX_SHIFT) && 
            (abs (shift2) < MAX_SHIFT)) 
    { 
        if (!shift1_found) 
            shift1 = shift1 + step1; 
 
        if (!shift2_found) 
            shift2 = shift2 + step2; 
 
        memcpy (sinogram1, origional_sinogram1, sizeof(float)*sinogram_x_dim*sinogram_y_dim); 
        memcpy (sinogram2, origional_sinogram2, sizeof(float)*sinogram_x_dim*sinogram_y_dim); 
 
        OffCenterCorrSingleManual (sinogram1, shift1); 
        RingCorrectionSingle (sinogram1, ring_coeff); 
 
        OffCenterCorrSingleManual (sinogram2, shift2); 
        RingCorrectionSingle (sinogram2, ring_coeff); 
 
        //reconstruct both sinograms 
        Gridrec (); 
 
        count=0; 
        for (j=0;j<(sinogram_x_dim*sinogram_x_dim);j++) 
            if (recon1[j] != 0) 
                count+=1; 
 
        if (count >= best_count1) 
        { 
            shift1 = shift1 - step1; 
            shift1_found = true; 
        } 
        else 
            best_count1 = count; 
 
        count=0; 
        for (j=0;j<(sinogram_x_dim*sinogram_x_dim);j++) 
            if (recon2[j] != 0) 
                count+=1; 
 
        if (count >= best_count2) 
        { 
            shift2 = shift2 - step2; 
            shift2_found = true; 
        } 
        else 
            best_count2 = count; 
    } 
 
    *shift_1 = shift1; 
    *shift_2 = shift2; 
} 

 
//--------------------------------------------------------------------------- 
 
CenteringClass::~CenteringClass ()
{
    if (sinogram1 != NULL)
        free (sinogram1);
    if (sinogram2 != NULL)
        free (sinogram2);

    if (recon1 != NULL)
        free (recon1);
    if (recon2 != NULL)
        free (recon2);

    if (shifted_data != NULL)
        free (shifted_data);

    if (wtemp != NULL)
        free(wtemp); 
 
    if (mean_vect != NULL)
        free (mean_vect);

    if (mean_sino_line_data != NULL)
        free (mean_sino_line_data);
    if (low_pass_sino_lines_data != NULL) 
        free (low_pass_sino_lines_data); 
 
} 

//--------------------------------------------------------------------------- 
//--------------------------------------------------------------------------- 
//Private Methods 
//--------------------------------------------------------------------------- 
//--------------------------------------------------------------------------- 
 
void CenteringClass::OffCenterCorrSingleManual (float *data, int shift) 
{ 
  int         j, m; 
 
  // performs logarithm of the projections 
  LogProj (data); 
 
  if (shift >= 0) {
    for (j=0;j<sinogram_y_dim;j++) 
      memcpy (&shifted_data[j*sinogram_x_dim+shift], &data[j*sinogram_x_dim], sizeof(float)*(sinogram_x_dim-shift)); 
  }
  else {
    for (j=0;j<sinogram_y_dim;j++) 
      memcpy (&shifted_data[j*sinogram_x_dim], &data[j*sinogram_x_dim+abs (shift)], sizeof(float)*(sinogram_x_dim-abs (shift))); 
  }
  memcpy (data, shifted_data, sizeof(float)*sinogram_x_dim*sinogram_y_dim); 
 
} 
 
//--------------------------------------------------------------------------- 
 
void CenteringClass::RingCorrectionSingle (float *data, float ring_coeff) 
{ 
  int         i, j, m; 
  float       mean_total; 
  float       tmp; 
 
    for (m=0;m<20;m++) 
    { 
        // normalization of each projection: mean values estimation 
        for (i=0;i<sinogram_y_dim;i++) 
            mean_vect[i] = 0.0; 
        mean_total = 0.0; 
 
        for (i=0;i<sinogram_y_dim;i++) 
        { 
            for (j=0;j<sinogram_x_dim;j++) 
                mean_vect[i] += data[i*sinogram_x_dim+j]; 
 
            mean_vect[i] /= sinogram_x_dim; 
            mean_total += mean_vect[i]; 
        } 
        mean_total /= sinogram_y_dim; 
 
        // renormalization of each projection to the global mean 
        for (i=0;i<sinogram_y_dim;i++) 
            for (j=0;j<sinogram_x_dim;j++) 
                if (mean_vect[i] != 0.0) 
                   data[i*sinogram_x_dim+j] = data[i*sinogram_x_dim+j]*mean_total/mean_vect[i];        // ring filtering: sum of projection and low-pass filter of the result 
 
        for (i=0;i<sinogram_x_dim;i++) 
            mean_sino_line_data[i] = 0.0; 
 
        for (i=0;i<sinogram_y_dim;i++) 
            for (j=0;j<sinogram_x_dim;j++) 
                mean_sino_line_data[j] += data[i*sinogram_x_dim+j]; 
 
        for (i=0;i<sinogram_x_dim;i++) 
            mean_sino_line_data[i] /= sinogram_y_dim; 
 
        for (j=1;j<sinogram_x_dim-1;j++) 
            low_pass_sino_lines_data[j] = (mean_sino_line_data[j-1]+mean_sino_line_data[j]+mean_sino_line_data[j+1])/3.0; 
        low_pass_sino_lines_data[0] = mean_sino_line_data[0]; 
        low_pass_sino_lines_data[sinogram_x_dim-1] = mean_sino_line_data[sinogram_x_dim-1]; 
 
        // ring corrections 
        for (i=0;i<sinogram_y_dim;i++) 
            for (j=0;j<sinogram_x_dim;j++) 
            { 
                tmp = mean_sino_line_data[j]-low_pass_sino_lines_data[j]; 
		if ((data[i*sinogram_x_dim+j] - (tmp * ring_coeff) ) > 0.0) 
		  data[i*sinogram_x_dim+j] -= (tmp * ring_coeff); 
		else 
		  data[i*sinogram_x_dim+j] = 0.0; 
            } 
    } 
 
} 
 
//--------------------------------------------------------------------------- 
 
void CenteringClass::LogProj(float *data) 
{ 
  int     i, k; 
  float   mean, max; 
  
  for (i=0;i<sinogram_y_dim;i++) {
	
    max = data[i*sinogram_x_dim]; 
    for (k=0;k<sinogram_x_dim;k++) {
      if (data[i*sinogram_x_dim+k] > max) 
	max = data[i*sinogram_x_dim+k]; 
    }

    for (k=0;k<sinogram_x_dim;k++) { 
      if (data[i*sinogram_x_dim+k] <= 0.0) 
	data[i*sinogram_x_dim+k] = 1.0; // this is really only to check if is == 0 
 
      data[i*sinogram_x_dim+k] = log (max/data[i*sinogram_x_dim+k]); 
    } 
  } 
} 
 
//---------------------------------------------------------------------------
/*
//This is a broken routine--the data is being normalized to the average value of the sum of the first and 
//last 10 pixels in each line over the entire frame.  This leads to drop outs of data that look cause
//problems in the reconstruction that look very much like missing slices.
void CenteringClass::LogProj (float *data)
{
int     i, 
        k; 
float   mean; 
 
    mean=0.0; 
    for (i=0;i<sinogram_y_dim;i++) 
        for (k=0;k<10;k++) 
            mean += data[i*sinogram_x_dim+k]+data[(i+1)*sinogram_x_dim-k-1]; 
 
    mean = mean/20/(sinogram_y_dim); 
 
    if (mean<=0.0) 
        mean=1.0; 
 
    for (i=0;i<sinogram_y_dim;i++) 
        for (k=0;k<sinogram_x_dim;k++) 
        { 
            if (data[i*sinogram_x_dim+k] > mean) 
                data[i*sinogram_x_dim+k] = mean; 
            if (data[i*sinogram_x_dim+k] <= 0.0) 
                data[i*sinogram_x_dim+k] = 1.0; 
 
            data[i*sinogram_x_dim+k] = log (mean/data[i*sinogram_x_dim+k]); 
        } 
} 
*/
//---------------------------------------------------------------------------

void CenteringClass::Gridrec ()
{
int     i, 
        j; 
float   min_image_value1=0, 
        min_image_value2 = 0, 
        max_image_value1 = -1e30, 
        max_image_value2 = -1e30; 
 
    recon_algorithm->setSinoAndReconBuffers (1, sinogram1, recon1); 
    recon_algorithm->setSinoAndReconBuffers (2, sinogram2, recon2); 
 
    recon_algorithm->reconstruct (); 
 
    for (i=0;i<sinogram_x_dim;i++) 
        for (j=0;j<sinogram_x_dim;j++) 
        { 
            if (max_image_value1 < recon1[i*sinogram_x_dim+j]) 
                max_image_value1 = recon1[i*sinogram_x_dim+j]; 
 
            if (max_image_value2 < recon2[i*sinogram_x_dim+j]) 
                max_image_value2 = recon2[i*sinogram_x_dim+j]; 
        } 
 
    for (i=0;i<sinogram_x_dim;i++) 
        for (j=0;j<sinogram_x_dim;j++) 
        { 
            if ((recon1[i*sinogram_x_dim+j]-min_image_value1) > 0.0) 
                recon1[i*sinogram_x_dim+j] = ((recon1[i*sinogram_x_dim+j]-min_image_value1)/(max_image_value1-min_image_value1)*65535); 
            else 
                recon1[i*sinogram_x_dim+j] = 0; 
 
            if ((recon2[i*sinogram_x_dim+j]-min_image_value2) > 0.0) 
                recon2[i*sinogram_x_dim+j] = ((recon2[i*sinogram_x_dim+j]-min_image_value2)/(max_image_value2-min_image_value2)*65535); 
            else 
                recon2[i*sinogram_x_dim+j] = 0; 
        } 
} 
 
//---------------------------------------------------------------------------
