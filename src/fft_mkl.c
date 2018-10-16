/*** File fft_intel.c 5/5/2000
     This file calls routines from the Intel Math Kernel Library.  It
	 emulates the 1-D and n-D routines from Numerical Recipes, so that
     gridrec requires essentially no modification

    Written by Mark Rivers
 **/

#include <stdlib.h>
#include <stdio.h>
//#include <mkl_fftc_ln.h>
//#include <mkl_fft.h>

static float    *wsave,
                *real_data_1D,
	            *imaginary_data_1D,
                *real_data_ND,
	            *imaginary_data_ND;

//_______________________________________________________________________________________________________________

void  MCFFT1D (float data[], unsigned long n, int zero, float *wsave)
{
int	    loop;

    if (real_data_1D == NULL)
    	real_data_1D = (float *) malloc (sizeof (float) * n);
    if (imaginary_data_1D == NULL)
    	imaginary_data_1D = (float *) malloc (sizeof (float) * n);

	for (loop=0;loop<n;loop++)
	{
		real_data_1D[loop] = data[loop*2];
		imaginary_data_1D[loop] = data[loop*2+1];
	}

	cfft1dc (real_data_1D, imaginary_data_1D, n, zero, wsave);

	for (loop=0;loop<n;loop++)
	{
		data[loop*2] = real_data_1D[loop];
		data[loop*2+1] = imaginary_data_1D[loop];
	}
}

//_______________________________________________________________________________________________________________

void  MCFFT2D (float data[], unsigned long nx, unsigned long ny, int isign)
{
int	    loop;

    if (real_data_ND == NULL)
    	real_data_ND = (float *) malloc (sizeof (float) * nx * ny);
    if (imaginary_data_ND == NULL)
    	imaginary_data_ND = (float *) malloc (sizeof (float) * nx * ny);

	for (loop=0;loop<nx*ny;loop++)
	{
		real_data_ND[loop] = data[loop*2];
		imaginary_data_ND[loop] = data[loop*2+1];
	}

	cfft2dc (real_data_ND, imaginary_data_ND, nx, ny, isign);

	for (loop=0;loop<nx*ny;loop++)
	{
		data[loop*2] = real_data_ND[loop];
		data[loop*2+1] = imaginary_data_ND[loop];
	}
}

//_______________________________________________________________________________________________________________

void initFFTMemoryStructures (void)
{
    wsave = NULL;

    real_data_1D = NULL;
    imaginary_data_1D = NULL;

    real_data_ND = NULL;
    imaginary_data_ND = NULL;
}

//_______________________________________________________________________________________________________________

void destroyFFTMemoryStructures (void)
{
char msg[256];
FILE *output;

    if (wsave != NULL)
      	free(wsave);

    if (real_data_1D != NULL)
    	free (real_data_1D);
    if (imaginary_data_1D != NULL)
    	free (imaginary_data_1D);

    if (real_data_ND != NULL)
    	free (real_data_ND);
    if (imaginary_data_ND != NULL)
    	free (imaginary_data_ND);


output = fopen ("output.log", "wt");
fprintf (output, "Done!\n");
fclose (output);
}

//_______________________________________________________________________________________________________________

void four1(float data[], unsigned long nn, int isign)
{
   static int n_prev;
   float scale, *p;
   int n = nn;
   int i;
   int zero = 0;

   /* Call the Intel Math Kernel Library routine */
	if ((isign == 0) || (n != n_prev))
    {
   	    /* The required storage is (3*N)/2 complex elements = 3*N floats */
        if (wsave != NULL);
            free (wsave);
        wsave = malloc(3 * n * sizeof(float));

		n_prev = n;
		MCFFT1D(&data[1], n, zero, wsave);
	}

	/* The Numerical Recipes routines are passed a pointer to one element
	   before the start of the array - add one */
	MCFFT1D(&data[1], n, isign, wsave);

	/* Must rescale data if isign is 1, since Numerical Recipes output is scaled by N, CFFT1D is not */
	if (isign == 1)
    {
		scale = (float)n;
	    for (i=0, p=data+1; i<2*n; i++)
            *p++ *= scale;
	}

}

//_______________________________________________________________________________________________________________

void fourn(float data[], unsigned long nn[], int ndim, int isign)
{

	int nx = nn[2];
	int ny = nn[1];
	int i;
	float scale, *p;

	/* NOTE: This function only works for ndim=2 */
	if (ndim != 2) return;
	/* Call Intel Math Kernel Library routine */
	/* The Numerical Recipes routines are passed a pointer to one element
	   before the start of the array - add one */
	MCFFT2D(&data[1], nx, ny, isign);

	/* Must rescale data if isign is 1, since Numerical Recipes output is scaled by N, CFFT2D is not */
	if (isign == 1)
    {
		scale = (float)(nx*ny);
	    for (i=0, p=data+1; i<2*nx*ny; i++)
            *p++ *= scale;
	}
}

//_______________________________________________________________________________________________________________

