/* tomosolo.cpp
 *
 * Basic example of how to call the tomography reconstruction classes.
 *
 * Author: Brian Tieman
 * Created: 9/8/2009
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "recon_algorithm.h"
#include "filteredbackprojection.h"
#include "gridrec.h"

int				sinogram_dimx,
				sinogram_dimy,
				num_sinogram_elements,
				num_thetas;
unsigned short	*sinogram;
float			*normalized_sinogram,
				*reconstruction,
	 			*theta_list;

ReconAlgorithm  *recon_algorithm;

//_____________________________________________________________________________________

/*
 * Input data is already formated as a "sinogram".  A sinogram is a 2 dimensional
 * array where each row is a single line from the projected image.  Each row in the
 * sinogram image is the same line taken from all the projections in the angle order
 * of the thetaList array.  IE: sinogram1 would be the first line from each projection.
 * Row 1 would be the first like from projection angle 0, row 2 the first line from
 * projection angle 0.25, row 3 from projection angle 0.5, etc...
 *
 * The angle list is in degrees and ordered numerically.
 *
 */

void InputData (char *sinogramFile, char *thetaListFile)
{
FILE	*inputFile;
char	fileName[256];
int		recordsRead;

	//First read sinogram
	//Sinogram is stored as x<int>, y<int>, elements<int>, data<short int>
	inputFile = fopen (sinogramFile, "rb");

	recordsRead = fread ((void *) &sinogram_dimx, sizeof (int), 1, inputFile);
	recordsRead = fread ((void *) &sinogram_dimy, sizeof (int), 1, inputFile);
	recordsRead = fread ((void *) &num_sinogram_elements, sizeof (int), 1, inputFile);

	sinogram = (unsigned short *) malloc (sizeof (unsigned short) * num_sinogram_elements);
	normalized_sinogram = (float *) malloc (sizeof (float) * num_sinogram_elements);
	reconstruction = (float *) malloc (sizeof (float) * sinogram_dimx * sinogram_dimx);
	recordsRead = fread ((void *) sinogram, sizeof (unsigned short), num_sinogram_elements, inputFile);

	fclose (inputFile);

	//Then read the angle list
	//Angle list is stored as num_angles<int>, angles<float>
	inputFile = fopen (thetaListFile, "rb");

	recordsRead = fread ((void *) &num_thetas, sizeof (int), 1, inputFile);

	theta_list = (float *) malloc (sizeof (float) * num_thetas);
	recordsRead = fread ((void *) theta_list, sizeof (float), num_thetas, inputFile);

	fclose (inputFile);

}

//_____________________________________________________________________________________

/*
 * The reconstruction class defines an API that all reconstruction algorithms use.
 * This routine specifically initializes a Filtered Back Projection instance of the
 * class. The FBP class only reconstructs 1 sinogram at a time.  The optimization
 * levels are defined as follows:
 *
 * FBP () -- no optimization
 * OptimizedFBP () -- hand optimized version of the code.  Harder to read, but faster.
 * CircleFBP () -- same code as the optimized version but ignores the "corners" of the image
 *			that don't contribute to the quality of the reconstruction
 *
 * Theta list is in degrees.
 *
 * Filters are applied to the data during the reconstruction process.  Valid values are:
 *
 * #define         FILTER_NONE		               0
 * #define         FILTER_SHEPP_LOGAN              1
 * #define         FILTER_HANN                     2
 * #define         FILTER_HAMMING                  3
 * #define         FILTER_RAMP                     4
 * #define         FILTER_FBP                      5
 *
 * and are defined in recon_algorithm.h
 *
 */

void InitFBP (int optimizationLevel)
{

	switch (optimizationLevel)
	{
		case 1 : recon_algorithm = new FBP (); break;
		case 2 : recon_algorithm = new OptimizedFBP (); break;
		case 3 : recon_algorithm = new CircleFBP (); break;
		default : recon_algorithm = new OptimizedFBP (); break;
	}

	recon_algorithm->setSinogramDimensions(sinogram_dimx, sinogram_dimy);
	recon_algorithm->setThetaList (theta_list, num_thetas);
	recon_algorithm->setFilter (FILTER_NONE);

	recon_algorithm->init();

}

//_____________________________________________________________________________________

/*
 * The reconstruction class defines an API that all reconstruction algorithms use.
 * This routine specifically initializes a Gridrec instance of the class.  Gridrec
 * is capable of reconstructing 2 singrams at the same time.
 *
 * Theta list is in degrees.
 *
 * Filters are applied to the data during the reconstruction process.  Valid values are:
 *
 * #define         FILTER_NONE		               0
 * #define         FILTER_SHEPP_LOGAN              1
 * #define         FILTER_HANN                     2
 * #define         FILTER_HAMMING                  3
 * #define         FILTER_RAMP                     4
 * #define         FILTER_FBP                      5
 *
 * and are defined in recon_algorithm.h
 *
 */

void InitGridrec (void)
{
int				loop;

	recon_algorithm = new GridRec ();

	recon_algorithm->setSinogramDimensions(sinogram_dimx, sinogram_dimy);
	recon_algorithm->setThetaList (theta_list, num_thetas);
	recon_algorithm->setFilter (FILTER_NONE);

	recon_algorithm->init();

}

//_____________________________________________________________________________________

/* This is where normalization would be done.  In the parallel code, the normalization
 * is done on the sinograms in the client nodes for performance reasons.  Normalization
 * on the sinogram is tricky because different backgrounds/write fields may need to
 * be applied to different rows in the sinogram.  If possible, it is generally best to
 * normalize the projections prior to making sinograms.
 *
 * This version of "normalize" performas no normalization and simply converts the data
 * to floating point for the reconstruction class.
 */

void Normalize (unsigned short *short_sino, float *norm_sino)
{

	for (int loopx=0;loopx<sinogram_dimx;loopx++)
		for (int loopy=0;loopy<sinogram_dimy;loopy++)
			norm_sino[loopy*sinogram_dimx+loopx] = 1.0 * short_sino[loopy*sinogram_dimx+loopx];

}

//---------------------------------------------------------------------------

/* Takes the log of the data.
 */

void LogSinogram (float *data)
{
int     i,
        k;

    for (i=0;i<sinogram_dimy;i++)
        for (k=0;k<sinogram_dimx;k++)
		{
			if (data[i*sinogram_dimx+k] > 0)
				data[i*sinogram_dimx+k] = -1 * log (data[i*sinogram_dimx+k]);
			else
				data[i*sinogram_dimx+k] = 0;
		}

}

//_____________________________________________________________________________________

/* Clean up our mess...
 */

void CleanUp ()
{
	if (sinogram != NULL)
		free (sinogram);

	if (normalized_sinogram != NULL)
		free (normalized_sinogram);

	if (reconstruction != NULL)
		free (reconstruction);

	if (theta_list != NULL)
		free (theta_list);

	if (recon_algorithm != NULL)
		delete (recon_algorithm);

}

//_____________________________________________________________________________________

/* Be nice to the poor user--tell them what they are doing wrong.
 */

void PrintUsage ()
{
	printf ("Usage:\n");
	printf (">tomosolo <full_path_to_sinogram_file> <full_path_to_theta_file> <fbp|gridrec>\n\n");
	exit (0);
}

//_____________________________________________________________________________________

int main (int argc, char* argv[])
{
char	sinogramFile[256],
		thetaListFile[256],
		reconstructionFile[256];
FILE	*outputFile;

	sinogram = NULL;
	normalized_sinogram = NULL;
	reconstruction = NULL;
	theta_list = NULL;
	recon_algorithm = NULL;

	if (argc != 4)
		PrintUsage ();
	if ((strcmp (argv[3], "fbp") != 0) && (strcmp (argv[3], "gridrec") != 0))
		PrintUsage ();

    strcpy (sinogramFile, argv[1]);
    strcpy (thetaListFile, argv[2]);

	InputData (sinogramFile, thetaListFile);

	if (strcmp (argv[3], "fbp") == 0)
	{
		printf ("Using Filtered Back Projection...\n");
		InitFBP(1);
	}
	if (strcmp (argv[3], "gridrec") == 0)
	{
		printf ("Using gridrec...\n");
		InitGridrec();
	}

	Normalize (sinogram, normalized_sinogram);
	LogSinogram (normalized_sinogram);

	//recon_algorithm->setSinoAndReconBuffers(slot, normalized_sinogram, reconstruction);
	//slot = slot to put the sinogram in--gridrec has slots 1 and 2, FBP has slot 1 only
	//normalized_sinogram = normalized sinogram
	//reconstruction = pre-allocated array of size sinogram_dimx x sinogram_dimx x sizeof (float)
	//When reconstruction is done, the memory pointer reconstruction will contain the reconstructed frame.

	//if gridrec, we can work on 2 sinograms at a time--we'll just feed it the same one twice for demonstration
	//else if fbp we can use a single sinogram
	if (strcmp (argv[3], "gridrec") == 0)
	{
		recon_algorithm->setSinoAndReconBuffers(1, normalized_sinogram, reconstruction);
		recon_algorithm->setSinoAndReconBuffers(2, normalized_sinogram, reconstruction);
	}
	else
		recon_algorithm->setSinoAndReconBuffers(1, normalized_sinogram, reconstruction);


	recon_algorithm->reconstruct ();


	//write the reconstructed image
	sprintf (reconstructionFile, "./reconstructed.bin");
	outputFile = fopen (reconstructionFile, "wb");

	fwrite ((void *) &sinogram_dimx, sizeof (int), 1, outputFile);
	fwrite ((void *) &sinogram_dimx, sizeof (int), 1, outputFile);
	fwrite ((void *) reconstruction, sizeof (float), sinogram_dimx*sinogram_dimx, outputFile);

	fclose (outputFile);

	CleanUp ();

    return 0;
}

//_____________________________________________________________________________________
