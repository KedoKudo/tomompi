#ifndef gridrecH
#define gridrecH
//---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h> 
#include <string.h> 
#include <math.h> 
#include <stddef.h> 
#include <time.h> 
#include <sys/stat.h> 

//---------------------------------------------------------------------------

#include "recon_algorithm.h"

//---------------------------------------------------------------------------

/**** Macros and typedefs ****/
#define max(A,B) ((A)>(B)?(A):(B)) 
 
#define min(A,B) ((A)<(B)?(A):(B)) 
 
//#define free_matrix(A) (free(*(A)),free(A)) 
 
#define abs(A) ((A)>0 ?(A):-(A)) 
 
#define Cmult(A,B,C) {(A).r=(B).r*(C).r-(B).i*(C).i; (A).i=(B).r*(C).i+(B).i*(C).r;} /** A=B*C for complex A, B, C. A must be distinct, and an lvalue */ 
 
#ifdef INTERP 
#define Cnvlvnt(X) (wtbl[(int)X]+(X-(int)X)*dwtbl[(int)X]) /* Linear interpolation version */ 
#else 
#define Cnvlvnt(X) (wtbl[(int)(X+0.5)]) /* Nearest nbr version - no interpolation */ 
#endif 
 
#define TOLERANCE 0.1	/* For comparing centers of two sinugrams */ 
#define LTBL_DEF 512	/* Default lookup table length */ 
 
//---------------------------------------------------------------------------
#define NO_PSWFS 5
 
typedef struct PSWF_STRUCT { 
/*Prolate spheroidal wave fcn (PSWF) data */ 
    float   C,          /* Parameter for particular 0th order pswf being used*/ 
            lmbda;   	/* Eigenvalue */ 
    int     nt;         /* Degree of Legendre polynomial expansion */ 
    float   coefs[15];	/* Coeffs for Legendre polynomial expansion */ 
} pswf_struct; 
 
//--------------------------------------------------------------------------- 

class GridRec : public ReconAlgorithm
{
public:

    GridRec (void);

    void setSinoAndReconBuffers (int   number, 
                                float *sinogram_address, 
                                float *reconstruction_address);

    void init (void);
    void reconstruct (void);
    void destroy (void);

    static void acknowledgements (LogFileClass *acknowledge_file);

private:
    int     flag;

    long    pdim, 
            M, 
            M0, 
            M02, 
            ltbl, 
            imgsiz; 

    float   sampl, 
            scale, 
            L, 
            X0, 
            Y0; 

    float   *SINE, 
            *COSE, 
            *wtbl, 
            *dwtbl, 
            *work, 
            *winv;

    complex *cproj, 
            *filphase, 
            *H;

    float   **G1, 
            **G2, 
            **S1, 
            **S2;

    pswf_struct pswf_db[NO_PSWFS]; 
 
    void phase1 (void);
    void phase2 (void);
    void phase3 (void);

    void trig_su (int geom, int n_ang);
    void filphase_su (long pd, float center, complex *A);
    void pswf_su (pswf_struct *pswf, long ltbl, long linv, float* wtbl, float* dwtbl, float* winv);
	float legendre (int n, float *coefs, float x);
    void get_pswf (float C, pswf_struct **P);

};

//---------------------------------------------------------------------------
#endif
