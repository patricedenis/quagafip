#ifndef _FOURIER_QUATERNION_H_
#define _FOURIER_QUATERNION_H_

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

//fourier quaternion
#include "quaternion.h"
#include "matrix_double.h"
#include "common.h"
#include "couleur.h"

//fourier quaternion rapide
#include "complex.h"
#include "MinMax.h"
#include "matrix.h"
#include "matrix_double.h"

#include "chaine.h"
#include "MagickIO.h"


void Quaternion_Matrix_Fourier_Transform(Quaternion ** QMatrix,Quaternion *** pQMatrixFT,
    Quaternion QNu,int intHeight,int intWidth);
void Quaternion_Matrix_Fourier_Transform_ExpLeft(Quaternion ** QMatrix,Quaternion *** pQMatrixFT,
    Quaternion QNu,int intHeight,int intWidth);
void Quaternion_Init_Freq(int intHeight,int intWidth,int intNbStripe, Quaternion QNu,double K,char * strFile);
void QFFT(Quaternion **QIn,int * dim,Quaternion *** pQOut,Quaternion mu1,Quaternion mu2,int intSens);
int symetrie_partieR_OK(Quaternion ** QMat,int *dim,double dblPrecision,double dblSeuilTolerance);
int symetrie_partieImI_OK(Quaternion ** QMat,int *dim,double dblPrecision,double dblSeuilTolerance);
void QFFT_Verifie_Symetries(char * strFile, Quaternion QMu1,Quaternion QMu2);
int main_fourier_quaternion(int argc, char *argv[]);

#endif

