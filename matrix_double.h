#ifndef _MATRIX_DOUBLE_H
#define _MATRIX_DOUBLE_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "bool.h"
#include "MinMax.h"
#include "MagickIO.h"

int Matrix_Allocate_Double(int intHeight,int intWidth,double *** pMatrix);
void Matrix_Free_Double(int intHeight,double *** pMatrix);
void Matrix_Init_Double(double dblValue,int *intDim,double *** pdblMatrix);
//cette procedure initialise la matrice avec des valeurs nulles
int DblMatInitZero(double ***pdblMat, int * dim);

void DblLogMatrix1D(const double * dblData, int * dim,double ** dblLogData);
void DblLogMatrix(const double ** dblMat, int * dim,double *** pdblLogMat);
void DblLogTab(const double * dblTab, int dim,double ** pdblLogTab);
void DblLogP1Matrix(const double ** dblMat, int * dim,double *** pdblLogMat);
void Dbl_Matrix_Min_Max2(const double * dblMat,int * dim,double * dblMin,double * dblMax);
void DblMatMaxed(double *** pdblMatMaxed,const double ** dblMat2,int * dim);
void DblMatMixed(double *** pdblMatMaxed,const double ** dblMat2,double dblCoeffDiv, int * dim);
int Dbl_Mat_ChgScale_IntLvlMax(const double ** dblMat, int * dim,int intLvlMax,
	double dblMin,double dblMax, double *** pdblMatScaled);
int Dbl_Mat_ChgScaleLog_IntLvlMax(const double ** dblMat, int * dim,int intLvlMax,
	double dblMin,double dblMax, double *** dblMatScaled);
int Dbl_Mat_ChgScale_IntLvlMax2(const double * dblMat, int * dim,int intLvlMax,
	double dblMin,double dblMax, double ** pdblMatScaled);

void Dbl_Tab_Min_Max(const double * dblTab,int dim,double * dblMin,double * dblMax);
int Dbl_Tab_ChgScale_IntLvlMax(const double * dblTab, int dim,int intLvlMax,
	double dblMin,double dblMax, double ** pdblTabScaled);
    
void SetDblMatrixFromTab(double tabSet [],int intHeight,int intWidth,double *** pMatrix);
double VectMean_Double(double * tab, int * dim);
int IsVectDoubleNullEpsilon(double * tab, int dim,double dblEpsilon);
int DblIsMatNull(double ** dblMat, int * intDim);
void DblImageThreshold(const double ** dblMatFrom, int * intDim,
  double *** pdblMatTo, double dblThreshold, int bln255);
void DblImageProductTBT(const double ** dblMatFrom1,const double ** dblMatFrom2,
  int * intDim, double *** pdblMatTo);
void DblImageInverse(const double ** dblMatFrom, int * intDim,
  double *** pdblMatTo);

void DblFiltrePrewittH(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim);
void DblFiltrePrewittV(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim);
void DblFiltrePrewittD1(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim);
void DblFiltrePrewittD2(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim);

void DblMatL2Normalization(const double ** dblMatH1,const double ** dblMatH2,
  const double ** dblMatV1,const double ** dblMatV2,
  const double ** dblMatD1,const double ** dblMatD2,
  double ***pdblMatNormalized,int * intDim);
void DblMatL2Norme(double *** pdblMatDistMaxL2,const double ** dblMatDistH,const double ** dblMatDistV,
  const double ** dblMatDistD1,const double ** dblMatDistD2,int * intDim);
  
void DblFiltrePrewittTotal(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim);

void DblGaussienne2D(double dblAmp, double dblSigma, double * dblCoordCentre, int * intDim, double *** pDblTab,int blnInvert);
void DblCopySubMatrix(double ** DblMatFrom,double *** pDblMatTo,int intXFrom,int intYFrom,
  int intXTo,int intYTo,int intHeight,int intWidth);
void DblMatrixShift(double ** DblMatFrom,double *** pDblMatShifted,int * dim);
	
#endif
