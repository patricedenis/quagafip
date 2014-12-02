#ifndef _COMPLEX_H_
#define _COMPLEX_H_

#include "bool.h"
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "MinMax.h"

struct stComplex
{
  double r;
  double i;
};

typedef struct stComplex Complex;

///////////////////////////////////////////////////////////
//       OPERATIONS SUR COMPLEXES SIMPLES
///////////////////////////////////////////////////////////

Complex Complex_InitXY(double x,double y);
double CReal(Complex cpx);
double CIm(Complex cpx);
Complex Complex_InitRoTheta(double dblRo,double dblTheta);
Complex Complex_InitExpTheta(double dblTheta);
Complex Complex_Sum(const Complex c1,const Complex c2);
Complex Complex_Diff(const Complex c1,const Complex c2);
Complex Complex_Mult(const Complex c1,const Complex c2);
Complex Complex_Scal_Mult(const Complex c,const double d);
Complex Complex_Div(const Complex c1,const Complex c2);
Complex Complex_Conjugate(const Complex c);
double Complex_Modulus(const Complex c);
double Complex_Phase(const Complex c);
double Complex_Phase2(const Complex c);
void Complex_Display(const Complex c,FILE * stream);

///////////////////////////////////////////////////////////
//       OPERATIONS SUR MATRICES COMPLEXES
///////////////////////////////////////////////////////////


int Complex_Matrix_Allocate(int intHeight,int intWidth,Complex *** pcpxMatrix);
void Complex_Matrix_Free(int intHeight,Complex *** cpxMatrix);
int Get_Complex_Matrix_From_Image_23(double ** MatYUV,int intHeight,int intWidth,Complex *** pcpxMatrix);
void Complex_Init_Matrix_Double(int intHeight,int intWidth,const double dblValue, Complex *** pcpxMatrix);
void CMatrixCopy(const Complex ** cpxIn, Complex *** pcpxOut, int * dim);
void Complex_Init_Matrix(int intHeight,int intWidth,const Complex cpxValue, Complex *** cpxMatrix);
void Set_Complex_Matrix_To_Image_23(double *** MatYUV,int intHeight,int intWidth,const Complex ** cpxMatrix);
void Complex_Copy_SubMatrix(const Complex ** MatFrom,Complex *** MatTo,int intXFrom,int intYFrom,
    int intXTo,int intYTo,int intHeight,int intWidth);
void Complex_Matrix_Shift(const Complex ** MatFrom,Complex *** pMatShifted,int intHeight,int intWidth);
int CSetRMatrix(double tabRSet [],int intHeight,int intWidth,Complex *** pCMatrix);
void CSetMatrix(double tabRSet [],double tabISet [],int intHeight,int intWidth,	Complex *** pCMatrix);
void CSetMatrixFromDblTab(double * dlbTabRealPart, double * dblTabImaginaryPart,int * dim,Complex *** pcpxMat);
void CSetMatrixFromDblMats(double ** dlbMatRealPart, double ** dblMatImaginaryPart,int * dim,Complex *** pcpxMat);
void CGetMatrixDoubleTab(double * pr[],double * pi[],int intHeight,int intWidth,const Complex ** CMatrix);
void CGetMatrixTab(Complex * pcpx[],int intHeight,int intWidth,const Complex ** CMatrix);
void CMatrixMod(double * dblTabModulus [],int intHeight,int intWidth,const Complex ** CMat);
int IsCMatRealInt(int * intDim,const Complex ** cpxMat);
int IsCMatImInt(int * intDim,const Complex ** cpxMat);
int IsCMatRealEpsilon(int * intDim,const Complex ** cpxMat,double dblEpsilon);
void CMatDisp(const Complex **CMat,int * intDim,FILE * stream);
void CMatMultTBT(const Complex ** CMatIn1,const Complex ** CMatIn2, Complex *** pCMatOut,int * dim);
void CMatConj(const Complex ** CMatIn, int * dim,	Complex *** pCMatConj);
void CMatScalMult(const Complex ** CMatIn, int * dim,	Complex *** pCMatOut,double dblScal);


void CSubImageUnion(const Complex ** cpxll,const Complex ** cpxlh,const Complex ** cpxhl,
  const Complex ** cpxhh, int *dim,Complex *** pcpxComp); 
void CThumbnailsConstruct(Complex *** pcpxMatOut,int * intDimOut, const Complex **cpxLow,
 const Complex **cpxHighV,const Complex **cpxHighH,const Complex **cpxHighD);
void CAdd4(Complex *** pcpxCMat,int *dim,Complex ** cpxCLow,
  Complex ** cpxCHighV,Complex ** cpxCHighH,Complex ** cpxCHighD);
void CThumbnailsExtract(const Complex ** cpxMat,int * intDim, Complex ***pcpxLow,
 Complex ***pcpxHighV,Complex ***pcpxHighH,Complex ***pcpxHighD);
void CNThumbnailsExtract(const Complex **cpxMat,int *intDim,int intN,Complex ***pcpxll,Complex ***pcpxlh,
Complex ***pcpxhl,Complex ***pcpxhh);
void CThumbnailsConstructWithoutLow(Complex *** pcpxMatOut,int * intDimOut,
 const Complex **cpxHighV,const Complex **cpxHighH,const Complex **cpxHighD);
void CThumbnailsConstructLow(Complex *** pcpxMatOut,int * intDimLow, const Complex **cpxLow);
void CMatCopy(Complex ***pcpxMatCpy,int *intDim,const Complex **cpxMat);
void CThumbnailsExtractLow(const Complex ** cpxMat,int * intDim, Complex ***pcpxLow);
void CMatAdd4(Complex *** pcpxCMat,int *dim,const Complex ** cpxCLow,
  const Complex ** cpxCHighV,const Complex ** cpxCHighH,const Complex ** cpxCHighD);
//void CMatScaleIm0255(const Complex **cpxMatIn,int *dim, Complex *** pcpxMatOut);
void CMatScaleRe0255(const Complex** cpxlh,const Complex** cpxhl,const Complex** cpxhh,
  int *dim, Complex *** pcpxlhOut,Complex *** pcpxhlOut,Complex *** pcpxhhOut);
void CMat_Min_Max(const Complex ** CMat,int * dim,double * dblRMin,double * dblRMax,
double * dblIMin,double * dblIMax);
void CMat_Min_Max2(const Complex ** CMatPure,int * dim,double * dblMin,double * dblMax);
int CMat_ChgScale_IntLvlMax(const Complex ** CMat, int * dim,int intLvlMax,
	double dblRMin,double dblRMax,double dblIMin,double dblIMax,Complex *** pCMatScaled);
int CReMat_ChgScale_IntLvlMax(const Complex ** CMat, int * dim,int intLvlMax,
	double dblIMin,double dblIMax,Complex *** pCMatScaled);
int CMat_ChgScale_IntLvlMax_PerChannel(const Complex ** CMat, int * dim,int intLvlMax,
	double dblRMin,double dblRMax,double dblIMin,double dblIMax,Complex *** pCMatScaled);

///////////////////////////////////////////////////////////
//       OPERATIONS SUR VECTEURS COMPLEXES
///////////////////////////////////////////////////////////

void CVectMultTBT(const Complex * CVect1,const Complex * CVect2,Complex ** pCMult,int dim);
int CVectAlloc(Complex **pCVect,int dim);
void CVectFree(Complex **pCVect);
void CInitVect(Complex **pCVect,int dim,Complex cpxValue);
void CVectScalMult(Complex **pVect,int dim,double dblScal);
void CVectConj(Complex **pVect,Complex *Vect,int dim);
void CVectDisp(const Complex *Vect,int dim,FILE * stream);
void CGetMatIRow(const Complex ** Mat,Complex **pVect,int * dim,int n);
void CSetMatIRow(Complex *** pMat,const Complex *Vect,int * dim,int n);
void CGetMatIColumn(const Complex ** Mat,Complex **pVect,int * dim,int n);
void CSetMatIColumn(Complex *** pMat,const Complex *Vect,int * dim,int n);

#endif
