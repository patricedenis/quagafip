#ifndef _QUATERNION_H_
#define _QUATERNION_H_

#include <stdio.h>
#include <math.h>
#include <float.h>

#include "bool.h"
#include "MinMax.h"

struct stQuaternion
{
  double a;
  double b;
  double c;
  double d;
};

typedef struct stQuaternion Quaternion;

//Quaternion Base Functions

Quaternion QInit(const double a,const double b,const double c,const double d);
Quaternion QInitZero();

Quaternion QConj(const Quaternion q1);

Quaternion QOpp(const Quaternion q1);

double QReal(const Quaternion q);

double QImI(const Quaternion q);

double QImJ(const Quaternion q);

double QImK(const Quaternion q);

Quaternion QInitExp(const Quaternion Nu,const double theta);

Quaternion QImag(const Quaternion q);

double QNorm(const Quaternion q);

Quaternion QInv(const Quaternion q);

Quaternion QDiv(const Quaternion q1, const Quaternion q2);

void QDisp(const Quaternion q,FILE * stream);

Quaternion QAdd(const Quaternion q1,const Quaternion q2);

Quaternion QScalAdd(const Quaternion q1,const double dbl1);

Quaternion QDiff(const Quaternion q1,const Quaternion q2);

Quaternion QMult(const Quaternion q1,const Quaternion q2);

double QScalProd(const Quaternion q1,const Quaternion q2);

Quaternion QScalMult(const Quaternion q1,const double dbl1);

Quaternion QScalDiv(const Quaternion q1,const double dbl1);

void QExp_LogRo_Mu_Phi(const Quaternion q, double * LogRo, Quaternion * qMu, double * Phi);

void QParOrtho(const Quaternion q,const Quaternion p, Quaternion * qPar, Quaternion * qOrtho);

//Quaternion Matrix Functions

int QMatrixAllocate(int intHeight,int intWidth,Quaternion *** pQMatrix);

void QMatrixFree(int intHeight,Quaternion *** pQMatrix);

void QMatrixDisp(const Quaternion **qMat,int intHeight,int intWidth,int intPart,FILE * stream);

//cette procedure initialise la matrice avec des valeurs nulles
int QMatInitZero(Quaternion ***pqtnMat, int * dim);

//cette procedure initialise la matrice avec les tableaux lineaires de doubles passés en paramatre
int QSetMatrix(double r[],double g[],double b[],int intHeight,int intWidth,Quaternion *** pQMatrix);
// cette procedure remplit la matrice avec les matrices doubles passées en parametres
int QSetMatrixFromDblMat(const double ** dblCompR,const double ** dblCompI, const double ** dblCompJ,
  const double ** dblCompK,int * intDim, Quaternion *** pQMatrix);


// this procedure take the composant from r, i, j and k arrays
// to fill the Quaternion Matrix QMatrix.
int QSetMatrixComp(double r [],double i[],double j[],double k[],int intHeight,int intWidth,
	Quaternion *** pQMatrix);

void QGetMatrixImagPart(double * i[],double * j[],double * k[],int intHeight,int intWidth,
	Quaternion ** QMatrix);

void QGetMatrixCompl(double * mr[],double * mi[],double * mj[],double * mk[],int intHeight,int intWidth,
	Quaternion ** QMatrix);	

void QInitMatrix(int intHeight,int intWidth, Quaternion QInitValue, Quaternion *** pQMatrix);

void QMatrixShift(Quaternion ** QMatFrom,Quaternion *** pQMatShifted,int intHeight,int intWidth);

void QGetMatrixComponent(double * comp[],int intHeight,int intWidth,Quaternion ** QMatrix,char * strComp);	

void QSetMatrixComponent(double  comp[],int intHeight,int intWidth,Quaternion *** pQMatrix,char * strComp);

void QCopyMatrixPartR(const Quaternion ** MatrixFrom , Quaternion *** pQMatrixTo, int * intDim);
void QCopyMatrixPartI(const Quaternion ** MatrixFrom , Quaternion *** pQMatrixTo, int * intDim);
void QCopyMatrixPartJ(const Quaternion ** MatrixFrom , Quaternion *** pQMatrixTo, int * intDim);
void QCopyMatrixPartK(const Quaternion ** MatrixFrom , Quaternion *** pQMatrixTo, int * intDim);

int IsQMatImEpsilon(int * intDim,Quaternion ** qtnMat,double dblEpsilon);
int IsQMatReEpsilon(int * intDim,Quaternion ** qtnMat,double dblEpsilon);

void QMatDisp(const Quaternion **QMat,int * intDim,FILE * stream);

void QMatMultTBT(Quaternion ** QMatIn1, Quaternion ** QMatIn2, Quaternion *** pQMatOut,int * dim);

//cette fonction effectue la multiplication terme par terme avec une matrice de reels
void QMatMultTBTWithDblMat(Quaternion ** QMatIn1, double ** dblMatIn2, Quaternion *** pQMatOut,int * dim);

void QMatConj(Quaternion ** QMatIn, int * dim,	Quaternion *** pQMatConj);

void QMatScalMult(Quaternion ** QMatIn, int * dim,	Quaternion *** pQMatOut,double dblScal);

void QSubImageUnion(Quaternion ** qtnll,Quaternion ** qtnlh,Quaternion ** qtnhl,Quaternion ** qtnhh,
  int *dim,Quaternion *** pqtnComp); 

void QThumbnailsConstruct(Quaternion *** pqtnMatOut,int * intDimOut, Quaternion **qtnLow,
 Quaternion **qtnHighV,Quaternion **qtnHighH,Quaternion **qtnHighD);

void QThumbnailsExtract(const Quaternion ** qtnMat,int * intDim, Quaternion ***pqtnLow,
 Quaternion ***pqtnHighV,Quaternion ***pqtnHighH,Quaternion ***pqtnHighD);

void QNThumbnailsExtract(const Quaternion **qtnMat,int *intDim,int intN,Quaternion ***pqtnll,Quaternion ***pqtnlh,
Quaternion ***pqtnhl,Quaternion ***pqtnhh);

void QThumbnailsConstructWithoutLow(Quaternion *** pqtnMatOut,int * intDimOut,
 Quaternion **qtnHighV,Quaternion **qtnHighH,Quaternion **qtnHighD);

void QThumbnailsConstructLow(Quaternion *** pqtnMatOut,int * intDimLow, Quaternion **qtnLow);

void QMatCopy(Quaternion ***pqtnMatCpy,int *intDim,const Quaternion **qtnMat);

void QThumbnailsExtractLow(Quaternion ** qtnMat,int * intDim, Quaternion ***pqtnLow);

//cette fonction somme les deux matrices passées en parametre
void QMatAdd(Quaternion *** pqtnMatAdd,Quaternion ** qtnMat1,Quaternion ** qtnMat2,int *dim);

void QAdd4(Quaternion *** pqtnQMat,int *dim,Quaternion ** qtnQLow,
  Quaternion ** qtnQHighV,Quaternion ** qtnQHighH,Quaternion ** qtnQHighD);

//void QMatScaleIm0255(const Quaternion **qtnMatIn,int *dim, Quaternion *** pqtnMatOut);
void QMatScaleIm0255(Quaternion** qtnlh,Quaternion** qtnhl,Quaternion** qtnhh,
  int *dim, Quaternion *** pqtnlhOut,Quaternion *** pqtnhlOut,Quaternion *** pqtnhhOut);

void QMat_Min_Max(const Quaternion ** QMat,int * dim,double * dblRMin,double * dblRMax,
double * dblIMin,double * dblIMax,double * dblJMin,double * dblJMax,double * dblKMin,double * dblKMax);

void QMat_Min_Max2(const Quaternion ** QMatPure,int * dim,double * dblMin,double * dblMax);

int QMat_ChgScale_IntLvlMax(const Quaternion ** QMat, int * dim,int intLvlMax,
	double dblRMin,double dblRMax,double dblIMin,double dblIMax,double dblJMin,
    double dblJMax,double dblKMin,double dblKMax,Quaternion *** pQMatScaled);

int QImMat_ChgScale_IntLvlMax(const Quaternion ** QMat, int * dim,int intLvlMax,
	double dblIMin,double dblIMax,double dblJMin,
    double dblJMax,double dblKMin,double dblKMax,Quaternion *** pQMatScaled);

int QMat_ChgScale_IntLvlMax_PerChannel(const Quaternion ** QMat, int * dim,int intLvlMax,
	double dblRMin,double dblRMax,double dblIMin,double dblIMax,double dblJMin,
    double dblJMax,double dblKMin,double dblKMax,Quaternion *** pQMatScaled);
    
void QMatAbs(const Quaternion ** qtnqMat,int * dim,Quaternion *** pqtnMatAbs);

void QMatMax(const Quaternion ** qtnqMat,int * dim,double *** pdblMatMax);
void QMatMax2(const Quaternion ** qtnqMat,int * dim,double ** pdblMatMax);
    
//////////////////////////////////////////////////////////////////////////////////////////
//                OPERATION SUR DES VECTEURS QUATERNIONIQUES
//////////////////////////////////////////////////////////////////////////////////////////

//multiplication term by term
void QVectMultTBT(const Quaternion * QVect1,const Quaternion * QVect2,Quaternion ** pQMult,int dim);

int QVectAlloc(Quaternion **pQVect,int dim);

void QVectFree(Quaternion **pQVect);

void QInitVect(Quaternion **pQVect,int dim,Quaternion qtnValue);

void QVectScalMult(Quaternion **pVect,int dim,double dblScal);

void QVectConj(Quaternion **pVect,Quaternion *Vect,int dim);

void QVectDisp(Quaternion *Vect,int dim,FILE * stream);

#endif
