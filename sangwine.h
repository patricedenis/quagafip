#ifndef _SANGWINE_H_
#define _SANGWINE_H_

#include "quaternion.h" 
#include "MagickIO.h" //pour l'ouverture des images

void qtnSympleticDecomposition(const Quaternion ** qtnMatrixInQtnBasis,
  double *** pdblMatrixProjR,  double *** pdblMatrixProjMu1,
  double *** pdblMatrixProjMu2, double *** pdblMatrixProjMu1Mu2,
  int * intDim, Quaternion qtnMu1, Quaternion qtnMu2);
void qtnChangementDeBaseSurLesImag(const Quaternion ** qtnMatrixInQtnBasis,
  double *** pdblMatProjE1,double *** pdblMatProjE2, double *** pdblMatProjE3,
  int * intDim, Quaternion qtnE1, Quaternion qtnE2, Quaternion qtnE3);
void qtnConstructParaMat(const double ** dblMatMu1, Quaternion *** pqtnMatPara, int * intDim, Quaternion qtnMu1);
void qtnConstructPerpMat(const double ** dblMatMu2, const double ** dblMatMu3, Quaternion *** pqtnMatPerp,
  int * intDim, Quaternion qtnMu1, Quaternion qtnMu2);
void qtnRecomposition(const double ** dblMatMu1,const double ** dblMatMu2,const double ** dblMatMu3,
  int * intDim, Quaternion qtnMu1, Quaternion qtnMu2,Quaternion *** pQMatrix);
void qtnConstructSympleticDecomposition(char * strFile)

#endif

