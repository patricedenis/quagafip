#ifndef _MAGICKIO_H_
#define _MAGICKIO_H_

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <magick\api.h>

#include "complex.h"
#include "quaternion.h"
#include "matrix_double.h"

//methodes utilisant simplement les tableaux de doubles
//lecture
void MgkTypeImage(char *nom_file,int dim[],int type[]);
void MgkLireImgGris(char *nom_file,double intensite[]);
void MgkLireImgCouleur(char *nom_file,double r[],double g[],double b[]);
//ecriture
void MgkEcrireImgGris(char *nom_file,double intensite[],int dim[]);
void MgkEcrireImgCouleur(char *nom_file,double r[],double g[],double b[],int dim[]);

//matrice double
void MgkDblGetMatrixAndEcrireImgGris(int * dim,const double ** dblMat,
  char * strFileNameOut,FILE *stream);
void MgkLireImgGrisAndAllocAndDblSetMatrix(char * strInFileName,int intDim[],double *** pdblMat);
//complexes
void MgkLireImgGrisAndAllocAndCSetMatrix(char * strInFileName,int intDim[],Complex *** pcpxMat);
void MgkCGetMatrixAndEcrireImgGris(int * dim,const Complex ** cpxMat,char * strFileNameOut);

void MgkLireImgCouleurAndAllocAndCSetMatrixMarg(char * strInFileName,int intDim[],
  Complex *** pcpxMat1,Complex *** pcpxMat2,Complex *** pcpxMat3);
  
void MgkCGetMatrixAndEcrireImgCouleurMarg(int * dim,const Complex ** cpxMat1,const Complex ** cpxMat2,
  const Complex ** cpxMat3,char * strFileNameOut);  

//quaternions
void MgkLireImgCouleurAndAllocAndQSetMatrix(char * strInFileName,int intDim[],Quaternion *** pqtnMat);
void MgkQGetMatrixAndEcrireImgCouleur(int * dim,Quaternion ** qtnMat,char * strFileNameOut);

#endif

