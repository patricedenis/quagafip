#include "complex.h"

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Fundamental complex definition functions

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

Complex Complex_InitXY(double x,double y)
{
  Complex c;
  c.r = x;
  c.i = y;
  return c;
}

double CReal(Complex cpx)
{
  return cpx.r;
}

double CIm(Complex cpx)
{
  return cpx.i;
}

//build a complex with his modulus Ro and phase Theta
Complex Complex_InitRoTheta(double dblRo,double dblTheta)
{
  Complex c;
  
  if(dblRo<0)
  {
    dblRo = -dblRo;//modulus is always positive;
    printf("modulus is always positive\n");
  }  
  c.r = dblRo * cos(dblTheta);
  c.i = dblRo * sin(dblTheta);
  return c;
}

Complex Complex_InitExpTheta(double dblTheta)
{
  return Complex_InitRoTheta(1,dblTheta);
}

Complex Complex_Sum(const Complex c1,const Complex c2)
{
  Complex Csum;
  
  Csum.r = c1.r + c2.r;
  Csum.i = c1.i + c2.i;
  return Csum;
}

Complex Complex_Diff(const Complex c1,const Complex c2)
{
  Complex Cdiff;
  
  Cdiff.r = c1.r - c2.r;
  Cdiff.i = c1.i - c2.i;
  return Cdiff;
}

Complex Complex_Mult(const Complex c1,const Complex c2)
{
  Complex Cmult;
  
  Cmult.r = c1.r*c2.r - c1.i*c2.i;
  Cmult.i = c1.r*c2.i + c1.i*c2.r;
  return Cmult;
}

Complex Complex_Scal_Mult(const Complex c,const double d)
{
  Complex Cmult;
  
  Cmult.r = c.r*d;
  Cmult.i = c.i*d;
  return Cmult;
}


Complex Complex_Div(const Complex c1,const Complex c2)
{
  return Complex_InitRoTheta( Complex_Modulus(c1)/Complex_Modulus(c2),
  		Complex_Phase(c1)-Complex_Phase(c2));
}

Complex Complex_Conjugate(const Complex c)
{
  Complex z;
  
  z.r = c.r;
  z.i = - c.i;
  return z;
}

double Complex_Modulus(const Complex c)
{
  return sqrt(c.r*c.r + c.i*c.i);
}

double Complex_Phase(const Complex c)
{
  if(c.r == 0)
    return M_PI/2;
  if(c.i == 0)
    return 0;
  if(c.r < 0)
    return atan(c.i/c.r);
  else
    return atan(c.i/c.r) + M_PI;
}

double Complex_Phase2(const Complex c)
{
  if (-DBL_EPSILON<c.r && c.r<DBL_EPSILON)
  	return M_PI;
  else
  {
    if (c.r>0 && c.i>0)
      return atan(c.i/c.r);
    else
    {
      if (c.r>0 && c.i<0)
        return atan(c.i/c.r) + 2*M_PI;
      else
        return atan(c.i/c.r) + M_PI;
    }
  }
}


/*
void ComplexTrigCoordinates(Complex z, double *r, double *phi)
{
    *r = sqrt(z.re*z.re+z.im*z.im);
    if (-DBL_EPSILON<z.re && z.re<DBL_EPSILON) {
      // On est dans un cas limite, par continuité on rend phi = PI/2 
      *phi = PI/2;
    }
    // Traitement des cas "normaux" 
    else if (z.re>0 && z.im>0) {
        *phi = atan(z.im/z.re);
				
    } else if (z.re>0 && z.im<0) {
        *phi = atan(z.im/z.re) + 2*PI;
    } else {
        *phi = atan(z.im/z.re) + PI;
    }
}


Complex ComplexDiv(Complex z1, Complex z2)
{
    Complex result;
    double z2ModSquare = (z2.re*z2.re + z2.im*z2.im);
    
    if (z2ModSquare < DBL_EPSILON) {
        result.re = result.im = DBL_MAX;
    }

    // Traitement des cas "normaux" 
    else {
        result.re = (z1.re*z2.re + z1.im*z2.im)/z2ModSquare;
        result.im = (-z1.re*z2.im + z1.im*z2.re)/z2ModSquare;
    }
    
    return result;
}

*/


void Complex_Display(const Complex c,FILE * stream)
{
  fprintf(stream,"z = %g + (%g) i\n",c.r,c.i);
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Matrix operations on complex

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

int Complex_Matrix_Allocate(int intHeight,int intWidth,Complex *** pcpxMatrix)
{
  int i,intCount,j;
	
	(*pcpxMatrix)=(struct stComplex **)malloc(sizeof(struct stComplex *)*intHeight);
  if ( (*pcpxMatrix) != NULL) /*allocation successful*/
  {
		for(i=0;i<intHeight;i++)
		{
			(*pcpxMatrix)[i]=(struct stComplex *)malloc(sizeof(struct stComplex)*intWidth);
			if ((*pcpxMatrix)[i] == NULL) //allocation error
    	{
				printf("allocation memory error\n");
				//we need to free the matrix memory that is already allocated
				for(j=i-1;j<=0;j--)
	  			free((*pcpxMatrix)[j]);
				free(*pcpxMatrix);
				return FALSE;
    	}
		}
	}
	else /*allocation error*/
  {
    printf("allocation memory error\n");
    return FALSE;
  }
  return TRUE;
}

/* we need to free the matrix memory */
void Complex_Matrix_Free(int intHeight,Complex *** cpxMatrix)
{
  int i;
  
  for(i=0;i<=intHeight-1;i++)
  {
      free((*cpxMatrix)[i]);
      (*cpxMatrix)[i]=NULL;
  }
  free(*cpxMatrix);
  (*cpxMatrix)=NULL;
  printf("complex memory free\n");
}

// this procedure take the 2nd and 3rd composant from the image MatYUV
// to fill a complex Matrix cpxMatrix.
int Get_Complex_Matrix_From_Image_23(double ** MatYUV,int intHeight,int intWidth,
	Complex *** pcpxMatrix)
{
	int i,j;
	
	if( Complex_Matrix_Allocate(intHeight,intWidth,pcpxMatrix) == TRUE)
	{
		for(i=0;i<intHeight;i++)
			for(j=0;j<intHeight;j++)
			{
				(*pcpxMatrix)[i][j] = Complex_InitXY(MatYUV[i][3*j+1],MatYUV[i][3*j+2]);
			}
	}
	else return FALSE;
	return TRUE;
}

void CMatrixCopy(const Complex ** cpxIn, Complex *** pcpxOut, int * dim)
{
  int i,j;
  for(i=0;i<dim[0];i++)
    for(j=0;j<dim[1];j++)
      (*pcpxOut)[i][j]=cpxIn[i][j];
}

// this fonction initialise the matrix given in "cpxMatrix"
//each reel and imaginary componant of the matrix will be 
// set with the "dblValue" value.
void Complex_Init_Matrix_Double(int intHeight,int intWidth,const double dblValue, Complex *** pcpxMatrix)
{
	int i,j;
	for(i=0;i<intHeight;i++)
		for(j=0;j<intHeight;j++)
		{
			((*pcpxMatrix)[i][j]).r = dblValue;
			((*pcpxMatrix)[i][j]).i = dblValue;
		}
}


// this fonction initialise the matrix given in "cpxMatrix"
// each componant of the matrix will be 
// set with the "CpxValue" value.
void Complex_Init_Matrix(int intHeight,int intWidth,const Complex cpxValue, Complex *** cpxMatrix)
{
  int i,j;
  for(i=0;i<intHeight;i++)
    for(j=0;j<intHeight;j++)
	  ((*cpxMatrix)[i][j]) = cpxValue;
}


//this procedure copy the real part of the complex matrix to fill the 2nd
//componant of the Image MatYUV, and the imaginary part for the 3rd componant
void Set_Complex_Matrix_To_Image_23(double *** MatYUV,int intHeight,int intWidth,
	const Complex ** cpxMatrix)
{
	int i,j;
	
	for(i=0;i<intHeight;i++)
		for(j=0;j<intHeight;j++)
		{
			(*MatYUV)[i][3*j+1] = (cpxMatrix[i][j]).r;
			(*MatYUV)[i][3*j+2] = (cpxMatrix[i][j]).i;
		}
}


void Complex_Copy_SubMatrix(const Complex ** MatFrom,Complex *** MatTo,int intXFrom,int intYFrom,
int intXTo,int intYTo,int intHeight,int intWidth)
{
	int i,j;
	for(i=0;i<intHeight;i++)
		for(j=0;j<intWidth;j++)
		{
			(*MatTo)[i+intXTo][j+intYTo].r = MatFrom[i+intXFrom][j+intYFrom].r;
			(*MatTo)[i+intXTo][j+intYTo].i = MatFrom[i+intXFrom][j+intYFrom].i;
		}
}
//  This function will shift the contents of the two sub-matrix pointed by
//  MatFrom and MatTo as describe in the following scheme
//		|-----------|				|-----------|
//		|  1  |  2  |				|  4  |  3  |
//		|-----------|  --->         |-----------|
//		|  3  |  4  |				|  2  |  1  |
//		|-----------|				|-----------|
void Complex_Matrix_Shift(const Complex ** MatFrom,Complex *** pMatShifted,int intHeight,int intWidth)
{
	int intMid;
	
	intMid = intHeight/2;
	Complex_Copy_SubMatrix(MatFrom,pMatShifted,0,0,intMid,intMid,intMid,intMid);	//1->4
	Complex_Copy_SubMatrix(MatFrom,pMatShifted,intMid,intMid,0,0,intMid,intMid);	//4->1
	Complex_Copy_SubMatrix(MatFrom,pMatShifted,intMid,0,0,intMid,intMid,intMid);	//2->3
	Complex_Copy_SubMatrix(MatFrom,pMatShifted,0,intMid,intMid,0,intMid,intMid);	//3->2
}


// this procedure take the composant from the array
// to fill the Complex Matrix CMatrix real part.
// caution the pCMatrix should be allocated before !
int CSetRMatrix(double tabRSet [],int intHeight,int intWidth,
	Complex *** pCMatrix)
{
  int i,j,ind=0;
  for(i=0;i<intHeight;i++)
    for(j=0;j<intHeight;j++)
	{
	  (*pCMatrix)[i][j].r = tabRSet[ind];
	  ind++;
    }
}

// this procedure take the composant from R and I arrays
// to fill the Complex Matrix CMatrix.
// caution the pCMatrix should be allocated before !
void CSetMatrix(double tabRSet [],double tabISet [],int intHeight,int intWidth,
	Complex *** pCMatrix)
{
  int i,j,ind=0;
  for(i=0;i<intHeight;i++)
	for(j=0;j<intHeight;j++)
	{
	  (*pCMatrix)[i][j].r = tabRSet[ind];
	  (*pCMatrix)[i][j].i = tabISet[ind];
	  ind++;
	}
}

//cette fonction permet de remplir la matrice avec les tableaux donnés pour sa partie reelle et sa partie imaginaire.
void CSetMatrixFromDblTabs(double * dblTabRealPart, double * dblTabImaginaryPart,int * dim,Complex *** pcpxMat)
{
  int i,j,ind=0;
  for(i=0;i<dim[0];i++)
    for(j=0;j<dim[1];j++)
    {
	    ((*pcpxMat)[i][j]) = Complex_InitXY(dblTabRealPart[ind],dblTabImaginaryPart[ind]);
	    ind++;
    }
}

//cette fonction permet de remplir la matrice avec les matrices données pour sa partie reelle et sa partie imaginaire.
void CSetMatrixFromDblMats(double ** dblMatRealPart, double ** dblMatImaginaryPart,int * dim,Complex *** pcpxMat)
{
  int i,j;
  for(i=0;i<dim[0];i++)
    for(j=0;j<dim[1];j++)
	    ((*pcpxMat)[i][j]) = Complex_InitXY(dblMatRealPart[i][j],dblMatImaginaryPart[i][j]);
}


// this procedure fill the r, and i arrays
// from the Complex Matrix QMatrix imaginary composants.
void CGetMatrixDoubleTab(double * pr[],double * pi[],int intHeight,int intWidth,const Complex ** CMatrix)
{
	int i,j,ind=0;
	 
	for(i=0;i<intHeight;i++)
		for(j=0;j<intWidth;j++)
		{
		  (*pr)[ind] = CReal(CMatrix[i][j]);
      (*pi)[ind] = CIm(CMatrix[i][j]);
		  ind++;
		}
}

// this procedure fill the r, and i arrays
// from the Complex Matrix QMatrix imaginary composants.
void CGetMatrixTab(Complex * pcpx[],int intHeight,int intWidth,const Complex ** CMatrix)
{
	int i,j,ind=0;
	 
	for(i=0;i<intHeight;i++)
		for(j=0;j<intWidth;j++)
		{
		  (*pcpx)[ind] = CMatrix[i][j];
		  ind++;
		}
}

void CMatrixMod(double * dblTabModulus [],int intHeight,int intWidth,const Complex ** CMat)
{
  int i,j,ind=0;
  
  for(i=0;i<intHeight;i++)
	for(j=0;j<intWidth;j++)
	{
	  (*dblTabModulus)[ind] = Complex_Modulus(CMat[i][j]);
	  ind++;
	}
}


void Complex_Matrix_Fourier_Transform( const Complex ** Matrix,
Complex *** pMatrix_Fourier_Transformed ,int intHeight,int intWidth)
{
	int p,q,P,Q;
	Complex cpxExp,cpxSum;
	double dblFactor;
	
	//dblFactor = -2*M_PI/((double)(sqrt(intHeight)*sqrt(intWidth)));
	dblFactor = -2*M_PI;
	
	// I'm looking for the frequency image, its coordinates will be describe
	// by the P and Q index.
	// U and V frequency components are the real and imaginary parts of the
	// Matrix_Fourier_Transformed Matrix respectively.
	// In the spacial image, the pixel will be pointed by the spacial 
	// coordinates p and q.
	// u and v are the real and imaginary component of the 
	// spatial Complex Matrix "Matrix".

	for (P=0;P<intHeight;P++)
		for (Q=0;Q<intWidth;Q++)
		{
			cpxSum = Complex_InitXY(0,0);
			for(p=0;p<intHeight;p++)
			{
				for(q=0;q<intWidth;q++)
				{
					cpxExp = Complex_InitExpTheta(dblFactor*(P*p/((double)intHeight)+Q*q/((double)intWidth)));
					cpxSum = Complex_Sum(cpxSum,Complex_Mult(Matrix[p][q],cpxExp));
				}
			}
			(*pMatrix_Fourier_Transformed)[P][Q].r = cpxSum.r;	//	R(P,Q)
			(*pMatrix_Fourier_Transformed)[P][Q].i = cpxSum.i;	//	I(P,Q)
		}
}

int IsCMatRealInt(int * intDim,const Complex ** cpxMat)
{
  int i,j;
  for(i=0;i<intDim[0];i++)
    for(j=0;j<intDim[1];j++)
      if(((int)(CIm(cpxMat[i][j])))!=0) return FALSE;
  return TRUE;
}

int IsCMatImInt(int * intDim,const Complex ** cpxMat)
{
  int i,j;
  for(i=0;i<intDim[0];i++)
    for(j=0;j<intDim[1];j++)
      if(((int)(CReal(cpxMat[i][j])))!=0) return FALSE;
  return TRUE;
}

int IsCMatRealEpsilon(int * intDim,const Complex ** cpxMat,double dblEpsilon)
{
  int i,j;
  for(i=0;i<intDim[0];i++)
    for(j=0;j<intDim[1];j++)
      if(fabs(CIm(cpxMat[i][j]))>=dblEpsilon) return FALSE;
  return TRUE;
}

void CMatDisp(const Complex **CMat,int * intDim,FILE * stream)
{
  int i,j;
  for(i=0;i<intDim[0];i++)
  {
    for(j=0;j<intDim[1];j++)
    {
      fprintf(stream,"%4.2lf+(%4.2lf)i  ",CMat[i][j].r,CMat[i][j].i);
    }
  fprintf(stream,"\n");
  }
}
//MatrixIn1 * MatrixIn2 -> MatrixOut multiplication Term by Term
void CMatMultTBT(const Complex ** CMatIn1, const Complex ** CMatIn2, Complex *** pCMatOut,int * dim)
{
  int i,j;
  
  Complex_Init_Matrix(dim[0],dim[1],Complex_InitXY(1.0,0.),pCMatOut);
  for(i=0;i<dim[0];i++)
    for(j=0;j<dim[1];j++)
      (*pCMatOut)[i][j]=Complex_Mult(CMatIn1[i][j],CMatIn2[i][j]);
}

void CMatConj(const Complex ** CMatIn, int * dim,	Complex *** pCMatConj)
{
  int i,j;
  
  Complex_Init_Matrix(dim[0],dim[1],Complex_InitXY(1.0,0.),pCMatConj);
  for(i=0;i<dim[0];i++)
    for(j=0;j<dim[1];j++)
      (*pCMatConj)[i][j]=Complex_Conjugate(CMatIn[i][j]);
}

void CMatScalMult(const Complex ** CMatIn, int * dim,	Complex *** pCMatOut,double dblScal)
{
  int i,j;
  
  Complex_Init_Matrix(dim[0],dim[1],Complex_InitXY(0.0,0.),pCMatOut);
  for(i=0;i<dim[0];i++)
    for(j=0;j<dim[1];j++)
      (*pCMatOut)[i][j]=Complex_Scal_Mult(CMatIn[i][j],dblScal);
}




void CSubImageUnion(const Complex ** cpxll,const Complex ** cpxlh,const Complex ** cpxhl,
  const Complex ** cpxhh, int *dim,Complex *** pcpxComp)
{
  int i,j;
  
  Complex_Init_Matrix(dim[0],dim[1],Complex_InitXY(0.,0.),pcpxComp);
  for(i=0;i<dim[0]/2;i++)
  {
    for(j=0;j<dim[1]/2;j++)
    {
      (*pcpxComp)[i][j]= cpxll[i][j];
      (*pcpxComp)[i+dim[0]/2][j]= cpxlh[i][j];
      (*pcpxComp)[i][j+dim[1]/2]= cpxhl[i][j];
      (*pcpxComp)[i+dim[0]/2][j+dim[1]/2]= cpxhh[i][j];
    }
  }
}


// ici on contruit la matrice qui contient toutes les informations : l'imagette basse frequence
// plus les details horizontaux, verticaux et diagonaux
void CThumbnailsConstruct(Complex *** pcpxMatOut,int * intDimOut, const Complex **cpxLow,
 const Complex **cpxHighV,const Complex **cpxHighH,const Complex **cpxHighD)
{
  int i,j;
  
  //ici on remplit la matrice des sortie avec uniCuement un point sur deux de la matrice d'origine.
  for(i=0;i<intDimOut[0]/2;i++)  
    for(j=0;j<intDimOut[1]/2;j++)
    {
      (*pcpxMatOut)[i][j]=cpxLow[i][j];
      (*pcpxMatOut)[i][j+intDimOut[1]/2]=cpxHighV[i][j];
      (*pcpxMatOut)[i+intDimOut[0]/2][j]=cpxHighH[i][j];
      (*pcpxMatOut)[i+intDimOut[0]/2][j+intDimOut[1]/2]=cpxHighD[i][j];
    }
}


void CAdd4(Complex *** pcpxCMat,int *dim,Complex ** cpxCLow,Complex ** cpxCHighV,Complex ** cpxCHighH,Complex ** cpxCHighD)
{
  int i,j;
  
  for(i=0;i<dim[0];i++)  
    for(j=0;j<dim[1];j++)
    {
      (*pcpxCMat)[i][j]=Complex_Sum(cpxCLow[i][j],cpxCHighV[i][j]);
      (*pcpxCMat)[i][j]=Complex_Sum((*pcpxCMat)[i][j],cpxCHighH[i][j]);
      (*pcpxCMat)[i][j]=Complex_Sum((*pcpxCMat)[i][j],cpxCHighD[i][j]);
    }
}

void CThumbnailsConstructWithoutLow(Complex *** pcpxMatOut,int * intDimOut,
 const Complex **cpxHighV,const Complex **cpxHighH,const Complex **cpxHighD)
{
  int i,j;
  
  //ici on remplit la matrice des sortie avec uniCuement un point sur deux de la matrice d'origine.
  for(i=0;i<intDimOut[0]/2;i++)  
    for(j=0;j<intDimOut[1]/2;j++)
    {
      (*pcpxMatOut)[i][j+intDimOut[1]/2]=cpxHighV[i][j];
      (*pcpxMatOut)[i+intDimOut[0]/2][j]=cpxHighH[i][j];
      (*pcpxMatOut)[i+intDimOut[0]/2][j+intDimOut[1]/2]=cpxHighD[i][j];
    }
}

void CThumbnailsConstructLow(Complex *** pcpxMatOut,int * intDimLow, const Complex **cpxLow)
{
  int i,j;
  
  //ici on remplit la matrice des sortie avec uniCuement un point sur deux de la matrice d'origine.
  for(i=0;i<intDimLow[0];i++)  
    for(j=0;j<intDimLow[1];j++)
      (*pcpxMatOut)[i][j]=cpxLow[i][j];
}

void CMatCopy(Complex ***pcpxMatCpy,int *intDim,const Complex **cpxMat)
{
  int i,j;
  for(i=0;i<intDim[0];i++)
    for(j=0;j<intDim[1];j++)
      (*pcpxMatCpy)[i][j]=cpxMat[i][j];
}

void CThumbnailsExtract(const Complex ** cpxMat,int * intDim, Complex ***pcpxLow,
 Complex ***pcpxHighV,Complex ***pcpxHighH,Complex ***pcpxHighD)
{
int i,j;
  
  //les sous matrices complexes sont deja allouées
  for(i=0;i<intDim[0]/2;i++)  
    for(j=0;j<intDim[1]/2;j++)
    {
      (*pcpxLow)[i][j]=cpxMat[i][j];
      (*pcpxHighV)[i][j]=cpxMat[i][j+intDim[1]/2];
      (*pcpxHighH)[i][j]=cpxMat[i+intDim[0]/2][j];
      (*pcpxHighD)[i][j]=cpxMat[i+intDim[0]/2][j+intDim[1]/2];
    }
}

void CNThumbnailsExtract(const Complex **cpxMat,int *intDim,int intN,Complex ***pcpxll,Complex ***pcpxlh,
Complex ***pcpxhl,Complex ***pcpxhh)
{
  int i,j;
  int dimTempB[2];
  
  dimTempB[0]=intDim[0]>>intN-1;
  dimTempB[1]=intDim[1]>>intN-1;
  
  //les sous matrices complexes sont deja allouées
  for(i=0;i<dimTempB[0]/2;i++)  
    for(j=0;j<dimTempB[1]/2;j++)
    {
      (*pcpxll)[i][j]=cpxMat[i][j];
      (*pcpxlh)[i][j]=cpxMat[i][j+dimTempB[1]/2];
      (*pcpxhl)[i][j]=cpxMat[i+dimTempB[0]/2][j];
      (*pcpxhh)[i][j]=cpxMat[i+dimTempB[0]/2][j+dimTempB[1]/2];
    }
}


void CThumbnailsExtractLow(const Complex ** cpxMat,int * intDim, Complex ***pcpxLow)
{
  int i,j;
  
  //les sous matrices complexes sont deja allouées
  for(i=0;i<intDim[0]/2;i++)  
    for(j=0;j<intDim[1]/2;j++)
      (*pcpxLow)[i][j]=cpxMat[i][j];
}


void CMatAdd4(Complex *** pcpxCMat,int *dim,const Complex ** cpxCLow,
  const Complex ** cpxCHighV,const Complex ** cpxCHighH,const Complex ** cpxCHighD)
{
  int i,j;
  
  for(i=0;i<dim[0];i++)  
    for(j=0;j<dim[1];j++)
    {
      (*pcpxCMat)[i][j]=Complex_Sum(cpxCLow[i][j],cpxCHighV[i][j]);
      (*pcpxCMat)[i][j]=Complex_Sum((*pcpxCMat)[i][j],cpxCHighH[i][j]);
      (*pcpxCMat)[i][j]=Complex_Sum((*pcpxCMat)[i][j],cpxCHighD[i][j]);
    }
}


//void CMatScaleIm0255(const Complex **cpxMatIn,int *dim, Complex *** pcpxMatOut)
//{
//  double dblMin,dblMax;
//  double dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax;
//  
//  
//  //on va chercher les extrema
//  CMat_Min_Max(cpxMatIn,dim,&dblRMin,&dblRMax,&dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
//  //ensuite on modifie l'echelle des parties imaginaires
//  CImMat_ChgScale_IntLvlMax(cpxMatIn,dim,255,
//	dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax,pcpxMatOut);
//}

void CMatScaleRe0255(const Complex** cpxlh,const Complex** cpxhl,const Complex** cpxhh,
  int *dim, Complex *** pcpxlhOut,Complex *** pcpxhlOut,Complex *** pcpxhhOut)
{
  double dblMin,dblMax;
  double dblRMin,dblRMax,dblIMin,dblIMax;
  
  //LH
  //on va chercher les extrema
  CMat_Min_Max((const Complex **)cpxlh,dim,&dblRMin,&dblRMax,&dblIMin,&dblIMax);
  //ensuite on modifie l'echelle des parties imaginaires
  CReMat_ChgScale_IntLvlMax((const Complex **)cpxlh,dim,255,dblRMin,dblRMax,pcpxlhOut);
  //HL
  CMat_Min_Max((const Complex **)cpxhl,dim,&dblRMin,&dblRMax,&dblIMin,&dblIMax);
  //ensuite on modifie l'echelle des parties imaginaires
  CReMat_ChgScale_IntLvlMax((const Complex **)cpxhl,dim,255,dblRMin,dblRMax,pcpxhlOut);
  //HH
  CMat_Min_Max((const Complex **)cpxhh,dim,&dblRMin,&dblRMax,&dblIMin,&dblIMax);
  //ensuite on modifie l'echelle des parties imaginaires
  CReMat_ChgScale_IntLvlMax((const Complex **)cpxhh,dim,255,dblRMin,dblRMax,pcpxhhOut);
}

void CMat_Min_Max(const Complex ** CMat,int * dim,double * dblRMin,double * dblRMax,
double * dblIMin,double * dblIMax)
{
int bx,by;
    
    *dblRMin=DBL_MAX;
    *dblRMax=DBL_MIN;
    *dblIMin=DBL_MAX;
    *dblIMax=DBL_MIN;
    
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        *dblRMin = Denis_Min(*dblRMin,CMat[bx][by].r);
        *dblRMax = Denis_Max(*dblRMax,CMat[bx][by].r);
        *dblIMin = Denis_Min(*dblIMin,CMat[bx][by].i);
        *dblIMax = Denis_Max(*dblIMax,CMat[bx][by].i);
      }
}

void CMat_Min_Max2(const Complex ** CMatPure,int * dim,double * dblMin,double * dblMax)
{
int bx,by;
    
    *dblMin=DBL_MAX;
    *dblMax=DBL_MIN;
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        *dblMin = Denis_Min(*dblMin,CMatPure[bx][by].r);
        *dblMax = Denis_Max(*dblMax,CMatPure[bx][by].r);
        *dblMin = Denis_Min(*dblMin,CMatPure[bx][by].i);
        *dblMax = Denis_Max(*dblMax,CMatPure[bx][by].i);
      }
}
//////////////////////////////////////////////////////////////////////////////

//              C H A N G E M E N T        D ' E C H E L L E 

//////////////////////////////////////////////////////////////////////////////



int CMat_ChgScale_IntLvlMax(const Complex ** CMat, int * dim,int intLvlMax,
	double dblRMin,double dblRMax,double dblIMin,double dblIMax,Complex *** pCMatScaled)
{
  int bx,by;
  double Min,Max;
  
  Min=Denis_Min(dblRMin,dblIMin);
    
  Max=Denis_Max(dblRMax,dblIMax);
    
  if(((dblRMax-dblRMin) != 0.)&&((dblIMax-dblIMin) != 0.))
  {
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        ((*pCMatScaled)[bx][by]).r = (double)((int)(((CMat[bx][by].r - Min)*intLvlMax/(Max - Min))+ 0.5));
        ((*pCMatScaled)[bx][by]).i = (double)((int)(((CMat[bx][by].i - Min)*intLvlMax/(Max - Min))+ 0.5));     
      }
      return TRUE;
  }
  else
  {
    printf("division by zero\n");
    return FALSE;
  }
}

int CReMat_ChgScale_IntLvlMax(const Complex ** CMat, int * dim,int intLvlMax,
	double dblRMin,double dblRMax,Complex *** pCMatScaled)
{
  int bx,by;
  double Min,Max;
      
  if(((dblRMax-dblRMin) != 0.))
  {
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        ((*pCMatScaled)[bx][by]).r = (double)((int)(((CMat[bx][by].r - dblRMin)*intLvlMax/(dblRMax - dblRMin))+ 0.5));
      }
      return TRUE;
  }
  else
  {
    printf("division by zero\n");
    return FALSE;
  }
}


int CMat_ChgScale_IntLvlMax_PerChannel(const Complex ** CMat, int * dim,int intLvlMax,
	double dblRMin,double dblRMax,double dblIMin,double dblIMax,Complex *** pCMatScaled)
{
  int bx,by;

  if(((dblRMax-dblRMin) != 0.)&&((dblIMax-dblIMin) != 0.))
  {
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        ((*pCMatScaled)[bx][by]).r = (double)((int)(((CMat[bx][by].r - dblRMin)*intLvlMax/(dblRMax - dblRMin))+ 0.5));
        ((*pCMatScaled)[bx][by]).i = (double)((int)(((CMat[bx][by].i - dblIMin)*intLvlMax/(dblIMax - dblIMin))+ 0.5));
      }
      return TRUE;
  }
  else
  {
    printf("division by zero\n");
    return FALSE;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////
//                OPERATION SUR DES VECTEURS COMPLEXES
//////////////////////////////////////////////////////////////////////////////////////////

//multiplication term by term
void CVectMultTBT(const Complex * CVect1,const Complex * CVect2,Complex ** pCMult,int dim)
{
  int i;
  for(i=0;i<dim;i++)
      (*pCMult)[i]=Complex_Mult(CVect1[i],CVect2[i]);
}

int CVectAlloc(Complex **pCVect,int dim)
{
  (*pCVect)=(Complex *)malloc(dim*sizeof(Complex));
  if(*pCVect)  return TRUE;
  else {printf("complex vector allocation error\n"); return FALSE;}
}

void CVectFree(Complex **pCVect)
{
  free(*pCVect);
  (*pCVect)=NULL;
}

void CInitVect(Complex **pCVect,int dim,Complex cpxValue)
{
  int i;
  for(i=0;i<dim;i++)
    (*pCVect)[i]=cpxValue;
}

void CVectScalMult(Complex **pVect,int dim,double dblScal)
{
  int i;
  for(i=0;i<dim;i++)
    (*pVect)[i]=Complex_Scal_Mult((*pVect)[i],dblScal);
}

void CVectConj(Complex **pVect,Complex *Vect,int dim)
{
  int i;
  for(i=0;i<dim;i++)
    (*pVect)[i]=Complex_Conjugate(Vect[i]);
}

void CVectDisp(const Complex *Vect,int dim,FILE * stream)
{
  int i,j;
  for(i=0;i<dim;i++)
    fprintf(stream,"%4.2lf+(%4.2lf)i  ",Vect[i].r,Vect[i].i);
  fprintf(stream,"\n");
}

//cette fonction permet de recuperer la nieme ligne de la matrice Mat
void CGetMatIRow(const Complex ** Mat,Complex **pVect,int * dim,int n)
{
  int i;
  for(i=0;i<dim[1];i++)
    (*pVect)[i]=Mat[n][i];
}

void CSetMatIRow(Complex *** pMat,const Complex *Vect,int * dim,int n)
{
  int i;
  for(i=0;i<dim[0];i++)
    (*pMat)[n][i]=Vect[i];
}

void CGetMatIColumn(const Complex ** Mat,Complex **pVect,int * dim,int n)
{
  int i;
  for(i=0;i<dim[0];i++)
    (*pVect)[i]=Mat[i][n];
}

void CSetMatIColumn(Complex *** pMat,const Complex *Vect,int * dim,int n)
{
  int i;
  for(i=0;i<dim[0];i++)
    (*pMat)[i][n]=Vect[i];
}


