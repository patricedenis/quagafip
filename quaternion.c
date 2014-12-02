#include "quaternion.h"

//creates and returns a quaternion with its Cartesian form
Quaternion QInit(const double a,const double b,const double c,const double d)
{
  Quaternion q;
  q.a = a;
	q.b = b;
	q.c = c;
	q.d = d;
	return q;
}

Quaternion QInitZero()
{
    Quaternion q;	
    q.a = 0.;
	q.b = 0.;
	q.c = 0.;
	q.d = 0.;
	return q;
}

//returns the conjugate quaternion
Quaternion QConj(const Quaternion q1)
{
    Quaternion q;
    q.a = q1.a;
    q.b = - q1.b;
    q.c = - q1.c;
    q.d = - q1.d;
    return q;
}

//returns the opponent quaternion
Quaternion QOpp(const Quaternion q1)
{
    Quaternion q;
    q.a = - q1.a;
    q.b = - q1.b;
    q.c = - q1.c;
    q.d = - q1.d;
    return q;
}

//returns the real part
double QReal(const Quaternion q)
{
    return q.a;
}

//returns the first imaginary part
double QImI(const Quaternion q)
{
    return q.b;
}

//returns the second imaginary part
double QImJ(const Quaternion q)
{
    return q.c;
}

//returns the third imaginary part
double QImK(const Quaternion q)
{
    return q.d;
}

//Creates and returns a quaternion from its exponential form
//Nu must be a single unit quaternion
Quaternion QInitExp(const Quaternion Nu,const double theta)
{
    double sinus;
    sinus=sin(theta);
    return QInit(cos(theta),QImI(Nu)*sinus,QImJ(Nu)*sinus,QImK(Nu)*sinus);
}

//returns the complete imaginary part
Quaternion QImag(const Quaternion q)
{
    Quaternion q1;
    q1.a = 0;
    q1.b = q.b;
    q1.c = q.c;
    q1.d = q.d;
    return q1;
}

//returns the quaternion norm
double QNorm(const Quaternion q)
{
    return sqrt(q.a*q.a + q.b*q.b + q.c*q.c + q.d*q.d);
}

//return the quaternion inverse
Quaternion QInv(const Quaternion q)
{
    return QScalDiv(QConj(q),QNorm(q)*QNorm(q));
}

//returns the division between quaternion q1 and q2 
//q2 must be different than 0
Quaternion QDiv(const Quaternion q1, const Quaternion q2)
{
    return QMult(q1,QInv(q2));
}

//displays the quaternion q on the input stream
void QDisp(const Quaternion q,FILE * stream)
{
  fprintf(stream,"q = %g + (%g) i + (%g) j + (%g) k",q.a,q.b,q.c,q.d);
}

//returns the sum between q1 & q2
Quaternion QAdd(const Quaternion q1,const Quaternion q2)
{
	Quaternion q;
    q.a = q1.a + q2.a;
    q.b = q1.b + q2.b;
    q.c = q1.c + q2.c;
    q.d = q1.d + q2.d;
    return q;
}


//returns the sum between q1 and a real
Quaternion QScalAdd(const Quaternion q1,const double dbl1)
{
    Quaternion q;
    q.a = q1.a + dbl1;
    q.b = q1.b;
    q.c = q1.c;
    q.d = q1.d;
    return q;
}

//returns the difference between q1 & q2
Quaternion QDiff(const Quaternion q1,const Quaternion q2)
{
	Quaternion q;
    q.a = q1.a - q2.a;
    q.b = q1.b - q2.b;
    q.c = q1.c - q2.c;
    q.d = q1.d - q2.d;
    return q;
}

//returns the multiplication between q1 & q2
Quaternion QMult(const Quaternion q1,const Quaternion q2)
{
    Quaternion q;
    q.a = q1.a*q2.a - q1.b*q2.b - q1.c*q2.c - q1.d*q2.d;
    q.b = q1.a*q2.b + q1.b*q2.a + q1.c*q2.d - q1.d*q2.c;
    q.c = q1.a*q2.c - q1.b*q2.d + q1.c*q2.a + q1.d*q2.b;
    q.d = q1.a*q2.d + q1.b*q2.c - q1.c*q2.b + q1.d*q2.a;
    return q;
}

//returns the scalar product of q1 & q2
double QScalProd(const Quaternion q1,const Quaternion q2)
{
    return q1.a*q2.a + q1.b*q2.b + q1.c*q2.c + q1.d*q2.d;
}

//returns the multiplication of quaternion q1 with a scalar
Quaternion QScalMult(const Quaternion q1,const double dbl1)
{
    Quaternion q;
    q.a = q1.a * dbl1;
    q.b = q1.b * dbl1;
    q.c = q1.c * dbl1;
    q.d = q1.d * dbl1;
    return q;
}

//returns the division of quaternion q by a scalar
Quaternion QScalDiv(const Quaternion q1,const double dbl1)
{
    Quaternion q;
    q.a = q1.a / dbl1;
    q.b = q1.b / dbl1;
    q.c = q1.c / dbl1;
    q.d = q1.d / dbl1;
    return q;
}

void QExp_LogRo_Mu_Phi(const Quaternion q, double * LogRo, Quaternion * qMu, double * Phi)
{
    *LogRo = log(QNorm(q));
    *qMu = QScalDiv(QImag(q),QNorm(QImag(q)));//qMu est donc un quaternion pur
    *Phi = atan(QNorm(QImag(q))/QReal(q));
}

void QExp_Ro_Mu_Phi(const Quaternion q, double * Ro, Quaternion * qMu, double * Phi)
{
    *Ro = QNorm(q);
    *qMu = QScalDiv(QImag(q),QNorm(QImag(q)));
    *Phi = atan(QNorm(QImag(q))/QReal(q));
}

void QParOrtho(const Quaternion q,const Quaternion p, Quaternion * qPar, Quaternion * qOrtho)
{
    *qPar = QScalMult(QAdd(q,QMult(p,QMult(q,p))),0.5);
    *qOrtho = QScalMult(QAdd(q,QOpp(QMult(p,QMult(q,p)))),0.5);
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Matrix operations on quaternion

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

int QMatrixAllocate(int intHeight,int intWidth,Quaternion *** pQMatrix)
{
    int i,intCount,j;
	
    (*pQMatrix)=(struct stQuaternion **)malloc(sizeof(struct stQuaternion *)*intHeight);
    if ( (*pQMatrix) != NULL) /*allocation successful*/
    {
		for(i=0;i<intHeight;i++)
		{
			(*pQMatrix)[i]=(struct stQuaternion *)malloc(sizeof(struct stQuaternion)*intWidth);
			if ((*pQMatrix)[i] == NULL) //allocation error
  			{
				printf("allocation memory error\n");
				//we need to free the matrix memory that is already allocated
				for(j=i-1;j<=0;j--)
	  			free((*pQMatrix)[j]);
				free(*pQMatrix);
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
void QMatrixFree(int intHeight,Quaternion *** pQMatrix)
{
  int i;
  
  for(i=0;i<=intHeight-1;i++)
      free((*pQMatrix)[i]);
  free(*pQMatrix);
  printf("quaternion matrix memory free\n");
}

void QMatrixDisp(const Quaternion **qMat,int intHeight,int intWidth,int intPart,FILE * stream)
{
  int i,j;

  switch(intPart)
  {
    case 1 : printf("partie reelle\n");
             for(i=0;i<intHeight;i++)
             {
               for(j=0;j<intWidth;j++)
                 printf("%6.2lf ",qMat[i][j].a);
               printf("\n");
             }
             break;
    case 2 : printf("partie imaginaire I\n");
             for(i=0;i<intHeight;i++)
             {
               for(j=0;j<intWidth;j++)
                 printf("%6.2lf ",qMat[i][j].b);
               printf("\n");
             }
             break;
    case 3 : printf("partie imaginaire J\n");
             for(i=0;i<intHeight;i++)
             {
               for(j=0;j<intWidth;j++)
                 printf("%6.2lf ",qMat[i][j].c);
               printf("\n");
             }
             break;
    case 4 : printf("partie imaginaire K\n");
             for(i=0;i<intHeight;i++)
             {
               for(j=0;j<intWidth;j++)
                 printf("%6.2lf ",qMat[i][j].d);
               printf("\n");
             }
             break;
    default : printf("partie reelle\n");
             for(i=0;i<intHeight;i++)
             {
               for(j=0;j<intWidth;j++)
                 printf("%6.2lf ",qMat[i][j].a);
               printf("\n");
             }
             printf("partie imaginaire I\n");
             for(i=0;i<intHeight;i++)
             {
               for(j=0;j<intWidth;j++)
                 printf("%6.2lf ",qMat[i][j].b);
               printf("\n");
             }
             printf("partie imaginaire J\n");
             for(i=0;i<intHeight;i++)
             {
               for(j=0;j<intWidth;j++)
                 printf("%6.2lf ",qMat[i][j].c);
               printf("\n");
             }
             printf("partie imaginaire K\n");
             for(i=0;i<intHeight;i++)
             {
               for(j=0;j<intWidth;j++)
                 printf("%6.2lf ",qMat[i][j].d);
               printf("\n");
             }
  }
}

//cette procedure initialise la matrice avec des valeurs nulles
int QMatInitZero(Quaternion ***pqtnMat, int * dim)
{
  int i,j,ind=0;
  
  if(pqtnMat)
  {
    for(i=0;i<dim[0];i++)
	  for(j=0;j<dim[0];j++)
	  (*pqtnMat)[i][j] = QInitZero();
  }
  else return FALSE;
  return TRUE;
}

// this procedure take the composant from r, g, and b arrays
// to fill the Quaternion Matrix QMatrix.
int QSetMatrix(double r [],double g[],double b[],int intHeight,int intWidth,
	Quaternion *** pQMatrix)
{
	int i,j,ind=0;
	
	if(pQMatrix)
	{
		for(i=0;i<intHeight;i++)
			for(j=0;j<intWidth;j++)
			{
				(*pQMatrix)[i][j] = QInit(0.,r[ind],g[ind],b[ind]);
				ind++;
			}
	}
	else return FALSE;
	return TRUE;
}

// cette procedure remplit la matrice avec les matrices doubles passées en parametres
int QSetMatrixFromDblMat(const double ** dblCompR,const double ** dblCompI, const double ** dblCompJ,
  const double ** dblCompK,int * intDim, Quaternion *** pQMatrix)
{
	int i,j;
	
	if(pQMatrix)
	{
		for(i=0;i<intDim[0];i++)
			for(j=0;j<intDim[1];j++)
				(*pQMatrix)[i][j] = QInit(dblCompR[i][j],dblCompI[i][j],dblCompJ[i][j],dblCompK[i][j]);
	}
	else return FALSE;
	return TRUE;
}


// this procedure take the composant from r, i, j and k arrays
// to fill the Quaternion Matrix QMatrix.
int QSetMatrixComp(double r [],double i[],double j[],double k[],int intHeight,int intWidth,
	Quaternion *** pQMatrix)
{
	int bx,by,ind=0;
	
	if( QMatrixAllocate(intHeight,intWidth,pQMatrix) == TRUE)
	{
		for(bx=0;bx<intHeight;bx++)
			for(by=0;by<intWidth;by++)
			{
				(*pQMatrix)[bx][by] = QInit(r[ind],i[ind],j[ind],k[ind]);
				ind++;
			}
	}
	else return FALSE;
	return TRUE;
}

// this procedure fill the r, g, and b arrays
// from the Quaternion Matrix QMatrix imaginary composants.
void QGetMatrixCompl(double * mr[],double * mi[],double * mj[],double * mk[],int intHeight,int intWidth,
	Quaternion ** QMatrix)
{
	int i,j,ind=0;
	 
	for(i=0;i<intHeight;i++)
		for(j=0;j<intWidth;j++)
		{
		  (*mr)[ind] = QReal(QMatrix[i][j]);
          (*mi)[ind] = QImI(QMatrix[i][j]);
		  (*mj)[ind] = QImJ(QMatrix[i][j]);
		  (*mk)[ind] = QImK(QMatrix[i][j]);
		  ind++;
		}
}

// this procedure fill the i, j, and k arrays
// from the Quaternion Matrix QMatrix imaginary composants.
void QGetMatrixImagPart(double * i[],double * j[],double * k[],int intHeight,int intWidth,
	Quaternion ** QMatrix)
{
	int x,y,ind=0;
	 
	for(x=0;x<intHeight;x++)
		for(y=0;y<intWidth;y++)
		{
		  (*i)[ind] = QImI(QMatrix[x][y]);
		  (*j)[ind] = QImJ(QMatrix[x][y]);
		  (*k)[ind] = QImK(QMatrix[x][y]);
		  ind++;
		}
}


void QCopyMatrixPartR(const Quaternion ** QMatrixFrom , Quaternion *** pQMatrixTo, int * intDim)
{
	int x,y,ind=0;
	 
	for(x=0;x<intDim[0];x++)
		for(y=0;y<intDim[1];y++)
		{
		  (*pQMatrixTo)[x][y].a = QMatrixFrom[x][y].a;
		}
}

void QCopyMatrixPartI(const Quaternion ** QMatrixFrom , Quaternion *** pQMatrixTo, int * intDim)
{
	int x,y,ind=0;
	 
	for(x=0;x<intDim[0];x++)
		for(y=0;y<intDim[1];y++)
		{
		  (*pQMatrixTo)[x][y].b = QMatrixFrom[x][y].b;
		}
}

void QCopyMatrixPartJ(const Quaternion ** QMatrixFrom , Quaternion *** pQMatrixTo, int * intDim)
{
	int x,y,ind=0;
	 
	for(x=0;x<intDim[0];x++)
		for(y=0;y<intDim[1];y++)
		{
		  (*pQMatrixTo)[x][y].c = QMatrixFrom[x][y].c;
		}
}

void QCopyMatrixPartK(const Quaternion ** QMatrixFrom , Quaternion *** pQMatrixTo, int * intDim)
{
	int x,y,ind=0;
	 
	for(x=0;x<intDim[0];x++)
		for(y=0;y<intDim[1];y++)
		{
		  (*pQMatrixTo)[x][y].d = QMatrixFrom[x][y].d;
		}
}




void QGetMatrixComponent(double * comp[],int intHeight,int intWidth,Quaternion ** QMatrix,char * strComp)
{
  int x,y,ind=0;

  for(x=0;x<intHeight;x++)
	for(y=0;y<intWidth;y++)
	{
	  if(strcmp(strComp,"R")==0)
        (*comp)[ind] = QReal(QMatrix[x][y]);
      if(strcmp(strComp,"I")==0)
        (*comp)[ind] = QImI(QMatrix[x][y]);
      if(strcmp(strComp,"J")==0)
        (*comp)[ind] = QImJ(QMatrix[x][y]);
      if(strcmp(strComp,"K")==0)
        (*comp)[ind] = QImK(QMatrix[x][y]);
	  ind++;
	}
}

void QSetMatrixComponent(double  comp[],int intHeight,int intWidth,Quaternion *** pQMatrix,char * strComp)
{
  int i,j,ind=0;
	
	for(i=0;i<intHeight;i++)
  	  for(j=0;j<intWidth;j++)
	  {
	    if(strcmp(strComp,"R")==0)
          (*pQMatrix)[i][j].a = comp[ind];
        if(strcmp(strComp,"I")==0)
          (*pQMatrix)[i][j].b = comp[ind];
        if(strcmp(strComp,"J")==0)
          (*pQMatrix)[i][j].c = comp[ind];
        if(strcmp(strComp,"K")==0)
          (*pQMatrix)[i][j].d = comp[ind];
		ind++;
	  }
}

// this fonction initialise the matrix given in "pQMatrix"
//each componant of the matrix will be 
// set with the "QInitValue" quaternion.
void QInitMatrix(int intHeight,int intWidth, Quaternion QInitValue, Quaternion *** pQMatrix)
{
	int i,j;
	for(i=0;i<intHeight;i++)
		for(j=0;j<intWidth;j++)
		{
			(*pQMatrix)[i][j] = QInitValue;
		}
}

//  This function will shift the contents of the two sub-matrix pointed by
//  MatFrom and MatTo as describe in the following scheme
//		|-----------|				|-----------|
//		|  1  |  2  |				|  4  |  3  |
//		|-----------|     --->      |-----------|
//		|  3  |  4  |				|  2  |  1  |
//		|-----------|				|-----------|
void QCopySubMatrix(Quaternion ** QMatFrom,Quaternion *** pQMatTo,int intXFrom,int intYFrom,
int intXTo,int intYTo,int intHeight,int intWidth)
{
	int i,j;
	for(i=0;i<intHeight;i++)
		for(j=0;j<intWidth;j++)
			(*pQMatTo)[i+intXTo][j+intYTo] = QMatFrom[i+intXFrom][j+intYFrom];
}

void QMatrixShift(Quaternion ** QMatFrom,Quaternion *** pQMatShifted,int intHeight,int intWidth)
{
	int intMid;
	
	intMid = intHeight/2;
	QCopySubMatrix(QMatFrom,pQMatShifted,0,0,intMid,intMid,intMid,intMid);	//1->4
	QCopySubMatrix(QMatFrom,pQMatShifted,intMid,intMid,0,0,intMid,intMid);	//4->1
	QCopySubMatrix(QMatFrom,pQMatShifted,intMid,0,0,intMid,intMid,intMid);	//2->3
	QCopySubMatrix(QMatFrom,pQMatShifted,0,intMid,intMid,0,intMid,intMid);	//3->2

}

int IsQMatImEpsilon(int * intDim,Quaternion ** qtnMat,double dblEpsilon)
{
  int i,j;
  for(i=0;i<intDim[0];i++)
    for(j=0;j<intDim[1];j++)
      if(fabs(QReal(qtnMat[i][j]))>=dblEpsilon) return FALSE;
  return TRUE;
}

int IsQMatReEpsilon(int * intDim,Quaternion ** qtnMat,double dblEpsilon)
{
  int i,j;
  for(i=0;i<intDim[0];i++)
    for(j=0;j<intDim[1];j++)
      if((fabs(QImI(qtnMat[i][j]))>=dblEpsilon)||(fabs(QImJ(qtnMat[i][j]))>=dblEpsilon)
         ||(fabs(QImK(qtnMat[i][j]))>=dblEpsilon)) return FALSE;
  return TRUE;
}


void QMatDisp(const Quaternion **QMat,int * intDim,FILE * stream)
{
  int i,j;
  for(i=0;i<intDim[0];i++)
  {
    for(j=0;j<intDim[1];j++)
    {
      fprintf(stream,"%4.2lf+(%4.2lf)i+(%4.2lf)j+(%4.2lf)k  ",QReal(QMat[i][j]),QImI(QMat[i][j]),QImJ(QMat[i][j]),QImK(QMat[i][j]));
    }
  fprintf(stream,"\n");
  }
}
//MatrixIn1 * MatrixIn2 -> MatrixOut multiplication Term by Term
void QMatMultTBT(Quaternion ** QMatIn1, Quaternion ** QMatIn2, Quaternion *** pQMatOut,int * dim)
{
  int i,j;
  
  QInitMatrix(dim[0],dim[1],QInit(1.0,0.,0.,0.),pQMatOut);
  for(i=0;i<dim[0];i++)
    for(j=0;j<dim[1];j++)
      (*pQMatOut)[i][j]=QMult(QMatIn1[i][j],QMatIn2[i][j]);
}

//cette fonction effectue la multiplication terme par terme avec une matrice de reels
void QMatMultTBTWithDblMat(Quaternion ** QMatIn1, double ** dblMatIn2, Quaternion *** pQMatOut,int * dim)
{
  int i,j;
  
  QInitMatrix(dim[0],dim[1],QInit(1.0,0.,0.,0.),pQMatOut);
  for(i=0;i<dim[0];i++)
    for(j=0;j<dim[1];j++)
      (*pQMatOut)[i][j]=QScalMult(QMatIn1[i][j],dblMatIn2[i][j]);
}


void QMatConj(Quaternion ** QMatIn, int * dim,	Quaternion *** pQMatConj)
{
  int i,j;
  
  QInitMatrix(dim[0],dim[1],QInit(1.0,0.,0.,0.),pQMatConj);
  for(i=0;i<dim[0];i++)
    for(j=0;j<dim[1];j++)
      (*pQMatConj)[i][j]=QConj(QMatIn[i][j]);
}


void QMatScalMult(Quaternion ** QMatIn, int * dim,	Quaternion *** pQMatOut,double dblScal)
{
  int i,j;
  
  QInitMatrix(dim[0],dim[1],QInit(0.0,0.,0.,0.),pQMatOut);
  for(i=0;i<dim[0];i++)
    for(j=0;j<dim[1];j++)
      (*pQMatOut)[i][j]=QScalMult(QMatIn[i][j],dblScal);
}


void QSubImageUnion(Quaternion ** qtnll,Quaternion ** qtnlh,Quaternion ** qtnhl,Quaternion ** qtnhh,
  int *dim,Quaternion *** pqtnComp)
{
  int i,j;
  
  QInitMatrix(dim[0],dim[1],QInit(0.0,0.,0.,0.),pqtnComp);
  for(i=0;i<dim[0]/2;i++)
  {
    for(j=0;j<dim[1]/2;j++)
    {
      (*pqtnComp)[i][j]= qtnll[i][j];
      (*pqtnComp)[i+dim[0]/2][j]= qtnlh[i][j];
      (*pqtnComp)[i][j+dim[1]/2]= qtnhl[i][j];
      (*pqtnComp)[i+dim[0]/2][j+dim[1]/2]= qtnhh[i][j];
    }
  }
}

void QThumbnailsConstruct(Quaternion *** pqtnMatOut,int * intDimOut, Quaternion **qtnLow,
 Quaternion **qtnHighV,Quaternion **qtnHighH,Quaternion **qtnHighD)
{
  int i,j;
  
  //ici on remplit la matrice des sortie avec uniquement un point sur deux de la matrice d'origine.
  for(i=0;i<intDimOut[0]/2;i++)  
    for(j=0;j<intDimOut[1]/2;j++)
    {
      (*pqtnMatOut)[i][j]=qtnLow[i][j];
      (*pqtnMatOut)[i][j+intDimOut[1]/2]=qtnHighV[i][j];
      (*pqtnMatOut)[i+intDimOut[0]/2][j]=qtnHighH[i][j];
      (*pqtnMatOut)[i+intDimOut[0]/2][j+intDimOut[1]/2]=qtnHighD[i][j];
    }
}

void QThumbnailsConstructWithoutLow(Quaternion *** pqtnMatOut,int * intDimOut,
 Quaternion **qtnHighV,Quaternion **qtnHighH,Quaternion **qtnHighD)
{
  int i,j;
  
  //ici on remplit la matrice des sortie avec uniquement un point sur deux de la matrice d'origine.
  for(i=0;i<intDimOut[0]/2;i++)  
    for(j=0;j<intDimOut[1]/2;j++)
    {
      (*pqtnMatOut)[i][j+intDimOut[1]/2]=qtnHighV[i][j];
      (*pqtnMatOut)[i+intDimOut[0]/2][j]=qtnHighH[i][j];
      (*pqtnMatOut)[i+intDimOut[0]/2][j+intDimOut[1]/2]=qtnHighD[i][j];
    }
}

void QThumbnailsConstructLow(Quaternion *** pqtnMatOut,int * intDimLow, Quaternion **qtnLow)
{
  int i,j;
  
  //ici on remplit la matrice des sortie avec uniquement un point sur deux de la matrice d'origine.
  for(i=0;i<intDimLow[0];i++)  
    for(j=0;j<intDimLow[1];j++)
      (*pqtnMatOut)[i][j]=qtnLow[i][j];
}

void QMatCopy(Quaternion ***pqtnMatCpy,int *intDim,const Quaternion **qtnMat)
{
  int i,j;
  for(i=0;i<intDim[0];i++)
    for(j=0;j<intDim[1];j++)
      (*pqtnMatCpy)[i][j]=qtnMat[i][j];
}

void QThumbnailsExtract(const Quaternion ** qtnMat,int * intDim, Quaternion ***pqtnLow,
 Quaternion ***pqtnHighV,Quaternion ***pqtnHighH,Quaternion ***pqtnHighD)
{
int i,j;
  
  //les sous matrices complexes sont deja allouées
  for(i=0;i<intDim[0]/2;i++)  
    for(j=0;j<intDim[1]/2;j++)
    {
      (*pqtnLow)[i][j]=qtnMat[i][j];
      (*pqtnHighV)[i][j]=qtnMat[i][j+intDim[1]/2];
      (*pqtnHighH)[i][j]=qtnMat[i+intDim[0]/2][j];
      (*pqtnHighD)[i][j]=qtnMat[i+intDim[0]/2][j+intDim[1]/2];
    }
}

void QNThumbnailsExtract(const Quaternion **qtnMat,int *intDim,int intN,Quaternion ***pqtnll,Quaternion ***pqtnlh,
Quaternion ***pqtnhl,Quaternion ***pqtnhh)
{
  int i,j;
  int dimTempB[2];
  
  dimTempB[0]=intDim[0]>>intN-1;
  dimTempB[1]=intDim[1]>>intN-1;
  
  //les sous matrices complexes sont deja allouées
  for(i=0;i<dimTempB[0]/2;i++)  
    for(j=0;j<dimTempB[1]/2;j++)
    {
      (*pqtnll)[i][j]=qtnMat[i][j];
      (*pqtnlh)[i][j]=qtnMat[i][j+dimTempB[1]/2];
      (*pqtnhl)[i][j]=qtnMat[i+dimTempB[0]/2][j];
      (*pqtnhh)[i][j]=qtnMat[i+dimTempB[0]/2][j+dimTempB[1]/2];
    }
}


void QThumbnailsExtractLow(Quaternion ** qtnMat,int * intDim, Quaternion ***pqtnLow)
{
  int i,j;
  
  //les sous matrices complexes sont deja allouées
  for(i=0;i<intDim[0]/2;i++)  
    for(j=0;j<intDim[1]/2;j++)
      (*pqtnLow)[i][j]=qtnMat[i][j];
}

//cette fonction somme les deux matrices passées en parametre
void QMatAdd(Quaternion *** pqtnMatAdd,Quaternion ** qtnMat1,Quaternion ** qtnMat2,int *dim)
{
  int i,j;
 
  for(i=0;i<dim[0];i++)  
    for(j=0;j<dim[1];j++)
      (*pqtnMatAdd)[i][j]=QAdd(qtnMat1[i][j],qtnMat1[i][j]);
}


void QAdd4(Quaternion *** pqtnQMat,int *dim,Quaternion ** qtnQLow,
  Quaternion ** qtnQHighV,Quaternion ** qtnQHighH,Quaternion ** qtnQHighD)
{
  int i,j;
  
  for(i=0;i<dim[0];i++)  
    for(j=0;j<dim[1];j++)
    {
      (*pqtnQMat)[i][j]=QAdd(qtnQLow[i][j],qtnQHighV[i][j]);
      (*pqtnQMat)[i][j]=QAdd((*pqtnQMat)[i][j],qtnQHighH[i][j]);
      (*pqtnQMat)[i][j]=QAdd((*pqtnQMat)[i][j],qtnQHighD[i][j]);
    }
}


//void QMatScaleIm0255(const Quaternion **qtnMatIn,int *dim, Quaternion *** pqtnMatOut)
//{
//  double dblMin,dblMax;
//  double dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax;
//  
//  
//  //on va chercher les extrema
//  QMat_Min_Max(qtnMatIn,dim,&dblRMin,&dblRMax,&dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
//  //ensuite on modifie l'echelle des parties imaginaires
//  QImMat_ChgScale_IntLvlMax(qtnMatIn,dim,255,
//	dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax,pqtnMatOut);
//}

void QMatScaleIm0255(Quaternion** qtnlh,Quaternion** qtnhl,Quaternion** qtnhh,
  int *dim, Quaternion *** pqtnlhOut,Quaternion *** pqtnhlOut,Quaternion *** pqtnhhOut)
{
  double dblMin,dblMax;
  double dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax;
  
  //LH
  //on va chercher les extrema
  QMat_Min_Max((const Quaternion **)qtnlh,dim,&dblRMin,&dblRMax,&dblIMin,&dblIMax,
    &dblJMin,&dblJMax,&dblKMin,&dblKMax);
  //ensuite on modifie l'echelle des parties imaginaires
  QImMat_ChgScale_IntLvlMax((const Quaternion **)qtnlh,dim,255,dblIMin,dblIMax,
    dblJMin,dblJMax,dblKMin,dblKMax,pqtnlhOut);
  //HL
  QMat_Min_Max((const Quaternion **)qtnhl,dim,&dblRMin,&dblRMax,&dblIMin,&dblIMax,
    &dblJMin,&dblJMax,&dblKMin,&dblKMax);
  //ensuite on modifie l'echelle des parties imaginaires
  QImMat_ChgScale_IntLvlMax((const Quaternion **)qtnhl,dim,255,dblIMin,dblIMax,
    dblJMin,dblJMax,dblKMin,dblKMax,pqtnhlOut);
  //HH
  QMat_Min_Max((const Quaternion **)qtnhh,dim,&dblRMin,&dblRMax,&dblIMin,&dblIMax,
    &dblJMin,&dblJMax,&dblKMin,&dblKMax);
  //ensuite on modifie l'echelle des parties imaginaires
  QImMat_ChgScale_IntLvlMax((const Quaternion **)qtnhh,dim,255,dblIMin,dblIMax,
    dblJMin,dblJMax,dblKMin,dblKMax,pqtnhhOut);
}

void QMat_Min_Max(const Quaternion ** QMat,int * dim,double * dblRMin,double * dblRMax,
double * dblIMin,double * dblIMax,double * dblJMin,double * dblJMax,double * dblKMin,double * dblKMax)
{
int bx,by;
    
    *dblRMin=DBL_MAX;
    *dblRMax=DBL_MIN;
    *dblIMin=DBL_MAX;
    *dblIMax=DBL_MIN;
    *dblJMin=DBL_MAX;
    *dblJMax=DBL_MIN;
    *dblKMin=DBL_MAX;
    *dblKMax=DBL_MIN;
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        *dblRMin = Denis_Min(*dblRMin,QMat[bx][by].a);
        *dblRMax = Denis_Max(*dblRMax,QMat[bx][by].a);
        *dblIMin = Denis_Min(*dblIMin,QMat[bx][by].b);
        *dblIMax = Denis_Max(*dblIMax,QMat[bx][by].b);
        *dblJMin = Denis_Min(*dblJMin,QMat[bx][by].c);
        *dblJMax = Denis_Max(*dblJMax,QMat[bx][by].c);
        *dblKMin = Denis_Min(*dblKMin,QMat[bx][by].d);
        *dblKMax = Denis_Max(*dblKMax,QMat[bx][by].d);
      }
}

void QMat_Min_Max2(const Quaternion ** QMatPure,int * dim,double * dblMin,double * dblMax)
{
int bx,by;
    
    *dblMin=DBL_MAX;
    *dblMax=DBL_MIN;
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        *dblMin = Denis_Min(*dblMin,QMatPure[bx][by].a);
        *dblMax = Denis_Max(*dblMax,QMatPure[bx][by].a);
        *dblMin = Denis_Min(*dblMin,QMatPure[bx][by].b);
        *dblMax = Denis_Max(*dblMax,QMatPure[bx][by].b);
        *dblMin = Denis_Min(*dblMin,QMatPure[bx][by].c);
        *dblMax = Denis_Max(*dblMax,QMatPure[bx][by].c);
        *dblMin = Denis_Min(*dblMin,QMatPure[bx][by].d);
        *dblMax = Denis_Max(*dblMax,QMatPure[bx][by].d);
      }
}
//////////////////////////////////////////////////////////////////////////////

//              C H A N G E M E N T        D ' E C H E L L E 

//////////////////////////////////////////////////////////////////////////////



int QMat_ChgScale_IntLvlMax(const Quaternion ** QMat, int * dim,int intLvlMax,
	double dblRMin,double dblRMax,double dblIMin,double dblIMax,double dblJMin,
    double dblJMax,double dblKMin,double dblKMax,Quaternion *** pQMatScaled)
{
  int bx,by;
  double Min,Max;
  
  Min=Denis_Min(dblRMin,dblIMin);
  Min=Denis_Min(Min,dblJMin);
  Min=Denis_Min(Min,dblKMin);
  
  Max=Denis_Max(dblRMax,dblIMax);
  Max=Denis_Max(Max,dblJMax);
  Max=Denis_Max(Max,dblKMax);
  
    
  
  if(((dblRMax-dblRMin) != 0.)&&((dblIMax-dblIMin) != 0.)&&((dblJMax-dblJMin) != 0.)&&((dblKMax-dblKMin) != 0.))
  {
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        ((*pQMatScaled)[bx][by]).a = (double)((int)(((QMat[bx][by].a - Min)*intLvlMax/(Max - Min))+ 0.5));
        ((*pQMatScaled)[bx][by]).b = (double)((int)(((QMat[bx][by].b - Min)*intLvlMax/(Max - Min))+ 0.5));
        ((*pQMatScaled)[bx][by]).c = (double)((int)(((QMat[bx][by].c - Min)*intLvlMax/(Max - Min))+ 0.5));
        ((*pQMatScaled)[bx][by]).d = (double)((int)(((QMat[bx][by].d - Min)*intLvlMax/(Max - Min))+ 0.5));       
      }
      return TRUE;
  }
  /*if(((dblIMax-dblIMin) != 0.)&&((dblJMax-dblJMin) != 0.)&&((dblKMax-dblKMin) != 0.))
  {
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        ((*QMatScaled)[bx][by]).a = (double)((int)(((QMat[bx][by].a - dblRMin)*intLvlMax/(dblRMax - dblRMin))+ 0.5));
        ((*QMatScaled)[bx][by]).b = (double)((int)(((QMat[bx][by].b - dblIMin)*intLvlMax/(dblIMax - dblIMin))+ 0.5));
        ((*QMatScaled)[bx][by]).c = (double)((int)(((QMat[bx][by].c - dblJMin)*intLvlMax/(dblJMax - dblJMin))+ 0.5));
        ((*QMatScaled)[bx][by]).d = (double)((int)(((QMat[bx][by].d - dblKMin)*intLvlMax/(dblKMax - dblKMin))+ 0.5));       
      }
      return TRUE;
  }*/
  else
  {
    printf("division by zero\n");
    return FALSE;
  }
}

int QImMat_ChgScale_IntLvlMax(const Quaternion ** QMat, int * dim,int intLvlMax,
	double dblIMin,double dblIMax,double dblJMin,
    double dblJMax,double dblKMin,double dblKMax,Quaternion *** pQMatScaled)
{
  int bx,by;
  double Min,Max;
  
  Min=Denis_Min(dblIMin,dblJMin);
  Min=Denis_Min(Min,dblKMin);
  
  Max=Denis_Max(dblIMax,dblJMax);
  Max=Denis_Max(Max,dblKMax);
  
    
  
  if(((dblIMax-dblIMin) != 0.)&&((dblJMax-dblJMin) != 0.)&&((dblKMax-dblKMin) != 0.))
  {
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        ((*pQMatScaled)[bx][by]).b = (double)((int)(((QMat[bx][by].b - Min)*intLvlMax/(Max - Min))+ 0.5));
        ((*pQMatScaled)[bx][by]).c = (double)((int)(((QMat[bx][by].c - Min)*intLvlMax/(Max - Min))+ 0.5));
        ((*pQMatScaled)[bx][by]).d = (double)((int)(((QMat[bx][by].d - Min)*intLvlMax/(Max - Min))+ 0.5));       
      }
      return TRUE;
  }
  else
  {
    printf("division by zero\n");
    return FALSE;
  }
}


int QMat_ChgScale_IntLvlMax_PerChannel(const Quaternion ** QMat, int * dim,int intLvlMax,
	double dblRMin,double dblRMax,double dblIMin,double dblIMax,double dblJMin,
    double dblJMax,double dblKMin,double dblKMax,Quaternion *** pQMatScaled)
{
  int bx,by;

  if(((dblRMax-dblRMin) != 0.)&&((dblIMax-dblIMin) != 0.)&&((dblJMax-dblJMin) != 0.)&&((dblKMax-dblKMin) != 0.))
  {
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        ((*pQMatScaled)[bx][by]).a = (double)((int)(((QMat[bx][by].a - dblRMin)*intLvlMax/(dblRMax - dblRMin))+ 0.5));
        ((*pQMatScaled)[bx][by]).b = (double)((int)(((QMat[bx][by].b - dblIMin)*intLvlMax/(dblIMax - dblIMin))+ 0.5));
        ((*pQMatScaled)[bx][by]).c = (double)((int)(((QMat[bx][by].c - dblJMin)*intLvlMax/(dblJMax - dblJMin))+ 0.5));
        ((*pQMatScaled)[bx][by]).d = (double)((int)(((QMat[bx][by].d - dblKMin)*intLvlMax/(dblKMax - dblKMin))+ 0.5));       
      }
      return TRUE;
  }
  else
  {
    printf("division by zero\n");
    return FALSE;
  }
}


void QMatAbs(const Quaternion ** qtnqMat,int * dim,Quaternion *** pqtnMatAbs)
{
  int bx,by;

  for (bx=0;bx<dim[0];bx++)
    for (by=0;by<dim[1];by++)
    {
        ((*pqtnMatAbs)[bx][by]).a = fabs(qtnqMat[bx][by].a);
        ((*pqtnMatAbs)[bx][by]).b = fabs(qtnqMat[bx][by].b);
        ((*pqtnMatAbs)[bx][by]).c = fabs(qtnqMat[bx][by].c);
        ((*pqtnMatAbs)[bx][by]).d = fabs(qtnqMat[bx][by].d);       
    }
}


void QMatMax(const Quaternion ** qtnqMat,int * dim,double *** pdblMatMax)
{
  int bx,by;
  double temp;

  for (bx=0;bx<dim[0];bx++)
    for (by=0;by<dim[1];by++)
    {
        temp = Denis_Max(qtnqMat[bx][by].a,qtnqMat[bx][by].b);
        temp = Denis_Max(temp,qtnqMat[bx][by].c);
        (*pdblMatMax)[bx][by] = Denis_Max(temp,qtnqMat[bx][by].d);       
    }
}



void QMatMax2(const Quaternion ** qtnqMat,int * dim,double ** pdblMatMax)
{
  int bx,by,ind=0;
  double temp;

  for (bx=0;bx<dim[0];bx++)
    for (by=0;by<dim[1];by++)
    {
        temp = Denis_Max(qtnqMat[bx][by].a,qtnqMat[bx][by].b);
        temp = Denis_Max(temp,qtnqMat[bx][by].c);
        (*pdblMatMax)[ind] = Denis_Max(temp,qtnqMat[bx][by].d);
        ind++;
    }
}


//////////////////////////////////////////////////////////////////////////////////////////
//                OPERATION SUR DES VECTEURS QUATERNIONIQUES
//////////////////////////////////////////////////////////////////////////////////////////

//multiplication term by term
void QVectMultTBT(const Quaternion * QVect1,const Quaternion * QVect2,Quaternion ** pQMult,int dim)
{
  int i;
  for(i=0;i<dim;i++)
      (*pQMult)[i]=QMult(QVect1[i],QVect2[i]);
}

int QVectAlloc(Quaternion **pQVect,int dim)
{
  (*pQVect)=(Quaternion *)malloc(dim*sizeof(Quaternion));
  if(*pQVect)  return TRUE;
  else {printf("quaternionic vector allocation error\n"); return FALSE;}
}

void QVectFree(Quaternion **pQVect)
{
  free(*pQVect);
  (*pQVect)=NULL;
}

void QInitVect(Quaternion **pQVect,int dim,Quaternion qtnValue)
{
  int i;
  for(i=0;i<dim;i++)
    (*pQVect)[i]=qtnValue;
}

void QVectScalMult(Quaternion **pVect,int dim,double dblScal)
{
  int i;
  for(i=0;i<dim;i++)
    (*pVect)[i]=QScalMult((*pVect)[i],dblScal);
}

void QVectConj(Quaternion **pVect,Quaternion *Vect,int dim)
{
  int i;
  for(i=0;i<dim;i++)
    (*pVect)[i]=QConj(Vect[i]);
}

void QVectDisp(Quaternion *Vect,int dim,FILE * stream)
{
  int i,j;
  for(i=0;i<dim;i++)
    fprintf(stream,"%4.2lf+(%4.2lf)i+(%4.2lf)j+(%4.2lf)k  ",QReal(Vect[i]),QImI(Vect[i]),QImJ(Vect[i]),QImK(Vect[i]));
  fprintf(stream,"\n");
}






