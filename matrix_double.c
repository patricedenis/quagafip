#include "matrix_double.h"


//memory allocation to place the Image Matrix 
int Matrix_Allocate_Double(int intHeight,int intWidth,double *** pMatrix)
{
  int i,intCount,j;
  
  *pMatrix = (double **)malloc(intHeight*sizeof(double *));
  if ( (*pMatrix) != NULL) /*allocation successful*/
  {
    for(i=0;i<=intHeight-1;i++)
    {
      // place reservation for the three componants RGB
      //(*pMatrix)[i]=(double *)malloc(3*intWidth*sizeof(double));
      // je ne vois pas pourquoi je reservais de la place pour les trois composantes alors je ne laisse
      // qu'une dimension comme on le comprend normalement en utilisant la notion de matrice de double
      
      (*pMatrix)[i]=(double *)malloc(intWidth*sizeof(double));
      if ((*pMatrix)[i] == NULL) //allocation error
      {
				printf("allocation memory error\n");
				//we need to free the matrix memory that is already allocated
				for(j=i-1;j<=0;j--)
				{
	 			  free(((*pMatrix)[j]));
	 			  (*pMatrix)[j]=NULL;
			  }
				free(*pMatrix);
				(*pMatrix)=NULL;
				return FALSE;
      }
    }
  }
  else //allocation error
  {
    printf("allocation memory error\n");
    return FALSE;
  }
  return TRUE;
}

// we need to free the matrix memory 
void Matrix_Free_Double(int intHeight,double *** pMatrix)
{
  int i;
  
  for(i=0;i<=intHeight-1;i++)
      free((*pMatrix)[i]);
  free(*pMatrix);
  (*pMatrix)=NULL;
  printf("double matrix memory free\n");
}

void Matrix_Init_Double(double dblValue,int *intDim,double *** pdblMatrix)
{
    int bx,by;
    
    for (bx=0;bx<intDim[0];bx++)
      for (by=0;by<intDim[1];by++)
        (*pdblMatrix)[bx][by]=dblValue;
}

//cette procedure initialise la matrice avec des valeurs nulles
int DblMatInitZero(double ***pdblMat, int * dim)
{
  int i,j;
  
  if(pdblMat)
  {
    for(i=0;i<dim[0];i++)
	  for(j=0;j<dim[0];j++)
	  (*pdblMat)[i][j] = 0.;
  }
  else return FALSE;
  return TRUE;
}


void DblLogMatrix1D(const double * dblData, int * dim,double ** dblLogData)
{
    int bx,by,ind=0;
    
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        (*dblLogData)[ind]=log(dblData[ind]);
        ind++;
      }  
}

void DblLogMatrix(const double ** dblMat, int * dim,double *** pdblLogMat)
{
    int bx,by,ind=0;
    
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
        (*pdblLogMat)[bx][by]=log(dblMat[bx][by]);
}

void DblLogTab(const double * dblTab, int dim,double ** pdblLogTab)
{
    int bx;
    
    for (bx=0;bx<dim;bx++)
        (*pdblLogTab)[bx]=log(dblTab[bx]);
}

void DblLogP1Matrix(const double ** dblMat, int * dim,double *** pdblLogMat)
{
    int bx,by,ind=0;
    
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
        (*pdblLogMat)[bx][by]=log(dblMat[bx][by]+1)/log(255.);
}
//////////////////////////////////////////////////////////////////////////////

//                  M  I  N   et  M  A  X

//////////////////////////////////////////////////////////////////////////////
 
Dbl_Matrix_Min_Max(const double ** dblMat,int * dim,double * dblMin,double * dblMax)
{
int bx,by;
    
    *dblMin=DBL_MAX;
    *dblMax=DBL_MIN;
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        *dblMin = Denis_Min(*dblMin,dblMat[bx][by]);
        *dblMax = Denis_Max(*dblMax,dblMat[bx][by]);
      }

}


void Dbl_Matrix_Min_Max2(const double * dblMat,int * dim,double * dblMin,double * dblMax)
{
  int bx,by,ind=0;
    
    *dblMin=DBL_MAX;
    *dblMax=DBL_MIN;
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        *dblMin = Denis_Min(*dblMin,dblMat[ind]);
        *dblMax = Denis_Max(*dblMax,dblMat[ind]);
        ind++;
      }

}

//renvoie la matrice dont chaque element est le max des elements des deux matrices en entree
void DblMatMaxed(double *** pdblMatMaxed,const double ** dblMat2,int * dim)
{
  int bx,by;
  for (bx=0;bx<dim[0];bx++)
    for (by=0;by<dim[1];by++)
      (*pdblMatMaxed)[bx][by] = Denis_Max((*pdblMatMaxed)[bx][by],dblMat2[bx][by]);
}

//renvoie la matrice dont chaque element est le max de la matrice d'entree
//et de la moitier de la deuxieme matrice
void DblMatMixed(double *** pdblMatMaxed,const double ** dblMat2,double dblCoeffDiv, int * dim)
{
  int bx,by;
  for (bx=0;bx<dim[0];bx++)
    for (by=0;by<dim[1];by++)
      (*pdblMatMaxed)[bx][by] = Denis_Max((*pdblMatMaxed)[bx][by],dblMat2[bx][by]/dblCoeffDiv);
}


//////////////////////////////////////////////////////////////////////////////

//              C H A N G E M E N T        D ' E C H E L L E 

//////////////////////////////////////////////////////////////////////////////

int Dbl_Mat_ChgScale_IntLvlMax(const double ** dblMat, int * dim,int intLvlMax,
	double dblMin,double dblMax, double *** pdblMatScaled)
{
  int bx,by;
    
  if((dblMax-dblMin) != 0.)
  {
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
        (*pdblMatScaled)[bx][by]= (double)((int)(((dblMat[bx][by] - dblMin)*intLvlMax/(dblMax - dblMin))+ 0.5));
      return TRUE;
  }
  else
  {
    printf("division by zero\n");
    return FALSE;
  }
}

int Dbl_Mat_ChgScaleLog_IntLvlMax(const double ** dblMat, int * dim,int intLvlMax,
	double dblMin,double dblMax, double *** dblMatScaled)
{
  int bx,by;
    
  if((dblMax-dblMin) != 0.)
  {
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
        (*dblMatScaled)[bx][by]= (double)((int)(((log(dblMat[bx][by]+1) - log(dblMin+1))*intLvlMax/(log(dblMax+1) - log(dblMin+1)))+ 0.5));
      return TRUE;
  }
  else
  {
    printf("division by zero\n");
    return FALSE;
  }
}


int Dbl_Mat_ChgScale_IntLvlMax2(const double * dblMat, int * dim,int intLvlMax,
	double dblMin,double dblMax, double ** pdblMatScaled)
{
  int bx,by,ind=0;
    
  if((dblMax-dblMin) != 0.)
  {
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        (*pdblMatScaled)[ind]= (double)((int)(((dblMat[ind] - dblMin)*intLvlMax/(dblMax - dblMin))+ 0.5));
        ind++;
      }
      return TRUE;
  }
  else
  {
    printf("division by zero\n");
    return FALSE;
  }
}


int Dbl_ChgScale_LvlMin_LvlMax(const double dblA,double LvlMin,double LvlMax,
	double dblMin,double dblMax, double * dblScaled)
{
    
  if(((dblMax-dblMin) != 0.)&&((LvlMax-LvlMin) != 0.))
  {
     *dblScaled = ((dblA - dblMin)*(LvlMax-LvlMin)/(dblMax - dblMin));
      return TRUE;
  }
  else
  {
    printf("division by zero\n");
    return FALSE;
  }
}


int Dbl_Mat_ChgScale_LvlMin_LvlMax(const double ** dblMat, int * dim,double LvlMin,double LvlMax,
	double dblMin,double dblMax, double *** dblMatScaled)
{
  int bx,by;
    
  if(((dblMax-dblMin) != 0.)&&((LvlMax-LvlMin) != 0.))
  {
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
        (*dblMatScaled)[bx][by]= ((dblMat[bx][by] - dblMin)*(LvlMax-LvlMin)/(dblMax - dblMin));
      return TRUE;
  }
  else
  {
    printf("division by zero\n");
    return FALSE;
  }
}

void Dbl_Tab_Min_Max(const double * dblTab,int dim,double * dblMin,double * dblMax)
{
  int bx;
    
  *dblMin=DBL_MAX;
  *dblMax=DBL_MIN;
  for (bx=0;bx<dim;bx++)
  {
    *dblMin = Denis_Min(*dblMin,dblTab[bx]);
    *dblMax = Denis_Max(*dblMax,dblTab[bx]);
  }
}

int Dbl_Tab_ChgScale_IntLvlMax(const double * dblTab, int dim,int intLvlMax,
	double dblMin,double dblMax, double ** pdblTabScaled)
{
  int ind;
    
  if((dblMax-dblMin) != 0.)
  {
    for (ind=0;ind<dim;ind++)
      (*pdblTabScaled)[ind]= (double)((int)(((dblTab[ind] - dblMin)*intLvlMax/(dblMax - dblMin))+ 0.5));
    return TRUE;
  }
  else
  {
    printf("division by zero\n");
    return FALSE;
  }
}

//////////////////////////////////////////////////////////////////////////////

//                    G   E   T

//////////////////////////////////////////////////////////////////////////////


void GetDblTabMatrix(double ** dblMat,double ** dblDataTab, int * dim)
{
    int bx,by,ind=0;
    
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        (*dblDataTab)[ind]=dblMat[bx][by];
        ind++;
      }

}

//////////////////////////////////////////////////////////////////////////////

//                    S   E   T

//////////////////////////////////////////////////////////////////////////////

void SetDblMatrixFromTab(double tabSet [],int intHeight,int intWidth,double *** pMatrix)
{
  int i,j,ind=0;
  for(i=0;i<intHeight;i++)
	for(j=0;j<intHeight;j++)
	{
	  (*pMatrix)[i][j] = tabSet[ind];
	  ind++;
	}
}


double VectMean_Double(double * tab, int * dim)
{
  int bx,by,ind=0;
  double sum=0.;
  for (bx=0;bx<dim[0];bx++)
    for (by=0;by<dim[1];by++)
    {
      sum += tab[ind];
      ind++;
    }
  return (double)sum / ((double)(ind));
}

int IsVectDoubleNullEpsilon(double * tab, int dim,double dblEpsilon)
{
  int i;
  for (i=0;i<dim;i++)
    if (tab[i]>=dblEpsilon) return FALSE;
  return TRUE;
}

int DblIsMatNull(double ** dblMat, int * intDim)
{
  int bx,by;

  for (bx=0;bx<intDim[0];bx++)
    for (by=0;by<intDim[1];by++)
    {
      if(dblMat[bx][by]!=0) return FALSE;
    }
  return TRUE;
}

//this function is to apply a threshold on a double matrix
void DblImageThreshold(const double ** dblMatFrom, int * intDim,
  double *** pdblMatTo, double dblThreshold, int bln255)
{
  int bx,by;

  for (bx=0;bx<intDim[0];bx++)
    for (by=0;by<intDim[1];by++)
    {
      if(dblMatFrom[bx][by]>=dblThreshold)
      {
        if (bln255)
          (*pdblMatTo)[bx][by] = 255.;
        else
          (*pdblMatTo)[bx][by] = 1.;
      }
      else (*pdblMatTo)[bx][by] = 0.;
    }
}

void DblImageProductTBT(const double ** dblMatFrom1,const double ** dblMatFrom2,
  int * intDim, double *** pdblMatTo)
{
  int bx,by;

  for (bx=0;bx<intDim[0];bx++)
    for (by=0;by<intDim[1];by++)
    {
      (*pdblMatTo)[bx][by] = dblMatFrom1[bx][by] * dblMatFrom2[bx][by];
    }  
}

//cette fonction inverse les niveaux de gris d'une image,
//le blanc devient noir et le noir devient blanc
void DblImageInverse(const double ** dblMatFrom, int * intDim,
  double *** pdblMatTo)
{
int bx,by;

  for (bx=0;bx<intDim[0];bx++)
    for (by=0;by<intDim[1];by++)
    {
      (*pdblMatTo)[bx][by] = fabs(dblMatFrom[bx][by] - 255.);
    }  
}

//cette fonction effectue un filtrage par la methode de Prewitt horizontale
void DblFiltrePrewittH1(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim)
{
  int bx,by;
  //on gère les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<intDim[0]-1;bx++)
  {
    for(by=1;by<intDim[1]-1;by++)
    {
      (*pdblMatFiltre)[bx][by]= 
      ( dblMatFrom[bx-1][by-1] + dblMatFrom[bx-1][by] + dblMatFrom[bx-1][by+1]
      - dblMatFrom[bx+1][by-1] - dblMatFrom[bx+1][by] - dblMatFrom[bx+1][by+1]);
    }
  }
}

//cette fonction effectue un filtrage par la methode de Prewitt horizontale
void DblFiltrePrewittH2(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim)
{
  int bx,by;
  //on gère les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<intDim[0]-1;bx++)
  {
    for(by=1;by<intDim[1]-1;by++)
    {
      (*pdblMatFiltre)[bx][by]= 
      -( dblMatFrom[bx-1][by-1] - dblMatFrom[bx-1][by] - dblMatFrom[bx-1][by+1]
      + dblMatFrom[bx+1][by-1] + dblMatFrom[bx+1][by] + dblMatFrom[bx+1][by+1]);
    }
  }
}

void DblFiltrePrewittV1(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim)
{
  int bx,by;
  //on gère les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<intDim[0]-1;bx++)
  {
    for(by=1;by<intDim[1]-1;by++)
    {
      (*pdblMatFiltre)[bx][by]= 
      ( dblMatFrom[bx-1][by-1] + dblMatFrom[bx][by-1] + dblMatFrom[bx+1][by-1]
      - dblMatFrom[bx-1][by+1] - dblMatFrom[bx][by+1] - dblMatFrom[bx+1][by+1]);
    }
  }
}

void DblFiltrePrewittV2(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim)
{
  int bx,by;
  //on gère les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<intDim[0]-1;bx++)
  {
    for(by=1;by<intDim[1]-1;by++)
    {
      (*pdblMatFiltre)[bx][by]= 
      -( dblMatFrom[bx-1][by-1] - dblMatFrom[bx][by-1] - dblMatFrom[bx+1][by-1]
      + dblMatFrom[bx-1][by+1] + dblMatFrom[bx][by+1] + dblMatFrom[bx+1][by+1]);
    }
  }
}


void DblFiltrePrewittD1(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim)
{
  int bx,by;
  //on gère les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<intDim[0]-1;bx++)
  {
    for(by=1;by<intDim[1]-1;by++)
    {
      (*pdblMatFiltre)[bx][by]= 
      ( dblMatFrom[bx-1][by] + dblMatFrom[bx-1][by-1] + dblMatFrom[bx][by-1]
      - dblMatFrom[bx][by+1] - dblMatFrom[bx+1][by+1] - dblMatFrom[bx+1][by]);
    }
  }
}

void DblFiltrePrewittD2(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim)
{
  int bx,by;
  //on gère les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<intDim[0]-1;bx++)
  {
    for(by=1;by<intDim[1]-1;by++)
    {
      (*pdblMatFiltre)[bx][by]= 
      ( dblMatFrom[bx][by-1] + dblMatFrom[bx+1][by-1] + dblMatFrom[bx+1][by]
      - dblMatFrom[bx][by+1] - dblMatFrom[bx-1][by+1] - dblMatFrom[bx-1][by]);
    }
  }
}

void DblFiltrePrewittD3(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim)
{
  int bx,by;
  //on gère les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<intDim[0]-1;bx++)
  {
    for(by=1;by<intDim[1]-1;by++)
    {
      (*pdblMatFiltre)[bx][by]= 
      ( - dblMatFrom[bx-1][by] - dblMatFrom[bx-1][by-1] - dblMatFrom[bx][by-1]
      + dblMatFrom[bx][by+1] + dblMatFrom[bx+1][by+1] + dblMatFrom[bx+1][by]);
    }
  }
}


void DblFiltrePrewittD4(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim)
{
  int bx,by;
  //on gère les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<intDim[0]-1;bx++)
  {
    for(by=1;by<intDim[1]-1;by++)
    {
      (*pdblMatFiltre)[bx][by]= 
      ( -dblMatFrom[bx][by-1] - dblMatFrom[bx+1][by-1] - dblMatFrom[bx+1][by]
      + dblMatFrom[bx][by+1] + dblMatFrom[bx-1][by+1] + dblMatFrom[bx-1][by]);
    }
  }
}

//cette procedure effectue la normalisation avec la norme L2
void DblMatL2Normalization(const double ** dblMatH1,const double ** dblMatH2,
  const double ** dblMatV1,const double ** dblMatV2,
  const double ** dblMatD1,const double ** dblMatD2,
  double ***pdblMatNormalized,int * intDim)
{
  int bx,by;
  //la nomre L2 c'est la racine carre de tous les elements au carre
  for(bx=1;bx<intDim[0]-1;bx++)
    for(by=1;by<intDim[1]-1;by++)
      (*pdblMatNormalized)[bx][by] = sqrt(
      pow(dblMatH1[bx][by],2) + pow(dblMatH2[bx][by],2)
      + pow(dblMatV1[bx][by],2) + pow(dblMatV2[bx][by],2)
      + pow(dblMatD1[bx][by],2) + pow(dblMatD2[bx][by],2)
      );
}

//cette fonction est un clone de la précédente mais contient moins de matrices d'entree
void DblMatL2Norme(double *** pdblMatNormeL2,const double ** dblMatH,const double ** dblMatV,
  const double ** dblMatD1,const double ** dblMatD2,int * intDim)
{
  int bx,by;
  //la nomre L2 c'est la racine carre de tous les elements au carre
  for(bx=1;bx<intDim[0]-1;bx++)
    for(by=1;by<intDim[1]-1;by++)
      (*pdblMatNormeL2)[bx][by] = sqrt(
      pow(dblMatH[bx][by],2) + pow(dblMatV[bx][by],2)
      + pow(dblMatD1[bx][by],2) + pow(dblMatD2[bx][by],2)
      );
}

void DblFiltrePrewittTotal(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim)
{
  int bx,by;
  //on gère les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  double ** dblMatPrewittH,**dblMatPrewittV,**dblMatPrewittD1,**dblMatPrewittD2,
    ** dblMatPrewittH2,**dblMatPrewittV2,**dblMatPrewittD3,**dblMatPrewittD4;
  
  Matrix_Allocate_Double(intDim[0],intDim[1],&dblMatPrewittH);
  Matrix_Allocate_Double(intDim[0],intDim[1],&dblMatPrewittV);
  Matrix_Allocate_Double(intDim[0],intDim[1],&dblMatPrewittH2);
  Matrix_Allocate_Double(intDim[0],intDim[1],&dblMatPrewittV2);
  Matrix_Allocate_Double(intDim[0],intDim[1],&dblMatPrewittD1);
  Matrix_Allocate_Double(intDim[0],intDim[1],&dblMatPrewittD2);
  Matrix_Allocate_Double(intDim[0],intDim[1],&dblMatPrewittD4);
  Matrix_Allocate_Double(intDim[0],intDim[1],&dblMatPrewittD3);
  
  //on clacule pour chacune des direction un gradient par filtre de Prewitt
  DblFiltrePrewittH1((const double **)dblMatFrom,&dblMatPrewittH,intDim);
  DblFiltrePrewittV1((const double **)dblMatFrom,&dblMatPrewittV,intDim);
  DblFiltrePrewittH2((const double **)dblMatFrom,&dblMatPrewittH2,intDim);
  DblFiltrePrewittV2((const double **)dblMatFrom,&dblMatPrewittV2,intDim);
  DblFiltrePrewittD1((const double **)dblMatFrom,&dblMatPrewittD1,intDim);
  DblFiltrePrewittD2((const double **)dblMatFrom,&dblMatPrewittD2,intDim);
  DblFiltrePrewittD3((const double **)dblMatFrom,&dblMatPrewittD3,intDim);
  DblFiltrePrewittD4((const double **)dblMatFrom,&dblMatPrewittD4,intDim);
  
  //on calcule ensuite le gradient total avec la norme Linf : le max
  Matrix_Init_Double(0.,intDim,pdblMatFiltre);
  DblMatMaxed(&(*pdblMatFiltre),(const double **)dblMatPrewittH,intDim);
  DblMatMaxed(&(*pdblMatFiltre),(const double **)dblMatPrewittH2,intDim);
  DblMatMaxed(&(*pdblMatFiltre),(const double **)dblMatPrewittV,intDim);
  DblMatMaxed(&(*pdblMatFiltre),(const double **)dblMatPrewittV2,intDim);
  DblMatMaxed(&(*pdblMatFiltre),(const double **)dblMatPrewittD1,intDim);
  DblMatMaxed(&(*pdblMatFiltre),(const double **)dblMatPrewittD2,intDim);
  DblMatMaxed(&(*pdblMatFiltre),(const double **)dblMatPrewittD3,intDim);
  DblMatMaxed(&(*pdblMatFiltre),(const double **)dblMatPrewittD4,intDim);
  
  /*DblMatL2Normalization((const double **)dblMatPrewittH,(const double **)dblMatPrewittH2,
    (const double **)dblMatPrewittV,(const double **)dblMatPrewittV2,
    (const double **)dblMatPrewittD1,(const double **)dblMatPrewittD2,
    &(*pdblMatFiltre),intDim);*/

  
  //ensuite on libere les matrice temporaires
  Matrix_Free_Double(intDim[0],&dblMatPrewittH);
  Matrix_Free_Double(intDim[0],&dblMatPrewittV);
  Matrix_Free_Double(intDim[0],&dblMatPrewittH2);
  Matrix_Free_Double(intDim[0],&dblMatPrewittV2);
  Matrix_Free_Double(intDim[0],&dblMatPrewittD1);
  Matrix_Free_Double(intDim[0],&dblMatPrewittD2);
  Matrix_Free_Double(intDim[0],&dblMatPrewittD3);
  Matrix_Free_Double(intDim[0],&dblMatPrewittD4);}


//cette fonction effectue un filtrage par la methode de Prewitt horizontale
void DblFiltreSobelH(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim)
{
  int bx,by;
  //on gère les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Sobel 3*3 puisse se faire sans probleme.
  for(bx=1;bx<intDim[0]-1;bx++)
  {
    for(by=1;by<intDim[1]-1;by++)
    {
      (*pdblMatFiltre)[bx][by]= 1./8. * 
      ( dblMatFrom[bx-1][by-1] + 2*dblMatFrom[bx-1][by] + dblMatFrom[bx-1][by+1]
      - dblMatFrom[bx+1][by-1] - 2*dblMatFrom[bx+1][by] - dblMatFrom[bx+1][by+1]);
    }
  }
}

void DblFiltreSobelV(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim)
{
  int bx,by;
  //on gère les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Sobel 3*3 puisse se faire sans probleme.
  for(bx=1;bx<intDim[0]-1;bx++)
  {
    for(by=1;by<intDim[1]-1;by++)
    {
      (*pdblMatFiltre)[bx][by]= 1./8. * 
      ( dblMatFrom[bx-1][by-1] + 2*dblMatFrom[bx][by-1] + dblMatFrom[bx+1][by-1]
      - dblMatFrom[bx-1][by+1] - 2*dblMatFrom[bx][by+1] - dblMatFrom[bx+1][by+1]);
    }
  }
}

void DblFiltreSobelD1(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim)
{
  int bx,by;
  //on gère les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Sobel 3*3 puisse se faire sans probleme.
  for(bx=1;bx<intDim[0]-1;bx++)
  {
    for(by=1;by<intDim[1]-1;by++)
    {
      (*pdblMatFiltre)[bx][by]= 1./8. * 
      ( dblMatFrom[bx-1][by] + 2*dblMatFrom[bx-1][by-1] + dblMatFrom[bx][by-1]
      - dblMatFrom[bx][by+1] - 2*dblMatFrom[bx+1][by+1] - dblMatFrom[bx+1][by]);
    }
  }
}

void DblFiltreSobelD2(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim)
{
  int bx,by;
  //on gère les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Sobel 3*3 puisse se faire sans probleme.
  for(bx=1;bx<intDim[0]-1;bx++)
  {
    for(by=1;by<intDim[1]-1;by++)
    {
      (*pdblMatFiltre)[bx][by]= 1./8. * 
      ( dblMatFrom[bx][by-1] + 2*dblMatFrom[bx+1][by-1] + dblMatFrom[bx+1][by]
      - dblMatFrom[bx][by+1] - 2*dblMatFrom[bx+1][by+1] - dblMatFrom[bx+1][by]);
    }
  }
}



void DblFiltreSobelTotal(const double ** dblMatFrom,double ***pdblMatFiltre,int * intDim)
{
  int bx,by;
  //on gère les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  double ** dblMatPrewittH,**dblMatPrewittV,**dblMatPrewittD1,**dblMatPrewittD2;
  
  Matrix_Allocate_Double(intDim[0],intDim[1],&dblMatPrewittH);
  Matrix_Allocate_Double(intDim[0],intDim[1],&dblMatPrewittV);
  Matrix_Allocate_Double(intDim[0],intDim[1],&dblMatPrewittD1);
  Matrix_Allocate_Double(intDim[0],intDim[1],&dblMatPrewittD2);
  //on clacule pour chacune des direction un gradient par filtre de Prewitt
  DblFiltreSobelH((const double **)dblMatFrom,&dblMatPrewittH,intDim);
  DblFiltreSobelV((const double **)dblMatFrom,&dblMatPrewittV,intDim);
  DblFiltreSobelD1((const double **)dblMatFrom,&dblMatPrewittD1,intDim);
  DblFiltreSobelD2((const double **)dblMatFrom,&dblMatPrewittD2,intDim);

  //on calcule ensuite le gradient total avec la norme Linf : le max
  Matrix_Init_Double(0.,intDim,pdblMatFiltre);
  DblMatMaxed(&(*pdblMatFiltre),(const double **)dblMatPrewittH,intDim);
  DblMatMaxed(&(*pdblMatFiltre),(const double **)dblMatPrewittV,intDim);
  DblMatMaxed(&(*pdblMatFiltre),(const double **)dblMatPrewittD1,intDim);
  DblMatMaxed(&(*pdblMatFiltre),(const double **)dblMatPrewittD2,intDim);
  
  //ensuite on libere les matrice temporaires
  Matrix_Free_Double(intDim[0],&dblMatPrewittH);
  Matrix_Free_Double(intDim[0],&dblMatPrewittV);
  Matrix_Free_Double(intDim[0],&dblMatPrewittD1);
  Matrix_Free_Double(intDim[0],&dblMatPrewittD2);
}



void DblGaussienne2D(double dblAmp, double dblSigma, double * dblCoordCentre, int * intDim, double *** pDblTab,int blnInvert)
{
  int i,j;
  for(i=0;i<intDim[0];i++)
    for(j=0;j<intDim[1];j++)
    {
      if(blnInvert==FALSE)
      {
        (*pDblTab)[i][j] = dblAmp * exp( - 
        (  ( (i-dblCoordCentre[0])*(i-dblCoordCentre[0])
            +(j-dblCoordCentre[1])*(j-dblCoordCentre[1]) )
        / 2*dblSigma*dblSigma ) );
        }
        else {
          (*pDblTab)[i][j] = dblAmp - dblAmp * exp( - 
        (  ( (i-dblCoordCentre[0])*(i-dblCoordCentre[0])
            +(j-dblCoordCentre[1])*(j-dblCoordCentre[1]) )
        / 2*dblSigma*dblSigma ) );
        }
    }
}

        
        
//  This function will shift the contents of the two sub-matrix pointed by
//  MatFrom and MatTo as describe in the following scheme
//		|-----------|				|-----------|
//		|  1  |  2  |				|  4  |  3  |
//		|-----------| --->  |-----------|
//		|  3  |  4  |				|  2  |  1  |
//		|-----------|				|-----------|
void DblCopySubMatrix(double ** DblMatFrom,double *** pDblMatTo,int intXFrom,int intYFrom,
int intXTo,int intYTo,int intHeight,int intWidth)
{
	int i,j;
	for(i=0;i<intHeight;i++)
		for(j=0;j<intWidth;j++)
			(*pDblMatTo)[i+intXTo][j+intYTo] = DblMatFrom[i+intXFrom][j+intYFrom];
}

void DblMatrixShift(double ** DblMatFrom,double *** pDblMatShifted,int * dim)
{
	int intMid;
	
	intMid = dim[0]/2;
	DblCopySubMatrix(DblMatFrom,pDblMatShifted,0,0,intMid,intMid,intMid,intMid);	//1->4
	DblCopySubMatrix(DblMatFrom,pDblMatShifted,intMid,intMid,0,0,intMid,intMid);	//4->1
	DblCopySubMatrix(DblMatFrom,pDblMatShifted,intMid,0,0,intMid,intMid,intMid);	//2->3
	DblCopySubMatrix(DblMatFrom,pDblMatShifted,0,intMid,intMid,0,intMid,intMid);	//3->2
}

void test_Prewitt(char * strFile)
{
  double ** dblMatImage,**dblGradient;
  int dim[2];
char * strGradient;
   MgkLireImgGrisAndAllocAndDblSetMatrix(strFile,dim,&dblMatImage);
   
    //ensuite il faut creer un gradient sur cette partie
    Matrix_Allocate_Double(dim[0],dim[1],&dblGradient);
    DblFiltrePrewittTotal((const double **)dblMatImage,&dblGradient,dim);
    Matrix_Free_Double(dim[0],&dblMatImage);
    //on enregistre ce gradient
    strGradient = (char *)malloc((30+strlen(strFile))*sizeof(char));
    strcpybtw(&strGradient,strFile,"-GradientPrewitt",'.');
    MgkDblGetMatrixAndEcrireImgGris(dim,(const double **)dblGradient,strGradient,stdout);
    free(strGradient);strGradient=NULL;  
    Matrix_Free_Double(dim[0],&dblGradient);
   
}

int main_matrix_double(int argc, char *argv[])
{
    test_Prewitt("C:/images/G3Algebra/filtre/DenisBerthier/synthese2-1024NG.bmp");
    return 1;
}
