#include "filtrage_spectral_quaternionique.h"

//pour toute nos fenetres on considerera que la dimmension est un nombre pair

//cette procedure cree une fenetre de Hamming de la dimmension souhaitee
//la fenetre de Hamming s'etendra donc sur toute la largeur
void dblHammingWindow1D(double * pdblHammingVect,int intDim)
{
  int i;
  
  //comme la fenetre de Hamming sera utilisée pour etre multiplié à un signal
  //frequentiel ou spatial, on va la coder sur un tableau de reel.
  // la formule est donnée par h(t) = 0.54 - 0.46 * cos(2* MATH_PI t/T);
  //avec un t compris dans [0 , T]
  for(i=0;i<intDim;i++)
  {
    (pdblHammingVect)[i]= 0.54 - 0.46 * cos(2*M_PI * i /(double)(intDim-1));
  }
}

//construction d'une fenêtre rectangulaire
//la fenetre s'etendra donc sur toute la largeur
void dblRectangularWindow1D(double * pdblRectVect,int intDim)
{
  int i;
  
  for(i=0;i<intDim;i++)
  {
    (pdblRectVect)[i]= 1.0;
  }
}

//cette procedure cree une fenetre de Hamming sur un intervalle donné 
//en fonction de la dimension du vecteur
//intNbCoupe represente le nombre de fois que l'on coupe la tranche dim
void dblHammingWindow1DIntervalle(double * pdblHammingVect,int intDim, int intNbCoupe)
{
  int i;
  int intHammingDim,intHammingDebut;
  
  //on calcule sur quel intervalle on appliquera la fenetre de Hamming sur le vecteur
  intHammingDim = intDim / intNbCoupe + 1;
  intHammingDebut = intDim/2 - intHammingDim/2 ;
  for(i=0;i<intDim;i++)
  {
    (pdblHammingVect)[i]= 0.;
  }
  for(i=intHammingDebut;i<intHammingDebut+intHammingDim;i++)
  {
    (pdblHammingVect)[i]= 1.-(0.54 - 0.46 * cos(2*M_PI * i /(double)(intHammingDim-1)));
  }
}

//cree une fenetre carrée de Hamming 2d 
void dblHammingSquareWindow2DIntervalle(double *** pdblHamming2D,int intDim, int intNbCoupe)
{
  int i,j;
  double * dblHamming1D;
  
  //on cree le vecteur comportant la porte de Hamming
  dblHamming1D = (double *)malloc(intDim*sizeof(double));
  dblHammingWindow1DIntervalle(dblHamming1D,intDim,intNbCoupe);     

  //on recopie separement
  //sur les lignes
  for(i=0;i<intDim;i++)
  {
    //et sur les colonnes
    for(j=0;j<intDim;j++)
      (*pdblHamming2D)[i][j] = dblHamming1D[j];
    (*pdblHamming2D)[i][j] = (*pdblHamming2D)[i][j] * dblHamming1D[i];  
  }
  free(dblHamming1D);
  dblHamming1D = NULL;
}


//on veut pouvoir effectuer l'operation de switch sur les vecteurs 1D qui consiste a faire:
//		|-----------|				|-----------|
//		|  1  |  2  |				|  2  |  1  |
//		|-----------| --->  |-----------|
void dblCopySubVect(double * DblVectFrom,double ** pDblVectTo,int intXFrom,int intXTo,int intDim)
{
	int i,j;
	for(i=0;i<intDim;i++)
			(*pDblVectTo)[i+intXTo] = DblVectFrom[i+intXFrom];
}

void DblVectShift(double * DblVectFrom,double ** pDblVectShifted,int dim)
{
	int intMid;
	
	intMid = dim/2;
	dblCopySubVect(DblVectFrom,pDblVectShifted,0,intMid,intMid);	//1->2
	dblCopySubVect(DblVectFrom,pDblVectShifted,intMid,0,intMid);	  //2->1
}



// nbcoupe permet de savoir comment sera initialisé le spectre
// en effet si nbcoupe = 4 il y aura un quart du spectre qui sera 
// initialiser a une certaine valeur, le reste du spectre restant nul.
void QFiltreFreqPasseBasHamming(char * strFileIn,char * strFileOut,int nbcoupe)
{
  int dim[2],type;
  //Complex ** Mat;
  Quaternion ** qtnImage,** qtnQMat;
  Quaternion ** qtnMatLL, ** qtnFiltrePasseBas, ** qtnFll,**qtnFllC, **qtnFLL;

  Quaternion QMu1,QMu2;
    
  double dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax;
  
  double ** dblHamming2D;

  QMu1 = QInit(0.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)); //Mu luminance
  QMu2 = QInit(0.,0.,1./sqrt(2.),-1./sqrt(2.)); //mu2 donné dans l'article de Sangwine et Ell
    
  //allocation matrice Quaternion puis copie
  MgkLireImgCouleurAndAllocAndQSetMatrix(strFileIn,dim,&qtnImage);

  //allocation d'un filtre reel pour la fenetre de Hamming
  Matrix_Allocate_Double(dim[0],dim[1],&dblHamming2D);
  dblHammingSquareWindow2DIntervalle(&dblHamming2D,dim[0],nbcoupe);
  
  //on calcule la TF
  QMatrixAllocate(dim[0],dim[1],&qtnQMat);
  QFFT(qtnImage,dim,&qtnQMat,QMu1,QMu2,1);
  QMatrixFree(dim[0],&qtnImage);

  //il faut maintenant multiplier le spectre par la fenetre de Hamming
  //on alloue pour filtrer le spectre
  QMatrixAllocate(dim[0],dim[1],&qtnFLL);
  QMatMultTBTWithDblMat(qtnQMat,dblHamming2D,&qtnFLL,dim);  
  printf("%lf\n",dblHamming2D[dim[0]/2][dim[1]/2]);
  //ca peche ici ...
  QMatrixFree(dim[0],&qtnQMat);
  printf("%lf\n",dblHamming2D[dim[0]/2][dim[1]/2]);
  Matrix_Free_Double(dim[0],&dblHamming2D);
  printf("%lf\n",dblHamming2D[dim[0]/2][dim[1]/2]);    

  QMatrixAllocate(dim[0],dim[1],&qtnFll);
  QFFT(qtnFLL,dim,&qtnFll,QMu1,QMu2,-1);
  
  QMat_Min_Max((const Quaternion **) qtnFll,dim,&dblRMin,&dblRMax,
    &dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
  printf("rmin%lf rmax%lf\nimin%lf imax%lf\njmin%lf jmax%lf\nkmin%lf kmax%lf\n",
      dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax);
  
  //on caste
  QMatrixAllocate(dim[0],dim[1],&qtnFllC);
  QMatCast0255((const Quaternion **)qtnFll,&qtnFllC,dim);
  QMatrixFree(dim[0],&qtnFll);
  
  MgkQGetMatrixAndEcrireImgCouleur(dim,qtnFllC,strFileOut);
  QMatrixFree(dim[0],&qtnFllC);
  
  //QMatrixFree(dim[0],&qtnQMat);
}


/*
void QFiltreFreqPasseHaut(char * strFileIn,char * strFileOut,int nbcoupe)
{
  int dim[2],type;
  //Complex ** Mat;
  Quaternion ** qtnqMat,** qtnQMat=NULL;
  Quaternion ** qtnMatHH, ** qtnFHH, ** qtnFhh,**qtnFhhC;

  Quaternion QMu1,QMu2;
  
  double dblEpsilon = 0.000000001;
  double dblRMin=0,dblRMax=0,dblIMin=0,dblIMax=0,dblJMin=0,dblJMax=0,dblKMin=0,dblKMax=0;
  double dblMin=0,dblMax=0;

  QMu1 = QInit(0.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)); //Mu luminance
  QMu2 = QInit(0.,0.,1./sqrt(2.),-1./sqrt(2.)); //mu2 donné dans l'article de Sangwine et Ell
    
  //allocation matrice Quaternion puis copie
  MgkLireImgCouleurAndAllocAndQSetMatrix(strFileIn,dim,&qtnqMat);

//QMat_Min_Max((const Quaternion **) qtnqMat,dim,&dblRMin,&dblRMax,
//    &dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
//  printf("rmin%lf rmax%lf\nimin%lf imax%lf\njmin%lf jmax%lf\nkmin%lf kmax%lf\n",
//      dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax);

  //allocation des differents filtres
  QMatrixAllocate(dim[0],dim[1],&qtnMatHH);
  QInitMatrix(dim[0],dim[1],QInit(0.0,0.,0.,0.),&qtnMatHH);
  //remplissage
  QMatHHFilter2(&qtnMatHH,dim,nbcoupe);
  
//  QMat_Min_Max((const Quaternion **) qtnMatHH,dim,&dblRMin,&dblRMax,
//    &dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
//  printf("3rmin%lf rmax%lf\nimin%lf imax%lf\njmin%lf jmax%lf\nkmin%lf kmax%lf\n",
//      dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax);
  
  //on calcule la TF
  QMatrixAllocate(dim[0],dim[1],&qtnQMat);
  QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&qtnQMat);
  //fft_2d(qtncMat,dim,qtnQMat,1);
  QFFT(qtnqMat,dim,&qtnQMat,QMu1,QMu2,1);
  QMatrixFree(dim[0],&qtnqMat);
  
//  QMat_Min_Max((const Quaternion **) qtnQMat,dim,&dblRMin,&dblRMax,
//    &dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
//  printf("four rmin%lf rmax%lf\nimin%lf imax%lf\njmin%lf jmax%lf\nkmin%lf kmax%lf\n",
//      dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax);

  //on alloue pour filtrer le spectre
  QMatrixAllocate(dim[0],dim[1],&qtnFHH);
  QInitMatrix(dim[0],dim[1],QInit(0.0,0.,0.,0.),&qtnFHH);
  QMatMultTBT(qtnQMat,qtnMatHH,&qtnFHH,dim);
  
//  QMat_Min_Max((const Quaternion **) qtnFHH,dim,&dblRMin,&dblRMax,
//    &dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
//  printf("rmin%lf rmax%lf\nimin%lf imax%lf\njmin%lf jmax%lf\nkmin%lf kmax%lf\n",
//      dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax);
      
  QMatrixFree(dim[0],&qtnQMat);
  QMatrixFree(dim[0],&qtnMatHH);
  
  QMatrixAllocate(dim[0],dim[1],&qtnFhh);
  QInitMatrix(dim[0],dim[1],QInit(0.0,0.,0.,0.),&qtnFhh);
  QFFT(qtnFHH,dim,&qtnFhh,QMu1,QMu2,-1);
  
  QMat_Min_Max((const Quaternion **) qtnFhh,dim,&dblRMin,&dblRMax,
    &dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
  printf("rmin%lf rmax%lf\nimin%lf imax%lf\njmin%lf jmax%lf\nkmin%lf kmax%lf\n",
      dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax);
  
  QMat_Min_Max2((const Quaternion **) qtnFhh,dim,&dblMin,&dblMax);
  QMatrixAllocate(dim[0],dim[1],&qtnFhhC);
  QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&qtnFhhC);
  QImMat_ChgScale_IntLvlMax((const Quaternion **) qtnFhh,dim,255,
    dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax,&qtnFhhC);
  QMatrixFree(dim[0],&qtnFhh);
  
  QMat_Min_Max((const Quaternion **) qtnFhhC,dim,&dblRMin,&dblRMax,
    &dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
  printf("scale rmin%lf rmax%lf\nimin%lf imax%lf\njmin%lf jmax%lf\nkmin%lf kmax%lf\n",
      dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax);
  
  MgkQGetMatrixAndEcrireImgCouleur(dim,qtnFhhC,strFileOut);
  QMatrixFree(dim[0],&qtnFhhC);

} 


void QAbsolue(char * strFile,char * strFileOut)
{
  int dim[2],type;
  //Complex ** Mat;
  Quaternion ** qtnqMat=NULL,** qtnMatAbs=NULL;
  double *dblMatMax=NULL, *dblMatScaled=NULL;
  double dblMin=0,dblMax=0;

  //allocation matrice Quaternion puis copie
  MgkLireImgCouleurAndAllocAndQSetMatrix(strFile,dim,&qtnqMat);
  //on recupere les valeurs absolues
  QMatrixAllocate(dim[0],dim[1],&qtnMatAbs);
  QMatAbs((const Quaternion **)qtnqMat,dim,&qtnMatAbs);
  QMatrixFree(dim[0],&qtnqMat);
  
  //on recupere une matrice de double max des valeurs absolues
  dblMatMax=(double *)malloc(dim[0]*dim[1]*sizeof(double));
  QMatMax2((const Quaternion **)qtnMatAbs,dim,&dblMatMax);
  QMatrixFree(dim[0],&qtnMatAbs);
  
  //on remet a l'echelle entre 0 et 255
  
  dblMatScaled=(double *)malloc(dim[0]*dim[1]*sizeof(double));
  Dbl_Matrix_Min_Max2(dblMatMax,dim,&dblMin,&dblMax);
  Dbl_Mat_ChgScale_IntLvlMax2((const double *)dblMatMax,dim,255,dblMin,dblMax,&dblMatScaled);
	free(dblMatMax);dblMatMax=NULL;

  //on enregistre l'image
  MgkEcrireImgGris(strFileOut,dblMatScaled,dim);
  free(dblMatScaled);dblMatScaled=NULL;
}
*/

void affichevect(double * vect,int dim)
{
  int i;
  for(i=0;i<dim;i++)
    printf("%lf ",vect[i]);
   printf("\n");     
}

int main_filtrage_spectral_quaternionique(int argc, char *argv[])
{
  #define DIM 16
  int nb=8;//on divise le spectre en huit parties
  double hamming[DIM];
  double rect[DIM];  
 

  //on teste les fenetres -> OK fonctionnent
  //dblRectangularWindow(rect,DIM);
  //affichevect(rect,DIM);
  //dblHammingWindow1D(hamming,DIM);  
  //affichevect(hamming,DIM);
  
  //ensuite on creer des filtres 1D
  //dblHammingWindow1DIntervalle(hamming,DIM,2);
  //affichevect(hamming,DIM);
  

  QFiltreFreqPasseBasHamming("D:/images/base/formes.bmp","D:/images/filtre/quaternions/spectral/formesFiltreHaut.bmp",nb);
  //QAbsolue("D:/images/filtre/freq/formesFiltreHaut.bmp","D:/images/filtre/freq/formesFiltreAbs.bmp");

  //ensuite on fait le meme filtrage avec une fenetre de Hamming
  return 0;
}

