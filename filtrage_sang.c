#include "filtrage_sang.h"

//
//QIMAGE QIMAGE::SgwH (quaternion mu)
//{
//	QIMAGE qiTemp(n,m);
//	unsigned int x,y;
//	
//	for (y=1;y<m-1;y++)
//	for (x=1;x<n-1;x++)
//	{
//		qiTemp.image[x+y*n] = image[x-1+(y-1)*n] + image[x-1+y*n] + image[x-1+(y+1)*n] +
//				     mu * (image[x+1+(y-1)*n] + image[x+1+y*n] + image[x+1+(y+1)*n]) * (~mu);
//	}
//		
//	return qiTemp;
//}
//
//QIMAGE QIMAGE::SgwV (quaternion mu)
//{
//	QIMAGE qiTemp(n,m);
//	unsigned int x,y;
//	
//	for (y=1;y<m-1;y++)
//	for (x=1;x<n-1;x++)
//	{
//		qiTemp.image[x+y*n] = image[x-1+(y-1)*n] + image[x+(y-1)*n] + image[x+1+(y-1)*n] +
//				     mu * (image[x-1+(y+1)*n] + image[x+(y+1)*n] + image[x+1+(y+1)*n]) * (~mu);
//	}
//		
//	return qiTemp;
//}

//cette fonction permet d'appliquer un filtre de prewitt quaternionique horizontal
//comme decrit dans l'article de S Sangwine.
// attention a ne pas oublier d'alouer et d'initialiser
//la matrice QMatFiltre à 0 avant de l'appeler dans cette fonction
void filtre_Sang_H(Quaternion ** QMat,Quaternion ***pQMatFiltre, Quaternion QMu,int * dim)
{
  int bx,by,ind=0;
  //on gere les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<dim[0]-1;bx++)
  {
    for(by=1;by<dim[1]-1;by++)
    {
      (*pQMatFiltre)[bx][by]= QScalMult( QAdd(
          QAdd(QAdd(QMat[bx-1][by-1],QMat[bx-1][by]),QMat[bx-1][by+1]) ,
          QMult(  QMult(QMu,(QAdd(QAdd(QMat[bx+1][by-1],QMat[bx+1][by]),QMat[bx+1][by+1])) ),QInv(QMu))
            ),1./6.);
    }
  }
}

void filtre_Sang_Hbis(Quaternion ** QMat,Quaternion ***pQMatFiltre, Quaternion QMu,int * dim)
{
  int bx,by,ind=0;
  //on gere les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<dim[0]-1;bx++)
  {
    for(by=1;by<dim[1]-1;by++)
    {
      (*pQMatFiltre)[bx][by]= QScalMult( QAdd(
          QAdd(QAdd(QMat[bx+1][by-1],QMat[bx+1][by]),QMat[bx+1][by+1]) ,
          QMult(  QMult(QMu,(QAdd(QAdd(QMat[bx-1][by-1],QMat[bx-1][by]),QMat[bx-1][by+1])) ),QInv(QMu))
            ),1./6.);
    }
  }
}

//cette fonction permet d'appliquer un filtre de prewitt quaternionique vertical
//comme decrit dans l'article de S Sangwine.
// attention a ne pas oublier d'alouer et d'initialiser
//la matrice QMatFiltre à 0 avant de l'appeler dans cette fonction
void filtre_Sang_V(Quaternion ** QMat,Quaternion ***pQMatFiltre, Quaternion QMu,int * dim)
{
  int bx,by,ind=0;
  //on gere les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<dim[0]-1;bx++)
  {
    for(by=1;by<dim[1]-1;by++)
    {
      (*pQMatFiltre)[bx][by]= QScalMult( QAdd(
          QAdd(QAdd(QMat[bx-1][by-1],QMat[bx][by-1]),QMat[bx+1][by-1]) ,
          QMult(  QMult(QMu,(QAdd(QAdd(QMat[bx-1][by+1],QMat[bx][by+1]),QMat[bx+1][by+1])) ),QInv(QMu))
            ),1./6.);
    }
  }
}

void filtre_Sang_Vbis(Quaternion ** QMat,Quaternion ***pQMatFiltre, Quaternion QMu,int * dim)
{
  int bx,by,ind=0;
  //on gere les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<dim[0]-1;bx++)
  {
    for(by=1;by<dim[1]-1;by++)
    {
      (*pQMatFiltre)[bx][by]= QScalMult( QAdd(
          QAdd(QAdd(QMat[bx-1][by+1],QMat[bx][by+1]),QMat[bx+1][by+1]) ,
          QMult(  QMult(QMu,(QAdd(QAdd(QMat[bx-1][by-1],QMat[bx][by-1]),QMat[bx+1][by-1])) ),QInv(QMu))
            ),1./6.);
    }
  }
}

void filtre_Sang_D1(Quaternion ** QMat,Quaternion ***pQMatFiltre, Quaternion QMu,int * dim)
{
  int bx,by,ind=0;
  //on gere les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<dim[0]-1;bx++)
  {
    for(by=1;by<dim[1]-1;by++)
    {
      (*pQMatFiltre)[bx][by]= QScalMult( QAdd(
          QAdd(QAdd(QMat[bx-1][by],QMat[bx-1][by-1]),QMat[bx][by-1]) ,
          QMult(  QMult(QMu,(QAdd(QAdd(QMat[bx][by+1],QMat[bx+1][by+1]),QMat[bx+1][by])) ),QInv(QMu))
            ),1./6.);
    }
  }
}

void filtre_Sang_D1bis(Quaternion ** QMat,Quaternion ***pQMatFiltre, Quaternion QMu,int * dim)
{
  int bx,by,ind=0;
  //on gere les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<dim[0]-1;bx++)
  {
    for(by=1;by<dim[1]-1;by++)
    {
      (*pQMatFiltre)[bx][by]= QScalMult( QAdd(
          QAdd(QAdd(QMat[bx][by+1],QMat[bx+1][by+1]),QMat[bx+1][by]) ,
          QMult(  QMult(QMu,(QAdd(QAdd(QMat[bx-1][by],QMat[bx-1][by-1]),QMat[bx][by-1])) ),QInv(QMu))
            ),1./6.);
    }
  }
}

void filtre_Sang_D2(Quaternion ** QMat,Quaternion ***pQMatFiltre, Quaternion QMu,int * dim)
{
  int bx,by,ind=0;
  //on gere les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<dim[0]-1;bx++)
  {
    for(by=1;by<dim[1]-1;by++)
    {
      (*pQMatFiltre)[bx][by]= QScalMult( QAdd(
          QAdd(QAdd(QMat[bx][by-1],QMat[bx+1][by-1]),QMat[bx+1][by]) ,
          QMult(  QMult(QMu,(QAdd(QAdd(QMat[bx][by+1],QMat[bx+1][by+1]),QMat[bx+1][by])) ),QInv(QMu))
            ),1./6.);
    }
  }
}

void filtre_Sang_D2bis(Quaternion ** QMat,Quaternion ***pQMatFiltre, Quaternion QMu,int * dim)
{
  int bx,by,ind=0;
  //on gere les effets de bords
  //les premieres et dernieres lignes et colonnes de l'image
  //ne sont pas prises en compte pour que la convolution de
  //notre filtre de Prewitt 3*3 puisse se faire sans probleme.
  for(bx=1;bx<dim[0]-1;bx++)
  {
    for(by=1;by<dim[1]-1;by++)
    {
      (*pQMatFiltre)[bx][by]= QScalMult( QAdd(
          QAdd(QAdd(QMat[bx][by+1],QMat[bx+1][by+1]),QMat[bx+1][by]) ,
          QMult(  QMult(QMu,(QAdd(QAdd(QMat[bx][by-1],QMat[bx+1][by-1]),QMat[bx+1][by])) ),QInv(QMu))
            ),1./6.);
    }
  }
}

void distance_couleur(Quaternion ** QMat,double ***pdblMatDist, Quaternion QMu,int * dim)
{
  int bx,by,ind=0;
  Quaternion QTemp;
  double dblEpsilon = 0.000000001;
  for(bx=1;bx<dim[0]-1;bx++)
  {
    for(by=1;by<dim[1]-1;by++)
    {
      QTemp = QScalMult(QAdd(QMat[bx][by],QMult(QMult(QMu,QMat[bx][by]),QMu)),1./2.);
      //verifier que la matrice est strictement reelle
      //if((QImI(QTemp)>=dblEpsilon)||(QImJ(QTemp)>=dblEpsilon)||(QImK(QTemp)>=dblEpsilon))
//        printf("QTemp non reel\n");
      (*pdblMatDist)[bx][by] = QNorm(QTemp);
    }
  }
}

//effectue un filtrage de Prewitt quaternionique sur l'image passée en paramatre.
//le paramatre intSens permet d'indiquer si le filtrage de Prewitt est horizontal:1
//ou vertical:2
void filtrage_image(char * strFile,int intSens, char * strFileOut)
{
  int dim[2],type;
  double *r=NULL,*g=NULL,*b=NULL;
  double dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax;
  Quaternion **QMatrix=NULL,**QMatFiltre=NULL,**QMatScaled;

  Quaternion QMu,QExpMuFoisPiSur2;
  
  printf("Filtrage de Prewitt Quaternionique (Sangwine)\n");
  if ((intSens == 1)||(intSens == 2)||(intSens == 3)||(intSens == 4)) // bon parametre de sens
  {
    //on remplit la matrice quaternionique correspond a l'image du chemin passé en parametre
    //lecture de l'en-tete
    MgkTypeImage(strFile,dim,&type);
  
    //allocation des tableaux r,g et b
    r = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    g = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    b = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  
    //remplissage des tableaux
    MgkLireImgCouleur(strFile,r,g,b);
  
    //copie dans matrice quaternionique
    QMatrixAllocate(dim[0],dim[1],&QMatrix);
    QSetMatrix(r,g,b,dim[0],dim[1],&QMatrix);
    free(r);r=NULL;
    free(g);g=NULL;
    free(b);b=NULL;
  
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
  
    QMu = QInit(0.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.));
    QExpMuFoisPiSur2 = QInitExp(QMu,M_PI/2.);
    
    //on filtre
    switch(intSens)
    {
      case 1 : filtre_Sang_H(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
               break;
      case 2 : filtre_Sang_V(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
               break;
      case 3 : filtre_Sang_Hbis(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
               break;
      case 4 : filtre_Sang_Vbis(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
               break;
    }
    
    QMatrixFree(dim[0],&QMatrix);
    //on verifie que la martie reelle est bien nulle.
    printf("la matrice filtre est quaternionique pure (1 oui,0 non):%d\n",Is_Pure_Image(QMatFiltre,dim));
    
    //à la sortie du filtre, certaines couleurs ne sont pas dans le range 0-255 donc il faut faire
    //un rescale pour eviter les "fausses" couleurs.
    QMatrixAllocate(dim[0],dim[1],&QMatScaled);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatScaled);
    QMat_Min_Max((const Quaternion **)QMatFiltre,dim,&dblRMin,&dblRMax,&dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
    QImMat_ChgScale_IntLvlMax((const Quaternion **)QMatFiltre,dim,255,dblIMin,dblIMax,
      dblJMin,dblJMax,dblKMin,dblKMax,&QMatScaled);
    QMatrixFree(dim[0],&QMatFiltre);
    
    //on enregistre l'image filtre
    r = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    g = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    b = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    QGetMatrixImagPart(&r,&g,&b,dim[0],dim[1],QMatScaled);
    QMatrixFree(dim[0],&QMatScaled);
    MgkEcrireImgCouleur(strFileOut,r,g,b,dim);
    printf("\nimage %s creee\n\n",strFileOut);
    free(r);r=NULL;
    free(g);g=NULL;
    free(b);b=NULL;
  }
  else
    printf("le sens est soit horizontal:1,3, soit vertical :2,4\n");
}

void filtrage_distance_couleur(char * strFile,int intSens)
{
  int dim[2],type;
  double *r=NULL,*g=NULL,*b=NULL;
  Quaternion **QMatrix=NULL,**QMatFiltre=NULL;
  char * result=NULL; //pour la chaine de caractere modifiant le nom du fichier
  Quaternion QMu,QExpMuFoisPiSur2;
  double **dblMatDist,**dblMatDistScaled; //les matrices distances couleur
  double * dblDist; //la meme mais version 1D
  double dblMin,dblMax;
  
  printf("Filtrage Distance Couleur\n");
  if ((intSens == 1)||(intSens == 2)||(intSens == 3)||(intSens == 4)) // bon parametre de sens
  {
    //on remplit la matrice quaternionique qui correspond a l'image du chemin passé en parametre
    //lecture de l'en-tete
    MgkTypeImage(strFile,dim,&type);
  
    //allocation des tableaux r,g et b
    r = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    g = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    b = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  
    //remplissage des tableaux
    MgkLireImgCouleur(strFile,r,g,b);
  
    //copie dans matrice quaternionique
    QMatrixAllocate(dim[0],dim[1],&QMatrix);
    QSetMatrix(r,g,b,dim[0],dim[1],&QMatrix);
    free(r);r=NULL;
    free(g);g=NULL;
    free(b);b=NULL;
  
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
  
    QMu = QInit(0.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.));
    QExpMuFoisPiSur2 = QInitExp(QMu,M_PI/2.);
    
    //on filtre
    switch(intSens)
    {
      case 1 : filtre_Sang_H(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
               break;
      case 2 : filtre_Sang_V(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
               break;
      case 3 : filtre_Sang_Hbis(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
               break;
      case 4 : filtre_Sang_Vbis(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
               break;
    }
    QMatrixFree(dim[0],&QMatrix);
    
    //on verifie que la martie reelle est bien nulle.
    printf("la matrice filtre est quaternionique pure (1 oui,0 non):%d\n",Is_Pure_Image(QMatFiltre,dim));
    
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDist);
    distance_couleur(QMatFiltre,&dblMatDist,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    //il faut faire un rescale
    Dbl_Matrix_Min_Max((const double **)dblMatDist,dim,&dblMin,&dblMax);
    
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistScaled);
    Dbl_Mat_ChgScale_IntLvlMax((const double **)dblMatDist,dim,255,dblMin,dblMax,&dblMatDistScaled);
    Matrix_Free_Double(dim[0],&dblMatDist);
    
    //il faut convertir la matrice 2D en matrice 1D
    dblDist = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    GetDblTabMatrix(dblMatDistScaled,&dblDist,dim);
    Matrix_Free_Double(dim[0],&dblMatDistScaled);

    
    result = (char *)malloc((10+strlen(strFile))*sizeof(char));
    strcpybtw(&result,strFile,"-dist",'.');
    MgkEcrireImgGris(result,dblDist,dim);
    printf("\nimage %s creee\n\n",result);
    free(result);result=NULL;
    free(dblDist);dblDist=NULL;
  }
  else
    printf("le sens est soit horizontal:1, soit vertical :2\n");
}

void filtrage_max_distances_couleur_L2(char * strFile)
{
  int dim[2],type;
  double *r=NULL,*g=NULL,*b=NULL;
  Quaternion **QMatrix=NULL,**QMatFiltre=NULL;
  char * result=NULL; //pour la chaine de caractere modifiant le nom du fichier
  Quaternion QMu,QExpMuFoisPiSur2;
  double **dblMatDistMaxL2,**dblMatDistMaxScaled; //les matrices distances couleur
  double **dblMatDistH=NULL,**dblMatDistV=NULL,**dblMatDistD1=NULL,**dblMatDistD2=NULL;
  double * dblDistMax1D; //la meme mais version 1D
  double dblMin=0,dblMax=0;
  
  printf("Filtrage Max Distances Couleur\n");
    //on remplit la matrice quaternionique correspond a l'image du chemin passé en parametre
    //lecture de l'en-tete
    MgkTypeImage(strFile,dim,&type);
  
    //allocation des tableaux r,g et b
    r = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    g = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    b = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  
    //remplissage des tableaux
    MgkLireImgCouleur(strFile,r,g,b);
  
    //copie dans matrice quaternionique
    QMatrixAllocate(dim[0],dim[1],&QMatrix);
    QSetMatrix(r,g,b,dim[0],dim[1],&QMatrix);
    free(r);r=NULL;
    free(g);g=NULL;
    free(b);b=NULL;
  
    QMu = QInit(0.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.));
    QExpMuFoisPiSur2 = QInitExp(QMu,M_PI/2.);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_H(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistH);
    Matrix_Init_Double(0.,dim,&dblMatDistH);
    distance_couleur(QMatFiltre,&dblMatDistH,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_V(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistV);
    Matrix_Init_Double(0.,dim,&dblMatDistV);
    distance_couleur(QMatFiltre,&dblMatDistV,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_D1(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistD1);
    Matrix_Init_Double(0.,dim,&dblMatDistD1);
    distance_couleur(QMatFiltre,&dblMatDistD1,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_D2(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistD2);
    Matrix_Init_Double(0.,dim,&dblMatDistD2);
    distance_couleur(QMatFiltre,&dblMatDistD2,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixFree(dim[0],&QMatrix);
    
    //on calcule la matrice de max de distance avec la norme L2
    //autrement dit la norme euclidienne
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistMaxL2);
    Matrix_Init_Double(0.,dim,&dblMatDistMaxL2);
    DblMatL2Norme(&dblMatDistMaxL2,dblMatDistH,dblMatDistV,dblMatDistD1,dblMatDistD2,dim);
    
    Matrix_Free_Double(dim[0],&dblMatDistH);
    Matrix_Free_Double(dim[0],&dblMatDistV);
    Matrix_Free_Double(dim[0],&dblMatDistD1);
    Matrix_Free_Double(dim[0],&dblMatDistD2);
    
    //ici changement d'echelle
    Dbl_Matrix_Min_Max((const double **)dblMatDistMaxL2,dim,&dblMin,&dblMax);
    printf("%lf %lf\n",dblMin,dblMax);
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistMaxScaled);
    Dbl_Mat_ChgScale_IntLvlMax((const double **)dblMatDistMaxL2,dim,255,dblMin,dblMax,&dblMatDistMaxScaled);
    Matrix_Free_Double(dim[0],&dblMatDistMaxL2);

    //il faut convertir la matrice 2D en matrice 1D
    dblDistMax1D = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    GetDblTabMatrix(dblMatDistMaxScaled,&dblDistMax1D,dim);
    Matrix_Free_Double(dim[0],&dblMatDistMaxScaled);
    
    result = (char *)malloc((20+strlen(strFile))*sizeof(char));
    strcpybtw(&result,strFile,"-distMaxL2",'.');
    MgkEcrireImgGris(result,dblDistMax1D,dim);
    printf("\nimage %s creee\n\n",result);
    free(result);result=NULL;
    free(dblDistMax1D);dblDistMax1D=NULL;
}

void filtrage_max_distances_couleur_Linf(char * strFile)
{
  int dim[2],type;
  double *r=NULL,*g=NULL,*b=NULL;
  Quaternion **QMatrix=NULL,**QMatFiltre=NULL;
  char * result=NULL; //pour la chaine de caractere modifiant le nom du fichier
  Quaternion QMu,QExpMuFoisPiSur2;
  double **dblMatDistMax,**dblMatDistMaxScaled; //les matrices distances couleur
  double **dblMatDistH=NULL,**dblMatDistV=NULL,**dblMatDistD1=NULL,**dblMatDistD2=NULL;
  double * dblDistMax1D; //la meme mais version 1D
  double dblMin=0,dblMax=0;
  
  printf("Filtrage Max Distances Couleur\n");
    //on remplit la matrice quaternionique correspond a l'image du chemin passé en parametre
    //lecture de l'en-tete
    MgkTypeImage(strFile,dim,&type);
  
    //allocation des tableaux r,g et b
    r = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    g = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    b = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  
    //remplissage des tableaux
    MgkLireImgCouleur(strFile,r,g,b);
  
    //copie dans matrice quaternionique
    QMatrixAllocate(dim[0],dim[1],&QMatrix);
    QSetMatrix(r,g,b,dim[0],dim[1],&QMatrix);
    free(r);r=NULL;
    free(g);g=NULL;
    free(b);b=NULL;
  
    QMu = QInit(0.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.));
    QExpMuFoisPiSur2 = QInitExp(QMu,M_PI/2.);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_H(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistH);
    Matrix_Init_Double(0.,dim,&dblMatDistH);
    distance_couleur(QMatFiltre,&dblMatDistH,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_V(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistV);
    Matrix_Init_Double(0.,dim,&dblMatDistV);
    distance_couleur(QMatFiltre,&dblMatDistV,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_D1(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistD1);
    Matrix_Init_Double(0.,dim,&dblMatDistD1);
    distance_couleur(QMatFiltre,&dblMatDistD1,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_D2(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistD2);
    Matrix_Init_Double(0.,dim,&dblMatDistD2);
    distance_couleur(QMatFiltre,&dblMatDistD2,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixFree(dim[0],&QMatrix);
    
    //on calcule la matrice de max de distance avec la norme Linfini
    //autrement dit la norme maximale
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistMax);
    Matrix_Init_Double(0.,dim,&dblMatDistMax);
    DblMatMaxed(&dblMatDistMax,dblMatDistH,dim);
    DblMatMaxed(&dblMatDistMax,dblMatDistV,dim);
    DblMatMaxed(&dblMatDistMax,dblMatDistD1,dim);
    DblMatMaxed(&dblMatDistMax,dblMatDistD2,dim);
    
    Matrix_Free_Double(dim[0],&dblMatDistH);
    Matrix_Free_Double(dim[0],&dblMatDistV);
    Matrix_Free_Double(dim[0],&dblMatDistD1);
    Matrix_Free_Double(dim[0],&dblMatDistD2);
    
    //ici changement d'echelle
    Dbl_Matrix_Min_Max((const double **)dblMatDistMax,dim,&dblMin,&dblMax);
    printf("%lf %lf\n",dblMin,dblMax);
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistMaxScaled);
    Dbl_Mat_ChgScale_IntLvlMax((const double **)dblMatDistMax,dim,255,dblMin,dblMax,&dblMatDistMaxScaled);
    //Matrix_Free_Double(dim[0],&dblMatDistMax);

    //il faut convertir la matrice 2D en matrice 1D
    dblDistMax1D = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    GetDblTabMatrix(dblMatDistMaxScaled,&dblDistMax1D,dim);
    Matrix_Free_Double(dim[0],&dblMatDistMaxScaled);
    
    result = (char *)malloc((20+strlen(strFile))*sizeof(char));
    strcpybtw(&result,strFile,"-distMaxLinf",'.');
    MgkEcrireImgGris(result,dblDistMax1D,dim);
    printf("\nimage %s creee\n\n",result);
    free(result);result=NULL;
    free(dblDistMax1D);dblDistMax1D=NULL;
}

void filtrage_logmax_distances_couleur(char * strFile)
{
  int dim[2],type;
  double *r=NULL,*g=NULL,*b=NULL;
  Quaternion **QMatrix=NULL,**QMatFiltre=NULL;
  char * result=NULL; //pour la chaine de caractere modifiant le nom du fichier
  Quaternion QMu,QExpMuFoisPiSur2;
  double **dblMatDistMax,**dblMatDistMaxScaled; //les matrices distances couleur
  double **dblMatDistH=NULL,**dblMatDistV=NULL,**dblMatDistD1=NULL,**dblMatDistD2=NULL;
  double **dblLogMat=NULL;
  double * dblDistMax1D; //la meme mais version 1D
  double dblMin=0,dblMax=0;
  
  printf("Filtrage Log Max Distances Couleur\n");
    //on remplit la matrice quaternionique correspond a l'image du chemin passé en parametre
    //lecture de l'en-tete
    MgkTypeImage(strFile,dim,&type);
  
    //allocation des tableaux r,g et b
    r = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    g = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    b = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  
    //remplissage des tableaux
    MgkLireImgCouleur(strFile,r,g,b);
  
    //copie dans matrice quaternionique
    QMatrixAllocate(dim[0],dim[1],&QMatrix);
    QSetMatrix(r,g,b,dim[0],dim[1],&QMatrix);
    free(r);r=NULL;
    free(g);g=NULL;
    free(b);b=NULL;
  
    QMu = QInit(0.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.));
    QExpMuFoisPiSur2 = QInitExp(QMu,M_PI/2.);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_H(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistH);
    Matrix_Init_Double(0.,dim,&dblMatDistH);
    distance_couleur(QMatFiltre,&dblMatDistH,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_V(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistV);
    Matrix_Init_Double(0.,dim,&dblMatDistV);
    distance_couleur(QMatFiltre,&dblMatDistV,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_D1(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistD1);
    Matrix_Init_Double(0.,dim,&dblMatDistD1);
    distance_couleur(QMatFiltre,&dblMatDistD1,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_D2(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistD2);
    Matrix_Init_Double(0.,dim,&dblMatDistD2);
    distance_couleur(QMatFiltre,&dblMatDistD2,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixFree(dim[0],&QMatrix);
    
    //on calcule la matrice de max de distance
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistMax);
    Matrix_Init_Double(0.,dim,&dblMatDistMax);
    DblMatMaxed(&dblMatDistMax,dblMatDistH,dim);
    DblMatMaxed(&dblMatDistMax,dblMatDistV,dim);
    DblMatMaxed(&dblMatDistMax,dblMatDistD1,dim);
    DblMatMaxed(&dblMatDistMax,dblMatDistD2,dim);
    
    Matrix_Free_Double(dim[0],&dblMatDistH);
    Matrix_Free_Double(dim[0],&dblMatDistV);
    Matrix_Free_Double(dim[0],&dblMatDistD1);
    Matrix_Free_Double(dim[0],&dblMatDistD2);

    //ici changement d'echelle logarithmique
    Dbl_Matrix_Min_Max((const double **)dblMatDistMax,dim,&dblMin,&dblMax);
    printf("%lf %lf\n",dblMin,dblMax);
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistMaxScaled);
    Dbl_Mat_ChgScaleLog_IntLvlMax((const double **)dblMatDistMax,dim,255,dblMin,dblMax,&dblMatDistMaxScaled);
    Matrix_Free_Double(dim[0],&dblMatDistMax);

    //il faut convertir la matrice 2D en matrice 1D
    dblDistMax1D = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    GetDblTabMatrix(dblMatDistMaxScaled,&dblDistMax1D,dim);
    Matrix_Free_Double(dim[0],&dblMatDistMaxScaled);
    
    result = (char *)malloc((10+strlen(strFile))*sizeof(char));
    strcpybtw(&result,strFile,"-distLogMax",'.');
    MgkEcrireImgGris(result,dblDistMax1D,dim);
    printf("\nimage %s creee\n\n",result);
    free(result);result=NULL;
    free(dblDistMax1D);dblDistMax1D=NULL;
}

void RGBtoHSVTabPat(double *r,double *g,double *b,double **ph,double **ps,double **pv,int dim)
{
  double *r2,*g2,*b2,*h2,*s2,*v2;
  int i;
  r2 = (double*)malloc(sizeof(double)*dim);
  g2 = (double*)malloc(sizeof(double)*dim);
  b2 = (double*)malloc(sizeof(double)*dim);
  h2 = (double*)malloc(sizeof(double)*dim);
  s2 = (double*)malloc(sizeof(double)*dim);
  v2 = (double*)malloc(sizeof(double)*dim);
  for(i=0;i<dim;i++)
  {
    RGB0255toRGB01(r[i],g[i],b[i],&r2[i],&g2[i],&b2[i]);
    RGBtoHSVPat(r2[i],g2[i],b2[i],&h2[i],&s2[i],&v2[i]);
    H360S1V1toHSV255(h2[i],s2[i],v2[i],&((*ph)[i]),&((*ps)[i]),&((*pv)[i]));
  }
  free(r2);
  free(g2);
  free(b2);
  free(h2);
  free(s2);
  free(v2);
}

void HSVtoRGBTabPat(double *h,double *s,double *v,double **pr,double **pg,double **pb,int dim)
{
  double *r2,*g2,*b2,*h2,*s2,*v2;
  int i;
  r2 = (double*)malloc(sizeof(double)*dim);
  g2 = (double*)malloc(sizeof(double)*dim);
  b2 = (double*)malloc(sizeof(double)*dim);
  h2 = (double*)malloc(sizeof(double)*dim);
  s2 = (double*)malloc(sizeof(double)*dim);
  v2 = (double*)malloc(sizeof(double)*dim);
  for(i=0;i<dim;i++)
  {
    HSV255toH360S1V1(h[i],s[i],v[i],&h2[i],&s2[i],&v2[i]);
    HSVtoRGBPat(h2[i],s2[i],v2[i],&r2[i],&g2[i],&b2[i]);
    RGB01toRGB0255(r2[i],g2[i],b2[i],&((*pr)[i]),&((*pg)[i]),&((*pb)[i]));
  }
  free(r2);
  free(g2);
  free(b2);
  free(h2);
  free(s2);
  free(v2);
}



void filtrage_max_distances_couleurPIMHAI(const char * strFile1,
    const char * strFile2,const char * strFile3,const char * strFileOut)
{
  int dim[2],type;
  double *c1=NULL,*c2=NULL,*c3=NULL;
  Quaternion **QMatrix=NULL,**QMatFiltre=NULL;
  char * result=NULL; //pour la chaine de caractere modifiant le nom du fichier
  Quaternion QMu,QExpMuFoisPiSur2;
  double **dblMatDistMax,**dblMatDistMaxScaled; //les matrices distances couleur
  double **dblMatDistH=NULL,**dblMatDistV=NULL,**dblMatDistD1=NULL,**dblMatDistD2=NULL;
  double * dblDistMax1D; //la meme mais version 1D
  double dblMin=0,dblMax=0;
  
  printf("Filtrage Max Distances Couleur PIMHAI\n");
    //on remplit la matrice quaternionique correspond a l'image du chemin passé en parametre
    //lecture de l'en-tete
    
    MgkTypeImage(strFile1,dim,&type);
    c1 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    MgkLireImgGris(strFile1,c1);

    MgkTypeImage(strFile2,dim,&type);
    c2 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    MgkLireImgGris(strFile2,c2);
    
    MgkTypeImage(strFile3,dim,&type);
    c3 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    MgkLireImgGris(strFile3,c3);


    //on ecrit une image couleur pour voir le resultat
    //MgkEcrireImgCouleur("D:/images/filtre/PIMHAI/meandre001band05-08-01.bmp",c1,c2,c3,dim);
    
    //comme les composantes c1,c2 et c3 sont stockées sous forme de tableaux
    //nous pouvons nous en servir pour allouer une matrice de quaternions
    //et faire un traitement identique a un traitement couleur. Chaque
    //composante sera représentée par une composante couleur.
    //copie dans matrice quaternionique
    QMatrixAllocate(dim[0],dim[1],&QMatrix);
    QSetMatrix(c1,c2,c3,dim[0],dim[1],&QMatrix);
    free(c1);c1=NULL;
    free(c2);c2=NULL;
    free(c3);c3=NULL;
  
    QMu = QInit(0.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.));
    QExpMuFoisPiSur2 = QInitExp(QMu,M_PI/2.);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_H(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistH);
    Matrix_Init_Double(0.,dim,&dblMatDistH);
    distance_couleur(QMatFiltre,&dblMatDistH,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_V(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistV);
    Matrix_Init_Double(0.,dim,&dblMatDistV);
    distance_couleur(QMatFiltre,&dblMatDistV,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_D1(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistD1);
    Matrix_Init_Double(0.,dim,&dblMatDistD1);
    distance_couleur(QMatFiltre,&dblMatDistD1,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_D2(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistD2);
    Matrix_Init_Double(0.,dim,&dblMatDistD2);
    distance_couleur(QMatFiltre,&dblMatDistD2,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixFree(dim[0],&QMatrix);
    
    //on calcule la matrice de max de distance
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistMax);
    Matrix_Init_Double(0.,dim,&dblMatDistMax);
    DblMatMaxed(&dblMatDistMax,dblMatDistH,dim);
    DblMatMaxed(&dblMatDistMax,dblMatDistV,dim);
    DblMatMaxed(&dblMatDistMax,dblMatDistD1,dim);
    DblMatMaxed(&dblMatDistMax,dblMatDistD2,dim);
    
    Matrix_Free_Double(dim[0],&dblMatDistH);
    Matrix_Free_Double(dim[0],&dblMatDistV);
    Matrix_Free_Double(dim[0],&dblMatDistD1);
    Matrix_Free_Double(dim[0],&dblMatDistD2);
    
    //ici changement d'echelle
    Dbl_Matrix_Min_Max((const double **)dblMatDistMax,dim,&dblMin,&dblMax);
    printf("%lf %lf\n",dblMin,dblMax);
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistMaxScaled);
    Dbl_Mat_ChgScale_IntLvlMax((const double **)dblMatDistMax,dim,255,dblMin,dblMax,&dblMatDistMaxScaled);
    //Matrix_Free_Double(dim[0],&dblMatDistMax);

    //il faut convertir la matrice 2D en matrice 1D
    dblDistMax1D = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    GetDblTabMatrix(dblMatDistMaxScaled,&dblDistMax1D,dim);
    Matrix_Free_Double(dim[0],&dblMatDistMaxScaled);
    
    result = (char *)malloc((10+strlen(strFileOut))*sizeof(char));
    strcpybtw(&result,strFileOut,"-distMax",'.');
    MgkEcrireImgGris(result,dblDistMax1D,dim);
    printf("\nimage %s creee\n\n",result);
    free(result);result=NULL;
    free(dblDistMax1D);dblDistMax1D=NULL;
}


void filtrage_logmax_distances_couleurPIMHAI(const char * strFile1,
    const char * strFile2,const char * strFile3,const char * strFileOut)
{
  int dim[2],type;
  double *c1=NULL,*c2=NULL,*c3=NULL;
  Quaternion **QMatrix=NULL,**QMatFiltre=NULL;
  char * result=NULL; //pour la chaine de caractere modifiant le nom du fichier
  Quaternion QMu,QExpMuFoisPiSur2;
  double **dblMatDistMax,**dblMatDistMaxScaled; //les matrices distances couleur
  double **dblMatDistH=NULL,**dblMatDistV=NULL,**dblMatDistD1=NULL,**dblMatDistD2=NULL;
  double **dblLogMat=NULL;
  double * dblDistMax1D; //la meme mais version 1D
  double dblMin=0,dblMax=0;
  
  printf("Filtrage PIMHAI Log Max Distances Couleur\n");
    //on remplit la matrice quaternionique correspond a l'image du chemin passé en parametre
    //lecture de l'en-tete
    
    MgkTypeImage(strFile1,dim,&type);
    c1 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    MgkLireImgGris(strFile1,c1);

    MgkTypeImage(strFile2,dim,&type);
    c2 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    MgkLireImgGris(strFile2,c2);
    
    MgkTypeImage(strFile3,dim,&type);
    c3 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    MgkLireImgGris(strFile3,c3);
  
    //comme les composantes c1,c2 et c3 sont stockées sous forme de tableaux
    //nous pouvons nous en servir pour allouer une matrice de quaternions
    //et faire un traitement identique a un traitement couleur. Chaque
    //composante sera représentée par une composante couleur.
    //copie dans matrice quaternionique
    QMatrixAllocate(dim[0],dim[1],&QMatrix);
    QSetMatrix(c1,c2,c3,dim[0],dim[1],&QMatrix);
    free(c1);c1=NULL;
    free(c2);c2=NULL;
    free(c3);c3=NULL;
  
    QMu = QInit(0.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.));
    QExpMuFoisPiSur2 = QInitExp(QMu,M_PI/2.);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_H(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistH);
    Matrix_Init_Double(0.,dim,&dblMatDistH);
    distance_couleur(QMatFiltre,&dblMatDistH,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_V(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistV);
    Matrix_Init_Double(0.,dim,&dblMatDistV);
    distance_couleur(QMatFiltre,&dblMatDistV,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_D1(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistD1);
    Matrix_Init_Double(0.,dim,&dblMatDistD1);
    distance_couleur(QMatFiltre,&dblMatDistD1,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixAllocate(dim[0],dim[1],&QMatFiltre);
    QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatFiltre);
    filtre_Sang_D2(QMatrix,&QMatFiltre,QExpMuFoisPiSur2,dim);
    //on calcule la distance couleur
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistD2);
    Matrix_Init_Double(0.,dim,&dblMatDistD2);
    distance_couleur(QMatFiltre,&dblMatDistD2,QMu,dim);
    QMatrixFree(dim[0],&QMatFiltre);
    
    QMatrixFree(dim[0],&QMatrix);
    
    //on calcule la matrice de max de distance
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistMax);
    Matrix_Init_Double(0.,dim,&dblMatDistMax);
    DblMatMaxed(&dblMatDistMax,dblMatDistH,dim);
    DblMatMaxed(&dblMatDistMax,dblMatDistV,dim);
    DblMatMaxed(&dblMatDistMax,dblMatDistD1,dim);
    DblMatMaxed(&dblMatDistMax,dblMatDistD2,dim);
    
    Matrix_Free_Double(dim[0],&dblMatDistH);
    Matrix_Free_Double(dim[0],&dblMatDistV);
    Matrix_Free_Double(dim[0],&dblMatDistD1);
    Matrix_Free_Double(dim[0],&dblMatDistD2);

    //ici changement d'echelle logarithmique
    Dbl_Matrix_Min_Max((const double **)dblMatDistMax,dim,&dblMin,&dblMax);
    printf("%lf %lf\n",dblMin,dblMax);
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatDistMaxScaled);
    Dbl_Mat_ChgScaleLog_IntLvlMax((const double **)dblMatDistMax,dim,255,dblMin,dblMax,&dblMatDistMaxScaled);
    Matrix_Free_Double(dim[0],&dblMatDistMax);

    //il faut convertir la matrice 2D en matrice 1D
    dblDistMax1D = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    GetDblTabMatrix(dblMatDistMaxScaled,&dblDistMax1D,dim);
    Matrix_Free_Double(dim[0],&dblMatDistMaxScaled);
    
    result = (char *)malloc((10+strlen(strFileOut))*sizeof(char));
    strcpybtw(&result,strFileOut,"-distLogMax",'.');
    MgkEcrireImgGris(result,dblDistMax1D,dim);
    printf("\nimage %s creee\n\n",result);
    free(result);result=NULL;
    free(dblDistMax1D);dblDistMax1D=NULL;
}

//on compare le filtrage de Sangwine suivant la même direction mais en changeant l'ordre de convolution
int main_filtrage_sangwine(int argc, char *argv[])
{
  int sens;
  char strFile[200];
  char * result=NULL; //pour la chaine de caractere modifiant le nom du fichier
  
  //sens = 1;//filtre horizontal 1
  //sens = 2;//filtre vertical 1
  //sens = 3;//filtre horizontal 2
  //sens = 4;//filtre vertical 2

  //on copie la chaine de caractere
  strcpy (strFile,(const char *)"D:/images/filtre/Sangwine/house.bmp");
  
  //on définie la chaine de caractere du premier fichier de sortie
  result = (char *)malloc((10+strlen(strFile))*sizeof(char));
  strcpybtw(&result,strFile,"-filtreV1",'.');
  filtrage_image(strFile,2,result); 
  free(result);result=NULL;

  //on définie la chaine de caractere du deuxieme fichier de sortie
  result = (char *)malloc((10+strlen(strFile))*sizeof(char));
  strcpybtw(&result,strFile,"-filtreV2",'.');
  filtrage_image(strFile,4,result); 
  free(result);result=NULL;

  return -1;
}

int main_filtrage_sang(int argc, char *argv[])
{
  int sens;
  char strFile[200];  

  char * result=NULL; //pour la chaine de caractere modifiant le nom du fichier
  
  //printf("Image d'origine: %s\n",argv[1]);

  //sens = 1;//filtre horizontal 1
  //sens = 2;//filtre vertical 1
  sens = 3;//filtre horizontal 2
  //sens = 4;//filtre vertical 2

  //on copie la chaine de caractere
  //strcpy (strFile,(const char *)"D:/images/filtre/PIMHAI/meandre001band05-08-01.bmp");
  
  //on définie la chaine de caractere du premier fichier de sortie
  //result = (char *)malloc((10+strlen(strFile))*sizeof(char));
  //strcpybtw(&result,strFile,"-filtreH",'.');

  //filtrage_image(strFile,sens,result);
  //filtrage_image(strFile,1,result); //sens horizontal
  //filtrage_image(strFile,3,result); //sens horizontal bis
  //free(result);result=NULL;


  //on effecttue le filtrage par distance couleur quaternionique
  // pour obtenir un gradient de saturation 
  filtrage_max_distances_couleur_Linf("D:/images/quaternions/filtrage/spatial/zebre2.jpg");
  filtrage_max_distances_couleur_L2("D:/images/quaternions/filtrage/spatial/zebre2.jpg");
  
  filtrage_max_distances_couleur_Linf("D:/images/quaternions/filtrage/spatial/zebre1.jpg");
  filtrage_max_distances_couleur_L2("D:/images/quaternions/filtrage/spatial/zebre1.jpg");

  filtrage_max_distances_couleur_Linf("D:/images/quaternions/filtrage/spatial/perroquet.jpg");
  filtrage_max_distances_couleur_L2("D:/images/quaternions/filtrage/spatial/perroquet.jpg");

  
  //filtrage_logmax_distances_couleur("D:/images/quaternions/filtrage/spatial/formes.bmp");
  
  //filtrage_max_distances_couleurHSV(argv[1]);
  //filtrage_logmax_distances_couleurHSV(argv[1]);
  
  //ici on cherche a effectuer un filtrage sur les images de PIMHAI
  //filtrage_max_distances_couleurPIMHAI("D:/images/filtre/PIMHAI/meandre004-band2.ppm",
  //  "D:/images/filtre/PIMHAI/meandre004-band8.ppm",
  //  "D:/images/filtre/PIMHAI/meandre004-band13.ppm",
  //  "D:/images/filtre/PIMHAI/PIMHAI.bmp");
  
  //filtrage_max_distances_couleurPIMHAI("D:/images/filtre/PIMHAI/multi/meandre-est-001.05.ppm",
  //  "D:/images/filtre/PIMHAI/multi/meandre-est-001.08.ppm",
  //  "D:/images/filtre/PIMHAI/multi/meandre-est-001.01.ppm",
  //  "D:/images/filtre/PIMHAI/PIMHAIdistmax.bmp");
    
  //filtrage_logmax_distances_couleurPIMHAI("D:/images/filtre/PIMHAI/meandre004-band2.ppm",
  //  "D:/images/filtre/PIMHAI/meandre004-band8.ppm",
  //  "D:/images/filtre/PIMHAI/meandre004-band13.ppm",
  //  "D:/images/filtre/PIMHAI/PIMHAI.bmp");
  
  //filtrage_logmax_distances_couleurPIMHAI("D:/images/filtre/PIMHAI/multi/meandre-est-001.05.ppm",
  //  "D:/images/filtre/PIMHAI/multi/meandre-est-001.08.ppm",
  //  "D:/images/filtre/PIMHAI/multi/meandre-est-001.01.ppm",
  //  "D:/images/filtre/PIMHAI/PIMHAIlogdistmax.bmp");
  return -1;
}





