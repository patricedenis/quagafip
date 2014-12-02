#include "fourier_quaternion.h"
  

/*void QFFT_InDouble(double *r,double *g,double *b,int * dim,Quaternion *** pqfft)
{
  int i,j,bx,by,ind;
  Quaternion *q,mu1,mu2,mu3;
  Complex *c1,*c2,*C1,*C2;
  
  //pour l'instant on va faire une fft avec mu1=mulum et mu2=(j-k)/racine(2)  
  mu1 = QInit(0.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)); //Mu luminance
  mu2 = QInit(0.,0.,1./sqrt(2.),-1./sqrt(2.)); //mu2 donné dans l'article de Sangwine et Ell
  mu3 = QMult(mu1,mu2);
  q = (Quaternion *)malloc(dim[0]*dim[1]*sizeof(Quaternion ));
  c1 = (Complex *)malloc(dim[0]*dim[1]*sizeof(Complex ));
  c2 = (Complex *)malloc(dim[0]*dim[1]*sizeof(Complex ));
  ind=0;
  //initialisation des tableaux
  // on cherche a obtenir les deux matrices de complexes pour effectuer ensuite
  // une fft sur chacune d'entre elles.
  for(bx=0;bx<dim[0];bx++)
    for(by=0;by<dim[1];by++)
    {
      q[ind]=QInit(0.,r[ind],g[ind],b[ind]);
      // ici on initialise les tableaux de complexes pour utiliser ensuite deux fft
      //c1[ind] = 0. + x'*mu1
      //c2[ind] = y' + z'*mu1
      // x',y' et z' sont obtenus en utilisant la relation de projection sur un axe
      c1[ind] = Complex_InitXY(0.,(-1./2.)*(QAdd(QMult(q[ind],mu1),QMult(mu1,q[ind]))).a);  
      c2[ind] = Complex_InitXY((-1./2.)*(QAdd(QMult(q[ind],mu2),QMult(mu2,q[ind]))).a,      
            (-1./2.)*(QAdd(QMult(q[ind],mu3),QMult(mu3,q[ind]))).a);
      ind++;
    }
  //on effectue maintenant les deux ffts
  //fft_2d(Complex **f,int dim[2],Complex **fhat,int sens)
  printf("fast fourier transform performing ...\n");
  fft_2d(c1,dim,C1,1);
  printf("done.\n");
  free(c1);
  printf("fast fourier transform performing ...\n");
  fft_2d(c2,dim,C2,1);
  free(c2);
  printf("done.\n");
  // on reconstruit ensuite le resultat dans la matrice de quaternions
  // F(v,u) = C1(v,u) + C2(v,u) mu2 = W(v,u) + X(v,u)mu1 + Y(v,u)mu2 + Z(v,u)mu3
  // = C1.r(v,u) + C1.i(v,u)mu1 + C2.r(v,u)mu2 + C2.i(v,u)mu3
  ind=0;
  for(bx=0;bx<dim[0];bx++)
    for(by=0;by<dim[1];by++)
    {
      (*pqfft)[bx][by] = QScalAdd(QScalMult(mu1,C1[ind].i),C1[ind].r);
      (*pqfft)[bx][by] = QAdd((*pqfft)[bx][by],QScalMult(mu2,C2[ind].r));
      (*pqfft)[bx][by] = QAdd((*pqfft)[bx][by],QScalMult(mu3,C2[ind].i));
      ind ++;
    }
  free(C1);
  free(C2);

}*/


void QFFT(Quaternion **QIn,int * dim,Quaternion *** pQOut,Quaternion mu1,Quaternion mu2,int intSens)
{
  int bx,by;
  Quaternion mu3;
  Complex **c1,**c2,**C1,**C2;
  
  mu3 = QMult(mu1,mu2);
  Complex_Matrix_Allocate(dim[0],dim[1],&c1);
  Complex_Matrix_Allocate(dim[0],dim[1],&c2);
  //on est à l'origine dans la base (1,i,j,k)  
  // initialisation des tableaux
  // on cherche a obtenir les deux matrices de complexes pour effectuer ensuite
  // une fft sur chacune d'entre elles.
  // on effectue donc un changement de base pour faire les calculs dans (1,mu1,mu2,mu1mu2).
  for(bx=0;bx<dim[0];bx++)
    for(by=0;by<dim[1];by++)
    {
      // ici on initialise les tableaux de complexes pour utiliser ensuite deux ffts
      //c1[ind] = w' + x'*mu1
      //c2[ind] = y' + z'*mu1
      // w', x',y' et z' sont obtenus en utilisant la relation de projection sur un axe
      c1[bx][by] = Complex_InitXY(QReal(QIn[bx][by]),
      (-1./2.)*(QAdd(QMult(QIn[bx][by],mu1),QMult(mu1,QIn[bx][by]))).a);  
      c2[bx][by] = Complex_InitXY((-1./2.)*(QAdd(QMult(QIn[bx][by],mu2),QMult(mu2,QIn[bx][by]))).a,      
            (-1./2.)*(QAdd(QMult(QIn[bx][by],mu3),QMult(mu3,QIn[bx][by]))).a);
    }
  Complex_Matrix_Allocate(dim[0],dim[1],&C1);
  Complex_Matrix_Allocate(dim[0],dim[1],&C2);

  //on effectue maintenant les deux ffts dans (1,mu1,mu2,mu1mu2).
  //fft_2d(Complex **f,int dim[2],Complex **fhat,int sens)
  printf("quaternionic fast fourier transform performing ...\n");
  fft_2d(c1,dim,C1,intSens);
  Complex_Matrix_Free(dim[0],&c1);
  fft_2d(c2,dim,C2,intSens);
  Complex_Matrix_Free(dim[0],&c2);  
  
  // on reconstruit ensuite le resultat dans la matrice de quaternions
  // F(v,u) = C1(v,u) + C2(v,u) mu2 = W(v,u) + X(v,u)mu1 + Y(v,u)mu2 + Z(v,u)mu3
  // = C1.r(v,u) + C1.i(v,u)mu1 + C2.r(v,u)mu2 + C2.i(v,u)mu3
  for(bx=0;bx<dim[0];bx++)
    for(by=0;by<dim[1];by++)
    {
      (*pQOut)[bx][by] = QScalAdd(QScalMult(mu1,C1[bx][by].i),C1[bx][by].r);
      (*pQOut)[bx][by] = QAdd((*pQOut)[bx][by],QScalMult(mu2,C2[bx][by].r));
      (*pQOut)[bx][by] = QAdd((*pQOut)[bx][by],QScalMult(mu3,C2[bx][by].i));
    }
  //on est donc revenu ici à la base d'origine en (1,i,j,k)
  Complex_Matrix_Free(dim[0],&C1);
  Complex_Matrix_Free(dim[0],&C2);
  printf("done.\n");

}


//////////////////////////////////////////////////////////////////////////////

//         PROCEDURES DE TEST PERMETTANT DE SAVOIR 
//         SI LES PROPRIETES DE SYMETRIE DU SPECTRE
//       QUATERNIONIQUE DES IMAGES COULEURS EST RESPECTE         

//////////////////////////////////////////////////////////////////////////////

//on passe une matrice de quaternion et on verifie sa partie reelle
//on considere ici que la matrice a deja ete shiftee et donc nous sommes dans le repere
// ou l'origine est au milieu de l'image
//on assume aussi que les dimaensions de l'image sont multiples de deux
//etant donne que lors des calculs il y a des erreurs, on acceptera quelques erreurs dans la symetrie
//ce nombre etant fixe par un seuil
int symetrie_partieR_OK(Quaternion ** QMat,int *dim,double dblPrecision,double dblSeuilTolerance)
{
  int bx,by;
  int mid0 = dim[0]/2;
  int nberreur = 0;
  double dblRapport;
  
  //
  //les proprietes de symetrie pour la partie reelle sont données ici:
  //Qr(0,0)=Qr(N/2,M/2)=0
  //et Qr(-S,-T)=-Qr(S,T) pout tout S,T dans la matrice
  //
  //elements isolés
  if(QMat[0][0].a != 0.0)
    nberreur++;
  if(QMat[0][mid0].a != 0.0)
    nberreur++;
  if(QMat[mid0][0].a != 0.0)
    nberreur++;
  if(QMat[mid0][mid0].a != 0.0)
    nberreur++;
  //premiere ligne et colonne + ligne milieu
  //on considere ici que l'image est carrée
  for(bx=1;bx<mid0;bx++)
  {
    if(QMat[0][bx].a + QMat[0][dim[1]-bx].a > dblPrecision)
      nberreur++;
    if(QMat[bx][0].a + QMat[dim[1]-bx][0].a > dblPrecision)
      nberreur++;
    if(QMat[mid0][bx].a + QMat[mid0][dim[1]-bx].a > dblPrecision)
      nberreur++;
  }
  //le reste de l'image
  for (bx=1;bx<dim[0];bx++)
      for (by=1;by<mid0;by++)
      {
          //on procede par ligne pour ne pas verifier en doublon
          if(QMat[bx][by].a + QMat[dim[0]-bx][dim[1]-by].a > dblPrecision)
            nberreur++;
      }
  dblRapport=((double)nberreur/(dim[0]*dim[1]));
  printf("precision %lf, nb de points %d\n, imprecis:%d, %% d'erreur:%0.2lf, seuil tolerance%lf\n",
      dblPrecision,dim[0]*dim[1],nberreur,dblRapport,dblSeuilTolerance);
  if(dblRapport <= dblSeuilTolerance)
  {
    return TRUE;
  }
  else return FALSE;
}

//on passe une matrice de quaternion et on verifie sa partie imagainaire I
//etant donne que lors des calculs il y a des erreurs, on acceptera quelques erreurs dans la symetrie
//ce nombre etant fixe par un seuil
//on assume aussi que les dimaensions de l'image sont multiples de deux
int symetrie_partieImI_OK(Quaternion ** QMat,int *dim,double dblPrecision,double dblSeuilTolerance)
{
  int bx,by;
  int mid0 = dim[0]/2;
  int nberreur = 0;
  double dblRapport;
  
  //
  //les proprietes de symetrie pour la partie reelle sont données ici:
  //et Qr(-S,-T)=Qi(S,T) pout tout S,T dans la matrice
  //
  //premiere ligne et colonne + ligne milieu
  //on considere ici que l'image est carrée
  for(bx=1;bx<mid0;bx++)
  {
    if(QMat[0][bx].b - QMat[0][dim[1]-bx].b > dblPrecision)
      nberreur++;
    if(QMat[bx][0].b - QMat[dim[1]-bx][0].b > dblPrecision)
      nberreur++;
    if(QMat[mid0][bx].b - QMat[mid0][dim[1]-bx].b > dblPrecision)
      nberreur++;
  }
  //le reste de l'image
  for (bx=1;bx<dim[0];bx++)
      for (by=1;by<mid0;by++)
      {
          //on procede par ligne pour ne pas verifier en doublon
          if(QMat[bx][by].b - QMat[dim[0]-bx][dim[1]-by].b > dblPrecision)
            nberreur++;
      }
  dblRapport=((double)nberreur/(dim[0]*dim[1]));
  printf("precision %lf, nb d'erreur:%d, nb de points %d\n, %% d'erreur:%0.2lf, seuil tolerance%lf\n",
      dblPrecision,dim[0]*dim[1],nberreur,dblRapport,dblSeuilTolerance);
  if(dblRapport <= dblSeuilTolerance)
  {
    return TRUE;
  }
  else return FALSE;
}


void QFFT_Verifie_Symetries(char * strFile, Quaternion QMu1,Quaternion QMu2)
{
    
    int i,j,dim[2],type,blnReal,blnSym;
    double *r,*g,*b, *dblData;
    Quaternion **QMatrix,**QMatFourier,**QMatrixShifted;
    
    double dblMin,dblMax,dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax;
    char * result; //pour la chaine de caractere modifiant le nom du fichier
    double dblPrecision,dblTolerance;
      
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
    free(r);
	free(g);
	free(b);
	
    // traitement dans l'espace de fourier
    QMatrixAllocate(dim[0],dim[1],&QMatFourier);
    
	printf("fast quaternionic fourier transform is processing ...\n");
    QFFT(QMatrix,dim,&QMatFourier,QMu1,QMu2,1);
    printf("done.\n");
    QMatrixFree(dim[0],&QMatrix);
       
    //QMatrixDisp(QMatFourier,dim[0],dim[1],1,stdout);
    //QMatrixDisp(QMatFourier,dim[0],dim[1],2,stdout);
    //QMatrixDisp(QMatFourier,dim[0],dim[1],5,stdout);
    //il faut aussi verifier les symetries dans le spectre frequentiel
    dblPrecision = 0.001; //deux nombres différents d'au moins ce seuil compteront comme erreur de calcul
    dblTolerance = 0.0; // 0% d'erreur autorisée pour dire que la symetrie est respectée
    blnSym=symetrie_partieR_OK(QMatFourier,dim,dblPrecision,dblTolerance);
    printf("la partie reelle presente les symetries desirees? %d (oui 1, non 0)\n",blnSym);
    blnSym=symetrie_partieImI_OK(QMatFourier,dim,dblPrecision,dblTolerance);
    printf("la partie imaginaire I presente les symetries desirees? %d (oui 1, non 0)\n",blnSym);
 
	QMatrixFree(dim[0],&QMatFourier);
}





// strFile pointe vers l'image à ouvrir
// QMu indique la direction de la transformée de Fourier Utilisée dans cette fonction
void QFFT_Construct_Vues_Frequentielles(char * strFile, Quaternion QMu1,Quaternion QMu2, double seuil)
{
    
    double **dblMatAngle,**dblMatLogModule,**dblMatLogModuleScaled,
        **dblMatAngleScaled,**dblMatAngleMod,**dblMatAngle02PI;
    int i,j,dim[2],type,blnReal,blnSym;
    unsigned char ** blnLogModulusMaskMat;
    double *r,*g,*b, *dblData;
    Quaternion **QMatrix,**QMatFourier,**QMatrixShifted,**QMatAxe,
        **QAngleColourMasked,**QMatAxeMasked,**QMatAxeMaskedScaled;
    
    double dblMin,dblMax,dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax;
    char * result; //pour la chaine de caractere modifiant le nom du fichier
  
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
    printf("Moyenne de l'image d'origine :R%lf G%lf B%lf\n",VectMean_Double(r,dim),VectMean_Double(g,dim),VectMean_Double(b,dim));
    free(r);
	  free(g);
	  free(b);
	
  // traitement dans l'espace de fourier
  QMatrixAllocate(dim[0],dim[1],&QMatFourier);
    
	printf("fast quaternionic fourier transform is processing ...\n");
  QFFT(QMatrix,dim,&QMatFourier,QMu1,QMu2,1);
  printf("done.\n");
  QMatrixFree(dim[0],&QMatrix);
        
  QMatrixAllocate(dim[0],dim[1],&QMatrixShifted);
	QMatrixShift(QMatFourier,&QMatrixShifted,dim[0],dim[1]);
	QMatrixFree(dim[0],&QMatFourier);

  //maintenant on a la matrice QMatFourier qui comprend les informations fréquentielles de l'image

  
  // soucis avec cette méthode, on utilise des quaternions purs dans le domaine statial
  // faire attention qu'il n'y ait pas de partie réelle nulle dans le domaine fréquentiel
  // car problème après avec S(q)=0 pour calculer le phi=atan(|V(q)|/S(q))
  
  // donc ici on vérifie pas de partie réelle nulle
  // à faire
  blnReal=IsRealPresent(QMatrixShifted,dim);
  printf("La partie reelle de la transformee de Fourier Quaternionique est non nulle (1 oui, 0 non):%d\n",blnReal);
  
  
  
  
  // on alloue les différentes matrices qui vont nous servir juste aprrs.
  Matrix_Allocate_Double(dim[0],dim[1],&dblMatLogModule);
  Matrix_Allocate_Double(dim[0],dim[1],&dblMatAngle);
  QMatrixAllocate(dim[0],dim[1],&QMatAxe);
  
  //puis on calcule pour chaque point de la matrice, son module, son angle et son axe.
  Calcul_Matrice_LogRo_Mu_Phi((const Quaternion **)QMatrixShifted,&dblMatLogModule,&dblMatAngle,&QMatAxe,dim[0],dim[1]);
  //Calcul_Matrice_Ro_Mu_Phi((const Quaternion **)QMatrixShifted,&dblMatLogModule,&dblMatAngle,&QMatAxe,dim[0],dim[1]);
  QMatrixFree(dim[0],&QMatrixShifted);


  ///////////////////////////////////////////////////////////////////////////////////////////////
  //  MODULE
  printf("construction de la vue du module...\n");
  //le log du module est calculé suivant cette formule: *LogRo = log(QNorm(q));
  
  Dbl_Matrix_Min_Max((const double **) dblMatLogModule,dim,&dblMin,&dblMax);
  Matrix_Allocate_Double(dim[0],dim[1],&dblMatLogModuleScaled);
  
  Dbl_Mat_ChgScale_IntLvlMax((const double **) dblMatLogModule,dim,255,dblMin,dblMax,&dblMatLogModuleScaled);
  Matrix_Free_Double(dim[0],&dblMatLogModule);

  //on construit le masque a partir du log du module pour appliquer sur l'angle et l'axe
  Matrix_Allocate(dim[0],dim[1],&blnLogModulusMaskMat);
  Mask_From_LogModulus(seuil,dim,dblMatLogModuleScaled,&blnLogModulusMaskMat);
  
  //allocation sous forme de tableau monodimensionnel pour passe en argument a la fonction MgkEcrireImgGris
  dblData = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  
  
  //on recupere donc les données matricielles sous forme de tableau 1D
  GetDblTabMatrix(dblMatLogModuleScaled,&dblData,dim);
  Matrix_Free_Double(dim[0],&dblMatLogModuleScaled);
  
  
  result = (char *)malloc((10+strlen(strFile))*sizeof(char));
  strcpybtw(&result,strFile,"-qfft-mod",'.');
  MgkEcrireImgGris(result,dblData,dim);
  printf("\nimage %s creee\n\n",result);
  free(result);
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // PHASE

  printf("construction de la vue de l'angle...\n");
  //la phase est calculée suivant cette formule:  *Phi = atan(QNorm(QImag(q))/QReal(q));
  // -pi/2 < phi < pi/2
  
  Dbl_Matrix_Min_Max((const double **) dblMatAngle,dim,&dblMin,&dblMax);
  printf("Angle min %lf Angle max %lf\n",dblMin,dblMax);
  
  //Matrix_Allocate_Double(dim[0],dim[1],&dblMatAngle02PI);
  //on replace les valeurs de la matrice d'angle entre 0 et 2PI
  //ici peu importe le dim min et max on les considere egaux a -pi/2 et pi/2
  
  //Dbl_Mat_ChgScale_LvlMin_LvlMax((const double **) dblMatAngle,dim,0.,2*M_PI,-M_PI/2.,M_PI/2.,&dblMatAngle02PI);
  
  //Matrix_Free_Double(dim[0],&dblMatAngle);
  
  QMatrixAllocate(dim[0],dim[1],&QAngleColourMasked);
  //QInitMatrix(dim[0],dim[1], QInit(0.,0.,0.,0.),&QAngleColourMasked);
  // avant de pouvoir construire la matrice d'angle entre 0 et 2PI rad ou 0 et 360 degre
  //QMat_Angle_Construct_WithLogModMask(dim,blnLogModulusMaskMat,dblMatAngle02PI,&QAngleColourMasked);
  QMat_Angle_Construct_WithLogModMask(dim,blnLogModulusMaskMat,dblMatAngle,&QAngleColourMasked);
  Matrix_Free_Double(dim[0],&dblMatAngle);
  
  r = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  g = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  b = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  QGetMatrixImagPart(&r,&g,&b,dim[0],dim[1],QAngleColourMasked);
  QMatrixFree(dim[0],&QAngleColourMasked);
  
  result = (char *)malloc((10+strlen(strFile))*sizeof(char));
  strcpybtw(&result,strFile,"-qfft-angle",'.');
  MgkEcrireImgCouleur(result,r,g,b,dim);
  printf("\nimage %s creee\n\n",result);
  free(result);
      
  printf("Moyenne de l'image de l'angle :R%lf G%lf B%lf\n",VectMean_Double(r,dim),VectMean_Double(g,dim),VectMean_Double(b,dim));
  free(r);
	free(g);
	free(b);
   
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // AXE

  printf("construction de la vue de l'axe...\n");
  
  //l'axe est calculé suivant cette formule: *qMu = QScalDiv(QImag(q),QNorm(QImag(q)));
  //l'axe est donc un quaternion pur    

  QMat_Min_Max2((const Quaternion **)QMatAxe,dim,&dblMin,&dblMax);
  
  /*QMatrixAllocate(dim[0],dim[1],&QMatAxeScaled2);
  //bug
  QMat_ChgScale_IntLvlMax((const Quaternion **)QMatAxe,dim,(dblMax-dblMin),dblRMin,dblRMax,
      dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax,&QMatAxeScaled2);
  QMatrixFree(dim[0],&QMatAxe);
  
  QMatrixAllocate(dim[0],dim[1],&QPureAxeMatMasked);
  QMat_Axe_Construct_WithLogModMask(dim,blnLogModulusMaskMat,QMatAxeScaled2,&QPureAxeMatMasked);
  //QPureMat_Mask_Product(dim,blnModulusMaskMat,QMatAxe,&QPureAddMat);
  Matrix_Free(dim[0],&blnLogModulusMaskMat);
  //QMatrixFree(dim[0],&QMatAxe);
  QMatrixFree(dim[0],&QMatAxeScaled2);
  
  QMat_Min_Max((const Quaternion **)QPureAxeMatMasked,dim,&dblRMin,&dblRMax,&dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
  //printf("minI %lf maxI %lf minJ %lf maxJ %lf minK %lf maxK %lf\n",dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax);
  QMatrixAllocate(dim[0],dim[1],&QMatAxeScaled);
  QMat_ChgScale_IntLvlMax((const Quaternion **)QPureAxeMatMasked,dim,255,dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,
  dblKMin,dblKMax,&QMatAxeScaled);
  QMatrixFree(dim[0],&QPureAxeMatMasked);

  r = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  g = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  b = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  QGetMatrixImagPart(&r,&g,&b,dim[0],dim[1],QMatAxeScaled);
  QMatrixFree(dim[0],&QMatAxeScaled);*/

  QMatrixAllocate(dim[0],dim[1],&QMatAxeMasked);
  QMat_Axe_Construct_WithLogModMask(dim,blnLogModulusMaskMat,QMatAxe,&QMatAxeMasked);
  Matrix_Free(dim[0],&blnLogModulusMaskMat);
  QMatrixFree(dim[0],&QMatAxe);
  
  QMat_Min_Max((const Quaternion **)QMatAxeMasked,dim,&dblRMin,&dblRMax,&dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);

  QMatrixAllocate(dim[0],dim[1],&QMatAxeMaskedScaled);
  QMat_ChgScale_IntLvlMax((const Quaternion **)QMatAxeMasked,dim,255,dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,
  dblKMin,dblKMax,&QMatAxeMaskedScaled);
  QMatrixFree(dim[0],&QMatAxeMasked);

  r = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  g = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  b = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  QGetMatrixImagPart(&r,&g,&b,dim[0],dim[1],QMatAxeMaskedScaled);
  QMatrixFree(dim[0],&QMatAxeMaskedScaled);
  
  result = malloc((10+strlen(strFile))*sizeof(char));
  strcpybtw(&result,strFile,"-qfft-axe",'.');
  MgkEcrireImgCouleur(result,r,g,b,dim);
  printf("\nimage %s creee\n\n",result);
  free(result);
  free(r);r=NULL;
  free(g);g=NULL;
  free(b);b=NULL;
}


int QFFT_Init_Freq_RGB()
{
  Quaternion Qmu1,Qmu2;
  int intNbStripe;
  double InvRacine3;
  double K=2000.;
  int dim[2],i;
  double *dblCmp1,*dblCmp2,*dblCmp3, *rs,*gs,*bs;
  double dblMinR, dblMaxR, dblMinG, dblMaxG, dblMinB, dblMaxB, dblMin, dblMax;
 
  Quaternion **QMatrix,**QMatrixShifted, **QMatOUT;
    
  /*64 64 2 0.3333333 0.3333333 0.3333333*/
    
  intNbStripe = 2;
  dim[0]=256;dim[1]=256;
  Qmu1 = QInit(0.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)); //Mu luminance
  Qmu2 = QInit(0.,0.,1./sqrt(2.),-1./sqrt(2.)); //mu2 donné dans l'article de Sangwine et Ell
 
  //Qmu1 = QInit(0.,1./sqrt(2.),1./sqrt(2.),0.); //Mu jaune
  //Qmu2 = QInit(0.,0.,0.,1.); 
  
  //Qmu1 = QInit(0.,0.,0.,1.); //mu bleu
  //Qmu2 = QInit(0.,1.,0.,0.);
  
  //Qmu1 = QInit(0.,1./sqrt(2.),0.,1./sqrt(2.)); //Mu magenta
  //Qmu2 = QInit(0.,0.,1.,0.);
  
  // Allocation for a Quaterniotic matrix
   QMatrixAllocate(dim[0],dim[1],&QMatrix);
  //Initialisation of this matrix with zeros
  QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatrix);
	            
  //j'initialise sur la partie réelle
  //QMatrix[dim[0]/2 + intNbStripe][dim[1]/2 - intNbStripe].a =K; 
  //QMatrix[dim[0]/2 - intNbStripe][dim[1]/2 + intNbStripe].a =-K;
	
  //QMatrix[dim[0]/2 + intNbStripe][dim[1]/2 + intNbStripe].a =K; 
  //QMatrix[dim[0]/2 - intNbStripe][dim[1]/2 - intNbStripe].a =-K;
	
  //initialisation sur la prmeiere partie imaginaire		
  QMatrix[dim[0]/2 + intNbStripe][dim[1]/2 + intNbStripe].b =K; 
  QMatrix[dim[0]/2 - intNbStripe][dim[1]/2 - intNbStripe].b =K;
 		
	//QMatrix[dim[0]/2 + intNbStripe][dim[1]/2 + intNbStripe].c =K; 
 	//QMatrix[dim[0]/2 - intNbStripe][dim[1]/2 - intNbStripe].c =K;
 	
 	//QMatrix[dim[0]/2 + intNbStripe][dim[1]/2 + intNbStripe].d =K; 
 	//QMatrix[dim[0]/2 - intNbStripe][dim[1]/2 - intNbStripe].d =K;
 		
  //allocation de la matrice shift
  QMatrixAllocate(dim[0],dim[1],&QMatrixShifted);
  //we need to shift the matrix before performing the inverse Fourier transform
			
  QMatrixShift(QMatrix,&QMatrixShifted,dim[0],dim[1]);
  QMatrixFree(dim[0],&QMatrix);
			
  //initialisation of the matrix that will receive the insverse Fourier image
  QMatrixAllocate(dim[0],dim[1],&QMatOUT);
  // perform the inverse fourier transform
  //on initialise la matrice de sortie avec des zeros
  QInitMatrix(dim[0],dim[1],QInit(0.,0.,0.,0.),&QMatOUT);
  QFFT(QMatrixShifted,dim,&QMatOUT,Qmu1,Qmu2,1);
  QMatrixFree(dim[0],&QMatrixShifted);
  
  //on alloue les tableaux qui serviront pour l'ecriture sur le fichier  
  dblCmp1 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  dblCmp2 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  dblCmp3 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  
  //on recupere les données couleurs de la matrice inverse (donc spatiale) dans les tableaux  
  QGetMatrixImagPart(&dblCmp1,&dblCmp2,&dblCmp3,dim[0],dim[1],QMatOUT);

  //for(i=0;i<dim[0]*dim[1];i++)
  //{
    //dblCmp1[i] = dblCmp1[i]*10+128;
    //dblCmp2[i] = dblCmp2[i]*10+128;
    //dblCmp3[i] = dblCmp3[i]*8+128;
  //}

  QMatrixFree(dim[0],&QMatOUT);
  
  //il faut maintenant effectuer un rescale de la matrice statiale pour obtenir 
  //un résultat affichable correctement
  Dbl_Tab_Min_Max((const double *)dblCmp1,dim[0]*dim[1],&dblMinR,&dblMaxR);
  Dbl_Tab_Min_Max((const double *)dblCmp2,dim[0]*dim[1],&dblMinG,&dblMaxG);
  Dbl_Tab_Min_Max((const double *)dblCmp3,dim[0]*dim[1],&dblMinB,&dblMaxB);
  dblMin = Denis_MIN3(dblMinR,dblMinG,dblMinB);
  dblMax = Denis_MAX3(dblMaxR,dblMaxG,dblMaxB);
  rs = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  gs = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  bs = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  
  //changement d'échelle
  Dbl_Tab_ChgScale_IntLvlMax((const double *)dblCmp1,dim[0]*dim[1],255,dblMin,dblMax,&rs);
  Dbl_Tab_ChgScale_IntLvlMax((const double *)dblCmp2,dim[0]*dim[1],255,dblMin,dblMax,&gs);
  Dbl_Tab_ChgScale_IntLvlMax((const double *)dblCmp3,dim[0]*dim[1],255,dblMin,dblMax,&bs);

  
  //on ecrit le fichier de sortie
  MgkEcrireImgCouleur("C:/images/variations/test.bmp",rs,gs,bs,dim);
  free(dblCmp1);dblCmp1=NULL;
  free(dblCmp2);dblCmp2=NULL;
  free(dblCmp3);dblCmp3=NULL;
  free(rs);rs=NULL;
  free(gs);gs=NULL;
  free(bs);bs=NULL;
}





int main_qfft(int argc, char *argv[])
{
    double InvRacine3;
    double seuil;//seuil définissant un masque de module
    int dim[2];//dimension de l'image
    Quaternion Qmu1,Qmu2;//Quaternion pur utilisé comme direction de l'analyse de Fourier

    /*
    //Initialisation des paramètres par défaut
    Qmu1 = QInit(0.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)); //Mu luminance
    Qmu2 = QInit(0.,0.,1./sqrt(2.),-1./sqrt(2.)); //mu2 donné dans l'article de Sangwine et Ell
    
    //Qmu1 = QInit(0.,0.,1.,0.); //Mu j
    //Qmu1 = QInit(0.,1.,0.,0.); //Mu j
    //Qmu2 = QInit(0.,0.,0.,1.); //mu k
    
    //le seuil pour le masque de module
    seuil = 0.23;
    //D:/images/fft/couleur/lenna.png
    //"D:/images/fft/couleur/disques.ppm"
    printf("image 3 debug qfft, fichier :%s\t + chaine en argument : %s\n",argv[1],argv[2]);
    QFFT_Construct_Vues_Frequentielles(argv[1],Qmu1,Qmu2,seuil);
    QFFT_Verifie_Symetries(argv[1],Qmu1,Qmu2);
    */
    
    QFFT_Init_Freq_RGB();
}



