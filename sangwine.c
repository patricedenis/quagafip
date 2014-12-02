#include "sangwine.h"

//la procedure suivante effectue la decomposition sympletique sur la matrice
//de quaternions passée en exemple. la base de départ est donc
// (1 i j k) et la base d'arrivée est (1 mu1 mu2 mu1mu2) avec mu1 et mu2 passés en parametre.
// on obtient en résultat les trois matrices de double associées à chaque vecteur de la nouvelle base
// Sachant que la composante réelle de l'image d'origine ne change pas son equivalent dans la nouvelle base
// est obtenu directement.
// L'image d'origine peut se décomposer comme suit:
// f = a + bi + cj + dk
// la decomposition sympletique permet de separer l'image d'origine en une partie parallele (f'1 = simplexe)
// à l'axe mu1 et une partie perpendiculaire (f'2 = perplexe)
// f' = f'1 + f'2 = (a' + b'mu1) + (c' + d'mu1)mu2
void qtnSympleticDecomposition(const Quaternion ** qtnMatrixInQtnBasis,
  double *** pdblMatrixProjR,  double *** pdblMatrixProjMu1,
  double *** pdblMatrixProjMu2, double *** pdblMatrixProjMu1Mu2,
  int * intDim, Quaternion qtnMu1, Quaternion qtnMu2)
{
  int intCpt1,intCpt2;  //les compteurs de boucles
  Quaternion qtnMu1Mu2; //dernier elt de la base
  
  //initialisation du troisième element de la base d'arrivée
  qtnMu1Mu2 = QMult(qtnMu1,qtnMu2);
  
  // on effectue donc un changement de base pour faire les calculs dans (1,mu1,mu2,mu1mu2).
  for(intCpt1=0;intCpt1<intDim[0];intCpt1++)
    for(intCpt2=0;intCpt2<intDim[1];intCpt2++)
    {
      //on calcule chaque coefficient de la nouvelle base qu'on place dans les matrices de double
      //la premiere partie est la partie reelle (a')
      (*pdblMatrixProjR)[intCpt1][intCpt2] = QReal(qtnMatrixInQtnBasis[intCpt1][intCpt2]);
      //la deuxieme partie est la projection sur mu1 (b')
      (*pdblMatrixProjMu1)[intCpt1][intCpt2] = (-1./2.)*QReal(QAdd(QMult(qtnMatrixInQtnBasis[intCpt1][intCpt2],qtnMu1),
            QMult(qtnMu1,qtnMatrixInQtnBasis[intCpt1][intCpt2]))); 
      //la troisieme partie est la projection sur mu2 (c')
      (*pdblMatrixProjMu2)[intCpt1][intCpt2] = (-1./2.)*QReal(QAdd(QMult(qtnMatrixInQtnBasis[intCpt1][intCpt2],qtnMu2),
            QMult(qtnMu2,qtnMatrixInQtnBasis[intCpt1][intCpt2])));
       //la derniere partie est la projection sur mu1mu2 (d')
      (*pdblMatrixProjMu1Mu2)[intCpt1][intCpt2] = (-1./2.)*QReal(QAdd(QMult(qtnMatrixInQtnBasis[intCpt1][intCpt2],qtnMu1Mu2),
            QMult(qtnMu1Mu2,qtnMatrixInQtnBasis[intCpt1][intCpt2])) );
    }
}


//cette procedure effectue le changement de base entre la base actuelle et celle
//specifier par les vecteurs de la base (E1,E2,E3)
//il est a notre cependant qu'aucune modification n'est apportée a la composante
//réelle des quaternions de la matrice d'origine. En effet on suppose qu'ils sont nuls.
void qtnChangementDeBaseSurLesImag(const Quaternion ** qtnMatrixInQtnBasis,
  double *** pdblMatProjE1,double *** pdblMatProjE2, double *** pdblMatProjE3,
  int * intDim, Quaternion qtnE1, Quaternion qtnE2, Quaternion qtnE3)
{
  int intCpt1,intCpt2;  //les compteurs de boucles
  
  // on effectue donc un changement de base pour faire les calculs dans (1,E1,E2,E3).
  for(intCpt1=0;intCpt1<intDim[0];intCpt1++)
    for(intCpt2=0;intCpt2<intDim[1];intCpt2++)
    {
      //on calcule chaque coefficient de la nouvelle base qu'on place dans les matrices de double
      //la deuxieme partie est la projection sur E1 (b')
      (*pdblMatProjE1)[intCpt1][intCpt2] = (-1./2.)*QReal(QAdd(QMult(qtnMatrixInQtnBasis[intCpt1][intCpt2],qtnE1),
            QMult(qtnE1,qtnMatrixInQtnBasis[intCpt1][intCpt2]))); 
      //la troisieme partie est la projection sur mu2 (c')
      (*pdblMatProjE2)[intCpt1][intCpt2] = (-1./2.)*QReal(QAdd(QMult(qtnMatrixInQtnBasis[intCpt1][intCpt2],qtnE2),
            QMult(qtnE2,qtnMatrixInQtnBasis[intCpt1][intCpt2])));
       //la derniere partie est la projection sur mu1mu2 (d')
      (*pdblMatProjE3)[intCpt1][intCpt2] = (-1./2.)*QReal(QAdd(QMult(qtnMatrixInQtnBasis[intCpt1][intCpt2],qtnE3),
            QMult(qtnE3,qtnMatrixInQtnBasis[intCpt1][intCpt2])) );
    }
}

// la procedure qui suit permet de construire la matrice
// de la partie parallele résultant de la décomposition sympletique
// à partir des coefficients b' contenus dans le matrice
// de doubles passée en parametre. Sachant que le coefficient 
// a' ne sert a rien nous ne le traiterons pas ici
// elle est calculée comme suit :
// f'1 = b'Mu1
void qtnConstructParaMat(const double ** dblMatMu1, Quaternion *** pqtnMatPara, int * intDim, Quaternion qtnMu1)
{
  int intCpt1,intCpt2;  //les compteurs de boucles

  //on obtient le résultat en multipliant la matrice de double par mu1
  for(intCpt1=0;intCpt1<intDim[0];intCpt1++)
    for(intCpt2=0;intCpt2<intDim[1];intCpt2++)
    {
      (*pqtnMatPara)[intCpt1][intCpt2] = QScalMult((const Quaternion)qtnMu1,
        dblMatMu1[intCpt1][intCpt2]);
    }     
}


void qtnConstructPerpMat(const double ** dblMatMu2, const double ** dblMatMu3, Quaternion *** pqtnMatPerp,
  int * intDim, Quaternion qtnMu1, Quaternion qtnMu2)
{
  int intCpt1,intCpt2;  //les compteurs de boucles
  Quaternion qtnMu3; //dernier elt de la base
  
  //initialisation du troisième element de la base d'arrivée
  qtnMu3 = QMult(qtnMu1,qtnMu2);
  
  //on obtient la partie perplexe en multipliant la matrice de double par mu2 et mu3
  for(intCpt1=0;intCpt1<intDim[0];intCpt1++)
    for(intCpt2=0;intCpt2<intDim[1];intCpt2++)
    {
      (*pqtnMatPerp)[intCpt1][intCpt2] = QAdd( QInit(0.,dblMatMu2[intCpt1][intCpt2],0.,0.),
        QScalMult((const Quaternion)qtnMu2,dblMatMu3[intCpt1][intCpt2]) );
    }     
}

void qtnRecomposition(const double ** dblMatMu1,const double ** dblMatMu2,const double ** dblMatMu3,
  int * intDim, Quaternion qtnMu1, Quaternion qtnMu2,Quaternion *** pQMatrix)
{
  Quaternion qtnMu3;
  int intCpt1,intCpt2;  //les compteurs de boucles
  
  qtnMu3 = QMult(qtnMu1,qtnMu2);
  for(intCpt1=0;intCpt1<intDim[0];intCpt1++)
    for(intCpt2=0;intCpt2<intDim[1];intCpt2++)
      (*pQMatrix)[intCpt1][intCpt2] = QAdd(QAdd(      
        QScalMult((const Quaternion)qtnMu1,dblMatMu1[intCpt1][intCpt2]), //b'Mu1
        QScalMult((const Quaternion)qtnMu2,dblMatMu2[intCpt1][intCpt2])), //c'Mu2 
        QScalMult((const Quaternion)qtnMu3,dblMatMu3[intCpt1][intCpt2])); //d'Mu3 
}

// strFile pointe vers l'image à ouvrir
// on va décomposer l'image RGB de maniere simpletique en utilisant la procédure SympleticDecomposition
void qtnConstructSympleticDecomposition(char * strFile)
{
    
  Quaternion qtnI,qtnJ,qtnK; //les element de la base d'origine
  Quaternion qtnMu1, qtnMu2, qtnMu1Mu2; //les elements de la nouvelle base
    
  double *r,*g,*b; //les tableaux de double pour ouvrir les images couleur
  double *dblData1,*dblData2,*dblData3; //idem pour image en NdG
  int dim[2]; //pour les dimensions de l'image
  int type; //le type de l'image
  //les différentes matrices de quaternions nécessaires.
  Quaternion **QMatrix,**QMatrixPara, **QMatrixPerp,
    **QMatScaled,**QMatrixParaIJK,**QMatrixPerpIJK,**QMatrixSomme;
  //les différentes matrices de double résultant de la décomposition symplétique
  double **dblMatR, **dblMatE1, **dblMatE2, **dblMatE3,**dblMatScal, **dblMatMu1, **dblMatMu2, **dblMatMu3;
  //les differents extrema pour les rescale
  double dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax,dblMin,dblMax;
  char * result; //pour la chaine de caractères modifiant le nom du fichier
  
  //initialisation de i,j et K
  qtnI = QInit(0.,1.,0.,0.); //i
  qtnJ = QInit(0.,0.,1.,0.); //j
  qtnK = QInit(0.,0.,0.,1.); //k  
  //initialisation de mu1 et mu2
  qtnMu1 = QInit(0.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)); //Mu luminance ou Mugris
  //qtnMu2 = QInit(0.,0.,1./sqrt(2.),-1./sqrt(2.)); //mu2 donné dans l'article de Sangwine et Ell
  qtnMu2 = QInit(0.,sqrt(2./3.),-sqrt(2./3.)/2.,-sqrt(2./3.)/2.); //mu2 rejection du rouge (i) par rapport à MuGris
  qtnMu1Mu2 = QMult(qtnMu1,qtnMu2);
  
        
  //OUVERTURE DE L'IMAGE
  //lecture de l'en-tete
  MgkTypeImage(strFile,dim,&type);

  //allocations
  //allocation des tableaux r,g et b
  r = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  g = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  b = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  QMatrixAllocate(dim[0],dim[1],&QMatrix);  
  Matrix_Allocate_Double(dim[0],dim[1],&dblMatR);
  Matrix_Allocate_Double(dim[0],dim[1],&dblMatMu1);
  Matrix_Allocate_Double(dim[0],dim[1],&dblMatMu2);
  Matrix_Allocate_Double(dim[0],dim[1],&dblMatMu3);
  Matrix_Allocate_Double(dim[0],dim[1],&dblMatE1);
  Matrix_Allocate_Double(dim[0],dim[1],&dblMatE2);
  Matrix_Allocate_Double(dim[0],dim[1],&dblMatE3);
  result = malloc((20+strlen(strFile))*sizeof(char));
  QMatrixAllocate(dim[0],dim[1],&QMatrixPara);
  QMatrixAllocate(dim[0],dim[1],&QMatrixPerp);    
  QMatrixAllocate(dim[0],dim[1],&QMatrixParaIJK);
  QMatrixAllocate(dim[0],dim[1],&QMatrixPerpIJK);    
  QMatrixAllocate(dim[0],dim[1],&QMatScaled);
  QMatrixAllocate(dim[0],dim[1],&QMatrixSomme);
  //remplissage des tableaux
  MgkLireImgCouleur(strFile,r,g,b);
  //copie dans matrice quaternionique
  QSetMatrix(r,g,b,dim[0],dim[1],&QMatrix);
  	
  //CHANGEMENT DE BASE
  qtnSympleticDecomposition((const Quaternion **)QMatrix,&dblMatR,&dblMatMu1,
    &dblMatMu2,&dblMatMu3, dim, qtnMu1, qtnMu2);   
  //on construit la matrice de la partie paralelle
  QMatInitZero(&QMatrixPara,dim);
  qtnConstructParaMat((const double **)dblMatMu1,&QMatrixPara,dim,qtnMu1);
  //on construit la matrice de la partie perpendiculaire
  QMatInitZero(&QMatrixPerp,dim);
  qtnConstructPerpMat((const double **)dblMatMu2,(const double **)dblMatMu3,&QMatrixPerp,dim,qtnMu2,qtnMu1Mu2);

  //remarque les changements de bases qui suivent ne sont pas necessaires, on obtient les memes resultats
  // sans changer de base car de toute maniere nos quaternions sont codés dans la base (1,i,j,k)
  // et lorsqu'on effectue des multiplications elles se font dans cette base.
  
  //on veut afficher la matrice de quaternions "parallele"
  //pour cela il faut la repasser dans la base (i,j,k) qui est utilisée de coder les couleurs
  //des images en rvb telles que nous l'avons defini pour enregistrer 
  qtnChangementDeBaseSurLesImag((const Quaternion **)QMatrixPara,
    &dblMatE1,&dblMatE2,&dblMatE3,dim,qtnI,qtnJ,qtnK);
  QSetMatrixFromDblMat((const double **)dblMatR,(const double **)dblMatE1,(const double **)dblMatE2,
    (const double **)dblMatE3,dim, &QMatrixParaIJK);
  QMat_Min_Max((const Quaternion **)QMatrixParaIJK,dim,&dblRMin,&dblRMax,
    &dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
  QMat_ChgScale_IntLvlMax((const Quaternion **)QMatrixParaIJK,dim,255,dblRMin,dblRMax,
    dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax,&QMatScaled);
  QGetMatrixImagPart(&r,&g,&b,dim[0],dim[1],QMatScaled);
  strcpybtw(&result,strFile,"-para",'.');
  MgkEcrireImgCouleur(result,r,g,b,dim);
  printf("\nimage %s creee\n\n",result);

  //on veut afficher la matrice de quaternions "perpendiculaire" ou "perplexe"
  //pour cela il faut la repasser dans la base (i,j,k) qui est utilisée de coder les couleurs
  //des images en rvb telles que nous l'avons defini pour enregistrer 
  qtnChangementDeBaseSurLesImag((const Quaternion **)QMatrixPerp,
    &dblMatE1,&dblMatE2,&dblMatE3,dim,qtnI,qtnJ,qtnK);
  QSetMatrixFromDblMat((const double **)dblMatR,(const double **)dblMatE1,(const double **)dblMatE2,
    (const double **)dblMatE3,dim, &QMatrixPerpIJK);
  QMat_Min_Max((const Quaternion **)QMatrixPerpIJK,dim,&dblRMin,&dblRMax,
    &dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
  QMat_ChgScale_IntLvlMax((const Quaternion **)QMatrixPerpIJK,dim,255,dblRMin,dblRMax,
    dblIMin,dblIMax,dblJMin,dblJMax,dblKMin,dblKMax,&QMatScaled);
  QGetMatrixImagPart(&r,&g,&b,dim[0],dim[1],QMatScaled);
  strcpybtw(&result,strFile,"-perp",'.');
  MgkEcrireImgCouleur(result,r,g,b,dim);
  printf("\nimage %s creee\n\n",result);

  //SOMME
  //on effectue la somme et on fait le changement de base inverse
  qtnRecomposition((const double **)dblMatMu1,(const double **)dblMatMu2,
    (const double **)dblMatMu3,dim,qtnMu1,qtnMu2,&QMatrixSomme);
  QMat_Min_Max((const Quaternion **)QMatrixSomme,dim,&dblRMin,&dblRMax,&dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
  QMat_ChgScale_IntLvlMax((const Quaternion **)QMatrixSomme,dim,255,dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,
    dblKMin,dblKMax,&QMatScaled);
  QGetMatrixImagPart(&r,&g,&b,dim[0],dim[1],QMatScaled);
  strcpybtw(&result,strFile,"-somme",'.');
  MgkEcrireImgCouleur(result,r,g,b,dim);
  printf("\nimage %s creee\n\n",result);
   
  //on libere
  QMatrixFree(dim[0],&QMatrix);
  QMatrixFree(dim[0],&QMatScaled);
  QMatrixFree(dim[0],&QMatrixPara);
  QMatrixFree(dim[0],&QMatrixPerp);
  QMatrixFree(dim[0],&QMatrixParaIJK);
  QMatrixFree(dim[0],&QMatrixPerpIJK);
  Matrix_Free_Double(dim[0],&dblMatR);
  Matrix_Free_Double(dim[0],&dblMatMu1);
  Matrix_Free_Double(dim[0],&dblMatMu2);
  Matrix_Free_Double(dim[0],&dblMatMu3);
  Matrix_Free_Double(dim[0],&dblMatE1);
  Matrix_Free_Double(dim[0],&dblMatE2);
  Matrix_Free_Double(dim[0],&dblMatE3);
  QMatrixFree(dim[0],&QMatrixSomme);  
  free(result);
  free(r);
  free(g);
  free(b);
}



int main_sangwine(int argc, char *argv[])
{
  qtnConstructSympleticDecomposition("C:/images/quaternions/SympleticDecomposition/lenna.bmp");
  qtnConstructSympleticDecomposition("C:/images/quaternions/SympleticDecomposition/baboon.bmp");
  qtnConstructSympleticDecomposition("C:/images/quaternions/SympleticDecomposition/tiffany.bmp");
  qtnConstructSympleticDecomposition("C:/images/quaternions/SympleticDecomposition/parrots1024.bmp"); //ok ca marche !!
  
}

