#include "fourier_quaternion.h"


//////////////////////////////////////////////////////////////////////////////

//                  Q   F   T

//////////////////////////////////////////////////////////////////////////////
void Quaternion_Matrix_Fourier_Transform_ExpLeft(Quaternion ** QMatrix,Quaternion *** pQMatrixFT,
    Quaternion QMu,int intHeight,int intWidth)
{
	int s,t,S,T;
	Quaternion QExp,QSum;
	double dblFactor,dblQuotient;
	int intInverse;

	dblFactor = -2. * M_PI ;
	dblQuotient = sqrt((double)intHeight*(double)intWidth);
 	
	// the frequency image coordinates will be describe
	// by the S and T index.
	// In the spatial image, the pixel will be pointed by the spatial 
	// coordinates s and t.

	for (S=0;S<intHeight;S++)
		for (T=0;T<intWidth;T++)
		{
			QSum = QInit(0.,0.,0.,0.);
			//QDisp(QSum,stdout); 
			for(s=0;s<intHeight;s++)
			{
				for(t=0;t<intWidth;t++)
				{
					QExp = QInitExp(QMu,dblFactor*(S*s/((double)intHeight)+T*t/((double)intWidth)));
					//QDisp(QExp,stdout);
					QSum = QAdd(QSum,QMult(QExp,QMatrix[s][t]));
				}
			}
			(*pQMatrixFT)[S][T].a = QSum.a/dblQuotient;		
			(*pQMatrixFT)[S][T].b = QSum.b/dblQuotient;	//	R(S,T)
			(*pQMatrixFT)[S][T].c = QSum.c/dblQuotient;	//	G(S,T)
			(*pQMatrixFT)[S][T].d = QSum.d/dblQuotient;	//	B(S,T)
		}
}

void Quaternion_Matrix_Fourier_Transform_ExpLeft_Inv(Quaternion ** QMatrix,Quaternion *** pQMatrixFT,
    Quaternion QMu,int intHeight,int intWidth)
{
	int s,t,S,T;
	Quaternion QExp,QSum;
	double dblFactor,dblQuotient;
	
	dblFactor = 2. * M_PI ;
	dblQuotient = sqrt((double)intHeight*(double)intWidth);
	
	// the frequency image coordinates will be describe
	// by the S and T index.
	// In the spatial image, the pixel will be pointed by the spatial 
	// coordinates s and t.

	for (s=0;s<intHeight;s++)
		for (t=0;t<intWidth;t++)
		{
			QSum = QInit(0.,0.,0.,0.);
			//QDisp(QSum,stdout);  
			for(S=0;S<intHeight;S++)
			{
				for(T=0;T<intWidth;T++)
				{
					QExp = QInitExp(QMu,dblFactor*(S*s/((double)intHeight)+T*t/((double)intWidth)));
					//QDisp(QExp,stdout);
					QSum = QAdd(QSum,QMult(QExp,QMatrix[S][T]));
				}
			}
			(*pQMatrixFT)[s][t].a = QSum.a/dblQuotient;		
			(*pQMatrixFT)[s][t].b = QSum.b/dblQuotient;	//	R(s,t)
			(*pQMatrixFT)[s][t].c = QSum.c/dblQuotient;	//	G(s,t)
			(*pQMatrixFT)[s][t].d = QSum.d/dblQuotient;	//	B(s,t)
		}
}


//////////////////////////////////////////////////////////////////////////////

//              FONCTIONS TEST    is real?   is pure?

//////////////////////////////////////////////////////////////////////////////

int IsRealPresent(Quaternion ** QMat,int *dim,double epsilon)
{
  int bx,by;
    
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        //on teste avec un epsilon car il subsiste souvent des elements non nuls dus aux calculs numeriques
        if(QReal((const Quaternion)QMat[bx][by]) > epsilon )
                return TRUE;
      }
    return FALSE;
}

//permet de savoir si la matrice passée en argument est une matrice de quaternions
//purs : TRUE ou s'il existe des éléments possédant une partie réelle
int Is_Pure_Image(Quaternion ** mat,int dim[])
{
    int bx,by;
    for (bx=0;bx<dim[0];bx++)
       for (by=0;by<dim[1];by++)
         if((int)(mat[bx][by].a + 0.5)!=0 )
           return FALSE;
    return TRUE;
}

int Is_Null_Mat(double * dblMat,int dim[])
{
    int bx,by,ind=0;
    for (bx=0;bx<dim[0];bx++)
       for (by=0;by<dim[1];by++)
       {
         if(dblMat[ind] !=0 )
           return FALSE;
         ind++;
       }
    return TRUE;
}




void QLogMatrix(const Quaternion ** mat, int * dim,Quaternion *** pLogmat)
{
    int bx,by;
    
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        ((*pLogmat)[bx][by]).a=log((mat[bx][by]).a);
        ((*pLogmat)[bx][by]).b=log((mat[bx][by]).b);
        ((*pLogmat)[bx][by]).c=log((mat[bx][by]).c);
        ((*pLogmat)[bx][by]).d=log((mat[bx][by]).d);
      }  
}





//////////////////////////////////////////////////////////////////////////////

//                  FONCTIONS SUR IMAGES

//////////////////////////////////////////////////////////////////////////////


int CompareImage(char * Im1,char * Im2)
{
    int dim1[2],dim2[2],type;
    double *r1,*r2,*g1,*g2,*b1,*b2;
    int bx,by,ind=0;
    

    MgkTypeImage(Im1,dim1,&type);
    MgkTypeImage(Im2,dim2,&type);

    if((dim1[0]!=dim2[0])||(dim1[1]!=dim2[1]))
        return FALSE;
        
    //allocation des tableaux r,g et b
    r1 = (double *)malloc(dim1[0]*dim1[1]*sizeof(double));
    g1 = (double *)malloc(dim1[0]*dim1[1]*sizeof(double));
    b1 = (double *)malloc(dim1[0]*dim1[1]*sizeof(double));
    //remplissage des tableaux
    MgkLireImgCouleur(Im1,r1,g1,b1);
    


    //allocation des tableaux r,g et b
    r2 = (double *)malloc(dim2[0]*dim2[1]*sizeof(double));
    g2 = (double *)malloc(dim2[0]*dim2[1]*sizeof(double));
    b2 = (double *)malloc(dim2[0]*dim2[1]*sizeof(double));
    //remplissage des tableaux
    MgkLireImgCouleur(Im2,r2,g2,b2);
    

    for (bx=0;bx<dim1[0];bx++)
      for (by=0;by<dim1[1];by++)
      {
         if((r1[ind]!=r2[ind])||(g1[ind]!=g2[ind])||(b1[ind]!=b2[ind]))
         {
           free(r1);
           free(g1);
           free(b1);
           free(r2);
           free(g2);
           free(b2);
           return FALSE;
         }
         ind++;
      }
    free(r1);
    free(g1);
    free(b1);
    free(r2);
    free(g2);
    free(b2);
    return TRUE;
}




//////////////////////////////////////////////////////////////////////////////

//                  

//////////////////////////////////////////////////////////////////////////////

Matrix_Init_Freq(Quaternion ***QMatrix,double amplitude,int offsetX,int offsetY,int colour,int intHeight,int intWidth)
{
    switch (colour) 
    {
        case 1 : (*QMatrix)[intHeight/2+offsetY][intWidth/2+offsetX].a = amplitude;
        case 2 : (*QMatrix)[intHeight/2+offsetY][intWidth/2+offsetX].b = amplitude;
        case 3 : (*QMatrix)[intHeight/2+offsetY][intWidth/2+offsetX].c = amplitude;
        case 4 : (*QMatrix)[intHeight/2+offsetY][intWidth/2+offsetX].d = amplitude;
        default : ;
    }
}

//on calcule ici pour chaque point de la matrice, son module, son angle et son axe.
Calcul_Matrice_LogRo_Mu_Phi(const Quaternion ** QMatFourier,double *** dblMatLogModule,
double *** dblMatAngle,Quaternion *** QMatAxe,int intHeight,int intWidth)
{
    int dim[2],bx,by;
    double LogRo,Phi;
    Quaternion Mu;
    
    dim[0]=intHeight;
    dim[1]=intWidth;
    
    //la matrice QFourier contient les informations fréquentielles de l'image
    for (bx=0;bx<dim[0];bx++)
          for (by=0;by<dim[1];by++)
              {
                  QExp_LogRo_Mu_Phi(QMatFourier[bx][by],&LogRo,&Mu,&Phi);
                  //on veut construire une vue du module de cette matrice
                  (*dblMatLogModule)[bx][by]=LogRo;
                  //on veut enfin une vue de l'axe de la matrice
                  (*QMatAxe)[bx][by]=Mu; //Mu quaternion pur
                  //on veut construire aussi une vue de l'angle de la matrice
                  (*dblMatAngle)[bx][by]=Phi;
              }
}

//on calcule ici pour chaque point de la matrice, son module, son angle et son axe.
Calcul_Matrice_Ro_Mu_Phi(const Quaternion ** QMatFourier,double *** dblMatModule,
double *** dblMatAngle,Quaternion *** QMatAxe,int intHeight,int intWidth)
{
    int dim[2],bx,by;
    double Ro,Phi;
    Quaternion Mu;
    
    dim[0]=intHeight;
    dim[1]=intWidth;
    
    //la matrice QFourier contient les informations fréquentielles de l'image
    for (bx=0;bx<dim[0];bx++)
          for (by=0;by<dim[1];by++)
              {
                  QExp_Ro_Mu_Phi(QMatFourier[bx][by],&Ro,&Mu,&Phi);
                  //on veut construire une vue du module de cette matrice
                  (*dblMatModule)[bx][by]=Ro;
                  //on veut enfin une vue de l'axe de la matrice
                  (*QMatAxe)[bx][by]=Mu; //Mu quaternion pur
                  //on veut construire aussi une vue de l'angle de la matrice
                  (*dblMatAngle)[bx][by]=Phi;
              }
}


//////////////////////////////////////////////////////////////////////////////

//                  M A S Q U E S

//////////////////////////////////////////////////////////////////////////////

// On construit un masque "modulusMaskMat" a partir de la matrice de module "modulusMat"
// en selectionnant seulement les valeurs du module au dessus du seuil "dblthreshold" qui
// correspond a un pourcentage de la valeur maximale du module.
// donc 0 <= dblThreshold <= 1
void Mask_From_LogModulus(double dblThreshold,int * dim,double ** dblLogModulusMat,unsigned char *** pblnLogModulusMaskMat)
{
    double dblMin,dblMax,dblSeuil;
    int bx,by;
    
    Dbl_Matrix_Min_Max((const double **)dblLogModulusMat,dim,&dblMin,&dblMax);
    dblSeuil=dblThreshold*dblMax;
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
            if(dblLogModulusMat[bx][by]>=dblSeuil)
              (*pblnLogModulusMaskMat)[bx][by]=TRUE;
            else
              (*pblnLogModulusMaskMat)[bx][by]=FALSE;
      }
}


void dblMat_Mask_Add(int * dim,double ** dblMaskMat,double ** dblMat,double *** pdblAddMat)
{
    int bx,by;
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
        (*pdblAddMat)[bx][by]=dblMaskMat[bx][by]+dblMat[bx][by];
}

void QMat_Angle_Construct_WithLogModMask(int * dim,unsigned char ** blnLogMaskMat,double ** dblMat,
Quaternion *** pQAngleColourMasked)
{
    int bx,by;
    double r,g,b;
    
    // on regarde la ou l'information est pertinente
    // et ensuite on construit la matrice d'angle représentée dans l'espace couleur HSV
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        if(blnLogMaskMat[bx][by]==TRUE)
        {
          //la teinte sera representée sur la premiere composante de l'espace HSV : Hue comme teinte
          // avec une saturation ainsi qu'une clarté maximales = 1
          // on fera ensuite la conversion vers l'espace couleur RGB
          // pour pouvoir afficher correctement les couleurs de l'angle
          
          //h = MOD360(((dblMat[bx][by]*180)/M_PI));
          //HSVtoRGBPat(&r,&g,&b,h,1.,1.);//ici r,g et b sont entre 0 et 1
          
          HSVtoRGBPat(&r,&g,&b,MOD360(((dblMat[bx][by]*180)/M_PI)),1.,1.);
          (*pQAngleColourMasked)[bx][by].a=0.;
          
          RGB01toRGB0255(r,g,b,&(*pQAngleColourMasked)[bx][by].b,
          &(*pQAngleColourMasked)[bx][by].c,&(*pQAngleColourMasked)[bx][by].d); //ici les sorties sont entre 0 et 255     
        }
        else
        {
          //(*pQAngleColourMasked)[bx][by] = QInit(0.,0.,0.,0.);
          (*pQAngleColourMasked)[bx][by].a=0.;
          (*pQAngleColourMasked)[bx][by].b=0.;
          (*pQAngleColourMasked)[bx][by].c=0.;
          (*pQAngleColourMasked)[bx][by].d=0.;
        }
      }
}

/////////////////////////////////////////////////////////////////////////////////////
void QMat_Angle_Construct_WithModMask(int * dim,unsigned char ** blnMaskMat,double ** dblMat,
Quaternion *** pQAngleColourMasked)
{
    int bx,by;
    double r,g,b;
    
    // on regarde la ou l'information est pertinente
    // et ensuite on construit la matrice d'angle représentée dans l'espace couleur HSV
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        if(blnMaskMat[bx][by]==1)
        {
          //la teinte sera representée sur la premiere composante de l'espace HSV
          // avec une saturation ainsi qu'une clarté maximales
          // on fera ensuite la conversion vers l'espace couleur RGB
          // pour pouvoir afficher correctement les couleurs de l'angle
          
          
          HSVtoRGBPat(&r,&g,&b,MOD360(((dblMat[bx][by]*180)/M_PI)),1.,1.);
          
          (*pQAngleColourMasked)[bx][by].a=0.;
          
          RGB01toRGB0255(r,g,b,&(*pQAngleColourMasked)[bx][by].b,
          &(*pQAngleColourMasked)[bx][by].c,&(*pQAngleColourMasked)[bx][by].d);        
        }
        else
        {
          (*pQAngleColourMasked)[bx][by].a=0.;
          (*pQAngleColourMasked)[bx][by].b=0.;
          (*pQAngleColourMasked)[bx][by].c=0.;
          (*pQAngleColourMasked)[bx][by].d=0.;                
        }
      }
}
/////////////////////////////////////////////////////////////////////////////////////


void dblMat_Mask_Product(int * dim,double ** dblMaskMat,double ** dblMat,double *** pdblProductMat)
{
    int bx,by;
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
        (*pdblProductMat)[bx][by]=dblMaskMat[bx][by]*dblMat[bx][by];
}

void QMat_Axe_Construct_WithLogModMask(int * dim,unsigned char ** blnLogMaskMat,Quaternion ** QPureMat,Quaternion *** pQPureAxeMatMasked)
{
    int bx,by;
    for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
        if(blnLogMaskMat[bx][by]==TRUE)
        {
          (*pQPureAxeMatMasked)[bx][by].a=QPureMat[bx][by].a;
          (*pQPureAxeMatMasked)[bx][by].b=QPureMat[bx][by].b;
          (*pQPureAxeMatMasked)[bx][by].c=QPureMat[bx][by].c;
          (*pQPureAxeMatMasked)[bx][by].d=QPureMat[bx][by].d;
        }
        else
        {
          (*pQPureAxeMatMasked)[bx][by].a=0.;
          (*pQPureAxeMatMasked)[bx][by].b=0.;
          (*pQPureAxeMatMasked)[bx][by].c=0.;
          (*pQPureAxeMatMasked)[bx][by].d=0.;
        }
      }
}
//////////////////////////////////////////////////////////////////////////////

//                           A F F I C H A G E

//////////////////////////////////////////////////////////////////////////////
void affiche(Quaternion ** mat,int dim[],char composante)
{
    int bx,by;
    switch(composante)
    {
      case 'a' : printf("composante r\n");
               for (bx=0;bx<dim[0];bx++)
               {
                 for (by=0;by<dim[1];by++)
                   printf("| %lf ",mat[bx][by].a);
                 printf("|\n");
               }
               break;
      case 'b' : printf("composante i\n");
               for (bx=0;bx<dim[0];bx++)
               {
                 for (by=0;by<dim[1];by++)
                   printf("| %lf ",mat[bx][by].b);
                 printf("|\n");
               }
               break;
      case 'c' : printf("composante j\n");
               for (bx=0;bx<dim[0];bx++)
               {
                 for (by=0;by<dim[1];by++)
                   printf("| %lf ",mat[bx][by].c);
                 printf("|\n");
               }
               break;
      case 'd' : printf("composante k\n");
               for (bx=0;bx<dim[0];bx++)
               {
                 for (by=0;by<dim[1];by++)
                   printf("| %lf ",mat[bx][by].d);
                 printf("|\n");
               }
               break;
      default : printf("erreur\n");
    }
}

void Verifie_Symetries(char * strFile, Quaternion QMu)
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
    
	printf("quaternionic fourier transform is processing ...\n");
    Quaternion_Matrix_Fourier_Transform_ExpLeft(QMatrix,&QMatFourier,QMu,dim[0],dim[1]);
    printf("done.\n");
    QMatrixFree(dim[0],&QMatrix);
       
    //QMatrixDisp(QMatFourier,dim[0],dim[1],1,stdout);
    //QMatrixDisp(QMatFourier,dim[0],dim[1],2,stdout);
    QMatrixDisp((const Quaternion **)QMatFourier,dim[0],dim[1],5,stdout);
    //il faut aussi verifier les symetries dans le spectre frequentiel
    dblPrecision = 0.001; //deux nombres différents d'au moins ce seuil compteront comme erreur de calcul
    dblTolerance = 0.0; // 0% d'erreur autorisée pour dire que la symetrie est respectée
    //blnSym=symetrie_partieR_OK(QMatFourier,dim,dblPrecision,dblTolerance);
    //printf("la partie reelle presente les symetries desirees? %d (oui 1, non 0)\n",blnSym);
    //blnSym=symetrie_partieImI_OK(QMatFourier,dim,dblPrecision,dblTolerance);
    //printf("la partie imaginaire I presente les symetries desirees? %d (oui 1, non 0)\n",blnSym);
 
	QMatrixFree(dim[0],&QMatFourier);
}

//////////////////////////////////////////////////////////////////////////////

//            P R O C E D U R E S       P R I N C I P A L E S 

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


int Variation_Couleur_RGB()
{
    Quaternion QMu;
    int intHeight,intWidth,intNbStripe;
    double InvRacine3;
    double K=2000.;
    int dim[2];
 
    Quaternion **QMatrix,**QMatrixShifted, **QMatOUT;
    
    /*64 64 2 0.3333333 0.3333333 0.3333333*/
    
        intHeight = 8;
        intWidth = 8;
        intNbStripe = 2;
        InvRacine3 = 1./sqrt(3.);
	dim[0]=intHeight;
	dim[1]=intWidth;
	
        
   	//the Analysis will follow the axis given by QMu
   	//QMu = QInit(0.,1./3.,1./3.,1./3.);
    //QMu = QInit(0.,InvRacine3,InvRacine3,InvRacine3); //Nu luminance
    //QMu2 = QInit(0.,5./sqrt(50.),3./sqrt(50.),4./sqrt(50.)); //Nu quelconque imag pur
    //QMu = QInit(0.,1.,0.,0.);//rouge
    //QMu2 = QInit(0.,0.,1.,0.);//vert
    //QMu2 = QInit(0.,0.,0.,1.);//bleu
    //QMu2 = QInit(0.,1./sqrt(2.),0.,1./sqrt(2.));//rouge-bleu
    //QMu = QInit(0.,1./sqrt(2.),1./sqrt(2.),0.);//rouge-vert : jaune
    

	// Allocation for a Quaterniotic matrix
	if(QMatrixAllocate(intHeight,intWidth,&QMatrix)==TRUE)
	{
		//Initialisation of this matrix with zeros
		QInitMatrix(intHeight,intWidth,QInit(0.,0.,0.,0.),&QMatrix);
	            
 		//j'initialise sur la partie réelle
 		QMatrix[intHeight/2 + intNbStripe][intHeight/2 + intNbStripe].a =K; 
 		QMatrix[intHeight/2 - intNbStripe][intHeight/2 - intNbStripe].a =-K;
 		
 		/*QMatrix[intHeight/2 + intNbStripe][intHeight/2 + intNbStripe].b =K; 
 		QMatrix[intHeight/2 - intNbStripe][intHeight/2 - intNbStripe].b =K;
 		
 		QMatrix[intHeight/2 + intNbStripe][intHeight/2 + intNbStripe].c =K; 
 		QMatrix[intHeight/2 - intNbStripe][intHeight/2 - intNbStripe].c =K;
 		
 		QMatrix[intHeight/2 + intNbStripe][intHeight/2 + intNbStripe].d =K; 
 		QMatrix[intHeight/2 - intNbStripe][intHeight/2 - intNbStripe].d =K;*/
		
		//initialization of another quaternionic matrix
		if(QMatrixAllocate(intHeight,intWidth,&QMatrixShifted)==TRUE)
		{
			//we need to shift the matrix before performing te inverse Fourier transform
			
			QMatrixShift(QMatrix,&QMatrixShifted,intHeight,intWidth);
			QMatrixFree(intHeight,&QMatrix);
			
			//initialisation of the matrix that will receive the insverse Fourier image
			if (QMatrixAllocate(intHeight,intWidth,&QMatOUT)==TRUE)
			{
				// perform the inverse fourier transform
				printf("fourier transform is processing ...\n");
               Quaternion_Matrix_Fourier_Transform_ExpLeft_Inv(QMatrixShifted,&QMatOUT,QMu,intHeight,intWidth);
               QMatrixFree(intHeight,&QMatrixShifted);
               QMatrixDisp((const Quaternion **)QMatOUT,dim[0],dim[1],5,stdout);
               printf("element reel (non nul) present ? %d (0 non , 1 oui)\n",IsRealPresent(QMatOUT,dim,0.00001));
            } 
            QMatrixFree(intHeight,&QMatOUT);
        }
    }
}


// strFile pointe vers l'image à ouvrir
// QMu indique la direction de la transformée de Fourier Utilisée dans cette fonction
void Construct_Vues_Frequentielles(char * strFile, Quaternion QMu, double seuil)
{
    
    double **dblMatAngle,**dblMatLogModule,**dblMatLogModuleScaled,
        **dblMatAngleScaled,**dblMatAngleMod;
    int i,j,dim[2],type,blnReal;
    unsigned char ** blnLogModulusMaskMat;
    double *r,*g,*b, *dblData;
    Quaternion **QMatrix,**QMatFourier,**QMatrixShifted,**QMatAxe,
        **QMatAxeScaled,**QPureAxeMatMasked,**QAngleColourMasked,**QMatAxeScaled2;
    
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

	printf("quaternionic fourier transform is processing ...\n");
    Quaternion_Matrix_Fourier_Transform_ExpLeft(QMatrix,&QMatFourier,QMu,dim[0],dim[1]);
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
    blnReal=IsRealPresent(QMatrixShifted,dim,0.00001);
    printf("La partie réelle de la transformée de Fourier Quaternionique est non nulle (1 oui, 0 non):%d\n",blnReal);
    
    
    // on alloue les différentes matrices qui vont nous servir juste après.
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatLogModule);
    Matrix_Allocate_Double(dim[0],dim[1],&dblMatAngle);
    QMatrixAllocate(dim[0],dim[1],&QMatAxe);
    
    //puis on calcule pour chaque point de la matrice, son module, son angle et son axe.
    Calcul_Matrice_LogRo_Mu_Phi((const Quaternion **)QMatrixShifted,&dblMatLogModule,&dblMatAngle,&QMatAxe,dim[0],dim[1]);
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
    
    
    //GetUsgCharMatrix(Spectre255,&dblData,dim);
    GetDblTabMatrix(dblMatLogModuleScaled,&dblData,dim);
    Matrix_Free_Double(dim[0],&dblMatLogModuleScaled);
    result = (char *)malloc((10+strlen(strFile))*sizeof(char));
    strcpybtw(&result,strFile,"-qft-mod",'.');
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
    
    //Matrix_Allocate_Double(dim[0],dim[1],&dblMatAngleScaled);
    
    //Dbl_Mat_ChgScale_LvlMin_LvlMax((const double **) dblMatAngle,dim,0.,2*M_PI,dblMin,dblMax,&dblMatAngleScaled);
    //Matrix_Free_Double(dim[0],&dblMatAngle);
    
    QMatrixAllocate(dim[0],dim[1],&QAngleColourMasked);
    QMat_Angle_Construct_WithLogModMask(dim,blnLogModulusMaskMat,dblMatAngle,&QAngleColourMasked);    
    Matrix_Free_Double(dim[0],&dblMatAngle);
    
    r = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    g = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    b = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    QGetMatrixImagPart(&r,&g,&b,dim[0],dim[1],QAngleColourMasked);
    QMatrixFree(dim[0],&QAngleColourMasked);

    result = (char *)malloc((10+strlen(strFile))*sizeof(char));
    strcpybtw(&result,strFile,"-qft-angle",'.');
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
    QMat_Min_Max2((const Quaternion **)QMatAxe,dim,&dblMin,&dblMax);


    QMatrixAllocate(dim[0],dim[1],&QPureAxeMatMasked);
    QMat_Axe_Construct_WithLogModMask(dim,blnLogModulusMaskMat,QMatAxe,&QPureAxeMatMasked);
    Matrix_Free(dim[0],&blnLogModulusMaskMat);
    QMatrixFree(dim[0],&QMatAxe);
    
    
    QMat_Min_Max((const Quaternion **)QPureAxeMatMasked,dim,&dblRMin,&dblRMax,&dblIMin,&dblIMax,&dblJMin,&dblJMax,&dblKMin,&dblKMax);
    
    QMatrixAllocate(dim[0],dim[1],&QMatAxeScaled);
    
    QMat_ChgScale_IntLvlMax((const Quaternion **)QPureAxeMatMasked,dim,255,dblRMin,dblRMax,dblIMin,dblIMax,dblJMin,dblJMax,
    dblKMin,dblKMax,&QMatAxeScaled);
    
    QMatrixFree(dim[0],&QPureAxeMatMasked);

    r = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    g = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    b = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    QGetMatrixImagPart(&r,&g,&b,dim[0],dim[1],QMatAxeScaled);
    QMatrixFree(dim[0],&QMatAxeScaled);
    result = malloc((10+strlen(strFile))*sizeof(char));
    strcpybtw(&result,strFile,"-qft-axe",'.');
    MgkEcrireImgCouleur(result,r,g,b,dim);
    printf("\nimage %s creee\n\n",result);
    free(result);
    free(r);
	free(g);
	free(b); 
 	
	
}


//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


int main_fourier_quaternion(int argc, char *argv[])
{
    double dblQMuRedPart,dblQMuGreenPart,dblQMuBluePart;
    double InvRacine3;
    double seuil;//seuil définissant un masque de module
    int dim[2]={64,64};//dimension de l'image
    
    Quaternion QMu;//Quaternion pur utilisé comme direction de l'analyse de Fourier
    int intNbStripe;//nbre de raies souhaitées apres une IQFT
    //double K=2000.;//cste d'initialisation dans le domaine de Fourier bien pour RGB
    //double K=20.; //K correspondant a la bonne valeur pour Variation couleur Yuv
    
    //Initialisation des paramètres par défaut
    InvRacine3 = 1./sqrt(3.);

    //Différentes valeurs de Mu possibles avec Mu Quaternion unitaire pur
    //QMu = QInit(0.,InvRacine3,InvRacine3,InvRacine3); //Mu luminance
    QMu = QInit(0.,1./1.,0,0); //Mu rouge
    //QMu = QInit(0.,0,1,0); //Mu vert
    //QMu = QInit(0.,0,0,1); //Mu bleu
    //QMu = QInit(0.,1./sqrt(2.),1./sqrt(2.),0); //Mu jaune

    //le seuil pour le masque de module
    seuil = 0.23;

    /*printf("test min max\n");
    printf("test min(30.5,-6,6)=%lf\n",MIN3(30.5,-6,6.));
    printf("test max(30.5,-6,6)=%lf\n",MAX3(30.5,-6,6.));*/


    printf("fichier :%s\t + chaine en argument : %s\n",argv[1],argv[2]);

    Construct_Vues_Frequentielles(argv[1],QMu,seuil);
    //Verifie_Symetries(argv[1],QMu);
    
    //Variation_Couleur_RGB();

}


