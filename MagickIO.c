#include "MagickIO.h"

//////////////////////////////////////////////////////////////////////////////

//                  I M A G E    M A G I C K

//////////////////////////////////////////////////////////////////////////////

void MgkEcrireImgCouleur(char *nom_file,double r[],double g[],double b[],int dim[])
{
    ImageInfo *Entete;
    Image *image;
    unsigned char *pixels_r,*pixels_g,*pixels_b; 
    int bx,by,ind;
    
    pixels_r = (unsigned char *)malloc(sizeof(unsigned char)*dim[0]*dim[1]);
    pixels_g = (unsigned char *)malloc(sizeof(unsigned char)*dim[0]*dim[1]);
    pixels_b = (unsigned char *)malloc(sizeof(unsigned char)*dim[0]*dim[1]);

    Entete = CloneImageInfo (NULL);
    image=AllocateImage (Entete);
    SetImageType(image,TrueColorType);
    image->rows=(unsigned int)dim[0];
    image->columns=(unsigned int)dim[1];
    ind=0;
      for (bx=0;bx<dim[0];bx++)
          for (by=0;by<dim[1];by++)
              {
                  pixels_r[ind]=(unsigned char)(r[ind]+0.5);
                  pixels_g[ind]=(unsigned char)(g[ind]+0.5);
                  pixels_b[ind]=(unsigned char)(b[ind]+0.5);
                  ind++;
              }
   ImportImagePixels (image, 0,0, (const unsigned long)dim[1],(const unsigned long)dim[0],"R", CharPixel,pixels_r);
   ImportImagePixels (image, 0,0, (const unsigned long)dim[1],(const unsigned long)dim[0],"G", CharPixel,pixels_g);
   ImportImagePixels (image, 0,0, (const unsigned long)dim[1],(const unsigned long)dim[0],"B", CharPixel,pixels_b);
   free(pixels_r);
   free(pixels_g);
   free(pixels_b); 

   (void) strcpy(image->filename,nom_file);
     WriteImage(Entete,image);
     DestroyImage(image);
     DestroyImageInfo(Entete);

 }
 
void MgkEcrireImgGris(char *nom_file,double intensite[],int dim[])
{
    ImageInfo *Entete;
    Image *image;
    unsigned char *pixels_intensite; 
    int bx,by,ind;
    
    pixels_intensite = (unsigned char *)malloc(sizeof(unsigned char)*dim[0]*dim[1]);
    
    Entete = CloneImageInfo (NULL);
    image=AllocateImage (Entete);
    SetImageType(image,TrueColorType);
    image->rows=(unsigned int)dim[0];
    image->columns=(unsigned int)dim[1];
    ind=0;
      for (bx=0;bx<dim[0];bx++)
          for (by=0;by<dim[1];by++)
              {
                  pixels_intensite[ind]=(unsigned char)(intensite[ind]+0.5);
                  ind++;
              }
   ImportImagePixels (image, 0,0, (const unsigned long)dim[1],(const unsigned long)dim[0],"I", CharPixel,pixels_intensite);
   free(pixels_intensite);

   (void) strcpy(image->filename,nom_file);
     WriteImage(Entete,image);
     DestroyImage(image);
     DestroyImageInfo(Entete);
}

void MgkTypeImage(char *nom_file,int dim[],int type[])
{
 ExceptionInfo  exception;
 Image *image;
 ImageInfo *image_info;

 GetExceptionInfo(&exception);
 image_info=CloneImageInfo((ImageInfo *) NULL);
 (void) strcpy(image_info->filename,nom_file);
 image=ReadImage(image_info,&exception);
  if (exception.severity != UndefinedException)
        CatchException(&exception);
   if (image == (Image *) NULL)
        exit(1);
  dim[0]=(int)image->rows;dim[1]=(int)image->columns;    
  *type =  IsGrayImage (image, &exception );
  if (exception.severity != UndefinedException)
        CatchException(&exception);
  DestroyImage(image);
  DestroyImageInfo(image_info);
  DestroyExceptionInfo(&exception);
}
  
void MgkLireImgCouleur(char *nom_file,double r[],double g[],double b[])
{
 ExceptionInfo  exception;
 Image *image;
 ImageInfo *image_info;
 int dim[2],bx,by,ind=0;
 unsigned char *t_r,*t_g,*t_b;
 
 GetExceptionInfo(&exception);
 image_info=CloneImageInfo((ImageInfo *) NULL);
 (void) strcpy(image_info->filename,nom_file);
 image=ReadImage(image_info,&exception);
  if (exception.severity != UndefinedException)
        CatchException(&exception);
   if (image == (Image *) NULL)
        exit(1);
   dim[0]=(int)image->rows;dim[1]=(int)image->columns;    

   t_r = (unsigned char *)malloc(sizeof(unsigned char)*dim[0]*dim[1]);
   t_g = (unsigned char *)malloc(sizeof(unsigned char)*dim[0]*dim[1]);
   t_b = (unsigned char *)malloc(sizeof(unsigned char)*dim[0]*dim[1]);


   ExportImagePixels ( image, 0,0, image->columns,image->rows, "R", CharPixel,t_r,&exception );
   if (exception.severity != UndefinedException)
      CatchException(&exception);
   ExportImagePixels ( image, 0,0, image->columns,image->rows, "G", CharPixel,t_g,&exception );
   if (exception.severity != UndefinedException)
      CatchException(&exception);
   ExportImagePixels ( image, 0,0, image->columns,image->rows, "B", CharPixel,t_b,&exception );
   if (exception.severity != UndefinedException)
      CatchException(&exception);
   DestroyImage(image);
   DestroyImageInfo(image_info);
   DestroyExceptionInfo(&exception);
  
   for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
         r[ind]=(double)t_r[ind];
         g[ind]=(double)t_g[ind];
         b[ind]=(double)t_b[ind];
         ind++;
      }
      free(t_r);
      free(t_g);
      free(t_b); 
}

void MgkLireImgGris(char *nom_file,double intensite[])
{
 ExceptionInfo  exception;
 Image *image;
 ImageInfo *image_info;
 int dim[2],bx,by,ind=0;
 unsigned char *t_intensite;
 
 GetExceptionInfo(&exception);
 image_info=CloneImageInfo((ImageInfo *) NULL);
 (void) strcpy(image_info->filename,nom_file);
 image=ReadImage(image_info,&exception);
  if (exception.severity != UndefinedException)
        CatchException(&exception);
   if (image == (Image *) NULL)
        exit(1);
   dim[0]=(int)image->rows;dim[1]=(int)image->columns;    

   t_intensite = (unsigned char *)malloc(sizeof(unsigned char)*dim[0]*dim[1]);

   ExportImagePixels ( image, 0,0, image->columns,image->rows, "I", CharPixel,t_intensite,&exception );
   if (exception.severity != UndefinedException)
      CatchException(&exception);
   DestroyImage(image);
   DestroyImageInfo(image_info);
   DestroyExceptionInfo(&exception);
  
   for (bx=0;bx<dim[0];bx++)
      for (by=0;by<dim[1];by++)
      {
         intensite[ind]=(double)t_intensite[ind];
         ind++;
      }
      free(t_intensite);
}


//Complexes

//cette fonction permet de remplir une matrice complexe directement en passant en argument
//le nom du fichier image en "niveau de gris" en entree.
void MgkLireImgGrisAndAllocAndCSetMatrix(char * strInFileName,int intDim[],Complex *** pcpxMat)
{
  int type,intInc;
  double *dblCmp1=NULL,*dblCmp2=NULL;
  
  //lecture de l'en-tete
  MgkTypeImage(strInFileName,intDim,&type);

  //allocation
  dblCmp1 = (double *)malloc(intDim[0]*intDim[1]*sizeof(double));
  dblCmp2 = (double *)malloc(intDim[0]*intDim[1]*sizeof(double));
  //remplissage
  MgkLireImgGris(strInFileName,dblCmp1);
  //on initialise la partie imaginaire a zero car c'est une matrice complexe "reelle"
  for(intInc=0;intInc<intDim[0]*intDim[1];intInc++)
    dblCmp2[intInc]=0.;
  
  //on alloue
  Complex_Matrix_Allocate(intDim[0],intDim[1],pcpxMat);
  CSetMatrix(dblCmp1,dblCmp2,intDim[0],intDim[1],pcpxMat);

  //on libere
  free(dblCmp1);dblCmp1=NULL;
  free(dblCmp2);dblCmp2=NULL;
}


void MgkCGetMatrixAndEcrireImgGris(int * dim,const Complex ** cpxMat,char * strFileNameOut)
{
  double * dblCmp1=NULL,* dblCmp2=NULL;
  int blnMatReal;
  double dblEpsilon=0.00000001;//rajouter epsilon en parametre
  
  //on verifie que la martie imaginaire est bien nulle pour ecrire une image NdG
  blnMatReal=IsCMatRealEpsilon(dim,(const Complex **)cpxMat,dblEpsilon);
  if(blnMatReal)
  {  
    printf("Image niveaux de gris valide (partie imaginaire de la matrice complexe nulle)\n");
    dblCmp1 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    dblCmp2 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    CGetMatrixDoubleTab(&dblCmp1,&dblCmp2,dim[0],dim[1],(const Complex **)cpxMat);
    free(dblCmp2);dblCmp2=NULL;
    //on ecrit le fichier de sortie
    MgkEcrireImgGris(strFileNameOut,dblCmp1,dim);
    printf("\nimage %s creee\n\n",strFileNameOut);
    free(dblCmp1);dblCmp1=NULL;
  }
  else
    printf("La matrice de complexes passee en argument n'est pas une matrice reelle donc l'image n'a pas ete creee\n");
}


//Doubles

void MgkLireImgGrisAndAllocAndDblSetMatrix(char * strInFileName,int intDim[],double *** pdblMat)
{
  int type,intInc;
  double *dblCmp1=NULL;
  
  //lecture de l'en-tete
  MgkTypeImage(strInFileName,intDim,&type);

  //allocation
  dblCmp1 = (double *)malloc(intDim[0]*intDim[1]*sizeof(double));
  //remplissage
  MgkLireImgGris(strInFileName,dblCmp1);

  //on alloue
  Matrix_Allocate_Double(intDim[0],intDim[1],pdblMat);
  SetDblMatrixFromTab(dblCmp1,intDim[0],intDim[1],pdblMat);
  //on libere
  free(dblCmp1);dblCmp1=NULL;
}


void MgkDblGetMatrixAndEcrireImgGris(int * dim,const double ** dblMat,char * strFileNameOut,FILE *stream)
{
  double * dblDataTab=NULL;
  
  //passer le tableau 2D de doubles sous forme de tableau 1D de doubles
  dblDataTab = (double *)malloc(dim[0]*dim[1]*sizeof(double));
  GetDblTabMatrix(dblMat,&dblDataTab,dim);
  
  //on ecrit le fichier de sortie
  MgkEcrireImgGris(strFileNameOut,dblDataTab,dim);
  fprintf(stream,"\nimage %s creee\n\n",strFileNameOut);
  free(dblDataTab);dblDataTab=NULL;
}



void MgkLireImgCouleurAndAllocAndCSetMatrixMarg(char * strInFileName,int intDim[],
  Complex *** pcpxMat1,Complex *** pcpxMat2,Complex *** pcpxMat3)
{
  int type,intInc;
  double *dblCmp1=NULL,*dblCmp2=NULL,*dblCmp3=NULL;
  double *dblCmp1i=NULL,*dblCmp2i=NULL,*dblCmp3i=NULL;
  
  //lecture de l'en-tete
  MgkTypeImage(strInFileName,intDim,&type);

  //allocation
  dblCmp1 = (double *)malloc(intDim[0]*intDim[1]*sizeof(double));
  dblCmp2 = (double *)malloc(intDim[0]*intDim[1]*sizeof(double));
  dblCmp3 = (double *)malloc(intDim[0]*intDim[1]*sizeof(double));
  
  //remplissage des tableaux
  MgkLireImgCouleur(strInFileName,dblCmp1,dblCmp2,dblCmp3);
 
  dblCmp1i = (double *)malloc(intDim[0]*intDim[1]*sizeof(double));
  dblCmp2i = (double *)malloc(intDim[0]*intDim[1]*sizeof(double));
  dblCmp3i = (double *)malloc(intDim[0]*intDim[1]*sizeof(double));
 
  //on initialise la partie imaginaire a zero car c'est une matrice complexe "reelle"
  for(intInc=0;intInc<intDim[0]*intDim[1];intInc++)
  {
    dblCmp1i[intInc]=0.;
    dblCmp2i[intInc]=0.;
    dblCmp3i[intInc]=0.;
  }
  
  //on alloue
  Complex_Matrix_Allocate(intDim[0],intDim[1],pcpxMat1);
  Complex_Init_Matrix(intDim[0],intDim[1],Complex_InitXY(0.,0.),pcpxMat1);
  //on remplit
  CSetMatrix(dblCmp1,dblCmp1i,intDim[0],intDim[1],pcpxMat1);
  //on alloue
  Complex_Matrix_Allocate(intDim[0],intDim[1],pcpxMat2);
  Complex_Init_Matrix(intDim[0],intDim[1],Complex_InitXY(0.,0.),pcpxMat2);
  //on remplit
  CSetMatrix(dblCmp2,dblCmp2i,intDim[0],intDim[1],pcpxMat2);
  //on alloue
  Complex_Matrix_Allocate(intDim[0],intDim[1],pcpxMat3);
  Complex_Init_Matrix(intDim[0],intDim[1],Complex_InitXY(0.,0.),pcpxMat3);
  //on remplit
  CSetMatrix(dblCmp3,dblCmp3i,intDim[0],intDim[1],pcpxMat3);

  //on libere
  free(dblCmp1);dblCmp1=NULL;
  free(dblCmp1i);dblCmp1i=NULL;
  free(dblCmp2);dblCmp2=NULL;
  free(dblCmp2i);dblCmp2i=NULL;
  free(dblCmp3);dblCmp3=NULL;
  free(dblCmp3i);dblCmp3i=NULL;
}

void MgkCGetMatrixAndEcrireImgCouleurMarg(int * dim,const Complex ** cpxMat1,const Complex ** cpxMat2,
  const Complex ** cpxMat3,char * strFileNameOut)
{
  double * dblCmp1=NULL,* dblCmp1i=NULL;
  double * dblCmp2=NULL,* dblCmp2i=NULL;
  double * dblCmp3=NULL,* dblCmp3i=NULL;

  int blnMatReal1,blnMatReal2,blnMatReal3;
  double dblEpsilon=0.00000001;//rajouter epsilon en parametre
  
  //on verifie que la martie imaginaire est bien nulle pour ecrire une image NdG
  blnMatReal1=IsCMatRealEpsilon(dim,cpxMat1,dblEpsilon);
  blnMatReal2=IsCMatRealEpsilon(dim,cpxMat2,dblEpsilon);
  blnMatReal3=IsCMatRealEpsilon(dim,cpxMat3,dblEpsilon);
  if(blnMatReal1&&blnMatReal2&&blnMatReal3)
  {  
    printf("Canaux couleur valides (parties imaginaires des matrices complexes nulles)\n");
    dblCmp1 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    dblCmp1i = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    CGetMatrixDoubleTab(&dblCmp1,&dblCmp1i,dim[0],dim[1],cpxMat1);
    dblCmp2 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    dblCmp2i = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    CGetMatrixDoubleTab(&dblCmp2,&dblCmp2i,dim[0],dim[1],cpxMat2);
    dblCmp3 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    dblCmp3i = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    CGetMatrixDoubleTab(&dblCmp3,&dblCmp3i,dim[0],dim[1],cpxMat3);

    //verifier que les tableaux de doubles imaginaires sont nuls.

//    printf("Partie imaginaire 1 nulle? 1(oui) 0(non): %d\n",

//      IsVectDoubleNullEpsilon(dblCmp1i,dim[0]*dim[1],dblEpsilon));
//    printf("Partie imaginaire 2 nulle? 1(oui) 0(non): %d\n",
//      IsVectDoubleNullEpsilon(dblCmp2i,dim[0]*dim[1],dblEpsilon));
//    printf("Partie imaginaire 3 nulle? 1(oui) 0(non): %d\n",
//      IsVectDoubleNullEpsilon(dblCmp3i,dim[0]*dim[1],dblEpsilon));
    free(dblCmp1i);dblCmp1i=NULL;
    free(dblCmp2i);dblCmp2i=NULL;
    free(dblCmp3i);dblCmp3i=NULL;        
  
    //on ecrit le fichier de sortie
    MgkEcrireImgCouleur(strFileNameOut,dblCmp1,dblCmp2,dblCmp3,dim);    
    printf("\nimage %s creee\n\n",strFileNameOut);
    free(dblCmp1);dblCmp1=NULL;
    free(dblCmp2);dblCmp2=NULL;
    free(dblCmp3);dblCmp3=NULL;
  }
  else
    printf("Au moins une des matrices de complexes passees en argument n'est pas une matrice reelle donc l'image n'a pas ete creee\n");
}


//l'allocation dynamique est comprise dans la fonction
void MgkLireImgCouleurAndAllocAndQSetMatrix(char * strInFileName,int intDim[],Quaternion *** pqtnMat)
{
  int type,intInc;
  double *dblCmp1=NULL,*dblCmp2=NULL,*dblCmp3=NULL;
  
  //lecture de l'en-tete
  MgkTypeImage(strInFileName,intDim,&type);

  //allocation
  dblCmp1 = (double *)malloc(intDim[0]*intDim[1]*sizeof(double));
  dblCmp2 = (double *)malloc(intDim[0]*intDim[1]*sizeof(double));
  dblCmp3 = (double *)malloc(intDim[0]*intDim[1]*sizeof(double));
  
  //remplissage des tableaux
  MgkLireImgCouleur(strInFileName,dblCmp1,dblCmp2,dblCmp3);
 
  //on alloue la matrice
  QMatrixAllocate(intDim[0],intDim[1],pqtnMat);
  QSetMatrix(dblCmp1,dblCmp2,dblCmp3,intDim[0],intDim[1],pqtnMat);

  //on libere
  free(dblCmp1);dblCmp1=NULL;
  free(dblCmp2);dblCmp2=NULL;
  free(dblCmp3);dblCmp3=NULL;
}

void MgkQGetMatrixAndEcrireImgCouleur(int * dim,Quaternion ** qtnMat,char * strFileNameOut)
{
  double * dblCmp1=NULL,* dblCmp2=NULL,* dblCmp3=NULL;
  int blnMatIm;
  double dblEpsilon=0.00000001;//rajouter epsilon en parametre
  
  
  //on verifie que la partie reelle est bien nulle pour ecrire une image couleur
  blnMatIm=IsQMatImEpsilon(dim,qtnMat,dblEpsilon);
  if(blnMatIm)
  {  
    printf("Image couleur valide (partie reelle de la matrice complexe nulle)\n");
    dblCmp1 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    dblCmp2 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    dblCmp3 = (double *)malloc(dim[0]*dim[1]*sizeof(double));
    
    QGetMatrixImagPart(&dblCmp1,&dblCmp2,&dblCmp3,dim[0],dim[1],qtnMat);

    //on ecrit le fichier de sortie
    MgkEcrireImgCouleur(strFileNameOut,dblCmp1,dblCmp2,dblCmp3,dim);
    printf("\nimage %s creee\n\n",strFileNameOut);
    free(dblCmp1);dblCmp1=NULL;
    free(dblCmp2);dblCmp2=NULL;
    free(dblCmp3);dblCmp3=NULL;
  }
  else
    printf("La matrice de quaternions passee en argument n'est pas une matrice imaginaire donc l'image n'a pas ete creee\n");
}
