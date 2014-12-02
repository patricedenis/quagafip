#include "matrix.h"

/*memory allocation to place the Image Matrix */
int Matrix_Allocate_Old(int intHeight,int intWidth,unsigned char *** pMat)
{
  int i,intCount,j;
  
  *pMat = (unsigned char **)malloc(intHeight*sizeof(unsigned char *));
  if ( (*pMat) != NULL) /*allocation successful*/
  {
    for(i=0;i<=intHeight-1;i++)
    {
      
      // place reservation for the three componants RGB
      (*pMat)[i]=(unsigned char *)malloc(3*intWidth*sizeof(unsigned char));
      if ((*pMat)[i] == NULL) //allocation error
      {
	printf("allocation simple matrix old memory error\n");
	//we need to free the matrix memory that is already allocated
	for(j=i-1;j<=0;j--)
	  free(((*pMat)[j]));
	free(*pMat);
	return FALSE;
      }
    }
  }
  else /*allocation error*/
  {
    printf("allocation simple matrix old memory error\n");
    return FALSE;
  }
  return TRUE;
}


/* we need to free the matrix memory */
void Matrix_Free_Old(int intHeight,unsigned char *** pMat)
{
  int i;
  
  for(i=0;i<=intHeight-1;i++)
      free((*pMat)[i]);
  free(*pMat);
  printf("simple matrix old memory free\n");
}


/*memory allocation to place the Matrix */
int Matrix_Allocate(int intHeight,int intWidth,unsigned char *** pMat)
{
  int i,intCount,j;
  
  *pMat = (unsigned char **)malloc(intHeight*sizeof(unsigned char *));
  if ( (*pMat) != NULL) /*allocation successful*/
  {
    for(i=0;i<=intHeight-1;i++)
    {
      
      // place reservation for one componant
      (*pMat)[i]=(unsigned char *)malloc(intWidth*sizeof(unsigned char));
      if ((*pMat)[i] == NULL) //allocation error
      {
	printf("allocation simple matrix memory error\n");
	//we need to free the matrix memory that is already allocated
	for(j=i-1;j<=0;j--)
	  free(((*pMat)[j]));
	free(*pMat);
	return FALSE;
      }
    }
  }
  else /*allocation error*/
  {
    printf("allocation simple matrix memory error\n");
    return FALSE;
  }
  return TRUE;
}


/* we need to free the matrix memory */
void Matrix_Free(int intHeight,unsigned char *** pMat)
{
  int i;
  
  for(i=0;i<=intHeight-1;i++)
      free((*pMat)[i]);
  free(*pMat);
  printf("simple matrix memory free\n");
}


