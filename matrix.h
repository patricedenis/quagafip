#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdio.h>
#include <stdlib.h>

#include "bool.h"

int Matrix_Allocate(int intHeight,int intWidth,unsigned char *** pImage);
void Matrix_Free(int intHeight,unsigned char *** pImage);

#endif

