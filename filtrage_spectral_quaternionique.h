#ifndef _FILTRAGE_SPECTRAL_QUATERNIONIQUE_H_
#define _FILTRAGE_SPECTRAL_QUATERNIONIQUE_H_

#include <math.h>


#include "quaternion.h"
#include "fourier_quaternion.h"


void dblHammingWindow1D(double * pdblHammingVect,int intDim);
int main_filtrage_spectral_quaternionique(int argc, char *argv[]);

#endif
