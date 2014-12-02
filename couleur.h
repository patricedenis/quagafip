#ifndef _COULEUR_H_
#define _COULEUR_H_

#include "MinMax.h"
#include "bool.h"

void RGB01toRGB0255( double r, double g, double b, double *r2, double *g2, double *b2);
void RGB0255toRGB01( double r, double g, double b, double *r2, double *g2, double *b2);
void H360S1V1toHSV255(double h1, double s1, double v1, double *h2, double *s2, double *v2);
void HSV255toH360S1V1(double h1, double s1, double v1, double *h2, double *s2, double *v2);
void RGBtoHSVPat( double r, double g, double b, double *h, double *s, double *v );
void HSVtoRGBPat( double *r, double *g, double *b, double h, double s, double v );
void rgb2yuv(double *y,double *u,double *v,double r,double g,double b);
void yuv2rgb(double *r,double *g,double *b,double y,double u,double v);
#endif

