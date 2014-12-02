#include "couleur.h"

// r,g,b values are from 0 to 1
// r2,g2,b2 values are from 0 to 255
void RGB01toRGB0255( double r, double g, double b, double *r2, double *g2, double *b2)
{
    *r2 = r * 255.;
    *g2 = g * 255.;
    *b2 = b * 255.;
}

// r,g,b values are from 0 to 255
// r2,g2,b2 values are from 0 to 1
void RGB0255toRGB01( double r, double g, double b, double *r2, double *g2, double *b2)
{
    *r2 = r / 255.;
    *g2 = g / 255.;
    *b2 = b / 255.;
}


// h1 = [0,360], s1 = [0,1], v1 = [0,1]
// h2 = [0,255], s2 = [0,255], v2 = [0,255]
void H360S1V1toHSV255(double h1, double s1, double v1, double *h2, double *s2, double *v2)
{
    *h2 = (h1*255.)/360.;
    *s2 = s1*255.;
    *v2 = v1*255.;
}


// h1 = [0,255], s2 = [0,255], v2 = [0,255]
// h2 = [0,360], s1 = [0,1], v1 = [0,1]
void HSV255toH360S1V1(double h1, double s1, double v1, double *h2, double *s2, double *v2)
{
    *h2 = (h1*360.)/255.;
    *s2 = s1/255.;
    *v2 = v1/255.;
}

// r,g,b values are from 0 to 1
// h = [0,360], s = [0,1], v = [0,1]
//		if s == 0, then h = -1 (undefined)
void RGBtoHSVPat( double r, double g, double b, double *h, double *s, double *v )
{
	double min, max, delta;

	min = Denis_MIN3( r, g, b );
	max = Denis_MAX3( r, g, b );
	*v = max;				// v

	delta = max - min;

	if( max != 0 )
		*s = delta / max;		// s
	else {
		// r = g = b = 0		// s = 0, v is undefined
		*s = 0;
		*h = -1;
		return;
	}

	if( r == max )
		*h = ( g - b ) / delta;		// between yellow & magenta
	else if( g == max )
		*h = 2 + ( b - r ) / delta;	// between cyan & yellow
	else
		*h = 4 + ( r - g ) / delta;	// between magenta & cyan

	*h *= 60;				// degrees
	if( *h < 0 )
		*h += 360;

}

// r,g,b values are from 0 to 1
// h = [0,360], s = [0,1], v = [0,1]
//		if s == 0, then h = -1 (undefined)
void HSVtoRGBPat( double *r, double *g, double *b, double h, double s, double v )
{
	int i;
	double f, p, q, t;

	if( s == 0 ) {
		// achromatic (grey)
		*r = *g = *b = v;
		return;
	}
	//printf("\nh%lf\n",h);
	
	h /= 60;			// sector 0 to 5
	i = (int)( h );
	f = h - i;			// factorial part of h
	p = v * ( 1 - s );
	q = v * ( 1 - s * f );
	t = v * ( 1 - s * ( 1 - f ) );
	
    /*printf("\nh%lf,i%lf,f%lf\n",h,i,f);
	printf("\nh%lf,s%lf,v%lf\n",h,s,v);
	printf("\nv%lf,p%lf,q%lf,t%lf\n",v,p,q,t);*/
	switch( i ) {
		case 0:
			*r = v;
			*g = t;
			*b = p;
			break;
		case 1:
			*r = q;
			*g = v;
			*b = p;
			break;
		case 2:
			*r = p;
			*g = v;
			*b = t;
			break;
		case 3:
			*r = p;
			*g = q;
			*b = v;
			break;
		case 4:
			*r = t;
			*g = p;
			*b = v;
			break;
		default:		// case 5:
			*r = v;
			*g = p;
			*b = q;
			break;
	}

}


//intervalles de définition de Y, u et v pour r,g et compris entre 0 et 255
// Ymin=-11.215155		Ymax=483.577155
// umin=-446.125050		umax=193.522050
// vmin=-100.289460		vmax=160.775460

void rgb2yuv(double *y,double *u,double *v,double r,double g,double b)
{
	double trans[9] = {  0.47285400000000 ,  1.42352700000000 , -0.04398100000000 , -0.60773500000000 , -1.14177500000000 ,  0.75891000000000 , 0.58584600000000 , -0.39329200000000 ,  0.04464600000000};

	*y = trans[0]* r + trans[1] * g + trans[2] * b;
	*u = trans[3]* r + trans[4] * g + trans[5] * b;
	*v = trans[6]* r + trans[7] * g + trans[8] * b;
}

void yuv2rgb(double *r,double *g,double *b,double y,double u,double v)
{
	double trans[9] = { 0.33060074687364 , -0.06178943963832 ,  1.37599729167600 ,   0.63013437898936  , 0.06261724756816 , -0.44364366862929 ,  1.21277796507073  , 1.36240555895106  , 0.43443947807711 };

	*r = trans[0]* y + trans[1] * u + trans[2] * v;
	*g = trans[3]* y + trans[4] * u + trans[5] * v;
	*b = trans[6]* y + trans[7] * u + trans[8] * v;
}



