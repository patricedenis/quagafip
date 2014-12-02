#include "couleur.h"

void couleur_test1()
{
	double c1,c2,c3,c4,c5,c6;
	
    /*RGB01toRGB0255(0.,1.,0.,&c1,&c2,&c3);
	printf("RGB01 %lf %lf %lf\n",0.,1.,0.);*/
	
	/*H360S1V1toHSV255(double h1, double s1, double v1, double *h2, double *s2, double *v2);
	HSV255toH360S1V1(double h1, double s1, double v1, double *h2, double *s2, double *v2);*/
	
    printf("RGB0255 %lf %lf %lf\n",0.,0.,255.);
	
    RGB0255toRGB01(0.,0.,255.,&c1,&c2,&c3);
    printf("RGB01 %lf %lf %lf\n",c1,c2,c3);
    
    
    // r,g,b values are from 0 to 1
    // h = [0,360], s = [0,1], v = [0,1]
    //		if s == 0, then h = -1 (undefined)
    RGBtoHSVPat(c1,c2,c3,&c4,&c5,&c6);
	printf("RGBtoHSV %lf %lf %lf\n",c4,c5,c6);
	
	H360S1V1toHSV255(c4,c5,c6,&c1,&c2,&c3);
	printf("H360S1V1toHSV255 %lf %lf %lf\n",c1,c2,c3);
	 
    HSV255toH360S1V1(c1,c2,c3,&c4,&c5,&c6);
	printf("HSV255toH360S1V1 %lf %lf %lf\n",c4,c5,c6);
	
	
	// r,g,b values are from 0 to 1
	// h = [0,360], s = [0,1], v = [0,1]
	//		if s == 0, then h = -1 (undefined)
	HSVtoRGBPat(&c1,&c2,&c3,c4,c5,c6);
	printf("HSVtoRGB01 %lf %lf %lf\n",c1,c2,c3);
	

}

void couleur_test2()
{
    double H360,S1,V1,R1,G1,B1;
   	// r,g,b values are from 0 to 1
	// h = [0,360], s = [0,1], v = [0,1]
	//		if s == 0, then h = -1 (undefined)
    //couleur donnée par Gimp en bmp
    /*printf("R255G255B255 %lf %lf %lf\n",128.0,255.0,0.0);
    RGB0255toRGB01(128.,255.,0.,&R1,&G1,&B1);
    printf("R1G1B1 %lf %lf %lf\n",R1,G1,B1);*/
    //printf("R255G255B255 %lf %lf %lf\n",127.5,255.0,0.0);
    //RGB0255toRGB01(127.5,255.,0.,&R1,&G1,&B1);
    printf("R255G255B255 %lf %lf %lf\n",255.,12.0,0.0);
    RGB0255toRGB01(255.,12.,0.,&R1,&G1,&B1);
    RGBtoHSVPat(R1,G1,B1,&H360,&S1,&V1);
    printf("angle correspondant en degre %lf %lf %lf\n",H360,S1,V1);

    
}

void conversion()
{
	int choix;
	double r,g,b,y,u,v;
	printf("Conversion rgb->Yuv 1;Yuv->rgb 2\nChoix ?=(1 ou 2)");
	scanf("%d",&choix);
	printf("\n");
	if(choix == 1)
	{
		printf("Entrez les valeurs r,g,b\n");
		scanf("%lf %lf %lf",&r,&g,&b);
		printf("\n");
		rgb2yuv(&y,&u,&v,r,g,b);
		printf("r=%lf g=%lf b=%lf -> Y=%lf u=%lf v=%lf\n",r,g,b,y,u,v);
	}
	else
	{
		printf("Entrez les valeurs Y,u,v\n");
		scanf("%lf %lf %lf",&y,&u,&v);
		printf("\n");
		yuv2rgb(&r,&g,&b,y,u,v);
		printf("Y=%lf u=%lf v=%lf -> r=%lf g=%lf b=%lf\n",y,u,v,r,g,b);
	}
}


int main_couleur(int argc, char *argv[])
{
    //test1(); //test général
    
    //test particulier pour connaitre les couleurs rgb de l'angle
    //d'une TFQ "vert" cf interpretation graphique de l'angle 
    //cf Caractérisation de l'espage quaternionique numérique
	couleur_test2();
	
    return 0;
}

