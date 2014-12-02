#include "quaternion.h"

void test_operations_base()
{
    Quaternion a,b,c,d,e,f,h,i,j,k,m,o,p;
    double g,l,n;
    
    //Init
    a = QInit(1.,2.,3.,4.);
    b = QInit(0.,1.,2.,-1.);
    c = QAdd(a,b);
    d = QMult(a,b);
    e = QDiff(a,b);
    f = QScalAdd(a,2.);
    g = QScalProd(a,b);
    h = QScalMult(a,2.);
    i = QScalDiv(a,2.);
    j = QConj(a);
    k = QOpp(a);
    l = QReal(a);
    m = QImag(a);
    n = QNorm(a);
    o = QInv(a);
    p = QDiv(a,b);
    
    //Display
    QDisp(a,stdout);  
    QDisp(b,stdout);
    printf("somme\n");
    QDisp(c,stdout);
    printf("produit\n");
    QDisp(d,stdout);
    printf("diff\n");
    QDisp(e,stdout);
    printf("add scal\n");
    QDisp(f,stdout);
    printf("produit scalaire %lf\n",g);
    printf("multi / scal\n");
    QDisp(h,stdout);
    printf("div / scal\n");
    QDisp(i,stdout);
    printf("conjugué\n");
    QDisp(j,stdout);
    printf("opposé\n");
    QDisp(k,stdout);
    printf("partie reelle %lf\n",l);
    printf("partie imaginaire\n");
    QDisp(m,stdout);
    printf("Norme %lf\n",n);
    printf("Inverse\n");
    QDisp(o,stdout);    
    printf("Division\n");
    QDisp(p,stdout);
}

void produit_quaternionique(int argc, char *argv[])
{
  Quaternion q1,q2,q3;
  double dbl1,dbl2,dbl3,dbl4,dbl5,dbl6,dbl7,dbl8;
  
  if(argc!=9)
    printf("Arguments : q1=$1+$2i+$3j+$4k, q2=$5+$6i+$7j+$8k\n");
  else
  {
    sscanf(argv[1],"%lf",&dbl1);
    sscanf(argv[2],"%lf",&dbl2);
    sscanf(argv[3],"%lf",&dbl3);
    sscanf(argv[4],"%lf",&dbl4);
    sscanf(argv[5],"%lf",&dbl5);
    sscanf(argv[6],"%lf",&dbl6);
    sscanf(argv[7],"%lf",&dbl7);
    sscanf(argv[8],"%lf",&dbl8);
    
    printf("Produit Quaternionique\n");
    q1=QInit(dbl1,dbl2,dbl3,dbl4);
    q2=QInit(dbl5,dbl6,dbl7,dbl8);
    printf("Q1: ");QDisp(q1,stdout);printf("\n");
    printf("Q2: ");QDisp(q2,stdout);printf("\n");
    printf("Q3: ");QDisp(QMult(q1,q2),stdout);printf("\n");
  }
}

void calcul()
{
  Quaternion q,mu1,mu2,mu3;
  Quaternion x,y,z;

  
  q=QInit(0.,2.,-3.,1.);
  mu1=QInit(0.,1/sqrt(3.),1/sqrt(3.),1/sqrt(3.));
  mu2=QInit(0.,1/sqrt(2.),-1/sqrt(2.),0.);
  mu3=QMult(mu1,mu2);
  
  printf("Q: ");QDisp(q,stdout);printf("\n");
  printf("Mu1: ");QDisp(mu1,stdout);printf("\n");
  printf("Mu2: ");QDisp(mu2,stdout);printf("\n");
  printf("Mu3: ");QDisp(mu3,stdout);printf("\n");
  
  x = QScalMult(QAdd(QMult(q,mu1),QMult(mu1,q)),-1./2.);
  y = QScalMult(QAdd(QMult(q,mu2),QMult(mu2,q)),-1./2.);
  z = QScalMult(QAdd(QMult(q,mu3),QMult(mu3,q)),-1./2.);
  printf("x: ");QDisp(x,stdout);printf("\n");
  printf("y: ");QDisp(y,stdout);printf("\n");
  printf("z: ");QDisp(z,stdout);printf("\n"); 
}

void changement_base()
{
  Quaternion q1,q2;
  Quaternion mu1,mu2,mu3;
  double e1,e2,e3;
  
  q1=QInit(1.,5.,2.,3.);
  QDisp((const Quaternion) q1,stdout);
  printf("norme de q : %lf\n",QNorm((const Quaternion)q1));
  mu1=QInit(0.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.));
  //mu2=mu1=QInit(0.,-1./sqrt(2.),1./sqrt(2.),0.);
  mu2 = QInit(0.,0.,1./sqrt(2.),-1./sqrt(2.));

  mu3=QMult(mu1,mu2);
  printf("\n");
  QDisp((const Quaternion) mu1,stdout);
  printf("norme de mu1 : %lf\n",QNorm((const Quaternion)mu1));
  printf("\n");
  QDisp((const Quaternion) mu2,stdout);
  printf("norme de mu2 : %lf\n",QNorm((const Quaternion)mu2));
  printf("\n");
  QDisp((const Quaternion) mu3,stdout);
  printf("norme de mu3q : %lf\n",QNorm((const Quaternion)mu3));
  printf("\n");
  
  
  e1= (-1./2.)*(QAdd(QMult(q1,mu1),QMult(mu1,q1))).a;
  printf("%lf\n",e1);
  e2= (-1./2.)*(QAdd(QMult(q1,mu2),QMult(mu2,q1))).a;
  printf("%lf\n",e2);
  e3= (-1./2.)*(QAdd(QMult(q1,mu3),QMult(mu3,q1))).a;
  printf("%lf\n",e3);
  q2=QInit(1.,e1,e2,e3);
  QDisp(q2,stdout);
  printf("norme de q2 : %lf\n",QNorm((const Quaternion)q2));
  printf("\n");
  
}

int main_quaternion(int argc, char *argv[])
{
    //test_operations_base();
    
    //produit_quaternionique(argc,argv);
    //calcul();
    //changement_base();
    
    return 0;
}

