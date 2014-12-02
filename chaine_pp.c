#include "chaine.h"


void chaine_test1()
{
  char source[20]={"coucou.bmp"};
  char separateur='3';
  char entre[10]={"-machin"};
  char * resultat;
  int i,ind=-1;
  printf("%d\n",strlen(source));
  for(i=0;i<strlen(source);i++)
  {
   if (source[i]==separateur)
    ind = i;
  }
  printf("%d\n",ind);
  printf("strcharind : %d\n",strchrind("coucou.bmp",'.'));
  resultat = (char *)malloc((strlen("coucou.bmp")+strlen("coucou.bmp"))*sizeof(char));
  strcpybtw(&resultat,"coucou.bmp","-angle",'.');
  printf("test concat between : %s\n",resultat);
}


int main_chaine_pp(int argc, char *argv[])
{
  int a;
  char * result;
  char c[20] = {"coucou.bmp"};
  
  chaine_test1();
  return 1;
}


