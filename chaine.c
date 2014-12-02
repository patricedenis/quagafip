#include "chaine.h"

/*//substring between file extention
void sbstr_btw_filext(const char * str,char * between, char ** result)
{
    char* p;
    char* buffer;
    //on cherche a ajouter la chaine "between" avant le "." de l'extention du fichier
    char* separateurs = {"."};
    
    buffer = malloc(strlen(str)*sizeof(char));
    
    strcpy(buffer,str);
    //premier appel, renvoie le premier mot : le début du nom de fichier
    p = strtok(buffer,separateurs);
    strcpy(buffer,p);
    strcpy(*result,buffer);
    strcat(*result,between);
    strcat(*result,separateurs);
    //on récupère l'extention du nom de fichier
    p = strtok(NULL,separateurs);
    strcat(*result,p);
    free(buffer);
}*/

//retourne la position du caractère chrSep dans la chaine strSource
//ou -1 si celui-ci n'est pas présent dans la chaine
int strchrind(char * strSource,char chrSep)
{
  int i;
  //on parcourt la chaine jusqu'a trouver le separateur
  for(i=0;i<strlen(strSource);i++)
    if(strSource[i]==chrSep)
      return i; //pour renvoyer l'indice
  return -1; // ou -1 quand pas de separateur present dans la chaine
}


// permet d'insérer la chaine pointée par "strBetween" juste avant le premier séparateur "chrSep"
// dans la chaine strSource. Le tout est pointé par "pstrDestination" une autre chaine
// qui aura du être allouée avant d'être appelée dans cette fonction.
void strcpybtw(char ** pstrDestination, char * strSource,char * strBetween,char chrSep)
{
  char  *fin;
  int indice;
   
  //on récupère la fin de la chaine
  fin = strchr (strSource,chrSep);
  //l'indice du séparateur dans la chaine source
  indice = strchrind(strSource,chrSep);
  //construction de la chaine résultat
  *pstrDestination[0]='\0';
  strncat (*pstrDestination,strSource,indice);
  strcat(*pstrDestination,strBetween);
  strcat(*pstrDestination,fin);
}






