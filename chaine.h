#ifndef _CHAINE_H_
#define _CHAINE_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//retourne la position du caractère chrSep dans la chaine strSource
//ou -1 si celui-ci n'est pas présent dans la chaine
int strchrind(char * strSource,char chrSep);

// permet d'insérer la chaine pointée par "strBetween" juste avant le premier séparateur "chrSep"
// dans la chaine strSource. Le tout est pointé par "pstrDestination" une autre chaine
// qui aura du être allouée avant d'être appelée dans cette fonction.
void strcpybtw(char ** pstrDestination, char * strSource,char * strBetween,char chrSep);

#endif

