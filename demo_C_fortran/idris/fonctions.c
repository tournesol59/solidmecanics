#include <stdio.h>
#include <string.h>
#include <malloc.h>

void function(char *chaine,
              int *entier,
              double *reel )
{
  strcat( chaine, "et chaine_c");
  *entier = strlen( chaine );
  *reel = 100.0;
}

typedef struct cel
{
  double r;
  int n;
} Cellule;

/*
 * Fonctions  retournant l adresse d'un objet
 * de type Cellule.
 */

Cellule *creation_cel( double *r, 
                       int *n)
{
  Cellule *c = malloc( sizeof(Cellule) );
  c->r = *r;
  c->n = *n;
  
  return c;
}

void modif_cel(Cellule **c,
               double * r,
               int *n)
{
   (*c)->r = *r;
   (*c)->n= *n;

   return;
}

void imp_cel( Cellule **c )
{
   printf(" Cellule: %f, %d \n", (*c)->r, (*c)->n );
   
   return;
}

void libere_cel( Cellule **c )  
{
  free( *c );

  return;
}

