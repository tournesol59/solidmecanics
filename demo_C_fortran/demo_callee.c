// Simple test mixed C-fortran
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

typedef struct INFO {
        double x[1];
} *tINFO;

typedef struct TABLE {
        double* x;
} *tTABLE;

void democallee_(void * forfortran) {

  int i,j; // n,l
  tINFO myinfo = malloc(sizeof(struct INFO));

  printf("just after calling C\n");
  
 // double mydouble[10] = {3.14, 2.21, 1.72, 1.41, 1., 1.57, 1.10, 0.788, 0.707, 0.5};

  for (i=0; i<1; i++) {
     myinfo->x[i]=i+1.5;

  }
  printf("after call before typecasting and return to fortran\n");
  forfortran = (void*) myinfo;
}
