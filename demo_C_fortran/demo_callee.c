// Simple test mixed C-fortran
#include <stdio.h>
#include <stdlib.h>

typedef struct INFO {
        double x[10];
} *tINFO;

typedef struct TABLE {
        double* x;
} *tTABLE;

void democallee(void * forfortran) {

  int i,j; // n,l

  
  float * forfortran = {3.14, 2.21, 1.72, 1.41, 1., 1.57, 1.10, 0.788, 0.707, 0.5};

}
