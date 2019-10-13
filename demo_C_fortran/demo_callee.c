// Simple test mixed C-fortran, only for scalar or simple arrays
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

typedef struct INFO {
        double x[1];
} *tINFO;

typedef struct TABLE {
        double* x;
} *tTABLE;

extern struct
{
  int ii,jj,kk;
} ijk_;

void democallee(double * forfortran) {

  int i,j; // n,l
  tINFO myinfo = malloc(sizeof(struct INFO));
  double myscalar=10.0;

  printf("just after calling C\n");
  
 // double mydouble[10] = {3.14, 2.21, 1.72, 1.41, 1., 1.57, 1.10, 0.788, 0.707, 0.5};

  for (i=0; i<1; i++) {
     myinfo->x[i]=((float) i)+1.5;

  }
  printf("after call before return to fortran\n");
 // forfortran = (void*) myptr; // type casting do not work for me
  myscalar=forfortran[0];
  printf("before return to fortran, value copied in C: %f\n", myscalar);
}


int doubleijk_(char *cc, int ll)
{

  cc[ll--]='\0';  // NULL terminate the string
  printf("From doubleIJK: %s \n",cc);
  ijk_.ii *=2;
  ijk_.jj *=2;
  ijk_.kk *=2;

  return (1);
}
