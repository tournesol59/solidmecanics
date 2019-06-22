// Simple test mixed C-fortran
#include <stdio.h>
#include <stdlib.h>

typedef struct INFO {
        double x[3];
} *tINFO;

typedef struct TABLE {
        double* x;
} *tTABLE;

extern void subtest(void* ptr1);


int main() {
 
// void* pt1;
 void* pt1;
 double v=5.0;
 tTABLE pt5;
// struct INFO info2;
 
// pt3=(tINFO) &info2;
// pt3->x[0]=1.0;  pt3->x[1]=0.0;  pt3->x[2]=-1.0; 

 pt1=NULL;
 pt5->x=&v;

  subtest(pt1); // fortran procedure
/*
 if (pt1 != NULL) {
    printf("OK smth for 1st pointer to alloc fileds!\n");
  }  
  else printf("NOK, 1st pointer NULL\n");
*/  
//  pt5->x = (double*) pt1;

  if (pt5 != NULL) {
     printf("1st value is: %f\n", pt5->x[0]);
  }

 return 0;
}

