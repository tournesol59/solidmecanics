// Simple test mixed C-fortran
#include <stdio.h>
#include <stdlib.h>

typedef struct INFO {
        double x[10];
} *tINFO;

typedef struct TABLE {
        double* x;
} *tTABLE;

extern void subtest(void* ptr1);

extern void subcopy(void* ptr1);

int main() {
 
// void* pt1;
 void* pt1;
 double v=5.0;
//TABLE pt5;
 struct INFO info2;


 pt1=NULL;
// pt5->x=&v;

//  subtest(pt1); // fortran procedure
/*
 if (pt1 != NULL) {
    printf("OK smth for 1st pointer to alloc fileds!\n");
  }  
  else printf("NOK, 1st pointer NULL\n");
*/  
//  pt5->x = (double*) pt1;

/////////////////// 1rst part : //////////////////////

  info2.x[0]=3.14;
  pt1 = (void*) &(info2.x[0]);
  printf("BEFORE FORTRAN SUB test struct content pointed by pt: %f \n",*((double*)pt1));

  subtest(pt1);
  
  if (pt1 != NULL) {
    printf("OK smth for 1st pointer to alloc AFTER FORTRAN SUB: %f \n",*((double*)pt1)); //: 0.000000 always whichever fortran value is affected
  }  
  else printf("NOK, 1st pointer NULL\n");
 

 return 0;
}

