/*
  Interface: geometry,curvature, numeric, right-hand-side and boundary
  structures that contains the elements of the FE shell program.
  This sweeps the Fortran types.f90 corresponding structures.
*/
//#include <stdio.h>
#include <stdbool.h>

typedef float Number;

//  for rectangles/squares elements-meshs
typedef struct tMESHINFO
   {
   int        nraumx, nraumy;
   Number     startx, starty;
   Number     endx, endy;
} tMESHINFO*;

typedef struct tMESH
   {
   int        nraumx, nraumy;
   Number     startx, starty;
   Number     endx, endy;
   Number     dx,dy,dxq,dyq;
   Number     dxx,dyy,dxxq,dyyq;
   Number     *x,*y,*z;        //  nodes coordinates
   Number     **Coefficients;  //  
   } tMESH;

typedef struct tMESHGEN
   {
   int        nodes;
   Number     *x,*y,*z;
   int        elmts;
   int        *quad;            // indices for rectangle elmt
   int        ntop, nbottom, nleft, nright;
   int        *quadtop, *quadbottom,*quadleft, *quadright;// boundary indices of 
                                    //quadrangles
   } tMESHGEN;

typedef struct tCURV
   {
   Number    **chicoeff;
   Number  *der_a1xx_x,*der_a1xx_y,*der_a1xy_x,*der_a1xy_y,*der_a2yx_x,*der_a2yx_y,*der_a2yy_x,*der_a2yy_y;
   } tCURV;

typedef struct tCONSTANTS
   {
   Number    Pi, EY, nu, h;
   int       automesh;
//  int       controlAus, imax;;
//  char*     sheme, loeser, Abbruch, RHS_Typ, RB_Typ, AW_Typ, RHS_Datei;

   } tCONSTANTS;


typedef struct tEXACT
   {
   Number    *loesung;
   int        kmaxex;
   } tEXACT;

