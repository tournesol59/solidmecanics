/* Subroutined swept from Fortran */

#ifndef DEMO_CINTERFACE_H
#define DEMO_CINTERFACE_H
#endif

#include "demo_ctypes.h"

/********************************************************
   Fortran module allocateFields.f90      
*********************************************************/
void __allocatefields_mod_MOD_allocatefields(tCONSTANTS *Const,tMESH *Mesh,tRANDBEDINGUNGEN *RB,
                   Number *Uvar, Number *rhs, tEXACT *Exakt);

void __allocatefields_mod_MOD_deallocatefields(tCONSTANTS *Const, tMESH *Mesh, tRANDBEDINGUNGEN *RB,
                  tEXACT *Exakt,Number *Uvar,Number *rhs);

void __allocatefields_mod_MOD_allocatefieldscurv(tMESH *Mesh, tCURV *CurvStruct);

void __allocatefields_mod_MOD_deallocatefieldscurv(tCURV *CurvStruct);

void __allocatefields_mod_MOD_allocatefieldsvarnum(tMESH *Mesh, tNUMERIC *VarNum);

void __allocatefields_mod_MOD_deallocatefieldsvarnum(tNUMERIC *VarNum);

void __allocatefields_mod_MOD_allocategitter2d(tMESHGEN *Mesh2D); 

void __allocatefields_mod_MOD_deallocategitter2d(tMESHGEN *Mesh2D);


/********************************************************
   Fortran module createMesh.f90      
*********************************************************/

void __createmesh_mod_MOD_createmesh(tMESH *Mesh);
void __createmesh_mod_MOD_createmesh2(tMESH *Mesh);
void __createmesh_mod_MOD_createmeshgen(tMESHGEN *Mesh2D);
void __createmesh_mod_MOD_createmeshgen2(tMESHGEN *Mesh2D);

/********************************************************
   Fortran module Input_mod.f90      
*********************************************************/

void __input_mod_MOD_input(tMESH *Mesh, tEXACT *Exakt, tCONSTANTS *Const, tFILEIO *FileIO);
void __input_mod_MOD_input_file(tMESHGEN *Mesh2D,tEXACT *Exakt, tCONSTANTS *Const, tFILEIO *FileIO);

/********************************************************
   Fortran module Interpol2D.f90      
*********************************************************/

void __interpol2d_mod_MOD_basisdifferenz(tMESH *Gitter,Number *der_a1xx_x,Number *der_a1xx_y,Number *der_a1xy_x,Number *der_a1xy_y, Number *der_a2yy_x,Number *der_a2yy_y, Number *der_a2yx_x,Number *der_a2yx_y);

void __interpol2d_mod_MOD_christoffei(tMESH *Gitter, Number **chicoeff, Number *der_a1xx_x,Number *der_a1xx_y, Number *der_a1xy_x,Number *der_a1xy_y, Number *der_a2yy_x,Number *der_a2yy_y, Number *der_a2yx_x,Number *der_a2yx_y);

void __interpol2d_mod_MOD_trychristoffei(tMESH *Gitter, Number **chicoeff);

/********************************************************
   Fortran module Lagrange1D.f90      
*********************************************************/

void __lagrange1d_mod_MOD_evaluation(tPOLYNOM *p, float *x, Number *res);

void __lagrange1d_mod_MOD_expand2nom(tPOLYNOM *pol_1,tPOLYNOM *pol);

void __lagrange1d_mod_MOD_calcderivative(tPOLYNOM pol, Number *x, Number *der);

void __lagrange1d_mod_MOD_expandfromroots(tPOLYNOM *pol_roots,tPOLYNOM *pol_nom);

void __lagrange1d_mod_MOD_expandlagrange(tPOLYNOM *pol_pts,tPOLYNOM *pol_value,tPOLYNOM *pol_lagr);

void __lagrange1d_mod_MOD_expandsquare(tPOLYNOM *pol_pts, int *i, tPOLYNOM *polsq);

void __lagrange1d_mod_MOD_squarecalcderivative(tPOLYNOM *pol_pts, tPOLYNOM *pol_roots, int *i, Number *der);


/********************************************************
   Fortran module Hermite1D.f90      
*********************************************************/

void __hermite1d_mod_MOD_expandhermite(tPOLYNOM *pol_pts, tPOLYNOM *pol_value, tPOLYNOM *pol_lagr);

/********************************************************
   Fortran module Hermite2D.f90      
*********************************************************/

void __hermite2d_mod_MOD_linhm2allocatesdns(tMESH *Gitter, Number **SkalarProdMatrixQ, Number **SkalarProdMatrixT);

void __hermite2d_mod_MOD_computegaussint(int *tt, int *ii, tPOLYNOM *pol1, tPOLYNOM *pol2, tPOLYNOM *pol3, Number *valend);
    // horreur programmatique

void __hermite2d_mod_MOD_linhm2dfillmatrixq(tMESH *Gitter, Number **chicoeff, tNUMERIC VarNum, Number **SkalarProdMatrixT);
    // semble bien complet mais jamais testé

//void linlg2ddeallocatedns_(tMESH *Gitter, Number **SkalarProdMatrixQ, Number **SkalarProdMatrixT);

/********************************************************
   Fortran module LagrangeGrad2D.f90      
*********************************************************/

void __lagrangegrad2d_mod_MOD_linlgallocatesparse(tMESH *Gitter, Number *SparseSkalarProdRow1, Number *SparseSkalarProdCol, Number *SparseSkalarProdVal);

void __lagrangegrad2d_mod_MOD_linlgfillmatrix(tMESH *Gitter, tCONSTANTS *Const, Number *SparseSkalarProdRow1, Number *SparseSkalarProdCol, Number *SparseSkalarProdVal);

void __lagrangegrad2d_mod_MOD_linlgfillbiegung(tMESH *Gitter, tCONSTANTS *Const, Number *VirtSkalarProdM);

/********************************************************
   Fortran module verifyquad2D.f90      
*********************************************************/

void __verifyquad2d_mod_MOD_checkderivcoord(tMESHGEN *Mesh2D,int *nel,Number *condit);

void __verifyquad2d_mod_MOD_checkalledgesofelmts(tMESHGEN *Mesh2D);
