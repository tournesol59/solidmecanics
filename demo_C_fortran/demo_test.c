#include <stdio.h>
#include <malloc.h>
#include "demo_cinterface.h"


int main() {

  //!--------------------------------------------------------------------------! 
  //! Variablendeklarationen                                                   ! 
  //!--------------------------------------------------------------------------! 

     tMESHINFO*                  NTestGitter; // Test Gitter Dimensions for std input !
     tMESH*                      Gitter;     // Gitterwerte in matriz           ! 
     tMESHGEN*                   Gitter2D;   // Gitterwerte in Generic inputfile! 
     tCURV*                      KurvVect;   // Numerische Felder fuer Kurve, Christoffei Koeffs abhangig, Rechnung ! 
     tRANDBEDINGUNGEN*           RB;         // Randbedingungen Code            ! 
     tCONSTANTS*                 Const;      // Konstanten                      ! 
     tFILEIO*                    FileIO;     // Ausgabesteuerung                ! 
     tEXACT*                     Exakt;      // Exakte Loesung                  ! 
     Number                     *Uvar;      // Feld mit der numerischen Loesung!!=Tr(spanng) ! 
     Number                     *rhs;       // Feld für die Rechte Seite       ! 
     tNUMERIC*                   VarNum;     // Numerische Felder mit Loesungen ! 
     //!--------------------------------------------------------------------------!    
//     Number                     **SkalarProdMatrixQ;   // full dense Matrix fuer Erfassung des 2. Sys. Gleichung, 1-Dim 
//     Number                     ***SkalarProdMatrixT;  // full dense Matrix fuer Erfassung des 1. Sys. Gleichung, 4-Dims 
    //int errn1, errn2;

//     Number                    *SparseSkalarProdRow1;  // sparse matrix csr format: three vectors (1)
//! size=c.a. 500
//     Number                    *SparseSkalarProdCol;   // sparse matrix csr format: three vectors (2)
//     Number                    *SparseSkalarProdVal;   // sparse matrix csr format: three vectors (3)
     
     printf("Einlesen der Werte:\n");
 
     NTestGitter = (tMESHINFO*) malloc (sizeof(tMESHINFO));
     Const = (tCONSTANTS*) malloc (sizeof(tCONSTANTS));

     // --- Rechengebiet: 
     scanf( "%f", &NTestGitter->startx );
     scanf( "%f", &NTestGitter->starty ); 
     scanf( "%f", &NTestGitter->endx ); 
     scanf( "%f", &NTestGitter->endy ); 
     scanf( "%d", &NTestGitter->nraumx );
     scanf( "%d", &NTestGitter->nraumy );
 //    scanf( "%f", &Const->EY ); 
 //    scanf( "%f", &Const->nu ); 

     Const->EY = 270000.0;
     Const->nu = 0.3;
     Const->h = 0.15 ;  // Could not explain why read fails for EY, nu, h 
 
    // ----------------------------------<  Kontroll-Ausgabe der Groessen >-- 
     printf("\n");    
     printf(" Rechengebiet: \n"); 
     printf("    startx   = %f \n" ,NTestGitter->startx );
     printf("    starty   = %f \n" ,NTestGitter->starty );
     printf("    endx     = %f \n" ,NTestGitter->endx );
     printf("    endy     = %f \n" ,NTestGitter->endy );
 
     printf("    nraumx   = %d \n" ,NTestGitter->nraumx );
     printf("    nraumy   = %d \n" ,NTestGitter->nraumy );
     printf("    platten material parameter \n");  
     printf("    h        = %f \n" ,Const->h ); 
     printf("    EY       = %f \n" ,Const->EY );
     printf("    nu       = %f \n" ,Const->nu );
     printf("\n");
 
     printf("confirm (1) or press (0) to choose a predefined curved square meshregion\n"); 
     scanf("%d ", &Const->automesh); 
     printf("    You have chosen   = %d \n" ,Const->automesh );
//     __allocatefields_mod_MOD_allocatefields(Const,Gitter,RB,Uvar,rhs,Exakt);
//     if (errn1 != 0) printf("Error allocate \n");
//     printf("Success allocate \n");
     __allocatefields_mod_MOD_allocatefieldscurv(Gitter,KurvVect);
//     if (errn1 != 0) printf("Error allocate Kurvatur \n");     
     printf("Success allocate KurvVect \n");
     __allocatefields_mod_MOD_allocatefieldsvarnum(Gitter,VarNum);
//     if (errn1 != 0) printf("Error allocate Numeric \n");
     printf("Success allocate Numeric \n");
//     __lagrangegrad2d_mod_MOD_linlgallocatesparse(Gitter, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal);
//     if (errn1 != 0) printf("Error allocate ProdMatrices \n");

     if (Const->automesh == 1) {
  //       errn2=createMesh_(Gitter);
  //       if (errn2 != 0) printf("Error create Mesh \n");         
         }
     else {
         printf("Error no procedure available for Generic Mesh \n");
         
     }
/*
    __interpol2d_mod_MOD_basisdifferenz(Gitter,KurvVect->der_a1xx_x,KurvVect->der_a1xx_y,
                                KurvVect->der_a1xy_x,KurvVect->der_a1xy_y, 
                                KurvVect->der_a2yy_x, KurvVect->der_a2yy_y, 
                                KurvVect->der_a2yx_x, KurvVect->der_a2yx_y); 

    __interpol2d_mod_MOD_christoffei(Gitter,KurvVect->chicoeff, 
                                KurvVect->der_a1xx_x,KurvVect->der_a1xx_y,
                                KurvVect->der_a1xy_x,KurvVect->der_a1xy_y, 
                                KurvVect->der_a2yy_x, KurvVect->der_a2yy_y, 
                                KurvVect->der_a2yx_x, KurvVect->der_a2yx_y); 
*/

//    __allocatefields_mod_MOD_deallocatefields(Const,Gitter,RB,Exakt,Uvar,rhs) ;
    __allocatefields_mod_MOD_deallocatefieldscurv(KurvVect) ;
    __allocatefields_mod_MOD_deallocatefieldsvarnum(VarNum) ;
    printf("Success END of the program \n");
   return 0;
}
