module Transplanar2D
! This module contains the necessary functions to assembly, !
! extend a stiffness matrix from a 2D local basis to anoth  !
! er 2D global basis, assembly two matrices of different el.!
! ements to a global matrix in sparse CSR format. The assum.!
! tions are the plane stress thin plate elements and contin.!
! uous flexion, with 3-nodes integration.                   !
! The format and dimension of the global matrix is the numb.!
! er of nodes of the triangular elements times the size of a!
! local stiffness matrix, (9+6)**2, in the present case. It !
! should also vary with the addition of boundary conditions.!

 implicit none
 private
 
 integer,dimension(9)  :: Indblk = /(2,2,2,2,1,2,1,2,1/)
                   ! corresponds to size of a stiffness     !
                   ! matrix blocks in the plane trig element!
                   !(2,2,2. correspond to the membrane sress!
                   !.2,1,2,1,2,1) to the flexion stress:2 is!
                   ! for the plane rot. (x,y) deg of freedom!
                   !1 is for the Wz displacemnt for 3 points!
                   !sum of each degree=15                   !
 integer               :: SumDeg = 15

 public:: calcblockrotmatrix, transblockstiffmatrix, &
          assemblystiffmatrix

 contains

!****************************************************************************************!
!  procedure assemblystiffmatrix                                                        !
!     csc sparse matrix format compressed rows consists in  :                            !
!     + table of Reduced Row Indexes and number of nonzeros-Element                      !
!     + table of Col Indexes of non-zero Elements
!     + table of values  (real matrix nonzero entries on line i are on columns Col[Row[i],Row[i+1]-1]  !
!        and their values are the matching entries in the vector values.
!                                                                                        !
!****************************************************************************************!
 subroutine assemblystiffmatrix(Gitter2D, elmt, SingleStiffMat, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal)
 use types
 implicit none
  type(tMeshGen) :: Gitter2D
  integer        :: elmt
  real,pointer   :: SingleStiffMat 
  real,pointer   :: SparseSkalarProdRow1(:)
  real,pointer   :: SparseSkalarProdCol(:)
  real,pointer   :: SparseSkalarProdVal(:)
  intent(in)          :: Gitter2D, elmt, SingleStiffmat
  intent(inout)       :: SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal
  
  integer        :: rowcount, colcount, i,j


 end subroutine assemblystiffmatrix

end module Transplanar2D

