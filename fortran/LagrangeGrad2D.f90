module LagrangeGrad2D_mod

! Use: in Gesetz(1) f_ = Sum,k-i mu,k-i*Phi,i(x,y)

implicit none
private

!  interface
!    module procedure linlgallocatesparse
!  end interface

!  interface
!    module procedure linlgfillmatrix
!  end interface

!  interface
!    module procedure linlgfillbiegung
!  end interface

public :: linlgallocatesparse, linlgfillmatrix, linlgfillbiegung

contains 

!****************************************************************************************!
!  procedure linlgallocatesparse                                                         !
!     csc sparse matrix format compressed rows consists in  :                            !
!     + table of Reduced Row Indexes and number of nonzeros-Element                      !
!     + table of Col Indexes of non-zero Elements
!     + table of values  (real matrix nonzero entries on line i are on columns Col[Row[i],Row[i+1]-1]  !
!        and their values are the matching entries in the vector values.
!                                                                                        !
!****************************************************************************************!
 SUBROUTINE linlgallocatesparse(Gitter, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal)
 use types
 use interpol2D_mod
 implicit none
  type(tMesh)         :: Gitter
  integer             :: allocStat  
  integer,dimension(4) :: SparseDiagLength=(/ 99, 99, 90, 90 /)
  real,pointer   :: SparseSkalarProdRow1(:)
  real,pointer   :: SparseSkalarProdCol(:)
  real,pointer   :: SparseSkalarProdVal(:)
  intent(in)          :: Gitter
  intent(inout)       :: SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal

 allocate(SparseSkalarProdRow1(1:(Gitter%nraumx*Gitter%nraumy)), &
          SparseSkalarProdCol(1:(Gitter%nraumx*Gitter%nraumy)), &
          SparseSkalarProdVal(1:(Gitter%nraumx*Gitter%nraumy)+SparseDiagLength(1)+SparseDiagLength(2)+ &
                      SparseDiagLength(3)+SparseDiagLength(4)), &
         STAT = allocStat)
 END SUBROUTINE linlgallocatesparse

!*************************!
! and deallocate in main()!
!*************************!


!********************************************************************!
!  procedure linlgfillmatrix: fill matrix mit skalar Produkten von   !
!     (2x2) Gradienten einfach                                       !
!     wird benutzt werden, um alle Termen der Variation-Formula zu bestimmen !
!********************************************************************!
 SUBROUTINE linlgfillmatrix(Gitter, Const, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal)
 use types
 use interpol2D_mod
 implicit none
 type(tMesh)             :: Gitter
 type(tConstants)        :: Const
  real,pointer   :: SparseSkalarProdRow1(:)
  real,pointer   :: SparseSkalarProdCol(:)
  real,pointer   :: SparseSkalarProdVal(:)
  intent(in)          :: Gitter, Const
  intent(inout)       :: SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal

 ! lokal Variablen !
 integer                 :: rowcount, colcount, i
 
  rowcount=0
  colcount=0

! Limited to double gradient calculation, the matrix of the finite elements has size
! (Gitter%nraumx*Gitter%nraumy)^2, which is 10000 if disretization is 10 in x and y
! thus a sparse representation (CSR Compressed Sparse Row, see in Wikipedia) of the matrix
! is taken.
! Remind that the full matrix shall have this form:
! X  y  0  .  .  z  0  .  .  .  .  .  0   ! (0)
! y  X  y  0  .  .  z  0  .  .  .  .  .   ! (1)
! 0  y  X  y  0  .  .  z  0  .  .  .  .   ! (2)
! 0  0  y  X  y  0  .  .  z  0  .  .  .   ! (3) 
! 0  .  0  y  X  y  0  .  .  z  0  .  .   ! (4)
! z  .  .  0  y  X  y  0  .  .  z  0  .   ! (5)
! 0  z  .  .  0  y  X  y  0  .  .  z  0   ! (6)
! .  0  z  .  .  0  y  X  y  0  .  .  z   ! (7)
! .  .  0  z  0  .  .  y  X  y  0  .  0   ! (8)
! .  .  .  0  z  0  .  .  y  X  y  0  0   ! (9)
! .  .  .  .  0  z  0  .  .  y  X  y  0   ! (10)
! .  .  .  .  .  0  z  0  .  .  y  X  y   ! (11)
! 0  .  .  .  .  .  0  z  .  .  0  y  X   ! (12)

! where X,y,z shall be calculated from tabulated values in the global variable GradientScalarProd
! but here, we will calculate it as constants
!   X= 16 = skalar Produkt of an hexagonal 2D Element with itself = 2*(1+1+2+1+1+2)= 16 !
!   y= -2 = skalar Product of an hexagonal 2D Element with a +1 (x or y direction) decay of itself = -1+(-1)
!   z= -2 = skalar Product of an hexagonal 2D Element with a +Nraum (x or y direction) decay of itself = -1+(-1)

  do while ((rowcount).le.(Gitter%nraumx*Gitter%nraumy-1))
     if ((rowcount).le.(10-1)) then  !(Gitter%nraumx-1))               ! fill lines (0..4) in the illustration
         if (rowcount.le.0) then                            ! fill line (0)
             SparseSkalarProdCol(colcount+1)=0
             SparseSkalarProdVal(colcount+1)=16.0
             colcount=colcount+1
             SparseSkalarProdCol(colcount+1)=1
             SparseSkalarProdVal(colcount+1)=-2.0
             colcount=colcount+1
             SparseSkalarProdCol(colcount+1)=10
             SparseSkalarProdVal(colcount+10)=-2.0
             colcount=colcount+1
        else                                          ! fill other lines (1..4)
             SparseSkalarProdCol(colcount+1)=rowcount-1
             SparseSkalarProdVal(colcount+1)=-2.0
             colcount=colcount+1
             SparseSkalarProdCol(colcount+1)=rowcount
             SparseSkalarProdVal(colcount+1)=16.0
             colcount=colcount+1
             SparseSkalarProdCol(colcount+1)=rowcount+1
             SparseSkalarProdVal(colcount+1)=-2.0
             colcount=colcount+1
             SparseSkalarProdCol(colcount+1)=rowcount+10
             SparseSkalarProdVal(colcount+10)=-2.0
             colcount=colcount+1
        endif
     elseif ((rowcount).le.(Gitter%nraumx*(Gitter%nraumy-1)-1)) then  ! fill lines (5..7) in the illustration
             SparseSkalarProdCol(colcount+1)= rowcount-10
             SparseSkalarProdVal(colcount+1)=-2.0
             colcount=colcount+1
             SparseSkalarProdCol(colcount+1)=rowcount-1
             SparseSkalarProdVal(colcount+10)=-2.0
             colcount=colcount+1
             SparseSkalarProdCol(colcount+1)=rowcount
             SparseSkalarProdVal(colcount+1)=16.0
             colcount=colcount+1
             SparseSkalarProdCol(colcount+1)=rowcount+1
             SparseSkalarProdVal(colcount+1)=-2.0
             colcount=colcount+1
             SparseSkalarProdCol(colcount+1)=rowcount+10
             SparseSkalarProdVal(colcount+10)=-2.0
             colcount=colcount+1
     else                                                 ! fill lines (7..12) in the illustration
             SparseSkalarProdCol(colcount+1)= rowcount-10
             SparseSkalarProdVal(colcount+1)=-2.0
             colcount=colcount+1
             SparseSkalarProdCol(colcount+1)=rowcount-1
             SparseSkalarProdVal(colcount+10)=-2.0
             colcount=colcount+1
             SparseSkalarProdCol(colcount+1)=rowcount
             if ((rowcount).eq.(Gitter%nraumx*Gitter%nraumy-1)) then  ! fill last line (12)
                 SparseSkalarProdVal(colcount+1)=16.0
                 colcount=colcount+1
             endif
     endif
     SparseSkalarProdRow1(rowcount)=colcount
     rowcount=rowcount+1
  end do
  !SparseSkalarProdRow1(rowcount)=colcount ! normally already done by the algorithm: test it!

 END SUBROUTINE linlgfillmatrix

!********************************************************************!
!  procedure linlgfillmatrix: fill matrix mit skalar Produkten von   !
!     (2x2) Gradienten einfach                                       !
!     wird benutzt werden, um alle Termen der Variation-Formula zu bestimmen !
!********************************************************************!
! SUBROUTINE linlgfillmatrix(Gitter, Const, VirtSkalarProdM)
! use types
! use interpol2D_mod
! implicit none
!  type(tMesh)           :: Gitter
!  type(tConstants)      :: Const
!  real,pointer          :: VirtSkalarProdM(:,:)   ! (nraumx*nraumy, nraumx*nraumy) muss spaeter als Diagonal BlockMatrix !implementiert werden !
!  ! lokal Variablen  !
!  real                  :: produkt
!  integer               :: i,j,k, m,n,l, p
!
!  do i=0,Gitter%nraumx-1
!    do j=0,Gitter%nraumy-1
!      k=j*Gitter%nraumx+i+1
!
!      do m=0,Gitter%nraumx-1
!        do n=0,Gitter%nraumy-1
!          l=n*Gitter%nraumx+m+1
!          p= l-k
!          ! fill Matrix zero by Default
!          VirtSkalarProdM(k,l)=0.0    
!          select case(p) 
!          case (0)
!             VirtSkalarProdM(k,l)=16.0  ! berechnet von Hand mit Tabelle der Skalar Produkte
!          case (1)
!             VirtSkalarProdM(k,l)=-2.0  ! berechnet
!          case (-1)
!             VirtSkalarProdM(k,l)=-2.0  ! berechnet  
!          case (10)   !Gitter%nraumx ist nicht konstant aber fuer Praesent
!             VirtSkalarProdM(k,l)=-2.0  ! berechnet
!          case (-10)   !-Gitter%nraumx ist nicht konstant aber fuer Praesent
!             VirtSkalarProdM(k,l)=-2.0  ! berechnet
!          end select
!! end, normalws. die berechnete Werte wuerden besser detailliert !
!        enddo
!      enddo
!    enddo
!  enddo
! END SUBROUTINE linlgfillmatrix

!********************************************************************!
!  procedure linlgfillbiegung                                        !
!     wird benutzt werden, um Termen der Biegung in der Variation-Formula zu bestimmen 
!      int,S{ curvature cdot (n vdot BiegMoment) }                   !
!********************************************************************!
 SUBROUTINE linlgfillbiegung(Gitter, Const, VirtSkalarProdM)
 use types
 use interpol2D_mod
 implicit none
  type(tMesh)           :: Gitter
  type(tConstants)      :: Const
  real,pointer          :: VirtSkalarProdM(:,:)   ! (nraumx*nraumy, nraumx*nraumy) muss spaeter als Diagonal BlockMatrix implementiert werden !
  ! lokal Variablen  !
!  real                  :: biegung_1, gradbiegung_1, 
  integer               :: i,j,k, m,n,l, p
  do i=0,Gitter%nraumx-1
    do j=0,Gitter%nraumy-1
      k=j*Gitter%nraumx+i+1

      do m=0,Gitter%nraumx-1
        do n=0,Gitter%nraumy-1
          l=n*Gitter%nraumx+m+1
          p= l-k
          ! search the Kurvatures        
          select case(p) 
          case (0)
             VirtSkalarProdM(k,l)=16.0  ! berechnet 
          case (1)
             VirtSkalarProdM(k,l)=-2.0  ! berechnet
          case (-1)
             VirtSkalarProdM(k,l)=-2.0  ! berechnet  
          case (10)   !Gitter%nraumx ist nicht konstant aber fuer Praesent
             VirtSkalarProdM(k,l)=-2.0  ! berechnet
          case (-10)   !-Gitter%nraumx ist nicht konstant aber fuer Praesent
             VirtSkalarProdM(k,l)=-2.0  ! berechnet
          end select
        enddo
      enddo
    enddo
  enddo
 END SUBROUTINE linlgfillbiegung


end module LagrangeGrad2D_mod
