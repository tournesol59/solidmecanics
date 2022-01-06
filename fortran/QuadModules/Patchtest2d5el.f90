module Patchtest2d5el_mod
!******************************************************!
! The Rigidity Matrix is directly assembled from a     ! 
! 5 elements assumed patch test and not from a grid    !
! mesh treated by a from input data algorithm          !
!******************************************************!
  use types
 ! use createMesh_mod
 
  implicit none

  private
!<  x(4)--------------x(3) >
!< '\               / '    >
!< '  x(8)------x(7)  '    >
!< '  '         '     '    >
!< '  '        '      '    >
!< '  x(5)-----x(6)   '    >
!< ' /             \  '    >
! x(1)-------------x(2)

! or equivalent :

!<  ------------------^B   >
!< '                  '    >
!< '->E\ ^F    |->G\H'    >
!< '  '  '      '  ^  '    >
!< '  '  '     '   '  '    >
!< -->C\^D----' L^ '  '    >
!< '--------->K  '    '    >
!  ------------------>A

! SEE createMesh.f90 to define the coordinates of this mesh
! Assumptions: y(4)=y(3) and x(4)=x(1)

#define PAT5GITTER_A patchXe(2)
#define PAT5GITTER_B patchYe(3)
#define PAT5GITTER_C patchXe(5)
#define PAT5GITTER_D patchYe(5) 
#define PAT5GITTER_E (patchXe(8)-patchXe(5))
#define PAT5GITTER_F (patchYe(8)-patchYe(5))
#define PAT5GITTER_G (patchXe(7)-patchXe(6))
#define PAT5GITTER_H (patchYe(7)-patchYe(6))
#define PAT5GITTER_K patchXe(6)
#define PAT5GITTER_L patchYe(6)

! global variables constants
  real,dimension(8) :: patchXe = (/ 0.0, 100.0, 100.0, 0.0, 25.0, 75.0, 75.0, 25.0 /)
  real,dimension(8) :: patchYe = (/ 0.0, 0.0, 50.0, 50.0, 12.5, 12.5, 37.5, 37.5 /)
  real, dimension(8) :: patchForces = (/ -1.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0 /)

  real,dimension(9) :: Planestress = (/ (/200999., 66333., 0. /), & 
                                              (/ 66333., 200999., 0. /), &
                                              (/ 0. ,0. , 290000.0 /) /)
  real,dimension(2)  :: gaussP = (/-0.57735, 0.57735  /)
  real,dimension(96) :: QuadMatrixNConst = (/ (/ (/ -1.57735, 0.0, 1.57735, 0.0, 0.42265, 0.0, -1.57735, 0.0 /), &
                           (/ 0.0, -1.57735, 0.0, -0.42265, 0.0, 1.57735, 0.0, 1.57735 /), &
          (/ -1.57735, -1.57735, -0.42265, 1.57735, 1.57735, 0.42265, 1.57735,-1.57735 /) /), &
                        (/ (/ -1.57735, 0.0, 1.57735, 0.0, 0.42265, 0.0, -0.42265, 0.0 /), &
                           (/ 0.0, -0.42265, 0.0, -1.57735, 0.0, 1.57735, 0.0, 0.42265 /), &
          (/ -0.42265, -1.57735, -1.57735, 1.57735, 1.57735, 0.42265, 0.42265, 0.42265 /) /), &
                        (/ (/ -0.42265, 0.0, 0.42265, 0.0, 1.57735, 0.0, -1.57735, 0.0 /), &
                           (/ 0.0, -0.42265, 0.0, -1.57735, 0.0, 1.57735, 0.0, 0.42265 /), &
          (/ -0.42265, -0.42265, -1.57735, 0.42265, 1.57735, 1.57735, 0.42265, -1.57735 /) /), &
                        (/ (/ -0.42265, 0.0, 0.42265, 0.0, 1.57735, 0.0, -1.57735, 0.0 /), &
                           (/ 0.0, -1.57735, 0.0, -0.42265, 0.0, 0.42265, 0.0, 1.57735 /), &
          (/ -1.57735, -0.42265, -0.42265, 0.42265, 0.42265, 1.57735, 1.57735, -1.57735 /) /) /) 

! global variable definitions in this package
  real,pointer            :: matrixJ(:,:)  !dimension(1:4,1:2)
  real,pointer            :: matrixX(:,:)  ! dimension(1:2,1:4) 
  real,pointer            :: matrixInvJ(:,:)! dimension(1:4,1:2)
  real,pointer            :: matrixN(:,:)  ! dimension(1:8,1:3)
  real,pointer            :: matrixB(:,:)  ! dimension(1:3,1:8)
  real,pointer            :: matrixKimd1(:,:)  ! dimension(1:3,1:8)
  real,pointer            :: matrixKtrB(:,:)  ! dimension(1:8,1:3)
  real,pointer            :: matrixKinc(:,:)  ! dimension(1:8,1:8)*
  real,pointer            :: matrixKsum(:,:)  ! dimension(1:8,1:8)


  public :: linlgallocmatrixpatch5el, linlgfillmatrixXpatch5elX, linlgfillmatrixJpatch5elJ, & 
             linlgfillmatrixNpatch5elN, linlgintegmatrixpatch5el,linlgprodmatrixBpatch5elB, & 
             linlgdeallocmatrixpatch5el

contains
!****************************************!
! linlgallocmatrixpatch5el               !
!
  subroutine linlgallocmatrixpatch5el(Mesh2D, RB, Uvar, rhs, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal)
 use types
 implicit none
  type(tMeshGen)         :: Mesh2D
  type(tRandbedingungen) :: RB
  real,pointer        :: Uvar(:)
  real,pointer        :: rhs(:)
  integer             :: allocStat  
  integer,dimension(4) :: SparseDiagLength=(/ 8, 8, 8, 8 /)
  real,pointer        :: SparseSkalarProdRow1(:)
  real,pointer        :: SparseSkalarProdCol(:)
  real,pointer        :: SparseSkalarProdVal(:)
  intent(in)          :: Mesh2D
  intent(inout)       :: SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal

     integer          :: randltype      ! RandTyp links: 1=Fixed, 2=Momentum, 3=Freibiegung
     real,pointer     :: randl(:)       ! Randwerte links

! rando(1:4) because 2xvalues (Fx,Fy) for each point (2) in sides left, above, right, bottom
  allocate(RB%randl(1:4), RB%rando(1:4), RB%randr(1:4), RB%randu(1:4), &
          Uvar(1:8*2), rhs(1:8*2), &
          SparseSkalarProdRow1(1:8), &
          SparseSkalarProdCol(1:8), &
          SparseSkalarProdVal(1:(8+SparseDiagLength(1)+SparseDiagLength(2)+ &
                      SparseDiagLength(3)+SparseDiagLength(4)) ) , STAT = allocStat)

  allocate(matrixJ(1:4,1:2), matrixX(1:2,1:4), matrixInvJ(1:4,1:2), & 
           matrixN(1:8,1:3), matrixKinc(1:8,1:8), matrixKsum(1:8,1:8), & 
           matrixB(1:3,1:8), matrixKimd1(1:3,1:8), matrixKtrB(1:8,1:3), STAT=allocStat)

  end subroutine linlgallocmatrixpatch5el

!****************************************
! linlgcreatestructpatch5el
!
  subroutine linlgcreatestructpatch5el(Mesh2D, RB)
  use types
  implicit none
  type(tMeshGen)         :: Mesh2D
  type(tRandbedingungen) :: RB
    ! Timewatch revolution sense
  RB%randltype = 0  ! tensile force
  RB%randl(1)  = patchForces(1)
  RB%randl(2)  = patchForces(2)
  RB%randl(3)  = patchForces(7)
  RB%randl(4)  = patchForces(8)

  RB%rando(1)  = patchForces(7)
  RB%rando(2)  = patchForces(8)
  RB%rando(3)  = patchForces(5)
  RB%rando(4)  = patchForces(6)

  RB%randr(1)  = patchForces(5)
  RB%randr(2)  = patchForces(5)
  RB%randr(3)  = patchForces(3)
  RB%randr(4)  = patchForces(4)

  RB%randu(1)  = patchForces(3)
  RB%randu(2)  = patchForces(4)
  RB%randu(3)  = patchForces(1)
  RB%randu(1)  = patchForces(2)
  end subroutine linlgcreatestructpatch5el

!****************************************!
! linlgfillmatrixXpatch5elX                !
!
  subroutine linlgfillmatrixXpatch5elX(Mesh2D,ne)
  use types
  implicit none

  type(tMeshGen)        :: Mesh2D
  integer               :: ne
  integer               :: i,j,k

  select case(ne) ! fill matrix coord: fortran convention column, 
                  ! case of 5 elements of the form shown above
    case(1)  ! the element at the bottom
      matrixX(1,1)= 0.
      matrixX(2,1)= 0.
      matrixX(1,2)= PAT5GITTER_A 
      matrixX(2,2)= 0.
      matrixX(1,3)= PAT5GITTER_K
      matrixX(2,3)= PAT5GITTER_L
      matrixX(1,4)= PAT5GITTER_C
      matrixX(2,4)= PAT5GITTER_D
    case(2)  ! the element at the right
      matrixX(1,1)= PAT5GITTER_K
      matrixX(2,1)= PAT5GITTER_L
      matrixX(1,2)= PAT5GITTER_A
      matrixX(2,2)= 0.
      matrixX(1,3)= PAT5GITTER_A
      matrixX(2,3)= PAT5GITTER_B
      matrixX(1,4)= PAT5GITTER_K+PAT5GITTER_G
      matrixX(2,4)= PAT5GITTER_L+PAT5GITTER_H

    case(3)  ! the element at the top
      matrixX(1,1)=  PAT5GITTER_C + PAT5GITTER_E 
      matrixX(2,1)=  PAT5GITTER_D + PAT5GITTER_F
      matrixX(1,2)=  PAT5GITTER_K + PAT5GITTER_G
      matrixX(2,2)=  PAT5GITTER_L + PAT5GITTER_H
      matrixX(1,3)= PAT5GITTER_A
      matrixX(2,3)= PAT5GITTER_B
      matrixX(1,4)= 0.
      matrixX(2,4)= PAT5GITTER_B

    case(4)  ! the element at the left
      matrixX(1,1)= 0.
      matrixX(2,1)= 0.
      matrixX(1,2)= PAT5GITTER_C 
      matrixX(2,2)= PAT5GITTER_D
      matrixX(1,3)=  PAT5GITTER_C + PAT5GITTER_E 
      matrixX(2,3)=  PAT5GITTER_D + PAT5GITTER_F
      matrixX(1,4)= 0.
      matrixX(2,4)= PAT5GITTER_B

    case(5)  ! the center element
      matrixX(1,1)= PAT5GITTER_C 
      matrixX(2,1)= PAT5GITTER_D
      matrixX(1,2)= PAT5GITTER_K
      matrixX(2,2)= PAT5GITTER_L
      matrixX(1,3)=  PAT5GITTER_K + PAT5GITTER_G
      matrixX(2,3)=  PAT5GITTER_L + PAT5GITTER_H
      matrixX(1,4)=  PAT5GITTER_C + PAT5GITTER_E 
      matrixX(2,4)=  PAT5GITTER_D + PAT5GITTER_F

  end select

  end subroutine linlgfillmatrixXpatch5elX

!****************************************!
! linlgfillmatrixJinvpatch5elJ        !
! fill matrix of the determinant of the geometry by reference 
! to an ideal square and the inverse to have deriv(idealxy)/d(realX,Y)
  subroutine linlgfillmatrixJpatch5elJ(Mesh2D,ne,ni)
  use types
  implicit none

  type(tMeshGen)        ::Mesh2D
  integer               :: ne,ni
  real                  :: detJ          
  integer i,j,k

  select case(ne)
    case(1)
      select case(ni)
        case(1)
          detJ=1/16.0*(3*PAT5GITTER_A*PAT5GITTER_B/4.0)*(1+gaussP(1))
          matrixInvJ(1,1)=(1/2.0*PAT5GITTER_B)/detJ
          matrixInvJ(1,2)=(0)/detJ
          matrixInvJ(2,1)=gaussP(1)*PAT5GITTER_A/2.0/detJ
          matrixInvJ(2,2)=3.0*PAT5GITTER_A/2.0/detJ
        case(2)
          detJ=1/16.0*(3*PAT5GITTER_A*PAT5GITTER_B/4.0)*(1+1/3.0*gaussP(2))
          matrixInvJ(1,1)=(1/2.0*PAT5GITTER_B)/detJ  ! Simplification _L=_D=_B/4
          matrixInvJ(1,2)=(PAT5GITTER_B/2.0)/detJ
          matrixInvJ(2,1)=PAT5GITTER_A/2.0*gaussP(2)/detJ ! Simplification _K=_A*3/4, _C=_A/4
          matrixInvJ(2,2)=3.0/2.0*PAT5GITTER_A+1/2.0*PAT5GITTER_A*gaussP(2)/detJ
        case(3)
          detJ=1/16.0*(3*PAT5GITTER_A*PAT5GITTER_B/4.0)*(1-1/3.0*gaussP(2))   
          matrixInvJ(1,1)=(PAT5GITTER_B/2.0)/detJ
          matrixInvJ(1,2)=(0.0)/detJ
          matrixInvJ(2,1)=PAT5GITTER_A/2.0*gaussP(2)/detJ
          matrixInvJ(2,2)=3.0/2.0*PAT5GITTER_A*(1-1/3.0*gaussP(2))/detJ
        case(4)
          detJ=1/16.0*(3*PAT5GITTER_A*PAT5GITTER_B/4.0*(1-1/3.0*gaussP(2)))
          matrixInvJ(1,1)=PAT5GITTER_B/2.0/detJ
          matrixInvJ(1,2)=(0.0)/detJ
          matrixInvJ(2,1)=(-PAT5GITTER_A/2.0)*gaussP(2)/detJ
          matrixInvJ(2,2)=3.0/2.0*PAT5GITTER_A*(1-1/3.0*gaussP(2))/detJ
     end select
    case(2)


    case(3)


    case(4)


    case(5)


  end select

  end subroutine linlgfillmatrixJpatch5elJ

!****************************************!
! linlgfillmatrixNinvpatch5elN        !
!
  subroutine linlgfillmatrixNpatch5elN(Mesh2D,ne,ni)
  use types
  implicit none

  type(tMeshGen)        :: Mesh2D
  integer               :: ne, ni ! element-Zahl und integrationspunkt-Zahl
! Form des Matrix N hier beschreibt in Referenzkoordinaten:
!-------------------------------------------------!
! N1,x  O     N2,x  0    N3,x   0     N4,x  0     !
! 0     N1,y  0     N2,y  0     N3,y  0     N4,y  !
! N1,y  N1,x  N2,y  N2,x  N3,y  N3,x  N4,y  N4,x  !
!-------------------------------------------------!
  integer i,j,k

  select case(ne)
    case(1)
      do i=1,3
        do j=1,8
          MatrixN(i,j)=1/4.* QuadMatrixNConst(ni*24+i*8+j) ! TBC with Fortran col convention
        enddo
      enddo

    case(2)
      do i=1,3
        do j=1,8
          MatrixN(i,j)=1/4.* QuadMatrixNConst(ni*24+i*8+j)
        enddo
      enddo

    case(3)
      do i=1,3
        do j=1,8
          MatrixN(i,j)=1/4.* QuadMatrixNConst(ni*24+i*8+j)
        enddo
      enddo

    case(4)
      do i=1,3
        do j=1,8
          MatrixN(i,j)=1/4.* QuadMatrixNConst(ni*24+i*8+j)
        enddo
      enddo

    case(5)
      do i=1,3
        do j=1,8
          MatrixN(i,j)=1/4.* QuadMatrixNConst(ni*24+i*8+j)
        enddo
      enddo

  end select
  end subroutine linlgfillmatrixNpatch5elN

!****************************************
! linlgprodmatrixBpatch5elB
  subroutine linlgprodmatrixBpatch5elB
  use types
  implicit none
  ! no input data since the matrix are "static" vars at the 
  ! top of this package (see allocate)
  integer              :: i,j,k


  i=1
  do while (i.ne.4)
    ! 1. first step, N,x=N,u*Jinv(1,1)+N,v*Jinv(2,1)
    matrixB(1,2*i-1) = matrixN(1,2*i-1)*matrixInvJ(1,1) + matrixN(2,2*i)*matrixInvJ(2,1)
    ! 2. second step, N,y=N,u*Jinv(1,2)+N,v*Jinv(2,2)
    matrixB(2,2*i) = matrixN(1,2*i-1)*matrixInvJ(1,2) + matrixN(2,2*i)*matrixInvJ(2,2)
    ! 3. third step, copy these derivative in third line
    matrixB(3,2*i-1) = matrixB(2,2*i)
    matrixB(3,2*i)   = matrixB(1,2*i-1)
    
    i=i+1  ! next group of columns
  end do

  end subroutine linlgprodmatrixBpatch5elB


!****************************************!
! linlgintegmatrixpatch5el                !
! four calls to the fillmatrixprocedure
! to have [N,x],i  i=integ point
! and use matrix product to compute Ke
  subroutine linlgintegmatrixpatch5el(Mesh2D,ne)
  use types
  implicit none

  type(tMeshGen)        :: Mesh2D
  integer               :: ne
  integer               :: ni,i,j,k
  real                  :: sumup

 ! init
  do i=1,8
      do j=1,8
         matrixKsum(i,j)=0.0
      enddo
  enddo
  ni=1
  do while (ni.le.4) ! Perform a quadrature sum over the quadrangle
    ! call the subroutines "linlgfillmatrix.." which set the coefficients
    ! depending from element ne, integration-point ni
    call linlgfillmatrixXpatch5elX(Mesh2D,ne)
    call linlgfillmatrixJpatch5elJ(Mesh2D,ne,ni) ! TBC
    call linlgfillmatrixNpatch5elN(Mesh2D,ne,ni)
    call linlgprodmatrixBpatch5elB  ! TBC
    ! calculate K=tr(B)*D*B in several steps:
    do i=1,3
      do j=1,8
        sumup=0.0
        do k=1,3
           sumup=sumup+ Planestress(i*3+k)*matrixB(k,j)
        enddo
        matrixKimd1(i,j)=sumup  ! intermediary product
      enddo
    enddo
    do i=1,8
      do j=1,3
        matrixKtrB(i,j)=matrixB(j,i)
      enddo
    enddo
    do i=1,8
      do j=1,8
        sumup=0.0
        do k=1,3
           sumup=sumup+ matrixKtrB(i,k)*matrixKimd1(k,j)
        enddo
        matrixKinc(i,j)=sumup
      enddo
    enddo
   ! add the contribution of the particular integr point to the result
    do i=1,8
      do j=1,8
         matrixKsum(i,j)=matrixKsum(i,j)+matrixKinc(i,j)
      enddo
    enddo
    ni=ni+1 ! next point

  end do


  end subroutine linlgintegmatrixpatch5el


!****************************************!
! linlgdeallocmatrixpatch5el               !
!
  subroutine linlgdeallocmatrixpatch5el(Mesh2D, RB, Uvar, rhs, SparseSkalarProdRow1, &
               SparseSkalarProdCol, SparseSkalarProdVal)
  use types

  implicit none
  type(tMeshGen)      :: Mesh2D
  type(tRandbedingungen) :: RB
  real,pointer        :: Uvar(:)
  real,pointer        :: rhs(:)
  integer             :: allocStat  
  integer,dimension(4) :: SparseDiagLength=(/ 8, 8, 8, 8 /)
  real,pointer        :: SparseSkalarProdRow1(:)
  real,pointer        :: SparseSkalarProdCol(:)
  real,pointer        :: SparseSkalarProdVal(:)
  intent(in)          :: Mesh2D
  intent(inout)       :: SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal

  deallocate(RB%randu, RB%rando, RB%randl, RB%randr, & 
          Uvar, rhs, SparseSkalarProdRow1, SparseSkalarProdCol, &
          SparseSkalarProdVal, STAT = allocStat)

  deallocate(matrixJ, matrixX, matrixInvJ, matrixN, & 
           matrixKinc, matrixKsum, matrixB, matrixKimd1, matrixKtrB, & 
           STAT=allocStat)
  end subroutine linlgdeallocmatrixpatch5el

end module Patchtest2d5el_mod
