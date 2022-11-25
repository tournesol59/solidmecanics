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

! important
#define DEF_QUAD_PATCHTEST 0  ! 1=only for patch test, 0=use for generic mesh
#define FILEOUTPUT 27

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
  real,dimension(96) :: QuadMatrixNConst = (/ (/ (/ -1.57735, 0.0, 1.57735, 0.0, 0.42265, 0.0, -0.42265, 0.0 /), &
                           (/ 0.0, -1.57735, 0.0, -0.42265, 0.0, 0.42265, 0.0, 1.57735 /), &
          (/ -1.57735, -1.57735, -0.42265, 1.57735, 0.42265, 0.42265, 1.57735,-0.42265 /) /), &
                        (/ (/ -1.57735, 0.0, 1.57735, 0.0, 0.42265, 0.0, -0.42265, 0.0 /), &
                           (/ 0.0, -0.42265, 0.0, -1.57735, 0.0, 1.57735, 0.0, 0.42265 /), &
          (/ -0.42265, -1.57735, -1.57735, 1.57735, 1.57735, 0.42265, 0.42265, -0.42265 /) /), &
                        (/ (/ -0.42265, 0.0, 0.42265, 0.0, 1.57735, 0.0, -1.57735, 0.0 /), &
                           (/ 0.0, -0.42265, 0.0, -1.57735, 0.0, 1.57735, 0.0, 0.42265 /), &
          (/ -0.42265, -0.42265, -1.57735, 0.42265, 1.57735, 1.57735, 0.42265, -1.57735 /) /), &
                        (/ (/ -0.42265, 0.0, 0.42265, 0.0, 1.57735, 0.0, -1.57735, 0.0 /), &
                           (/ 0.0, -1.57735, 0.0, -0.42265, 0.0, 0.42265, 0.0, 1.57735 /), &
          (/ -1.57735, -0.42265, -0.42265, 0.42265, 0.42265, 1.57735, 1.57735, -1.57735 /) /) /) 

! global variable definitions in this package
  real,pointer            :: matrixJ(:,:)  !dimension(1:4,1:2)
  real,pointer            :: matrixX(:,:)  ! dimension(1:2,1:4) 
  real,pointer            :: matrixDX(:,:)  ! dimension(1:2,1:4)
  real,pointer            :: matrixInvJ(:,:)! dimension(1:4,1:2)
  real,pointer            :: matrixN(:,:)  ! dimension(1:8,1:3)
  real,pointer            :: matrixB(:,:)  ! dimension(1:3,1:8)
  real,pointer            :: matrixKimd1(:,:)  ! dimension(1:3,1:8)
  real,pointer            :: matrixKtrB(:,:)  ! dimension(1:8,1:3)
  real,pointer            :: matrixKinc(:,:)  ! dimension(1:8,1:8)*
!  real,pointer            :: matrixKsum(:,:)  ! dimension(1:8,1:8) this matrix must be accessible by ref from main and package LagrangePlan2D


  public :: linlgallocmatrixpatch5el, linlgdeallocmatrixpatch5el, linlgfillmatrixXpatch5elX, & 
          linlgfillmatrixJpatch5elJ, linlgfillmatrixJinvquadgen, & 
          linlgfillmatrixNpatch5elN, linlgintegmatrixpatch5el,linlgprodmatrixBpatch5elB, & 
             linlgassemblymatrixKpatch5el, linlgprintmatoneelmtpatch5el

contains
!****************************************!
! linlgallocmatrixpatch5el               !
!
  subroutine linlgallocmatrixpatch5el(matrixKsum)
  use types
  implicit none
  real,pointer        :: matrixKsum(:,:)  ! dimension(1:8,1:8)

  integer             :: allocStat  
  integer             :: i,j,randltype      ! RandTyp links: 1=Fixed, 2=Momentum, 3=Freibiegung
     real,pointer     :: randl(:)       ! Randwerte links

   allocate(matrixJ(1:4,1:2), matrixX(1:2,1:4), matrixInvJ(1:4,1:2), & 
           matrixN(1:3,1:8), matrixKinc(1:8,1:8), matrixKsum(1:8,1:8), & 
           matrixB(1:3,1:8), matrixKimd1(1:3,1:8), matrixKtrB(1:8,1:3), &
           matrixDX(1:2,1:4), STAT=allocStat)
   if (allocStat.ne.0) then 
      write(*,*) 'Error linlgallocmatrixpatch5el: could not allocate all elements'
   endif
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
!   fills the position coordinates of an element
!   since there are 5 elements: 5 times called
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
! why so much code? the idea is to validate the numeric sheme
! by focusing on the assembly and solving. Of course so much code is not generalizable
  subroutine linlgfillmatrixJpatch5elJ(Mesh2D,ne,ni)
  use types
  implicit none

  type(tMeshGen)        ::Mesh2D
  integer               :: ne,ni
  real                  :: Xki,Yeta, detJ  
  real,dimension(1:2,1:4) :: QuadderGauss  
  integer i,j,k

  ! *** initialization
  select case(ni)  ! select integration point
    case(1) 
       Xki=gaussP(1)
       Yeta=gaussP(1)
    case(2)
       Xki=gaussP(2)
       Yeta=gaussP(1)
    case(3)
       Xki=gaussP(2)
       Yeta=gaussP(2)
    case(4)
       Xki=gaussP(1)
       Yeta=gaussP(2)
  end select
  ! complete derivation of Quad polynom position based matrix
  QuadderGauss(1,1)=-(1-Yeta) 
  QuadderGauss(1,2)=(1-Yeta) 
  QuadderGauss(1,3)=(1+Yeta)
  QuadderGauss(1,4)=-(1+Yeta) 
  QuadderGauss(2,1)=-(1-Xki) 
  QuadderGauss(2,2)=(1-Xki)
  QuadderGauss(2,3)=(1+Xki)
  QuadderGauss(2,4)=(1+Xki)

  !*** calculus
  select case(ne)
    case(1)
 ! element 1 (A) below
!      select case(ni)
!        case(1)
!          detJ=1/16.0*(3*PAT5GITTER_A*PAT5GITTER_B/4.0)*(1+1/3.0*gaussP(1))
!          matrixInvJ(1,1)=(1/2.0*PAT5GITTER_B)/detJ
!          matrixInvJ(1,2)=(0)/detJ
!          matrixInvJ(2,1)=gaussP(1)*PAT5GITTER_A/2.0/detJ
!          matrixInvJ(2,2)=3.0*PAT5GITTER_A/2.0/detJ
!        case(2)
!          detJ=1/16.0*(3*PAT5GITTER_A*PAT5GITTER_B/4.0)*(1+1/3.0*gaussP(2))
!          matrixInvJ(1,1)=(1/2.0*PAT5GITTER_B)/detJ  ! Simplification _L=_D=_B/4
!          matrixInvJ(1,2)=(PAT5GITTER_B/2.0)/detJ
!          matrixInvJ(2,1)=PAT5GITTER_A/2.0*gaussP(2)/detJ ! Simplification _K=_A*3/4, _C=_A/4
!          matrixInvJ(2,2)=3.0/2.0*PAT5GITTER_A+1/2.0*PAT5GITTER_A*gaussP(2)/detJ
!        case(3)
!          detJ=1/16.0*(3*PAT5GITTER_A*PAT5GITTER_B/4.0)*(1-1/3.0*gaussP(2))   
!          matrixInvJ(1,1)=(PAT5GITTER_B/2.0)/detJ
!          matrixInvJ(1,2)=(0.0)/detJ
!          matrixInvJ(2,1)=PAT5GITTER_A/2.0*gaussP(2)/detJ
!          matrixInvJ(2,2)=3.0/2.0*PAT5GITTER_A*(1-1/3.0*gaussP(2))/detJ
!        case(4)
!          detJ=1/16.0*(3*PAT5GITTER_A*PAT5GITTER_B/4.0*(1-1/3.0*gaussP(2)))
!          matrixInvJ(1,1)=PAT5GITTER_B/2.0/detJ
!          matrixInvJ(1,2)=(0.0)/detJ
!          matrixInvJ(2,1)=(-PAT5GITTER_A/2.0)*gaussP(2)/detJ
!          matrixInvJ(2,2)=3.0/2.0*PAT5GITTER_A*(1-1/3.0*gaussP(2))/detJ
!     end select

      matrixInvJ(1,1)=(PAT5GITTER_L+PAT5GITTER_D)+(PAT5GITTER_L-PAT5GITTER_D)*Xki
      matrixInvJ(1,2)=-(PAT5GITTER_L-PAT5GITTER_D)-(PAT5GITTER_L-PAT5GITTER_D)*Yeta
      matrixInvJ(2,1)=-(-PAT5GITTER_A+PAT5GITTER_K+PAT5GITTER_C)-(PAT5GITTER_A-PAT5GITTER_K-PAT5GITTER_C)*Xki
      matrixInvJ(2,2)=(PAT5GITTER_A+PAT5GITTER_K-PAT5GITTER_C)+(PAT5GITTER_A-PAT5GITTER_K+PAT5GITTER_C)*Xki
      detJ=matrixInvJ(1,1)*matrixInvJ(2,2)-matrixInvJ(1,2)*matrixInvJ(2,1)
      matrixInvJ(1,1)=matrixInvJ(1,1)/detJ
      matrixInvJ(1,2)=matrixInvJ(1,2)/detJ
      matrixInvJ(2,1)=matrixInvJ(2,1)/detJ
      matrixInvJ(2,2)=matrixInvJ(2,2)/detJ

!    case(2)
!      select case(ni)
! Element 2 (B): right
!        case(1)
!          detJ=1/16.0*(3*PAT5GITTER_A*PAT5GITTER_B/4.0)*(1+1/3.0*gaussP(1))
!          matrixInvJ(1,1)=(3/2.0*PAT5GITTER_B*(1.0+1/3.0*gaussP(1)))/detJ
!          matrixInvJ(1,2)=gaussP(1)*PAT5GITTER_B/2.0/detJ
!          matrixInvJ(2,1)=(0)/detJ
!          matrixInvJ(2,2)=PAT5GITTER_A/2.0/detJ

!      end select
! Element 3 (C): top
!    case(3)
!      select case(ni)
!        case(1)
!          detJ=1/16.0*(3*PAT5GITTER_A*PAT5GITTER_B/4.0)*(1+1/3.0*gaussP(1))
!          matrixInvJ(1,1)=(0)/detJ
!          matrixInvJ(1,2)= -PAT5GITTER_B/2.0/detJ
!          matrixInvJ(2,1)=(3/2.0*PAT5GITTER_A*(1.0+1/3.0*gaussP(1)))/detJ
!          matrixInvJ(2,2)=PAT5GITTER_A/2.0*gaussP(2)/detJ
 
!      end select
! Element 4 (D): left
!    case(4)
!      select case(ni)
!        case(1)
!        
        !        case(2)
        !          detJ=1/16.0*(3*PAT5GITTER_A*PAT5GITTER_B/4.0)*(1+1/3.0*gaussP(1))
        !          matrixInvJ(1,1)=(3/2.0*PAT5GITTER_B*(1.0+1/3.0*gaussP(1)))/detJ
        !          matrixInvJ(1,2)=gaussP(1)*PAT5GITTER_B/2.0/detJ
        !          matrixInvJ(2,1)=(0)/detJ
        !          matrixInvJ(2,2)=PAT5GITTER_A/2.0/detJ
        ! 
        !      end select
        ! Element 5 (Z): center
        !    case(5)

        !      select case(ni)
        !        case(1)
        ! 
        !      end select
           case(2:5)
           ! element A was decomposed in order to test and validate the code later
           ! Element B to Z : general code
              ! not completely necessary but improves the understandability: calc. DX
              matrixDX(1,1)=0.
              matrixDX(2,1)=0.
              matrixDX(1,2)=matrixX(1,2)-matrixX(1,1)
              matrixDX(2,2)=matrixX(2,2)-matrixX(2,1)
              matrixDX(1,3)=matrixX(1,3)-matrixX(1,1)
              matrixDX(2,3)=matrixX(2,3)-matrixX(2,1)
              matrixDX(1,4)=matrixX(1,4)-matrixX(1,1)
              matrixDX(2,4)=matrixX(2,4)-matrixX(2,1)
              matrixInvJ(1,1)=0.
              matrixInvJ(1,2)=0.
              matrixInvJ(2,1)=0. 
              matrixInvJ(2,2)=0. 
              do j=1,4
              ! prepare inverse: the 2x2 ;atrix of cofactors is assigned insteadof J
                matrixInvJ(1,1)=matrixInvJ(1,1)+QuadderGauss(2,j)*matrixDX(2,j)!co11=J22
                matrixInvJ(1,2)=matrixInvJ(1,2)-QuadderGauss(1,j)*matrixDX(2,j)! -J12
                matrixInvJ(2,1)=matrixInvJ(2,1)-QuadderGauss(2,j)*matrixDX(1,j)! -J21
                matrixInvJ(2,2)=matrixInvJ(2,2)+QuadderGauss(1,j)*matrixDX(1,j)! J11
              enddo
              detJ=matrixInvJ(1,1)*matrixInvJ(2,2)-matrixInvJ(1,2)*matrixInvJ(2,1)
              matrixInvJ(1,1)=matrixInvJ(1,1)/detJ
              matrixInvJ(1,2)=matrixInvJ(1,2)/detJ
              matrixInvJ(2,1)=matrixInvJ(2,1)/detJ
              matrixInvJ(2,2)=matrixInvJ(2,2)/detJ
          end select ! end select element
          ! recopy bottom part of the matrix, because it will operate on 
          matrixInvJ(3,1)=matrixInvJ(1,1)
          matrixInvJ(3,2)=matrixInvJ(1,2)
          matrixInvJ(4,1)=matrixInvJ(2,1)
          matrixInvJ(4,2)=matrixInvJ(2,2)

  end subroutine linlgfillmatrixJpatch5elJ

!****************************************!
! linlgfillmatrixJinvquadgen             !
!
  subroutine linlgfillmatrixJinvquadgen(Mesh2D,ne,ni)
  use types
  implicit none

  type(tMeshGen)        :: Mesh2D
  integer               :: ne,ni
  real                  :: detJ  
  real,dimension(1:2,1:4) :: QuadderGauss  
  integer i,j,k

  ! copy derivatives Quad polynom position matrix from Patchtest5el2d packg
  QuadderGauss(1,1)=QuadMatrixNConst((ni-1)*24+1)
  QuadderGauss(1,2)=QuadMatrixNConst((ni-1)*24+3)
  QuadderGauss(1,3)=QuadMatrixNConst((ni-1)*24+5)
  QuadderGauss(1,4)=QuadMatrixNConst((ni-1)*24+7)
  QuadderGauss(2,1)=QuadMatrixNConst((ni-1)*24+8+2)
  QuadderGauss(2,2)=QuadMatrixNConst((ni-1)*24+8+4)
  QuadderGauss(2,3)=QuadMatrixNConst((ni-1)*24+8+6)
  QuadderGauss(2,4)=QuadMatrixNConst((ni-1)*24+8+8)

  matrixDX(1,1)=0.
  matrixDX(2,1)=0.
  matrixDX(1,2)=Mesh2D%x(ELMTi_NODEj( ne,2 )) - Mesh2D%x(ELMTi_NODEj( ne,1 )) 
  matrixDX(2,2)=Mesh2D%y(ELMTi_NODEj( ne,2 )) - Mesh2D%y(ELMTi_NODEj( ne,1 ))
  matrixDX(1,3)=Mesh2D%x(ELMTi_NODEj( ne,3 )) - Mesh2D%x(ELMTi_NODEj( ne,1 )) 
  matrixDX(2,3)=Mesh2D%y(ELMTi_NODEj( ne,3 )) - Mesh2D%y(ELMTi_NODEj( ne,1 )) 
  matrixDX(1,4)=Mesh2D%x(ELMTi_NODEj( ne,4 )) - Mesh2D%x(ELMTi_NODEj( ne,1 )) 
  matrixDX(2,4)=Mesh2D%y(ELMTi_NODEj( ne,4 )) - Mesh2D%y(ELMTi_NODEj( ne,1 )) 
  matrixInvJ(1,1)=0.
  matrixInvJ(1,2)=0.
  matrixInvJ(2,1)=0. 
  matrixInvJ(2,2)=0. 
 
  do j=1,4
              ! prepare inverse: the 2x2 ;atrix of cofactors is assigned insteadof J           
     matrixInvJ(1,1)=matrixInvJ(1,1)+QuadderGauss(2,j)*matrixDX(2,j)!co11=J22
     matrixInvJ(1,2)=matrixInvJ(1,2)-QuadderGauss(1,j)*matrixDX(2,j)! -J12
     matrixInvJ(2,1)=matrixInvJ(2,1)-QuadderGauss(2,j)*matrixDX(1,j)! -J21
     matrixInvJ(2,2)=matrixInvJ(2,2)+QuadderGauss(1,j)*matrixDX(1,j)! J11
  enddo
  detJ=matrixInvJ(1,1)*matrixInvJ(2,2)-matrixInvJ(1,2)*matrixInvJ(2,1)
  matrixInvJ(1,1)=matrixInvJ(1,1)/detJ
  matrixInvJ(1,2)=matrixInvJ(1,2)/detJ
  matrixInvJ(2,1)=matrixInvJ(2,1)/detJ
  matrixInvJ(2,2)=matrixInvJ(2,2)/detJ
  ! recopy bottom part of the matrix, because it will operate on 
  matrixInvJ(3,1)=matrixInvJ(1,1)
  matrixInvJ(3,2)=matrixInvJ(1,2)
  matrixInvJ(4,1)=matrixInvJ(2,1)
  matrixInvJ(4,2)=matrixInvJ(2,2)

  end subroutine linlgfillmatrixJinvquadgen

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

  select case(ne) ! why no comments? because the following is repetitive
    case(1)
      do i=1,3
        do j=1,8
          matrixN(i,j)=1/4.* QuadMatrixNConst((ni-1)*24+(i-1)*8+j) ! TBC with Fortran col convention
        enddo
      enddo

    case(2)
      do i=1,3
        do j=1,8
          matrixN(i,j)=1/4.* QuadMatrixNConst((ni-1)*24+(i-1)*8+j)
        enddo
      enddo

    case(3)
      do i=1,3
        do j=1,8
          matrixN(i,j)=1/4.* QuadMatrixNConst((ni-1)*24+(i-1)*8+j)
        enddo
      enddo

    case(4)
      do i=1,3
        do j=1,8
          matrixN(i,j)=1/4.* QuadMatrixNConst((ni-1)*24+(i-1)*8+j)
        enddo
      enddo

    case(5)
      do i=1,3
        do j=1,8
          matrixN(i,j)=1/4.* QuadMatrixNConst((ni-1)*24+(i-1)*8+j)
        enddo
      enddo

  end select
  end subroutine linlgfillmatrixNpatch5elN

!****************************************
! linlgprodmatrixBpatch5elB: derivative matrix
  subroutine linlgprodmatrixBpatch5elB()
  use types
  implicit none
  ! if a matrix wise product would have been adopted this were:
  !                        
  ! ! J1,k  J1,l  J1,k  J1,l   !  !Ni,k  0   !   !Ni,x  0  !
  ! ! J2,k  J2,l  J2,k  J2,l   !x !Ni,l  0   ! = !0   Ni,y !
  !                               !0     Ni,k!   
  !                               !0     Ni,l!
  !
  ! where J is the inverse of the jacobian. k and l
  ! but a functional form is preferred here
  ! no input data since the matrix are "static" vars at the 
  ! top of this package (see allocate)
  integer              :: i,j,k

  i=1
  do while (i.le.4)
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
  subroutine linlgintegmatrixpatch5el(Mesh2D,ne,matrixKsum)
  use types
  implicit none

  type(tMeshGen)        :: Mesh2D
  integer               :: ne
  real,pointer          :: matrixKsum(:,:)  ! dimension(1:8,1:8)

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
#ifdef DEF_QUAD_PATCHTEST
    call linlgfillmatrixXpatch5elX(Mesh2D,ne)
    call linlgfillmatrixJpatch5elJ(Mesh2D,ne,ni) ! TBC
    call linlgfillmatrixNpatch5elN(Mesh2D,ne,ni)
    call linlgprodmatrixBpatch5elB()  ! TBC
#else

    call linlgfillmatrixJinvquadgen(Mesh2D,ne,ni) ! TBC
    call linlgfillmatrixNpatch5elN(Mesh2D,ne,ni)
    call linlgprodmatrixBpatch5elB()  ! TBC
#endif
    ! calculate K=tr(B)*D*B in several steps:
    do i=1,3
      do j=1,8
        sumup=0.0
        do k=1,3
        sumup=sumup+ Planestress((i-1)*3+k)*matrixB(k,j) ! are there any errors of type (i-1) ?
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
    ! log into file the intermediary results
    call linlgprintmatoneelmtpatch5el(ne)
    ! add the contribution of the particular integr point to the result
    do i=1,8
      do j=1,8
         matrixKsum(i,j)=matrixKsum(i,j)+matrixKinc(i,j)
      enddo
    enddo
    ni=ni+1 ! next point
  end do
! log into file the end rigidity matrix   
  write(FILEOUTPUT,10) "----- ", adjustl("--------matrixKsum")
  do i=1,8
     write(FILEOUTPUT,18) (matrixKsum(i,j), j=1,8)
  enddo
10 format (a)
11 format (i10)
12 format (2f10.5)
13 format (3f10.5)
14 format (4f10.5)
18 format (8f10.5)

  end subroutine linlgintegmatrixpatch5el

!***************************************
! linlgassemblymatrixKpatch5el
 subroutine linlgassemblymatrixKpatch5el(Mesh2D, matrixKsum, ne, SparseSkalarProdRow1, &
                SparseSkalarProdCol, SparseSkalarProdVal) 
 use types
 
 implicit none
  type(tMeshGen)      :: Mesh2D
  real,pointer        :: matrixKsum(:,:)  ! dimension(1:8,1:8)
  integer             :: ne
  real,pointer        :: SparseSkalarProdRow1(:)
  real,pointer        :: SparseSkalarProdCol(:)
  real,pointer        :: SparseSkalarProdVal(:)
  intent(in)          :: Mesh2D, ne
  intent(inout)       :: SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal
  integer             :: i,j,k,l,m
  integer,dimension(8) :: indexA1=(/1,2,3,4,9,10,11,12/) ! from (1,2,5,6)
  integer,dimension(8) :: indexB2=(/9,10,3,4,5,6,13,14/) ! from (5,2,3,7)
  integer,dimension(8) :: indexC3=(/15,16,13,14,5,6,7,8/) ! from (8,7,3,4)
  integer,dimension(8) :: indexD4=(/1,2,9,10,15,16,7,8/) ! from (1,5,8,4)
  integer,dimension(8) :: indexZ5=(/9,10,11,12,13,14,15,16/) ! from (5,6,7,8)
  
  select case(ne)

    case(1)
       do i=1,8
         k=indexA1(i)
         do j=1,8
           l=indexA1(j)
           m=1 !(k-1)*16+l
           SparseSkalarProdVal(m)=SparseSkalarProdVal(m)+matrixKsum(i,j)
         enddo
       enddo
    case(2)
       do i=1,8
         k=indexB2(i)
         do j=1,8
           l=indexB2(j)
           m=(k-1)*16+l
           SparseSkalarProdVal(m)=SparseSkalarProdVal(m)+matrixKsum(i,j)
         enddo
       enddo
    case(3)
       do i=1,8
         k=2*indexC3(i)
         do j=1,8
           l=2*indexC3(j)
           m=(k-1)*16+l
           SparseSkalarProdVal(m)=SparseSkalarProdVal(m)+matrixKsum(i,j)
         enddo
       enddo
     case(4)
          do i=1,8
            k=2*indexD4(i)
            do j=1,8
              l=2*indexA1(j)
              m=(k-1)*16+l
              SparseSkalarProdVal(m)=SparseSkalarProdVal(m)+matrixKsum(i,j)
            enddo
          enddo
    case(5)
       do i=1,8
         k=2*indexZ5(i)
         do j=1,8
           l=2*indexA1(j)
           m=(k-1)*16+l
           SparseSkalarProdVal(m)=SparseSkalarProdVal(m)+matrixKsum(i,j)
         enddo
       enddo
 end select
 ! print to output with 8x8 blocks
 write(*,501) "------------ linlgassemblymatrixKpatch5el"

 do k=1,2
 do l=1,2
    do i=1,8
       write(*,504) "patch matrix block ", (2*(k-1)+l)
       write(*,501) "|   x1   |   x2   |   x3   |   x4   |   x5   |   x6   |   x7   |   x8   |"
 
       write(*,503) (SparseSkalarProdVal(128*(k-1)+64*(l-1)+8*(i-1)+j), j=1,8)

   enddo
 enddo
 enddo
501 format(a)
502 format(8i10)
503 format(8f10.5)
504 format(a,i10)

 end subroutine linlgassemblymatrixKpatch5el
 
!*******************************************
! linlgasssemblysparsematKgen
! to do:

!****************************************!
! linlgdeallocmatrixpatch5el               !
!
  subroutine linlgdeallocmatrixpatch5el()
  use types

  implicit none
  integer             :: lli, allocStat  

  deallocate(matrixJ, matrixX, matrixDX, matrixInvJ, matrixN, & 
           matrixKinc, matrixKsum, matrixB, matrixKimd1, matrixKtrB, & 
           STAT=allocStat)

  end subroutine linlgdeallocmatrixpatch5el

!*******************************************
! subroutine linlgprintmatoneelmtpatch5el
subroutine linlgprintmatoneelmtpatch5el(ne)
  use types
  implicit none
  integer             :: ne
!  real,pointer        :: SparseSkalarProdRow1(:)
!  real,pointer        :: SparseSkalarProdCol(:)
!  real,pointer        :: SparseSkalarProdVal(:)
!  intent(in)          :: ne, SparseSkalarProdRow1, & 
!         SparseSkalarProdCol, SparseSkalarProdVal
  integer             :: i,j,k

  write(FILEOUTPUT,100) "----- ", adjustl("----print mat one element patch")
  ! matrix tr Jinv
  write(FILEOUTPUT,100) "----- ", adjustl("--------matrixJinv")
  do i=1,4
    write(FILEOUTPUT,102) (matrixinvJ(i,j), j=1,2)
  enddo
  write(FILEOUTPUT,100) "----- ", adjustl("--------matrixX")
  do i=1,2
    write(FILEOUTPUT,104) (matrixX(i,j), j=1,4)
  enddo
  write(FILEOUTPUT,100) "----- ", adjustl("--------matrixKtrB")
  do i=1,8
    write(FILEOUTPUT,103) (matrixKtrB(i,j), j=1,3)
  enddo
  write(FILEOUTPUT,100) "----- ", adjustl("--------matrixD*B")
  do i=1,3
    write(FILEOUTPUT,108) (matrixKimd1(i,j), j=1,8)
  enddo
  write(FILEOUTPUT,100) "----- ", adjustl("--------matrixKinc")
  do i=1,8
    write(FILEOUTPUT,108) (matrixKinc(i,j), j=1,8)
  enddo

100 format (a)
101 format (i10)
102 format (2f10.5)
103 format (3f10.5)
104 format (4f10.5)
108 format (8f10.5)
end subroutine linlgprintmatoneelmtpatch5el


end module Patchtest2d5el_mod
