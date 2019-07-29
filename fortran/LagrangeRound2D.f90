! Use: interpolate with a0+a1*X+a2*Y+a3*X*Y
module LagrangeRound2D_mod
!---
implicit none
private
  interface interpol3PSolve
     module procedure interpol3PSolve
  end interface
  interface interpol3PCoord
     module procedure interpol3PCoord
  end interface

public :: interpol3PSolve, interpol3PCoord  
! , interpol3PEvalDerivX, interpol3PEvalDerivY

contains 
!********************************************************************!
!  procedure interpol3PSolve !
!********************************************************************!
 SUBROUTINE interpol3PSolve(p1, p2, p3, coef)
 use typesround
 implicit none
  ! Argumenten
  type(tNode)         :: p1, p2 , p3
  real,dimension(6)   :: coef
  intent(in)          :: p1, p2, p3
  intent(inout)       :: coef
  !-----------------------------------------------------
  ! lokales Vars
  integer             :: k 
 
  coef(5)=0.0
  coef(6)=0.0

 END SUBROUTINE interpol3PSolve
!********************************************************************!
!  procedure interpol3PCoord !
!********************************************************************!
 SUBROUTINE interpol3PCoord( NTestMeshR, MeshR, MeshCoeffR)
 use typesround
 implicit none
  ! Argumenten
  type(tRoundMeshInfo):: NTestMeshR
  type(tRoundMesh)    :: MeshR
  type(tRoundCoeff)   :: MeshCoeffR
  intent(in)          :: NTestMeshR, MeshR
  intent(inout)       :: MeshCoeffR
  !-----------------------------------------------------
  ! lokales Vars
  integer             :: nn, ne, k ,l, allocStat
  integer             :: n1, n2, n3
  type(tNode)         :: p1, p2 ,p3, p4
  real,dimension(:),allocatable     :: coef
!  real                :: a1, a2, a3, a4 

  ne = NTestMeshR%nElem
  do k=1,ne 
     n1 = MeshR%elems(k)%numeros(1)
     n2 = MeshR%elems(k)%numeros(2)
     n3 = MeshR%elems(k)%numeros(3)
     do l=1,3
       p1%vects(l) = MeshR%nodes(n1)%vects(l)
       p2%vects(l) = MeshR%nodes(n2)%vects(l)
       p3%vects(l) = MeshR%nodes(n3)%vects(l)
     enddo
     allocate(coef(1:6),       &
         STAT = allocStat )
     call interpol3PSolve(p1,p2,p3, coef) 
     ! save result
     do l=1,6
       MeshCoeffR%Coefficients(k*6+l)=coef(l)
     ! next iteration
     enddo
  enddo

END SUBROUTINE interpol3PCoord 

end module LagrangeRound2D_mod

