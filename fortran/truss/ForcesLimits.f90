module ForceLimits_mod
 use TypesBalken
 implicit none

 private

  interface ExtractTrussUnknowns
    module procedure ExtractTrussUnknowns
  end interface

 public :: ExtractTrussUnknowns

 contains

!****************************************
! Sub ExtractTrussUnknowns
!   accepts a whole rigidity matrix and copies
!   only a smaller submatrix for the unknown
!   displacement (columns->rows in Fortran) 
!   and knowns forces (rows->columns in Fortran)
 subroutine ExtractTrussUnknowns(MeshInfo, Dunknowns, Fmovemt, Ke_H, A_H, B_H)
 use TypesBalken
 implicit none
 type(tMeshInfo)                   :: MeshInfo
 integer,pointer                   :: Dunknowns(:,:)     ! 1: disp imposed, 0: disp to determine
 integer,pointer                   :: Fmovemt(:,:)       ! why this, if symetrical to Dunknonwns?
 type(tRigidFullMat)               :: Ke_H
 type(tRigidFullMat)               :: A_H, B_H
 ! Vorhabe
 intent(in)                        :: MeshInfo,Dunknowns,Fmovemt,Ke_H
 intent(inout)                     :: A_H, B_H
 ! lokale
 integer                           :: nn,nt,i,j,ni
 
 nt=MeshInfo%nt
 nn=MeshInfo%nn
 
 do i=1,nt
! if
 enddo 
 end subroutine ExtractTrussUnknowns


end module ForceLimits_mod

