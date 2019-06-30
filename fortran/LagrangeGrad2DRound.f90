module LagrangeGrad2DRound_mod
! MODULE NOT DEVELOPPED FULLY NOR INTEGRATED 
! should implement the sprse matrix version
! Use: in Gesetz(1) f_ = Sum,k-i mu,k-i*Phi,i(x,y)

implicit none
! use TypesRound

private

  interface linlgallocatesparse
    module procedure linlgallocatesparse
  end interface

!  interface
!    module procedure linlgfillmatrix
!  end interface

!  interface
!    module procedure linlgfillbiegung
!  end interface

public :: linlgallocatesparse


contains 

!***********************************************
subroutine linlgallocatesparse(NTestMeshR,MeshR)
!----------------------------------------------
use TypesRound
implicit none

type(tRoundMeshInfo) :: NTestMeshR
type(tRoundMesh)     :: MeshR

intent(inout)        :: NTestMeshR, MeshR

  ! do nothing

end subroutine linlgallocatesparse


end module LagrangeGrad2DRound_mod
