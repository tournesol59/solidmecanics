module RigidTrussContacts

!**********************************************************************!
! The Balken have 2 nodes, each of whom have 3 displacements and 3 rotations
!**********************************************************************!
  use typesbalken
  implicit none

  private
    interface Matrice_Ke_H_local
      module procedure Matrice_Ke_H_local
    end interface

    interface Matrice_Ke_H_Truss
      module procedure Matrice_Ke_H_Truss      
    end interface


   public::  Matrice_Ke_H_local, Matrice_Ke_H_Truss

  contains
!**************************************************************!
!							       !
!  Subroutine Matrice_Ke_H_local                               !
!   creates and fills the matrix from 6 ddl to 6 force components!
!**************************************************************!

SUBROUTINE Matrice_Ke_H_local(MeshT, VarT, Keii, Keij, Kejj)
  use typesbalken
  implicit none

  type(tMeshElmt), intent(in)  :: MeshT
  type(tVarElmt), intent(in)   :: VarT
  type(tRigidMat), intent(inout) :: Keii, Keij, Kejj
!  integer, dimension(6)         :: connect
! lokalen Variablen
  real                        :: EIsL, EIsL2, EIsL3, ESsL
  integer                     :: i,j

02 format (f10.5)
END SUBROUTINE Matrice_Ke_H_local

!**************************************************************!
!							       !
!  Subroutine Matrice_Ke_H_Truss                               !
!   shall be called in the main Program as there are elements  !
!**************************************************************!

SUBROUTINE Matrice_Ke_H_Truss(MeshT1, VarT1, Ke_H, j,k, n, ne, connectelmt, Dunknowns, Dimposed, Fmovemt, Fapplied)
  use typesbalken
  implicit none

  type(tMeshElmt), intent(in)  :: MeshT1  ! One bar or bulk element
  type(tVarElmt), intent(in)   :: VarT1   ! Unknown variables (not used in this subroutine)
                                          ! and applied forces (used)
  real,pointer,intent(inout)   :: Ke_H(:,:)   ! Whole global rigidity matrix,
  integer,intent(in)           :: n, ne   ! Max number of nodes, max number of deg freedom (3 or 6)
  integer,intent(in)           :: j,k       ! index of connected nodes to be assemblied in Ke_H
  integer,pointer,intent(in)   :: connectelmt(:,:)  ! These indexes (j,k) are also found if connectelmt(j,k)=1
  integer,pointer,intent(in)   :: Dunknowns(:,:)    ! Dunkowns(i,j)=1 means j-th of 4 deg of 
                                                    ! freedom of i-th of n nodes is unknown, 0 by default
  real,pointer,intent(in)      :: Dimposed(:,:)  ! imponierte Bewegungen at move elsewhere Dunknowns(i,j)=1
  integer,pointer,intent(in)   :: Fmovemt(:,:)      ! Fmovemt(i,j)=1 means there one Force to calculate at move j at node i
  real,pointer,intent(in)      :: Fapplied(:,:)  ! Force applied at move j at node i (provided Fmovemt(i,j)=0)

! lokalen Variablen

  integer                       :: i,l,p,allocStat
  type(tRigidMat)               :: Keii, Keij, Kejj

 ! real                          :: smallEIsL = (270000*0.001/0.5)*1e-4
 ! real                          :: smallEIsL2 = (270000*0.001/0.25)*1e-4
 ! real                          :: smallEIsL3 = (270000*0.001/0.125)*1e-4



  call Matrice_Ke_H_local(MeshT1, VarT1, Keii, Keij, Kejj)

END SUBROUTINE Matrice_Ke_H_Truss


end module RigidTrussContacts

