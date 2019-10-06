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

SUBROUTINE Matrice_Ke_H_local(MeshT, VarT, Keii, Keij, Kejj, connect)
  use typesbalken
  implicit none

  type(tMeshElmt), intent(in)  :: MeshT
  type(tVarElmt), intent(in)   :: VarT
  type(tRigidMat), intent(inout) :: Keii, Keij, Kejj
  integer, dimension(6)         :: connect
! lokalen Variablen
   real                        :: EIsL, EIsL2, EIsL3, ESsL
  integer                     :: i,j

02 format (f10.5)
END SUBROUTINE Matrice_Ke_H_local

!**************************************************************!
!							       !
!  Subroutine Matrice_Ke_H_Truss                               !
!   creates a global matrix of 6 elements Balken with boundary conds!
!**************************************************************!
SUBROUTINE Matrice_Ke_H_Truss(MeshT1,VarT1,VarT2,VarT3,VarT4,VarT5,VarT6, Ke_H, connect6)
  use typesbalken
  implicit none

  type(tMeshElmt), intent(in)  :: MeshT1
  type(tVarElmt), intent(in)   :: VarT1, VarT2,VarT3,VarT4,VarT5,VarT6
  real,pointer,intent(inout) :: Ke_H(:) ! 24*24, intent(inout) ::

  integer, dimension(6,6),intent(in)   :: connect6
! lokalen Variablen
  integer, dimension(6)         :: connect
  integer                       :: i,j,k ,allocStat
  type(tRigidMat)               :: Keii, Keij, Kejj
  integer, dimension(6)         ::  nodeone= (/ 0, 2, 0, 0, 0, 0 /) !! Assume only one has connect >0
  integer, dimension(6)         :: nodead= (/ 0, 0, 0, 4, 0, 6 /)
  integer                       :: node_1, node_2

  real                          :: smallEIsL = (270000*0.001/0.5)*1e-4
  real                          :: smallEIsL2 = (270000*0.001/0.25)*1e-4
  real                          :: smallEIsL3 = (270000*0.001/0.125)*1e-4

  allocate(Ke_H(1:24*24), STAT = allocStat)
  if (allocStat.ne.0) then
    print *," Errro Allocation global matrix !"
  endif
  !1 let us find a mean to exploit better the conectitivity info of connect6
  ! remind connect(i2)=1 means 'nodes i is connected by articulation of axis 2, if 0 by fixed (closet)
    
  i=1

  do j=1,6
    connect(j)=connect6((i-1),j)
  enddo
  call Matrice_Ke_H_local(MeshT1, VarT1, Keii, Keij, Kejj, connect)
!**************************************************************************************
!! let an example: articulation at node 1 in y, one otation %vector ez => [2,4],  [4,4]
!!         u v w ph,th,ps


  do while(i<=6)
    !!! 1. check connectivity single Element  
    i=0
    node_1=0
    do while ((node_1.eq.(0)).and.(i.le.6))
         i=i+1
         node_1=nodeone(i)
    enddo
! node1(i) has non zero element
    j=0
    node_2=0
    do while ((node_2.eq.(0)).and.(j.le.6))
         j=j+1
         node_2=nodead(j)      
    enddo
  if(.not.((node_1.le.0).or.(node_2.le.0)) ) then
    if ( abs( Keij%Ke(node_1,node_2) ).le.(smallEIsL3) ) then
		!do nothing korrrekt
    else 
      print *,"ERROR lacks free Articulation mechanism"
    endif
  endif   
 enddo

END SUBROUTINE Matrice_Ke_H_Truss


end module RigidTrussContacts

