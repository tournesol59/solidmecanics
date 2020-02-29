module RigidTrussContacts_mod

!**********************************************************************!
! The Balken have 2 nodes, each of whom have 3 displacements and 3 rotations
!**********************************************************************!
  use typesbalken
  ! use printabstractstruct_mod

  implicit none

#define mgetI(n) connectelmt(n,2)

#define mgetJ(n) connectelmt(n,3)

#define mcantileverI(n) ((Dunknowns(n,1).eq.0).and.(Dunknowns(n,2).eq.0).and.(Dunknowns(n,3).eq.0))

#define mcantileverJ(n) ((Dunknowns(n,4).eq.0).and.(Dunknowns(n,5).eq.0).and.(Dunknowns(n,6).eq.0))

#define marticulatedI(n) ((Fmovemt(n,3).eq.1))

#define marticulatedJ(n) ((Fmovemt(n,6).eq.1))

  private

    interface MatrixRotMultiply
      module procedure MatrixRotMultiply
    end interface

    interface Matrice_Ke_H_local_dfix
      module procedure Matrice_Ke_H_local_dfix
    end interface

    interface Matrice_Ke_H_Truss
      module procedure Matrice_Ke_H_Truss      
    end interface

    interface Vector_F_RHS_local
      module procedure Vector_F_RHS_local
    end interface

    interface Vector_F_RHS_Truss
      module procedure Vector_F_RHS_Truss
    end interface

    integer    :: TYP_BAR=1
    integer    :: TYP_BEAM=2

   public::  MatrixRotMultiply, Matrice_Ke_H_local_dfix, Matrice_Ke_H_Truss, Vector_F_RHS_local, Vector_F_RHS_Truss

  contains

!**************************************************************!
!							       !
!  Subroutine MatrixRotMultiply                                !
! operate left then right multiplication matrix-wise (rotation)!
! assumes square matrices order n, Rigid matrix is containes   !
! in a structure, while Rortho is a generatable
!**************************************************************!

SUBROUTINE MatrixRotMultiply(Krigid, Rortho, n)
  use typesbalken
  implicit none

  type(tRigidMat), intent(inout)  :: Krigid !!assumed 3x3
  real,dimension(9),intent(in)    :: Rortho
  integer,intent(in)              :: n
  
! lokale Variablen
  real,dimension(9)               :: Kcopy
  integer                         :: row,col,k

  do row=1,n
    do col=1,n      
      Kcopy(3*(col-1)+row)=0.0
      do k=1,n
         Kcopy(3*(col-1)+row)= Kcopy(3*(col-1)+row) + Krigid%Ke(row,k)*Rortho(3*(col-1)+k)
      enddo
    enddo
  enddo

  do row=1,n
    do col=1,n      
      Krigid%Ke(row,col)=0.0
      do k=1,n
!!        Krigid%Ke(row,col)= Krigid%Ke(row,col) + Rortho(3*(row-1)+k)*Kcopy(3*(col-1)+k) !! error in first ver
        Krigid%Ke(row,col)= Krigid%Ke(row,col) + Rortho(3*(k-1)+row)*Kcopy(3*(col-1)+k) ! multiply left by transpose of Rortho
      enddo
    enddo
  enddo

END SUBROUTINE MatrixRotMultiply

!**************************************************************!
!							       !
!  Subroutine Matrice_Ke_H_local                               !
!   creates and fills the matrix from 6 ddl to 6 force components!
!**************************************************************!

SUBROUTINE Matrice_Ke_H_local_dfix(MeshT, VarT, Keii, Keij, Kejj)
  use typesbalken
  use printabstractstruct_mod

  implicit none

  type(tMeshElmt), intent(in)  :: MeshT
  type(tVarElmt), intent(in)   :: VarT
  type(tRigidMat), intent(inout) :: Keii, Keij, Kejj
!  integer,intent(in)             :: code ! 1 fest eingespannt x2; 2 l-fest eingespannt r-freies Moment, 
!                                         ! 3 l-freies Moment r-fest eingespannt, 4 freies Moment x2
!  integer, dimension(6)         :: connect
! lokalen Variablen
  character(len=30)           :: SubName="Matrice_Ke_H_local_dfix       "
!  type(tRigidMat)             :: Ke_preii, Ke_preij, Ke_prejj
  real,pointer                :: M_rot(:)
  real                        :: EIsL, EIsL2, EIsL3, ESsL
  real                        :: costh, sinth
  integer                     :: i,j, allocstat
  integer,pointer             :: inull(:)
  
  allocate(inull(1), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"Matrice_Ke_H_local_dfix: error allocate null"
  endif
  allocate(M_rot(9), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"Matrice_Ke_H_local_dfix: error allocate M_rot"
  endif
  do i=1,9
    M_rot(i)=0.0
  enddo

  costh=(MeshT%endx-MeshT%startx)/MeshT%dlen
  sinth=(MeshT%endy-MeshT%starty)/MeshT%dlen

  ESsL=MeshT%EY*MeshT%SArea/MeshT%dlen
  EIsL=MeshT%EY*MeshT%CI/(MeshT%dlen)
  EIsL2=MeshT%EY*MeshT%CI/(MeshT%dlen**2)
  EIsL3=MeshT%EY*MeshT%CI/(MeshT%dlen**3)
  ! initialisierung
  do i=1,3
    do j=1,3
      Keii%Ke(i,j)=0.
      Keij%Ke(i,j)=0.
      Kejj%Ke(i,j)=0.
      M_rot(3*(i-1)+j)=0.
    enddo
  enddo
  print *, MeshT%typ
  ! ausfuellen
  if ((MeshT%typ).eq.TYP_BAR) then ! elmt is a bar, no flexion
    Keii%Ke(1,1) = ESsL    
    Keij%Ke(1,1) = -ESsL
    Kejj%Ke(1,1) = ESsL     
  elseif ((MeshT%typ).eq.TYP_BEAM) then  ! elmt is a beam with flexion around axis z only.
      Keii%Ke(1,1) = ESsL
      Keii%Ke(2,2) = 12*EIsL3
      Keii%Ke(2,3) = -6*EIsL2
      Keii%Ke(3,2) = -6*EIsL2
      Keii%Ke(3,3) = 4*EIsL

      Keij%Ke(1,1) = -ESsL
      Keij%Ke(2,2) = -12*EIsL3
      Keij%Ke(2,3) = 6*EIsL2
      Keij%Ke(3,2) = -6*EIsL2  !! is negative
      Keij%Ke(3,3) = 2*EIsL

      Kejj%Ke(1,1) = ESsL
      Kejj%Ke(2,2) = 12*EIsL3
      Kejj%Ke(2,3) = -6*EIsL2
      Kejj%Ke(3,2) = -6*EIsL2
      Kejj%Ke(3,3) = 4*EIsL

  endif

!  write(*,401) "----- ", adjustl(SubName)
! 401  format (a,a)

  M_rot(1)=costh
  M_rot(2)=sinth
  M_rot(4)=-sinth
  M_rot(5)=costh
  M_rot(9)=1.0
  
  call prntsimplearray(9, 2, inull, M_rot, SubName)

  call prntsimplestruct(3,3,Keii,SubName)

  call prntsimplestruct(3,3,Keij,SubName)

  call prntsimplestruct(3,3,Kejj,SubName)


  call MatrixRotMultiply(Keii, M_rot, 3)
  call prntsimplestruct(3,3,Keii,SubName)

  call MatrixRotMultiply(Keij, M_rot, 3)
  call prntsimplestruct(3,3,Keij,SubName)

  call MatrixRotMultiply(Kejj, M_rot, 3)
  call prntsimplestruct(3,3,Kejj,SubName)

  deallocate(inull, M_rot, STAT = allocstat)
  if (allocstat.ne.0) then
    print *,"Matrice_Ke_H_Truss: error deallocate M_rot"
  endif
END SUBROUTINE Matrice_Ke_H_local_dfix

!**************************************************************!
! Subroutine Matrice_Ke_H_assembly                             !
!   shall be called after each element submatrix               !
!**************************************************************!
SUBROUTINE Matrice_Ke_K_assembly(Ke_H,Keii,Keij,Kejj,jj,kk,n,ne)

  use typesbalken
!  use printabstractstruct_mod

  implicit none

  type(tRigidFullMat),intent(inout) :: Ke_H   ! Whole global rigidity matrix
  type(tRigidMat),intent(in)   :: Keii, Keij, Kejj
  integer,intent(in)           :: n,ne   ! Max number of nodes, max number of deg freedom (3 or 6)
  integer,intent(in)           :: jj,kk       ! index of connected nodes to be assemblied in Ke_H
  character(len=30)            :: SubName="Matrice_Ke_H_assembly         "
! lokalen Variablen
  integer                       :: i,j,l,p,allocstat


  do i=1,3
     do j=1,3
       Ke_H%Ke((jj-1)*ne+i, (jj-1)*ne+j) = Ke_H%Ke((jj-1)*ne+i, (jj-1)*ne+j) + Keii%Ke(i,j)
       Ke_H%Ke((jj-1)*ne+i, (kk-1)*ne+j) = Ke_H%Ke((jj-1)*ne+i, (kk-1)*ne+j) + Keij%Ke(i,j)
       Ke_H%Ke((kk-1)*ne+i, (jj-1)*ne+j) = Ke_H%Ke((kk-1)*ne+i, (jj-1)*ne+j) + Keij%Ke(j,i)
       Ke_H%Ke((kk-1)*ne+i, (kk-1)*ne+j) = Ke_H%Ke((kk-1)*ne+i, (kk-1)*ne+j) + Kejj%Ke(i,j)
     enddo
  enddo

END SUBROUTINE  Matrice_Ke_K_assembly

!*************************************************************!
!							       !
!  Subroutine Matrice_Ke_H_Truss                               !
!   shall be called in the main Program as there are elements  !
! not terminated
!**************************************************************!

SUBROUTINE Matrice_Ke_H_Truss(MeshT1, VarT1, Ke_H, jj, kk, n, ne, connectelmt, &
                              Dunknowns, Dimposed, &
                              Fmovemt, Fapplied)
  use typesbalken
  use printabstractstruct_mod

  implicit none

  type(tMeshElmt), intent(in)  :: MeshT1  ! One bar or beam element
  type(tVarElmt), intent(in)   :: VarT1   ! Unknown variables (not used in this subroutine)
                                          ! and applied forces (used)
  type(tRigidFullMat)          :: Ke_H   ! Whole global rigidity matrix
  integer,intent(in)           :: n, ne   ! Max number of nodes, max number of deg freedom (3 or 6)
  integer,intent(in)           :: jj,kk       ! index of connected nodes to be assemblied in Ke_H
  integer,pointer,intent(in)   :: connectelmt(:,:)  ! These indexes (j,k) are also found if connectelmt(j,k)=1
  integer,pointer,intent(in)   :: Dunknowns(:,:)    ! Dunkowns(i,j)=1 means j-th of 4 deg of 
                                                    ! freedom of i-th of n nodes is unknown, 0 by default
  real,pointer,intent(in)      :: Dimposed(:,:)  ! imponierte Bewegungen at move elsewhere Dunknowns(i,j)=-1
  integer,pointer,intent(in)   :: Fmovemt(:,:)      ! Fmovemt(i,j)=1 means there one Force to calculate at move j at node i
  real,pointer,intent(in)      :: Fapplied(:,:)  ! Force applied at move j at node i (provided Fmovemt(i,j)=0)
  character(len=30)            :: SubName="Matrice_Ke_H_Truss            "

! lokalen Variablen

  integer                       :: i,j,l,p,allocstat
  type(tRigidMat)               :: Keii, Keij, Kejj

!  allocate(Ke_H%Ke(1:n*ne,1:n*ne), STAT=allocstat) !! BY NO WAY !! (left as comment to not do that again!)
!  if (allocstat.ne.0) then
!    print *,"Matrice_Ke_H_Truss: error allocate global Ke_H"
!  endif
  allocate(Keii%Ke(1:3,1:3), Keij%Ke(1:3,1:3), Kejj%Ke(1:3,1:3), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"Matrice_Ke_H_Truss: error allocate Keii,ij,jj"
  endif
 
  call Matrice_Ke_H_local_dfix(MeshT1, VarT1, Keii, Keij, Kejj)  
  call Matrice_Ke_K_assembly(Ke_H,Keii,Keij,Kejj,jj,kk,n,ne)

  deallocate(Keii%Ke, Keij%Ke, Kejj%Ke, STAT = allocstat)
  if (allocstat.ne.0) then
    print *,"Matrice_Ke_H_Truss: error deallocate Keii,Keij,Kejj"
  endif
END SUBROUTINE Matrice_Ke_H_Truss

!**************************************************************!
!							       !
!  Subroutine Vector_F_RHS_local                               !
!   creates and fills the vector from 3 ddl corresp forces and Moments!
!**************************************************************!

SUBROUTINE Vector_F_RHS_local(MeshT, VarT, code)
  use typesbalken
  implicit none

  type(tMeshElmt), intent(in)  :: MeshT
  type(tVarElmt), intent(inout)   :: VarT  ! contains the force at nodes i and j
  integer,intent(in)              :: code ! 1 fest eingespannt x2; 2 l-fest eingespannt r-freies Moment, 
                                         ! 3 l-freies Moment r-fest eingespannt, 4 freies Moment x2
  ! TO COMPLETE !

END SUBROUTINE Vector_F_RHS_local

!***************************************************************!
!  Subroutine Vector_F_RHS _Truss                               !
!   shall be called in the main Program as there are nodes      !
!***************************************************************!

SUBROUTINE Vector_F_RHS_Truss(MeshT1, VarT1, V_RHS, j, k, n, ne, connectelmt)
  use typesbalken
  implicit none

  type(tMeshElmt), intent(in)  :: MeshT1  ! One bar or bulk element
  type(tVarElmt), intent(in)   :: VarT1   ! Unknown variables (not used in this subroutine)
                                          ! and applied forces (used)
  real,pointer,intent(inout)   :: V_RHS(:)   ! Whole global rigidity matrix,
  integer,intent(in)           :: n, ne   ! Max number of nodes, max number of deg freedom (3 or 6)
  integer,intent(in)           :: j,k       ! index of connected nodes to be assemblied in Ke_H
  integer,pointer,intent(in)   :: connectelmt(:,:)  ! These indexes (j,k) are also found if connectelmt(j,k)=1

! TO COMPLETE

END SUBROUTINE Vector_F_RHS_Truss


end module RigidTrussContacts_mod



