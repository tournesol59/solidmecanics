module RigidTrussContacts_mod

!**********************************************************************!
! The Balken have 2 nodes, each of whom have 3 displacements and 3 rotations
!**********************************************************************!
  use typesbalken
  use printabstractstruct_mod
  use CreateTruss_mod

  implicit none
 ! caution: duplicate preprocessor but could not find other mean (see CreateTruss.f)
#define LINK_BAREND 1
#define LINK_BARSLIDE 2
#define LINK_BEAMEND 3
#define LINK_BEAMSLIDE 4
#define LINK_BEAMFREE 5
#define LINK_BARTOBAR 6
#define LINK_BEAMTOBEAM 7

#define mgetI(n) connectelmt(n,2)

#define mgetJ(n) connectelmt(n,3)

#define mcantileverI(n,i) (Dunknowns(n,i).eq.1)  
! I have changed it recently to the inverse (1 insof 0) because intuitively a cantilever at the wall side is the determination of
! moment and force at the wall where D=forced zero
#define mcantileverJ(n,i) (Dunknowns(n,i+3).eq.1)
#define marticulatedI(n) (Fmovemt(n,3).eq.1)
#define marticulatedJ(n) (Fmovemt(n,6).eq.1)

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
! caution: linkage typ codes are defined by proprocessor in File: CreateTruss.f

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
  type(tPassageMat),intent(in)    :: Rortho
  integer,intent(in)              :: n
  
! lokale Variablen
  real,dimension(9)               :: Kcopy
  integer                         :: row,col,k

  do row=1,n
    do col=1,n      
      Kcopy(3*(col-1)+row)=0.0
      do k=1,n
         Kcopy(3*(col-1)+row)= Kcopy(3*(col-1)+row) + Krigid%Ke(row,k)*Rortho%Mp(3*(col-1)+k)
      enddo
    enddo
  enddo

  do row=1,n
    do col=1,n      
      Krigid%Ke(row,col)=0.0
      do k=1,n
!!        Krigid%Ke(row,col)= Krigid%Ke(row,col) + Rortho%Mp(3*(row-1)+k)*Kcopy(3*(col-1)+k) !! error in first ver
        Krigid%Ke(row,col)= Krigid%Ke(row,col) + Rortho%Mp(3*(k-1)+row)*Kcopy(3*(col-1)+k) ! multiply left by transpose of Rortho
      enddo
    enddo
  enddo

END SUBROUTINE MatrixRotMultiply

!**************************************************************!
!							       !
!  Subroutine Matrice_Ke_H_local                               !
!   creates and fills the matrix from 6 ddl to 6 force components!
!**************************************************************!

SUBROUTINE Matrice_Ke_H_local_dfix(MeshT, linktyp_1, linktyp_2, Keii, Keij, Kejj) !MeshT, VarT, Keii, Keij, Kejj)
  use typesbalken
  use printabstractstruct_mod

  implicit none
!  integer,intent(in)           :: iel ! think about passing the whole struct array tMeshElmt,pointer or just an Elmt
  type(tMeshElmt), intent(in)  :: MeshT
!  type(tVarElmt), intent(in)   :: VarT
  integer                      :: linktyp_1, linktyp_2
  type(tRigidMat), intent(inout) :: Keii, Keij, Kejj
!  integer,intent(in)             :: code ! 1 fest eingespannt x2; 2 l-fest eingespannt r-freies Moment, 
!                                         ! 3 l-freies Moment r-fest eingespannt, 4 freies Moment x2
!  integer, dimension(6)         :: connect
! lokalen Variablen
  character(len=30)           :: SubName="Matrice_Ke_H_local_dfix       "
!  type(tRigidMat)             :: Ke_preii, Ke_preij, Ke_prejj
  type(tPassagemat)           :: M_rot

  real                        :: costh, sinth
  integer                     :: i,j, allocstat
  integer,pointer             :: inull(:)
  
  allocate(inull(1), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"Matrice_Ke_H_local_dfix: error allocate null"
  endif
  allocate(M_rot%Mp(9), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"Matrice_Ke_H_local_dfix: error allocate M_rot"
  endif
  do i=1,9
    M_rot%Mp(i)=0.0
  enddo

  costh=(MeshT%endx-MeshT%startx)/MeshT%dlen
  sinth=(MeshT%endy-MeshT%starty)/MeshT%dlen
 !! OLD: 
!  if ((MeshT%typ).eq.TYP_BAR) then ! elmt is a bar, no flexion
!     call matrixrodcondition(MeshT, Keii)
!  elseif ((MeshT%typ).eq.TYP_BEAM) then  ! elmt is a beam with flexion around axis z only.
!     call matrixfixedcondition(MeshT, Keii, 1)
!     call matrixfixedcondition(MeshT, Keij, 2)
!     call matrixfixedcondition(MeshT, Kejj, 3)
!  endif

!  write(*,401) "----- ", adjustl(SubName)
! 401  format (a,a)

  ! actualized implementation: the logic of element type is performed by calling
  ! sub Matrice_Ke_H_Truss and a selection of the element is done here with one more parameter: linktyp
  
  select case(MeshT%typ)
    case(1)
      call matrixrodcondition(MeshT, Keii, 1)
    case(2)
     ! in this case we should distinguish bet. more linkage cases
      if ((linktyp_1.eq.LINK_BEAMEND).and.(linktyp_2.eq.LINK_BEAMEND)) then
        print *, "Beam fixed"
        call matrixfixedcondition(MeshT, Keii, 1)
        call matrixfixedcondition(MeshT, Keij, 2)
        call matrixfixedcondition(MeshT, Kejj, 3)
      else if ((linktyp_1.eq.LINK_BEAMEND).and.(linktyp_2.eq.LINK_BEAMFREE)) then
        print *, "Beam fixed free"
        call matrixfixedfreecondition(MeshT, Keii, 1)
        call matrixfixedfreecondition(MeshT, Keij, 2)
        call matrixfixedfreecondition(MeshT, Kejj, 3)
      else if ((linktyp_1.eq.LINK_BEAMFREE).and.(linktyp_2.eq.LINK_BEAMEND)) then
        print *, "Beam free fixed"
        call matrixfreefixedcondition(MeshT, Keii, 1)
        call matrixfreefixedcondition(MeshT, Keij, 2)
        call matrixfreefixedcondition(MeshT, Kejj, 3)
!      else if ((linktyp_1.eq.LINK_BEAMFREE).and.(linktyp_2.eq.LINKBEAMFREE)) then
      end if
  end select
  
  M_rot%Mp(1)=costh
  M_rot%Mp(2)=sinth
  M_rot%Mp(4)=-sinth
  M_rot%Mp(5)=costh
  M_rot%Mp(9)=1.0
  
  call prntsimplearray(9, 2, inull, M_rot%Mp, SubName)
  print *,"Keii"
  write(*,1002) Keii%Ke(2,2)
 1002  format (f10.5)
 
  call prntsimplestruct(3,3,Keii,SubName)

  call prntsimplestruct(3,3,Keij,SubName)

  call prntsimplestruct(3,3,Kejj,SubName)

  call MatrixRotMultiply(Keii, M_rot, 3)
  call prntsimplestruct(3,3,Keii,SubName)

  call MatrixRotMultiply(Keij, M_rot, 3)
  call prntsimplestruct(3,3,Keij,SubName)

  call MatrixRotMultiply(Kejj, M_rot, 3)
  call prntsimplestruct(3,3,Kejj,SubName)

  deallocate(inull, M_rot%Mp, STAT = allocstat)
  if (allocstat.ne.0) then
    print *,"Matrice_Ke_H_local_dfix: error deallocate M_rot"
  endif
END SUBROUTINE Matrice_Ke_H_local_dfix

!**************************************************************!
! Subroutine Matrice_Ke_H_assembly                             !
!   shall be called after each element submatrix               !
!**************************************************************!
SUBROUTINE Matrice_Ke_K_assembly(Ke_H,Keii,Keij,Kejj,jj,kk,n,ne)

  use typesbalken


  implicit none

  type(tRigidFullMat),intent(inout) :: Ke_H   ! Whole global rigidity matrix
  type(tRigidMat),intent(in)   :: Keii, Keij, Kejj
  integer,intent(in)           :: n,ne   ! Max number of nodes, max number of deg freedom (3 or 6)
  integer,intent(in)           :: jj,kk       ! index of connected nodes to be assemblied in Ke_H
! lokalen Variablen
  integer                     :: i,j,l,p,allocstat
  character(len=30)            :: SubName="Matrice_Ke_H_assembly         "

  do i=1,3
     do j=1,3
       Ke_H%Ke((jj-1)*ne+i, (jj-1)*ne+j) = Ke_H%Ke((jj-1)*ne+i, (jj-1)*ne+j) + Keii%Ke(i,j)
       Ke_H%Ke((jj-1)*ne+i, (kk-1)*ne+j) = Ke_H%Ke((jj-1)*ne+i, (kk-1)*ne+j) + Keij%Ke(i,j)
       Ke_H%Ke((kk-1)*ne+i, (jj-1)*ne+j) = Ke_H%Ke((kk-1)*ne+i, (jj-1)*ne+j) + Keij%Ke(j,i)
       Ke_H%Ke((kk-1)*ne+i, (kk-1)*ne+j) = Ke_H%Ke((kk-1)*ne+i, (kk-1)*ne+j) + Kejj%Ke(i,j)
     enddo
  enddo

END SUBROUTINE  Matrice_Ke_K_assembly

!**************************************************************!
!							       !
!  Subroutine Matrice_Ke_H_Truss                               !
!   shall be called in the main Program as frequent there      !
!   are many elements                                          !
! not terminated                                               !
!**************************************************************!

SUBROUTINE Matrice_Ke_H_Truss(MeshT1, VarT1, Ke_H, jj, kk, nn, ne, n, &
                connectelmt, Dunknowns, Dimposed, Fmovemt, Fapplied)
  use typesbalken
  use printabstractstruct_mod,only : prntadjuststruct

  implicit none

  type(tMeshElmt), intent(in)  :: MeshT1  ! One bar or beam element
  type(tVarElmt), intent(in)   :: VarT1   ! Unknown variables (not used in this subroutine)
                                          ! and applied forces (used)
  type(tRigidFullMat)          :: Ke_H   ! Whole global rigidity matrix
  integer,intent(in)           :: n, nn, ne   ! Current elmt, Max number of nodes, max number of deg freedom (3 or 6)
  integer,intent(in)           :: jj,kk       ! index of connected nodes to be assemblied in Ke_H
  integer,pointer,intent(in)   :: connectelmt(:,:)  ! These indexes (jj,kk) are also found if connectelmt(j,k)=1
  integer,pointer,intent(in)   :: Dunknowns(:,:)    ! Dunkowns(i,j)=1 means j-th of 4 deg of 
                                                    ! freedom of i-th of n nodes is unknown, 0 by default
  real,pointer,intent(in)      :: Dimposed(:,:)  ! imponierte Bewegungen at move elsewhere Dunknowns(i,j)=-1
  integer,pointer,intent(in)   :: Fmovemt(:,:)      ! Fmovemt(i,j)=1 means there one Force to calculate at move j at node i
  real,pointer,intent(in)      :: Fapplied(:,:)  ! Force applied at move j at node i (provided Fmovemt(i,j)=0)
  character(len=30)            :: SubName="Matrice_Ke_H_Truss            "

! lokalen Variablen

  integer                       :: i,j,l,p,n1,n2,allocstat
  integer                       :: linktyp_1, linktyp_2
  type(tRigidMat)               :: Keii, Keij, Kejj

!  allocate(Ke_H%Ke(1:n*ne,1:n*ne), STAT=allocstat) !! BY NO WAY !! (left as comment to not do that again!)
!  if (allocstat.ne.0) then
!    print *,"Matrice_Ke_H_Truss: error allocate global Ke_H"
!  endif
  allocate(Keii%Ke(1:3,1:3), Keij%Ke(1:3,1:3), Kejj%Ke(1:3,1:3), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"Matrice_Ke_H_Truss: error allocate Keii,ij,jj"
  endif
  n1=MeshT1%node_1
  n2=MeshT1%node_2
  ! use of macros to determine which linkage is to handle
  ! case 1: all forces have to be determined
  if ((mcantileverI(n,1)).and.(mcantileverI(n,2)) &
         .and.(mcantileverI(n,3)).and.(mcantileverJ(n,1)) &
         .and.(mcantileverJ(n,2)).and.(mcantileverJ(n,3))) then

     call Matrice_Ke_H_local_dfix(MeshT1, 3,3, Keii, Keij, Kejj)  
     ! tbd/tbc: include case 7 by inquiring the index of beam elmt in a discrete BEAM

     !to complete other cases
     ! case 2: only displacement at one (let's say second, convention) node tbd
  else if ((mcantileverI(n,1)).and.(mcantileverI(n,2)) &
         .and.(mcantileverI(n,3)).and.(marticulatedJ(n))) then
         
     call Matrice_Ke_H_local_dfix(MeshT1, 3,5, Keii, Keij, Kejj) ! 1fixed,2free

  else if ((marticulatedI(n)).and.(marticulatedJ(n))) then
     call Matrice_Ke_H_local_dfix(MeshT1, 1,1, Keii, Keij, Kejj)
  endif   
 
  call Matrice_Ke_K_assembly(Ke_H,Keii,Keij,Kejj,jj,kk,n,ne)

  call prntadjuststruct(nn, nn, Ke_H, SubName)

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

SUBROUTINE Vector_F_RHS_local(MeshT, VarT, linktyp_1, linktyp_2, fx1,fy1,mz1, fx2,fy2,mz2)
  use typesbalken
  implicit none

  type(tMeshElmt), intent(in)  :: MeshT
  type(tVarElmt), intent(inout) :: VarT  ! contains the force at nodes i and j
  integer                      :: linktyp_1, linktyp_2

!  integer,intent(in)              :: code ! 1 fest eingespannt x2; 2 l-fest eingespannt r-freies Moment, 
                                         ! 3 l-freies Moment r-fest eingespannt, 4 freies Moment x2
  real                         :: fx1,fy1,mz1, fx2,fy2,mz2
  real                         :: costh, sinth

  costh=(MeshT%endx-MeshT%startx)/MeshT%dlen
  sinth=(MeshT%endy-MeshT%starty)/MeshT%dlen

  ! TO COMPLETE !
  select case(MeshT%typ)
    case(1) ! bar elmt
      VarT%FeIi(1)=fx1
      VarT%FeIi(2)=fy1
      VarT%FeIi(3)=0.0 ! will not be used
      VarT%FeJi(1)=fx2
      VarT%FeJi(2)=fy1
      VarT%FeJi(3)=0.0 ! will not be used
    case(2) ! beam element
     ! in this case we should distinguish bet. more linkage cases
      if ((linktyp_1.eq.LINK_BEAMEND).and.(linktyp_2.eq.LINK_BEAMEND)) then
        print *, "Beam fixed"
        VarT%UeIi(1)=0.0
        VarT%UeIi(2)=0.0
        VarT%UeIi(3)=0.0
        VarT%UeJi(1)=0.0
        VarT%UeJi(2)=0.0
        VarT%UeJi(3)=0.0
        ! Nothing to do in Fe since displacement imposed 
      elseif ((linktyp_1.eq.LINK_BEAMEND).and.(linktyp_2.eq.LINK_BEAMFREE)) then
        print *, "Beam fixed and free"
        VarT%UeIi(1)=0.0
        VarT%UeIi(2)=0.0
        VarT%UeIi(3)=0.0

        VarT%FeJi(1)=MeshT%q1  / 2. * MeshT%dlen
! Interpolation betw q1 and q2 will be treated after
        VarT%FeJi(2)=MeshT%q1 / 2. * MeshT%dlen
        VarT%FeJi(3)=MeshT%q1 /12. * ( MeshT%dlen )**2

        ! Nothing to do in Fe since displacement imposed 
      elseif ((linktyp_1.eq.LINK_BEAMFREE).and.(linktyp_2.eq.LINK_BEAMEND)) then
        print *, "Beam free and fixed"
       VarT%FeIi(1)=MeshT%q1  / 2. * MeshT%dlen
! Interpolation betw q1 and q2 will be treated after
        VarT%FeIi(2)=MeshT%q1 / 2. * MeshT%dlen
        VarT%FeIi(3)=MeshT%q1 /12. * ( MeshT%dlen )**2

        VarT%UeJi(1)=0.0
        VarT%UeJi(2)=0.0
        VarT%UeJi(3)=0.0 
        ! Nothing to do in Fe since displacement imposed 
      endif
  end select

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



