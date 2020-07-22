! Interpolation means the calculation of displacements and rotations
! through the Element. It can be have a support of 6 nodes triangles
! (membrane effect with 3 integration points) ((or 3 nodes triangle))
!
! Use: interpolate with a0+a1*X+a2*Y+a3*X*Y
module LagrangeRound2D_mod
!---
implicit none
 private 

 real,dimension(18)      :: K_reference=(/ (/ 0.5, 0.0, 0.5, 0.0, 0.0, 0.0 /), &
                                           (/ 0.0, 0.5, 0.0, 0.0, 0.5, 0.5 /), &
                                           (/ 0.5, 0.5, 0.0, 0.5, 0.0, 0.0 /) /)

 interface multmatrixAcB
    module procedure multmatrixAcB
  end interface

  interface linprntmatrixMem
    module procedure linprntmatrixMem
  end interface

  interface linlgcreatematrixMem
    module procedure linlgcreatematrixMem
  end interface

  interface linlgintegratematrixMem
    module procedure linlgintegratematrixMem
  end interface

public :: multmatrixAcB, linprntmatrixMem, linlgcreatematrixMem, linlgintegratematrixMem, linlgassemblyPartAMem


contains 
#define ELnum_j_Xcoord(n,j)	 MeshR%nodes( MeshR%elems(n)%numeros(j+1) )%vects(1)
#define ELnum_j_Ycoord(n,j)	 MeshR%nodes( MeshR%elems(n)%numeros(j+1) )%vects(2)
#define ELnum_j_Zcoord(n,j)	 MeshR%nodes( MeshR%elems(n)%numeros(j+1) )%vects(3)

!***********************************************
! Matrix-Product very simple function
  subroutine multmatrixAcB(i,k,j,A,B,C)
  implicit none
  integer            :: i,k,j
  real,pointer       :: A(:,:),B(:,:)
  real, pointer      :: C(:,:)
  intent(in)         :: i,k,j,A,B
  intent(inout)      :: C
  
  integer            :: l,m,n
  real               :: sumnum

  do l=1,i
    do n=1,j
       sumnum=0.0
       do m=1,k   ! A(l,m) size{i,k} * B(m,n) size{k,j}
 !        sumnum = sumnum+A(i*(m-1)+l)*B(k*(n-1)+m)  ! m number of indiv row in A, k size of rows in B, n number of cols in B
          sumnum = sumnum+A(l,m)*B(m,n)
       enddo 
       C(l,n)=sumnum  ! C(l,n) C(i*(n-1)+l) size{i,j}
    enddo

  enddo
 
  end subroutine multmatrixAcB

!***********************************************
! Print a Matrix
  subroutine linprntmatrixMem(K_integrate) 

  use TypesRound
  implicit none
  ! uebergegebene
  real,pointer        :: K_integrate(:,:)    !  dimension(36) 
  intent(in)                 :: K_integrate
  integer i,j
  print *,"  -------- linprntmatrixMem:"
  do i=1,6
    write(*,101) "!--------!--------!--------!--------!--------!---------!"
    write(*,301) "!", (K_integrate(i,j), j=1,6)
  enddo
    write(*,101) "!--------!--------!--------!--------!--------!---------!"
101 format (a,a)
201 format (a,i10)
301 format (a,6f10.5,a)
  end subroutine linprntmatrixMem


!***********************************************
! create a matrix for membrane stress "binded" to an element
! from determination of triplet and passage matrix for
! a single trinagle
  subroutine linlgcreatematrixMem(NTestMeshR,MeshR,ni,K_eval,Pmat_xyz)
!------!------!------!------!------!------!
! y3   ! 0    !  y3  !  0   !  0   !  0   !
!------!------!------!------!------!------!
!  0   ! x2-x3!  0   !  0   ! x2   !  0   ! 
!------!------!------!------!------!------!
! x2-x3!  y3  !  0   !  y3  !  0   !  x2  ! 
!------!------!------!------!------!------!

  use TypesRound
  use verifymeshround_mod   ! important module to calculate the geometry of an elementar triangle elmt given its points
  implicit none

  type(tRoundMeshInfo) :: NTestMeshR
  type(tRoundMesh)     :: MeshR
  integer              :: ni
  type(tRoundCoeff)    :: K_eval
  type(tPassageMatrix) :: Pmat_xyz
  intent(in)           :: NTestMeshR, MeshR, ni
  intent(out)          :: K_eval, Pmat_xyz

! lokale Variablen for element local basis
  integer              :: i,j,allocStat
  type(tRoundlocal)    :: Pcharac
  type(tNode)              :: ux, uy, uz
  real,dimension(1:3)      :: dist
  real,dimension(1:3)      :: angle
  real,dimension(1:3)      :: triplet
  real,dimension(1:4)      :: area
  ! end vars
 
  allocate( ux%vects(3), uy%vects(3), uz%vects(3), &
          STAT=allocStat );  ! triangle elmt and only 6 d.o.f. (among total of 15) represents membrane effect
                                                 ! and x and y delta values in each of them
  if (allocStat.ne.0) then
     print *, "linlgcreatematrixMem allocation error, cannot allocate K_eval"
  endif
  ! copy of given element coordinates
  do j=1,3
     ux%vects(j) = ELnum_j_Xcoord(ni,j)  ! preprocessor macro
     uy%vects(j) = ELnum_j_Ycoord(ni,j)
     uz%vects(j) = ELnum_j_Zcoord(ni,j)
  enddo
  ! calculate elements charac in local basis
  call trigcharacteristics(ux,uy,uz, dist, angle, triplet, area)
  ! calculate passage matrix and recopy
  call trigcalcbasis(ux,uy,uz,Pmat_xyz)
  ! recopy result
  do j=1,3
    Pcharac%triplet(j)=triplet(j)
    do i=1,3
      Pcharac%axes_xyz%mat(3*(j-1)+i)=Pmat_xyz%mat(3*(j-1)+i)
    enddo
   !!! Pcharac%angleEuler(j)
  enddo

  ! fill coeffs=difference of local coord for !multiply! normalized reference elmt:
  ! this is a awkyard mean to program genralized element from a conform triangle rectangle
  ! (SEE the constant table K_reference at the top of this module) !

  K_eval%Coefficients(1,1) = triplet(3)                          !  Ux,%x(1) dof coeff corresp to y3
  K_eval%Coefficients(2,1) = 1.0                                 !  Uy,%y(1) dof coeff will mult zero value..
  K_eval%Coefficients(3,1) = triplet(1)-triplet(2)               !  Uy,%x(1) dof coeff correspond to x2-x3

  K_eval%Coefficients(1,2) = 1.0                                 !  Ux,%x(1) dof coeff will mult zero value..
  K_eval%Coefficients(2,2) = triplet(1)-triplet(2)               !  Uy,%y(1) dof coeff correspond to x2-x3 
  K_eval%Coefficients(3,2) = triplet(3)                          !  Ux,%y(1) dof coeff correspond to y3

  K_eval%Coefficients(1,3) = triplet(3)                          !  Ux,%x(2) dof coeff corresp to y3
  K_eval%Coefficients(2,3) = 1.0                                 !  Uy,%y(2) dof coeff will mult zero value..
  K_eval%Coefficients(3,3) = 1.0                                 !  Uy,%x(2) dof coeff will mult zero value..  

  K_eval%Coefficients(1,4) = 1.0                                 !  Uy,%x(2) dof coeff will mult zero value
  K_eval%Coefficients(2,4) = 1.0                                 !  Uy,%y(2) dof coeff will mult zero value
  K_eval%Coefficients(3,4) = triplet(3)                          !  Ux,%y(2) dof coeff correspond to y3

  K_eval%Coefficients(1,5) = 1.0                                 !  Ux,%x(3) dof coeff will mult zero value..
  K_eval%Coefficients(2,5) = 1.0                                 !  Uy,%y(3) dof coeff will mult zero value
  K_eval%Coefficients(3,5) = triplet(1)   ! ERROR HERE (see (2,6)                       !  Uy,%x(3) dof coeff correspond to x2

  K_eval%Coefficients(1,6) = 1.0                                 !  Uy,%x(3) dof coeff will mult zero value
  K_eval%Coefficients(2,6) = triplet(1)  ! ERROR HERE (see (3,5)                        !  Uy,%y(3) dof coeff correspond to x2..
  K_eval%Coefficients(3,6) = 1.0                                 !  Ux,%y(3) dof coeff will mult zero value

! deallocate nodes only
   deallocate(ux%vects, uy%vects, uz%vects, STAT=allocStat );
  if (allocStat.ne.0) then
     print *, " linlgcreatematrixMem deallocation error, cannot deallocate ux,uy,yz"
  endif


  end subroutine linlgcreatematrixMem

!***********************************************
! follows latest subroutine, and evaluate a matrix for membrane stress
! to integrate it
  subroutine linlgintegratematrixMem(K_membrane,x,y,K_eval,Pmat_xyz)

  use TypesRound
  implicit none
  type(tRoundCoeff)                 :: K_eval
  real,pointer                      :: K_membrane(:,:)
  real                              :: x,y
  type(tPassageMatrix)              :: Pmat_xyz
  intent(in)                        :: Pmat_xyz ,      x,y
  intent(inout)                     :: K_eval, K_membrane
 
! lokale Variable
  real,pointer                      :: K_trans(:,:), D_sig(:,:), K_interm(:,:), K_interm2(:,:)
  integer                           :: i,j,allocStat
  real                              :: sumnum, vu, EY
  vu=0.25
  EY=210000

  allocate(K_membrane(1:6,1:6), K_trans(1:6,1:3), D_sig(1:3,1:3), K_interm(1:6,1:3), &
          K_interm2(1:3,1:6), STAT=allocStat ); 
  if (allocStat.ne.0) then
     print *, "linlgintegratematrixMem allocation error, cannot allocate K_membrane, K_trans, D_sig"
  endif
  ! fills matrix of plane stress behaviour Dm
  D_sig(1,1)=EY/(1-vu*vu)*1
  D_sig(2,1)=EY/(1-vu*vu)*(vu) !! +nu
  D_sig(3,1)=0.0
  D_sig(1,3)=0.0
  D_sig(1,2)=EY/(1-vu*vu)*(vu)  !! +nu
  D_sig(2,2)=EY/(1-vu*vu)*1
  D_sig(2,3)=0.0
  D_sig(3,2)=0.0
  D_sig(3,3)=EY/(1-vu*vu)*(1-vu)/2

  ! fills local matrices trans(Bm) (K_trans, used in d=Bm*Ve) by means of a reference pattern and local coeffs (K_eval)
  !                and Bm itself (K_interm2)
  do i=1,3
    do j=1,6
      K_trans(j,i) = K_eval%Coefficients(i,j)*K_reference(3*(j-1)+i) ! i row minor
      K_interm2(i,j) = K_eval%Coefficients(i,j)*K_reference(3*(j-1)+i)
    enddo
  enddo
  call multmatrixAcB(6,3,3,K_trans,D_sig, K_interm)  ! result into K_interm: Dm*d= Dm*Bm
  call multmatrixAcB(3,3,6,K_interm,K_interm2,K_membrane)   ! End membrane matrix tBm*Dm*Bm

  deallocate(K_trans, D_sig, K_interm, STAT=allocStat );
  if (allocStat.ne.0) then
     print *, " linlgintegratematrixMem deallocation error, cannot deallocate K_eval, K_trans, D_sig"
  endif
  end subroutine linlgintegratematrixMem

!********************************************************
! follows the last subroutine in order to begin to assemble the rigitity matrix
! (membrane rigitity terms
 subroutine linlgassemblyPartAMem(NTestMeshR,K_mem,ie,iabooleans,iplacetoja)

  use TypesRound
  implicit none

  type(tRoundMeshInfo)              :: NTestMeshR
  real,pointer                      :: K_mem(:,:)
  integer                           :: ie
  type(tSparseYMat)                 :: RigidYMat
  integer,pointer                   :: iabooleans(:) ! dimension nn
  integer,pointer                   :: iplacetoja(:,:) ! dimension nn,1..2 two data foreach elmt
  intent(inout)                     :: iabooleans ! pass to modify it
  intent(inout)                     :: iplacetoja  ! has only be once created
  integer                           :: i,k,l,nn

  nn=NTestMeshR%nn
!  from iaboolean(i)=true(1) if and only if there already exists at least an integer at RigidYMat%ja(i) -> will have to update it
!  if there does not exist an integer at RigidYMat%ja(i) -> create it and for this set iabooleans(i)=new(2) (temporarly)

! RigidYMat%a(i) = nonzeros values as simple array
! RigidYMat%ia(i) = number of nnz values in each line, dim ne
! RigidYMat%ja(i) = index of columns of nnz values

 end subroutine linlgassemblyPartAMem

!********************************************************
! in this subroutine, the memory matrix will be really assemblied
 subroutine linlgassemblyPartBMem(NTestMeshR,K_mem,ie,RigidYMat,iplacetoja)

  use TypesRound
  implicit none

  type(tRoundMeshInfo) :: NTestMeshR
  integer,pointer                   :: iplacetoja(:,:) ! dimension 10*ne,1..2 two data foreach elmt
  real,pointer                      :: K_mem(:,:)
  integer                           :: ie
  type(tSparseYMat)                 :: RigidYMat
  integer,pointer                   :: iabooleans(:) ! dimension ne 
  integer,pointer                   :: iplacetoja(:,:) ! dimension ne,1..2 two data foreach elmt
  integer                           :: i,k,l,nn,ne,ibool,nbnew,nbentries,nbtotal

  intent(in)                     :: iabooleans,iplacetoja

  nn=NTestMeshR%nn
  ne=NTestMeshR%ne
  i=0
  nbnew=0
  nbtotal=0
  ibool=ibooleans(i)
! count new elements and elements alreay in (nbtotal)
  do while (i.le.nn)   ! there exists at least an integer at RigidYMat%ja(j(i)) j=1..RigidYMat%ia(i)
    if (ibool.eq.2) then
     nbnew=nbnew+1 
    elseif (ibool.ge.1) then
     nbtotal=nbtotal+1 
    endif
    i=i+1
    ibool=ibooleans(i)
  enddo
  nbentries=6*6  ! the size of an elementar matrix
! at this point we know that we will have to insert max nbnew*6 (pts connected) and consequently pushing the other existing values in RigidYMat%ja(:), how many of 1..6 that is the question

  ! pseudo code:
! if (iboolean(ie)=1)
  ! i=0
  ! while i<=36
      ! k=index line of (i)
      ! l=index line of (i)
      ! 
  ! end    
! else if (iboolean(ie)=2)



! endif
 end subroutine linlgassemblyPartBMem

end module LagrangeRound2D_mod

