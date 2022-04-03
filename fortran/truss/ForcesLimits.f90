module ForcesLimits_mod
 use TypesBalken
 implicit none

 private

  interface ExtractTrussUnknowns
    module procedure ExtractTrussUnknowns
  end interface

 public :: ExtractTrussUnknowns, ExtractAblock

 contains

!***************************************
! Help subroutine for extracting a 3x3 or 2x2 block

 subroutine ExtractAblock(Ke_H, V_H, A_H, B_H, nk, na, code)
 use TypesBalken
 type(tRigidFullMat)               :: Ke_H
 type(tVarFull)                    :: V_H
 type(tRigidFullMat)               :: A_H
 type(tVarFull)                    :: B_H
 integer                           :: na ! current position in matrixA 
! integer                           :: nb ! same in vector B
 integer                           :: nk ! current position from matrix K
 integer                           :: code ! =1.copy ux corresp coeffs 
                                           ! =2.copy uy corresp coeffs 
                                           ! =3.copy Oz corresp coeffs
                                          ! =4.copy ux+uy+Oz corresp coeff.
 ! Vorhabe
 intent(in)                        :: Ke_H, V_H, na, nk, code
 intent(inout)                     :: A_H, B_H
 ! lokale
 integer                           :: i,j
 if (code==4) then
  do i=1,3
    do j=1,3
       A_H%Ke(na+i, na+j)=Ke_H%Ke(nk+i,nk+j)
    enddo
    B_H%Fe(na+i)=V_H%Fe(nk+i)
  enddo
 elseif (code==1) then ! copy ux corresp coeff 
    A_H%Ke(na+1, na+1)=Ke_H%Ke(nk+1,nk+1)
    B_H%Fe(na+1)=V_H%Fe(nk+1)
 elseif (code==2) then ! copy ux corresp coeff 
    A_H%Ke(na+1, na+1)=Ke_H%Ke(nk+2,nk+2)
    B_H%Fe(na+1)=V_H%Fe(nk+2)
 elseif (code==3) then ! copy ux corresp coeff 
    A_H%Ke(na+1, na+1)=Ke_H%Ke(nk+3,nk+3)
    B_H%Fe(na+1)=V_H%Fe(na+3)
 endif
 end subroutine ExtractAblock


!****************************************
! Sub ExtractTrussUnknowns
!   accepts a whole rigidity matrix and copies
!   only a smaller submatrix for the unknown
!   displacement (columns->rows in Fortran) 
!   and knowns forces (rows->columns in Fortran)
 subroutine ExtractTrussUnknowns(MeshInfo, meshpattern, Dunknowns, Fmovemt, Ke_H, V_H, A_H, B_H, IndexesKe)
 use TypesBalken
 implicit none
 type(tMeshInfo)                   :: MeshInfo
 integer,pointer                   :: meshpattern(:,:) 
 integer,pointer                   :: Dunknowns(:,:)     ! 1: disp imposed, 0: disp to determine
 integer,pointer                   :: Fmovemt(:,:)       ! why this, if symetrical to Dunknonwns?
 type(tRigidFullMat)               :: Ke_H, A_H
 type(tVarFull)                    :: V_H, B_H
 integer,pointer                   :: IndexesKe(:)
 ! Vorhabe
 
 intent(in)                        :: Dunknowns,Fmovemt,Ke_H
 intent(inout)                     :: MeshInfo,A_H, B_H
 ! lokale
 integer                           :: nn,nt,la,na,nb,nk,i,j,n1,n2,ni,allocStat
 
 nt=MeshInfo%nt  ! nof element
 nn=MeshInfo%nn  ! nof nodes
 la=0  ! dimension of extracted matrix A_H and vector B_H
 na=0  ! index in A_H
 nb=0  ! index in B_H
 nk=0 ! position index in Ke_H
 allocate(IndexesKe(1:nn), stat=allocStat)
 if (allocStat.ne.0) then
    write(*,702) " -----ExtractTrussUnknowns error, cannot allocate IndexesKe"
 endif
 do i=1,nn
   IndexesKe(i)=-1  ! Initialization, -1 is a code to denote no belonging to pb end unknowns in the linear sys A_H*X=B_H
 enddo
! First loop to determine the dimension of matrices A_H, B_H
! do i=1,nn  
!    if ((Dunknowns(i,1).eq.1).and.(Dunknowns(i,2).eq.1).and.(Dunknowns(i,3).eq.1)) then
!      ! nothing
!      la=la+3
!    elseif ((Dunknowns(i,1).eq.0).and.(Dunknowns(i,2).eq.1).and.(Dunknowns(i,3).eq.1)) then ! Ux dof
!      la=la+1
!    elseif ((Dunknowns(i,1).eq.1).and.(Dunknowns(i,2).eq.0).and.(Dunknowns(i,3).eq.1)) then ! Uy dof
!      la=la+1
!    elseif ((Dunknowns(i,1).eq.1).and.(Dunknowns(i,2).eq.1).and.(Dunknowns(i,3).eq.0)) then ! Oz dof
!      la=la+1
!    elseif ((Dunknowns(i,1).eq.0).and.(Dunknowns(i,2).eq.0).and.(Dunknowns(i,3).eq.0)) then ! all dof
!      la=la+1
!    endif
! enddo
 do j=1,nt
   n1=meshpattern(j,2)
   n2=meshpattern(j,3)
   do i=1,nn
    if (((i.eq.n1).or.(i.eq.n2)).and.(IndexesKe(i).eq.(-1))) then 
      if (i.eq.n1) then
         ni=0
      elseif (i.eq.n2) then
         ni=3
      endif
      if ((Dunknowns(j,ni+1).eq.1).and.(Dunknowns(j,ni+2).eq.1).and.(Dunknowns(j,ni+3).eq.1)) then
         ! nothing
     elseif ((Dunknowns(j,ni+1).eq.0).and.(Dunknowns(j,ni+2).eq.1).and.(Dunknowns(j,ni+3).eq.1)) then ! Ux dof
        IndexesKe(i)=la+1
        la=la+1
      elseif ((Dunknowns(j,ni+1).eq.1).and.(Dunknowns(j,ni+2).eq.0).and.(Dunknowns(j,ni+3).eq.1)) then ! Uy dof
        IndexesKe(i)=la+1
        la=la+1
      elseif ((Dunknowns(j,ni+1).eq.1).and.(Dunknowns(j,ni+2).eq.1).and.(Dunknowns(j,ni+3).eq.0)) then ! Oz dof
        IndexesKe(i)=la+1
        la=la+1
      elseif ((Dunknowns(j,ni+1).eq.0).and.(Dunknowns(j,ni+2).eq.0).and.(Dunknowns(j,ni+3).eq.0)) then ! all dof
         IndexesKe(i)=la+1  ! index from whom rest of prgm will refer to have the calculated dof in extracted matrix
        la=la+3
      endif

    endif
   enddo
 enddo 

 write(*,702) " -----ExtractTrussUnknowns: counting of A_H : ", la
 MeshInfo%la = la
 allocate(A_H%Ke(1:la,1:la), B_H%Ue(1:la), B_H%Fe(1:la), stat=allocStat)
 if (allocStat.ne.0) then
    write(*,702) " -----ExtractTrussUnknowns error, cannot allocate A_H, B_H"
 endif
 
! Second loop to extract the real values from the rigidity matrix that corresponds to the real unknowns of the problem

  do j=1,nt
   n1=meshpattern(j,2)
   n2=meshpattern(j,3)
   do i=1,nn
    if (((i.eq.n1).or.(i.eq.n2)).and.(IndexesKe(i).eq.(-1))) then 
      if (i.eq.n1) then
         ni=0
      elseif (i.eq.n2) then
         ni=3
      endif
      if ((Dunknowns(i,ni+1).eq.1).and.(Dunknowns(i,ni+2).eq.1).and.(Dunknowns(i,ni+3).eq.1)) then
            ! nothing to do but
       nk=nk+3
      elseif ((Dunknowns(i,ni+1).eq.0).and.(Dunknowns(i,ni+2).eq.1).and.(Dunknowns(i,ni+3).eq.1)) then ! Ux dof
       call ExtractAblock(Ke_H,V_H,A_H,B_H,nk,na,1)
       na=na+1
       nk=nk+3
      elseif ((Dunknowns(i,ni+1).eq.1).and.(Dunknowns(i,ni+2).eq.0).and.(Dunknowns(i,ni+3).eq.1)) then ! Uy dof
       call ExtractAblock(Ke_H,V_H,A_H,B_H,nk,na,2)
       na=na+1
       nk=nk+3
      elseif ((Dunknowns(i,ni+1).eq.1).and.(Dunknowns(i,ni+2).eq.1).and.(Dunknowns(i,ni+3).eq.0)) then ! Oz dof
       call ExtractAblock(Ke_H,V_H,A_H,B_H,nk,na,3)
       na=na+1
       nk=nk+3
      elseif ((Dunknowns(i,ni+1).eq.0).and.(Dunknowns(i,ni+2).eq.0).and.(Dunknowns(i,ni+3).eq.0)) then ! all dof
       call ExtractAblock(Ke_H,V_H,A_H,B_H,nk,na,4)
       na=na+1
       nk=nk+3
      endif
    endif
   enddo 
  enddo
702 format(a,i10)
 end subroutine ExtractTrussUnknowns


end module ForcesLimits_mod

