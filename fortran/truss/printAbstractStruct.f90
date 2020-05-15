module printabstractstruct_mod
  
!**********************************************************************!
! should print various type of 1D array, or particular 2D (1:3,1:n) arrays
! at the end it should be reused for Triangular Mesh of Element
!**********************************************************************!
  use typesbalken
  implicit none

#define FILEOUTPUT 27

  private

!    interface prntsimplearray
!      module procedure prntsimplearray
!    end interface

 !  interface  prntmultisimplearray
 !    module procedure prntmultisimplearray
 !  end interface

   public::  prntsimplearray, prntsimplestruct, prntadjuststruct, &
               printCreateTrussGen

  contains

!*******************************************************************
! subroutine prntsimplearray: if length <7, print as a whole line
!   if length> 10, print 10 float as a line, as much lines are needed
!*******************************************************************
  subroutine prntsimplearray(nlen, typ, ptr_int, ptr_dou, SubName)
  use typesbalken
  implicit none
  ! args
  integer                         :: nlen, typ
  integer,pointer                 :: ptr_int(:)
  real,pointer                    :: ptr_dou(:)
  character(len=30)               :: SubName
  ! Vorhabe
  intent(in)                      :: nlen,typ,ptr_int,ptr_dou,Subname
  ! lokale vars
  integer                         :: i,j,q,r

  write(FILEOUTPUT,201) "----- ", adjustl(SubName)
  if (nlen<=7) then
     if (typ.eq.1) then  ! integer
        write(FILEOUTPUT,202) (ptr_int(j), j=1,nlen)
     elseif (typ.eq.2) then ! real
        write(FILEOUTPUT,203) (ptr_dou(j), j=1,nlen)
     endif
  elseif (nlen>=8) then
     q=nlen/7
     r=nlen-7*q
     do i=1,q
       write(FILEOUTPUT,204) "      ", i
       if (typ.eq.1) then  ! integer
          write(FILEOUTPUT,202) (ptr_int((i-1)*7+j), j=1,7)
       elseif (typ.eq.2) then ! real
          write(FILEOUTPUT,203) (ptr_dou((i-1)*7+j), j=1,7)
       endif

    enddo
    write(FILEOUTPUT,204) "      ", q+1
    if (typ.eq.1) then  ! integer
        write(FILEOUTPUT,202) (ptr_int(q*7+j), j=1,r)
     elseif (typ.eq.2) then ! real
        write(FILEOUTPUT,203) (ptr_dou(q*7+j), j=1,r)
     endif
  endif

 201  format (a,a)
 202  format (i10)
 203  format (f10.5)
 204  format (a,i10)
  end subroutine prntsimplearray

!*******************************************************************
! subroutine prntsimplestruct: a line per s%array values

  subroutine prntsimplestruct(n1, n2, Ms, SubName)
  use typesbalken
  implicit none
  ! args
  integer                         :: n1, n2
  type(tRigidMat)                 :: Ms
  character(len=30)               :: SubName
  ! Vorhabe
  intent(in)                      :: n1,n2,Ms,Subname
  ! lokale vars
  integer                         :: i,j,q,r
  real                            :: scalar

  write(FILEOUTPUT,301) "----- ", adjustl(SubName)

    do i=1,n1
      write(FILEOUTPUT,304) "      ", i
      do j=1,n2 
        scalar = Ms%Ke(i,j)
        write(FILEOUTPUT,303) i,j,scalar
      enddo
    enddo
!  endif

 301  format (a,a)
 302  format (i10)
 303  format (i10,i10,f10.1)
 304  format (a,i10)
  end subroutine prntsimplestruct

!*******************************************************************
! subroutine prntadjuststruct : needs two functions
  integer function indexIJ(j,n2)
  implicit none
  integer, intent(in)   :: j,n2
  integer               :: ij
    ij=j/n2
   indexIJ= ij
  end function indexIJ

  integer function indexJJ(j,n2)
  implicit none
  integer, intent(in)   :: j,n2
  integer               :: ij,jj
   ij=j/n2 
   jj=j-ij*n2
  indexJJ= jj
  end function indexJJ

  subroutine prntadjuststruct(n1, n2, FMs, SubName)
  use typesbalken
  implicit none
  ! args
  integer                         :: n1, n2
  type(tRigidFullMat)             :: FMs
  character(len=30)               :: SubName
  ! Vorhabe
  intent(in)                      :: n1,n2,FMs,Subname
  ! lokale vars
  integer                         :: i,j,q,r,l

  write(FILEOUTPUT,501) "----- ", adjustl(SubName)

  q=n2/7
  r=n2-q*7
  do i=1,n1
  write(FILEOUTPUT,504) "      ", i
    do l=1,q
      write(FILEOUTPUT,507) (FMs%Ke(i,(l-1)*7+j), j=1,7)
    enddo
    write(FILEOUTPUT,507) (FMs%Ke(i,q*7+j), j=1,r)
  enddo

 501  format (a,a)
 504  format (a,i10)
 507  format (7f12.1)
  end subroutine prntadjuststruct

!****************************************************
! subroutine prntinputmatrices, need to print a ne*nd matrix
!   for verification of input data coding the problem
  subroutine prntsimplematrix(nrow,ncol, typ, ptr_int, ptr_dou, SubName)
  use typesbalken
  implicit none
  ! args
  integer                         :: nrow,ncol,typ
  integer,pointer                 :: ptr_int(:,:)
  real,pointer                    :: ptr_dou(:,:)
  character(len=30)               :: SubName
  ! Vorhabe
  intent(in)                      :: nrow,ncol,typ,ptr_int,ptr_dou,Subname
  ! lokale vars
  integer                         :: i,j,k,l,q,r

  write(FILEOUTPUT,601) "----- ", adjustl(SubName)
  if (typ.eq.1) then  ! integer
    do i=1,nrow
       q=ncol/7
       r=ncol-q*7
       do k=1,q
         write(FILEOUTPUT,602) (ptr_int(i, (7*k+l)), l=1,7)  ! j=q*7+k
       enddo
         write(FILEOUTPUT,602) (ptr_int(i, (7*k+l)), l=1,r) 
   enddo

  elseif (typ.eq.2) then
    do i=1,nrow
       q=ncol/5
       r=ncol-q*5
       do k=1,q
         write(FILEOUTPUT,603) (ptr_dou(i, (5*k+l)), l=1,5)  ! j=q*7+k
       enddo
         write(FILEOUTPUT,603) (ptr_dou(i, (5*k+l)), l=1,r) 
   enddo

  endif
 601  format (a,a)
 602  format (i10)
 603  format (f10.5)
 604  format (a,i10)
 605  format (7i10)
 606  format (5f10.5)
 end subroutine prntsimplematrix


 subroutine prntinputmatrices(ne, connectelmt, &
                              Dunknowns, Dimposed, &
                              Fmovemt, Fapplied)

  use typesbalken
  implicit none

  integer,intent(in)           :: ne
  integer,pointer,intent(in)   :: connectelmt(:,:)  ! These indexes (j,k) are also found if connectelmt(j,k)=1
  integer,pointer,intent(in)   :: Dunknowns(:,:)    ! Dunkowns(i,j)=1 means j-th of 4 deg of 
                                                    ! freedom of i-th of n nodes is unknown, 0 by default
  real,pointer,intent(in)      :: Dimposed(:,:)  ! imponierte Bewegungen at move elsewhere Dunknowns(i,j)=-1
  integer,pointer,intent(in)   :: Fmovemt(:,:)      ! Fmovemt(i,j)=1 means there one Force to calculate at move j at node i
  real,pointer,intent(in)      :: Fapplied(:,:)  ! Force applied at move j at node i (provided Fmovemt(i,j)=0)
  ! lokale Variables
  integer, pointer             :: inull(:,:)
  real, pointer                :: fnull(:,:)
  integer                      :: allocstat
  character(len=30)            :: SubName="Input_Truss                   "

  allocate(inull(1,1), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"Input_Truss : error allocate null"
  endif
  allocate(inull(1,1), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"Input_Truss : error allocate null"
  endif
 call prntsimplematrix(ne, 3, 1, connectelmt, fnull, SubName)
 call prntsimplematrix(ne, 3, 1, Dunknowns, fnull, SubName)
 call prntsimplematrix(ne, 3, 2, inull, Dimposed, SubName)
 call prntsimplematrix(ne, 3, 1, Fmovemt, fnull, SubName)
 call prntsimplematrix(ne, 3, 2, inull, Fapplied, SubName)

 end subroutine prntinputmatrices


!************************************************************
! Fills the (in main prgm allocated table of struct tMeshElmt
 subroutine printCreateTrussGen(nn,nt,MeshGen,MeshPoints,meshconnect)
 use typesbalken
 implicit none

 integer                               :: nn,nt
 type(tMeshElmt),pointer               :: MeshGen(:)
 type(tmeshCoord)                      :: MeshPoints
 integer,pointer                       :: meshconnect(:,:)
 intent(in)                            :: nn,nt,MeshPoints,meshconnect,MeshGen
! lokal
 integer                               :: i,j
 character(len=30)     :: SubName="printCreateTrussGen           "

 write(FILEOUTPUT,701) "----- ", adjustl(SubName)
 do i=1,nt
   write(FILEOUTPUT, 701) "   ", ""
   write(FILEOUTPUT,702) " --- Element no. ", i, " with ", MeshGen(i)%nraumx, " nodes"
   if ((MeshGen(i)%typ).eq.1) then  !bar.eq.1)
      write(FILEOUTPUT,701) "  ", "rod element"
   elseif ((MeshGen(i)%typ).eq.2) then  !beam.eq.2)
      write(FILEOUTPUT,701) "  ", "beam element"
   endif
   write(FILEOUTPUT,702) "  Node 1: ", MeshGen(i)%node_1, " Node 2: ", MeshGen(i)%node_2, " "
   write(FILEOUTPUT,701) "  ","Coord x,y"
   write(FILEOUTPUT,703)    MeshGen(i)%startx, MeshGen(i)%endx, &
                   MeshGen(i)%starty, MeshGen(i)%endy   

   write(FILEOUTPUT,704)    MeshGen(i)%SArea, MeshGen(i)%CI
   write(FILEOUTPUT,704)    MeshGen(i)%dlen, MeshGen(i)%anglez
   do j=1,4
     write(FILEOUTPUT,705) "   Polynoms coeff for index: ", j
     write(FILEOUTPUT,703) MeshGen(i)%CoeffsH1(j), MeshGen(i)%CoeffsH2(j), &
                   MeshGen(i)%CoeffsH3(j), MeshGen(i)%CoeffsH4(j)
   enddo
 enddo
 write(FILEOUTPUT, 701) "   ", ""

 701 format(a,a)
 702 format(a,i10,a,i10,a)
 703 format(4f10.5)
 704 format(2f10.5)
 705 format(a,i10)
 end subroutine printCreateTrussGen

end module printabstractstruct_mod

