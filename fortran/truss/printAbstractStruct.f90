module printabstractstruct_mod
  
!**********************************************************************!
! should print various type of 1D array, or particular 2D (1:3,1:n) arrays
! at the end it should be reused for Triangular Mesh of Element
!**********************************************************************!
  use typesbalken
  implicit none

  private

!    interface prntsimplearray
!      module procedure prntsimplearray
!    end interface

 !  interface  prntmultisimplearray
 !    module procedure prntmultisimplearray
 !  end interface

   public::  prntsimplearray, prntsimplestruct, prntadjuststruct

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

  write(*,201) "----- ", adjustl(SubName)
  if (nlen<=7) then
     if (typ.eq.1) then  ! integer
        write(*,202) (ptr_int(j), j=1,nlen)
     elseif (typ.eq.2) then ! real
        write(*,203) (ptr_dou(j), j=1,nlen)
     endif
  elseif (nlen>=8) then
     q=nlen/7
     r=nlen-7*q
     do i=1,q
       write(*,204) "      ", i
       if (typ.eq.1) then  ! integer
          write(*,202) (ptr_int((i-1)*7+j), j=1,7)
       elseif (typ.eq.2) then ! real
          write(*,203) (ptr_dou((i-1)*7+j), j=1,7)
       endif

    enddo
    write(*,204) "      ", q+1
    if (typ.eq.1) then  ! integer
        write(*,202) (ptr_int(q*7+j), j=1,r)
     elseif (typ.eq.2) then ! real
        write(*,203) (ptr_dou(q*7+j), j=1,r)
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

  write(*,301) "----- ", adjustl(SubName)

    do i=1,n1
      write(*,304) "      ", i
      do j=1,n2 
        scalar = Ms%Ke(i,j)
        write(*,303) i,j,scalar
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

  integer 	function indexJJ(j,n2)
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

  write(*,501) "----- ", adjustl(SubName)

  q=n2/7
  r=n2-q*7
  do i=1,n1
  write(*,504) "      ", i
    do l=1,q
      write(*,507) (FMs%Ke(i,(l-1)*7+j), j=1,7)
    enddo
    write(*,507) (FMs%Ke(i,q*7+j), j=1,r)
  enddo

 501  format (a,a)
 504  format (a,i10)
 507  format (7f12.1)
  end subroutine prntadjuststruct

end module printabstractstruct_mod

