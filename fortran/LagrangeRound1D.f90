MODULE Lagrange1d_mod

!**********************************************************************!
!     ARRAY (1,n+1) to store the values of the coefficients
!       OR store the n roots and the last value be the leading coefficient
!     SUBROUTINES
!     evalution(x) = pol(1)*x^n + pol(2)*X^(n-1)+..+pol(n+1)
!     expand the roots reduced ARRAY -> see subroutine expand2nom
!     expand the roots reduced ARRAY to find the derivatives at
!       node X,j -> see subroutine calcderivative
!     expand the sum of Lagrange polynoms at nodes X,i to find the interpolation polynom with values Z,i -> see subroutine expandlagrange
!     expand the square of a polynom given its roots to find its coefficient
!        and actualize its degree
!     expand the square of a polynom given its root to calculate 
!        the derivatives at node X,j 
!     
!**********************************************************************!
  use typesround
  implicit none
  
  private
    interface evaluation
      module procedure evaluation
    end interface

    interface expand2nom
      module procedure expand2nom      
    end interface


   public::  evaluation, expand2nom

  contains

!**************************************************************!
!							       !
!  Subroutine evaluation                                       !
!                                                              !
!**************************************************************!
SUBROUTINE evaluation(pol, x, res)
  use typesround
  implicit none
  type(tPolynom), intent(in)  :: pol
  real, intent(in)            :: x
  real, intent(out)           :: res
! lokalen Variablen
  real                        :: term
  integer                     :: i

  term=1.0
  res=0.0
  do i=pol%n+1,1,-1
    res=res + pol%coeffs(i)*term
    term=term*x
  enddo

END SUBROUTINE

!**************************************************************!
!							       !
!  Subroutine expand2nom, product of  two polynoms             !
!        pol degree n *pol degree 1 => deg n+1!                !
!                                                              !
!**************************************************************!
SUBROUTINE expand2nom(pol_1, pol)
  use typesround
  implicit none
  type(tPolynom), intent(in)  :: pol_1
  type(tPolynom), intent(inout):: pol
! lokalen Variablen
  real                        :: sumup
  integer                     :: i
! verification
  write(*,*) '-------enter expand2nom:'
!  write(*,102) (pol%coeffs(i), i=1,pol%n+1)

  do i=pol%n+1,1,-1
    if (i.eq.(pol%n+1)) then
       sumup=pol%coeffs(i)*pol_1%coeffs(2) ! let's try and it will be zero if i>=n+1
       pol%coeffs(i+1)=0.0  ! ensures shifting towards lower expand
    else
      sumup=pol%coeffs(i)*pol_1%coeffs(2) 
    endif
    pol%coeffs(i+1)=pol%coeffs(i+1)*pol_1%coeffs(1)+sumup
  enddo
  pol%coeffs(1)=pol%coeffs(1)*pol_1%coeffs(1) ! highest expand
  pol%n = pol%n+1

! verification
  write(*,102) (pol%coeffs(i), i=1,pol%n+1)
!  write(*,*)
102 format (f10.5)
END SUBROUTINE expand2nom


END MODULE Lagrange1d_mod



