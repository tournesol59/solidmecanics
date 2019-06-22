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
  use types
  implicit none
  
  private
    interface evaluation
      module procedure evaluation
    end interface

    interface expand2nom
      module procedure expand2nom      
    end interface

    interface calcderivative
      module procedure calcderivative      
    end interface

    interface expandfromroots
      module procedure expandfromroots
    end interface 

    interface expandlagrange
      module procedure expandlagrange      
    end interface

   interface expandsquare
      module procedure expandsquare      
    end interface

    interface squarecalcderivative
      module procedure squarecalcderivative      
    end interface


   public::  evaluation, expand2nom, calcderivative, expandfromroots, &
             expandlagrange, expandsquare, squarecalcderivative

  contains

!**************************************************************!
!							       !
!  Subroutine evaluation                                       !
!                                                              !
!**************************************************************!
SUBROUTINE evaluation(pol, x, res)
  use types
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
  use types
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

!**************************************************************!
!							       !
!  Subroutine calcderivative: calc the deriv VALUE at point x  !
!                                                              !
!**************************************************************!
SUBROUTINE calcderivative(pol, x, der)
  use types
  implicit none
  type(tPolynom), intent(in)  :: pol
  real, intent(in)            :: x
  real, intent(out)           :: der
! lokalen Variablen
  real                        :: term
  integer                     :: i

  term=1.0
  der=0.0
  do i=pol%n,1,-1
    der=der + (pol%n-i+1)*pol%coeffs(i)*term
    term=term*x
  enddo
END SUBROUTINE calcderivative

!**************************************************************!
!							       !
!  Subroutine  expandfromroots:given the n roots of a polynom  !
!              determinate its coefficients and degree         !
!**************************************************************!
SUBROUTINE expandfromroots(pol_roots, pol_nom)
  use types
  implicit none
  type(tPolynom), intent(in)   :: pol_roots
  type(tPolynom), intent(inout):: pol_nom
! lokalen Variablen
  real                        :: sumup
  integer                     :: i
  type(tPolynom)              :: pol_1
! verification
  write(*,*) '-------enter expandfromroots:'
!  write(*,113) pol_roots%n
  write(*,103) (pol_roots%coeffs(i), i=1,pol_roots%n+1) !(n+1)is not a root,but leading coeff
  pol_nom%n=1                          !  Initialization to a first order polynom
  pol_nom%coeffs(2)= - pol_roots%coeffs(1)  ! minus size
  pol_nom%coeffs(1)= 1.0
  do i=3,pol_nom%n+2  ! rest iz zero
     pol_nom%coeffs(i)=0.0
  enddo
  do i=2,pol_roots%n
    pol_1%n=1
    pol_1%coeffs(1)=1.0
    pol_1%coeffs(2)= - pol_roots%coeffs(i)   ! pol_1: (X-xi)
    call expand2nom(pol_1, pol_nom)
  enddo

  do i=1,pol_roots%n+1  ! multiplikation per Major Term
    pol_nom%coeffs(i)= pol_nom%coeffs(i)*pol_roots%coeffs(pol_roots%n+1)
  enddo

! verification
  write(*,103) (pol_nom%coeffs(i), i=1,pol_nom%n+1)
103 format (f10.5)
113 format (i10)
END SUBROUTINE expandfromroots

!**************************************************************!
!							       !
!  Subroutine   expandlagrange: expand a lagrange polynom, given !
!               the points of interp, pol value is unused here !
!**************************************************************!
SUBROUTINE expandlagrange(pol_pts, pol_value, pol_lagr)
  use types
  implicit none
  type(tPolynom), intent(in)   :: pol_pts, pol_value
  type(tPolynom), intent(inout):: pol_lagr
! lokalen Variablen
  real                        :: term
  integer                     :: i,j
  type(tPolynom)              :: pol_roots, pol_nom
! verification
  write(*,104) (pol_pts%coeffs(i), i=1,pol_pts%n)
  write(*,104) (pol_value%coeffs(i), i=1,pol_value%n)

! Initialization
  pol_lagr%n=pol_value%n-1   ! degree of inteprolation = number of points -1
  do i=1,pol_lagr%n
     pol_lagr%coeffs(i)=0.0
  enddo
  pol_nom%n=pol_value%n-1
  do i=1,pol_nom%n
     pol_nom%coeffs(i)=0.0
  enddo

! Calculus of each Lagrange Polynom corresp. to a node
  do i=1,pol_value%n
     pol_roots%n=pol_value%n-1
      write(*,104) pol_value%coeffs(i)
     term=1.0
     do j=1,pol_value%n   ! n is important to review all points, even if n-1 are taken
       if (j.le.(i-1)) then
         pol_roots%coeffs(j)= - pol_pts%coeffs(j)
         term=term/(pol_pts%coeffs(i)-pol_pts%coeffs(j)) 
       elseif ((i+1).le.j) then
         pol_roots%coeffs(j-1)= - pol_pts%coeffs(j)
         term=term/(pol_pts%coeffs(i)-pol_pts%coeffs(j)) 
       else
      !! nothing 
       endif
     enddo
     term=term
     pol_roots%coeffs(pol_roots%n+1)=term !!!
     ! now expand pol_root and sum,up to lagrange polynom each coeff
     call expandfromroots(pol_roots, pol_nom)
     do j=1,pol_nom%n+1
         pol_lagr%coeffs(j) = pol_lagr%coeffs(j)+pol_nom%coeffs(j)*pol_value%coeffs(i)  ! the sum of all Lagrange polynoms mult by their value i
     enddo   
  enddo

! verification
  write(*,104) (pol_lagr%coeffs(i), i=1,pol_lagr%n+1)
  write(*,*)'------ end expand -------'
104 format (f10.5)
END SUBROUTINE expandlagrange
  

!**************************************************************!
!							       !
!  Subroutine   expandsquare: as lagrange, but square of each  !
!                base polynom (X-xi)^2 , it is used in Hermite interp !
!**************************************************************!
SUBROUTINE expandsquare(pol_pts, i, polsq)
  use types
  implicit none
  type(tPolynom), intent(in)   :: pol_pts
  integer, intent(in)          :: i
  type(tPolynom), intent(inout):: polsq
! lokalen Variablen
  real                        :: term
  integer                     :: j
  type(tPolynom)              :: pol_roots
! verification
  write(*,*) '-------enter expandsquare:'
  write(*,105) (pol_pts%coeffs(j), j=1,pol_pts%n)

! Calculus of each Lagrange Polynom corresp. to a node and take the square of it
     pol_roots%n=pol_pts%n-1
     term=1.0
     do j=1,pol_pts%n
       if (j.le.(i-1)) then
         pol_roots%coeffs(2*j-1)= pol_pts%coeffs(j)
         pol_roots%coeffs(2*j)= pol_pts%coeffs(j)
         term=term/(pol_pts%coeffs(i)-pol_pts%coeffs(j))/(pol_pts%coeffs(i)-pol_pts%coeffs(j)) 
       elseif ((i+1).le.j) then
         pol_roots%coeffs(2*j-3)= pol_pts%coeffs(j)
         pol_roots%coeffs(2*j-2)= pol_pts%coeffs(j)
         term=term/(pol_pts%coeffs(i)-pol_pts%coeffs(j))/(pol_pts%coeffs(i)-pol_pts%coeffs(j)) 
       else
      !! nothing
       endif
     enddo
     pol_roots%n=2*pol_roots%n  ! now that the loop is terminated update the degree
     pol_roots%coeffs(pol_roots%n+1)=term
     ! now expand pol_root 
     call expandfromroots(pol_roots, polsq)
! verification
  write(*,105) (polsq%coeffs(j), j=1,polsq%n+1)
  write(*,*) '-------quit expandsquare:'
105 format (f10.5)
END SUBROUTINE expandsquare


!**************************************************************!
!							       !
!  Subroutine   squarecalcderivative: as the previous expandsquare !
!               but this time, only the value of the derivative!
!**************************************************************!
SUBROUTINE squarecalcderivative(pol_pts, pol_roots, i, der)
  use types
  implicit none
  type(tPolynom), intent(in)   :: pol_pts, pol_roots
  integer, intent(in)          :: i
  real, intent(out)            :: der
! lokalen Variablen
  real                        :: term
  integer                     :: j

! verification
  write(*,106) (pol_pts%coeffs(j), j=1,pol_pts%n)

  term=0.0
  do j=1,pol_roots%n
     if (j.le.(i-1)) then
        term=term + term/(pol_pts%coeffs(i)-pol_pts%coeffs(j))
       elseif ((i+1).le.j) then
        term=term + term/(pol_pts%coeffs(i)-pol_pts%coeffs(j))
       else
        ! do nothing
     endif
  enddo
  der = term
  write(*,106) der

106 format (f10.5)
END SUBROUTINE squarecalcderivative

END MODULE Lagrange1d_mod



