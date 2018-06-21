MODULE hermite1d_mod
!***************************************************************!
!                                                               !
! Calculates the Hermite one-dimensional ( in one interval)     !
!    interpolation polynoms and subroutine to calculate         !
!    its derivative and so that it could also be used for 2d    !
!    The Lagrange1d module will also be called for intermediary !
!    subroutines because:
!     polynom = ARRAY(1,n+1) where n=degree
!     Caution: degree 3 will be used only here
!
!***************************************************************!

use types
use Lagrange1d_mod
implicit none

private

  interface expandhermite
     module procedure expandhermite
  end interface

  interface quadhermitecub00
     module procedure quadhermitecub00
  end interface

public:: expandhermite, quadhermitecub00 

contains
!***********************************************!
! subroutine  expandhermite
!  to complete afterwards if there need of order > 3
!***********************************************!
SUBROUTINE expandhermite(pol_pts, pol_value, pol_der, pol_herm)
  use types
  use lagrange1d_mod
  implicit none
  type(tPolynom),intent(in)             :: pol_pts, pol_value, pol_der
  type(tPolynom),intent(inout)          :: pol_herm

! lokal Variables
  type(tPolynom)            :: polsq_lagr, pol_1, pol_mul
  real                      :: pval_1, res_1, der_1
  integer                   :: i,j

! verification
  write(*,101) (pol_pts%coeffs(i), i=1,pol_pts%n)
  write(*,101) (pol_value%coeffs(i), i=1,pol_value%n)
  write(*,101) (pol_der%coeffs(i), i=1,pol_der%n)
! initialization
  pol_herm%n=2*pol_pts%n-1
  do i=1,pol_herm%n+1
     pol_herm%coeffs(i)=0.0;
  enddo

  do i=1,pol_pts%n
     ! initialization
     polsq_lagr%n = 2*pol_pts%n-2
     do j=1,polsq_lagr%n
        polsq_lagr%coeffs(j)=0.0;
     enddo

     call expandsquare(pol_pts, i, polsq_lagr)

     pol_mul%n=polsq_lagr%n
     do j=1,pol_mul%n
        pol_mul%coeffs(j) = polsq_lagr%coeffs(j)  ! copy values from square of ith Lagrange polynom
     enddo
     call evaluation(polsq_lagr, pol_pts%coeffs(i), res_1)
     call calcderivative(polsq_lagr, pol_pts%coeffs(i), der_1)
     pval_1=der_1*pol_value%coeffs(i) - res_1*pol_der%coeffs(i)
   write(*,101) (pval_1)
     pol_1%n=2
     pol_1%coeffs(1)=-pval_1
     pol_1%coeffs(2)=(-pval_1)*(-pol_pts%coeffs(i)) + pol_value%coeffs(i)
     call expand2nom(pol_1, pol_mul)   ! Result ist im pol_mul ueberschrieben

    ! Kopie der ith polynom im sum
!     polsq_herm%n = 2*pol_pts%n-1
     do j=1,polsq_lagr%n+2
        pol_herm%coeffs(j)= pol_herm%coeffs(j)+pol_mul%coeffs(j);
     enddo  

  enddo

!  Verifikation
  write(*,101) (pol_herm%coeffs(i), i=1,pol_herm%n+1)
  write(*,*)'------ end expand hermite -------'
101 format (f10.5)
END SUBROUTINE expandhermite

!***********************************************!
! subroutine  quadhermitecub00 even if unused at the moment
!  find the simpson-quadrature coeffs to calculate
!   the integral int(h00(x)*f(x)dx
!  where: h00(x)= 2*x^3 - 3*x^2 + 1
!***********************************************!
SUBROUTINE quadhermitecub00(pol_herm, quad, resint)
  use types
  use lagrange1d_mod
  implicit none
  type(tPolynom),intent(in)          :: pol_herm
  real, dimension(4) ,intent(out)    :: quad
  real,intent(out)                   :: resint


END SUBROUTINE quadhermitecub00



END MODULE hermite1d_mod
