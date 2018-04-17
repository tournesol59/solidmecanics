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

  pol_herm%n=2*pol_pts%n
  
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
