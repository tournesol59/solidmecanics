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

!***** calcderivative : see how it can be used from Lagrange1D *************!

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
  !--------------------------------------------
  real                               :: t, sum1, sum2, sum3, sum4
  real                               :: res1, res2, res3, res4
  real                               :: val1, val2, val3, val4
  real, dimension(4)                 :: valp, valh
  integer                            :: i,j
  type(tPolynom)                     :: pol_h00
  type(tPolynom)                     :: pol_lagr1,pol_lagr2,pol_lagr3,pol_lagr4
  type(tPolynom)                     :: roots_T4, pol_val1, pol_val2, pol_val3, pol_val4

  ! initiate h00(x)
  pol_h00%n=3
  pol_h00%coeffs(1)=2.0
  pol_h00%coeffs(2)=-3.0
  pol_h00%coeffs(3)=0.0
  pol_h00%coeffs(4)=1.0

  ! compute chebyshev roots tik, poly is used only for structure array coeff
  roots_T4%n=3
  roots_T4%coeffs(1) = sin((-3)*(3.1415)/(2*4))
  roots_T4%coeffs(2) = sin((-1)*(3.1415)/(2*4))
  roots_T4%coeffs(3) = sin((1)*(3.1415)/(2*4))
  roots_T4%coeffs(4) = sin((3)*(3.1415)/(2*4))
  pol_val1%n=3
  pol_val1%coeffs = (/ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  pol_val2%n=3
  pol_val2%coeffs = (/ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  pol_val3%n=3
  pol_val3%coeffs = (/ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  pol_val4%n=3
  pol_val4%coeffs = (/ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

  call expandlagrange(roots_T4, pol_val1, pol_lagr1)  
  sum1=0.0
  do i=1,11
    t = -1.0 + 0.2*(i-1)
    call evaluation(pol_lagr1, t, res1)  ! simple  Riemann Integration
    sum1 = sum1 + 0.1*(res1*sqrt(1-t*t));
  enddo

  call expandlagrange(roots_T4, pol_val2, pol_lagr2) 
  sum2=0.0 
  do i=1,11
    t = -1.0 + 0.2*(i-1)
    call evaluation(pol_lagr2, t, res2)  ! simple  Riemann Integration
    sum2 = sum2 + 0.1*(res2*sqrt(1-t*t));
  enddo

  call expandlagrange(roots_T4, pol_val3, pol_lagr3) 
  sum3=0.0 
  do i=1,11
    t = -1.0 + 0.2*(i-1)
    call evaluation(pol_lagr3, t, res3)  ! simple  Riemann Integration
    sum3 = sum3 + 0.1*(res3*sqrt(1-t*t));
  enddo
  
  call expandlagrange(roots_T4, pol_val4, pol_lagr4) 
  sum4=0.0 
  do i=1,11
    t = -1.0 + 0.2*(i-1)
    call evaluation(pol_lagr4, t, res4)  ! simple  Riemann Integration
    sum4 = sum4 + 0.1*(res4*sqrt(1-t*t));
  enddo

  do i=1,4   !* the computed integral are, if exact integration, christoffei-coeffs, used for quadrature intgratio,
     quad(1)=sum1
     quad(2)=sum2
     quad(3)=sum3
     quad(4)=sum4
  enddo
 
  resint=0.0
  do i=1,4
    call evaluation(pol_herm, roots_T4%coeffs(i), valp(i))
    call evaluation(pol_h00, roots_T4%coeffs(i), valh(i))
    resint =resint + quad(i)*valp(i)*valh(i)
  enddo

!  Verifikation
  write(*,*)'------ calc chebyshev T4 roots          -------'
  write(*,101) (roots_T4%coeffs(i), i=1,4)
  write(*,*)'------ quadhemitecub00 : expand hermite -------'
  write(*,101) (pol_herm%coeffs(i), i=1,pol_herm%n+1)
  write(*,*)'------ calc quad cristoffei coeff       -------'
  write(*,101) (quad(i), i=1,4)
  write(*,*)'------ calc integ resultat              -------'
  write(*,101) (resint)

101 format (f10.5)
END SUBROUTINE quadhermitecub00

END MODULE hermite1d_mod
