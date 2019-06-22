MODULE Hermite2d_mod
!***************************************************************!
!
! Berechnet die Koeffizienten des Interpolationspolynoms
!  , das die derivative-kontinuive Interpolation macht und
!  berechnet die Operations gradient, div(gradient) und moeglws. rot(gradient)
!  aus, damit die Resultate in das Finite-Element program schnell
!  implementiert werden.
!
!***************************************************************!
  use types
  implicit none
 
  private
 !********************************************************************!
!  Variable Tabelle der innere produkt von Komposition-Gradienten    !
!        je nach Indiz des Netzelements und der Position der         !
!        subelement  in der Platte                                   !
!                                                                    !
!********************************************************************!
! Tabelle der Gauss Points for approximation of an integral over a square region
! used in linlg2dfillmatrixq
  real,dimension(8)   :: GaussPoints = (/ (/ 0.57735, 0.57735 /), &
                             (/ 0.57735, -0.57735 /), &
                             (/ -0.57735, 0.57735 /), &
                             (/ -0.57735, -0.57735 /) /)

  interface evalhermitepol
    module procedure evalhermitepol
  end interface

 !-----------------------------------------
  public:: linhm2allocatesdns, computeGaussInt, linhm2dfillmatrixq, evalhermitepol
 !-----------------------------------------
  contains
!********************************************************************!
!  procedure linhm2allocatesdns DENSE MATRICES: n=10--> 10000 memory !
!********************************************************************!
 SUBROUTINE linhm2allocatesdns(Gitter, SkalarProdMatrixQ, SkalarProdMatrixT)
 use types
 use interpol2D_mod
 use Lagrange1d_mod
 implicit none
  type(tMesh)         :: Gitter
  integer             :: allocStat  
  real,pointer        :: SkalarProdMatrixQ(:,:)
  real,pointer        :: SkalarProdMatrixT(:,:,:)
  intent(in)          :: Gitter
  intent(inout)       :: SkalarProdMatrixQ, SkalarProdMatrixT
  
 allocate(SkalarProdMatrixQ(1:Gitter%nraumx*Gitter%nraumy,1:Gitter%nraumx*Gitter%nraumy), & 
     SkalarProdMatrixT(1:Gitter%nraumx*Gitter%nraumy,1:Gitter%nraumx*Gitter%nraumy,1:4), & 
     STAT = allocStat) 
      !    MatrixQ : Gleichung, 1-Dim
      !    MatrixT : Gleichung, 4-Dims
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not allocate all variables!'
       STOP
    end if
 END SUBROUTINE linhm2allocatesdns

!********************************************************************!
!  procedure computeGaussInt   : Perform the filling of the SkalarProdMatrix   !
!              ! tt=0: Prop*Unity, 1:Prop*DerX, 2:Prop*DerX2, 3: DerX*DerX, 4:DerX*DerX2
!                   5: DerX2*DerX2, 6: DerX*DerX3
!              ! ii: Art of the polynoms: position of one elemt to another
!********************************************************************!
 SUBROUTINE computeGaussInt(tt, ii, pol1, pol2, pol3, valend)
 use types
 use interpol2D_mod
 use Lagrange1d_mod
 implicit none
 integer,intent(in)	     :: tt, ii

 real,intent(out)                    :: valend
 type(tPolynom),intent(in)           :: pol1, pol2, pol3
 real                                :: val1, val2, val3, val4
 real                                :: der1, der2, der3, der4
 integer                             :: j
 select case(tt) ! Aim of the calculus: should we have Prop*Unity...

 case(0)  !Prop*Prop
  do j=0,3 
   select case(ii)  ! Art of polynoms in the element
   case(0)   
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     !call evaluation(pol1, GaussPoints(2*j+1), val3)
     !call evaluation(pol1, GaussPoints(2*j+2), val4)
     valend=valend+val1*val2*val1*val2*1.0;
   case(-1)
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     call evaluation(pol2, GaussPoints(2*j+1), val3)
     call evaluation(pol2, GaussPoints(2*j+2), val4)
     valend=valend+val1*val2*val3*val4*1.0;
  case(+1)
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     call evaluation(pol3, GaussPoints(2*j+1), val3)
     call evaluation(pol3, GaussPoints(2*j+2), val4)
     valend=valend+val1*val2*val3*val4*1.0;
  end select
  enddo
 end select
! return valend
 END SUBROUTINE computeGaussInt

!********************************************************************!
!  procedure linhm2dfillmatrixq:                                     !
!    depends on computeGaussInt(tt, ii, pol1, pol2, pol3, valuex)    !
!        , which is to be actualized to permit derivation up to fourth degree !
!******************************************************************!
 SUBROUTINE linhm2dfillmatrixq(Gitter, chicoeff, VarNum, SkalarProdMatrixQ)
 use types
 use interpol2D_mod
 use Lagrange1d_mod
 implicit none
 type(tMesh),intent(in)          :: Gitter
  real,intent(inout)             :: chicoeff(1:8,1:100)  ! 2nd Dimension ist
  type(tNumeric),intent(inout)   :: VarNum
  real,pointer,intent(inout)     :: SkalarProdMatrixQ(:,:)
  ! lokale Variablen
  integer                        :: i,j,k,r, l,m,n,p
  real                           :: K_x, K_y, costh, sinth, cosps, sinps, valuex, valuey
  type(tPolynom)                 :: pol1, pol2, pol3  !dies sind die polynoms von hermite ueber segment [0 1] !
   pol1%n=3
   pol1%coeffs= (/ 2.0, -3.0, 0.0 , 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
   pol2%n=3
   pol2%coeffs= (/ 1.0, -2.0, 1.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
   pol3%n=3
   pol3%coeffs= (/ -2.0, 3.0, 0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
   pol3%n=3
   pol3%coeffs= (/ 1.0, -1.0, 0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

   ! NICHT ZU ENDE VERFASSEN

 END SUBROUTINE linhm2dfillmatrixq


   SUBROUTINE evalhermitepol(n, l, Coefficients, Values)
    use types
    implicit none
    ! Variables
    integer,intent(in)       :: n  !Order==3
    integer,intent(in)       :: l  !(node)
    real,intent(in)          :: Coefficients(1:10,50)  ! Der Input: 10 Koefficients x NNodes
    real,intent(inout)       :: Values(:)  ! Der Resultat fuer ein Node
    ! lokalen Variablen
    integer                  :: i,j,k
    real                     :: x(9)
    real                     :: y(9)
   do i=0,8
     x(i)=0.0+i*0.125
     y(j)=0.0+j*0.125
   enddo
   do k=1,1 ! Nur eine Teste Node, alle sind gleich
     do i=1,9
       do j=1,9
          Values((j-1)*9+i)=Coefficients(1,k)+ &
                    Coefficients(2,k)*x(i)+ &
                    Coefficients(3,k)*y(j)+ &
                    Coefficients(4,k)*x(i)*y(j)+ &
                    Coefficients(5,k)*y(j)*y(j)+ &
                    Coefficients(6,k)*x(i)*x(i)+ &
                    Coefficients(7,k)*x(j)*x(i)*x(i)+ &
                    Coefficients(8,k)*x(j)*x(i)*y(j)+ &
                    Coefficients(9,k)*x(j)*y(j)*y(j)+ &
                    Coefficients(10,k)*y(j)*y(j)*y(j)
       enddo
     enddo
   enddo
    ! Pruefen das Resultat der erste Interpolation am Schirm:
    do i=1,9
      write(*,201) (Values((j-1)*9+i), j=0,9)
      write(*,*)
    enddo
201  format (f10.5)
  END SUBROUTINE evalhermitepol



END MODULE Hermite2d_mod

