MODULE Hermite2d_mod
!***************************************************************!
!
! Berechnet die Koeffizienten des Interpolationspolynoms
!  , das die derivative-kontinuive Interpolation macht und
!  berechnet die Operations gradient, div(gradient) und rot(gradient)
!  aus, damit die Resultate in das Finite-Element program schnell
!  implementiert werden.
!
!***************************************************************!
  use types
  implicit none
 
  private
 
  interface fillhermitepol
    module procedure fillhermitepol
  end interface

  interface calchermitepol
    module procedure calchermitepol
  end interface

  interface evalhermitepol
    module procedure evalhermitepol
  end interface

 !-----------------------------------------
  public:: fillhermitepol, calchermitepol, evalhermitepol
 !-----------------------------------------
  contains

!**************************************************************!
!							       !
!  Subroutine fillhermitepol                                   !
!                             !
!**************************************************************!
    SUBROUTINE fillhermitepol(n, quad, MatrixInterp)
    use types
    implicit none
    integer,intent(in)            :: n ! Order der Interpolation ==3
    logical,intent(in)            :: quad ! if true Rechteck (0,1)*(0,1) als Referenz Element
    real,dimension(1),intent(out) :: MatrixInterp(100) ! Matrix entsprechend zum Netzelement
    ! lokal 
    integer                       :: i,j
    real                          :: powx,powy,term,derx,dery
    real,dimension(2)              :: x_1=(/-1,1/)
    real,dimension(2)              :: y_1=(/-1,1/)
  !  if (n.eq.3)  then    ! only order 3 is supported
       ! Hier Identity Matrix zumm Start musst aktualisiert werden !
     do i=0,9
       do j=0,9
         MatrixInterp(1+i+j*10)=0.0
       enddo
       MatrixInterp(i+i*10+1)=1.0
     enddo
       
  !  first category, no algorithm here but a plenty of code!!!
  !  erste Equation:: p(x,y)=(a0*+a1*X+a2*Y+a3*XY+a4*X2+a5*Y2+a6*X3+a7*Y3+a8*X2Y+a9*XY2) eine VanderMonde Matrix muss gefuellt werden.

!--------
!
!Complicated
!
!--



    ! Pruefen das Resultat der Matrix am Schirm:
    do i=1,10
      write(*,201) (MatrixInterp(i+j*10), j=0,9)
    enddo
201  format (f10.5)
    END SUBROUTINE fillhermitepol


    SUBROUTINE calchermitepol
    use types
    implicit none
    ! Variables
    integer,intent(in)       :: n  !Order==3
    real,pointer,intent(in)          :: Coefficients(:,:)  ! Der Input: 10 Koefficients x NNodes
    ! lokalen Variablen
    real,dimension(100)      :: MatrixInterp
    real,dimension(10)       :: solution
    integer                  :: k,i
    integer ,dimension(10)   :: pivot
    integer                  :: ok     ! used by LAPACK FORTRAN subroutine SGESV

   k=1
      call fillhermitepol(n, .true. , MatrixInterp)  ! Fuellt die Matrixinterp aus
                                                   !  ,um die Koeffizient per linear Inversion
					           !  in "solution" zu berechnen
      solution(1)=1.0   ! First equation, pi(ai)=1
      solution(2)=-1.0  ! Second equation, dpi%dx(ai)=-1
      solution(3)=-1.0  ! Third equation, dpi%dy(ai)=-1
      solution(4)=0.0   ! Fourth equation, pi(ai+1)=0
      solution(5)=0.0  ! Fifth equation, dpi%dx(ai+1)=0
      solution(6)=0.0   ! Sixth equation, pi(ai+2)=0
      solution(7)=-0.707  ! Seventh equation, dpi%dx(ai+2)=-0.707
      solution(8)=-0.707  ! Eighth equation, dpi%dy(ai+2)=-0.707
      solution(9)=0.0   ! Ninth equation, pi(ai+3)=0
      solution(10)=0.0 ! Tenth equation, dpi%dy(ai+3)=0

   ! call Lapack procedure, rhs solution will be erased and filled by the result of linear
 !     call SGESV( 10, 1, MatrixInterp, 10, pivot, solution, 10, ok)
      do i=1,10
         Coefficients(i,k)=solution(i)
      enddo 

    ! Pruefen das Resultat der Coeffs am Schirm:
    write(*,*)
    write(*,303) '    Coefficients k   =', k 
    do i=1,10
      write(*,203) (Coefficients(i,1))
    enddo

203  format (f10.5)
303  format (a, i10)
    END SUBROUTINE calchermitepol
 

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

