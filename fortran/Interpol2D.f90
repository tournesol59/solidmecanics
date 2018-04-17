!****************************************************************************!
!                                                                            !
! module Interpol2D                                                          !
!                                                                            !
!      Funktion: Interpol Funktionen definieren, welcher Domain eines        ! 
!                Netzelement entspricht                                      !
!****************************************************************************!

module Interpol2D_mod

  implicit none
  private
  
!  interface fillinterpol
!     module procedure fillinterpol
!  end interface

!  interface interpolref
!     module procedure interpolref
!  end interface
   
!  interface basisdifferenz
!     module procedure basisdifferenz
!  end interface
 
!  interface christoffei
!     module procedure christoffei
!  end interface

!  interface lgeval
!     module procedure lgeval
!  end interface
 
!  interface trychristoffei
!     module procedure trychristoffei
!  end interface

!------------------------------------------------------------!
 public :: fillinterpol, interpolref, basisdifferenz, christoffei &
            , lgeval, trychristoffei  !, lgderivS, lgderivf_S
!------------------------------------------------------------!
 ! Global Variablen
 !!! type(tNode),pointer  :: derivbasis(:, :)  ! 2D array von real basis Vektor-derivativen  in der naturalen x,y,z Basis Richtungen!

!!!  real,pointer        :: chicoeff(:,:,:) ! 2D array of 3-array von Christoffei componenten der Vektor-derivativen in der (u,t,n) Basis Richtungen

 contains
  !**************************************************************************!
  !                                                                          !
  ! Subroutine fillinterpol                                                  !
  !                                                                          !
  !      procedure: Ausfuellt die Matrix (Array) als 2D Van der Monde Matrix !
  !       , um die Interpolation Polynomen zwischen Real/Refnetz u bestimmen ! 
  !                                                                          !
  !**************************************************************************!

  SUBROUTINE fillinterpol(n, k, x, y, MatrixInterp, Vektor)
    use types
    implicit none
    integer,intent(in)            :: n ! Order der Interpolation
    integer,intent(in)            :: k ! Indiz des NetzElements
    real,pointer                  :: x(:)
    real,pointer                  :: y(:)
    real,dimension(1:36),intent(out) :: MatrixInterp ! Matrix entsprechend zum Netzelement
    real,dimension(1:6)           :: Vektor
    ! lokal 
    integer                       :: i,j
    real                          :: x_1,y_1
   ! MatrixInterp
    do j=1,4
    select case(j)
         case (1)
                x_1=0.0
                y_1=0.0
         case (2)
                x_1=1.0
                y_1=0.0
         case (3)
                x_1=1.0
                y_1=1.0
         case (4)
                x_1=0.0
                y_1=1.0
      end select
      do i=1,6
      select case(i)
            case (1)
                  MatrixInterp(j*6+i)=1
            case (2)
                  MatrixInterp(j*6+i)=x_1
            case (3)
                  MatrixInterp(j*6+i)=y_1
            case (4)
                  MatrixInterp(j*6+i)=x_1*y_1
            case (5)
                  MatrixInterp(j*6+i)=0.0  ! Kein Termen in x^2
            case (6)
                  MatrixInterp(j*6+i)=0.0  ! Kein Termen in y^2
         end select
      end do
    end do
    do j=1,6
      do i=5,6
          MatrixInterp(j*6+i)=0.0! 2 Kompletende Gleichungen:
      end do
    end do
    MatrixInterp(4*6+5)=1.0    ! Gleichung 1: e5=gewuenschte Value, here zero
    MatrixInterp(5*6+6)=1.0    ! Gleichung 2: f6=0.0

  ! RHS Vektor:
    do j=1,4
     Vektor(1)=x(k-1)
     Vektor(2)=x(k)
     Vektor(3)=y(k-1)
     Vektor(4)=y(k)
   enddo
   Vektor(5)=0.0   ! Gleichung 1: e5=gewuenschte Value, here zero
   Vektor(6)=0.0   ! Gleichung 2: f6=0.0

  END SUBROUTINE fillinterpol

  !**************************************************************************!
  !                                                                          !
  ! Subroutine interpolref                                                   !
  !                                                                          !
  !      procedure: call fillinterpol fuer each Element und loest das linear ! 
  !        System aus und kopiert each element in dem Coefficient            ! 
  !                                                                          !
  !**************************************************************************!

  SUBROUTINE interpolref(n, Gitter, Coefficients)
    use types
    implicit none
    integer,intent(in)            :: n ! Order der Interpolation
    type(tMesh),intent(in)        :: Gitter
    real,dimension(1:6,1:100),intent(out)    :: Coefficients
    !lokal
    real,dimension(1:36)          :: MatrixInterp
    real,dimension(1:6)           :: solution
    integer                       :: i,j,k, pivot(6), ok
    
    do k=1,100
      ! solution is the right hand side . Das wird mit der Loesung ueberschrieben
      call fillinterpol(n, k, Gitter%x, Gitter%y, MatrixInterp, solution)
    ! call LAPACK procedure
    !  call SGESV( 6, 1, MatrixInterp, 6, pivot, solution, 6, ok)

    !  do i=1,6
    !     Coefficients(i,k)=solution(i)  ! Kopie der Loesung in dem richtigen Table
    !  enddo
    ! und geht weiter
    enddo

  END SUBROUTINE interpolref

  !**************************************************************************!
  !                                                                          !
  ! Subroutine basisdifferenz                                                !
  !                                                                          !
  !      procedure: Differenziert den lokal Basis fuer jeden Node des realen !
  !     Netzes und dann macht eine Projektion (naechste Prozedur)            !
  !                                                                          !
  !**************************************************************************!

  SUBROUTINE basisdifferenz(Gitter,der_a1xx_x,der_a1xx_y,der_a1xy_x,der_a1xy_y &
                                 ,der_a2yy_x, der_a2yy_y, der_a2yx_x, der_a2yx_y)
   use types 
   implicit none
    type(tMesh),intent(in)                   :: Gitter
    real,pointer,intent(inout)               :: der_a1xx_x(:)  !=der(ax)%dx, _x Koord
    real,pointer,intent(inout)               :: der_a1xx_y(:)  !           , _y Koord

    real,pointer,intent(inout)               :: der_a1xy_x(:)  !=der(ax)%dx, _x Koord
    real,pointer,intent(inout)               :: der_a1xy_y(:)  !           , _y Koord

    real,pointer,intent(inout)               :: der_a2yy_x(:)  !=der(ay)%dy, _x Koord
    real,pointer,intent(inout)               :: der_a2yy_y(:)  !           , _y Koord

    real,pointer,intent(inout)               :: der_a2yx_x(:)  !=der(ay)%dx, _x Koord
    real,pointer,intent(inout)               :: der_a2yx_y(:)  !=der(ay)%dx, _y Koord
    ! lokal
    integer                                  :: i, j

    ! Algorithm fuer das erste Vektor in x-Richtung, von x differenziert
    do j=0,Gitter%nraumy-2
      do i=0,Gitter%nraumx-2
        der_a1xx_x(j*Gitter%nraumx+i+2)=Gitter%dxxq*(Gitter%x(i+3+j*Gitter%nraumx)  &
                                         - 2*Gitter%x(i+2+j*Gitter%nraumx)  &
                                          + Gitter%x(i+1+j*Gitter%nraumx) ) 
        der_a1xx_y(j*Gitter%nraumx+i+2)=Gitter%dxxq*(Gitter%y(i+3+j*Gitter%nraumx)  &
                                         - 2*Gitter%y(i+2+j*Gitter%nraumx)  &
                                          + Gitter%y(i+1+j*Gitter%nraumx) ) 
      end do
      der_a1xx_x(j*Gitter%nraumx+1)=der_a1xx_x(j*Gitter%nraumx+2)
      der_a1xx_y(j*Gitter%nraumx+1)=der_a1xx_y(j*Gitter%nraumx+2)
    end do
     do j=0,Gitter%nraumy-2  ! Kopie des vorheriges index fuer das Rand, weil differenzieren unmoeglich ist !
      der_a1xx_x(j*Gitter%nraumx+Gitter%nraumy-1)=der_a1xx_x(j*Gitter%nraumx+Gitter%nraumy-2)
      der_a1xx_x(j*Gitter%nraumx+Gitter%nraumy)=der_a1xx_x(j*Gitter%nraumx+Gitter%nraumy-2)
      der_a1xx_y(j*Gitter%nraumx+Gitter%nraumy-1)=der_a1xx_y(j*Gitter%nraumx+Gitter%nraumy-2)
      der_a1xx_y(j*Gitter%nraumx+Gitter%nraumy)=der_a1xx_y(j*Gitter%nraumx+Gitter%nraumy-2)
     end do

    ! Algorithm fuer das erste Vektor in x-Richtung, von y differenziert
   do j=0,Gitter%nraumy-2
     do i=0,Gitter%nraumx-2
        der_a1xy_x((j+1)*Gitter%nraumx+i+1) =Gitter%dyq*(Gitter%dxq*(Gitter%x(i+2+(j+2)*Gitter%nraumx)  &
                                                                - Gitter%x(i+1+(j+2)*Gitter%nraumx))  &
                                                     - Gitter%dxq*(Gitter%x(i+2+(j+1)*Gitter%nraumx)  &
                                                                - Gitter%x(i+1+(j+1)*Gitter%nraumx)) ) 
        der_a1xy_y((j+1)*Gitter%nraumx+i+1) =Gitter%dyq*(Gitter%dxq*(Gitter%y(i+2+(j+2)*Gitter%nraumx)  &
                                                                - Gitter%y(i+1+(j+2)*Gitter%nraumx))  &
                                                     - Gitter%dxq*(Gitter%y(i+2+(j+1)*Gitter%nraumx)  &
                                                                - Gitter%y(i+1+(j+1)*Gitter%nraumx)) ) 
    end do
    der_a1xy_x((j+1)*Gitter%nraumx+Gitter%nraumx) =der_a1xy_x((j+1)*Gitter%nraumx+Gitter%nraumx-1)
    der_a1xy_y((j+1)*Gitter%nraumx+Gitter%nraumx) =der_a1xy_y((j+1)*Gitter%nraumx+Gitter%nraumx-1)
  end do    
  do i=0,Gitter%nraumx-2
     der_a1xy_x((Gitter%nraumy-1)*Gitter%nraumx+i+1) =der_a1xy_x((Gitter%nraumy-2)*Gitter%nraumx+i+1)
     der_a1xy_y((Gitter%nraumy-1)*Gitter%nraumx+i+1) =der_a1xy_y((Gitter%nraumy-2)*Gitter%nraumx+i+1)
  end do

    ! Algorithm fuer das zweite Vektor in y-Richtung, von x differenziert
   do j=0,Gitter%nraumy-2
     do i=0,Gitter%nraumx-2
        der_a2yx_x((j+1)*Gitter%nraumx+i+1) =Gitter%dxq*(Gitter%dyq*(Gitter%x(i+2+(j+2)*Gitter%nraumx)  &
                                           - Gitter%x(i+2+(j+1)*Gitter%nraumx))  &
                                      - Gitter%dyq*(Gitter%x(i+1+(j+2)*Gitter%nraumx)  &
                                            - Gitter%x(i+1+(j+1)*Gitter%nraumx)))
        der_a2yx_y((j+1)*Gitter%nraumx+i+1) =Gitter%dxq*(Gitter%dyq*(Gitter%y(i+2+(j+2)*Gitter%nraumx)  &
                                           - Gitter%y(i+2+(j+1)*Gitter%nraumx))  &
                                      - Gitter%dyq*(Gitter%y(i+1+(j+2)*Gitter%nraumx)  &
                                            - Gitter%y(i+1+(j+1)*Gitter%nraumx)))
    end do
    der_a2yx_x((j+1)*Gitter%nraumx+Gitter%nraumx) =der_a2yx_x((j+1)*Gitter%nraumx+Gitter%nraumx-1)
    der_a2yx_y((j+1)*Gitter%nraumx+Gitter%nraumx) =der_a2yx_y((j+1)*Gitter%nraumx+Gitter%nraumx-1)
  end do
  do i=0,Gitter%nraumx-2
     der_a2yx_x((Gitter%nraumy-1)*Gitter%nraumx+i+1) =der_a1xy_x((Gitter%nraumy-2)*Gitter%nraumx+i+1)
     der_a2yx_y((Gitter%nraumy-1)*Gitter%nraumx+i+1) =der_a1xy_y((Gitter%nraumy-2)*Gitter%nraumx+i+1)
  end do

    ! Algorithm fuer das zweite Vektor in y-Richtung, von y differenziert
    do j=0,Gitter%nraumy-3
      do i=0,Gitter%nraumx-2
        der_a2yy_y((j+1)*Gitter%nraumx+i+1)= Gitter%dyyq*(Gitter%y(i+1+(j+2)*Gitter%nraumx)  & 
                                          - 2*Gitter%y(i+1+(j+1)*Gitter%nraumx)  &
                                            + Gitter%y(i+1+j*Gitter%nraumx) ) 
        der_a2yy_x((j+1)*Gitter%nraumx+i+1)= Gitter%dyyq*(Gitter%x(i+1+(j+2)*Gitter%nraumx)  & 
                                          - 2*Gitter%x(i+1+(j+1)*Gitter%nraumx)  &
                                            + Gitter%x(i+1+j*Gitter%nraumx) ) 
     end do
      der_a2yy_x(j*Gitter%nraumx+1)=der_a2yy_x(j*Gitter%nraumx+2)
      der_a2yy_y(j*Gitter%nraumx+1)=der_a2yy_y(j*Gitter%nraumx+2)
   end do
   do i=0,Gitter%nraumx-2
       der_a2yy_x((Gitter%nraumy-2)*Gitter%nraumx+i+1)=der_a2yy_x((Gitter%nraumy-3)*Gitter%nraumx+i+1)
       der_a2yy_x((Gitter%nraumy-1)*Gitter%nraumx+i+1)=der_a2yy_x((Gitter%nraumy-3)*Gitter%nraumx+i+1)
       der_a2yy_y((Gitter%nraumy-2)*Gitter%nraumx+i+1)=der_a2yy_y((Gitter%nraumy-3)*Gitter%nraumx+i+1)
       der_a2yy_y((Gitter%nraumy-1)*Gitter%nraumx+i+1)=der_a2yy_y((Gitter%nraumy-3)*Gitter%nraumx+i+1)
   end do

  END SUBROUTINE basisdifferenz

  !**************************************************************************!
  !                                                                          !
  ! Subroutine christoffei                                                   !
  !                                                                          !
  !    nach basisdifferenz macht eine Projektion ueber das lokalen Basis zum !
  !     ermitteln des Christoffei Koeffizienten aber sie werden alle Null,   !
  !     da das realNetz planar ist.                                          !
  !                                                                          !
  !**************************************************************************!
 
  SUBROUTINE christoffei(Gitter, chicoeff, der_a1xx_x,der_a1xx_y,der_a1xy_x,der_a1xy_y &
                                 ,der_a2yy_x, der_a2yy_y, der_a2yx_x, der_a2yx_y)
   use types 
   implicit none
    type(tMesh),intent(in)                   :: Gitter
    real,intent(inout)               :: chicoeff(1:8,1:100)  ! 2nd Dimension ist nraumx*nraumy
    ! Ein Kristoffei coefficient K{g}{a,b} wird oft in dem Formeln benutzt
    ! der(Ua)%dxb = sum,g K{g}{a,b}*Ug,  wo die Ua,Ug sind das lokal BasisVektoren
    ! Hier die der(Ua)%dxb Koordinaten, aber planar!
    real,pointer,intent(inout)               :: der_a1xx_x(:)  !=der(ax)%dx, _x Koord
    real,pointer,intent(inout)               :: der_a1xx_y(:)  !           , _y Koord

    real,pointer,intent(inout)               :: der_a1xy_x(:)  !=der(ax)%dx, _x Koord
    real,pointer,intent(inout)               :: der_a1xy_y(:)  !           , _y Koord

    real,pointer,intent(inout)               :: der_a2yy_x(:)  !=der(ay)%dy, _x Koord
    real,pointer,intent(inout)               :: der_a2yy_y(:)  !           , _y Koord

    real,pointer,intent(inout)               :: der_a2yx_x(:)  !=der(ay)%dx, _x Koord
    real,pointer,intent(inout)               :: der_a2yx_y(:)  !=der(ay)%dx, _y Koord
    ! Lokal Variables
    real                                  :: chi, a1x,a1y,a2x,a2y, norm
    integer                               :: i,j,k

    do j=0,Gitter%nraumy-1
      do i=0,Gitter%nraumx-1
        do k=0,7
           chicoeff(k+1,j*Gitter%nraumx+i+1)=0.0  ! Init
        enddo
      enddo
   enddo
          ! Berechnen die lokal Basisvecktorskomponenten:
   do j=0,Gitter%nraumy-2
     do i=0,Gitter%nraumx-2
        k=j*Gitter%nraumx+i+1

        a1x=Gitter%x(i+2+j*Gitter%nraumx)-Gitter%x(i+1+j*Gitter%nraumx)
        a1y=Gitter%y(i+2+j*Gitter%nraumx)-Gitter%y(i+1+j*Gitter%nraumx)

        a2x=Gitter%x(i+1+(j+1)*Gitter%nraumx)-Gitter%x(i+1+j*Gitter%nraumx)
        a2y=Gitter%y(i+1+(j+1)*Gitter%nraumx)-Gitter%y(i+1+j*Gitter%nraumx)

        ! Projektion der der_aixx,xy,yx,yy ueber diese Basisvektors
        
        norm=sqrt(a1x*a1x+a1y*a1y)
        chi=der_a1xx_x(k)*a1x+der_a1xx_y(k)*a1y   ! Gamma{x}{xx}
        chicoeff(1,k)=chi/norm

        norm=sqrt(a2x*a2x+a2y*a2y)
        chi=der_a1xx_x(k)*a2x+der_a1xx_y(k)*a2y   ! Gamma{y}{xx}
        chicoeff(2,k)=chi/norm

        norm=sqrt(a1x*a1x+a1y*a1y)
        chi=der_a2yx_x(k)*a1x+der_a2yx_y(k)*a1y   ! Gamma{x}{xy}
        chicoeff(3,k)=chi/norm

        norm=sqrt(a2x*a2x+a2y*a2y)
        chi=der_a2yx_x(k)*a2x+der_a2yx_y(k)*a2y   ! Gamma{y}{xy}
        chicoeff(4,k)=chi/norm                    	

        norm=sqrt(a1x*a1x+a1y*a1y)
        chi=der_a1xy_x(k)*a1x+der_a1xy_y(k)*a1y   ! Gamma{x}{yx}
        chicoeff(5,k)=chi/norm

        norm=sqrt(a2x*a2x+a2y*a2y)
        chi=der_a1xy_x(k)*a2x+der_a1xy_y(k)*a2y   ! Gamma{y}{yx}
        chicoeff(6,k)=chi/norm

        norm=sqrt(a1x*a1x+a1y*a1y)
        chi=der_a2yy_x(k)*a1x+der_a2yy_y(k)*a1y   ! Gamma{x}{yy}
        chicoeff(7,k)=chi/norm

        norm=sqrt(a2x*a2x+a2y*a2y)
        chi=der_a2yy_x(k)*a2x+der_a2yy_y(k)*a2y   ! Gamma{y}{yy}
        chicoeff(8,k)=chi/norm

     enddo
   enddo
  END SUBROUTINE christoffei


  !**************************************************************************!
  !                                                                          !
  ! Funktion lgeval                                                          !
  !                                                                          !
  !      Funktion: berechnet die Werten von Polynoms an dem Node (xi,yi)     !
  !     generische Lagrange 2D Polynoms fuer Interpolation der realen Loesung! 
  !                                                                          !
  !**************************************************************************!
  subroutine lgeval(n,x,y,pi,fxy)
     use types
     implicit none
     integer,intent(in)           :: n   ! Order der Interpolation2D 2 = x,y,x*y,x**2, y**2
     real,intent(in)              :: x
     real,intent(in)              :: y ! Punkt der Evaluation   
     integer,intent(in)           :: pi  ! Indiz des Facets
     real,intent(out)             :: fxy  ! Resultat
     ! pi nimmt 6 differenten Forms auf einen Referenznetz von 6 facets
     ! facet a/ eqn: x  +z=1 
     ! facet b/ eqn:   y+z=1 
     ! facet c/ eqn:-x+y+z=1 
     ! facet d/ eqn:-x  +z=1 
     ! facet e/ eqn:  -y+z=1 
     ! facet f/ eqn: x-y+z=1 
     select case (pi)
       case (1)
          fxy = 1-x
       case (2)
          fxy = 1-y
       case (3)
          fxy = 1+x-y
       case (4)
          fxy = 1+x
       case (5)
          fxy = 1+y
       case (6)
          fxy = 1-x+y
     end select
        
!     WRITE(*,*)  'Value p(',pi,') =', fxy
  END subroutine lgeval


  !**************************************************************************!
  !                                                                          !
  ! Funktion lgeval3                                                         !
  !                                                                          !
  !      Funktion: berechnet die Werten von Polynoms an dem Node (xi,yi)     !
  !      generische Lagrange 2D Polynoms mit zweiter Ordnung (x*y)           ! 
  !      fuer Interpolation der realen Loesung                               !
  !**************************************************************************!

  subroutine lgeval3(n,x,y,pi,fxy)
     use types
     implicit none
     integer,intent(in)           :: n   ! Order der Interpolation2D 2 = x,y,x*y,x**2, y**2
     real,intent(in)              :: x
     real,intent(in)              :: y ! Punkt der Evaluation   
     integer,intent(in)           :: pi  ! Indiz des Polynoms
     real,intent(out)             :: fxy  ! Resultat
     ! pi nimmt 4 differenten Forms auf einen Referenznetz
     ! a/ p1(x,y)=0.25(x*y-y-x+1)
     ! b/ p2(x,y)=0.25(-x*y-y+x+1)
     ! c/ p3(x,y)=0.25(x*y+y+x+1)
     ! d/ p4(x,y)=0.25(-x*y+y-x+1)

     select case (pi)
       case (1)
          fxy = 0.25*(x*y-y-x+1)
       case (2)
          fxy = 0.25*(-x*y-y+x+1)
       case (3)
          fxy = 0.25*(x*y+y+x+1)
       case (4)
          fxy = 0.25*(-x*y+y-x+1)
     end select
        
!     WRITE(*,*)  'Value p(',pi,') =', fxy
  END subroutine lgeval3


  !**************************************************************************!
  !                                                                          !
  ! Funktion trychristoffei                                                  !
  !                                                                          !
  !      Funktion: nach basisdifferenz und christoffei call, prueft die Werte!
  !                von chicoeff, indem man sie am Bildschirm displayen       ! 
  !                                                                          !
  !**************************************************************************!

  SUBROUTINE trychristoffei(Gitter, chicoeff)
     use types
     implicit none
    type(tMesh),intent(in)           :: Gitter
    real,intent(inout)               :: chicoeff(1:8,1:100)  ! 2nd Dimension ist nraumx*nraumy
    ! lokal Variables
    integer                          :: i,j,k,l

  ! Display input
    write(*,*)
    do j=0,Gitter%nraumy-1
         write(*,302) '     x-data column j   =',j
         write(*,102) (Gitter%x(j*Gitter%nraumx+i+1), i=0,Gitter%nraumx-1)
    enddo
    write(*,*)
    do j=0,Gitter%nraumy-1
         write(*,302) '     y-data column j   =',j
         write(*,102) (Gitter%y(j*Gitter%nraumx+i+1), i=0,Gitter%nraumx-1)
    enddo
  ! Display output
    write(*,*)
    do j=0,Gitter%nraumy-1
      do i=0,Gitter%nraumx-1
         k=j*Gitter%nraumx+i+1
         write(*,302) '    Node k   =' ,k,'  {'

          write(*,202) '      gamma{x}{xx}   = ' ,chicoeff(1,k)
          write(*,202) '      gamma{y}{xx}   = ' ,chicoeff(2,k)
          write(*,202) '      gamma{x}{xy}   = ' ,chicoeff(3,k)
          write(*,202) '      gamma{y}{xy}   = ' ,chicoeff(4,k)
          write(*,202) '      gamma{x}{yx}   = ' ,chicoeff(5,k)
          write(*,202) '      gamma{y}{yx}   = ' ,chicoeff(6,k)
          write(*,202) '      gamma{x}{yy}   = ' ,chicoeff(7,k)
          write(*,202) '      gamma{y}{yy}   = ' ,chicoeff(8,k)

        write(*,*) '   }'
        write(*,*)

      enddo
    enddo

 102  format (f10.5)
 202  format (a,f10.5)
 302  format (a,i10)
  END SUBROUTINE trychristoffei

end module Interpol2D_mod

