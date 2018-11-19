module verifyquad2D_mod

  implicit none

!**************************************************************************!
!                                                                          !
! module verifyquad2D                                                      !
!                                                                          !
!      Funktion: Pruefen, ob jede 8-Punkten basierte Quadrangle Element    !
!                seine Derivatien fuer Koordinaten x,y als invertierbar    !
!                hat.                                                      !
!                                                                          !
!**************************************************************************!

  public :: calcpolycoord, calcpolycoordderx, calcpolycoorddery, checkderivcoord, checkalledgesofelmts

!**************************************************************************!
!  Variable Tabelle der Polynomen fuer ein 8-Punkten basierte Quadrangle   !
!*          Element. Ordnung der Koeff ist: 1, x, y, x**2, x*y, y**2,      !
!                                            x**2*y, x*y**2                !
!**************************************************************************!

  ! MUSS DIVIDIERT VON 4.0 SEIN
  real,dimension(64)        :: coeffpolycoord = (/ (/ -1.0/0.4, 0.0, 0.0, 1.0/4.0, 1.0/4.0, 1.0/4.0, -1.0/4.0, -1.0/4.0 /), &
                                                   (/ 2.0/4.0, 0.0, -2.0/4.0, -2.0/4.0, 0.0, 0.0, -2.0/4.0, 0.0 /), &
                                                   (/ -1.0/4.0, 0.0, 0.0, 1.0/4.0, -1.0/4.0, 1.0/4.0, -1.0/4.0, 1.0/4.0 /), &
                                                   (/ 2.0/4.0, 2.0/4.0, 0.0, 0.0, 0.0, -2.0/4.0, 0.0, -2.0/4.0 /), &
                                                   (/ -1.0/4.0, 0.0, 0.0, 1.0/4.0, -1.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/4.0 /), &
                                                   (/ 2.0/4.0, 0.0, 2.0/4.0, -2.0/4.0, 0.0, 0.0, -2.0/4.0, 0.0 /), &
                                                   (/ -1.0/4.0, 0.0, 0.0, 1.0/4.0, -1.0/4.0, 1.0/4.0, 1.0/4.0, -1.0/4.0 /), &
                                                   (/ 2.0/4.0, -2.0/4.0, 0.0, 0.0, 0.0, -2.0/4.0, 0.0, 2.0/4.0 /) /)

 ! und derivaten Polynoms je nach Variablen xoder y  (1,8)   KORREKTIERT ABER NICHT SICHER      !
  real,dimension(64)        :: coeffpolycoordderx = (/ (/ 0.0, 2.0/4.0, 1.0/4.0, 0.0, -2.0/4.0, -1.0/4.0, 0.0, 0.0 /), &
                                                   (/ 0.0, -4.0/4.0, 0.0, 0.0, 4.0/4.0, 0.0, 0.0, 0.0 /), &
                                                   (/ 0.0, 2.0/4.0, -1.0/4.0, 0.0, -2.0/4.0, 1.0/4.0, 0.0, 0.0 /), &
                                                   (/ 2.0/4.0, 0.0, 0.0, 0.0, 0.0, -2.0/4.0, 0.0, 0.0 /), &
                                                   (/ 0.0, 2.0/4.0, 1.0/4.0, 2.0/4.0, 1.0/4.0, 0.0, 0.0, 0.0 /), &
                                                   (/ 0.0, -4.0/4.0, 0.0, 0.0, -4.0/4.0, 0.0, 0.0, 0.0 /), &
                                                   (/ 0.0, 2.0/4.0, -1.0/4.0, 0.0, 2.0/4.0, -1.0/4.0, 0.0, 0.0 /), &
                                                   (/ -2.0/4.0, 0.0, 0.0, 0.0, 0.0, 2.0/4.0, 0.0, 0.0 /) /)
  ! y SIEHT KORREKT AUS
  real,dimension(64)       :: coeffpolycoorddery = (/ (/ 0.0, 1.0/4.0, 2.0/4.0, -1.0/4.0, -2.0/4.0, 0.0, 0.0, 0.0 /), &
                                                   (/ -2.0/4.0, 0.0, 0.0, 2.0/4.0, 0.0, 0.0, 0.0, 0.0 /), &
                                                   (/ 0.0, -1.0/4.0, 2.0/4.0, -1.0/4.0, 2.0/4.0, 0.0, 0.0, 0.0 /), &
                                                   (/ 0.0, 0.0, -4.0/4.0, 0.0, -4.0/4.0, 0.0, 0.0, 0.0 /), &
                                                   (/ 0.0, 1.0/4.0, 2.0/4.0, 1.0/4.0, 2.0/4.0, 0.0, 0.0, 0.0 /), &
                                                   (/ 2.0/4.0, 0.0, 0.0, -2.0/4.0, 0.0, 0.0, 0.0, 0.0 /), &
                                                   (/ 0.0, -1.0/4.0, 2.0/4.0, 1.0/4.0, -2.0/4.0, 0.0, 0.0, 0.0 /), &
                                                   (/ 0.0, 0.0, -4.0/4.0, 0.0, 4.0/4.0, 0.0, 0.0, 0.0 /) /)


  contains

  subroutine calcpolycoord(i,x,y,val)

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  integer,intent(in)          :: i     !  vom 1 bis 8
  real,intent(in)             :: x,y
  real,intent(out)            :: val
  !--------------------------------------------------------------------------!
  !                                                                          !
  ! Local variable declaration                                               !  
  integer                     :: k          ! Zählvariablen                  !
  real                        :: summe                                       !
  !--------------------------------------------------------------------------!
  summe = 0.0
  do k=0,7
    select case(k) ! wahl des Polynoms
      case (0)
       summe = summe + coeffpolycoord(8*(i-1)+k+1)*1
      case (1)
       summe = summe + coeffpolycoord(8*(i-1)+k+1)*x
      case (2)
       summe = summe + coeffpolycoord(8*(i-1)+k+1)*y
      case (3)
       summe = summe + coeffpolycoord(8*(i-1)+k+1)*x*x 
      case (4)
        summe = summe + coeffpolycoord(8*(i-1)+k+1)*x*y
      case (5)
       summe = summe + coeffpolycoord(8*(i-1)+k+1)*y*y
      case (6)
       summe = summe + coeffpolycoord(8*(i-1)+k+1)*x*x*y
      case (7)
       summe = summe + coeffpolycoord(8*(i-1)+k+1)*x*y*y 
    end select
  enddo
  val = summe
  end subroutine calcpolycoord


  subroutine calcpolycoordderx(i,x,y,val)

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  integer,intent(in)          :: i     !  vom 1 bis 8
  real,intent(in)             :: x,y
  real,intent(out)            :: val
  !--------------------------------------------------------------------------!
  !                                                                          !
  ! Local variable declaration                                               !  
  integer                     :: k          ! Zählvariablen                  !
  real                        :: summe                                       !
  !--------------------------------------------------------------------------!
  summe = 0.0
  do k=0,7
    select case(k) ! wahl des Polynoms
      case (0)
       summe = summe + coeffpolycoordderx(8*(i-1)+k+1)*1
      case (1)
       summe = summe + coeffpolycoordderx(8*(i-1)+k+1)*x
      case (2)
       summe = summe + coeffpolycoordderx(8*(i-1)+k+1)*y
      case (3)
       summe = summe + coeffpolycoordderx(8*(i-1)+k+1)*x*x 
      case (4)
        summe = summe + coeffpolycoordderx(8*(i-1)+k+1)*x*y
      case (5)
       summe = summe + coeffpolycoordderx(8*(i-1)+k+1)*y*y
      case (6)
       summe = summe + coeffpolycoordderx(8*(i-1)+k+1)*x*x*y
      case (7)
       summe = summe + coeffpolycoordderx(8*(i-1)+k+1)*x*y*y 
    end select
  enddo
  val = summe
  end subroutine calcpolycoordderx


  subroutine calcpolycoorddery(i,x,y,val)

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  integer,intent(in)          :: i     !  vom 1 bis 8
  real,intent(in)             :: x,y
  real,intent(out)            :: val
  !--------------------------------------------------------------------------!
  !                                                                          !
  ! Local variable declaration                                               !  
  integer                     :: k          ! Zählvariablen                  !
  real                        :: summe                                       !
  !--------------------------------------------------------------------------!
  summe = 0.0
  do k=0,7
    select case(k) ! wahl des Polynoms
      case (0)
       summe = summe + coeffpolycoorddery(8*(i-1)+k+1)*1
      case (1)
       summe = summe + coeffpolycoorddery(8*(i-1)+k+1)*x
      case (2)
       summe = summe + coeffpolycoorddery(8*(i-1)+k+1)*y
      case (3)
       summe = summe + coeffpolycoorddery(8*(i-1)+k+1)*x*x 
      case (4)
        summe = summe + coeffpolycoorddery(8*(i-1)+k+1)*x*y
      case (5)
       summe = summe + coeffpolycoorddery(8*(i-1)+k+1)*y*y
      case (6)
       summe = summe + coeffpolycoorddery(8*(i-1)+k+1)*x*x*y
      case (7)
       summe = summe + coeffpolycoorddery(8*(i-1)+k+1)*x*y*y 
    end select
  enddo
  val = summe
  end subroutine calcpolycoorddery

  !**************************************************************************!
  !  subroutine checkderivcoord:                                             !
  !                                                                          !
  !  Berechnet alle Gradient-Matrizen an 4 prinzipallen Punkten des          !
  !  8 Punkt-Elements und stellt eine globale Konditionsnummer fuer das      !
  !  Element dar, condit, der -1 ist, wenn ein null Gradient berechnet ist   !
  !                                                                          !
  !**************************************************************************!
  subroutine checkderivcoord(Mesh2D,nel,condit)

  use types

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tMeshGen)             :: Mesh2D       ! Gitterwerte                   !
  integer                    :: nel          ! Indiz des 8-Elements          !
  real                       :: condit       ! stat. Resultat fuer das Elemt !
  !                                                                          !
  ! Local variable declaration                                               !
  !                                                                          !
  real                       :: x,y, mn,mx
  real,dimension(4)          :: nodedet ! Resultat, Genau Konditionng.       !
  real,dimension(4)          :: Jnod_det !, Jnod_norm2,                      !
  real,dimension(16)         :: Jnod    ! 4 2x2 Matrizen Jakobian in Node i  !
  real,dimension(8)          :: xc, f_derx, yc, f_dery                       !
  integer                    :: i,j,k      ! Zählvariablen                   !
  integer,dimension(4)       :: i_array                                      !
  !--------------------------------------------------------------------------!
  intent(inout)              :: Mesh2D, condit 
  intent(in)                 :: nel 
  !--------------------------------------------------------------------------!

!  write(*,304) '      elmt   = ' , nel
  do i=1,8
    xc(i) = Mesh2D%x( Mesh2D%quad(i,nel) )
    yc(i) = Mesh2D%y( Mesh2D%quad(i,nel) )
!    write(*,204) '      xc{ }   = ' , xc(i)
!    write(*,204) '      yc{ }   = ' , yc(i)
  enddo

  i_array = (/ 1,3,5,7 /)
  do k=1,4
    i = i_array(k)  ! array subscripting
    select case(i) 
      case (1)
        x=-1.0
        y=-1.0
      case (3)
        x=1.0
        y=-1.0
      case (5)
        x=1.0
        y=1.0
      case (7)
        x=-1.0
        y=1.0
!      case (2)
!        x=-0.0
!        y=-1.0
!      case (4)
!        x=1.0
!        y=0.0
!      case (6)
!        x=0.0
!        y=1.0
!      case (8)
!        x=-1.0
!        y=0.0
    end select
    
    do j=1,8
      call calcpolycoordderx( j, x, y, f_derx(j) )    ! X-Derivative von Nj an dem Punkt i
 !     write(*,204) '      f_derx{ }   = ' , f_derx(j)
      call calcpolycoorddery( j, x, y, f_dery(j) )    ! Y-Derivative von Nj an dem Punkt i
 !     write(*,204) '      f_dery{ }   = ' , f_dery(j)
    enddo

    ! Rechnen der Derivative der Koordinaten-Interpolation mit Skalar Produkt
    ! Hier auch eine Transformation vom (1,3,5,7) zum (0,1,2,3)
    Jnod(4*(k-1)+1)=0.0
    do j=1,8
      Jnod(4*(k-1)+1) = Jnod(4*(k-1)+1) + f_derx(j)*xc(j)
    enddo

    Jnod(4*(k-1)+2)=0.0
    do j=1,8
      Jnod(4*(k-1)+2) = Jnod(4*(k-1)+2) + f_dery(j)*xc(j)
    enddo

    Jnod(4*(k-1)+3)=0.0
    do j=1,8
      Jnod(4*(k-1)+3) = Jnod(4*(k-1)+3) + f_derx(j)*yc(j)
    enddo

    Jnod(4*(k-1)+4)=0.0
    do j=1,8
      Jnod(4*(k-1)+4) = Jnod(4*(k-1)+4) + f_dery(j)*yc(i)
    enddo
    Jnod_det(k) = (Jnod(4*(k-1)+1)*Jnod(4*(k-1)+4)-Jnod(4*(k-1)+2)*Jnod(4*(k-1)+3))
    !Jnod_norm2(k) = (Jnod(4*(k-1)+1)*Jnod(4*(k-1)+1)+Jnod(4*(k-1)+2)*Jnod(4*(k-1)+2))
    !Jnod_norm2(k) = Jnod_norm2(k)*(Jnod(4*(k-1)+3)*Jnod(4*(k-1)+3)+Jnod(4*(k-1)+4)*Jnod(4*(k-1)+4))

    nodedet(k) = abs(Jnod_det(k))  !   !!

!    write(*,304) '      node   = ' , i
!    write(*,204) '      J_1{node}   = ' , Jnod(4*(k-1)+1)
!    write(*,204) '      J_2{node}   = ' , Jnod(4*(k-1)+2)
!    write(*,204) '      J_3{node}   = ' , Jnod(4*(k-1)+3)
!    write(*,204) '      J_4{node}   = ' , Jnod(4*(k-1)+4)

  enddo  ! index k


  ! Suche nach dem Min/Max von nodedet
  mn=nodedet(1)
  mx=nodedet(1)
  if (mn > nodedet(2)) then
     mn = nodedet(2)
  elseif (mx < nodedet(2)) then
     mx = nodedet(2)
  endif
  if (mn > nodedet(3)) then
     mn = nodedet(3)
  elseif (mx < nodedet(3)) then
     mx = nodedet(3)
  endif
  if (mn > nodedet(4)) then
     mn = nodedet(4)
  elseif (mx < nodedet(4)) then
     mx = nodedet(4)
  endif
  ! End Result: Konditionnierungsnummer fuer das Element i (nur 4 von Punkten gewertet):
  if (abs(mx) > 0) then
    condit = abs(mn/mx)
  else
    condit = -1.0
  endif
 104  format (f10.5)
 204  format (a,f10.5)
 304  format (a,i10)
  end subroutine checkderivcoord

  !**************************************************************************!
  !  subroutine checkalledgesofelmts:                                        !
  !                                                                          !
  !  Geht durch den generish Gitter von 8-Elementen und ruft die subroutine  !
  !  checkderivcoord individuell an. Stoppt, wenn nur eine Jakobian Matrix   !
  !  an einem Punkt des Elemts singulaer ist.                                !
  !                                                                          !
  !**************************************************************************!
  subroutine checkalledgesofelmts(Mesh2D)

  use types

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tMeshGen)             :: Mesh2D       ! Gitterwerte                   !
  !                                                                          !
  ! Local variable declaration                                               !
  !                                                                          !
  logical                    :: singular                                     !
  real                       :: resu       ! Konditionsnummer                !
  integer                    :: i          ! Zählvariablen                   !
  !--------------------------------------------------------------------------!
  intent(inout)              :: Mesh2D
  !--------------------------------------------------------------------------!
  i=1
  singular = .false.
  do while( (i<(Mesh2D%elmts+1)).AND.(.NOT.(singular)) )
    call checkderivcoord(Mesh2D, i, resu)
    if (resu <0) then
      write(*,302) '         BREAK, singular Jacobian of Geometry of Elmt i = ', i
      singular = .true.
      i=i+1  ! Sicherheit
    else
      write(*,302) '      elmt   = ' , i
      write(*,202) '      Cond Number{elmt}   = ' , resu
      i=i+1
    endif
  enddo

 102  format (f10.5)
 202  format (a,f10.5)
 302  format (a,i10)
 end subroutine checkalledgesofelmts

end module verifyquad2D_mod
