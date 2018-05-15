module Input_mod
use types
implicit none

private
  interface input
     module procedure input
  end interface


!------------------------------------------------------------!
 public :: input

contains
!**************************************************************************!
!                                                                          !
! subroutine input                                                         !
!                                                                          !
!      Funktion: Einlesen der Eingabedaten                                 !
!                                                                          !
!**************************************************************************!
!                                                                          !
!                                                                          !
!                                                                          !
!      EINGABE-PARAMETER:                                                  !
!                                                                          !
!      - RECHENGEBIET:                                                     !
!            startx,starty  Position der linken unteren Ecke               !
!            endx,endy      Position der rechten oberen Ecke               !
!            nraumx,nraumy  Anzahl Gitterpunkte in x- bzw. y-Richtung      !
!                                                                          !							   !                                                !
!      - FUER PLATTEN SYSTEM LOESER                                        !
!            EY             Young Module                                   !
!            mu             Poisson Koeffizient                            !
!                                                                          !
!**************************************************************************!

subroutine Input(Mesh,Exakt,Const,FileIO)

  use types
  use allocateFields_mod

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tMesh)                :: Mesh       ! Gitterwerte                     !
  type(tExakt)               :: Exakt      ! Exakte Loesung                  !
  type(tConstants)           :: Const      ! Konstanten                      !
  type(tFileIO)              :: FileIO     ! Ausgabesteuerung                !
  !--------------------------------------------------------------------------!
 
 ! OPEN(UNIT=25, FILE='solidinput.txt')

  write(*,*) 'Einlesen der Werte:'

  ! --- Rechengebiet:
  read (5,*)
  read (5,*)
  read (5,201) Mesh%startx
  read (5,201) Mesh%starty
  read (5,201) Mesh%endx
  read (5,201) Mesh%endy
  read (5,301) Mesh%nraumx
  read (5,301) Mesh%nraumy
  

  ! --- Parameters:
  read (5,*)
  read (5,*)       

  read (5,201) Const%EY
  read (5,201) Const%nu
!  Const%h = 0.15  ! Could not explain why read fails for h

  ! ----------------------------------<  Kontroll-Ausgabe der Groessen >--
  
  write(*,102) ' Rechengebiet: '
  write(*,202) '    startx   = ' ,Mesh%startx
  write(*,202) '    starty   = ' ,Mesh%starty
  write(*,202) '    endx     = ' ,Mesh%endx
  write(*,202) '    endy     = ' ,Mesh%endy

  write(*,302) '    nraumx   = ' ,Mesh%nraumx
  write(*,302) '    nraumy   = ' ,Mesh%nraumy
  write(*,102) ' '

  write(*,102) 'platten material parameter'  
!  write(*,202) '    h        = ' ,Const%h
!  write(*,102) ' '
  write(*,202) '    EY       = ' ,Const%EY
  write(*,102) ' '
  write(*,202) '    nu       = ' ,Const%nu
  write(*,102) ' '

  write(*,*) 'confirm (1) or press (0) to choose a predefined curved square meshregion'

  ! --- Rechengebiet:
  read (5,301) Const%auto

  ! ----------------------------------<  set the logical values  >--------


  Const%Pi     = 4.*atan(1.0)

  ! ----------------------------------<  check input data  >--------------


  ! ----------------------------------<  Ein/Ausgabe-Formate >------------
     
 101  format (a20)
 102  format (a,a)
 201  format (f10.5)
 202  format (a,f10.5)
 301  format (i10)
 302  format (a,i10)

end subroutine Input

end module Input_mod
