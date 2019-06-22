!**************************************************************************!
!                                                                          !
! subroutine createMesh                                                    !
!                                                                          !
!      Funktion: Festlegen der Gitterdaten                                 !
!                                                                          !
!**************************************************************************!
!                                                                          !
!      x,y        Koordinaten der Gitterpunkte                             !
!      dx,  dy    Gitterschrittweite (delta x, delta y)                    !
!      dxq, dyq   Kehrwert (1/delta x, 1/delta y)                          !
!      dxx, dyy   Quadrat  (delta x * delta x, delta y * delta y)          !
!      dxxq,dyyq  Kehrwerte der Quadrate (1/dxx, 1/dyy)                    !
!                                                                          !
!**************************************************************************!
module createmesh_mod
  use TypesRound

  implicit none

  private
!********************************************************************!
!  Variable Tabelle der Test-Plate fur die Koordinaten der 9 Punkten !
!********************************************************************!
 X2 = (/ 0.0, 0.05, 0.0, 0.05, 0.1, 0.1, 0.0, 0.05, 0.1 /)

 Y2 = (/ 0.0, 0.0, 0.05, 0.05, 0.0, 0.05, 0.1, 0.1, 0.1 /)

 Z2 = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

  public :: createRoundMesh2, createRoundMeshFile, &
            createRoundRand2, createRoundRandFile

  contains

!******************************************************
  subroutine createRoundMesh2(NTestMeshR,MeshR)

  use TypesRound

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tRoundMeshInfo)            :: NTestMeshR  ! ne,nn                     !
  type(tRoundMesh)                :: MeshR       ! Gitterwerte               !
  !                                                                          !
  ! Local variable declaration                                               !
  !                                                                          !
  integer                    :: ne,nn,i,j        ! Zählvariablen             !
  !--------------------------------------------------------------------------!
  intent(inout)              :: NTestMeshR, MeshR                            !
  !--------------------------------------------------------------------------!


  !-----------------------------------<  calculate constants  >-----------
  !---<  nn   = anzahl der inneren punkte                 >-----------
  !---<  ne = anzahl der zellen                           >-----------

  NTestMeshR%nn = 9
  NTestMeshR%ne = 8
  nn = NTestMeshR%nn
  ne = NTestMeshR%ne

  !---<  node(1..3..7)             = linker  rand             >-----------
  !---<  node(5..6..9)             = rechter rand             >-----------
  !---<  node(1..2..5)             = unten rand               >-----------
  !---<  node(7..8..9)             = obenn rand               >-----------
  do i= 1,nn
    MeshR%nodes(i)%vect_x = X2(i)
    MeshR%nodes(i)%vect_y = Y2(i)
    MeshR%nodes(i)%vect_z = Z2(i)
  enddo 

  
  end subroutine createRoundMesh2

!******************************************************
  subroutine createRoundRand2(NTestMeshR,BCU,BCS)

  use TypesRound

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tRoundMeshInfo)            :: NTestMeshR  ! ne,nn                     !
  type(tRoundRandU)               :: BCU       ! Randbedingungen Biegg       !
  type(tRoundRandS)               :: BCS       ! Randbedingungen Kraft       !
  !                                                                          !
  ! Local variable declaration                                               !
  !                                                                          !
  integer                    :: ne,nn,i,j        ! Zählvariablen             !
  !--------------------------------------------------------------------------!
  intent(in)                 ::  NTestMeshR
  intent(inout)              ::  BCU,BCS                                     !
  !--------------------------------------------------------------------------!


  end subroutine createRoundRand2


!******************************************************
  subroutine createRoundMeshFile(NTestMeshR,MeshR)

  use TypesRound

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tRoundMeshInfo)            :: NTestMeshR  ! ne,nn                     !
  type(tRoundMesh)                :: MeshR       ! Gitterwerte               !
  !                                                                          !
  ! Local variable declaration                                               !
  !                                                                          !
  integer                    :: ne,nn,i,j        ! Zählvariablen             !
  !--------------------------------------------------------------------------!
  intent(inout)              :: NTestMeshR, MeshR                            !
  !--------------------------------------------------------------------------!


  !-----------------------------------<  calculate constants  >-----------
  !---<  nn   = anzahl der inneren punkte                 >-----------
  !---<  ne = anzahl der zellen                           >-----------

!  NTestMeshR%nn = 9
!  NTestMeshR%ne = 8
!  nn = NTestMeshR%nn
!  ne = NTestMeshR%ne

 105  format (a20)
 106  format (a,a)
 205  format (f10.5)
 206  format (a,f10.5)
 305  format (i10)
 306  format (a,i10)

  end subroutine createRoundMeshFile

!******************************************************
  subroutine createRoundRandFile(NTestMeshR,BCU,BCS)

  use TypesRound

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tRoundMeshInfo)            :: NTestMeshR  ! ne,nn                     !
  type(tRoundRandU)               :: BCU       ! Randbedingungen Biegg       !
  type(tRoundRandS)               :: BCS       ! Randbedingungen Kraft       !
  !                                                                          !
  ! Local variable declaration                                               !
  !                                                                          !
  integer                    :: ne,nn,i,j        ! Zählvariablen             !
  !--------------------------------------------------------------------------!
  intent(in)                 ::  NTestMeshR
  intent(inout)              ::  BCU,BCS                                     !
  !--------------------------------------------------------------------------!


  end subroutine createRoundRandFile




end module createmesh_mod
