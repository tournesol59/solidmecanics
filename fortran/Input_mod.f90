module Input_mod
use types
implicit none

private
  interface input
     module procedure input
  end interface

  interface input_file
     module procedure input_file
  end interface

!------------------------------------------------------------!
 public :: input, input_file

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
!      EINGABE-PARAMETER - OPTION 1: A MATRIZ OF POINTS                    !
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
!      EINGABE-PARAMETER - OPTION 2: ARRAYS OF NODES,ELEMENTS,BOUNDARY     !
!                                                                          !
!      - RECHENGEBIET:                                                     !
!            filename      Name of Gmsh file, which is organized as opt 2  !
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



subroutine Input_file(Mesh2D,Exakt,Const,FileIO)

  use types
  use allocateFields_mod

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tMeshGen)             :: Mesh2D     ! Gitterwerte                     !
  type(tExakt)               :: Exakt      ! Exakte Loesung                  !
  type(tConstants)           :: Const      ! Konstanten                      !
  type(tFileIO)              :: FileIO     ! Ausgabesteuerung                !
  !--------------------------------------------------------------------------!
  character(LEN=13)          :: Mshfilename  ! Name of the input file (13 characters) in current dir !
  character(LEN=6)           :: bndname
  integer                    :: i, nbound, cbound, bnddim, bndind ! counter integers
  integer                    :: nodes, nodeind
  real                       :: nodex, nodey, nodez
  integer                    :: elmts, elind, eldim, eltype, elsize, eldim2, elnod1, elnod2, elnod3, elnod4
  !write(*,106) ' Mesh filename: (e.g. Untitled1.msh) '
  !read(5,106)  MshFilename
  
  OPEN(UNIT=25, FILE='Untitled1.msh', ACTION='READ')
  ! .. couts Nodes and Elements and Boundary Element...

  read (25,*)                ! $MeshFormat
  read (25,*)                ! 2.2 0 8
  read (25,*)                ! $EndMeshFormat
  read (25,*)                ! $PhysicalNames
  read (25,305), nbound       ! 
  cbound=0
  do i=1,nbound
     read(25,405), bnddim, bndind, bndname  ! e.g.  1 1 "Inlet1"
     if (bnddim.EQ.1) then
        if (bndname.EQ."FREEF") then             ! Eine Kraft wird gesetzt, Biegungfrei (w,s=0)
           Mesh2D%fronttype(1, cbound+1)=2
        elseif (bndname.EQ."FREED") then        ! Eine Verschiebung wird gesetzt, Biegungfrei (w,s=0)
           Mesh2D%fronttype(1, cbound+1)=4   
        elseif (bndname.EQ."SETAM") then        ! Ein Moment wird gesetzt, ohne Kraft
           Mesh2D%fronttype(1, cbound+1)=1
        elseif (bndname.EQ."SETAB") then        ! Ein Angle wird gesetzt, ohne Verschiebung
           Mesh2D%fronttype(1, cbound+1)=4
        elseif (bndname.EQ."SETFM") then        ! Eine Kraft und Moment werden gesetzt
           Mesh2D%fronttype(1, cbound+1)=5
        elseif (bndname.EQ."SETDB") then        ! Eine Verschiebung und Angle werden gesetzt
           Mesh2D%fronttype(1, cbound+1)=6
        endif
        cbound=cbound+1
     endif        
  ! elseif (bnddim.EQ.2) ! wird nicht hier verstanden
  enddo                      !
  read (25,*)                ! $EndPhysicalNames

  read (25,*)                ! $Nodes
  read (25, 305)  Mesh2D%nodes  ! Anzahl Punkten (e.g. 56)
  do nodes=1,Mesh2D%nodes
     read(25,406), nodeind, nodex, nodey, nodez              ! e.g. 1 0 0 0
     Mesh2D%x(nodeind)=nodex
     Mesh2D%y(nodeind)=nodey
     Mesh2D%z(nodeind)=nodez
  enddo
  read (25,*)                ! $EndNodes

  read (25,*)                ! $Elements
  read (25, 305)  Mesh2D%elmts ! Anzahl Elements (e.g. 56)
  do elmts=1,Mesh2D%elmts
                             ! wird noch  nicht komplett verstanden
    read(25,508), elind, eldim, eltype, elsize, eldim2, elnod1, elnod2, elnod3
               ! wenn <tag> 2: Drei Werten am Ende der Linie sind wichtig
    Mesh2D%allelmt(1, elind)= elnod1
    Mesh2D%allelmt(2, elind)= elnod2
    Mesh2D%allelmt(3, elind)= elnod3

  enddo
  read (25,*)                ! $EndElements
! SET ANZAHL RAND ELEMENTEN MANUALLY : KEINE GEEIGNETE ALGORITHMUS GEFUNDEN !


  CLOSE(UNIT=25)

 105  format (a20)
 106  format (a,a)
 205  format (f18.16)
 206  format (a,f18.16)
 305  format (i10)
 306  format (a,i10)
 405  format (i10, i10, a6)
 406  format (i10, f18.16, f18.16, f18.16)
 507  format (i10, i10, i10, i10, i10, i10, i10)
 508  format (i10, i10, i10, i10, i10, i10, i10, i10)
 509  format (i10, i10, i10, i10, i10, i10, i10, i10, i10)
end subroutine Input_file

end module Input_mod
