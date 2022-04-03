module Input_mod
use typesround
implicit none

private

  interface input_mshfile
     module procedure input_mshfile
  end interface

  interface input_inpfile
     module procedure input_inpfile
  end interface

!------------------------------------------------------------!
 public :: input_mshfile, input_inpfile

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
!                                                                          !
!      - FUER PLATTEN SYSTEM LOESER                                        !
!            EY             Young Module                                   !
!            mu             Poisson Koeffizient                            !
!                                                                          !
!      EINGABE-PARAMETER - OPTION 2: ARRAYS OF NODES,ELEMENTS,BOUNDARY     !
!                                                                          !
!      - RECHENGEBIET:                                                     !
!            filename      Name of Gmsh file, which is organized as opt 2  !
!                                                                          !
!      - FUER PLATTEN SYSTEM LOESER                                        !
!            EY             Young Module                                   !
!            mu             Poisson Koeffizient                            !
!                                                                          !
!**************************************************************************!

!**** DELETED
!subroutine Input(Mesh,Exakt,Const,FileIO)
!end subroutine Input

 subroutine Input_mshfile(NTestMeshR,MeshR)  ! FileIO for .msh file

  use typesround
!  use allocateFieldsRound_mod

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tRoundMeshInfo)            :: NTestMeshR  ! ne,nn                     !
  type(tRoundMesh)                :: MeshR       ! Gitterwerte               !
  !--------------------------------------------------------------------------!
 !  character(LEN=13)          :: Mshfilename  ! Name of the input file (13 characters) in current dir !
  integer                    :: node, elmt, bound ! counter integers

  !write(*,106) ' Mesh filename: (e.g. Untitled1.msh) '
  !read(5,106)  MshFilename
  
  OPEN(UNIT=25, FILE='untitled2.msh', ACTION='READ')
  ! .. couts Nodes and Elements and Boudary Element...

  read (25,*)                ! $MeshFormat
  read (25,*)                ! 2.2 0 8
  read (25,*)                ! $EndMeshFormat
  read (25,*)                ! $PhysicalNames

  read (25, 304)  NTestMeshR%nBound              ! for ex. here 4 boundary
!  read (25,*)                ! 1 1 "inlet"
!  read (25,*)                ! 1 2 "top"
!  read (25,*)                ! 1 3 "bottom"
!  read (25,*)                ! 1 4 "outlet"
  do bound=1,NTestMeshR%nBound
    read(25,*)
  enddo
  read (25,*)                ! $EndPhysicalNames
  read (25,*)                ! $Nodes

  read (25, 304)  NTestMeshR%nNode  ! Anzahl Punkten (e.g. 56)
  do node=1,NTestMeshR%nNode
     read(25,*)              ! Skip it, for further, go on createMeshRound.f90
  enddo
  read (25,*)                ! $EndNodes

  read (25,*)                ! $Elements
  read (25, 304)  NTestMeshR%nElem ! Anzahl Elements (e.g. 56)
  do elmt=1,NTestMeshR%nElem
     read(25,*)              ! Skip it, for further, go on createMeshRound.f90
  enddo

  CLOSE(UNIT=25)  ! close msh file

304  format (i10)

  end subroutine Input_mshfile
!**************************************************************************!
!                                                                          !
! subroutine Input_inpfile                                                 !
!                                                                          !
!      Function: Count Elements and Nodes in the input text file           !
!                                                                          !
!**************************************************************************!

 subroutine Input_inpfile(NTestMeshR,MeshR)  ! FileIO for .inp file

  use typesround
!  use allocateFieldsRound_mod

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tRoundMeshInfo)            :: NTestMeshR  ! ne,nn                     !
  type(tRoundMesh)                :: MeshR       ! Gitterwerte               !
  !--------------------------------------------------------------------------!
  character(LEN=40)          :: Rawtext  
  ! character(LEN=13)          :: Inpfilename  ! Name of the input file (13 characters) in current dir !
  integer                    :: nodes, elmts, bound ! counter integers

  !write(*,106) ' Mesh filename: (e.g. untitled2.inp) '
  !read(5,106)  InpFilename
  
  OPEN(UNIT=25, FILE='patch12equi.inp', ACTION='READ')
  ! .. couts Nodes and Elements and Boudary Element...

  read(25,*)   ! *Heading
  read(25,*)   ! /home/fredrik/Documents/gmsh/untitled.inp
  read(25,*)   ! *NODE
  nodes=0
  read(25,106) Rawtext
  do while ((Rawtext).ne.('******* E L E M E N T S *************'))
     nodes=nodes+1
     read(25, 106) Rawtext
  enddo
  NTestMeshR%nNode=nodes  ! save the counted nodes so that the content of the file can be imported in createMeshRound module
  write(*,305) nodes
  read(25, 106) Rawtext
  do while ((Rawtext).ne.('*ELEMENT, type=CPS3, ELSET=Surface1'))
     read(25, 106) Rawtext
  enddo
  elmts=0
  read(25, 106) Rawtext
  do while ((Rawtext).ne.('*ELSET,ELSET=inlet'))
     elmts=elmts+1
     read(25, 106) Rawtext
  enddo
  write(*,305) elmts
  NTestMeshR%nElem=elmts ! save the counted elements
  CLOSE(UNIT=25) ! close inp input file

  ! TRICKHERE SHOULD BE UPDATED
  NTestMeshR%nxmx = max(5,nodes)  ! frontier length should be optimized in the future ..
  NTestMeshR%nymx = max(5,nodes)  ! but essential in the present. If not: seg. fault!
 !

105  format (a20)
106  format (a,a)
205  format (f10.5)
206  format (a,f10.5)
305  format (i10)
306  format (a,i10)

  end subroutine Input_inpfile  !

end module Input_mod
