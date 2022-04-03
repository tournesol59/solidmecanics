!!*****************************************************************************!
!* 
!*      $RCSfile: TypesDef.f90,v $
!*
!*      $Revision: 1.3 $
!*
!*      $Date: 2004/02/13 08:16:09 $
!*
!*      $Author: iagklima $
!*
!*****************************************************************************!

MODULE Types
  !---------------------------------------------------------------------------!
  !  use ISO_C_BINDING
  implicit none
  private
  !---------------------------------------------------------------------------!

  type :: tMeshInfo 
   ! Nur die Anzahl Gittepunkte sind wichtig.
   !Es kann z.B. mit variable Laengen der Elemente hinzugefuegt
     integer      :: nraumx         ! Anzahl Gitterpunkte x-Richtung  !
     integer      :: nraumy         ! Anzahl Gitterpunkte y-Richtung  !
     real       :: startx, starty ! Position der linken unteren Ecke!
     real       :: endx, endy     ! Position der rechten oberen Ecke!
  end type tMeshInfo

  type :: tMesh ! use only in 4-rectangular plate meshed as a matrix nx*ny pts
     integer      :: nraumx         ! Copy of tMeshInfo after standard input  !
     integer      :: nraumy         ! idem  !
     real       :: startx, starty ! idem  !
     real       :: endx, endy     ! idem  !
     real       :: xlen, ylen     ! Laengen des Rechengebiets in x/y !
     real       :: dx             ! Gitterschrittweite (Delta x)    !
     real       :: dy             ! Gitterschrittweite (Delta y)    !
     real       :: dxq            ! Kehrwert ( 1/ Delta x )         !
     real       :: dyq            ! Kehrwert ( 1/ Delta y )         !
     real       :: dxx            ! dx^2 (Delta x * Delta x)        !
     real       :: dyy            ! dy^2 (Delta y * Delta y )       !
     real       :: dxxq           ! Kehrwert des Quadrats (1/dxx)   !
     real       :: dyyq           ! Kehrwert des Quadrats (1/dyy)   !
     ! real,pointer (c_ptr) :: x(:)        ! do NOT work
!     type(c_ptr)         :: x,y,z
     
    real,pointer :: x(:),y(:),z(:)

  end type tMesh

  type :: tCoeff
     real,pointer        :: Coefficients(:,:,:)  ! (ai,i,j) fuer jede Quadrangle Element, das am Punkt (i,j) beginnt
  end type tCoeff

  type tMeshGen
     integer              :: nodes          ! Anzahl der Gitterpunkte (Nodes) !
     real,pointer         :: x(:),y(:),z(:) ! Koordinaten der Gitterpunkte    !
     integer              :: elmts          ! Anzahl der Quadranglen          !
     integer,pointer      :: quad(:,:)      ! Indiz der Punte fuer Quadrangle !
     integer,pointer      :: neighbours(:,:)
     integer              :: ntop           ! Anzahl der Quadrangle an oberen Wand !
     integer,pointer      :: quadtop(:)     ! Indiz der Quadrangle an oberen Wand  !
     integer              :: nbottom        ! Anzahl der Quadrangle an untern Wand !
     integer,pointer      :: quadbottom(:)  ! Indiz der Quadrangle an untern Wand  !
     integer              :: nleft          ! Anzahl der Quadrangle an linken Wand !
     integer,pointer      :: quadleft(:)    ! Indiz der Quadrangle an linken Wand  !
     integer              :: nright         ! Anzahl der Quadrangle an rechten Wand!
     integer,pointer      :: quadright(:)   ! Indiz der Quadrangle an rechten Wand !
  end type tMeshGen 

  type tCurv              ! Abhaengig type tMesh fuer (nraumx, nraumy) !
     real,pointer         :: chicoeff(:,:)          ! 8 Values fuer jedes Punktes: Christoffel Koeffizienten
     real,pointer         :: der_a1xx_x(:),der_a1xx_y(:),der_a1xy_x(:),der_a1xy_y(:), &
                          der_a2yx_x(:),der_a2yx_y(:),der_a2yy_x(:),der_a2yy_y(:)  
  end type tCurv

  type tNumeric
    real,pointer     :: Verschiebung(:,:) ! 3 Verschiebungen in der Mitte der Platte !
    real,pointer     :: TensionTg(:,:)    ! array pointer: Spannung tengential Tensor!
    real,pointer     :: Tension3(:,:)     ! array pointer: Spannung in Richtung 3  !
    real,pointer     :: BiegMoment(:,:)   ! BiegMoment Tensor, nicht tangential    !

!     real,pointer         :: Coefficients(:,:) ! des Polynoms fuer vordraengige Interpolation 

  end type tNumeric

  type tZeiten  ! wird nicht beutzt fuer Platten Festigkeit  !
     real                 :: zeit           ! Endzeitpunkt                    !
     real                 :: sigma          ! Sicherheit gegen Stabilitaet    !
     real                 :: dtin           ! delta t     (fuer delta t fest) !
     real                 :: time           ! Aktueller Zeitpunkt             !
     real                 :: dt             ! Aktuelle Zeitschrittweite       !
     real                 :: outputTime     ! Naechster Ausgabezeitpunkt      !
     integer              :: nZeit          ! Max. Anzahl Zeitschritte        !
     integer              :: nTime          ! Nr des aktuellen Zeitschritts   !
     character(20)        :: ccfl           ! Zeitschritt variabel oder fixed?!
     logical              :: druck          ! Ausgabezeitpunkt erreicht?      !
     logical              :: ende           ! Endzeitpunkt erreicht?          !
  end type tZeiten

  type :: tRandbedingungen
! 'trandbedingungen' is not C interoperable --> may be add type, bind(C)

     integer          :: randltype      ! RandTyp links: 1=Fixed, 2=Force
     real,pointer     :: randl(:)       ! Randwerte links                 !
     integer          :: randrtype      ! RandTyp rechts: 1=Fixed, 2=Force
      real,pointer    :: randr(:)      ! Randwerte rechts                !
     integer          :: randotype      ! RandTyp oben: 1=Fixed, 2=Force
     real,pointer     :: rando(:)      ! Randwerte oben                  !
     integer          :: randutype      ! RandTyp unten: 1=Fixed, 2=Force
    real,pointer      :: randu(:)      ! Randwerte unten                 !
  end type tRandbedingungen
  
  type tConstants  ! Original for Waermetransportsproblem geplant !
     real                 :: Pi             ! Pi                              !
     real                 :: EY             ! Young Module                    !
     real                 :: nu             ! Poisson coefficient             !
     real                 :: h             ! Groesse (Height), Konstant      !
     integer              :: automesh       ! Automatishes Auswahl eines vordefinierten Mesh !
         ! REST WIRD NICHT BENUTZT !
!     real                 :: kappa          ! Eigenfrequenz                   !
!     real                 :: waerme         ! Waermeleitfaehigkeit            !
!     real                 :: eps            ! Genauigkeit                     !
!     real                 :: omega          ! Winkelgeschwindigkeit           !
!     real                 :: r              ! Radius                          !
!     real                 :: u0             !                                 !
!     real                 :: x0,y0          ! Mittelpunkt der Drehung         !
!     real                 :: x1,y1          ! Mittelpunkt des Kegels          !
!     integer              :: controlAus     ! Kontrollausgaben ein/aus        !
!     integer              :: imax           ! Maximale Anzahl Iterationen     !
 
!     character(20)        :: scheme         ! Welches numerische Schema       !
!     character(20)        :: loeser         ! Welcher Gleichungssystemloeser  !
!     character(20)        :: Abbruch        ! Abbruchkriterium fuer elliptic  !
!     character(20)        :: RHS_Typ        ! Art der rechten Seite           !
!     character(20)        :: RB_Typ         ! Art der Randbedingungen         !
!     character(20)        :: AW_Typ         ! Art der Anfangswerte            !
!     character(20)        :: RHS_Datei      ! Datei mit Beschreibung der Quell!

  end type tConstants

  type tLogicals
     logical              :: xflag          ! x-Richtung oder y-Richtung?     !
  end type tLogicals

  type :: tExakt
     real,pointer         :: loesung(:,:)   ! Feld mit der exakten Loesung     !
!      type(c_ptr)         :: loesung        ! loesung(*,*)
     integer              :: kmaxex         ! Genauigkeit der exakten L\F6sung  !
  end type tExakt

  type tFileIO
     CHARACTER(LEN=60)    :: filename       ! Filename                        !
     integer              :: nAusgabe       ! Anzahl Ausgabezeitpunkte        !
     integer              :: bildNr         ! Nr. der naechsten Ausgabezeit   !
     integer              :: visuProg       ! Tecplot or IBM DX               !
  end type tFileIO

  type tElliptic
     real,pointer         :: Ualt(:,:)      ! alte Iterierte                  !
     real,pointer         :: Res(:,:)       ! Residuum je Zelle               !
     integer              :: iter           ! Iterationszaehler               !
  end type tElliptic

  type tHyperbolic
     real,pointer         :: transA(:,:)    ! Transportgeschwindigkeit        !
     real,pointer         :: transB(:,:)    ! Transportgeschwindigkeit        !
  end type tHyperbolic

  type tHyper1D
     real,pointer         :: trans(:)       ! Transportgeschwindigkeit        !
     logical              :: noTrans_Flag   ! alle trans(:) = 0 ?             !
     real,pointer         :: RS(:)          ! Rechte Seite                    !
     real                 :: dx             ! Schrittweite in Raumrichtung    !
     real                 :: dt             ! Schrittweite in Zeitrichtung    !
     real                 :: dxq            ! Kehrweite von dx                !
     real                 :: nraum          ! Anzahl der Gitterpunkte im Raum !
  end type tHyper1D

  type tPolynom
    integer               :: n              ! degree =>  dimension n+1        !
    real,dimension(11)    :: coeffs(1:11)      ! Koeffizienten  order<10         !
  end type tPolynom

  !---------------------------------------------------------------------------!
  public  :: tMeshInfo, tMesh, tMeshGen, tCurv, tZeiten, tRandbedingungen, tConstants, & 
             tLogicals, tCoeff, tExakt, tFileIO, tElliptic, tHyperbolic, tHyper1D, & 
             tNumeric, tPolynom
  !---------------------------------------------------------------------------!

! contains

end module Types
