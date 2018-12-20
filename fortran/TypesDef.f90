!*****************************************************************************!
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
  implicit none
  private
  !---------------------------------------------------------------------------!
  type tMesh
     integer              :: nraumx         ! Anzahl Gitterpunkte x-Richtung  !
     integer              :: nraumy         ! Anzahl Gitterpunkte y-Richtung  !
     real                 :: startx, starty ! Position der linken unteren Ecke!
     real                 :: endx, endy     ! Position der rechten oberen Ecke!
     real                 :: xlen, ylen     ! Laengen des Rechengebiets in x/y !
     real                 :: dx             ! Gitterschrittweite (Delta x)    !
     real                 :: dy             ! Gitterschrittweite (Delta y)    !
!     real                 :: dxq            ! Kehrwert ( 1/ Delta x )         !
!     real                 :: dyq            ! Kehrwert ( 1/ Delta y )         !
!     real                 :: dxx            ! dx^2 (Delta x * Delta x)        !
!     real                 :: dyy            ! dy^2 (Delta y * Delta y )       !
!     real                 :: dxxq           ! Kehrwert des Quadrats (1/dxx)   !
!     real                 :: dyyq           ! Kehrwert des Quadrats (1/dyy)   !

     real,pointer         :: x(:),y(:),z(:) ! Koordinaten der Gitterpunkte    !

     real,pointer         :: Coefficients(:,:) ! des Polynoms fuer vordraengige Interpolation!   
  end type tMesh

  type tMeshGen  ! Hypothesis Nmax=5 different Randbedingungen
     character(6)         :: elmttype       ! Typ der Elementen, hier Dreiecke "BE3***" !
     integer              :: nodes          ! Anzahl der Gitterpunkte (Nodes) !
     real,pointer         :: x(:),y(:),z(:) ! Koordinaten der Gitterpunkte    !
     integer              :: elmts          ! Anzahl der Elemente             !
     integer,pointer      :: allelmt(:,:)   ! Indiz der Punten fuer Dreiecken (3,:) !
     integer, pointer     :: neighbours(:,:) ! Indiz der Nachbarn-Elementen fuer Elmt i !
                                             ! Max Anzahl von Nachbarn ist hier 10 !
     integer, pointer     :: frontelmt(:,:) ! Indiz des Elements an dem i-ten (1:5,:) Wand !
                                      ! (ein max. Anzahl von 5 Waende ist hier vorgesehen) !
     integer,pointer      :: frontvert(:,:,:) ! Indiz der 2xPunkten an dem i-ten Wand 
                                              ! (1:5,1:2,:) !
     integer,dimension(2,5) :: fronttype    ! Typ der Randbedingungen und Anzahl der Elmt !
                                            ! an dem i-ten (1:5) Wand !
                             ! frontype(1,i) = RB Typ(1-4, 1=Moment gesetzt, 2=Force gesetzt!
                             ! , 3=Biegung gesetzt, 4=Displacement gesetzt, 5=eine Kombination !
                             ! von Moment u . Force, 6=eine Kombination von Biegung u. ! 
                             ! Displacement sind vorgesehen, in dem Fall werden
                             !  die Biegung (resp. Moment) Werten in frontvalue(i,3)   !
                             ! (resp. frontvalue(i,1)) und die Displacement (resp. Force) !
                             ! Werten in frontvalue(i,4) (resp. frontvalue(i,2)) gespeichert. !
                             ! fronttype(2,i) = Anzahl von Elementen, die an dem Wand/Rand i sind!
     real,pointer         :: frontvalue(:,:,:)  
                             ! Randwerten, je nach Typ der RB in fronttype (1:5,1:6,:) !
  end type tMeshGen 

  type tCurv              ! Abhaengig type tMesh fuer (nraumx, nraumy) !
     real,pointer         :: chicoeff(:,:)          ! 8 Values fuer jedes Punktes: Christoffel Koeffizienten
     real,pointer         :: der_a1xx_x(:),der_a1xx_y(:),der_a1xy_x(:),der_a1xy_y(:), &
                          der_a2yx_x(:),der_a2yx_y(:),der_a2yy_x(:),der_a2yy_y(:)  
  end type tCurv

  type tNumeric
    real, pointer     :: Verschiebung(:,:) ! 3 Verschiebungen in der Mitte der Platte !
    real, pointer     :: TensionTg(:,:)    ! array pointer: Spannung tengential Tensor!
    real, pointer     :: Tension3(:,:)     ! array pointer: Spannung in Richtung 3  !
    real, pointer     :: BiegMoment(:,:)   ! BiegMoment Tensor, nicht tangential    !
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

  type tRandbedingungen
     integer              :: randltype      ! RandTyp links: 1=Fixed, 2=Momentum, 3=Freibiegung
     real,pointer         :: randl(:)       ! Randwerte links                 !
     integer              :: randrtype      ! RandTyp rechts: 1=Fixed, 2=Momentum, 3=Freibiegung
     real,pointer         :: randr(:)       ! Randwerte rechts                !
     integer              :: randotype      ! RandTyp oben: 1=Fixed, 2=Momentum, 3=Freibiegung
     real,pointer         :: rando(:)       ! Randwerte oben                  !
     integer              :: randutype      ! RandTyp unten: 1=Fixed, 2=Momentum, 3=Freibiegung
     real,pointer         :: randu(:)       ! Randwerte unten                 !
  end type tRandbedingungen
  
  type tConstants  ! Original for Waermetransportsproblem geplant !
     real                 :: Pi             ! Pi                              !
     real                 :: EY             ! Young Module                    !
     real                 :: nu             ! Poisson coefficient             !
!     real                 :: h              ! Groesse (Height), Konstant      !
     integer              :: auto           ! Automatishes Auswahl eines vordefinierten Mesh !
         ! REST WIRD NICHT BENUTZT !
     real                 :: kappa          ! Eigenfrequenz                   !
     real                 :: waerme         ! Waermeleitfaehigkeit            !
     real                 :: eps            ! Genauigkeit                     !
     real                 :: omega          ! Winkelgeschwindigkeit           !
     real                 :: r              ! Radius                          !
     real                 :: u0             !                                 !
     real                 :: x0,y0          ! Mittelpunkt der Drehung         !
     real                 :: x1,y1          ! Mittelpunkt des Kegels          !
     integer              :: controlAus     ! Kontrollausgaben ein/aus        !
     integer              :: imax           ! Maximale Anzahl Iterationen     !
     integer              :: type           ! Welches numerische Schema       !
     character(20)        :: scheme         ! Welches numerische Schema       !
     character(20)        :: loeser         ! Welcher Gleichungssystemloeser  !
     character(20)        :: Abbruch        ! Abbruchkriterium fuer elliptic  !
     character(20)        :: RHS_Typ        ! Art der rechten Seite           !
     character(20)        :: RB_Typ         ! Art der Randbedingungen         !
     character(20)        :: AW_Typ         ! Art der Anfangswerte            !
     character(20)        :: RHS_Datei      ! Datei mit Beschreibung der Quell!

  end type tConstants

  type tLogicals
     logical              :: xflag          ! x-Richtung oder y-Richtung?     !
  end type tLogicals

  type tExakt
     real,pointer         :: loesung(:,:)   ! Feld mit der exakten L\F6sung     !
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
  public  :: tMesh, tMeshGen, tCurv, tZeiten, tRandbedingungen, tConstants, tLogicals, &
             tExakt, tFileIO, tElliptic, tHyperbolic, tHyper1D, tNumeric, tPolynom
  !---------------------------------------------------------------------------!

! contains

end module Types
