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
  use ISO_C_BINDING
  implicit none
  private
  !---------------------------------------------------------------------------!

  type,bind(C) :: tMeshInfo 
     integer(c_int)      :: nraumx         ! Anzahl Gitterpunkte x-Richtung  !
     integer(c_int)      :: nraumy         ! Anzahl Gitterpunkte y-Richtung  !
     real(c_float)       :: startx, starty ! Position der linken unteren Ecke!
     real(c_float)       :: endx, endy     ! Position der rechten oberen Ecke!
  end type tMeshInfo

  type,bind(C) :: tMesh
     integer(c_int)      :: nraumx         ! Copy of tMeshInfo after standard input  !
     integer(c_int)      :: nraumy         ! idem  !
     real(c_float)       :: startx, starty ! idem  !
     real(c_float)       :: endx, endy     ! idem  !
     real(c_float)       :: xlen, ylen     ! Laengen des Rechengebiets in x/y !
     real(c_float)       :: dx             ! Gitterschrittweite (Delta x)    !
     real(c_float)       :: dy             ! Gitterschrittweite (Delta y)    !
     real(c_float)       :: dxq            ! Kehrwert ( 1/ Delta x )         !
     real(c_float)       :: dyq            ! Kehrwert ( 1/ Delta y )         !
     real(c_float)       :: dxx            ! dx^2 (Delta x * Delta x)        !
     real(c_float)       :: dyy            ! dy^2 (Delta y * Delta y )       !
     real(c_float)       :: dxxq           ! Kehrwert des Quadrats (1/dxx)   !
     real(c_float)       :: dyyq           ! Kehrwert des Quadrats (1/dyy)   !
     ! real,pointer (c_ptr) :: x(:)        ! do NOT work
     type(c_ptr)         :: x,y,z
!     type(*)             :: x,y,z          ! to void* in C to Koordinaten Array derGitterpunkte Assumed type at (1) is not allowed for components   !
     !real(c_double),pointer :: x(:),y(:),z(:)
     type(c_ptr)         :: Coefficients   ! to void* to Koeff. des Polynoms fuer vordraengige Interpolation!
!     real(c_double),pointer :: Coefficients(:,:) !Component ‘coefficients’ at (1) cannot have the POINTER attribute because it is a member of the BIND(C) derived type ‘tmesh’

  end type tMesh

  type tMeshGen
     integer              :: nodes          ! Anzahl der Gitterpunkte (Nodes) !
     real,pointer         :: x(:),y(:),z(:) ! Koordinaten der Gitterpunkte    !
     integer              :: elmts          ! Anzahl der Quadranglen          !
     integer,pointer      :: quad(:,:)      ! Indiz der Punte fuer Quadrangle !
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
     real(C_DOUBLE),pointer         :: chicoeff(:,:)          ! 8 Values fuer jedes Punktes: Christoffel Koeffizienten
     real(C_DOUBLE),pointer         :: der_a1xx_x(:),der_a1xx_y(:),der_a1xy_x(:),der_a1xy_y(:), &
                          der_a2yx_x(:),der_a2yx_y(:),der_a2yy_x(:),der_a2yy_y(:)  
  end type tCurv

  type tNumeric
    real(C_DOUBLE),pointer     :: Verschiebung(:,:) ! 3 Verschiebungen in der Mitte der Platte !
    real(C_DOUBLE),pointer     :: TensionTg(:,:)    ! array pointer: Spannung tengential Tensor!
    real(C_DOUBLE),pointer     :: Tension3(:,:)     ! array pointer: Spannung in Richtung 3  !
    real(C_DOUBLE),pointer     :: BiegMoment(:,:)   ! BiegMoment Tensor, nicht tangential    !

     real,pointer         :: Coefficients(:,:) ! des Polynoms fuer vordraengige Interpolation! 

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

  type,bind(C) :: tRandbedingungen
! 'trandbedingungen' is not C interoperable --> may be add type, bind(C)

     integer(c_int)          :: randltype      ! RandTyp links: 1=Fixed, 2=Momentum, 3=Freibiegung
     type(c_ptr)             :: randl       ! Randwerte links                 !
     integer                 :: randrtype      ! RandTyp rechts: 1=Fixed, 2=Momentum, 3=Freibiegung
     type(c_ptr)             :: randr       ! Randwerte rechts                !
     integer                 :: randotype      ! RandTyp oben: 1=Fixed, 2=Momentum, 3=Freibiegung
     type(c_ptr)             :: rando      ! Randwerte oben                  !
     integer                 :: randutype      ! RandTyp unten: 1=Fixed, 2=Momentum, 3=Freibiegung
     type(c_ptr)             :: randu       ! Randwerte unten                 !
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

  type,bind(C) :: tExakt
!     real(c_float),pointer         :: loesung(:,:)   ! Feld mit der exakten Loesung     !
      type(c_ptr)         :: loesung        ! loesung(*,*)
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
