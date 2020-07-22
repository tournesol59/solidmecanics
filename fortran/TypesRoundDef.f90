module TypesRound
! This module is the center of definition of all data types used by the program
! from mesh definition, to solution fields to be exported
  !---------------------------------------------------------------------------!
!  use ISO_C_BINDING
  implicit none
  private
  !---------------------------------------------------------------------------!
  type tNode
!     real,dimension(3)             :: vects
      real,dimension(:),allocatable   :: vects(:)   ! previewed for dim=3
  end type tNode

  type tTrigElem
!     integer           :: hasBoundary      !  not useful, see field neighbours on next data !
!     integer,dimension(4)      :: numeros ! Triangle elmt numero, node_1, node_2, node_3  ! 
!     integer       :: numero, node_1, node_2, node_3
      integer,dimension(:),allocatable  :: numeros(:)  ! previewed for dim=4 ^
  end type tTrigElem

  type tRoundMeshInfo
     integer            :: nNode         ! Number of Nodes !
     integer            :: nElem         ! Number of Elements, equal to num. neighbours data !
     integer            :: nBound        ! number of frontiers that may be boundary conditions
     integer            :: nxmx, nymx    ! Max Bound for Seiten Index of Vertex (Boundary Condition)
     integer,dimension(4) :: BoundVert   ! Number of elements faces for 4 vertex (previewed for square mesh, but can be used also for round mesh)
  end type tRoundMeshInfo

  type tRoundMesh        
!     type(tNode), pointer       :: nodes(:)   ! Coordinates of Nodes listed !
! Alternative
     type(tNode),dimension(:),allocatable   :: nodes(:)  ! Same coordinates but can be allocated, nNodes x3 values
     type(tTrigElem),dimension(:),allocatable   :: elems(:)      ! Nodes index Elements listed   !
     type(tTrigElem),dimension(:),allocatable   :: neighbours(:) ! Re-use of structure type tTrigElem !
                                                 ! to store neighbours, -1 = has boundary
  end type tRoundMesh

   type tRoundRandU
     integer, pointer     :: randFixedIndex(:)      ! RandTyp=feste Biegung/Kurve !
     real,pointer         :: randFixed(:)           ! entsprechende Randwerte   !
     integer, pointer     :: randDisplaceIndex(:)   ! RandTyp=frei/gewuenschte Biegung !
     real,pointer         :: randDisplace(:)        ! entsprechende Randwerte   !   
  end type tRoundRandU

  type tRoundRandS
     integer, pointer     :: randForceIndex(:)      ! RandTyp=Kraft in Richtung parallel Platen,Scheibe
     real, pointer        :: randForce(:)           ! entsprechende Randwerte
     integer, pointer     :: randMomentIndex(:)     ! RandTyp=Momentum !
     real,pointer         :: randMoment(:)          ! entsprechende Randwerte   !  
  end type tRoundRandS

  type tPolynom
    integer               :: n              ! degree =>  dimension n+1        !
    real,dimension(11)    :: coeffs(1:11)      ! Koeffizienten  order<10         !
  end type tPolynom

  type tXYPolynom ! definiert in TypesDef
!    integer               :: n              ! degree, angenommen 2 hier       !
    real,dimension(10)    :: coeffs(1:10)     ! Koeffizienten order 1,x,y,x2,xy,y2! 
                                              ! besetzt       !
  end type tXYPolynom

  type tPassageMatrix            ! dient als passage matrix
    real,dimension(9)     :: mat(1:9)
  end type tPassageMatrix

  type tRoundlocal
    integer,dimension(3)   :: permute  ! setzt aus der Punkten (j,k,l), 
                                       ! was Punkt 1, Punkt 2 oder Punkt 3 ist
    real,dimension(4)      :: triplet  !  x2,x3,y3 Werten
    real,dimension(3)      :: dist     ! dist(1) ist der Laenge des Segments ggseitig zum nodes 1
    real,dimension(3)      :: angleEuler ! welche benoetigt werden, um Passage aus lokal
                                         ! zu global Basis zu machen
    type(tPassageMatrix)   :: axes_xyz 
  end type tRoundlocal

  type tRoundCoeff        ! benutzbar aber muss allokiert werden und a priori 
                          ! interpol order 2 => 1:6 coeffs per Elements
     real,pointer     :: Coefficients(:,:) ! 1: Matrix coeffs von Biegung zur Deformationen !
                                           ! 2: oder Coeffs des Polynoms fuer Interpolation, 1 fuer jedes Elements! 
!     real,pointer                 :: Jacobian(:)    ! 2x2 real, nur diesen Test vorgesehen !
!     logical,pointer              :: Test(:)        ! Resultat des Tests an dem Element !
  end type tRoundCoeff

  type tRoundNumeric ! Numerische Loesung fuer Post treatment
    real,pointer     :: displacement(:,:) ! Default 12 Verschiebungen w in der Mitte des Elements !
                   ! IMPORTANT CORRECTION SINCE version d4cd9f5 : (not here but preview in allocatefields)
                   ! before: 12 Verschiebungen per ELEMENT: probably (uz=w, rx,ry, rz)x3 and ux,uy somewhere else 
                   ! after : 6 Verschiebungen per NODE:    (ux,uy,uz,rx,ry,rz) 6 Verschiebungen
                   ! note that rz is not part of the shell theory but only associated to a fictive lower rigidity

    real,pointer     :: rhs(:,:)          ! rechte Seite Vektor, gleiche Dimension als Verschiebung !
    real,pointer     :: tension(:,:)    ! array pointer: Spannung tengential Tensor in Element!
!    real,pointer     :: Tension3(:,:)     ! array pointer: Spannung in Richtung normal  !
    real,pointer     :: biegmoment(:,:)   ! BiegMoment Tensor, nicht tangential    !

  end type tRoundNumeric

  type tRoundExakt
     real,pointer         :: loesung(:,:)        ! loesung(*,*)
     integer              :: kmaxex         ! Genauigkeit der exakten Loesung  !
  end type tRoundExakt

  type tSparseYMat
     integer                        :: nnz, rows, lband
     real,dimension(:),allocatable  :: a(:)
     integer,dimension(:),allocatable :: ia(:)
     integer,dimension(:),allocatable :: ja(:)
  end type tSparseYMat

!************** Mixed Fortran-C programming data exchange are based on fixed size array here
!              (no pointer solution found)    
!*******************************************************************************************
! input to C:
! Array items for block passed to C (50 passed items to "c_testconnrigidmatrix" each to form a LL)

! output from C:
  type tLLSparseMat                ! Linked List from C-procedure transformed in a Array of 50 items tLLSparseMat
    integer             :: nEdgeNr   ! nEdgeNew new number if sparse renumerotation algorithms were used
    integer             :: node_1, node_2
    integer             :: IndexI1                !index i+1 or -1 if branch of the tree is ended
    integer             :: nextb1, nextb2         !Index nextb1, nextb2 (max two) of possible two other connections
  end type tLLSparseMat

  !---------------------------------------------------------------------------!
  public  :: tNode, tTrigElem, tRoundMeshInfo, tRoundMesh, tRoundRandU, tRoundRandS, &
             tRoundNumeric, tRoundCoeff, tRoundExakt, tPolynom, tXYPolynom, &
             tRoundlocal, tPassageMatrix, tSparseYMat, tLLSparseMat
  !---------------------------------------------------------------------------!

end module TypesRound
