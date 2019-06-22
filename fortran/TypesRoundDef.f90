module TypesRound

  !---------------------------------------------------------------------------!
!  use ISO_C_BINDING
  implicit none
  private
  !---------------------------------------------------------------------------!
  type tNode
     real             :: vect_x, vect_y, vect_z
  end type tNode

  type tTrigElem
!     integer           :: hasBoundary      !  not useful, see field neighbours on next data !
     integer       :: numero
!     type(tNode)        :: node_1, node_2, node_3  ! not useful, two refs data struct instead !
     integer       :: node_1, node_2, node_3
  end type tTrigElem

  type tRoundMeshInfo
     integer            :: nNode         ! Number of Nodes !
     integer            :: nElem         ! Number of Elements, equal to num. neighbours data !
  end type tRoundMeshInfo

  type tRoundMesh        
     type(tNode), pointer       :: nodes(:)   ! Coordinates of Nodes listed !
! Alternative
!    type(tNode),dimension(:),allocatable   :: nodes
     type(tTrigElem), pointer   :: elems(:)      ! Nodes index Elements listed   !
     type(tTrigElem), pointer   :: neighbours(:) ! Re-use of structure type tTrigElem !
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

  type tXYPolynom ! definiert in TypesDef
!    integer               :: n              ! degree, angenommen 2 hier       !
    real,dimension(10)    :: coeffs(1:10)      ! Koeffizienten order 1,x,y,x2,xy,y2! 
                                              ! besetzt       !
  end type tXYPolynom

  type tRoundCoeff        ! Nicht recht benutztbar
     type(tXYPolynom),pointer     :: Coefficients(:) ! des Polynoms fuer vordraengige !
                                                     ! Interpolation, 1 fuer jedes Elements! 
!     real,pointer                 :: Jacobian(:)    ! 2x2 real, nur diesen Test vorgesehen !
!     logical,pointer              :: Test(:)        ! Resultat des Tests an dem Element !
  end type tRoundCoeff

  type tRoundNumeric ! Numerische Loesung fuer Post treatment
    real,pointer     :: displacement(:,:) ! Default 12 Verschiebungen w in der Mitte des Elements !
    real,pointer     :: rhs(:,:)          ! rechte Seite Vektor, gleiche Dimension als Verschiebung !
    real,pointer     :: tension(:,:)    ! array pointer: Spannung tengential Tensor in Element!
!    real,pointer     :: Tension3(:,:)     ! array pointer: Spannung in Richtung normal  !
    real,pointer     :: biegmoment(:,:)   ! BiegMoment Tensor, nicht tangential    !

  end type tRoundNumeric

  type tRoundExakt
     real,pointer         :: loesung(:,:)        ! loesung(*,*)
     integer              :: kmaxex         ! Genauigkeit der exakten Loesung  !
  end type tRoundExakt

  !---------------------------------------------------------------------------!
  public  :: tNode, tTrigElem, tRoundMeshInfo, tRoundMesh, tRoundRandU, tRoundRandS, &
             tRoundNumeric, tRoundCoeff, tRoundExakt, tXYPolynom
  !---------------------------------------------------------------------------!

end module TypesRound
