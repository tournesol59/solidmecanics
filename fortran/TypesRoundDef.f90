module TypesRound

  !---------------------------------------------------------------------------!
  implicit none
  private
  !---------------------------------------------------------------------------!
  type tNode
     real             :: Vect_x, Vect_y, Vect_z
  end type tNode

  type tTrigElem
     integer           :: hasBoundary      !  (0=no, 1=yes)  !
     integer            :: numero
     type(tNode)        :: Node_1, Node_2, Node_3
  end type tTrigElem

  type tRoundMesh        
     integer              :: nElem         ! Number of Elements !
     type(tTrigElem), pointer   :: elems(:)      ! Coordinates of Elements listed   !
  end type tRoundMesh

  type tRandbedingungen
     integer, pointer     :: randFixedIndex(:)      ! RandTyp=feste Biegung/Kurve !
     real,pointer         :: randFixed(:)           ! entsprechende Randwerte   !
     integer, pointer     :: randDisplaceIndex(:)   ! RandTyp=frei/gewuenschte Biegung !
     real,pointer         :: randDisplace(:)        ! entsprechende Randwerte   !
     integer, pointer     :: randMomentIndex(:)     ! RandTyp=Momentum !
     real,pointer         :: randMoment(:)          ! entsprechende Randwerte   !     
  end type tRandbedingungen

end module TypesRound
