module TypesRoundDef

  !---------------------------------------------------------------------------!
  implicit none
  private
  !---------------------------------------------------------------------------!
  type tNode
     real             :: Vect_x, Vect_y,Vect_z
  end type tPoint

  type tTrigElem
     boolean            :: hasBoundary
     integer            :: numero
     real               :: Node_1, Node_2, Node_3
  end type tTrigElem

  type tRoundMesh        
     integer              :: nElem         ! Number of Elements
     tTrigElem, pointer   :: elems(:)      ! Coordinates of Elements listed   !
  end type tMesh

  type tRandbedingungen
     real,pointer         :: randl(:)       ! Randwerte links                 !
     real,pointer         :: randr(:)       ! Randwerte rechts                !
     real,pointer         :: rando(:)       ! Randwerte oben                  !
     real,pointer         :: randu(:)       ! Randwerte unten                 !
  end type tRandbedingungen

end TypesRoundDef
