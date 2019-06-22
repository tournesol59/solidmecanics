!****************************************************************************!
!                                                                            !
! subroutine allocateFields                                                  !
!                                                                            !
!      Funktion: Bereitstellen des Speicherplatzes                           !
!                                                                            !
!****************************************************************************!
!                                                                            !
!      Mesh%nodes,Mesh%elems,Mesh%neighbours:Koordinaten der Punkten und Ind.!
!      Uvar:                Feld mit u.a. der numerischen Lösung an Elems i  !
!        Displacement(1:12,:) 3 Verschiebungen in der Mitte der Platte       !
!        Tension(1:6,:)    array pointer: Spannung tengential                !
!        ! TensionT3(1,:)    ! array pointer: Spannung normal despite that E,xz=0 !
!        BiegMoment(1:6,:)   ! array pointer: contains as last arrays above    !
!                           the inital value of the problem and then re-calculated.!
!****************************************************************************!

module allocateFieldsRound_mod
!  use ISO_C_BINDING
!  use TypesRound

  !--------------------------------------------------------------------------!
  implicit none
  private
  !--------------------------------------------------------------------------!

  integer, parameter  :: NumberofDisp = 12       ! dof, so kann es geaendert werden !
  integer, parameter  :: NumberofTension = 6       ! tang. Kraft, so kann es geaendert werden !
  integer, parameter  :: NumberofMom = 6       ! bieg. Momentums, so kann es geaendert werden !

  !--------------------------------------------------------------------------!
  interface allocateFieldsRound
     module procedure allocateFieldsRound
  end interface
 !  interface allocateFieldsCurv
 !     module procedure allocateFieldsCurv
 !  end interface
!  interface allocateFieldsVarNum
!     module procedure allocateFieldsVarNum
!  end interface

  interface deallocateFieldsRound
     module procedure deallocateFieldsRound
  end interface
!  interface deallocateFieldsCurv
!     module procedure deallocateFieldsCurv
!  end interface
  interface deallocateFieldsVarRound
     module procedure deallocateFieldsVarRound
  end interface

  !--------------------------------------------------------------------------!
  public  :: allocateFieldsRound, allocateFieldsVarRound, &
             deallocateFieldsRound, deallocateFieldsVarRound
  !--------------------------------------------------------------------------!

contains


!******************************************************
  subroutine allocateFieldsRound(NTestMeshR,MeshR,ExaktR) 
! bind (C, name="allocateFieldsRound")
!    use ISO_C_BINDING, only: c_char
    use TypesRound
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tRoundMeshInfo)    :: NTestMeshR  ! Dimensions                      !
    type(tRoundMesh)        :: MeshR       ! Gitterwerte                     !
    type(tRoundExakt)       :: ExaktR      ! Exakte Loesung Triangle //Round !
!    real,pointer           :: Uvar(:,:)  ! Feld mit der numerischen Auswertg!
                                          !=Tr(spanng)                       !
    !                                                                        !
    ! Local variable declaration                                             !
    !                                                                        !
    integer                 :: nn,ne,i,allocStat,allocStati                  !
    !------------------------------------------------------------------------!
    intent(inout)           :: MeshR,ExaktR                                  !
    intent(in)              :: NTestMeshR                                    !
    !------------------------------------------------------------------------!

    nn = NTestMeshR%nNode
    ne = NTestMeshR%nElem

    allocate(MeshR%nodes(1:nn), MeshR%elems(1:nn), MeshR%neighbours(1:nn), &            
         ExaktR%loesung(1:nn,1:NumberofDisp),       &
         STAT = allocStat )        
!    do i=1,nn
!       allocate(MeshR%nodes(i))
!    enddo 
!    do i=1,ne
!       allocate(MeshR%elems(i))
!       allocate(MeshR%neighbours(i))
!   enddo 
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsRound: Could not allocate all variables!'
       STOP
    end if
    
  end subroutine allocateFieldsRound

!******************************************************
  subroutine allocateFieldsVarRound(NTestMeshR,VarNumR) 
! bind (C, name="allocateFieldsVarRound")
!    use ISO_C_BINDING, only: c_char 
  use TypesRound

  implicit none
  ! Variable Deklaration  ---------------------------------------------------!
    type(tRoundMeshInfo)   :: NTestMeshR  ! Dimensions                       !
    type(tRoundNumeric)    :: VarNumR     ! Numerische Felder mit Loesungen  !
    !------------------------------------------------------------------------!
    integer                 :: allocStat,ne,nn                               !
    !------------------------------------------------------------------------!
    intent(inout)           :: VarNumR                                       !
    intent(in)              :: NTestMeshR                                    !
    !------------------------------------------------------------------------!

    nn = NTestMeshR%nNode
    ne = NTestMeshR%nElem

    allocate(VarNumR%displacement(1:nn,1:NumberofDisp),       &
             VarNumR%rhs(1:nn,1:NumberofDisp),       &
             VarNumR%tension(1:nn,1:NumberofTension),        &
             VarNumR%biegmoment(1:nn,1:NumberofMom),        &
         STAT = allocStat )        
 
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsVarRound: Could not allocate all variables!'
       STOP
    end if

  end subroutine allocateFieldsVarRound


!******************************************************
  subroutine deallocateFieldsRound(NTestMeshR,MeshR,ExaktR) 
! bind (C, name="deallocateFields")
!    use ISO_C_BINDING, only: c_char  
    use TypesRound
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tRoundMeshInfo)    :: NTestMeshR  ! Dimensions                      !
    type(tRoundMesh)        :: MeshR       ! Gitterwerte                     !
    type(tRoundExakt)       :: ExaktR      ! Exakte Loesung Triangle //Round !
    !                                                                        !
    ! Local variable declaration                                             !
    integer                 :: nn,ne,i,allocStat,allocStati                  !
    !------------------------------------------------------------------------!
    intent(inout)           :: MeshR,ExaktR                                  !
    intent(in)              :: NTestMeshR                                    !
    !------------------------------------------------------------------------!

    nn = NTestMeshR%nNode
    ne = NTestMeshR%nElem

!    do i=1,nn
!       deallocate(MeshR%nodes(i))
!    enddo 
!    do i=1,ne
!       deallocate(MeshR%elems(i))
!       deallocate(MeshR%neighbours(i))
!   enddo 
    deallocate(MeshR%nodes, MeshR%elems, MeshR%neighbours, &
              ExaktR%loesung, &
               STAT=allocStat)
    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsRound: Could not deallocate correctly!'
       STOP
    end if

  end subroutine deallocateFieldsRound

!******************************************************
  subroutine deallocateFieldsVarRound(VarNumR) 
! bind (C, name="deallocateFieldsVarNum")
!    use ISO_C_BINDING, only: c_char  

    use TypesRound
    
    implicit none

  ! Variable Deklaration  ---------------------------------------------------!
    type(tRoundNumeric)    :: VarNumR     ! Numerische Felder mit Loesungen  !
    !------------------------------------------------------------------------!
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(inout)           :: VarNumR                                       !
    !------------------------------------------------------------------------!

    deallocate(VarNumR%displacement, VarNumR%rhs, &
               VarNumR%tension, VarNumR%biegmoment, &
               STAT=allocStat)

    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsVarRound: Could not deallocate correctly!'
       STOP
    end if

  end subroutine deallocateFieldsVarRound


end module allocateFieldsRound_mod
