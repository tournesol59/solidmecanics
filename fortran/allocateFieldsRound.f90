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
  interface allocateFieldsBCSRound
     module procedure allocateFieldsBCSRound
  end interface
  interface allocateFieldsVarRound
     module procedure allocateFieldsVarRound
  end interface

  interface deallocateFieldsRound
     module procedure deallocateFieldsRound
  end interface
!  interface deallocateFieldsCurv
!     module procedure deallocateFieldsCurv
!  end interface
  interface deallocateFieldsBCSRound
     module procedure deallocateFieldsBCSRound
  end interface
  interface deallocateFieldsVarRound
     module procedure deallocateFieldsVarRound
  end interface

  !--------------------------------------------------------------------------!
  public  :: allocateFieldsRound, allocateFieldsVarRound, allocateFieldsBCSRound, &
             deallocateFieldsRound, deallocateFieldsVarRound, deallocateFieldsBCSRound
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
                                          !=Tr(spanng)                       
    ! Local variable declaration                                             !
    !                                                                        !
    integer                 :: nn,ne,i,allocStat      !
    !------------------------------------------------------------------------!
    intent(inout)           :: MeshR,ExaktR                                  !
    intent(in)              :: NTestMeshR                                    !
    !------------------------------------------------------------------------!

    nn = NTestMeshR%nNode
    ne = NTestMeshR%nElem

    allocate(MeshR%nodes(1:nn), MeshR%elems(1:nn), MeshR%neighbours(1:nn), &            
         ExaktR%loesung(1:nn,1:NumberofDisp),       &
         STAT = allocStat )        
! Below may not be used with simple const dimensions, or declare allocatable
    do i=1,nn
       allocate(MeshR%nodes(i)%vects(1:3))
    enddo 
    do i=1,ne
       allocate(MeshR%elems(i)%numeros(1:4))
       allocate(MeshR%neighbours(i)%numeros(1:4))
   enddo 
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsRound: Could not allocate all variables!'
       STOP
    end if

  
  end subroutine allocateFieldsRound

!******************************************************
  subroutine allocateFieldsBCSRound(NTestMeshR,BCU,BCS) 
! bind (C, name="allocateFieldsBCSRound")
!    use ISO_C_BINDING, only: c_char
     use TypesRound
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tRoundMeshInfo)    :: NTestMeshR  ! Dimensions                      !                     !
    type(tRoundRandU)               :: BCU       ! Randbedingungen Biegg     !
    type(tRoundRandS)               :: BCS       ! Randbedingungen Kraft     !
    ! Local variable declaration                                             !
    !                                                                        !
    integer                 :: nxmx,nymx,i,allocStat                         !
    !------------------------------------------------------------------------!                                 !
    intent(inout)           ::  BCU,BCS                                      !
    intent(in)              :: NTestMeshR                                    !
    !-------------------------------------------------------------------------!

    nxmx = NTestMeshR%nxmx
    nymx = NTestMeshR%nymx

    allocate(BCU%randFixedIndex(1:nxmx),BCU%randFixed(1:nxmx), &            
            BCU%randDisplaceIndex(1:nxmx),BCU%randDisplace(1:nxmx), &
            BCS%randForceIndex(1:nxmx),BCS%randForce(1:nxmx), &            
            BCS%randMomentIndex(1:nxmx),BCS%randMoment(1:nxmx), &
         STAT = allocStat ) 
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsBCSRound: Could not allocate all variables!'
       STOP
    end if
  
  end subroutine allocateFieldsBCSRound

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
    ne = NTestMeshR%nElem ! A priori je nach Element geordnet, node 1(ux,vy,wz,Tx,Ty,.), node 2, node 3

    allocate(VarNumR%displacement(1:ne,1:NumberofDisp),       &
             VarNumR%rhs(1:ne,1:NumberofDisp),       &
             VarNumR%tension(1:ne,1:NumberofTension),        &
             VarNumR%biegmoment(1:ne,1:NumberofMom),        &
         STAT = allocStat )        
 
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsVarRound: Could not allocate all variables!'
       STOP
    end if

  end subroutine allocateFieldsVarRound

!******************************************************
  subroutine allocateFieldsCoeffRound(NTestMeshR,MeshCoeffR) 
! bind (C, name="allocateFieldsCoeffRound")
!    use ISO_C_BINDING, only: c_char 
  use TypesRound

  implicit none
  ! Variable Deklaration  ---------------------------------------------------!
    type(tRoundMeshInfo)   :: NTestMeshR  ! Dimensions                       !
    type(tRoundCoeff)      :: MeshCoeffR     ! Numerische Felder mit Loesungen  !
    !------------------------------------------------------------------------!
    integer                 :: allocStat,ne,nn                               !
    !------------------------------------------------------------------------!
    intent(inout)           :: MeshCoeffR                                       !
    intent(in)              :: NTestMeshR                                    !
    !------------------------------------------------------------------------!

    nn = NTestMeshR%nNode
    ne = NTestMeshR%nElem ! A priori interpol order 2 => 6 coeffs per Elements

    allocate(MeshCoeffR%Coefficients(1:ne*6),       &
         STAT = allocStat )        
 
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsCoeffRound: Could not allocate all variables!'
       STOP
    end if

  end subroutine allocateFieldsCoeffRound

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
! Below cannot be not used if simple const dimensions
    do i=1,nn
       deallocate(MeshR%nodes(i)%vects)
    enddo 
    do i=1,ne
       deallocate(MeshR%elems(i)%numeros)
       deallocate(MeshR%neighbours(i)%numeros)
   enddo 
    deallocate(MeshR%nodes, MeshR%elems, MeshR%neighbours, &
              ExaktR%loesung, &
               STAT=allocStat)
    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsRound: Could not deallocate correctly!'
       STOP
    else
       print *, 'Terminated AllocateFieldsRound: OK'
    end if

  end subroutine deallocateFieldsRound

!******************************************************
  subroutine deallocateFieldsBCSRound(NTestMeshR,BCU,BCS) 
! bind (C, name="deallocateFieldsBCSRound")
!    use ISO_C_BINDING, only: c_char
    use TypesRound
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tRoundMeshInfo)    :: NTestMeshR  ! Dimensions                      !                     !
    type(tRoundRandU)               :: BCU       ! Randbedingungen Biegg     !
    type(tRoundRandS)               :: BCS       ! Randbedingungen Kraft     !
    ! Local variable declaration                                             !
    !                                                                        !
    integer                 :: nxmx,nymx,i,allocStat                         !
    !------------------------------------------------------------------------!                                 !
    intent(inout)           ::  BCU,BCS                                      !
    intent(in)              :: NTestMeshR                                    !
    !-------------------------------------------------------------------------!

    nxmx = NTestMeshR%nxmx
    nymx = NTestMeshR%nymx

    deallocate(BCU%randFixedIndex,BCU%randFixed, &            
            BCU%randDisplaceIndex,BCU%randDisplace, &
            BCS%randForceIndex,BCS%randForce, &            
            BCS%randMomentIndex,BCS%randMoment, &
         STAT = allocStat ) 
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsBCSRound: Could not deallocate all variables!'
       STOP
    end if
  
  end subroutine deallocateFieldsBCSRound

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
    else
       print *, 'Terminated AllocateFieldsVARRound: OK'
    end if

  end subroutine deallocateFieldsVarRound

!******************************************************
  subroutine deallocateFieldsCoeffRound(MeshCoeffR) 
! bind (C, name="deallocateFieldsVarNum")
!    use ISO_C_BINDING, only: c_char  

    use TypesRound
    
    implicit none

  ! Variable Deklaration  ---------------------------------------------------!
    type(tRoundCoeff)       :: MeshCoeffR     ! Numerische Felder Interpol Coeff  !
    !------------------------------------------------------------------------!
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(inout)           :: MeshCoeffR                                    !
    !------------------------------------------------------------------------!

    deallocate(MeshCoeffR%Coefficients, &
               STAT=allocStat)

    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsCoeffRound: Could not deallocate correctly!'
       STOP
    end if

  end subroutine deallocateFieldsCoeffRound

end module allocateFieldsRound_mod
