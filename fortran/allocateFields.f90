!****************************************************************************!
!                                                                            !
! subroutine allocateFields                                                  !
!                                                                            !
!      Funktion: Bereitstellen des Speicherplatzes                           !
!                                                                            !
!****************************************************************************!
!                                                                            !
!      Mesh%x,Mesh%y            Koordinaten der Gitterpunkte                 !
!      Uvar                     Feld mit der numerischen Lösung              !
!                                                                            !
!****************************************************************************!

module allocateFields_mod
  !--------------------------------------------------------------------------!
  implicit none
  private
  !--------------------------------------------------------------------------!
  interface allocateFields
     module procedure allocateFields
  end interface
  interface deallocateFields
     module procedure deallocateFields
  end interface
  !--------------------------------------------------------------------------!
  public  :: allocateFields, deallocateFields
  !--------------------------------------------------------------------------!

contains

  subroutine allocateFields(Mesh,RB,Uvar,rhs,Exakt,chicoeff,Const)
    
    use TYPES
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tMesh)            :: Mesh       ! Gitterwerte                       !
    type(tRandbedingungen) :: RB         ! Randbedingungen                   !
    real,pointer           :: Uvar(:,:)  ! Feld mit der numerischen Loesung  !
    real,pointer           :: rhs(:,:)   ! Feld für die Rechte Seite         !
    type(tExakt)           :: Exakt      ! Exakte Loesung                    !
    real,pointer           :: chicoeff(:,:) ! Feld fuer derivat-BasisVektor !
    type(tConstants)       :: Const      ! Konstanten                        !
    !                                                                        !
    ! Local variable declaration                                             !
    !                                                                        !
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(in)              :: Const
    intent(inout)           :: Mesh,RB,Exakt,chicoeff
    !------------------------------------------------------------------------!


    allocate(Mesh%x(1:Mesh%nraumx),Mesh%y(1:Mesh%nraumy),  &
         RB%randl(1:Mesh%nraumx), RB%randr(1:Mesh%nraumx), &
         RB%rando(1:Mesh%nraumy), RB%randu(1:Mesh%nraumy), &
         Uvar(1:Mesh%nraumx,1:Mesh%nraumy),                &
         rhs(1:Mesh%nraumx,1:Mesh%nraumy),                 &
         Exakt%loesung(1:Mesh%nraumx,1:Mesh%nraumy),       &
         chicoeff(1:8,1:Mesh%nraumx*Mesh%nraumy), &
         STAT = allocStat )
    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not allocate all variables!'
       STOP
    end if
    
  end subroutine allocateFields



  subroutine deallocateFields(Mesh,RB,Uvar,rhs,Exakt,chicoeff)
    
    use TYPES
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tMesh)            :: Mesh       ! Gitterwerte                       !
    type(tRandbedingungen) :: RB         ! Randbedingungen                   !
    type(tExakt)           :: Exakt      ! Exakte Loesung                    !
    real,pointer           :: Uvar(:,:)  ! Feld mit der numerischen Loesung  !
    real,pointer           :: rhs(:,:)   ! Feld für die Rechte Seite         !
    real,pointer           :: chicoeff(:,:,:) ! Feld fuer derivat-BasisVektor !
    !                                                                        !
    ! Local variable declaration                                             !
    !                                                                        !
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(inout)           :: Mesh,RB,Exakt,chicoeff
    !------------------------------------------------------------------------!
  
    deallocate(Mesh%x, Mesh%y, RB%randl, RB%randr, RB%rando, RB%randu, &
         Exakt%loesung, Uvar, rhs, chicoeff, STAT = allocStat )
    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not deallocate correctly!'
       STOP
    end if

    write(*,*) 
    write(*,*) '============= Program terminated correctly ================ '
    
  end subroutine deallocateFields

end module allocateFields_mod
