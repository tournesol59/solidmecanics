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
!      Verschiebung(1:3,:) ! 3 Verschiebungen in der Mitte der Platte        !
!      TensionTg(1:4,:)    ! array pointer: Spannung tengential              !
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

  subroutine allocateFields(Mesh,RB,Uvar,rhs,Exakt,chicoeff,Const,VarNum)
    
    use TYPES
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tMesh)            :: Mesh       ! Gitterwerte                       !
    type(tRandbedingungen) :: RB         ! Randbedingungen                   !
    real,pointer           :: Uvar(:,:)  ! Feld mit der numerischen Auswertg !=Tr(spanng)!
    real,pointer           :: rhs(:,:)   ! Feld für die Rechte Seite         !
    type(tExakt)           :: Exakt      ! Exakte Loesung                    !
    real,pointer           :: chicoeff(:,:) ! Feld fuer derivat-BasisVektor  !
    type(tConstants)       :: Const      ! Konstanten                        !
    type(tNumeric)         :: VarNum     ! Numerische Felder mit Loesungen   !
    !                                                                        !
    ! Local variable declaration                                             !
    !                                                                        !
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(in)              :: Const
    intent(inout)           :: Mesh,RB,Exakt,chicoeff,VarNum !
    !------------------------------------------------------------------------!


    allocate(Mesh%x(1:Mesh%nraumx*Mesh%nraumy), &
	 Mesh%y(1:Mesh%nraumx*Mesh%nraumy),  &
	 Mesh%z(1:Mesh%nraumx*Mesh%nraumy),  &
         RB%randl(1:Mesh%nraumx), RB%randr(1:Mesh%nraumx), &
         RB%rando(1:Mesh%nraumy), RB%randu(1:Mesh%nraumy), &
         Uvar(1:Mesh%nraumx,1:Mesh%nraumy),                &
         rhs(1:Mesh%nraumx,1:Mesh%nraumy),                 &
         Exakt%loesung(1:Mesh%nraumx,1:Mesh%nraumy),       &
         chicoeff(1:8,1:Mesh%nraumx*Mesh%nraumy), &
         VarNum%Verschiebung(1:3,Mesh%nraumx*Mesh%nraumy), &  
         VarNum%TensionTg(1:3,Mesh%nraumx*Mesh%nraumy), &
         VarNum%Tension3(1:2,Mesh%nraumx*Mesh%nraumy), &
         VarNum%BiegMoment(1:4,Mesh%nraumx*Mesh%nraumy), &
         STAT = allocStat )
     ! TensionTg sigma_xx, sigma_xy=sigma_yx and sigma_yy  !
     ! Tension3 sigma_xz and sigma_yz  (sigma_zz=0)
    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not allocate all variables!'
       STOP
    end if
    
  end subroutine allocateFields



  subroutine deallocateFields(Mesh,RB,Uvar,rhs,Exakt,chicoeff,Const,VarNum)
    
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
    type(tConstants)       :: Const      ! Konstanten                        !
    type(tNumeric)         :: VarNum    ! Numerische Felder mit Loesungen    !
    !                                                                        !
    ! Local variable declaration                                             !
    !                                                                        !
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(in)              :: Const
    intent(inout)           :: Mesh,RB,Exakt,chicoeff,VarNum 
    !------------------------------------------------------------------------!
  
    deallocate(Mesh%x, Mesh%y, Mesh%z, RB%randl, RB%randr, RB%rando, RB%randu, &
         Exakt%loesung, Uvar, rhs, chicoeff, & 
         VarNum%Verschiebung, VarNum%TensionTg, VarNum%Tension3, VarNum%BiegMoment, &
       STAT = allocStat )


    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not deallocate correctly!'
       STOP
    end if

    write(*,*) 
    write(*,*) '============= Program terminated correctly ================ '
    
  end subroutine deallocateFields

end module allocateFields_mod
