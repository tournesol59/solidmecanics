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
!      TensionT3(1,:)    ! array pointer: Spannung normal despite that E,xz=0 !
!      BiegMoment(1:3,:)   !array pointer: contains as last arrays above    !
!                           the inital value of the problem and then re-calculated.!
!****************************************************************************!

module allocateFields_mod

  use types

  !--------------------------------------------------------------------------!
  implicit none
  private
  !--------------------------------------------------------------------------!
  interface allocateFields
     module procedure allocateFields
  end interface
  interface allocateFieldsCurv
     module procedure allocateFieldsCurv
  end interface
  interface allocateFieldsVarNum
     module procedure allocateFieldsVarNum
  end interface
  interface allocateGitter2D
     module procedure allocateGitter2D
  end interface

  interface deallocateFields
     module procedure deallocateFields
  end interface
  interface deallocateFieldsCurv
     module procedure deallocateFieldsCurv
  end interface
  interface deallocateFieldsVarNum
     module procedure deallocateFieldsVarNum
  end interface
  interface deallocateGitter2D
     module procedure deallocateGitter2D
  end interface
  !--------------------------------------------------------------------------!
  public  :: allocateFields, allocateFieldsCurv, allocateFieldsVarNum, allocateGitter2D, &
             deallocateFields, deallocateFieldsCurv, deallocateFieldsVarNum, deallocateGitter2D
  !--------------------------------------------------------------------------!

contains

  subroutine allocateFields(Const,Mesh,RB,Uvar,rhs,Exakt)
    
    use TYPES
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tConstants)       :: Const      ! Konstanten                        !
    type(tMesh)            :: Mesh       ! Gitterwerte                       !
    type(tRandbedingungen) :: RB         ! Randbedingungen                   !
    real,pointer           :: Uvar(:,:)  ! Feld mit der numerischen Auswertg !=Tr(spanng)!
    real,pointer           :: rhs(:,:)   ! Feld für die Rechte Seite         !
    type(tExakt)           :: Exakt      ! Exakte Loesung                    !
    !                                                                        !
    ! Local variable declaration                                             !
    !                                                                        !
    intent(in)              :: Const    
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(inout)           :: Mesh,RB,Exakt                                 !
    !------------------------------------------------------------------------!


    allocate(Mesh%x(1:Mesh%nraumx*Mesh%nraumy), &
	 Mesh%y(1:Mesh%nraumx*Mesh%nraumy),  &
	 Mesh%z(1:Mesh%nraumx*Mesh%nraumy),  &
         Mesh%Coefficients(1:6,Mesh%nraumx*Mesh%nraumy), &
         RB%randl(1:Mesh%nraumx), RB%randr(1:Mesh%nraumx), &
         RB%rando(1:Mesh%nraumy), RB%randu(1:Mesh%nraumy), &
         Uvar(1:Mesh%nraumx,1:Mesh%nraumy),                &
         rhs(1:Mesh%nraumx,1:Mesh%nraumy),                 &
         Exakt%loesung(1:Mesh%nraumx,1:Mesh%nraumy),       &
         STAT = allocStat )
     ! TensionTg sigma_xx, sigma_xy=sigma_yx and sigma_yy  !
     ! Tension3 sigma_xz and sigma_yz  (sigma_zz=0)
    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not allocate all variables!'
       STOP
    end if
    
  end subroutine allocateFields


  subroutine allocateFieldsCurv(Mesh,CurvStruct)
    
    use TYPES
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tMesh)            :: Mesh       ! Gitterwerte                       !
    type(tCurv)            :: CurvStruct ! Kurven Koeff                      !
    !                                                                        !
    ! Local variable declaration                                             !
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(in)              :: Mesh                                          !
    intent(inout)           :: CurvStruct                                    !
    !------------------------------------------------------------------------!

    allocate(CurvStruct%chicoeff(8,Mesh%nraumx*Mesh%nraumy), &
         CurvStruct%der_a1xx_x(Mesh%nraumx*Mesh%nraumy), &
         CurvStruct%der_a1xx_y(Mesh%nraumx*Mesh%nraumy), &
         CurvStruct%der_a1xy_x(Mesh%nraumx*Mesh%nraumy), &
         CurvStruct%der_a1xy_y(Mesh%nraumx*Mesh%nraumy), &
         CurvStruct%der_a2yx_x(Mesh%nraumx*Mesh%nraumy), &
         CurvStruct%der_a2yx_y(Mesh%nraumx*Mesh%nraumy), &
         CurvStruct%der_a2yy_x(Mesh%nraumx*Mesh%nraumy), &
         CurvStruct%der_a2yy_y(Mesh%nraumx*Mesh%nraumy), &
         STAT = allocStat )
  
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsCurv: Could not allocate all variables!'
       STOP
    end if

  end subroutine allocateFieldsCurv


  subroutine allocateFieldsVarNum(Mesh,VarNum)
 
  use TYPES

  implicit none
  ! Variable Deklaration  ---------------------------------------------------!
    type(tMesh)            :: Mesh       ! Gitterwerte                       !
    type(tNumeric)         :: VarNum     ! Numerische Felder mit Loesungen   !
    !------------------------------------------------------------------------!
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(inout)           :: VarNum                                        !
    !------------------------------------------------------------------------!

    allocate(VarNum%Verschiebung(1:3,Mesh%nraumx*Mesh%nraumy), &  
         VarNum%TensionTg(1:3,Mesh%nraumx*Mesh%nraumy), &
         VarNum%Tension3(1:2,Mesh%nraumx*Mesh%nraumy), &
         VarNum%BiegMoment(1:4,Mesh%nraumx*Mesh%nraumy), &
         STAT = allocStat )
  
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsVarNum: Could not allocate all variables!'
       STOP
    end if

  end subroutine allocateFieldsVarNum
  

  subroutine allocateGitter2D(Mesh2D)
  
  use TYPES

  implicit none    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    type(tMeshGen)         :: Mesh2D     ! Gitter Punkten und Elements (Quad)!
    ! -----------------------------------------------------------------------!
    integer                :: allocStat

    allocate(Mesh2D%x(1:Mesh2D%nodes), &
	 Mesh2D%y(1:Mesh2D%nodes),  &
	 Mesh2D%z(1:Mesh2D%nodes),  &
         Mesh2D%quad(1:8, 1:Mesh2D%elmts), &
         Mesh2D%quadtop(1:Mesh2D%ntop), & 
         Mesh2D%quadbottom(1:Mesh2D%nbottom), & 
         Mesh2D%quadleft(1:Mesh2D%nleft), & 
         Mesh2D%quadright(1:Mesh2D%nright), & 
         STAT = allocStat )

    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not allocate all variables!'
       STOP
    end if

  end subroutine allocateGitter2D





  subroutine deallocateFields(Const,Mesh,RB,Exakt,Uvar,rhs)

    use TYPES
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tConstants)       :: Const      ! Konstanten                        !
    type(tMesh)            :: Mesh       ! Gitterwerte                       !
    type(tRandbedingungen) :: RB         ! Randbedingungen                   !
    type(tExakt)           :: Exakt      ! Exakte Loesung                    !
    real,pointer           :: Uvar(:,:)  ! Feld mit der numerischen Loesung  !
    real,pointer           :: rhs(:,:)   ! Feld für die Rechte Seite         !
    !                                                                        !
    ! Local variable declaration                                             !
    !                                                                        !
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(in)              :: Const    
    intent(inout)           :: Mesh,RB,Exakt 
    !------------------------------------------------------------------------!
  
    deallocate(Mesh%x, Mesh%y, Mesh%z, Mesh%Coefficients, &
              RB%randl, RB%randr, RB%rando, RB%randu,  &
              Exakt%loesung, Uvar, rhs, &
               STAT=allocStat)

    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not deallocate correctly!'
       STOP
    end if

  end subroutine deallocateFields


  subroutine deallocateFieldsCurv(CurvStruct)

    use TYPES
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tCurv)            :: CurvStruct ! Gitterwerte                       !
    !                                                                        !
    ! Local variable declaration                                             !
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(inout)           :: CurvStruct                                    !
    !------------------------------------------------------------------------!

    deallocate(CurvStruct%chicoeff, &
               CurvStruct%der_a1xx_x, CurvStruct%der_a1xx_y, &
               CurvStruct%der_a1xy_x, CurvStruct%der_a1xy_y, &
               CurvStruct%der_a2yx_x, CurvStruct%der_a2yx_y, &
               CurvStruct%der_a2yy_x, CurvStruct%der_a2yy_y, &
               STAT=allocStat)

    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsCurv: Could not deallocate correctly!'
       STOP
    end if

  end subroutine deallocateFieldsCurv


  subroutine deallocateFieldsVarNum(VarNum)

    use TYPES
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tNumeric)          :: VarNum    ! Gitterwerte                       !
    !                                                                        !
    ! Local variable declaration                                             !
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(inout)           :: VarNum                                    !
    !------------------------------------------------------------------------!
    deallocate(VarNum%Verschiebung, VarNum%TensionTg, &
               VarNum%Tension3, VarNum%BiegMoment, &
               STAT=allocStat)

  end subroutine deallocateFieldsVarNum


  subroutine deallocateGitter2D(Mesh2D)
  
  use TYPES

  implicit none    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    type(tMeshGen)         :: Mesh2D     ! Gitter Punkten und Elements (Quad)!
    ! -----------------------------------------------------------------------!
    integer                :: allocStat

    deallocate(Mesh2D%x, Mesh2D%y, Mesh2D%z, Mesh2D%quad, &
               Mesh2D%quadtop, Mesh2D%quadbottom, Mesh2D%quadleft, Mesh2D%quadright, & 
         STAT = allocStat )

    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not deallocate correctly!'
       STOP
    end if

  end subroutine deallocateGitter2D

end module allocateFields_mod
