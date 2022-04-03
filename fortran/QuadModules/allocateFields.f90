module allocateFields_mod

  use TYPES

  implicit none

  private

  interface allocateFields
     module procedure allocateFields
  end interface

  interface allocateFieldsPatch5el
     module procedure allocateFieldsPatch5el
  end interface

  interface allocateFieldsCurv
     module procedure allocateFieldsCurv
  end interface

  interface allocateFieldsVarNum
     module procedure allocateFieldsVarNum
  end interface

  interface deallocateFields
     module procedure deallocateFields
  end interface

  interface deallocateFieldsPatch5el
     module procedure deallocateFieldsPatch5el
  end interface

  interface deallocateFieldsCurv
     module procedure deallocateFieldsCurv
  end interface

  interface deallocateFieldsVarNum
     module procedure deallocateFieldsVarNum
  end interface

  public :: allocateFields, allocateFieldsPatch5el, allocateFieldsCurv, allocateFieldsVarNum, &
             deallocateFields, deallocateFieldsPatch5el, deallocateFieldsCurv, deallocateFieldsVarNum

  contains

  subroutine allocateFields(NTestMesh,Mesh2D,RB,Uvar,rhs,Exakt,coeffic)
    
    use TYPES
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tMeshInfo)        :: NTestMesh  
!    type(tConstants)       :: Const      ! Konstanten                        !
    type(tMesh)            :: Mesh2D       ! Gitterwerte                       !
    type(tRandbedingungen) :: RB         ! Randbedingungen                   !
    real,pointer           :: Uvar(:,:)  ! Feld mit der numerischen Auswertg !=Tr(spanng)!
    real,pointer           :: rhs(:,:)   ! Feld für die Rechte Seite         !
    type(tExakt)           :: Exakt      ! Exakte Loesung                    !
    type(tCoeff)           :: coeffic    ! polynom coefficient if useful
    !                                                                        !
    ! Local variable declaration                                             !
    !                                                                        !
    intent(in)              :: NTestMesh   
    integer                 :: allocStat,nx,ny                         !
    !------------------------------------------------------------------------!
    intent(inout)           :: Mesh2D,RB,Exakt                                 !
    !------------------------------------------------------------------------!
    nx=NTestMesh%nraumx
    ny=NTestMesh%nraumy

    allocate(Mesh2D%x(1:nx*ny), &
         Mesh2D%y(1:nx*ny),  &
         Mesh2D%z(1:nx*ny),  &
         coeffic%Coefficients(1:6,nx,ny), &
         RB%randl(1:nx), RB%randr(1:nx), &
         RB%rando(1:ny), RB%randu(1:ny), &
         Uvar(1:nx,1:ny),                &
         rhs(1:nx,1:ny),                 &
         Exakt%loesung(1:nx,1:ny),       &
         STAT = allocStat )
     ! TensionTg sigma_xx, sigma_xy=sigma_yx and sigma_yy  !
     ! Tension3 sigma_xz and sigma_yz  (sigma_zz=0)
    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not allocate all variables!'
       STOP
    end if
    
  end subroutine allocateFields

 
  subroutine allocateFieldsPatch5el(Const,Mesh2D,RB,Uvar,rhs,Exakt)
    
    use TYPES
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tConstants)       :: Const      ! Konstanten                        !
  !  type(tMeshInfo)        :: NTestMesh  ! Dimensions Info structure 
    type(tMeshGen)         :: Mesh2D       ! Gitterwerte                       !
    type(tRandbedingungen) :: RB         ! Randbedingungen                   !
    real,pointer           :: Uvar(:,:)  ! Feld mit der numerischen Auswertg !=Tr(spanng)!
    real,pointer           :: rhs(:,:)   ! Feld für die Rechte Seite         !
    type(tExakt)           :: Exakt      ! Exakte Loesung                    !
    !                                                                        !
    ! Local variable declaration                                             !
    integer                :: nn,ne
    !                                                                        !
    intent(in)              :: Const    
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(inout)           :: Mesh2D,RB,Exakt                                 !
    !------------------------------------------------------------------------!
    ! FOR THIS ALLOCATION ALL STRUCTS ARE FIXED DIMENSION SINCE MANUAL PURPOSE
    nn=8
    ne=5
    
    allocate(Mesh2D%x(1:8), &
         Mesh2D%y(1:8),  &
         Mesh2D%z(1:8),  &
         Mesh2D%quad(1:5,1:4), &
         Mesh2D%neighbours(1:5,1:4), &
         Mesh2D%quadtop(1:1), &
         Mesh2D%quadbottom(1:1), &
         Mesh2D%quadleft(1:1), &
         Mesh2D%quadright(1:1), &
         RB%randl(1:2), RB%randr(1:2), &
         RB%rando(1:2), RB%randu(1:2), &
         Uvar(1,1:2*8),                &
         rhs(1,1:2*8),                 &
         Exakt%loesung(1,1:8),       &
         STAT = allocStat )
     ! TensionTg sigma_xx, sigma_xy=sigma_yx and sigma_yy  !
     ! Tension3 sigma_xz and sigma_yz  (sigma_zz=0)
    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not allocate all variables!'
       STOP
    end if
    
  end subroutine allocateFieldsPatch5el

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
  
  subroutine deallocateFields(Const,Mesh,RB,Exakt,Uvar,rhs,coeffic)

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
    type(tCoeff)           :: coeffic
    !                                                                        !
    ! Local variable declaration                                             !
    !                                                                        !
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(in)              :: Const    
    intent(inout)           :: Mesh,RB,Exakt,coeffic 
    !------------------------------------------------------------------------!
  
    deallocate(Mesh%x, Mesh%y, Mesh%z, coeffic%Coefficients, &
              RB%randl, RB%randr, RB%rando, RB%randu,  &
              Exakt%loesung, Uvar, rhs, &
               STAT=allocStat)

    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not deallocate correctly!'
       STOP
    end if

  end subroutine deallocateFields

 subroutine deallocateFieldsPatch5el(Const,Mesh2D,RB,Exakt,Uvar,rhs)

    use TYPES
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tConstants)       :: Const      ! Konstanten                        !
    type(tMeshGen)         :: Mesh2D       ! Gitterwerte                       !
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
    intent(inout)           :: Mesh2D,RB,Exakt 
    !------------------------------------------------------------------------!
  
    deallocate(Mesh2D%x, Mesh2D%y, Mesh2D%z, Mesh2D%quad, Mesh2D%neighbours, & 
              Mesh2D%quadtop, Mesh2D%quadbottom, Mesh2D%quadleft, Mesh2D%quadright, &
              RB%randl, RB%randr, RB%rando, RB%randu,  &
              Exakt%loesung, Uvar, rhs, &

               STAT=allocStat)

    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not deallocate correctly!'
       STOP
    end if

  end subroutine deallocateFieldsPatch5el


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

end module allocateFields_mod
