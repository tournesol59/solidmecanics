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
!        Verschiebung(1:3,:) ! 3 Verschiebungen in der Mitte der Platte        !
!        TensionTg(1:4,:)    ! array pointer: Spannung tengential              !
!        TensionT3(1,:)    ! array pointer: Spannung normal despite that E,xz=0 !
!        BiegMoment(1:3,:)   !array pointer: contains as last arrays above    !
!                           the inital value of the problem and then re-calculated.!

! Before testing further called from C Fortran allocation, try a simplified array
!****************************************************************************!

module allocateFields_mod
  use, INTRINSIC :: ISO_C_BINDING
  use types

  !--------------------------------------------------------------------------!
  implicit none

!   interface
!    subroutine allocateFields(NTestMesh,Mesh,RB,Uvar,rhs,Exakt) bind name="allocateFields")
!       end subroutine allocateFields
!  end interface

!   interface
!    subroutine deallocateFields(Mesh,RB,Exakt,Uvar,rhs) bind (C, name="deallocateFields")
!    end subroutine deallocateFields
!  end interface

  private
  !--------------------------------------------------------------------------!
  interface allocateFields
     module procedure allocateFields
  end interface

  interface deallocateFields
     module procedure deallocateFields
  end interface

!  interface demo_accessFields
!     module procedure demo_accessFields
!  end interface

  !--------------------------------------------------------------------------!
  public  :: allocateFields, deallocateFields
 !--------------------------------------------------------------------------!

contains

  subroutine allocateFields(NTestMesh_,Mesh_,RB_,Uvar_,rhs_,Exakt_) bind (C, name="allocateFields")
    use, INTRINSIC :: ISO_C_BINDING
    use TYPES    
    implicit none
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente u. lokale Vars                         !
    !                                                                        !
    type(tMeshInfo),bind(C)  :: NTestMesh  ! Ref: is allowed because tMeshInfo has 2 int
!    type(c_ptr)            :: NTestMesh_  ! pointer to convert by callee !
    type(tMesh)              :: Mesh       ! Ref: bind(C) is not allowed bec. 
                                           !tMesh has pointers to allocate
    type(c_ptr)            :: Mesh_      ! .. !
    type(tRandbedingungen),bind(C)  :: RB
!    type(c_ptr)            :: RB_        ! .. !
    real(c_float),pointer           :: Uvar
!    type(c_ptr)            :: Uvar_ 
    real(c_float),pointer           :: rhs
!    type(c_ptr)            :: rhs_
    type(tExakt),bind(C)     :: Exakt
!    type(c_ptr)            :: Exakt_
    !                                                                        !
    intent(in)              :: NTestMesh  
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(inout)           :: Mesh,RB,Exakt,Uvar,rhs,Exakt                  !
    !------------------------------------------------------------------------!                                                                      !
    ! Local variable declaration                                             !    



    allocate(Mesh%x(1:NTestMesh%nraumx*NTestMesh%nraumy), &
	 Mesh%y(1:NTestMesh%nraumx*NTestMesh%nraumy),  &
	 Mesh%z(1:NTestMesh%nraumx*NTestMesh%nraumy),  &
         Mesh%Coefficients(1:6,NTestMesh%nraumx*NTestMesh%nraumy), &
         RB%randl(1:NTestMesh%nraumx), RB%randr(1:NTestMesh%nraumx), &
         RB%rando(1:NTestMesh%nraumy), RB%randu(1:NTestMesh%nraumy), &
         Uvar(1:NTestMesh%nraumx,1:NTestMesh%nraumy),                &
         rhs(1:NTestMesh%nraumx,1:NTestMesh%nraumy),                 &
         Exakt%loesung(1:NTestMesh%nraumx,1:NTestMesh%nraumy),       &
         STAT = allocStat )
    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not allocate all variables!'
       STOP
    end if
    
  end subroutine allocateFields


  subroutine deallocateFields(Mesh,RB,Exakt,Uvar,rhs) bind (C, name="deallocateFields")
    use, INTRINSIC :: ISO_C_BINDING
    use TYPES
    
    implicit none
    
    !------------------------------------------------------------------------!
    ! Variablendeklarationen                                                 !
    !------------------------------------------------------------------------!
    ! Liste der übergebenen Argumente                                        !
    !                                                                        !
    type(tMeshInfo),bind(C)  :: NTestMesh
!    type(c_ptr)            :: NTestMesh_  ! pointer to convert by callee !
    type(tMesh),bind(C)      :: Mesh
!    type(c_ptr)            :: Mesh_      ! .. !
    type(tRandbedingungen),bind(C)  :: RB
!    type(c_ptr)            :: RB_        ! .. !
     real(c_float),pointer  :: Uvar
!    type(c_ptr)            :: Uvar_ 
    real(c_float),pointer   :: rhs
!    type(c_ptr)            :: rhs_
    type(tExakt),bind(C)     :: Exakt
 !   type(c_ptr)            :: Exakt_
    !
    ! lokale Variabledeklaration:
    integer                 :: allocStat                                     !
    !------------------------------------------------------------------------!
    intent(in)              :: NTestMesh  
    intent(inout)           :: Mesh,RB,Exakt,Uvar,rhs                        !
    !------------------------------------------------------------------------!
!    NTestMesh = c_loc(NTestMesh_)
!    Mesh = c_loc(Mesh_)
!    RB = c_loc(RB_)
!    Uvar = c_loc(Uvar_)
!    rhs = c_loc(rhs_)
!    Exakt = c_loc(Exakt_)

    deallocate(Mesh%x, Mesh%y, Mesh%z,        &
         STAT = allocStat )

    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not deallocate correctly!'
       STOP
    end if

  end subroutine deallocateFields

end module allocateFields_mod
