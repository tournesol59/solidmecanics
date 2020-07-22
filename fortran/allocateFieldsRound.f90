!****************************************************************************!
!                                                                            !
! subroutine allocateFields                                                  !
!                                                                            !
!      Funktion: Bereitstellen des Speicherplatzes                           !
!                                                                            !
!****************************************************************************!
!      Z.Bsp                                                                 !
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
  ! siehe auch Spezifikation in TypesRoundDef.f90
  ! Diese Parameter sind nur so als <Parameter-NO> als Variabel nur einmal geschrieben in 
  ! allocateForXXXType() sub, so dass sie als Information fuer die Algorithmus
  ! koennen aber nicht modifiziert werden.

  ! Dieser Parameter ist fuer Node gedacht
  !  integer, parameter :: NumberofDisp = 6  ! mit parameter muss es nicht mofiziert werden

  integer  :: NumberofDisp = 6       ! 6 dof per nn, so kann es geaendert werden !
  ! Diese zwei folg. Parameter sind fuer Element gedacht (in Post treatment)
  integer  :: NumberofTension = 6    ! 6 tang. Kraft per ne, so kann es geaendert werden !
  integer  :: NumberofMom = 6        ! 6 bieg. Momentums per ne, so kann es geaendert werden !

  !--------------------------------------------------------------------------!
  interface allocateForC0BelyuType
     module procedure allocateForC0BelyuType
  end interface
  interface allocateForDKTType
     module procedure allocateForDKTType
  end interface
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
!  interface allocateFieldsCoeffRound
!     module procedure allocateFieldsCoeffRound
!  end interface
  interface allocateSparseRound
     module procedure allocateSparseRound
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
!  interface deallocateFieldsCoeffRound
!     module procedure deallocateFieldsCoeffRound
!  end interface
  interface deallocateSparseRound
     module procedure deallocateSparseRound
  end interface

  !--------------------------------------------------------------------------!
  public  :: allocateForC0BelyuType, allocateForDKTType, &
             allocateFieldsRound, allocateFieldsVarRound, allocateFieldsBCSRound, & ! allocateFieldsCoeffRound,&
             allocateSparseRound, &
             deallocateFieldsRound, deallocateFieldsVarRound, deallocateFieldsBCSRound, & ! deallocateFieldsCoeffRound,
             deallocateSparseRound
  !--------------------------------------------------------------------------!

contains

!******************************************************
! Configuration for Triangle C0 (or Belyushko Element)
  subroutine allocateForC0BelyuType()
    use TypesRound
    implicit none
    !------------------------------------------------------------------------!
    ! no input variables
   NumberofDisp = 15     ! <u1x,u1y,u1z,r1x,r1y,u2x,u2y,u2z,r2x,r2y,u3x,u3y,u3z,r3x,r3y>
   NumberofTension = 6   ! <F1x,F1y,F2x,F2y,F3x,F3y>
   NumberofMom = 6       ! <M1xx,M1yy,M1xy,M2xx,M2yy,M2xy,M3xx,M3yy,M3xy>
  
  end subroutine allocateForC0BelyuType

!******************************************************
! Configuration for Triangle C0 (or Belyushko Element)
  subroutine allocateForDKTType()
    use TypesRound
    implicit none
    !------------------------------------------------------------------------!
    ! no input variables
   NumberofDisp = 15     ! <u1x,u1y,u1z,r1x,r1y,u2x,u2y,u2z,r2x,r2y,u3x,u3y,u3z,r3x,r3y>
   NumberofTension = 6   ! <F1x,F1y,F2x,F2y,F3x,F3y>
   NumberofMom = 6       ! <M1xx,M1yy,M1xy,M2xx,M2yy,M2xy,M3xx,M3yy,M3xy>
  
  end subroutine allocateForDKTType

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
    integer                 :: nn,ne,i,allocStat,allocStatnn,allocStatne     !
    !------------------------------------------------------------------------!
    intent(inout)           :: MeshR,ExaktR                                  !
    intent(in)              :: NTestMeshR                                    !
    !------------------------------------------------------------------------!

    nn = NTestMeshR%nNode
    ne = NTestMeshR%nElem

    allocate(MeshR%nodes(1:nn), MeshR%elems(1:ne), MeshR%neighbours(1:ne), &            
         ExaktR%loesung(1:nn,1:NumberofDisp),       &
         STAT = allocStat )        
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsRound: Could not allocate all variables!'
       STOP
    end if
! Below may not be used with simple const dimensions, or declare allocatable
    do i=1,nn
       allocate(MeshR%nodes(i)%vects(1:3), STAT=allocStatnn)
       if (allocStatnn.NE.0) then
          print *, 'ERROR AllocateFieldsRound: Could not allocate all component of node!'
          STOP
       end if
    enddo 
    do i=1,ne
       allocate(MeshR%elems(i)%numeros(1:4), STAT=allocStatne)
       allocate(MeshR%neighbours(i)%numeros(1:4), STAT=allocStatne)
       if (allocStatne.NE.0) then
          print *, 'ERROR AllocateFieldsRound: Could not allocate all component of elem!'
          STOP
       end if
   enddo 
  
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
    type(tRoundMeshInfo)    :: NTestMeshR  ! Dimensions                      !
    type(tRoundRandU)               :: BCU       ! Randbedingungen Biegg     !
    type(tRoundRandS)               :: BCS       ! Randbedingungen Kraft     !
    ! Local variable declaration                                             !
    !                                                                        !
    integer                 :: nxmx,nymx,i,allocStat                         !
    !------------------------------------------------------------------------!
    intent(inout)           ::  BCU,BCS                                      !
    intent(in)              :: NTestMeshR                                    !
    !------------------------------------------------------------------------!

    nxmx = NTestMeshR%nxmx
    nymx = NTestMeshR%nymx
! why shall be only nxmx taken care?
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
                          ! gut aber nicht je nach Element sondern nach Node und Tz ist dabei: 6
                          ! fuer tension und biegmoment je nach Element, da es im Post Prozessing ist.
    allocate(VarNumR%displacement(1:nn,1:NumberofDisp),       &
             VarNumR%rhs(1:nn,1:NumberofDisp),       &
             VarNumR%tension(1:ne,1:NumberofTension),        &
             VarNumR%biegmoment(1:ne,1:NumberofMom),        &
         STAT = allocStat )        
 
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFieldsVarRound: Could not allocate all variables!'
       STOP
    end if

  end subroutine allocateFieldsVarRound

!******************************************************
!  subroutine allocateFieldsCoeffRound(NTestMeshR,MeshCoeffR) 
! bind (C, name="allocateFieldsCoeffRound")
!    use ISO_C_BINDING, only: c_char 
!  use TypesRound

!  implicit none
!  ! Variable Deklaration  ---------------------------------------------------!
!    type(tRoundMeshInfo)   :: NTestMeshR  ! Dimensions                       !
!    type(tRoundCoeff)      :: MeshCoeffR     ! Numerische Felder mit Loesungen  !
!    !------------------------------------------------------------------------!
!    integer                 :: allocStat,ne,nn                               !
!    !------------------------------------------------------------------------!
!    intent(inout)           :: MeshCoeffR                                       !
!    intent(in)              :: NTestMeshR                                    !
!    !------------------------------------------------------------------------!!!
!
!    nn = NTestMeshR%nNode
!    ne = NTestMeshR%nElem ! A priori interpol order 2 => 6 coeffs per Elements!
!
!    allocate(MeshCoeffR%Coefficients(1:ne*6),       &
!         STAT = allocStat )        
! 
!    if (allocStat.NE.0) then
!       print *, 'ERROR AllocateFieldsCoeffRound: Could not allocate all variables!'
!       STOP
!    end if!
!
!  end subroutine allocateFieldsCoeffRound

!*****************************************************
  subroutine allocateSparseRound(NTestMeshR, RigidYMat)
  use TypesRound

  implicit none
  
    !------------------------------------------------------------------------!
    ! Variabledeklarationen
    !------------------------------------------------------------------------!
    ! Liste der uebergebenen Argumente                                       !
    !                                                                        !
    type(tRoundMeshInfo)   :: NTestMeshR
    type(tSparseYMat)      :: RigidYMat
    !                                                                        !
    ! Local variable declaration                                             !
    integer                 :: nn,ne,i,allocStat                             !
    !------------------------------------------------------------------------!
    intent(inout)           :: RigidYMat                                     !
    intent(in)              :: NTestMeshR                                    !
    !------------------------------------------------------------------------!

    nn = NTestMeshR%nNode
    ! Die Spezifikation der Sparse Matrix ist auf Element Zahl basiert. (Kompliziert aber eine 
    ! dense Matrix waere von Dimension (ne*ne)=(deltax)**2*(deltay)**2 fuer quadrangles, auch wenn wir
    ! hier Triangles nutzen.
    allocate(RigidYMat%a(1:ne*6*15), RigidYMat%ia(1:ne*6*15), RigidYMat%ja(1:ne*6*15),      &
         STAT = allocStat )        
 
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateSparseRound: Could not allocate all variables!'
       STOP
    end if

  end subroutine allocateSparseRound



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
    integer                 :: nn,ne,i,allocStat,allocStatnn,allocStatne     !
    !------------------------------------------------------------------------!
    intent(inout)           :: MeshR,ExaktR                                  !
    intent(in)              :: NTestMeshR                                    !
    !------------------------------------------------------------------------!

    nn = NTestMeshR%nNode
    ne = NTestMeshR%nElem
! Below cannot be not used if simple const dimensions
    do i=1,nn
       deallocate(MeshR%nodes(i)%vects,STAT=allocStatnn)
       if (allocStatnn.NE.0) then
          print *, 'ERROR DellocateFieldsRound: Could not deallocate node correctly!'
          STOP
       endif
    enddo 
    do i=1,ne
       deallocate(MeshR%elems(i)%numeros,STAT=allocStatne)
       deallocate(MeshR%neighbours(i)%numeros,STAT=allocStatne)
       if (allocStatne.NE.0) then
          print *, 'ERROR DellocateFieldsRound: Could not deallocate elem correctly!'
          STOP
       endif
   enddo 
    deallocate(MeshR%nodes, MeshR%elems, MeshR%neighbours, &
              ExaktR%loesung, &
               STAT=allocStat)
    
    if (allocStat.NE.0) then
       print *, 'ERROR DellocateFieldsRound: Could not deallocate correctly!'
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
       print *, 'ERROR DeallocateFieldsBCSRound: Could not deallocate all variables!'
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
       print *, 'ERROR DeallocateFieldsVarRound: Could not deallocate correctly!'
       STOP
    else
       print *, 'Terminated AllocateFieldsVARRound: OK'
    end if

  end subroutine deallocateFieldsVarRound

!******************************************************
!  subroutine deallocateFieldsCoeffRound(MeshCoeffR) 
! bind (C, name="deallocateFieldCoeffRound")
!    use ISO_C_BINDING, only: c_char  

!    use TypesRound
    
!    implicit none

!  ! Variable Deklaration  ---------------------------------------------------!
!    type(tRoundCoeff)       :: MeshCoeffR     ! Numerische Felder Interpol Coeff  !
!    !------------------------------------------------------------------------!
!    integer                 :: allocStat                                     !
!    !------------------------------------------------------------------------!
!    intent(inout)           :: MeshCoeffR                                    !
!    !------------------------------------------------------------------------!!
!
!    deallocate(MeshCoeffR%Coefficients, &
!               STAT=allocStat)
!
!    if (allocStat.NE.0) then
!       print *, 'ERROR DeallocateFieldsCoeffRound: Could not deallocate correctly!'
!       STOP
!    end if
!
!  end subroutine deallocateFieldsCoeffRound

!*****************************************************
  subroutine deallocateSparseRound(RigidYMat)
  use TypesRound

  implicit none
  
    !-----------------------------------------------------------------------!
    ! Variabledeklarationen
    !-----------------------------------------------------------------------!
    ! Liste der uebergebenen Argumente                                      !
    !                                                                       !
    type(tSparseYMat)      :: RigidYMat
    !                                                                        !
    intent(inout)           :: RigidYMat                                  !
    !                                                                        !
    ! Local variable declaration                                             !
    integer                 :: nn,allocStat
    !------------------------------------------------------------------------!

    deallocate(RigidYMat%a, RigidYMat%ia, RigidYMat%ja,      &
         STAT = allocStat )        
 
    if (allocStat.NE.0) then
       print *, 'ERROR DeallocateSparseRound: Could not deallocate all variables!'
       STOP
    end if

  end subroutine deallocateSparseRound

end module allocateFieldsRound_mod
