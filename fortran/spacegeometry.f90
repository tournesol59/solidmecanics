!*****************************************************************************!
!                                                                             !
! module spacegeometry                                                        !
!                                                                             !
!      Funktion: Im Allgemein berechnen, welche Winkel aus einem Ebene        !
!                 , der aus einem anderen Ebene von einem Rotation erzeugt    !
!                 wird,(Winkeln) nötig sind und andere Transformattableion macht   !
!                  Deswegen werden die Datenstrukturen von TypesRoundDef.f90  !
!                 (tNode, tTrigElem, !tPassagemattablerix, benutzt.                !
!                                                                             !
!*****************************************************************************!


module spacegeometry_mod
  use TypesRound
  implicit none

  private
    real      :: TrigRadianToDegree = 57.297

  interface fillmattableX
    module procedure fillmattableX
  end interface

  interface fillmattableY
    module procedure fillmattableY
  end interface

  interface fillmattableZ
    module procedure fillmattableZ
  end interface

  interface trigoprotation
     module procedure trigoprotation
  end interface

  public :: fillmattableX, fillmattableY, fillmattableZ, trigoprotation
  contains

subroutine fillmattableX(mattable,angle)
  implicit none
!  real,pointer  :: mattable(:)
  real,dimension(:),pointer :: mattable
  real          :: angle
  intent(in)    :: angle
  intent(inout) :: mattable
  real                     :: cosval, sinval

  cosval= cos(angle/TrigRadianToDegree)
  sinval= sin(angle/TrigRadianToDegree)

  mattable(1)=1.0
  mattable(2)=0.0
  mattable(3)=0.0
  mattable(4)=0.0
  mattable(5)=cosval
  mattable(6)=-sinval
  mattable(7)=0.0
  mattable(8)=sinval
  mattable(9)=cosval

end subroutine fillmattableX

subroutine fillmattableY(mattable,angle)
  implicit none
!  real,pointer  :: mattable(:)
  real,dimension(:),pointer   :: mattable
  real          :: angle
  intent(in)    :: angle
  intent(inout) :: mattable
  real                     :: cosval, sinval

  cosval= cos(angle/TrigRadianToDegree)
  sinval= sin(angle/TrigRadianToDegree)

  mattable(1)=cosval
  mattable(2)=0.0
  mattable(3)=-sinval
  mattable(4)=0.0
  mattable(5)=1.0
  mattable(6)=0.0
  mattable(7)=sinval
  mattable(8)=0.0
  mattable(9)=cosval

end subroutine fillmattableY

subroutine fillmattableZ(mattable,angle)
  implicit none
!  real,pointer  :: mattable(:)
  real,dimension(:),pointer   :: mattable
  real          :: angle
  intent(in)    :: angle
  intent(inout) :: mattable
  real                     :: cosval, sinval

  cosval= cos(angle/TrigRadianToDegree)
  sinval= sin(angle/TrigRadianToDegree)

  mattable(1)=cosval
  mattable(2)=-sinval
  mattable(3)=0.0
  mattable(4)=sinval
  mattable(5)=cosval
  mattable(6)=0.0
  mattable(7)=0.0
  mattable(8)=0.0
  mattable(9)=1.0

end subroutine fillmattableZ

!********************************************
!* subroutine trigoprotation
!   operate a rotation from one of the canonic axis (iu=1 for x)
!   with an angle "angle" in degree
!   note that the rotated axes are built in an array of nodes (remember to allocate %vects)
!   because of the concision of results (instead of individual elemts) while the input axes
!   are given individually (result from mesh data, individual vertices)
!*******************************************
subroutine trigoprotation(iu,ux,uy,uz,angle,axout_xyz)
   use typesround
   implicit none

! uebergegebene Variabels
    integer                  :: iu
    type(tNode)              :: ux, uy, uz
    real                     :: angle
    type(tPassagematrix)     :: axout_xyz
! Vorhabe
    intent(in)               :: iu, ux, uy, uz, angle
    intent(inout)            :: axout_xyz
! lokale Variable
    integer                  :: i,j, allocStat
    real,dimension(:),pointer   :: mattable

    allocate(mattable(1:9), STAT = allocStat )        
 
    if (allocStat.NE.0) then
       print *, 'ERROR Allocate trigoprotation: Could not allocate mattable(:) !'
       STOP
    end if

    select case (iu)
      case (1) 
         call fillmattableX(mattable,angle)
      case (2)
         call fillmattableY(mattable,angle)
      case (3)
         call fillmattableZ(mattable,angle)
    end select

    axout_xyz%mat(1)= mattable(1)*ux%vects(1) + mattable(4)*uy%vects(1) + mattable(7)*uz%vects(1)
    axout_xyz%mat(2)= mattable(2)*ux%vects(1) + mattable(5)*uy%vects(1) + mattable(8)*uz%vects(1)
    axout_xyz%mat(3)= mattable(3)*ux%vects(1) + mattable(6)*uy%vects(1) + mattable(9)*uz%vects(1)

    axout_xyz%mat(4)= mattable(1)*ux%vects(2) + mattable(4)*uy%vects(2) + mattable(7)*uz%vects(2)
    axout_xyz%mat(5)= mattable(2)*ux%vects(2) + mattable(5)*uy%vects(2) + mattable(8)*uz%vects(2)
    axout_xyz%mat(6)= mattable(3)*ux%vects(2) + mattable(6)*uy%vects(2) + mattable(9)*uz%vects(2)

    axout_xyz%mat(4)= mattable(1)*ux%vects(3) + mattable(4)*uy%vects(3) + mattable(7)*uz%vects(3)
    axout_xyz%mat(5)= mattable(2)*ux%vects(3) + mattable(5)*uy%vects(3) + mattable(8)*uz%vects(3)
    axout_xyz%mat(6)= mattable(3)*ux%vects(3) + mattable(6)*uy%vects(3) + mattable(9)*uz%vects(3)

    deallocate(mattable, STAT=allocStat)
end subroutine trigoprotation

end module spacegeometry_mod
