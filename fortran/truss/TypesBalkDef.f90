!*****************************************************************************!
!* 
!*      $RCSfile: TypesBalkDef.f90,v $
!*
!*      $Revision: 1.0 $
!*
!*      $Date: 2019/10/9 08:16:09 $
!*
!*      $Author: Fpeugny
!*
!*****************************************************************************!

MODULE TypesBalken
  !---------------------------------------------------------------------------!
  use ISO_C_BINDING
  implicit none
  private
  !---------------------------------------------------------------------------!
  type tMeshInfo
     integer        :: nn,nt          ! Anzahl Nodes and Element
     real           :: EY             ! YOUNG Modulus
     real           :: nu             ! Poisson Coeff
  end type tMeshInfo

  type tMeshCoord
     real,pointer   :: x(:),y(:),z(:)        ! Nodes Coord
     integer,pointer :: elmts(:,:)   ! (1:nelmt; ie c1 d1 f1 c2 d2 f2)
                                     ! 7 numeros per line since index und per Nodes in elements (2)
                                     !  + types number (1 to 7) of connection at node_i
                                     !     + 1= >o--  bar articulated
                                     !     + 2= :>--  bar in a rail
                                     !     + 3= //=   beam fixed
                                     !     + 4= :/=   beam in rail
                                     !     + 5= |o|=  beam articulated
                                     !     + 6= ->o<- bar to bar
                                     !     + 7= =| |= beam to beam
                                     !  + boolean 1=displacement imposed, 0=not imposed (unknown)
                                     !  + boolean 1=force applied known, 0= no force or local reaction 
                                     !         to calculate (depending of types number)
                                     !  + and similar index number for the second node: c2 d2 f2
  end type tMeshCoord

  type tMeshElmt
     integer     :: typ            ! 1=Bar, 2=Balken
     integer     :: nraumx         ! subdivsion, default=2
     real        :: startx         ! node 1 x-koord  !
     real        :: endx           ! node 2 x-koord  !
     real        :: starty         ! node 1 y-koord  !
     real        :: endy           ! node 2 y-koord  !
     real        :: startz         ! node 2 z-koord  !
     real        :: endz           ! node 2 z-koord  !
     real        :: dlen           ! Laengen des Rechengebiets (Balk) in x/y !
     real        :: anglez         ! Angle beim Kurve auf z-Achse in radians
     real        :: SArea          ! Section area
     real        :: EY             ! YOUNG Modulus
     real        :: vu             ! Poisson Coeff
     real        :: CI             ! kinetics inertia moment
     real        :: q1, q2         ! verteilte  Lasten q(x)=x/dlen*q1+(1-x/dlen)*q2
!     real        :: q?
     integer     :: node_1, node_2

     ! real,pointer (c_ptr) :: x(:)        ! do NOT work
     !real,dimension(2)         :: x,y,z

     real,pointer         :: CoeffsH1(:)  ! Coefficients des Hermites Polynoms H1 (H1(0)=1, H1'(0)=0
     real,pointer         :: CoeffsH2(:)  ! Coefficients des Hermites Polynoms H2 (H2(0)=0, H2'(0)=1
     real,pointer         :: CoeffsH3(:)  ! Coefficients des Hermites Polynoms H3 (H3(0)=1, H3'(0)=0
     real,pointer         :: CoeffsH4(:)  ! Coefficients des Hermites Polynoms H4 (H4(0)=0, H4'(0)=1
  end type tMeshElmt

  type tVarElmt
     real,dimension(3)         :: UeIi, UeJi
     real,dimension(3)         :: FeIi,FeJi
  end type tVarElmt

  type tRigidMat
!     real,dimension(3,3)       :: Ke
!      real,dimension(:,:),allocatable  :: Ke
      real,pointer             :: Ke(:,:)
  end type tRigidMat

  type tVarFull
     real,pointer              :: Ue(:)
     real,pointer              :: Fe(:) 
  end type tVarFull

  type tRigidFullMat
     integer                   :: nn,ne
!     real,dimension(:,:),allocatable    :: Ke
     real,pointer              :: Ke(:,:)
  end type tRigidFullMat

  type tExakt
    real,pointer             :: Ux(:)
    real,pointer             :: Fx(:)
  end type tExakt

  type tPolynom
    integer              :: n              ! degree =>  dimension n+1        !
    real,dimension(7)    :: coeffs(1:7)      ! Koeffizienten  order<10         !
  end type tPolynom

!  type tConnectBalk
!    integer,dimension(1:12)  :: tab       ! tab(2),tab(3)=0 tab(5),tab(6)=1 in eine Momentfreie Bewegung
!                                          ! tab(2),tab(3)=0 tab(5) tab(6)=0 in eine Fest eingeklammerte Abbildung
!  end type tConnectBalk

  type tFileIO
     CHARACTER(LEN=60)    :: filename       ! Filename                        !
  end type tFileIO

  !---------------------------------------------------------------------------!
  public  :: tMeshInfo, tMeshCoord, tMeshElmt, tVarElmt, tRigidMat, &
             tRigidFullMat, tVarFull, tExakt, tPolynom, tFileIO
  !---------------------------------------------------------------------------!

! contains

end module TypesBalken
