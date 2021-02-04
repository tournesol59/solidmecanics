!************************************************************
! procedure to fill the -in main prgm allocated tables of struct tMeshElmt
!  and rigidity elements (in the moment redundant with RigidTrussContacts
module CreateTruss_mod
 use typesbalken
 implicit none

 private

  interface createTrussGen
    module procedure createTrussGen
  end interface createTrussGen

  interface matrixfixedcondition
    module procedure matrixfixedcondition
  end interface matrixfixedcondition

#define LINK_BAREND 1
#define LINK_BARSLIDE 2
#define LINK_BEAMEND 3
#define LINK_BEAMSLIDE 4
#define LINK_BEAMFREE 5
#define LINK_BARTOBAR 6
#define LINK_BEAMTOBEAM 7

 public :: createTrussGen, matrixrodcondition, matrixfixedcondition

 contains
!************************************************************
! Fills the (in main prgm allocated table of struct tMeshElmt
 subroutine createTrussGen(nn,nt,MeshGen,MeshPoints,meshconnect)
 use typesbalken
 implicit none
 
 integer                               :: nn,nt
 type(tMeshElmt),pointer               :: MeshGen(:)
 type(tMeshCoord)                      :: MeshPoints
 integer,pointer                       :: meshconnect(:,:)
 intent(in)                            :: nn,nt,MeshPoints,meshconnect
 intent(inout)                         :: MeshGen
! lokal
 integer                               :: i,j,n1,n2,c1,c2,allocStat

  do i=1,nt    ! fuer jedes Element
    n1 = meshconnect(i,2)       ! Index of the node corresp. as one extremity
    n2 = meshconnect(i,3)
    MeshGen(i)%node_1=n1
    MeshGen(i)%node_2=n2

    c1 = MeshPoints%elmts(i,2)  ! type of connection/element type at node_1
    c2 = MeshPoints%elmts(i,5)
    if ((c1.eq.LINK_BAREND).or.(c1.eq.LINK_BARSLIDE).or.(c1.eq.LINK_BARTOBAR)) then
        if ((c2.eq.LINK_BAREND).or.(c2.eq.LINK_BARSLIDE).or.(c2.eq.LINK_BARTOBAR)) then
           MeshGen(i)%typ=1   !bar
        else
           MeshGen(i)%typ=2   !beam
        endif
    elseif ((c1.eq.LINK_BEAMEND).or.(c1.eq.LINK_BEAMSLIDE).or.(c1.eq.LINK_BEAMTOBEAM)) then
           MeshGen(i)%typ=2   !beam
    endif   

    MeshGen(i)%nraumx = 2        ! subdivision, default=2
    MeshGen(i)%startx =  MeshPoints%x(n1)        ! node 1 x-koord  !
    MeshGen(i)%endx =  MeshPoints%x(n2)          ! node 2 x-koord  !
    MeshGen(i)%starty =  MeshPoints%y(n1)        ! node 1 y-koord  !
    MeshGen(i)%endy =  MeshPoints%y(n2)          ! node 2 y-koord  !
    MeshGen(i)%startz =  MeshPoints%z(n1)        ! node 2 z-koord  !
    MeshGen(i)%endz =  MeshPoints%z(n2)          ! node 2 z-koord  !
    MeshGen(i)%dlen = sqrt((  MeshGen(i)%endx -   MeshGen(i)%startx )**2+ &
                    (  MeshGen(i)%endy -   MeshGen(i)%starty )**2+ &
                    (  MeshGen(i)%endz -   MeshGen(i)%startz )**2)           ! Laengen des Rechengebiets (Balk) in x/y !
    MeshGen(i)%anglez = acos( (MeshGen(i)%endx - MeshGen(i)%startx) / MeshGen(i)%dlen)

  allocate(   MeshGen(i)%CoeffsH1(1:4),   MeshGen(i)%CoeffsH2(1:4),   MeshGen(i)%CoeffsH3(1:4), &
               MeshGen(i)%CoeffsH4(1:4), STAT=allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation polynoms elmt 1 !"
  endif  

 ! Fill Hermite  Polynoms, this would be only for implementation of other type of polynoms but will not be used
  MeshGen(i)%CoeffsH1(1)=1.0
  MeshGen(i)%CoeffsH1(2)=0.0
  MeshGen(i)%CoeffsH1(3)=-3.0
  MeshGen(i)%CoeffsH1(4)=2.0

  MeshGen(i)%CoeffsH2(1)=0.0
  MeshGen(i)%CoeffsH2(1)=1.0
  MeshGen(i)%CoeffsH2(1)=-2.0
  MeshGen(i)%CoeffsH2(1)=1.0

  MeshGen(i)%CoeffsH3(1)=0.0
  MeshGen(i)%CoeffsH3(1)=0.0
  MeshGen(i)%CoeffsH3(1)=3.0
  MeshGen(i)%CoeffsH3(1)=-2.0

  MeshGen(i)%CoeffsH4(1)=0.0
  MeshGen(i)%CoeffsH4(1)=0.0
  MeshGen(i)%CoeffsH4(1)=-1.0
  MeshGen(i)%CoeffsH4(1)=1.0

 enddo
 end subroutine createTrussGen


!***********************************************************
! fills a previous allocated rod rigidity matrix by Integration
 subroutine matrixrodcondition(MeshGen, Mat, kpart)
 use typesbalken
 implicit none

 integer                               :: kpart
 type(tMeshElmt)                       :: MeshGen 
 type(tRigidMat)                       :: Mat
 intent(in)                            :: kpart, MeshGen
 intent(inout)                         :: Mat
 ! lokale
 real                                  :: ESsL


  ESsL=MeshGen%EY * MeshGen%SArea / MeshGen%dlen
  if (kpart.eq.1) then
    Mat%Ke(1,1)=ESsL
  else if (kpart.eq.2) then
    Mat%Ke(1,1)=-ESsL
  else if (kpart.eq.3) then
    Mat%Ke(1,1)=ESsL 
  endif
  Mat%Ke(1,2)=0.
  Mat%Ke(1,3)=0.
  Mat%Ke(2,1)=0.
  Mat%Ke(2,2)=0.
  Mat%Ke(2,3)=0.
  Mat%Ke(3,1)=0.
  Mat%Ke(3,2)=0.
  Mat%Ke(3,3)=0.

 end subroutine matrixrodcondition

!***********************************************************
! fills a previous allocated beam rigidity matrix by Integration
 subroutine matrixfixedcondition(MeshGen, Mat, kpart)
 use typesbalken
 implicit none
 
 integer                               :: kpart
 type(tMeshElmt)                       :: MeshGen
 type(tRigidMat)                       :: Mat
 intent(in)                            :: kpart, MeshGen
 intent(inout)                         :: Mat
 ! lokale
  real          :: EIsL, EIsL2, EIsL3, ESsL

  ESsL=MeshGen%EY * MeshGen%SArea / MeshGen%dlen
  EIsL=(MeshGen%EY) * (MeshGen%CI) / (MeshGen%dlen)
  EIsL2=(MeshGen%EY) * (MeshGen%CI) / (MeshGen%dlen) **2
  EIsL3= (MeshGen%EY) * (MeshGen%CI) / (MeshGen%dlen) **3

 if (kpart.eq.1) then
  Mat%Ke(1,1)= ESsL
  Mat%Ke(1,2)=0.0
  Mat%Ke(1,3)=0.0
  Mat%Ke(2,1)=0.0
  Mat%Ke(2,2)=12.0 * EIsL3
  Mat%Ke(2,3)=-6.0 * EIsL2
  Mat%Ke(3,2)= Mat%Ke(2,3)
  Mat%Ke(2,2)= 4.0 * EIsL  
 else if (kpart.eq.2) then
  Mat%Ke(1,1)= -ESsL
  Mat%Ke(1,2)=0.0
  Mat%Ke(1,3)=0.0
  Mat%Ke(2,1)=0.0
  Mat%Ke(2,2)=-12.0 * EIsL3   !***
  Mat%Ke(2,3)=6.0 * EIsL2
  Mat%Ke(3,2)=-6.0 * EIsL2    !***
  Mat%Ke(2,2)= 2.0 * EIsL  
 else if (kpart.eq.3) then
  Mat%Ke(1,1)= ESsL
  Mat%Ke(1,2)=0.0
  Mat%Ke(1,3)=0.0
  Mat%Ke(2,1)=0.0
  Mat%Ke(2,2)=12.0 * EIsL3
  Mat%Ke(2,3)=-6.0 * EIsL2
  Mat%Ke(3,2)= Mat%Ke(2,3)
  Mat%Ke(2,2)= 4.0 * EIsL  
 endif

 end subroutine matrixfixedcondition

end module CreateTruss_mod
