module allocateTruss_mod
 use typesbalken
 implicit none

  private

  interface allocateTrussDef
    module procedure allocateTrussDef
  end interface allocateTrussDef

  interface deallocateTrussDef
    module procedure deallocateTrussDef
  end interface deallocateTrussDef

  interface allocateGridElmt
    module procedure allocateGridElmt
  end interface allocateGridElmt

  interface deallocateGridElmt
    module procedure deallocateGridElmt
  end interface deallocateGridElmt

  public :: allocateTrussDef, deallocateTrussDef, &
            allocateGridElmt, deallocateGridElmt !, & 
        ! deallocateRigidMat, allocateRigidMat

  contains

!*******************************************************
! procedure to allocate memory of problem definition (further see
! InputTruss.f90
  subroutine allocateTrussDef(MeshInfo,Dunknowns,Fmovemt, &
             Dimposed,Fimposed)
  use typesbalken

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!


  type(tMeshInfo)            :: MeshInfo     ! GitterDims       
   ! integer,pointer           :: meshconnect(:,:) ! Definition von Verbindungen !
  integer,pointer           :: Dunknowns(:,:)
  real,pointer               :: Dimposed(:,:)
  integer,pointer           :: Fmovemt(:,:)
  real,pointer               :: Fimposed(:,:)
  intent(inout)              :: Dunknowns, Fmovemt,Dimposed,Fimposed
  integer                    :: infoStat,nn,nt

   nn=MeshInfo%nn

   allocate(Dunknowns(1:nn,1:3), & 
            Dimposed(1:nn,1:3), Fmovemt(1:nn,1:3), Fimposed(1:nn,1:3), &
            STAT=infoStat)
    if (infoStat.ne.0) then
       print *,"Error allocate Truss Definition"
    endif

  end subroutine allocateTrussDef


     !**********
     ! deallocate

  subroutine deallocateTrussDef(Dunknowns,Fmovemt,Dimposed,Fimposed)
  use typesbalken

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!

  integer,pointer           :: Dunknowns(:,:)
  real,pointer               :: Dimposed(:,:)
  integer,pointer           :: Fmovemt(:,:)
  real,pointer               :: Fimposed(:,:)
  integer                    :: infoStat,nn,nt

  deallocate(Dunknowns, Fmovemt, Dimposed, Fimposed, & 
             STAT=infoStat)
  if (infoStat.ne.0) then
       print *,"Error deallocate Truss Definition"
  endif

  end subroutine deallocateTrussDef

!*******************************************
!  Grid Element MeshPoints+MeshGen
  subroutine allocateGridElmt(MeshInfo,MeshPoints,meshpattern,MeshGen)
  use typesbalken

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  type(tMeshInfo)       :: MeshInfo
  type(tMeshElmt),pointer :: MeshGen(:)   ! fuer nt Elementen 1-2 und 1-3 THIS WORKS
  type(tMeshCoord)      :: MeshPoints
  integer,pointer       :: meshpattern(:,:)
  integer               :: allocStat,nn,nt

   nn=MeshInfo%nn
   nt=MeshInfo%nt
   allocate(MeshPoints%x(1:nn), MeshPoints%y(1:nn), &
            MeshPoints%z(1:nn), &
            MeshPoints%elmts(1:nt,1:7), &
            STAT=allocStat)
   if (allocstat.ne.0) then
     print *,"allocateTruss: error allocate MeshPunkte"
   endif

   allocate(meshpattern(1:nt,1:3), STAT = allocStat)
   if (allocStat.ne.0) then
      print *," Error allocation connection matrix !"
   endif

  allocate(MeshGen(1:nt), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"allocateTruss: error allocate MeshGen null"
  endif

  end subroutine allocateGridElmt
  
  !***********************************************
  ! deallocate
  subroutine deallocateGridElmt(MeshPoints,meshpattern,MeshGen)
  use typesbalken

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  type(tMeshInfo)       :: MeshInfo
  type(tMeshElmt),pointer :: MeshGen(:)   ! fuer nt Elementen 1-2 und 1-3 THIS WORKS
  type(tMeshCoord)      :: MeshPoints
  integer,pointer       :: meshpattern(:,:)
  integer               :: allocStat

   deallocate(MeshPoints%x, MeshPoints%y, MeshPoints%z, &
            MeshPoints%elmts, &
            STAT=allocStat)
   if (allocstat.ne.0) then
     print *,"deallocateTruss: error deallocate MeshPunkte"
   endif

   deallocate(meshpattern, STAT = allocStat)
   if (allocStat.ne.0) then
      print *," Error deallocation connection matrix !"
   endif

  deallocate(MeshGen, STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"deallocateTruss: error allocate MeshGen null"
  endif

  end subroutine deallocateGridElmt


end module allocateTruss_mod
