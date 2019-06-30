!**************************************************************************!
!                                                                          !
! subroutine createMesh                                                    !
!                                                                          !
!      Funktion: Festlegen der Gitterdaten                                 !
!                                                                          !
!**************************************************************************!
!                                                                          !
!      x,y        Koordinaten der Gitterpunkte                             !
!      dx,  dy    Gitterschrittweite (delta x, delta y)                    !
!      dxq, dyq   Kehrwert (1/delta x, 1/delta y)                          !
!      dxx, dyy   Quadrat  (delta x * delta x, delta y * delta y)          !
!      dxxq,dyyq  Kehrwerte der Quadrate (1/dxx, 1/dyy)                    !
!                                                                          !
!**************************************************************************!
module createmeshround_mod
  use TypesRound

  implicit none

  private
!********************************************************************!
!  Variable Tabelle der Test-Plate fur die Koordinaten der 9 Punkten !
!********************************************************************!
 integer,dimension(2) :: NTest2 = (/ 3, 3/)
 real,dimension(9)  ::  X2 = (/ 0.0, 0.05, 0.1, 0.0, 0.05, 0.1, 0.0, 0.05, 0.1 /)

 real,dimension(9)  ::  Y2 = (/ 0.0, 0.0, 0.0, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1 /)

 real,dimension(9)  ::  Z2 = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

  public :: createRoundMesh2, createRoundMeshFile, &
            createRoundRand2, createRoundRandFile

  contains

!******************************************************
  subroutine createRoundMesh2(NTestMeshR,MeshR)

  use TypesRound

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tRoundMeshInfo)            :: NTestMeshR  ! ne,nn                     !
  type(tRoundMesh)                :: MeshR       ! Gitterwerte               !
  !                                                                          !
  ! Local variable declaration                                               !
  !                                                                          !
  integer                    :: ne,nn,i,j,k,l        ! Zählvariablen             !
  !--------------------------------------------------------------------------!
  intent(inout)              :: NTestMeshR, MeshR                            !
  !--------------------------------------------------------------------------!


  !-----------------------------------<  calculate constants  >-----------
  !---<  nn   = anzahl der inneren punkte                 >-----------
  !---<  ne = anzahl der zellen                           >-----------

  NTestMeshR%nNode = NTest2(1)*NTest2(2)            ! 9
  NTestMeshR%nElem = (NTest2(2)-1)*2*(NTest2(1)-1)  ! 8
  nn = NTestMeshR%nNode
  ne = NTestMeshR%nElem

  !---<  node(1..3..7)             = linker  rand             >-----------
  !---<  node(5..6..9)             = rechter rand             >-----------
  !---<  node(1..2..5)             = unten rand               >-----------
  !---<  node(7..8..9)             = obenn rand               >-----------
  do i= 1,nn
    MeshR%nodes(i)%vects(1) = X2(i)
    MeshR%nodes(i)%vects(2) = Y2(i)
    MeshR%nodes(i)%vects(3) = Z2(i)
  enddo 
  do i= 1,(NTest2(1)-1) ! no overshoot of indexes
    do j=1,(NTest2(2)-1)*2
      k=(i-1)*NTest2(2) + j  ! 1..3 4..6 7..9 Linie von 3 Punkten. Sieht X1,Y1,Z1  !
      l=(i-1)*(NTest2(2)-1)*2+j  ! 1..4 5 ..8 Linie von 4 Elements Triangle
      MeshR%elems(l)%numeros(1) = l
      if ((j).eq.1) then
        MeshR%elems(l)%numeros(2) = k
        MeshR%elems(l)%numeros(3) = k+1
        MeshR%elems(l)%numeros(4) = k+3
        MeshR%neighbours(l)%numeros(1)=-1
        MeshR%neighbours(l)%numeros(2)=-1
        MeshR%neighbours(l)%numeros(3)=l+1
      else if ((j).eq.2) then
        MeshR%elems(l)%numeros(2) = k+4
        MeshR%elems(l)%numeros(3) = k+3
        MeshR%elems(l)%numeros(4) = k+1
        MeshR%neighbours(l)%numeros(1)=l+1
        MeshR%neighbours(l)%numeros(2)=l-1+(NTest2(2)-1)*2
        MeshR%neighbours(l)%numeros(1)=l-1
      else if ((j).eq.3) then
        MeshR%elems(l)%numeros(2) = k+1
        MeshR%elems(l)%numeros(3) = k+2   
        MeshR%elems(l)%numeros(4) = k+4
        MeshR%neighbours(l)%numeros(1)=l-1
        MeshR%neighbours(l)%numeros(2)=-1
        MeshR%neighbours(l)%numeros(3)=l+1   
      else if ((j).eq.4) then
        MeshR%elems(l)%numeros(2) = k+5      
        MeshR%elems(l)%numeros(3) = k+4   
        MeshR%elems(l)%numeros(4) = k+2   
        MeshR%neighbours(l)%numeros(1)=-1
        MeshR%neighbours(l)%numeros(2)=l-1+(NTest2(2)-1)*2
        MeshR%neighbours(l)%numeros(3)=l-1
     endif
    enddo
  enddo 

   write(*, 105) "--------- Koordinaten von Mesh: --------- "
   do i= 1,(NTest2(1)-1) 
    write(*, 207) "y: ", i, MeshR%nodes(i*NTest2(2))%vects(2)
    do j=1,(NTest2(2))
      k=(i-1)*NTest2(2) + j  ! 1..3 4..6 7..9 Linie von 3 Punkten. Sieht X1,Y1,Z1  !
      l=(i-1)*(NTest2(2)-1)*2+j
      write(*, 208) " ---------- ", k, ": ", MeshR%nodes(l)%vects(1) 
    enddo

    write(*, 307) "el:        "
    do j=1,(NTest2(2)-1)*2,2
      l=(i-1)*(NTest2(2)-1)*2+j
      write(*, 308) " -------- ", l, " -\- ", l +1
    enddo
  enddo 

 105  format (a20)
! 106  format (a,a)
! 205  format (f10.5)
 207  format (a,i10, f10.5)
 208  format(a, i10,a,f10.5)
! 305  format (i10)
 307  format (a)
 308  format (a, i10,a,i10) 
  end subroutine createRoundMesh2

!******************************************************
  subroutine createRoundRand2(NTestMeshR,BCU,BCS)

  use TypesRound

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tRoundMeshInfo)            :: NTestMeshR  ! ne,nn                     !
  type(tRoundRandU)               :: BCU       ! Randbedingungen Biegg       !
  type(tRoundRandS)               :: BCS       ! Randbedingungen Kraft       !
  !                                                                          !
  ! Local variable declaration                                               !
  !                                                                          !
  integer                    :: ne,nn,i,j,nx,ny        ! Zählvariablen       !
  !--------------------------------------------------------------------------!
  intent(in)                 ::  NTestMeshR
  intent(inout)              ::  BCU,BCS                                     !
  !--------------------------------------------------------------------------!
  
  nn = NTestMeshR%nNode
  ne = NTestMeshR%nElem


  ! RHSeite: Fest eingeklammerte
  
  do i= 0,NTest2(1)-1
    j=i*(NTest2(2))+1
    BCU%randFixedIndex(j)=1 
    BCU%randFixed(j) = 0.0 ! 0.1=u,normal  
  enddo

  ! LHSeite, Force
  do i= 0,NTest2(1)-1
    j=i*(NTest2(2))+3
    BCS%randForceIndex(j)=1   ! Fest
    BCS%randForce(j)=0.1  ! Newton
  enddo

  ! Unt. Seite, Moment
  do j= 1,NTest2(2)
    i=j
    BCS%randMomentIndex(i)=1
    BCS%randMoment(i)=-0.1 ! Newton.m
  enddo

  ! Obr. Seite, Moment
  do j= 1,NTest2(2)
    i=2*(NTest2(2))+j
    BCS%randMomentIndex(i)=1
    BCS%randMoment(i)=0.1 ! Newton.m
  enddo

  end subroutine createRoundRand2


!******************************************************
  subroutine createRoundMeshFile(NTestMeshR,MeshR)

  use TypesRound

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tRoundMeshInfo)            :: NTestMeshR  ! ne,nn                     !
  type(tRoundMesh)                :: MeshR       ! Gitterwerte               !
  !                                                                          !
  ! Local variable declaration                                               !
  !                                                                          !
  integer                    :: ne,nn,i,j        ! Zählvariablen             !
  !--------------------------------------------------------------------------!
  intent(inout)              :: NTestMeshR, MeshR                            !
  !--------------------------------------------------------------------------!


  !-----------------------------------<  calculate constants  >-----------
  !---<  nn   = anzahl der inneren punkte                 >-----------
  !---<  ne = anzahl der zellen                           >-----------

!  NTestMeshR%nn = 9
!  NTestMeshR%ne = 8
!  nn = NTestMeshR%nn
!  ne = NTestMeshR%ne

  end subroutine createRoundMeshFile

!******************************************************
  subroutine createRoundRandFile(NTestMeshR,BCU,BCS)

  use TypesRound

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tRoundMeshInfo)            :: NTestMeshR  ! ne,nn                     !
  type(tRoundRandU)               :: BCU       ! Randbedingungen Biegg       !
  type(tRoundRandS)               :: BCS       ! Randbedingungen Kraft       !
  !                                                                          !
  ! Local variable declaration                                               !
  !                                                                          !
  integer                    :: ne,nn,i,j        ! Zählvariablen             !
  !--------------------------------------------------------------------------!
  intent(in)                 ::  NTestMeshR
  intent(inout)              ::  BCU,BCS                                     !
  !--------------------------------------------------------------------------!


  end subroutine createRoundRandFile




end module createmeshround_mod
