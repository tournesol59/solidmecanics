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
module createmesh_mod
  implicit none

  public :: createMesh, createMesh2, createMesh5el, createMeshGen, createMeshGen2

!********************************************************************!
!  Variable Tabelle der Biegung-Plate fuer Einlesen der (10x10)Werte !
!********************************************************************!
  real,dimension(19)        :: X1 = (/ -0.8000, -0.715, -0.6337, -0.545, -0.4589, -0.370, -0.2778, &
            -0.180, -0.0930, 0.0, 0.0930, 0.180, 0.2778, 0.370, 0.4589, 0.545, 0.6337, 0.715, 0.8000 /) 
  real,dimension(19)        :: Z1 = (/ -0.2144, -0.1725, -0.1309, -0.0910, -0.0672, -0.0455, -0.0243, &
            -0.0125, -0.0027, -0.0, -0.0027, -0.0125, -0.0243, -0.0455, -0.0672, -0.0910, -0.1309, -0.1725, -0.2144 /) 


!********************************************************************!
!  Variable Tabelle der 5 Elmts Patch Test fuer Einlesen der (10x10)Werte !
!********************************************************************!
  real,dimension(8)        :: patchXe = (/ 0.0, 3.0, 3.0, 0.0, 0.75, 2.25, 2.25, 0.75 /)
  real,dimension(8)        :: patchYe = (/ 0.0, 0.0, 2.0, 2.0, 0.5 , 0.75, 1.5, 1.5 /)

  contains

  subroutine createMesh(Mesh)

  use types

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tMesh)                :: Mesh       ! Gitterwerte                     !
  !                                                                          !
  ! Local variable declaration                                               !
  !                                                                          !
  integer                    :: i,j        ! Zählvariablen                   !
  !--------------------------------------------------------------------------!
  !!! intent(in)                 :: Const
  intent(inout)              :: Mesh
  !--------------------------------------------------------------------------!


  !-----------------------------------<  calculate constants  >-----------
  !
  !---<  nraumx   = anzahl der inneren punkte                 >-----------
  !---<  nraumx+1 = anzahl der inneren zellen                 >-----------

  Mesh%dx     = (Mesh%endx-Mesh%startx)/(Mesh%nraumx-1)
  Mesh%dxq    = 1./Mesh%dx
  Mesh%dxx    = Mesh%dx*Mesh%dx
  Mesh%dxxq   = Mesh%dxq*Mesh%dxq
  
  !---<  x(0)             = linker  rand                      >-----------
  !---<  x(nraumx+1)      = rechter rand                      >-----------
  !---<  x(i) i=1,nraumx  = innere  punkte                    >-----------
  do i= 0,Mesh%nraumx-1 
    do j= 0,Mesh%nraumy-1
     Mesh%x(j*(Mesh%nraumx)+i+1)   = Mesh%startx + i*Mesh%dx
    enddo
  enddo

  !---<  nraumy   = anzahl der inneren punkte                 >-----------
  !---<  nraumy+1 = anzahl der inneren zellen                 >-----------
  
  Mesh%dy     = (Mesh%endy - Mesh%starty)/(Mesh%nraumy-1)
  Mesh%dyq    = 1./Mesh%dy
  Mesh%dyy    = Mesh%dy*Mesh%dy
  Mesh%dyyq   = Mesh%dyq*Mesh%dyq
  
  !---<  y(0)             = unterer rand                      >-----------
  !---<  y(nraumy+1)      = oberer  rand                      >-----------
  !---<  y(i) i=1,nraumy  = innere  punkte                    >-----------
  
   do j= 0,Mesh%nraumy-1 
    do i= 0,Mesh%nraumx-1
     Mesh%y(j*(Mesh%nraumx)+i+1)   = Mesh%starty + j*Mesh%dy
    enddo
  enddo

  !erase dx,dxq,dxxq, dy,dyq,dyyq
!  Mesh%dx=1.0
!  Mesh%dxq=1.0
!  Mesh%dxx=1.0
!  Mesh%dxxq=1.0
!  Mesh%dy=1.0
!  Mesh%dyq=1.0
!  Mesh%dyy=1.0
!  Mesh%dyyq=1.0
  ! -----------------------------------<  laenge x, laenge y   >-----------
  
  Mesh%xlen = Mesh%endx - Mesh%startx
  Mesh%ylen = Mesh%endy - Mesh%starty
  
  end subroutine createMesh



  !--------------------------------------------------------------------------!
  !  Automatic creation of the mesh with the data in the header              !
  !--------------------------------------------------------------------------!
  subroutine createMesh2(Mesh)

  use types

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tMesh)                :: Mesh       ! Gitterwerte                     !
  !                                                                          !
  ! Local variable declaration                                               !
  !    
  !
  integer                    :: i,j        ! Zählvariablen                   !
  !--------------------------------------------------------------------------!
  !!! intent(in)                 :: Const
  intent(inout)              :: Mesh
  !--------------------------------------------------------------------------!


  !-----------------------------------<  set constants  >-----------
  !---<  nraumx   = anzahl der inneren punkte                 >-----------
  !---<  nraumx+1 = anzahl der inneren zellen                 >-----------
  Mesh%nraumx = 10
  Mesh%nraumy = 10
   
  Mesh%startx = X1(1)
  Mesh%endx = X1(10)

  Mesh%dx     = (Mesh%endx-Mesh%startx)/(Mesh%nraumx-1)
  Mesh%dxq    = 1./Mesh%dx
  Mesh%dxx    = Mesh%dx*Mesh%dx
  Mesh%dxxq   = Mesh%dxq*Mesh%dxq

  do i= 0,Mesh%nraumx-1 
    do j= 0,Mesh%nraumy-1
     Mesh%x(j*(Mesh%nraumx)+i+1)   = X1(i+1)
    enddo
  enddo

  Mesh%starty = 0.0
  Mesh%endy = 1.0
  
  Mesh%dy     = (Mesh%endy - Mesh%starty)/(Mesh%nraumy-1)
  Mesh%dyq    = 1./Mesh%dy
  Mesh%dyy    = Mesh%dy*Mesh%dy
  Mesh%dyyq   = Mesh%dyq*Mesh%dyq

   do j= 0,Mesh%nraumy-1 
    do i= 0,Mesh%nraumx-1
     Mesh%y(j*(Mesh%nraumx)+i+1)   = Mesh%starty + j*Mesh%dy
    enddo
  enddo
  
  !------------ no information about dz,dzq,dzz,dzzq is used ------!

  do i= 0,Mesh%nraumx-1 
     do j= 0,Mesh%nraumy-1
     Mesh%z(j*(Mesh%nraumx)+i+1)   = Z1(i+1)
    enddo
  enddo
  
  Mesh%xlen = Mesh%endx - Mesh%startx
  Mesh%ylen = Mesh%endy - Mesh%starty
 end subroutine createMesh2


  subroutine createMeshGen(Mesh2D)

  use types

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tMeshGen)             :: Mesh2D     ! Gitterwerte                     !
  !                                                                          !
  ! Local variable declaration                                               !
  !    
  !
  integer                    :: i,j        ! Zählvariablen                   !
  !--------------------------------------------------------------------------!
  intent(inout)              :: Mesh2D
  !--------------------------------------------------------------------------!
  
  OPEN(UNIT=25, FILE='Untitled1.msh', ACTION='READ')
  ! OPERATIONEN HIER ...  !
  read (25,*)
  CLOSE(UNIT=25)

 105  format (a20)
 106  format (a,a)
 205  format (f10.5)
 206  format (a,f10.5)
 305  format (i10)
 306  format (a,i10)

  end subroutine createMeshGen


  !--------------------------------------------------------------------------!
  !  Automatic creation of the mesh fuer Patch 5 Elmts with the data in the header              !
  !--------------------------------------------------------------------------!
  subroutine createMesh5el(Mesh2D,RB)

  use types

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tMeshGen)                :: Mesh2D       ! Gitterwerte                !
  type(tRandbedingungen)        :: RB           ! RandB struktur             !
  !                                                                          !
  ! Local variable declaration                                               !
  !    
  !
  integer                    :: i,j,n        ! Zählvariablen                   !
  !--------------------------------------------------------------------------!
  !!! intent(in)                 :: Const
  intent(inout)              :: Mesh2D
  !--------------------------------------------------------------------------!
  n=0
  do i= 1,8
      Mesh2D%x(i)   = patchXe(i)
      Mesh2D%z(i)   = 0.0
      Mesh2D%y(i)   = patchYe(i)
      n=n+1
  enddo 
  Mesh2D%nodes = n
  Mesh2D%elmts = 5
! Indiz der Punkte fuer 1. Quadrangle !
  Mesh2D%quad(1,1)=1
  Mesh2D%quad(1,2)=2
  Mesh2D%quad(1,3)=6
  Mesh2D%quad(1,4)=5
! Indiz der Punkte fuer 2. Quadrangle !
  Mesh2D%quad(2,1)=6
  Mesh2D%quad(2,2)=2
  Mesh2D%quad(2,3)=3
  Mesh2D%quad(2,4)=7
! Indiz der Punkte fuer 3. Quadrangle !
  Mesh2D%quad(3,1)=7
  Mesh2D%quad(3,2)=3
  Mesh2D%quad(3,3)=4
  Mesh2D%quad(3,4)=8
! Indiz der Punkte fuer 4. Quadrangle !
  Mesh2D%quad(4,1)=1
  Mesh2D%quad(4,2)=5
  Mesh2D%quad(4,3)=8
  Mesh2D%quad(4,4)=4
! Indiz der Punkte fuer 5. Quadrangle (Mitte) !
  Mesh2D%quad(5,1)=5
  Mesh2D%quad(5,2)=6
  Mesh2D%quad(5,3)=7
  Mesh2D%quad(5,4)=8
! Indiz der Nachbar Elemente fuer 1. Quadrangle convention: -1 ist Randbedingung!
  Mesh2D%neighbours(1,1)=-1 ! unten
  Mesh2D%neighbours(1,2)=2 ! rechts
  Mesh2D%neighbours(1,3)=5 ! top
  Mesh2D%neighbours(1,4)=4 ! links: hier Vertiz (1,4) linken Nachbar
! Indiz der Nachbar Elemente fuer 2. Quadrangle !
  Mesh2D%neighbours(2,1)=1 ! unten
  Mesh2D%neighbours(2,2)=-1 ! rechts
  Mesh2D%neighbours(2,3)=3 ! top
  Mesh2D%neighbours(2,4)=5 ! links: hier Vertiz (1,4) linken Nachbar
! Indiz der Nachbar Elemente fuer 3. Quadrangle !
  Mesh2D%neighbours(2,1)=5 ! unten
  Mesh2D%neighbours(2,2)=3 ! rechts
  Mesh2D%neighbours(2,3)=-1 ! top
  Mesh2D%neighbours(2,4)=4 ! links: hier Vertiz (1,4) linken Nachbar
!! Indiz der Nachbar Elemente fuer 4. Quadrangle !
  Mesh2D%neighbours(2,1)=1 ! unten
  Mesh2D%neighbours(2,2)=5 ! rechts
  Mesh2D%neighbours(2,3)=3 ! top
  Mesh2D%neighbours(2,4)=-1 ! links: hier Vertiz (1,4) linken Nachbar
! Indiz der Nachbar Elemente fuer 5. Quadrangle !
  Mesh2D%neighbours(2,1)=1 ! unten
  Mesh2D%neighbours(2,2)=2 ! rechts
  Mesh2D%neighbours(2,3)=3 ! top
  Mesh2D%neighbours(2,4)=4 ! links: hier Vertiz (1,4) linken Nachbar
!
! Indiz der Elemente an boundaries
  Mesh2D%ntop=1
  Mesh2D%quadtop(1)=3   ! 3st PatchTest Elmt ist im Top
  Mesh2D%nbottom=1
  Mesh2D%quadbottom(1)=1 ! 1st PatchTest Elmt ist unten
  Mesh2D%nleft=1
  Mesh2D%quadleft(1)=4  ! 4st PatchTest Elmt hat links Rand
  Mesh2D%nright=1
  Mesh2D%quadright(1)=2 ! 2st PatchTest Elmt hat rechts Rand
  
  RB%randltype=1      ! RandTyp links: 1=Fixed, 2=Force
  RB%randl(1)=0.0       ! Randwerte Ux links, der nur ein Ecke ist
  RB%randl(2)=0.0       ! Randwerte Uy links, der nur ein Ecke..
  RB%randutype=2      ! RandTyp unten: 2=Force
  RB%randu(1)=0.5       ! Randwerte Ux 
  RB%randu(2)=0.0       ! Randwerte Uy 
  RB%randrtype=1      ! RandTyp rechts: 1=Fixed
  RB%randr(1)=0.0       ! Randwerte Ux
  RB%randr(2)=0.0       ! Randwerte Uy 
  RB%randotype=2      ! RandTyp oben: 2=Force
  RB%rando(1)=0.5       ! Randwerte Ux
  RB%rando(2)=0.0       ! Randwerte Uy 
 
 end subroutine createMesh5el

  !**************************************************************************!
  !  subroutine  createMeshGen2:                                           !
  !                                                                          !
  !  Berechnet die Gitter2D Punkten (Punkten-Element Format)                 !
  ! funktioniert nur aber, wenn 400 Punkten und 81 Elementen                 !
  ! (8-Punkt-Elementen) vorhanden sind.                                      !
  !**************************************************************************!
  subroutine createMeshGen2(Mesh2D)

  use types

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tMeshGen)             :: Mesh2D     ! Gitterwerte                     !
  !                                                                          !
  ! Local variable declaration                                               !
  !    
  !
  integer                    :: i,j,k,n,el ! Zählvariablen                   !
  !--------------------------------------------------------------------------!
  intent(inout)              :: Mesh2D
  !--------------------------------------------------------------------------!
   
  n=0
  do i= 0,18 ! Fuellt 361 Punkten aus aber 81 davon werden nicht benutzt
    do j= 0,18
      Mesh2D%x(j*19+i+1)   = X1(i+1)  ! Index min ist 1, max ist 361 OK 
      Mesh2D%z(j*19+i+1)   = Z1(i+1)
      Mesh2D%y(j*19+i+1)   = -1.0 + real(j,8)*0.11111
      n=n+1
    enddo
  enddo 
!  Mesh2D%nodes = n

  el=0
  do i= 0,8
    do j= 0,8
      k = 2*j*19+2*i+1    ! min ist 1, max ist 321
      Mesh2D%quad(1, j*9+i+1)   = k
      Mesh2D%quad(3, j*9+i+1)   = k+2
      Mesh2D%quad(5, j*9+i+1)   = k+40    ! Index min ist 41, max ist 361=19*19 OK
      Mesh2D%quad(7, j*9+i+1)   = k+38

      Mesh2D%quad(2, j*9+i+1)   = k+1
      Mesh2D%quad(4, j*9+i+1)   = k+21
      Mesh2D%quad(6, j*9+i+1)   = k+39
      Mesh2D%quad(8, j*9+i+1)   = k+19
      el=el+1
    enddo
  enddo
!  Mesh2D%elmts = el

  n=0
  do i= 0,8
    k = 8*9+i+1
    Mesh2D%quadtop(i+1)   = k
    n=n+1
  enddo
!  Mesh2D%ntop = n

  n=0
  do i= 0,8
    k = i+1
    Mesh2D%quadbottom(i+1)   = k
    n=n+1
  enddo
!  Mesh2D%nbottom = n

  n=0
  do j= 0,8
    k = j*9+9+1
    Mesh2D%quadright(j+1)   = k
    n=n+1
  enddo
!  Mesh2D%nright = n

  n=0
  do j= 0,8
    k = j*9+1
    Mesh2D%quadleft(j+1)   = k
    n=n+1
  enddo
!  Mesh2D%nleft = n


  end subroutine createMeshGen2


end module createmesh_mod
