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

  public :: createMesh, createMesh2, createMeshGen, createMeshGen2

!********************************************************************!
!  Variable Tabelle der Biegung-Plate fuer Einlesen der (10x10)Werte !
!********************************************************************!
  real,dimension(10)        :: X1 = (/ -0.8000, -0.6337, -0.4589, -0.2778, &
            -0.0930, 0.0930, 0.2778, 0.4589, 0.6337, 0.8000 /) 
  real,dimension(10)        :: Z1 = (/ -0.2144, -0.1309, -0.0672, -0.0243, &
            -0.0027, -0.0027, -0.0243, -0.0672, -0.1309, -0.2144 /) 

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
  do i= 0,9 
    do j= 0,9
     Mesh2D%x(j*10+i+1)   = X1(i+1)
     Mesh2D%z(j*10+i+1)   = Z1(i+1)
     Mesh2D%y(j*10+i+1)   = 0.0 + j*1.11111
     n=n+1
    enddo
  enddo 
!  Mesh2D%nodes = n

  el=0
  do i= 0,8
    do j= 0,8
      k = j*10+i+1
      Mesh2D%quad(1, j*9+i+1)   = k
      Mesh2D%quad(2, j*9+i+1)   = k+1
      Mesh2D%quad(3, j*9+i+1)   = k+10+1
      Mesh2D%quad(4, j*9+i+1)   = k+10
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
