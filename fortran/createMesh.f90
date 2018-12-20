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

  public :: createMesh, createMesh2, createMeshGen     !, createMeshGen2

!********************************************************************!
!  Variable Tabelle der Biegung-Plate fuer Einlesen der (10x10)Werte !
!********************************************************************!
  real,dimension(19)        :: X1 = (/ -0.8000, -0.715, -0.6337, -0.545, -0.4589, -0.370, -0.2778, &
            -0.180, -0.0930, 0.0, 0.0930, 0.180, 0.2778, 0.370, 0.4589, 0.545, 0.6337, 0.715, 0.8000 /) 
  real,dimension(19)        :: Z1 = (/ -0.2144, -0.1725, -0.1309, -0.0910, -0.0672, -0.0455, -0.0243, &
            -0.0125, -0.0027, -0.0, -0.0027, -0.0125, -0.0243, -0.0455, -0.0672, -0.0910, -0.1309, -0.1725, -0.2144 /) 

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
!  Mesh%dxq    = 1./Mesh%dx
!  Mesh%dxx    = Mesh%dx*Mesh%dx
!  Mesh%dxxq   = Mesh%dxq*Mesh%dxq
  
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
!  Mesh%dyq    = 1./Mesh%dy
!  Mesh%dyy    = Mesh%dy*Mesh%dy
!  Mesh%dyyq   = Mesh%dyq*Mesh%dyq
  
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
!  Mesh%dxq    = 1./Mesh%dx
!  Mesh%dxx    = Mesh%dx*Mesh%dx
!  Mesh%dxxq   = Mesh%dxq*Mesh%dxq

  do i= 0,Mesh%nraumx-1 
    do j= 0,Mesh%nraumy-1
     Mesh%x(j*(Mesh%nraumx)+i+1)   = X1(i+1)
    enddo
  enddo

  Mesh%starty = 0.0
  Mesh%endy = 1.0
  
  Mesh%dy     = (Mesh%endy - Mesh%starty)/(Mesh%nraumy-1)
!  Mesh%dyq    = 1./Mesh%dy
!  Mesh%dyy    = Mesh%dy*Mesh%dy
!  Mesh%dyyq   = Mesh%dyq*Mesh%dyq

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


  !**************************************************************************!
  !Completion of the mesh data (field neighbour) after imported nodes in .inp!
  !**************************************************************************!
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
  !                                                                          !
  integer                    :: node1, node2, node3                          !
  integer                    :: i,j,k,kinc       ! Zählvariablen             !
  logical                    :: dejavu                                       !
  !--------------------------------------------------------------------------!
  intent(inout)              :: Mesh2D
  !--------------------------------------------------------------------------!

  !** Initialisierung **!
  do i=1,Mesh2D%elmts      
    Mesh2D%neighbours(1,i) = -1
    Mesh2D%neighbours(2,i) = -1
    Mesh2D%neighbours(3,i) = -1
    Mesh2D%neighbours(4,i) = -1
    Mesh2D%neighbours(5,i) = -1
    Mesh2D%neighbours(6,i) = -1
    Mesh2D%neighbours(7,i) = -1
    Mesh2D%neighbours(8,i) = -1
    Mesh2D%neighbours(9,i) = -1
    Mesh2D%neighbours(10,i) = -1
  enddo
  do i=1,Mesh2D%elmts
    node1 = Mesh2D%allelmt(1,i)
    node2 = Mesh2D%allelmt(2,i)
    node3 = Mesh2D%allelmt(3,i)
    kinc=1
  !--- Loop fuer erste Indiz der Punkten auf alle Elementen except i (j=1..j-1,j+1,end)
    do j=1,i-1
       dejavu = .false.
       do k=1,10
         dejavu = (dejavu).or.( (Mesh2D%neighbours(k,i)).eq.(j) )
       enddo
       if ( (.not.(dejavu)).and.((Mesh2D%allelmt(1,j)).eq.(node1)).and.((kinc).le.(11)) ) then
         Mesh2D%neighbours(kinc,i) = j
         kinc=kinc+1
       endif
    enddo

    do j=i+1,Mesh2D%elmts
       dejavu = .false.
       do k=1,10
         dejavu = (dejavu).or.( (Mesh2D%neighbours(k,i)).eq.(j) )
       enddo
       if ( (.not.(dejavu)).and.((Mesh2D%allelmt(1,j)).eq.(node1)).and.((kinc).le.(11)) ) then
         Mesh2D%neighbours(kinc,i) = j
         kinc=kinc+1
       endif
    enddo

  !--- Loop fuer zweite Indiz der Punkten auf alle Elementen except i (j=1..j-1,j+1,end)
    do j=1,i-1
       dejavu = .false.
       do k=1,10
         dejavu = (dejavu).or.( (Mesh2D%neighbours(k,i)).eq.(j) )
       enddo
       if ( (.not.(dejavu)).and.((Mesh2D%allelmt(1,j)).eq.(node2)).and.((kinc).le.(11)) ) then
         Mesh2D%neighbours(kinc,i) = j
         kinc=kinc+1
       endif
    enddo

    do j=i+1,Mesh2D%elmts
       dejavu = .false.
       do k=1,10
         dejavu = (dejavu).or.( (Mesh2D%neighbours(k,i)).eq.(j) )
       enddo
       if ( (.not.(dejavu)).and.((Mesh2D%allelmt(1,j)).eq.(node2)).and.((kinc).le.(11)) ) then
         Mesh2D%neighbours(kinc,i) = j
         kinc=kinc+1
       endif
    enddo

  !--- Loop fuer dritte Indiz der Punkten auf alle Elementen except i (j=1..j-1,j+1,end)
    do j=1,i-1
       dejavu = .false.
       do k=1,10
         dejavu = (dejavu).or.( (Mesh2D%neighbours(k,i)).eq.(j) )
       enddo
       if ( (.not.(dejavu)).and.((Mesh2D%allelmt(1,j)).eq.(node3)).and.((kinc).le.(11)) ) then
         Mesh2D%neighbours(kinc,i) = j
         kinc=kinc+1
       endif
    enddo

    do j=i+1,Mesh2D%elmts
       dejavu = .false.
       do k=1,10
         dejavu = (dejavu).or.( (Mesh2D%neighbours(k,i)).eq.(j) )
       enddo
       if ( (.not.(dejavu)).and.((Mesh2D%allelmt(1,j)).eq.(node3)).and.((kinc).le.(11)) ) then
         Mesh2D%neighbours(kinc,i) = j
         kinc=kinc+1
       endif
    enddo

  enddo  ! Element i

  end subroutine createMeshGen


end module createmesh_mod
