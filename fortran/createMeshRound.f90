!**************************************************************************!
!                                                                          !
! subroutine createMesh                                                    !
!                                                                          !
!      Funktion: Festlegen der Gitterdaten                                 !
!                                                                          !
!**************************************************************************!

module createmeshround_mod
  use TypesRound

  implicit none

  private
!**************************************************************************!
!  Variable Tabelle der Patch Test-Plate fur die Koordinaten der 9 Punkten !
!**************************************************************************!
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
        MeshR%elems(l)%numeros(2) = k+3 !k+4
        MeshR%elems(l)%numeros(3) = k+2 !k+3
        MeshR%elems(l)%numeros(4) = k  !+1
        MeshR%neighbours(l)%numeros(1)=l+1
        MeshR%neighbours(l)%numeros(2)=l-1+(NTest2(2)-1)*2
        MeshR%neighbours(l)%numeros(1)=l-1
      else if ((j).eq.3) then
        MeshR%elems(l)%numeros(2) = k-1 ! k+1
        MeshR%elems(l)%numeros(3) = k !k+2   
        MeshR%elems(l)%numeros(4) = k+2 !+4
        MeshR%neighbours(l)%numeros(1)=l-1
        MeshR%neighbours(l)%numeros(2)=-1
        MeshR%neighbours(l)%numeros(3)=l+1   
      else if ((j).eq.4) then
        MeshR%elems(l)%numeros(2) = k+1 !k+5      
        MeshR%elems(l)%numeros(3) = k  !k+4   
        MeshR%elems(l)%numeros(4) = k-2 !k+2   
        MeshR%neighbours(l)%numeros(1)=-1
        MeshR%neighbours(l)%numeros(2)=l-1+(NTest2(2)-1)*2
        MeshR%neighbours(l)%numeros(3)=l-1
     endif
    enddo
  enddo 

   write(*, 105) "--------- Kordinaten von Mesh: --------- "
   do i= 1,(NTest2(1)-1) 
    write(*, 207) "y: ", i, MeshR%nodes(i*NTest2(2))%vects(2)
    do j=1,(NTest2(2))
      k=(i-1)*NTest2(2) + j  ! 1..3 4..6 7..9 Linie von 3 Punkten. Sieht X1,Y1,Z1  !
      l=(i-1)*(NTest2(2)-1)*2+j
      write(*, 208) " -------x-y ", k, ": ", MeshR%nodes(k)%vects(1), ": ", MeshR%nodes(k)%vects(2) 
    enddo

    write(*, 307) "elem:        "
    do j=1,(NTest2(2)-1)*2,2
      l=(i-1)*(NTest2(2)-1)*2+j
      write(*, 308) " -------- ", l, " -\- ", l +1
    enddo
  enddo 
  write(*, 207) "y: ", i, MeshR%nodes(i*NTest2(2))%vects(2)
  i=3
  do j=1,3
      k=(i-1)*NTest2(2) + j  ! 1..3 4..6 7..9 Linie von 3 Punkten. Sieht X1,Y1,Z1  !
      l=(i-1)*(NTest2(2)-1)*2+j
      write(*, 208) " -------x-y ", k, ": ", MeshR%nodes(k)%vects(1), ": ", MeshR%nodes(k)%vects(2) 
  enddo
 105  format (a20)
! 106  format (a,a)
! 205  format (f10.5)
 207  format (a,i10, f10.5)
 208  format(a, i10,a,f10.5,a,f10.5)
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
    BCU%randFixedIndex(j)=1 ! fixed diplacement type
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

   print *,"--------- Randbedingungen"
   print *," Typ--Seite(Axis)     Ind(i/j)   Ind     values"
  do i= 0,NTest2(1)-1
    j=i*(NTest2(2))+1
    write(*, 217) " BCU--x1-y1-3 ", j, ": ",     BCU%randFixedIndex(j), ": ", BCU%randFixed(j) 
  enddo
  do i= 0,NTest2(1)-1
    j=i*(NTest2(2))+3
   write(*, 217) " BCS--x3-y1-3 ", j, ": ",     BCS%randForceIndex(j), ": ", BCS%randForce(j) 
  enddo
  do j= 1,NTest2(2)
    i=j
   write(*, 217) " BCS--x1-3-y1 ", i, ": ",     BCS%randForceIndex(i), ": ", BCS%randForce(i) 
  enddo
  do j= 1,NTest2(2)
    i=2*(NTest2(2))+j
   write(*, 217) " BCS--x1-3-y3 ", i, ": ",     BCS%randMomentIndex(i), ": ", BCS%randMoment(i) 
  enddo

217  format(a, i10,a,i10,a,f10.5)


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
  character(LEN=36)          :: Rawtext
  integer                    :: no,el,elmn,elmx,i,j        ! Zählvariablen             !
  !--------------------------------------------------------------------------!
  intent(inout)              :: NTestMeshR, MeshR                            !
  !--------------------------------------------------------------------------!


  !-----------------------------------<  calculate constants  >-----------
  
  OPEN(UNIT=25, FILE='Untitled.inp', ACTION='READ')
  ! .. couts Nodes and Elements and Boudary Element...

  read (25,*)                ! *Heading
  read (25,*)                !  /home/fredrik/Documents/gmsh/untitled3.inp
  read (25,*)                ! 2.2 0 8
  read (25,*)                ! *NODE

  do i=1,NTestMeshR%nNode
     read(25,207) no, MeshR%nodes(i)%vects(1), MeshR%nodes(i)%vects(2), MeshR%nodes(i)%vects(3)    ! what to do if no!=i?
  enddo
  if (no.ne.NTestMeshR%nNode) then
    print *, "nodes spec in .msh is not equal to nodes lines read in .inp Input File"
  endif
  read (25,*)                ! ******* E L E M E N T S *************

!!! crucial: search for
! *ELEMENT, type=CPS4, ELSET=Surface1
  do while ((Rawtext).NE.("*ELEMENT, type=T3D2, ELSET=Line1"))
    read (25,106) 
  end do

  do j=1,NTestMeshR%nElem
     read(25,307) MeshR%elems(j)%numeros(1), MeshR%elems(j)%numeros(2), MeshR%elems(j)%numeros(3), &
                   MeshR%elems(j)%numeros(4)             ! Skip it, for further, go on createMesh.f90
     if (j.eq.1) then
       elmn=MeshR%elems(j)%numeros(1)
     elseif (j.eq.NTestMeshR%nElem) then
       elmx=MeshR%elems(j)%numeros(1)
     endif
  enddo
  if ((elmx-elmn+1).ne.NTestMeshR%nElem) then
    print *, "elements spec in .msh is not equal to elements lines read in .inp Input File"
  endif  
  
  CLOSE(UNIT=25)

 106  format (a)
 205  format (f10.5)
 206  format (a,f10.5)
 207  format (i10,3f10.5)
 305  format (i10)
 306  format (a,i10)
 307  format (4i10)

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
  character(LEN=33)          :: Rawtext
  integer                    :: ne,nn,i,j        ! Zählvariablen             !
  !--------------------------------------------------------------------------!
  intent(inout)              ::  NTestMeshR, BCU,BCS                         !
  !--------------------------------------------------------------------------!
  
  OPEN(UNIT=25, FILE='Untitled.inp', ACTION='READ')
  ! .. opens a first time to count boundary nodes...

  read (25,*)                ! *Heading
  read (25,*)                !  /home/fredrik/Documents/gmsh/untitled3.inp
  read (25,*)                ! 2.2 0 8
  read (25,*)                ! *NODE

  do i=1,NTestMeshR%nNode
     read(25,*)
  enddo
  read (25,*)                ! $EndNodes

!!! crucial: search for
! *ELEMENT, type=CPS4, ELSET=Surface1
  i=1
  read (25,406)
  do while ((Rawtext).NE.("*ELEMENT, type=T3D2, ELSET=Line1"))
     i=i+1
     read (25,406) 
  end do
  NTestMeshR%BoundVert(1)=i-1
  i=1
  read (25,406)
  do while ((Rawtext).NE.("*ELEMENT, type=T3D2, ELSET=Line2"))
     i=i+1
     read (25,406) 
  end do
  NTestMeshR%BoundVert(2)=i-1
  i=1
  read (25,406)
  do while ((Rawtext).NE.("*ELEMENT, type=T3D2, ELSET=Line3"))
     i=i+1
     read (25,406) 
  end do
  NTestMeshR%BoundVert(3)=i-1
  i=1
  read (25,406)
  do while ((Rawtext).NE.("*ELEMENT, type=T3D2, ELSET=Line4"))
     i=i+1
     read (25,406) 
  end do
  NTestMeshR%BoundVert(4)=i-1

  CLOSE(UNIT=25)

  
  OPEN(UNIT=25, FILE='Untitled.inp', ACTION='READ')
  ! .. opens a second time to read boundary nodes...


  CLOSE(UNIT=25)

 406  format (a)
 705  format (i10)
 706  format (a,i10)
 707  format (4i10)

  end subroutine createRoundRandFile




end module createmeshround_mod
