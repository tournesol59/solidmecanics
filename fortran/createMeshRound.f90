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
 integer,dimension(2) :: NTest2 = (/ 3, 3/) ! dimension of Patch Test Mesh, then coordinates:
 
 real,dimension(9)  ::  X2 = (/ 0.0, 0.05, 0.1, 0.0, 0.05, 0.1, 0.0, 0.05, 0.1 /)

 real,dimension(9)  ::  Y2 = (/ 0.0, 0.0, 0.0, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1 /)

 real,dimension(9)  ::  Z2 = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

 integer,pointer    ::  RenumTable(:)

#define nodeI_isalreadypaired_J(i,j) (RenumTable(i)==j)
#define elemInodeJ(i,j) MeshR%elems(i)%numeros(j)
#define neighbourIJK(i,j,k) MeshR%elems(MeshR%neighbours(i)%numeros(j))%numeros(k)
  public :: createRoundMesh2, createRoundMeshFile, &
            createRoundRand2, createRoundRandFile, &
            allocateRenum, neighbourVertice, & 
            renumerotate, deallocateRenum

  contains
!****************
! allocate Tables specific to this module
  subroutine allocateRenum(NTestMeshR)
  use TypesRound

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  type(tRoundMeshInfo)            :: NTestMeshR  ! ne,nn                     !
  integer                         :: infoStat

  ! comment. is nNode sufficent or should we have a multiple of 3*nElem?
  allocate(RenumTable(1:NTestMeshR%nNode), STAT=infoStat)
  if (infoStat.ne.0) then
          print *,"error allocateRenum"
  endif
  end subroutine allocateRenum

!****************
  subroutine deallocateRenum
  use TypesRound

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  integer                         :: infoStat

  deallocate(RenumTable, STAT=infoStat)
  if (infoStat.ne.0) then
          print *,"error deallocateRenum"
  endif
  end subroutine deallocateRenum

!******************************************************
! neighbourVertice: 
  subroutine neighbourVertice(MeshR, elem, elemnei, numvert)
        ! numvert is by convention the neighbour vertex: 1 if (1,2) is 
        ! the common vertice, 2 if (1,3) and 3 if (2,3) -1 by found no vertice
  use TypesRound
  
  implicit none
  type(tRoundMesh)                :: MeshR  ! ne,nn                     !
  integer                         :: elem, elemnei, numvert
  integer                         :: i,j,k,a,b,c
  intent(in)                      :: MeshR, elem, elemnei
  intent(inout)                   :: numvert
  
  i=MeshR%elems(elem)%numeros(1)
  j=MeshR%elems(elem)%numeros(2)
  k=MeshR%elems(elem)%numeros(3)
  a=MeshR%elems(elemnei)%numeros(1)
  b=MeshR%elems(elemnei)%numeros(2)
  c=MeshR%elems(elemnei)%numeros(3)
  
  if (((i.eq.a).and.(j.eq.b)).or.((i.eq.b).and.(j.eq.a))) then
     numvert=1
  elseif (((i.eq.a).and.(k.eq.b)).or.((i.eq.b).and.(k.eq.a))) then
     numvert=1
  elseif (((j.eq.a).and.(k.eq.b)).or.((j.eq.b).and.(k.eq.a))) then
     numvert=1
  elseif (((i.eq.a).and.(j.eq.c)).or.((i.eq.c).and.(j.eq.a))) then
     numvert=2
  elseif (((i.eq.a).and.(k.eq.c)).or.((i.eq.c).and.(k.eq.a))) then
     numvert=2
  elseif (((j.eq.a).and.(k.eq.c)).or.((j.eq.c).and.(k.eq.a))) then
     numvert=2
  elseif (((i.eq.b).and.(j.eq.c)).or.((i.eq.c).and.(j.eq.b))) then
     numvert=3
  elseif (((i.eq.b).and.(k.eq.c)).or.((i.eq.c).and.(k.eq.b))) then
     numvert=3
  elseif (((j.eq.b).and.(k.eq.c)).or.((j.eq.c).and.(k.eq.b))) then
     numvert=3
  endif

  end subroutine neighbourVertice

!******************************************************
! renumerotate
  subroutine renumerotate(NTestMeshR,MeshR)
 
  use TypesRound

  implicit none
  type(tRoundMeshInfo)         :: NTestMeshR
  type(tRoundMesh)             :: MeshR
  integer                      :: i,j,k,l,m,numvert
! and operate on global variable RenumTable
  m=1
  l=1
  do while ((l.le.NTestMeshR%nElem).and.(m.le.NTestMeshR%nElem)) ! end cond
       do i=1,3 ! three sides of the element
         k=0
         do while ((k.le.m).and.(nodeI_isalreadypaired_J(k,elemInodeJ(l,i)))) ! see #define at top
            k=k+1 ! test if the nodes of elem(l) are already stored in a pair in RenumTable
              ! because m goes forward
         enddo
         if ((k.eq.m).and.(MeshR%neighbours(l)%numeros(i).ne.(-1))) then ! no boundary
                 ! if not found k<m -> call to func neighbourVertice
            call neighbourVertice(MeshR,l,MeshR%neighbours(l)%numeros(i),numvert) ! return an int code in numvert:
            select case (numvert) 
               case (1)
                  RenumTable(m) = elemInodeJ(l,1)
                  RenumTable(m+1) = elemInodeJ(l,2) 
                  RenumTable(m+2) = neighbourIJK(l,i,3) ! opposite node to face
                  m=m+3 
               case (2)
                  RenumTable(m) = elemInodeJ(l,1)
                  RenumTable(m+1) = elemInodeJ(l,3)
                  RenumTable(m+2) = neighbourIJK(l,i,2)
                  m=m+3
               case (3)
                  RenumTable(m) = elemInodeJ(l,2)
                  RenumTable(m+1) = elemInodeJ(l,3)
                  RenumTable(m+2) = neighbourIJK(l,i,1)
                  m=m+3
               case (-1) ! boundary side
                  RenumTable(m) = elemInodeJ(l,j)
                  RenumTable(m+1) = elemInodeJ(l,mod((j+1),3))
                  m=m+2
             end select
         endif
       enddo
     l=l+1
  enddo
! PRINT RESULT (by line of six)
 print *, " -- createMeshRound renumerotate --"
 do i=0,NTestMeshR%nNode/6-1
    write(*,807) ( j, j=6*i,6*(i+1) )
    write(*,807) ( RenumTable(j), j=6*i,6*(i+1) )
    print *, "  "
 enddo
 write(*,807) ( j, j=6*(NTestMeshR%nNode/6)+1,(NTestMeshR%nNode) )
 write(*,807) ( RenumTable(j), j=6*(NTestMeshR%nNode/6)+1,(NTestMeshR%nNode) )
 807  format (6i10)

 end subroutine renumerotate

!******************************************************
!createRoundMesh2:  generates 9 nodes and 8 elements simple mesh for Patch Test
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
! create simple boundary condition structures for Patch Test
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
  
  OPEN(UNIT=25, FILE='patch12equi.inp', ACTION='READ')
  ! .. couts Nodes and Elements and Boudary Element...

  read (25,*)                ! *Heading
  read (25,*)                !  /home/fredrik/Documents/gmsh/untitled2.inp
  read (25,*)                ! *NODE
  write(*,305) NTestMeshR%nNode
  do i=1,NTestMeshR%nNode
     read(25,207) no, MeshR%nodes(i)%vects(1), MeshR%nodes(i)%vects(2), MeshR%nodes(i)%vects(3)    ! what to do if no!=i?
  enddo
  if (no.ne.NTestMeshR%nNode) then
    print *, "nodes spec in .msh is not equal to nodes lines read in .inp Input File"
  endif
  read (25,*)                ! ******* E L E M E N T S *************

!!! crucial: search for
    read (25,106) Rawtext
  do while ((Rawtext).NE.("*ELEMENT, type=CPS3, ELSET=Surface1"))
    read (25,106) Rawtext
  end do

  j=1
  elmx=1
  elmn=1
  do j=1,NTestMeshR%nElem
     read(25,307) MeshR%elems(j)%numeros(1), MeshR%elems(j)%numeros(2), MeshR%elems(j)%numeros(3), &
                   MeshR%elems(j)%numeros(4)          
  !   if (elmn.ge.( MeshR%elems(j)%numeros(1) )) then
  !     elmn=MeshR%elems(j)%numeros(1)
  !   elseif (elmx.le.( MeshR%elems(j)%numeros(1) )) then
  !     elmx=MeshR%elems(j)%numeros(1)
  !   endif

  enddo
!  if (j.ne.NTestMeshR%nElem) then
!    print *, "elements spec in .msh is not equal to elements lines read in .inp Input File"
!  endif  
  
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
