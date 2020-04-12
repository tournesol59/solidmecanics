module CreateBeamTriangles_mod

  use typesbalken
  implicit none

  private

  public:: allocMesh2DfromFile, createMesh2DfromFile, deallocMesh2DfromFile

  contains
  !****************************
  ! subroutine allocate a field of type Mesh2DTorsion
 subroutine allocMesh2DfromFile(MeshNTest,MeshSec)
  use typesbalken
  implicit none
  type(tMesh2DInfo)            :: MeshNTest
  type(tMesh2DSection)         :: MeshSec
  ! intent
  intent(inout)                :: MeshNTest, MeshSec
  ! lokal
  integer                      :: i,l,nnodes,nelmts,allocStat
  real                         :: x,y,z
  character(len=6)             :: nodetag

  OPEN(UNIT=25, FILE='Lbeamin.msh', ACTION='READ')
  i=1
  read(25,*) nodetag
  do while ((i.le.17).and.(nodetag.eq.("$Nodes")))
    read(25,101) nodetag  ! search for Tag beyond Number of elements line
    i=i+1
  enddo
  i=1
  read(25,102) nnodes
  do i=1,nnodes
    read(25,*) ! # skips since only dimensions interests us for allocation
  enddo
  
  read(25,*)  !  $EndNodes
  read(25,*)  !  $Elements

  read(25,102) nelmts ! number of Elements

  CLOSE(UNIT=25)

  MeshNTest%nnodes=nnodes
  MeshNTest%nelmts=nelmts
  MeshNTest%nboundary=-1  ! not used at least

  allocate(MeshSec%y(nnodes), MeshSec%z(nnodes), MeshSec%elements(nelmts,1:3), &
           MeshSec%neighbours(nelmts,1:3), stat=allocStat)

 101  format (a6)
 102  format (i10)
 103  format (i10 3f8.7)
 104  format (4i10)
 end subroutine allocMesh2DfromFile

  !*****************************
  ! subroutine create a 2D triangles for a rectangular beam section
 subroutine createMesh2DfromFile(MeshNTest,MeshSec)
  use typesbalken
  implicit none
  type(tMesh2DInfo)            :: MeshNTest
  type(tMesh2DSection)         :: MeshSec
  integer                      :: i,j,k,l,m,nnodes,nelmts
  real                         :: x,y,z
  character(len=6)             :: nodetag

  OPEN(UNIT=25, FILE='Lbeamin.msh', ACTION='READ')
  i=1
  read(25,*) nodetag
  do while ((i.le.17).and.(nodetag.eq.("$Nodes")))
    read(25,201) nodetag  ! search for Tag beyond Number of elements line
    i=i+1
  enddo

  do i=1,nnodes
    read(25,203) l,x,y,z
    MeshSec%y(i)=x
    MeshSec%z(i)=y
  enddo
 
  read(25,*)  !  $EndNodes
  read(25,*)  !  $Elements
  read(25,*)  !  ! number of Elements

  do i=1,nelmts
    read(25,204) j,k,l,m
    MeshSec%elements(i,1)=k
    MeshSec%elements(i,2)=l
    MeshSec%elements(i,3)=m
  enddo

  CLOSE(UNIT=25)

 201  format (a6)
! 202  format (i10)
 203  format (i10 3f8.7)
 204  format (4i10)
 end subroutine createMesh2DfromFile

  !****************************
  ! subroutine allocate a field of type Mesh2DTorsion
 subroutine deallocMesh2DfromFile(MeshNTest,MeshSec)
  use typesbalken
  implicit none
  type(tMesh2DInfo)            :: MeshNTest
  type(tMesh2DSection)         :: MeshSec
  
  integer                      :: allocStat
  deallocate(MeshSec%y,MeshSec%z,MeshSec%elements,MeshSec%neighbours, stat=allocStat)
  if (allocStat.ne.0) then
    print *," ---deallocMesh2DfromFile ERROR could not deallocate properly "
  endif
 end subroutine deallocMesh2DfromFile

end module CreateBeamTriangles_mod

