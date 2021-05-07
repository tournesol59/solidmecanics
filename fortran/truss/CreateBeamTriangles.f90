module CreateBeamTriangles_mod

  use typesbalken
  implicit none

  private

#define BeELnum_Ycoord(m,n)	 (MeshSec%dy*(m-1))  
! for a rectangular global form section only
#define BeELnum_Zcoord(m,n)      (MeshSec%dz*(n-1))

    interface allocMesh2DfromFile
      module procedure allocMesh2DfromFile
    end interface

    interface createMesh2DfromFile
      module procedure createMesh2DfromFile
    end interface

    interface deallocMesh2DfromFile
      module procedure deallocMesh2DfromFile
    end interface

  public:: allocMesh2DfromFile, createMesh2DfromFile, deallocMesh2DfromFile

  ! Macros definition for access to mesh2DSection
  #define BEAM2D_NAME_OF_STRUCTURE_index_plus_index_of_neighbou 

  contains
  !****************************
  ! subroutine allocate a field of type tMesh2DSection
 subroutine allocMesh2DfromFile(MeshNTest,MeshSec)
  use typesbalken
  implicit none
  type(tMesh2DInfo)            :: MeshNTest
  type(tMesh2DSection)         :: MeshSec
  ! intent
  intent(inout)                :: MeshNTest, MeshSec
  ! lokal
  integer                      :: i,l,nnodes,nelmts,nboundary,allocStat
  real                         :: x,y,z
  character(len=6)             :: nodetag

!  OPEN(UNIT=25, FILE='Lbeamin.msh', ACTION='READ')
  OPEN(UNIT=25, FILE='squareforprog.msh', ACTION='READ')
  i=1
  do l=1,4
    read(25,*) ! skip up to $PhysicalNames next line
    i=i+1 
  enddo
!  read(25,*) nodetag
!  do while ((i.le.17).and.(nodetag.eq.("$Nodes")))
!    read(25,101) nodetag  ! search for Tag beyond Number of elements line
!    i=i+1
!  enddo
  read(25,102) nboundary
  i=i+1
  do l=1,nboundary+2
    read(25,*)  ! skip upt to $Nodes next line
    i=i+1
  enddo
  read(25,102) nnodes
  i=i+1
  do l=1,nnodes+2
    read(25,*) ! # skips since only dimensions interests us for allocation up to $Elements next line
    i=i+1
  enddo
  read(25,102) nelmts ! number of Elements
  i=i+1
  CLOSE(UNIT=25)

  MeshNTest%nnodes=nnodes
  MeshNTest%nelmts=nelmts
  MeshNTest%nboundary=nboundary

  allocate(MeshSec%y(nnodes), MeshSec%z(nnodes), MeshSec%elements(nelmts,1:3), &
            stat=allocStat)

 101  format (a6)
 102  format (i10)
 103  format (i10, 3f8.7)
 104  format (4i10)
 end subroutine allocMesh2DfromFile

  !*****************************
  ! subroutine create a 2D triangles for a rectangular beam section
  ! reads a mesh file of a square or rectangular geom
 subroutine createMesh2DfromFile(MeshNTest,MeshSec)
  use typesbalken
  implicit none
  type(tMesh2DInfo)            :: MeshNTest
  type(tMesh2DSection)         :: MeshSec
  integer                      :: i,j,k,l,m,n,t,u,v,w,s,nnodes,nelmts,nboundary
  real                         :: x,y,z
  character(len=6)             :: nodetag


  nnodes=MeshNTest%nnodes
  nelmts=MeshNTest%nelmts
  nboundary=MeshNTest%nboundary

!  OPEN(UNIT=25, FILE='Lbeamin.msh', ACTION='READ')
  OPEN(UNIT=25, FILE='squareforprog.msh', ACTION='READ')
  do i=1, 5+nboundary+3
    read(25,*)  ! Skip these lines 
  enddo
  j=1
  do i=5+nboundary+4, 5+nboundary+4+nnodes-1
!    read(25,203) l,x,y,z
    x=0.0
    y=0.0
    read(25,*) 
!   not a field from struct,x
    MeshSec%y(j)=x ! yes
    MeshSec%z(j)=y !yes
    j=j+1
  enddo
  do i=5+nboundary+4+nnodes,  5+nboundary+4+nnodes+2
    read(25,*)  ! Skip these lines 
  enddo
  i=1
  t=1
  do while ((t.ne.3).and.(i.le.nnodes))
    read(25,204) n,t,u,v,w,j,k
    i=i+1
  enddo
  s=1 ! init counter of elements
  MeshSec%elements(s,1)=j  ! form two triangles regular based upon the quadrangle read
  MeshSec%elements(s,2)=k
  MeshSec%elements(s,3)=m
  s=s+1
  MeshSec%elements(s,1)=k
  MeshSec%elements(s,2)=l
  MeshSec%elements(s,3)=m
  s=s+1
  do while ((t.eq.3).and.(i.le.nnodes))
    read(25,204) n,t,u,v,w,j,k,l,m
    i=i+1
    MeshSec%elements(s,1)=j
    MeshSec%elements(s,2)=k
    MeshSec%elements(s,3)=m
    s=s+1
    MeshSec%elements(s,1)=k
    MeshSec%elements(s,2)=l
    MeshSec%elements(s,3)=m
    s=s+1
  enddo
  CLOSE(UNIT=25)

! set true number of elements: triangles not quadrangles
  MeshNTest%nelmts=s

 201  format (a6)
! 202  format (i10)
 203  format (i5, 3f8.7)
 204  format (4i5)
 205  format (7i5)
 end subroutine createMesh2DfromFile

  !****************************
  ! subroutine allocate a field of type Mesh2DTorsion
 subroutine deallocMesh2DfromFile(MeshNTest,MeshSec)
  use typesbalken
  implicit none
  type(tMesh2DInfo)            :: MeshNTest
  type(tMesh2DSection)         :: MeshSec
  
  integer                      :: allocStat
  deallocate(MeshSec%y,MeshSec%z,MeshSec%elements, stat=allocStat)
  if (allocStat.ne.0) then
    print *," ---deallocMesh2DfromFile ERROR could not deallocate properly "
  endif
 end subroutine deallocMesh2DfromFile

end module CreateBeamTriangles_mod

