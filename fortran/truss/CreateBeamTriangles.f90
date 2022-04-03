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

    interface  importGrad2DfromFile
      module procedure  importGrad2DfromFile
    end interface

  public:: allocMesh2DfromFile, createMesh2DfromFile, deallocMesh2DfromFile, importGrad2DfromFile


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
  read(25,102) nboundary  ! read after "$PhysicalNames"
  i=i+1
  do l=1,nboundary+2
    read(25,*)  ! skip upt to $Nodes next line
    i=i+1
  enddo
  read(25,102) nnodes  ! read after "$Nodes"
  i=i+1
  do l=1,nnodes+2
    read(25,*) ! # skips since only dimensions interests us for allocation up to $Elements next line
    i=i+1
  enddo
  read(25,102) nelmts ! read after "$Elements" number of Elements
  i=i+1
  CLOSE(UNIT=25)

  MeshNTest%nnodes=nnodes
  MeshNTest%nelmts=nelmts
  MeshNTest%nboundary=nboundary

  allocate(MeshSec%y(1:nnodes), MeshSec%z(1:nnodes), MeshSec%elements(1:nelmts,1:3), &
            stat=allocStat)

 101  format (a6)
 102  format (i4)
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
  integer                      :: ioreadstat
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
!    read(25,203) l,x,y,z  ! A READ ERROR HAS TO BE CORRECTED HERE BEFORE ANY GOING FURTHER
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
  s=1  ! init counter of elmts
  do while ((t.ne.3).and.(i.le.nnodes))
  ! see Gmsh online documentation to see the significance of the first
  ! two numbers, t u (n is elmt index), the others last four are the nodes
    read(25,204,IOSTAT=ioreadstat) n,t,u,v,w,j,k
    if( ioreadstat < 0 )then
       write(6,'(A)') 'Warning: File containts less than 65 entries'
       exit
    else if (ioreadstat >0) then
   ! form two triangles regular based upon the quadrangle read
    MeshSec%elements(s,1)=v
    MeshSec%elements(s,2)=w
    MeshSec%elements(s,3)=k
    s=s+1
    MeshSec%elements(s,1)=w
    MeshSec%elements(s,2)=j
    MeshSec%elements(s,3)=k
    s=s+1
    i=i+1
    endif    
  enddo
  ! note: there was some init code suppressed there, hope this works
  do while ((t.eq.3).and.(i.le.nnodes))
  ! again see Gmsh online Documentation
    read(25,205,IOSTAT=ioreadstat) n,t,u,v,w,j,k,l,m
    if( ioreadstat < 0 )then
       write(6,'(A)') 'Warning: File containts less than 65 entries'
       exit
    else if (ioreadstat >0) then
    MeshSec%elements(s,1)=j
    MeshSec%elements(s,2)=k
    MeshSec%elements(s,3)=m
    s=s+1
    MeshSec%elements(s,1)=k
    MeshSec%elements(s,2)=l
    MeshSec%elements(s,3)=m
    s=s+1
    i=i+1
    endif
  enddo
  CLOSE(UNIT=25)

! set true number of elements: triangles not quadrangles
  MeshNTest%nelmts=s

 201  format (a6)
! 202  format (i10)
 203  format (i5, 3f8.7)
 204  format (7i4)
 205  format (9i4)
 end subroutine createMesh2DfromFile

  !****************************
  ! subroutine allocate a field of type Mesh2DTorsion
 subroutine deallocMesh2DfromFile(MeshNTest,MeshSec)
  use typesbalken
  implicit none
  type(tMesh2DInfo)            :: MeshNTest
  type(tMesh2DSection)         :: MeshSec
  
  integer                      :: allocStat
!  deallocate(MeshSec%y,MeshSec%z,MeshSec%elements, stat=allocStat)
  deallocate(MeshSec%y, stat=allocStat)
  deallocate(MeshSec%z, stat=allocStat)
  deallocate(MeshSec%elements, stat=allocStat)

  if (allocStat.ne.0) then
    print *," ---deallocMesh2DfromFile ERROR could not deallocate properly "
  endif
 end subroutine deallocMesh2DfromFile

 !**************************************************
  ! subroutine create structure for copying the result of FEniCS
 subroutine importGrad2DfromFile(MeshNTest,MeshSec)
  use typesbalken
  implicit none
  type(tMesh2DInfo)            :: MeshNTest
  type(tMesh2DSection)         :: MeshSec

  integer                      :: i,j, nnodes, nelmts, nboundary

  nnodes=MeshNTest%nnodes
  nelmts=MeshNTest%nelmts
  nboundary=MeshNTest%nboundary

!  OPEN(UNIT=25, FILE='Lbeamin.msh', ACTION='READ')
  OPEN(UNIT=25, FILE='squareforprog.msh', ACTION='READ')
  
  CLOSE(UNIT=25)

  end subroutine importGrad2DfromFile


end module CreateBeamTriangles_mod

