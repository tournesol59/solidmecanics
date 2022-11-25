module LagrangePlan2D_mod
! this module shall be used to perform a Planestress model
! for Quad elements. The model uses elementary subroutines and polynom
! operation from Patchtest5el2d.f90 and data structure subroutines from
! createMesh.f90 and re use qnd re write as much as possible LagrangeGrad2D.f90

  implicit none

  use Patchtest5el2d_mod

! access macros to get node from elemt shall be in createMesh.f90
!#define ELEMi_NODEj(i,j)  Mesh2D%elmts(i)%nodes(j)

   public :: linlgallocmatrixsparse, linlgcreatestructquadgen
 contains
! ********************
! linlgallocmatrixsparse

  subroutine linlgallocmatrixsparse(MeshInfo2D,SparseMatProdRow1, SparseMatProdCol, SparseMatProVal)
  use types
  implicit none
  tMeshInfo         :: Mesh2DInfo
  real, pointer     :: SparseMatProdRow1(:)
  real, pointer     :: SparseMatProdCol(:)
  real, pointer     :: SparseMatProVald(:)
  integer           :: allocStat, nx, ny

  nx=Mesh2DInfo%nx
  ny=Mesh2DInfo%ny
  allocate(SparseMatProdRow1(1:2*nx*ny), SparseMatProdCol(1:2*nx*ny), &
             SparseMatProVal(1:((7*2+6*2))*nx*ny), INFO=allocStat)


  end subroutine linlgallocmatrixsparse

!****************************************
! linlgcreatestructpatch5el
!
subroutine linlgcreatestructquadgen(MeshInfo, RB, RBtyp, RBconst)
  use types
  implicit none
  type(tMeshInfo)        :: MeshInfo
  type(tRandbedingungen) :: RB
  integer                :: RBtyp
  real,dimension(1:4)    :: RBconst
! local variable
  integer                :: i

    ! constant values are copied
  RB%randltype = RBtyp  !0=fixed displace, 1= tensile force
  do i=1,MeshInfo%nraumy
    RB%randl(i)  = RBconst(1)
    RB%randr(i)  = RBconst(3)
  enddo
  do i=1,MeshInfo%nraumx
    RB%rando(i)  = RBconst(2)
    RB%randu(i)  = RBconst(4)
  enddo

  end subroutine linlgcreatestructquadgen

end module LagrangePlan2D_mod
