module TypesDefInterface

! make das Interface between old surface structured code for elliptic partial eqn problem
! and the new FEM program, which has to manage also unstructured surface blocks
use Types
! will imort tMesh, tZeiten, from a legacy structured mesh code 

use TypesRound
! will import Circular - triangles meshed geometries,

end module TypesDefInterface
