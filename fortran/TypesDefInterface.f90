module TypesDefInterface
implicit none
! make das Interface between old surface structured code for elliptic partial eqn problem
! and the new FEM program, which has to manage also unstructured surface blocks
use TypesDef
! will imort tMesh, tZeiten, from a legacy structured mesh code 

use TypesRoundDef
! will import Circular - triangles meshed geometries,

end TypesDefInterface
