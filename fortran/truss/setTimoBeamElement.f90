module setTimoBeamElement_mod

 use TypesBalken
 implicit none

  private
  
  interface linplcreatematrixBeam
    module procedure linplcreatematrixBeam
  end interface linplcreatematrixBeam

  interface linplassemblymatrixBeam
    module procedure linplassemblymatrixBeam
  end interface linplassemblymatrixBeam

  public :: linplcreatematrixBeam, linplassemblymatrixBeam

 contains

!********************************************************
! subroutine linplcreatematrixBeam (pl=piecewise linear
 subroutine linplcreatematrixBeam(MeshBSet,Ke_Bsub,ie)
 use TypesBalken
 implicit none
 
 type(tMeshBeam)       :: MeshBSet
 type(tRigidMat)       :: Ke_Bsub
 integer               :: ie

 ! Ke_Bsub should have been allocated a size 4*6*6=144 
 ! before calling this program

 end subroutine linplcreatematrixBeam

!********************************************************
! subroutine linplassemblymatrixBeam (pl=piecewise linear
 subroutine linplassemblymatrixBeam(MeshBSet,Ke_Bglb,Var_Bglb)
 use TypesBalken
 implicit none
 
 type(tMeshBeam)       :: MeshBSet
 type(tRigidFullMat)   :: Ke_Bglb
 type(tVarFull)        :: Var_Bglb
 integer               :: allocStat

 type(tRigidMat)       :: Ke_Bsub
 type(tVarFull)        :: Var_Bsub

 allocate(Ke_Bsub%Ke(12,12), Var_Bsub%Ue(12), Var_Bsub%Fe(12), stat=allocStat);
 if (allocStat.ne.0) then
     print *," Error allocation Beam sub matrix !"
 endif  


 deallocate(Ke_Bsub%Ke, Var_Bsub%Ue, Var_Bsub%Fe, stat=allocStat); 
 if (allocStat.ne.0) then
     print *," Error deallocation Beam sub matrix !"
 endif  
 end subroutine linplassemblymatrixBeam

end module setTimoBeamElement_mod
