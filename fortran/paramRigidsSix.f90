module paramRigidSix_mod
!***********************************************************************!
! both: Rigidity calculation with the hand-made precalculated integrated!
!          expression (general, not deep into membrane, flexion, shear  !
!          separated parameters). Valable for 6 nodes triangular elmt   !
!       at the present shall be used for mass-matrix calc. Reuse in the 
!       future for a DKT formulation

  use TYPESROUND
  implicit none


  private

  interface massMatrix_M_E_local
     module procedure massMatrix_M_E_local
  end interface

  interface massMatrix_M_H_global
     module procedure massMatrix_M_H_global
  end interface

  module 

  end module procedure


  public :: massMatrix_M_E_local, massMatrix_M_H_global

  contains


end paramRigidSix_mod
