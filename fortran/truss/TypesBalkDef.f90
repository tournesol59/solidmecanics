!*****************************************************************************!
!* 
!*      $RCSfile: TypesDef.f90,v $
!*
!*      $Revision: 1.3 $
!*
!*      $Date: 2004/02/13 08:16:09 $
!*
!*      $Author: iagklima $
!*
!*****************************************************************************!

MODULE TypesBalken
  !---------------------------------------------------------------------------!
  use ISO_C_BINDING
  implicit none
  private
  !---------------------------------------------------------------------------!

  type tMeshElmt
     integer     :: nraumx         ! subdivsion, default=2
     real        :: startx         ! idem  !
     real        :: endx           ! idem  !
     real        :: xlen           ! Laengen des Rechengebiets (Balk) in x/y !
     real        :: SArea          ! Section area
     real        :: EY             ! YOUNG Modulus
     real        :: vu             ! Poisson Coeff
     real        :: CI             ! kinetics inertia moment
     real        :: dx             ! Gitterschrittweite (Delta x)    !
     real        :: dxq            ! Kehrwert ( 1/ Delta x )         !
     real        :: dxx            ! dx^2 (Delta x * Delta x)        !
     real        :: dxxq           ! Kehrwert des Quadrats (1/dxx)   !
     ! real,pointer (c_ptr) :: x(:)        ! do NOT work
     real,dimension(2)         :: x,y,z

     real,dimension(:),allocatable         :: CoeffsH1(:)
     real,dimension(:),allocatable         :: CoeffsH2(:)
     real,dimension(:),allocatable         :: CoeffsH3(:)
     real,dimension(:),allocatable         :: CoeffsH4(:)
  end type tMeshElmt

  type tVarElmt
     real,dimension(6)         :: UeIi, UeJi
     real,dimension(6)         :: FeIi,FeJi
  end type tVarElmt

  type tRigidMat
     real,dimension(6,6)       :: Ke
  end type tRigidMat

  type tPolynom
    integer               :: n              ! degree =>  dimension n+1        !
    real,dimension(7)    :: coeffs(1:7)      ! Koeffizienten  order<10         !
  end type tPolynom

  !---------------------------------------------------------------------------!
  public  :: tMeshElmt, tVarElmt, tRigidMat, tPolynom
  !---------------------------------------------------------------------------!

! contains

end module TypesBalken
