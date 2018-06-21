program main

  use types
  use input_mod
  use allocatefields_mod
  use createmesh_mod
  use interpol2d_mod
  use Lagrange2d_mod
  use Lagrangegrad2d_mod
  use Lagrange1d_mod
  use hermite1d_mod
  use hermite2d_mod

  implicit none

 ! Parameter des Gitter,netz: Ebene auf eine Fläche
 ! integer      :: NRAUMX=20
 ! integer      :: NRAUMY=10
 ! real         :: STARTX=0.0
 ! real         :: ENDX=30.0
 ! real         :: STARTY=0.0
 ! real         :: ENDY=15.0
 ! Referenz Flaeche ist ein Rechteck-quadratnetz mit NRAUMX*NRAUMY Elementen

  ! Parameter der Interpolation
  integer      :: ORDER=3   ! 2 before

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  type(tMesh)                :: Gitter     ! Gitterwerte                     !
  type(tRandbedingungen)     :: RB         ! Randbedingungen                 !
  type(tConstants)           :: Const      ! Konstanten                      !
  type(tFileIO)              :: FileIO     ! Ausgabesteuerung                !
  type(tExakt)               :: Exakt      ! Exakte Loesung                  !
  real,pointer               :: Uvar(:,:)  ! Feld mit der numerischen Loesung!!=Tr(spanng) !
  real,pointer               :: rhs(:,:)   ! Feld für die Rechte Seite       !
  type(tNumeric)             :: VarNum     ! Numerische Felder mit Loesungen !
  !--------------------------------------------------------------------------!   
  real,pointer   :: SkalarProdMatrixQ(:,:)  ! full dense Matrix fuer Erfassung des 2. Sys. Gleichung, 1-Dim
  real,pointer   :: SkalarProdMatrixT(:,:,:)  ! full dense Matrix fuer Erfassung des 1. Sys. Gleichung, 4-Dims

  real,pointer   :: SparseSkalarProdRow1(:) ! sparse matrix csr format: three vectors (1)
! size=c.a. 500
  real,pointer   :: SparseSkalarProdCol(:) ! sparse matrix csr format: three vectors (2)
  real,pointer   :: SparseSkalarProdVal(:) ! sparse matrix csr format: three vectors (3)
  real,pointer   :: chicoeff(:,:)          ! 8 Values fuer jedes Punktes: Christoffel Koeffizienten
  real,pointer   :: der_a1xx_x(:),der_a1xx_y(:),der_a1xy_x(:),der_a1xy_y(:) &
                        ,der_a2yx_x(:),der_a2yx_y(:),der_a2yy_x(:),der_a2yy_y(:)  
                   ! Alle Punkten der Area 2-zeit differenzierbar!
  !--------------------------------------------------------------------------!  
!  Testen Lagrange1d/hermite1d functions
  type(tPolynom) :: pol_pts, pol_root1, pol_nom1, pol_value, pol_lagr, pol_der, pol_herm
  real,pointer   :: Coefficients(:,:)  ! Referenz interpolation polynoms coeffsX,Y, size=2*100
  real           :: Jac
  integer        :: i,j,n, allocStat

  ! -------------------------------------------------------------< BEGIN >---!
  ! ------------------------------------------< Eingabegroessen einlesen >---!

  call input(Gitter,Exakt,Const,FileIO)       

  ! ------------------------------------------< Speicherplatz allokieren >---!

  call allocatefields(Const,Gitter,RB,Uvar,rhs,Exakt,chicoeff,VarNum)  
  call linlg2allocatesdns(Gitter, SkalarProdMatrixQ, SkalarProdMatrixT)
  !call linlgallocatesparse(Gitter, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal) 

  allocate(der_a1xx_x(100))
  allocate(der_a1xx_y(100))

  allocate(der_a1xy_x(100))
  allocate(der_a1xy_y(100))

  allocate(der_a2yx_x(100))
  allocate(der_a2yx_y(100))

  allocate(der_a2yy_x(100))
  allocate(der_a2yy_y(100))

   ! --------------------------< Gitter festlegen >---------------------!
   if (Const%auto == 1) then
      call createMesh(Gitter)    
   else 
      call createMesh2(Gitter)
   endif
  ! call trychristoffei(Gitter, chicoeff)


  ! --------------------------< check Referenz Netz interpolieren >---------------------!
 !  call interpolref(2, Gitter);
 !  call lgevalderiv3(2,0.3,0.15,1,Coefficients,Jac)

  ! --------------------------< Temporary test of one module >---------------------!
  !write(*,*) '================    Verification module Lagrange1d starts now ======== '
  pol_pts%n=3
  pol_pts%coeffs=  (/ -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

!  pol_root1%n=3
!  pol_root1%coeffs=  (/ -1.0, -2.0, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!  call expandfromroots(pol_root1, pol_nom1)

  pol_value%n=3
  pol_value%coeffs=  (/ -1.0, 0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  !call expandlagrange(pol_pts, pol_value, pol_lagr)
  pol_value%n=3
  pol_value%coeffs=  (/ -1.0, 0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  !call expandlagrange(pol_pts, pol_value, pol_lagr)
  !write(*,*) '================    Verification module Hermite1d starts now ======== '

  pol_der%n=3
  pol_der%coeffs=  (/ 0.5, 0.1, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  call expandhermite(pol_pts, pol_value, pol_der, pol_herm)
  n=3

  ! call calchermitepol(n, Coefficients)

  ! -----------< Calculation of curvature, possibility of interpolation after>-----!

  write(*,*) '================    Interpolation starts now    ================ '
  write(*,*)
! call from Interpol2d
  call basisdifferenz(Gitter,der_a1xx_x,der_a1xx_y,der_a1xy_x,der_a1xy_y &
                                 ,der_a2yy_x, der_a2yy_y, der_a2yx_x, der_a2yx_y)
  
  write(*,*) '     differenz der lokalen Basis kompletiert    '
! call from Interpol2d
  !call trychristoffei(Gitter, chicoeff)

! call from Interpol2d
  call christoffei(Gitter, chicoeff, der_a1xx_x,der_a1xx_y,der_a1xy_x,der_a1xy_y &
                                 ,der_a2yy_x, der_a2yy_y, der_a2yx_x, der_a2yx_y)

  write(*,*) '     Christoffei Koeffizienten berechnet   '
! call from Interpol2d
!  call trychristoffei(Gitter, chicoeff)


  !----------< Matrix fuellen :: generic name common to Lagrange // Hermite>-----------!

 ! call linlgfillmatrix(Gitter, Const, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal)
 ! call linlg2dfillmatrixq(Gitter, chicoeff, VarNum, SkalarProdMatrixQ)
  ! ------------------------------------------< Speicherplatz allokieren >---!

    deallocate(der_a1xx_x)
    deallocate(der_a1xx_y)

    deallocate(der_a1xy_x)
    deallocate(der_a1xy_y)

    deallocate(der_a2yx_x)
    deallocate(der_a2yx_y)

    deallocate(der_a2yy_x)
    deallocate(der_a2yy_y)

 !   deallocate(SkalarProdMatrixQ, SkalarProdMatrixT, STAT = allocStat )
 !   deallocate(SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal, &
 !          STAT = allocStat);

 !  call deallocateFields(Gitter, RB, Exakt, Uvar, rhs, chicoeff, Const, VarNum )
    
   ! if (allocStat.NE.0) then
   !    print *, 'ERROR AllocateFields: Could not deallocate correctly!'
   !    STOP
   ! end if 

      
 !   if (allocStat.NE.0) then
 !      print *, 'ERROR AllocateSparseMatrix: Could not deallocate correctly!'
 !      STOP
 !   end if 

    write(*,*) 
    write(*,*) '============= Program terminated correctly ================ '

end program main
