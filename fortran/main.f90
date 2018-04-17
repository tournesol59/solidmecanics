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
  real,pointer   :: SparseSkalarProdRow1(:) ! sparse matrix csr format: three vectors (1)
  real,pointer   :: SparseSkalarProdCol(:) ! sparse matrix csr format: three vectors (2)
  real,pointer   :: SparseSkalarProdVal(:) ! sparse matrix csr format: three vectors (3)
  real,pointer   :: chicoeff(:,:)
  real,pointer   :: der_a1xx_x(:),der_a1xx_y(:),der_a1xy_x(:),der_a1xy_y(:) &
                        ,der_a2yx_x(:),der_a2yx_y(:),der_a2yy_x(:),der_a2yy_y(:)
!  Testen Lagrange1d functions
  type(tPolynom) :: pol_pts, pol_root1, pol_nom1, pol_value, pol_lagr, pol_der, pol_herm

  integer        :: i,j, allocStat
  ! ------------------------------------------< Eingabegroessen einlesen >---!

  call input(Gitter,Exakt,Const,FileIO)       

  ! ------------------------------------------< Speicherplatz allokieren >---!

  call allocatefields(Gitter,RB,Uvar,rhs,Exakt,chicoeff,Const,VarNum)  

  call linlgallocatesparse(Gitter, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal) 

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
   call trychristoffei(Gitter, chicoeff)

  ! --------------------------< Temporary test of one module >---------------------!
  write(*,*) '================    Verification module Lagrange1d starts now ======== '
  pol_pts%n=3
  pol_pts%coeffs=  (/ -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

  pol_root1%n=3
  pol_root1%coeffs=  (/ -1.0, -2.0, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
  call expandfromroots(pol_root1, pol_nom1)

  pol_value%n=3
  pol_value%coeffs=  (/ -1.0, 0.5, -0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  call expandlagrange(pol_pts, pol_value, pol_lagr)

  write(*,*) '================    Verification module Hermite1d starts now ======== '

  pol_der%n=3
  pol_der%coeffs=  (/ 0.5, 0.1, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  call expandhermite(pol_pts, pol_value, pol_der, pol_herm)
  ! -----------< Calculation of curvature, possibility of interpolation after>-----!

  write(*,*) '================    Interpolation starts now    ================ '
  write(*,*)
! call from Interpol2d
  call basisdifferenz(Gitter,der_a1xx_x,der_a1xx_y,der_a1xy_x,der_a1xy_y &
                                 ,der_a2yy_x, der_a2yy_y, der_a2yx_x, der_a2yx_y)

  write(*,*) '     differenz der lokalen Basis kompletiert    '
! call from Interpol2d
  call trychristoffei(Gitter, chicoeff)

! call from Interpol2d
  call christoffei(Gitter, chicoeff, der_a1xx_x,der_a1xx_y,der_a1xy_x,der_a1xy_y &
                                 ,der_a2yy_x, der_a2yy_y, der_a2yx_x, der_a2yx_y)

  write(*,*) '     Christoffei Koeffizienten berechnet   '
! call from Interpol2d
  call trychristoffei(Gitter, chicoeff)


  !----------< Matrix fuellen :: generic name common to Lagrange // Hermite>-----------!

  call linlgfillmatrix(Gitter, Const, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal)

  ! ------------------------------------------< Speicherplatz allokieren >---!

  !call deallocateFields(Gitter,RB,Uvar,rhs,Exakt,chicoeff,Const) 
    deallocate(Gitter%x, Gitter%y, RB%randl, RB%randr, RB%rando, RB%randu, &
         Exakt%loesung, Uvar, rhs, chicoeff, STAT = allocStat )
    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not deallocate correctly!'
       STOP
    end if 

   !call from LagrangeGrad2d:
    deallocate(SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal, &
           STAT = allocStat);
      
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateSparseMatrix: Could not deallocate correctly!'
       STOP
    end if 

  deallocate(der_a1xx_x)
  deallocate(der_a1xx_y)

  deallocate(der_a1xy_x)
  deallocate(der_a1xy_y)

  deallocate(der_a2yx_x)
  deallocate(der_a2yx_y)

  deallocate(der_a2yy_x)
  deallocate(der_a2yy_y)

    write(*,*) 
    write(*,*) '============= Program terminated correctly ================ '

end program main
