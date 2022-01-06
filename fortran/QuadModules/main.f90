program main

  use types
!  use input_mod
  use printabstractstruct_mod
  use allocatefields_mod
  use createmesh_mod
!!  use interpol2d_mod
!!  use Lagrangegrad2d_mod
  use Patchtest2d5el_mod


  implicit none

  ! logical: 0: Ebene Flaeche, 1: mit Kurvatur Flaeche
  logical      :: MeshPlanar
  ! Referenz Flaeche ist ein Rechteck-quadratnetz mit NRAUMX*NRAUMY Elementen

  ! Parameter der Interpolation
  integer      :: ORDER=3   ! 2 before
  integer      :: allocStat
  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  type(tMeshInfo)            :: NTestGitter ! Gitterwerte in matriz           !
  type(tMesh)                :: Gitter     ! Gitterwerte in matriz           !
  type(tMeshGen)             :: Gitter2D   ! Gitterwerte in Generic inputfile!
  type(tCurv)                :: KurvVect   ! Numerische Felder fuer Kurve, Christoffei Koeffs abhangig, Rechnung !
  type(tRandbedingungen)     :: RB         ! Randbedingungen Code            !
  type(tConstants)           :: Const      ! Konstanten                      !
  type(tFileIO)              :: FileIO     ! Ausgabesteuerung                !
  type(tExakt)               :: Exakt      ! Exakte Loesung                  !
  real,pointer               :: Uvar(:)  ! Feld mit der numerischen Loesung!!=Tr(spanng) !
  real,pointer               :: rhs(:)   ! Feld für die Rechte Seite       !
  type(tNumeric)             :: VarNum     ! Numerische Felder mit Loesungen !
  !--------------------------------------------------------------------------!   
  ! Elementen Matrizen
  real,pointer   :: matrixN(:,:)
  real,pointer   :: SkalarProdMatrixQ(:,:)  ! full dense Matrix fuer Erfassung des 2. Sys. Gleichung, 1-Dim
  real,pointer   :: SkalarProdMatrixT(:,:,:)  ! full dense Matrix fuer Erfassung des 1. Sys. Gleichung, 4-Dims

  real,pointer   :: SparseSkalarProdRow1(:) ! sparse matrix csr format: three vectors (1)
! size=c.a. 500
  real,pointer   :: SparseSkalarProdCol(:) ! sparse matrix csr format: three vectors (2)
  real,pointer   :: SparseSkalarProdVal(:) ! sparse matrix csr format: three vectors (3)
                   ! Alle Punkten der Area 2-zeit differenzierbar!
  !--------------------------------------------------------------------------!  
!  Testen Lagrange1d/hermite1d functions
  type(tPolynom) :: pol_pts, pol_root1, pol_nom1, pol_value, pol_lagr, pol_der, pol_herm
  type(tPolynom) :: pol1, pol2, pol3
  real,pointer   :: Coefficients(:,:)  ! Referenz interpolation polynoms coeffsX,Y, size=2*100
  real           :: Jac, val1, val2
  integer        :: i,j,n, allocStat
  real,pointer,intent(inout) :: matrixN(:,:)  ! (1:8,1:3)
  ! -------------------------------------------------------------< BEGIN >---!
  !read(5,*) MeshPlanar

  ! --------------------------< Eingabegroessen einlesen: 2 Alternativen >---!
  ! -- commented here since Patch Test data are constants into code
!  call input(Gitter,Exakt,Const,FileIO)       
!  ! OR
!  call input_file(Gitter2D,Exakt,Const,FileIO)
!  Gitter2D%nodes=361
!  Gitter2D%elmts=81
!  Gitter2D%ntop = 9
!  Gitter2D%nbottom = 9
!  Gitter2D%nright = 9
!  Gitter2D%nleft = 9
  ! ------------------------------------------< Speicherplatz allokieren >---!

  !call allocateFields(NTestGitter,Gitter,RB,Uvar,rhs,Exakt)
  !call allocateFieldsCurv(NTestGitter,KurvVect)
  !call allocateFieldsVarNum(NTestGitter,VarNum)
  !call allocateGitter2D(Gitter2D)  
  !call linlg2allocatesdns(Gitter, SkalarProdMatrixQ, SkalarProdMatrixT)
  linlgallocmatrixpatch5el(Gitter2D, Uvar, rhs, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal)
  !call linlgallocmatrixpatch5el(Gitter2D, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal)
  !allocate(matrixN(1:8,1:3), STAT=allocStat)
  ! ! --------------------------< Gitter festlegen >---------------------!
   !if (Const%auto == 1) then
   !   call createMesh(Gitter)    
   !else 
   !   call createMesh2(Gitter)
   !endif

   !if (Const%auto == 1) then  
   !  call createMeshGen(Gitter2D)
   !else 
   !   call createMeshGen2(Gitter2D)
   !endif

  ! call trychristoffei(Gitter, chicoeff)

  ! --------------------------< check Referenz Netz interpolieren >---------------------!
 !  call interpolref(2, Gitter);
 !  call lgevalderiv3(2,0.3,0.15,1,Coefficients,Jac)

  !write(*,*) '================    Verifikation Kondnrg starts now    ================ '
  !----------< Alle Element Geometrie verifizieren :: wenn nur eins ausfaellt, dann krasht es >-----------!
  !call  checkalledgesofelmts(Gitter2D)

  ! -----------< Calculation of curvature, possibility of interpolation after>-----!

  !write(*,*) '================    Interpolation starts now    ================ '
  !write(*,*)
! call from Interpol2d
!!!  call basisdifferenz(Gitter,der_a1xx_x,der_a1xx_y,der_a1xy_x,der_a1xy_y &
!!!                                ,der_a2yy_x, der_a2yy_y, der_a2yx_x, der_a2yx_y)
  
  !write(*,*) '     differenz der lokalen Basis kompletiert    '
! call from Interpol2d
!!!  call trychristoffei(Gitter, chicoeff)

! call from Interpol2d
!!!  call christoffei(Gitter, chicoeff, der_a1xx_x,der_a1xx_y,der_a1xy_x,der_a1xy_y &
!!!                                 ,der_a2yy_x, der_a2yy_y, der_a2yx_x, der_a2yx_y)

 ! write(*,*) '     Christoffei Koeffizienten berechnet   '
! call from Interpol2d
 ! call trychristoffei(Gitter, chicoeff)

  !write(*,*) '================    Matrix ausfuellen starts now    ================ '
  !----------< Matrix fuellen :: generic name common to Lagrange // Hermite>-----------!
 write(*,*) '================    Start Patch test 5 elements  ================ '
 call linlgallocmatrixpatch5el(Gitter2D, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal)

 call linlgfillmatrixpatch5el(Gitter2D,matrixN,ne)
 102  format (f10.5)

!!! call linlgintegmatrixpatch5el(Gitter, chicoeff, VarNum, SkalarProdMatrixQ)

  ! ------------------------------------------< Speicherplatz de-allokieren >---!

 !   deallocate(SkalarProdMatrixQ, SkalarProdMatrixT, STAT = allocStat )
 !   deallocate(SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal, &
 !          STAT = allocStat);

 ! call deallocateFields(Const,Gitter,RB,Exakt,Uvar,rhs)
 ! call deallocateFieldsCurv(KurvVect)
 ! call deallocateFieldsVarNum(VarNum)
 ! call deallocateGitter2D(Gitter2D)  
 ! deallocate(matrixN, STAT=allocStat)
  call linlgdeallocmatrixpatch5el(Gitter2D, Uvar, rhs, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal)      
 !   if (allocStat.NE.0) then
 !      print *, 'ERROR AllocateSparseMatrix: Could not deallocate correctly!'
 !      STOP
 !   end if 

    write(*,*) 
    write(*,*) '============= Program terminated correctly ================ '

end program main
