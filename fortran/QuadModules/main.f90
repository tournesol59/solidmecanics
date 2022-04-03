program main

  use types
!  use input_mod
!  use printabstractstruct_mod
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
!  type(tNumeric)             :: VarNum     ! Numerische Felder mit Loesungen !
  !--------------------------------------------------------------------------!   
 ! Elementen Matrizen
!  real,pointer   :: SkalarProdMatrixQ(:,:)  ! full dense Matrix fuer Erfassung des 2. Sys. Gleichung, 1-Dim
!  real,pointer   :: SkalarProdMatrixT(:,:,:)  ! full dense Matrix fuer Erfassung des 1. Sys. Gleichung, 4-Dims

  real,pointer   :: SparseSkalarProdRow1(:) ! sparse matrix csr format: three vectors (1)
! size=c.a. 500
  real,pointer   :: SparseSkalarProdCol(:) ! sparse matrix csr format: three vectors (2)
  real,pointer   :: SparseSkalarProdVal(:) ! sparse matrix csr format: three vectors (3)
                   ! Alle Punkten der Area 2-zeit differenzierbar!
  !--------------------------------------------------------------------------!  
!  Testen Lagrange1d/hermite1d functions
  integer        :: i,j,n, allocStat
  real,pointer,intent(inout) :: matrixN(:,:)  ! (1:8,1:3)
  ! -------------------------------------------------------------< BEGIN >---!
  !read(5,*) MeshPlanar

  ! --------------------------< Eingabegroessen einlesen: 2 Alternativen >---!
  ! ------------------------------------------< Speicherplatz allokieren >---!

  call allocateFields(NTestGitter,Gitter,RB,Uvar,rhs,Exakt)
  !call allocateFieldsCurv(NTestGitter,KurvVect)
  !call allocateFieldsVarNum(NTestGitter,VarNum)
  !call allocateGitter2D(Gitter2D)  
  !call linlg2allocatesdns(Gitter, SkalarProdMatrixQ, SkalarProdMatrixT)
  call linlgallocmatrixpatch5el(Gitter2D, Uvar, rhs, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal)

  !allocate(matrixN(1:8,1:3), STAT=allocStat)
  ! ! --------------------------< Gitter festlegen >---------------------!
   Const%auto = 1
   !if (Const%auto == 1) then
   !   call createMesh(Gitter)    
   !else 
   !   call createMesh2(Gitter)
   !endif

   if (Const%auto == 1) then  
     call createMesh5el(Gitter2D,RB)
   else 
      call createMeshGen2(Gitter2D)
   endif

  !write(*,*) '================    Verifikation Kondnrg starts now    ================ '

  !write(*,*) '================    Interpolation starts now    ================ '

  !write(*,*) '================    Matrix ausfuellen starts now    ================ '
  !----------< Matrix fuellen :: generic name common to Lagrange // Hermite>-----------!
 write(*,*) '================    Start Patch test 5 elements  ================ '

 do ne=1,5
    call linlgfillmatrixpatch5el(Gitter2D,matrixN,ne)
    call linlgintegmatrixpatch5el(Gitter2D,ne)
 enddo
 102  format (f10.5)

!!! call linlgintegmatrixpatch5el(Gitter, chicoeff, VarNum, SkalarProdMatrixQ)

  ! ------------------------------------------< Speicherplatz de-allokieren >---!

 !   deallocate(SkalarProdMatrixQ, SkalarProdMatrixT, STAT = allocStat )
    deallocate(SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal, &
           STAT = allocStat);

  call deallocateFields(Const,Gitter,RB,Exakt,Uvar,rhs)
 ! call deallocateFieldsCurv(KurvVect)
 ! call deallocateFieldsVarNum(VarNum)
 ! call deallocateGitter2D(Gitter2D)  
 ! deallocate(matrixN, STAT=allocStat)
  call linlgdeallocmatrixpatch5el(Gitter2D, Uvar, rhs, SparseSkalarProdRow1, SparseSkalarProdCol, SparseSkalarProdVal)      
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateSparseMatrix: Could not deallocate correctly!'
       STOP
    end if 

    write(*,*) 
    write(*,*) '============= Program terminated correctly ================ '

end program main
