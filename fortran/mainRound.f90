program main

  use typesround
  use input_mod
  use allocatefieldsround_mod
  use createmeshround_mod
  use verifymeshround_mod
  use spacegeometry_mod
!  use lagrangeGrad2DRound_mod
  use lagrangeround2d_mod
  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  type(tRoundMeshInfo)            :: NTestGitterR ! Gitterwerte in matriz           !
  type(tRoundMesh)                :: GitterR      ! Gitterwerte in matriz           !
  type(tRoundNumeric)             :: VarFelderR   ! Numerische Felder mit Loesungen !
  type(tRoundExakt)               :: ExaktLR      ! Exakte Loesung                  !
!  type(tRoundFileIO)              :: FileIOR      ! Ausgabesteuerung                !
  type(tRoundCoeff)               :: GitterCoeffR ! Interpolation
  type(tRoundRandU)               :: RBU          ! Randbedingungen Biegung
  type(tRoundRandS)               :: RBS          ! Randbedingungen Force
  type(tPassageMatrix)            :: axref_xyz, axout_xyz, Pmat_xyz    ! a 3x3 Matrix contains local axes 
  type(tSparseYMat)               :: RigidYMatrix
                                                  ! vs global axes coordinates!
 ! type(tNode),dimension(1:3)      :: vects        ! contains 3x 3 Vectors
  type(tNode)                     :: vx,vy,vz     ! Zum test of check in verifymeshRound!
  real,dimension(1:3)             :: dist, angle, triplet
  real,dimension(1:4)             :: area
  type(tRoundCoeff)               :: C_evalmembrane ! Alle K_ sind steifigkeitsmatrix, C_ sind nur Teile
  real,pointer                    :: K_membrane(:,:)

 !------------------------------------------------------------------------------
 ! Variable for Mixed C (callee)-Fortran (caller) Programming
 !--------------------------------------------------------------------------!  

! Other
   real           :: x,y
   integer        :: i,j,k,l,m,n,p,ne,ni, allocStat

  ! -------------------------------------------------------------< BEGIN >---!
 ! call input(Gitter,Exakt,Const,FileIO)
  call Input_inpfile(NTestGitterR,GitterR)   ! untitled.inp 
  ! ------------------------------------------< Speicherplatz allokieren >---!
!  NTestGitterR%nNode = 9
!  NTestGitterR%nElem = 8
!  NTestGitterR%nxmx = 9 ! frontier length should be optimized in the future ..
!  NTestGitterR%nymx = 9  ! but essential in the present. If not: seg. fault!
 ! allocate(vx%vects(3), vy%vects(3), vz%vects(3), STAT=allocStat) 
 ! if (allocStat.NE.0) then
 !    print *, 'ERROR Allocate vx vy vz: Could not allocate all variables!'
 !    STOP
 ! end if
  call allocateFieldsRound(NTestGitterR,GitterR,ExaktLR)
  call allocateRenum(NTestGitterR)
  call allocateFieldsBCSRound(NTestGitterR,RBU,RBS)
  call allocateFieldsVarRound(NTestGitterR,VarFelderR)
!  call allocateFieldsCoeffRound(NTestGitterR,GitterCoeffR) 
  call allocateSparseRound(NTestGitterR, RigidYMatrix)

   ! --------------------------< Gitter festlegen >---------------------!
   call createRoundMesh2(NTestGitterR,GitterR)    
!   call createRoundRand2(NTestGitterR,RBU,RBS)
!   call createRoundMeshFile(NTestGitterR,GitterR)   ! patch12equi.inp, untitled2.inp
   call renumerotate(NTestGitterR,GitterR)
  
  print *, 'Terminated createRoundMesh2 and createRoundRand: OK'
  
   ! --------------------< functional test of verifymesh module subroutines >-!
   ne=3

 ! copy Nodes coordinates in three vectors ux,uy,uz, ith-component owns ith-node 
 !  do i=1,3
 !    select case (i)
 !       case (1)
 !        vx%vects(1)= 0.0
 !        vx%vects(2)= 0.9659
 !        vx%vects(3)= 0.25882
 !      case (2)
 !        vy%vects(1)= 0.0
 !        vy%vects(2)= -0.19827
 !        vy%vects(3)= 0.7400
 !      case (3)
 !        vz%vects(1)= 0.0
 !        vz%vects(2)= 0.16637
 !        vz%vects(3)= 0.62089
 !   end select
 ! enddo
 ! deallocate(vx%vects, vy%vects, vz%vects, STAT=allocStat) 
 ! if (allocStat.NE.0) then
 !    print *, 'ERROR Deallocate vx vy vz !'
 !    STOP
 ! end if

  ! ------------------------------------------< Single element test calculation >---!
  allocate(C_evalmembrane%Coefficients(1:3,1:6), STAT=allocStat) 
  if (allocStat.NE.0) then
     print *, 'ERROR Allocate K_eval: Could not allocate all variables!'
     STOP
  end if
  ni=1
  call linlgcreatematrixMem(NTestGitterR,GitterR,ni,C_evalmembrane,Pmat_xyz)  ! membrane Teil Matrix Kurvature
  call linlgintegratematrixMem(K_membrane,x,y,C_evalmembrane,Pmat_xyz) ! membrane Steifigkeit Matrix

  call linprntmatrixMem(NTestGitterR,K_membrane) 

  deallocate(C_evalmembrane%Coefficients, STAT=allocStat)
 
  ! ------------------------------------------< Speicherplatz de-allokieren >---!

  call deallocateFieldsRound(NTestGitterR,GitterR,ExaktLR) 
  call deallocateRenum
  print *, 'Deallocate Gitter: OK'
  call deallocateFieldsBCSRound(NTestGitterR,RBU,RBS) 
  print *, 'Deallocate Randbedingungen: OK'
  call deallocateFieldsVarRound(VarFelderR)
  print *, 'Deallocate Numeric Variables: OK'
 ! call deallocateFieldsCoeffRound(GitterCoeffR)
  print *, 'Deallocate Coefficients: OK'

    write(*,*) 
    write(*,*) '============= Program terminated correctly ================ '

 103  format (a20)
 403  format (a, i10, i10, f10.5, i10, f10.5, i10, f10.5)

end program main
