program main

  use typesround
!  use inputround_mod
  use allocatefieldsround_mod
  use createmeshround_mod
  use verifymeshround_mod
  use spacegeometry_mod

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
  type(tPassageMatrix)            :: axref_xyz, axout_xyz    ! a 3x3 Matrix contains local axes 
                                                  ! vs global axes coordinates!
  type(tNode),dimension(1:3)      :: vects        ! contains 3x 3 Vectors
  !--------------------------------------------------------------------------!  
   integer        :: i,j,n, allocStat

  ! -------------------------------------------------------------< BEGIN >---!
 ! call input(Gitter,Exakt,Const,FileIO)

  ! ------------------------------------------< Speicherplatz allokieren >---!
  NTestGitterR%nNode = 9
  NTestGitterR%nElem = 8
  NTestGitterR%nxmx = 3
  NTestGitterR%nymx = 3
  call allocateFieldsRound(NTestGitterR,GitterR,ExaktLR)
  call allocateFieldsBCSRound(NTestGitterR,RBU,RBS)
  call allocateFieldsVarRound(NTestGitterR,VarFelderR)
  call allocateFieldsCoeffRound(NTestGitterR,GitterCoeffR) 

   ! --------------------------< Gitter festlegen >---------------------!
  call createRoundMesh2(NTestGitterR,GitterR)    
  print *, 'Terminated createRoundMesh2: OK'

   ! --------------------< functional test of verifymesh module subroutines >-!
  call trigrotationxcheck(1, GitterR, axref_xyz, axout_xyz, vects)
  print *, 'Terminated test verifymesh: OK'

  ! ------------------------------------------< Speicherplatz de-allokieren >---!

  call deallocateFieldsRound(NTestGitterR,GitterR,ExaktLR) 
  call deallocateFieldsBCSRound(NTestGitterR,RBU,RBS) 
  call deallocateFieldsVarRound(VarFelderR)
  call deallocateFieldsCoeffRound(GitterCoeffR)

    write(*,*) 
    write(*,*) '============= Program terminated correctly ================ '


end program main
