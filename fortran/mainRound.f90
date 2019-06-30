program main

  use typesround
!  use inputround_mod
  use allocatefieldsround_mod
  use createmeshround_mod

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  type(tRoundMeshInfo)            :: NTestGitterR ! Gitterwerte in matriz           !
  type(tRoundMesh)                :: GitterR     ! Gitterwerte in matriz           !
  type(tRoundNumeric)             :: VarFelderR     ! Numerische Felder mit Loesungen !
  type(tRoundExakt)               :: ExaktLR      ! Exakte Loesung                  !
!  type(tRoundFileIO)              :: FileIOR     ! Ausgabesteuerung                !
  type(tRoundRandU)               :: RBU          ! Randbedingungen Biegung
  type(tRoundRandS)               :: RBS          ! Randbedingungen Force
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

   ! --------------------------< Gitter festlegen >---------------------!
  call createRoundMesh2(NTestGitterR,GitterR)    
  print *, 'Terminated createRoundMesh2: OK'
  ! ------------------------------------------< Speicherplatz de-allokieren >---!

  call deallocateFieldsRound(NTestGitterR,GitterR,ExaktLR) 
  call deallocateFieldsBCSRound(NTestGitterR,RBU,RBS) 
  call deallocateFieldsVarRound(VarFelderR)

    write(*,*) 
    write(*,*) '============= Program terminated correctly ================ '


end program main
