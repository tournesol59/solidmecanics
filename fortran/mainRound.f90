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
 ! type(tNode),dimension(1:3)      :: vects        ! contains 3x 3 Vectors
  type(tNode)                     :: vx,vy,vz     ! Zum test of check in verifymeshRound!
  real,dimension(1:3)             :: dist, angle, triplet
  real,dimension(1:4)             :: area
  !--------------------------------------------------------------------------!  
   integer        :: i,j,k,l,m,n,p,ne, allocStat

  ! -------------------------------------------------------------< BEGIN >---!
 ! call input(Gitter,Exakt,Const,FileIO)

  ! ------------------------------------------< Speicherplatz allokieren >---!
  NTestGitterR%nNode = 9
  NTestGitterR%nElem = 8
  NTestGitterR%nxmx = 9 ! frontier length should be optimized in the future ..
  NTestGitterR%nymx = 9 
  allocate(vx%vects(3), STAT=allocStat)
  allocate(vy%vects(3), STAT=allocStat)
  allocate(vz%vects(3), STAT=allocStat) 
  if (allocStat.NE.0) then
     print *, 'ERROR Allocate vx vy vz: Could not allocate all variables!'
     STOP
  end if
  call allocateFieldsRound(NTestGitterR,GitterR,ExaktLR)
  call allocateFieldsBCSRound(NTestGitterR,RBU,RBS)
 ! call allocateFieldsVarRound(NTestGitterR,VarFelderR)
 ! call allocateFieldsCoeffRound(NTestGitterR,GitterCoeffR) 

   ! --------------------------< Gitter festlegen >---------------------!
   call createRoundMesh2(NTestGitterR,GitterR)    
   call createRoundRand2(NTestGitterR,RBU,RBS)
  print *, 'Terminated createRoundMesh2 and createRoundRand: OK'

   ! --------------------< functional test of verifymesh module subroutines >-!
   ne=3
   j=GitterR%elems(ne)%numeros(2)
   k=GitterR%elems(ne)%numeros(3)
   l=GitterR%elems(ne)%numeros(4) 

 ! copy Nodes coordinates in three vectors ux,uy,uz, ith-component owns ith-node 
   do i=1,3
     select case (i)
       case (1)
         vx%vects(1)= GitterR%nodes(j)%vects(1)
         vx%vects(2)= GitterR%nodes(k)%vects(1)
         vx%vects(3)= GitterR%nodes(l)%vects(1)
       case (2)
         vy%vects(1)= GitterR%nodes(j)%vects(2)
         vy%vects(2)= GitterR%nodes(k)%vects(2)
         vy%vects(3)= GitterR%nodes(l)%vects(2)
       case (3)
         vz%vects(1)= GitterR%nodes(j)%vects(3)
         vz%vects(2)= GitterR%nodes(k)%vects(3)
         vz%vects(3)= GitterR%nodes(l)%vects(3)

    end select
  enddo
   write(*, 403) "  initial value in vx: ", ne, j, vx%vects(1), k, vx%vects(2), l, vx%vects(3)
   print *, 'Terminated test verifymesh part 1: OK'

  call trigcharacteristics(vx,vy,vz, dist, angle, triplet, area)
  print *, 'Terminated test verifymesh part 2: OK'

  ! ------------------------------------------< Speicherplatz de-allokieren >---!

  call deallocateFieldsRound(NTestGitterR,GitterR,ExaktLR) 
  print *, 'Deallocate Gitter: OK'
  call deallocateFieldsBCSRound(NTestGitterR,RBU,RBS) 
  print *, 'Deallocate Randbedingungen: OK'
 ! call deallocateFieldsVarRound(VarFelderR)
  print *, 'Deallocate Numeric Variables: OK'
 ! call deallocateFieldsCoeffRound(GitterCoeffR)
  print *, 'Deallocate Coefficients: OK'

    write(*,*) 
    write(*,*) '============= Program terminated correctly ================ '

 103  format (a20)
 403  format (a, i10, i10, f10.5, i10, f10.5, i10, f10.5)

end program main
