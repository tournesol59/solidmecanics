program mainTruss
  use typesbalken
  use allocateTruss_mod
  use InputTruss_mod
  use CreateTruss_mod
  use RigidTrussContacts_mod
  use printabstractstruct_mod
  use CreateBeamTriangles_mod
  use ForcesLimits_mod

  implicit none

#define BEAM_ELMT_SIZE  6
#define NUM_DOFREEDOM_MAX  3
! <-------------------- INTERFACE --------------------->
!  interface 
!    subroutine Input_file(MeshI,MeshT,MeshG,meshpattern,Dimposed,Fimposed)
!    use typesbalken
!    type(tGitterInfo)            :: MeshI     ! Gitterwerte                   !
!    type(tMeshCoord)           :: MeshT     ! Gitterwerte                      !
!    type(tMeshElmt),pointer    :: MeshG(:)
!    integer,pointer            :: meshpattern(:,:) ! Definition von Elementen !
!    real,pointer               :: Dimposed(:,:)
!    real,pointer               :: Fimposed(:,:)
!    end subroutine Input_file
!  end interface 
! <-------------------- END INTERFACE ------------------>
  integer               :: ProblemSelect
!<----------------- MECH. STRUCTURE ELMTs DECLARATION -------->
!  integer               :: ne=3 ! dof Freiheitsgraden Anzahl
  integer               :: sizeDun=6 ! 3*2 (c1 d1 f1 c2 d2 f2
  type(tMeshInfo)       :: GitterInfo
  type(tMesh2DInfo)     :: Gitter2DInfo
  type(tMeshCoord)      :: MeshPunkte
!  type(tMeshElmt)       :: MeshT1, MeshT2
!  type(tMeshElmt),dimension(:),allocatable :: MeshGen(:)   ! Elementen 1-2 und 1-3 THIS CANNOT BE PASSSED AS POINTER
!  type(tVarElmt),dimension(:),allocatable   :: VarGen(:)    ! Werten (U,F) in Elementen
  type(tMeshElmt),pointer :: MeshGen(:)   ! fuer nt Elementen 1-2 und 1-3 THIS WORKS
  type(tMeshBeam),pointer :: Meshbalk(:) ! fuer nbeam Beam Elementen
  type(tVarElmt),pointer :: VarGen(:)    ! Werten (U,F) in Elementen
  integer,pointer       :: elmtabbdung(:,:) 
  type(tMesh2DSection)    :: Gitter2D
  type(tVar2DSection)     :: Werte2D

!<----------------- MECH. STRUCTURE PROBLEM DEFINITION -------->
  type(tRigidFullMat)   :: Ke_H   ! Whole global rigidity matrix - will not be used
  type(tRigidFullMat)   :: A_H  ! Sub matrices from Ke_H (f=A.U,unknown + B.U,applied)
  type(tRigidFullMat)   :: Ke_Be   ! whole Timoshenko elements rigidity matrix
  type(tVarFull)        :: Var_H, B_H, Var_Be  ! grouped solution in Ue,Fe vectors and also for Timoshenko element
  type(tExakt)          :: Exakt ! Struct fuer Exakte Loesung
  integer,pointer       :: Dunbekannt(:,:)  ! Unbekannte Bewegungen zu berechnen
  real,pointer          :: Dimponiert(:,:)  ! imponierte Bewegungen
  integer,pointer       :: Kraftbewgg(:,:)  ! Result of Calc: unbekannte Kraefte in bewegungsfrei Lenkung
  real,pointer          :: Kraftimponiert(:,:)  ! Bekannte imponiert Kraefte
  integer,pointer       :: IndizeKe(:)

!<----------------- INTEGERS and MISCs ----------------------->
  integer               :: jj,kk,km1,kd
  integer               :: i,j,k,l,p,allocStat

  integer,pointer       :: inull(:)
  character(len=30)     :: SubName="mainTruss                     "

!<------------------- BEGIN COMMON PARTS TO TRUSS AND SECTION -------------->
!  !********************* Enter Manually sections: Beam that are modeles by segment of doublets ***********!

  write(*,*) "Enter problems selection: truss(1) or section(2)"
  read(*,311) ProblemSelect

! ************** read only first line of input def file, the rest is parsed by Input_file procedure
  OPEN(UNIT=25, FILE='Untitled3.in', ACTION='READ')
  read(25,*) ! #
  read(25,705) GitterInfo%nn, GitterInfo%nt, GitterInfo%nbeam, GitterInfo%ngeobeam, GitterInfo%nsection 
  CLOSE(UNIT=25)


!----------------- PART 1 ALLOCATION FOR TRUSS ----------------------->
! ******************* allocate Geometry Strukturen ***************
  call allocateGridElmt(GitterInfo,MeshPunkte,elmtabbdung,MeshGen)

  allocate(inull(1), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate null"
  endif
  
  !--- Beam Section for 2D calculation (see PART 2) must be placed here
  ! because InputTruss procedure needs this be allocated
  allocate(Meshbalk(1:GitterInfo%nbeam), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate Meshbalk null"
  endif
  do k=1,GitterInfo%nbeam
    allocate( Meshbalk(k)%nodes(10), STAT=allocStat)  ! This is not a good choice of Nmax 10 nodes in code but how to know the number before diving in input file?
    if (allocStat.ne.0) then
       print *," Error allocation node array of Meshbalk !"
    endif
  enddo
 
! ******** allocate Elementen Strukturen (allocatable) *********
!--- structure data for 1D calculation
  call allocateTrussDef(GitterInfo,Dunbekannt,Kraftbewgg,Dimponiert, &
                       Kraftimponiert)
  allocate(VarGen(1:GitterInfo%nt), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate VarGen null"
  endif

  
! ********* allocate definition arrays in the main program **************
  allocate(Ke_H%Ke(1: NUM_DOFREEDOM_MAX *GitterInfo%nn,1: NUM_DOFREEDOM_MAX *GitterInfo%nn), Var_H%Ue(1: NUM_DOFREEDOM_MAX *GitterInfo%nn), Var_H%Fe(1: NUM_DOFREEDOM_MAX *GitterInfo%nn), &
           Exakt%Ux(1: NUM_DOFREEDOM_MAX *GitterInfo%nn), Exakt%Fx(1: NUM_DOFREEDOM_MAX *GitterInfo%nn), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation global matrix !"
  endif  
  do i=1, NUM_DOFREEDOM_MAX *GitterInfo%nn       ! initialisierung
    do j=1, NUM_DOFREEDOM_MAX *GitterInfo%nn
       Ke_H%Ke(i,j)=0.0
    enddo
    Var_H%Ue(i)=0.0
    Var_H%Fe(i)=0.0   
  enddo
! <-------------------- 2nd TIME: PARSE INPUT FILE FOR DEFINITION of points and imposed bc--------------->

  call Input_file(GitterInfo, MeshPunkte, MeshGen, Meshbalk, elmtabbdung, Dimponiert, Kraftimponiert)

  call InputSetUnknown(GitterInfo, MeshPunkte, elmtabbdung, Dunbekannt, & 
                       Dimponiert, Kraftbewgg, Kraftimponiert)

  OPEN(UNIT=27, FILE='Untitled3.log', ACTION='WRITE')
  write(27,701) "----inputTruss"
  
  do i=1,GitterInfo%nt
  write(*,703) elmtabbdung(i,1), elmtabbdung(i,2), elmtabbdung(i,3)
  enddo

  !verify and print output 
  call prntinputmatrices(GitterInfo%nn,GitterInfo%nt,elmtabbdung,Dunbekannt,Dimponiert,Kraftbewgg,Kraftimponiert,Subname)
 
! ****************************** Definitions von Punkten und Elementen Typen ***********************************
 
  call createTrussGen(GitterInfo,MeshGen,MeshPunkte,elmtabbdung)

  call printCreateTrussGen(GitterInfo%nn,GitterInfo%nt,MeshGen,MeshPunkte,elmtabbdung,Subname)

!************** First Part: Resolve the Laplacian for torsion stress  *****************************
! <-------------------- DO SOMETHING --------------------->
  write(27,701) "----- main2DSection"
! call Matrice_Ke_H_Torsion(Gitter2DInfo, Gitter2D, Ke_H_2D)

!************** Second (Main) Part: Build and Assembly Truss matrices                     *****************************
! <-------------------- DO SOMETHING --------------------->
  write(27,701) "----- mainTruss"
!  call Matrice_Ke_H_Truss(MeshGen(1), VarGen(1), Ke_H, 1,1,2, GitterInfo%nn, GitterInfo%nt, &
    !      elmtabbdung, Dunbekannt, Dimponiert, Kraftbewgg, Kraftimponiert)

!! THINK ABOUT TO GENERALIZE EXPRESSION BELOW: IT COULD HAVE ALSO A SPC condition

    !************* CALL ELEMENT BUILT MATRIZE **************
 
do km1=1,GitterInfo%nt
 do jj=1,GitterInfo%nn   ! it has been chosen to loop over the nodes and test if
 ! there exists an element which connect the two. a priori it works
  do kk=1,GitterInfo%nn 
   if ((jj.eq.elmtabbdung(km1,2)).and.(kk.eq.elmtabbdung(km1,5))) then
      write(27,702) "----- mainTruss: elmt typ= ", MeshGen(km1)%typ
      call Matrice_Ke_H_Truss(MeshGen(km1),VarGen(km1),Ke_H,jj,kk, & 
              GitterInfo%nn, GitterInfo%nt, km1 , &
              elmtabbdung, Dunbekannt, Dimponiert, Kraftbewgg, &
              Kraftimponiert)

      call prntadjuststruct(9, 9, Ke_H, SubName)
   endif
  enddo !! do kk
 enddo !! do jj
 enddo !! km1

! Extract real unknowns submatrix A_H and K_H, which are also allocated A_H%Ke ..
  call ExtractTrussUnknowns(GitterInfo, elmtabbdung, Dunbekannt, Kraftbewgg, Ke_H, Var_H, A_H, B_H, IndizeKe)

  ! copy a submatrix A=(K_PP,K_PL; K_LP, K_LL) from the full matrix Ke_H (L=imposed, S=null, P=unkowns)
  ! thus here for displacements: P={1:1, 1:2} L={} S={1:3, 2:1, 2:2, 2:3, 3:1, 3:2, 3:3}
  ! and forces: P={2:1, 2:2, 2:3 3:1, 3:2} L={1:2} S={1:1 1:3 3:3}

  ! Verify the extraction procedure
  write(27, 701) "----- mainTruss : original Ke_H: "
  call prntadjuststruct(GitterInfo%nn, GitterInfo%nn, Ke_H, SubName)
  write(27, 701) "----- mainTruss : extracted A_H: "
  call prntadjuststruct(GitterInfo%la, GitterInfo%la, A_H, SubName)
 


!****************** DEALLOCATE PART1  ******************

  call deallocateGridElmt(MeshPunkte,elmtabbdung,MeshGen)

  call deallocateTrussDef(Dunbekannt,Kraftbewgg,Dimponiert,Kraftimponiert)

  deallocate(Ke_H%Ke, Var_H%Ue, Var_H%Fe, Exakt%Ux, Exakt%Fx, &
             A_H%Ke, B_H%Ue, B_H%Fe, IndizeKe, STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error deallocation Matrizen !"
  endif  

  deallocate(VarGen, STAT=allocstat)
  if (allocstat.ne.0) then
    print *," Error deallocate VarGen null"
  endif

  
!-------------------------------------------------------------------------!
! PART 2: BEAM SECTION PROBLEM  PARSE INPUT FILE FOR 
!   DEFINITION of section 2S mesh"
!-------------------------------------------------------------------------!
 
 ! user interface
  write(*,*) 'Enter index of the Beam'
  read(*,311) Meshbalk(1)%section_index

  write(*,*) 'Enter number of elements in the Beam'
  read(*,311) Meshbalk(1)%ibnn

  write(*,*) 'Enter surface of the section of the Beam'
  read(*,205) Meshbalk(1)%SArea

  write(*,*) 'Enter moment of inertia of the Beam'
  read(*,205) Meshbalk(1)%CI

  write(*,*) 'Enter Young Modulus of the Beam'
  read(*,206) Meshbalk(1)%EY
 
  write(*,*) 'Enter Poisson number of the Beam'
  read(*,206) Meshbalk(1)%nu
 
  !  assume there is only one beam structure Meshbalk
  allocate(Ke_Be%Ke(1:((Meshbalk(1)%ibnn-1)*BEAM_ELMT_SIZE),1:((Meshbalk(1)%ibnn-1)*BEAM_ELMT_SIZE)), &
           Var_Be%Ue(1:((Meshbalk(1)%ibnn-1)*BEAM_ELMT_SIZE)), Var_Be%Fe(1:((Meshbalk(1)%ibnn-1)*BEAM_ELMT_SIZE)), &
             STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation Beam rigidity matrix !"
  endif  
  do i=1,(Meshbalk(1)%ibnn-1)*BEAM_ELMT_SIZE       ! initialisierung
    do j=1,(Meshbalk(1)%ibnn-1)*BEAM_ELMT_SIZE
      ! write(*,702) "-----mainTruss (i*ibnn*6)+j = ", i*(Meshbalk(1)%ibnn-1)*BEAM_ELMT_SIZE+j
       Ke_Be%Ke(i,j)=0.0
    enddo
    Var_Be%Ue(i)=0.0
    Var_Be%Fe(i)=0.0   
  enddo

!<---------------------- PARSE INPUT FILE FOR DEFINITION of section 2S mesh" ---------------------------->

  call allocMesh2DfromFile(Gitter2DInfo,Gitter2D)
  
  call createMesh2DfromFile(Gitter2DInfo,Gitter2D)

  call deallocMesh2DfromFile(Gitter2DInfo,Gitter2D)
 
!!  call createTrussGen(GitterInfo,MeshGen,MeshPunkte,elmtabbdung)

! ************ First Part: Resolve the Laplacian for torsion stress  *********
! <-------------------- DO SOMETHING --------------------->
  write(27,701) "----- main2DSection"
! call Matrice_Ke_H_Torsion(Gitter2DInfo, Gitter2D, Ke_H_2D)
!************** Second (Main) Part: Build and Assembly Truss matrices                     *****************************

  deallocate(Ke_Be%Ke, Var_Be%Ue, Var_Be%Fe, &
             STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error deallocation Matrizen !"
  endif

  do k=1,GitterInfo%nbeam
    deallocate(Meshbalk(k)%nodes, STAT=allocstat)
    if (allocstat.ne.0) then
      print *," Error deallocate Meshbalk(1)%nodes null"
    endif
  enddo
  deallocate(Meshbalk, STAT=allocstat)
  if (allocstat.ne.0) then
    print *," Error deallocate Meshbalk null"
  endif
! close log file
  CLOSE(UNIT=27)
  print *,"========== PROGRAM TERMINATED CORRECTLY =========="

205  format (e10.2)             
206  format (2e10.2)
311  format (i5)
701 format(a)
702 format(a,i5)
703 format(3i10)
705 format (5i5)

end program mainTruss

