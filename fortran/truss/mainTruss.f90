program mainTruss
  use typesbalken
  use InputTruss_mod
  use CreateTruss_mod
  use RigidTrussContacts_mod
  use printabstractstruct_mod
  use CreateBeamTriangles_mod

  implicit none

#define BEAM_ELMT_SIZE  6

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

!<----------------- MECH. STRUCTURE ELMTs DECLARATION -------->
  integer               :: ne=3 ! dof Freiheitsgraden Anzahl
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
  type(tRigidFullMat)   :: A_H, B_H  ! Sub matrices from Ke_H (f=A.U,unknown + B.U,applied)
  type(tRigidFullMat)   :: Ke_Be   ! whole Timoshenko elements rigidity matrix
  type(tVarFull)        :: Var_H, Var_Be  ! grouped solution in Ue,Fe vectors and also for Timoshenko element
  type(tExakt)          :: Exakt ! Struct fuer Exakte Loesung
  integer,pointer       :: Dunbekannt(:,:)  ! Unbekannte Bewegungen zu berechnen
  real,pointer          :: Dimponiert(:,:)  ! imponierte Bewegungen
  integer,pointer       :: Kraftbewgg(:,:)  ! Result of Calc: unbekannte Kraefte in bewegungsfrei Lenkung
  real,pointer          :: Kraftimponiert(:,:)  ! Bekannte imponiert Kraefte

!<----------------- INTEGERS and MISCs ----------------------->
  integer               :: jj,kk
  integer               :: i,j,k,km1,l,p,allocStat

  integer,pointer       :: inull(:)
  character(len=30)     :: SubName="mainTruss                     "

!<--------------------------------- BEGIN ---------------------------->

! ************** read only first line of input def file, the rest is parsed by Input_file procedure
  OPEN(UNIT=25, FILE='Untitled3.in', ACTION='READ')
  read(25,*) ! #
  read(25,705) GitterInfo%nn, GitterInfo%nt, GitterInfo%nbeam, GitterInfo%ngeobeam, GitterInfo%nsection 
  CLOSE(UNIT=25)


! **************************** allocate Geometry Strukturen ************************
  allocate(inull(1), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate null"
  endif
!  allocate(M_rot(9), STAT=allocstat)
!  if (allocstat.ne.0) then
!    print *,"mainTruss: error allocate M_rot"
!  endif
  allocate(MeshPunkte%x(1:GitterInfo%nn), MeshPunkte%y(1:GitterInfo%nn), &
            MeshPunkte%z(1:GitterInfo%nn), &
            MeshPunkte%elmts(1:GitterInfo%nt,1:7), &
            STAT=allocStat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate MeshPunkte"
  endif

  allocate(elmtabbdung(1:GitterInfo%nt,1:3), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation connection matrix !"
  endif

! ***************************** allocate Elementen Strukturen (allocatable) *****************************************
  allocate(Dunbekannt(1:GitterInfo%nt,1:6), Dimponiert(1:GitterInfo%nn,1:3), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation Bewegungs definitionsmatrizen !"
  endif 
  allocate(Kraftbewgg(1:GitterInfo%nt,1:6), Kraftimponiert(1:GitterInfo%nn,1:3), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation Matrizen zur berechnenden Bewegungskraefte matrix !"
  endif 
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
  allocate(MeshGen(1:GitterInfo%nt), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate MeshGen null"
  endif
! Coeffs1 not used
!  do k=1,GitterInfo%nt
!    allocate( MeshGen(k)%CoeffsH1(1:4), MeshGen(k)%CoeffsH2(1:4), &
!                 MeshGen(k)%CoeffsH3(1:4), MeshGen(k)%CoeffsH4(1:4), STAT=allocStat)
!    if (allocStat.ne.0) then
!       print *," Error allocation polynoms elmt 2 !"
!    endif
!  enddo
  allocate(VarGen(1:GitterInfo%nt), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate VarGen null"
  endif
  !********************* Enter Manually sections: Beam that are modeles by segment of doublets ***********!
  write(*,*) 'Enter index of the Beam'
  read(*,311) Meshbalk%section_index

  write(*,*) 'Enter number of elements in the Beam'
  read(*,311) Meshbalk%ibnn

  write(*,*) 'Enter surface of the section of the Beam'
  read(*,205) Meshbalk%SArea

  write(*,*) 'Enter moment of inertia of the Beam'
  read(*,205) Meshbalk%CI

  write(*,*) 'Enter Young Modulus of the Beam'
  read(*,206) Meshbalk%EY
 
  write(*,*) 'Enter Poisson number of the Beam'
  read(*,206) Meshbalk%vu
  
! ***************************** allocate definition arrays in the main program *********************************
  allocate(Ke_H%Ke(1:ne*GitterInfo%nn,1:ne*GitterInfo%nn), Var_H%Ue(1:ne*GitterInfo%nn), Var_H%Fe(1:ne*GitterInfo%nn), &
           Exakt%Ux(1:ne*GitterInfo%nn), Exakt%Fx(1:ne*GitterInfo%nn), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation global matrix !"
  endif  
  do i=1,ne*GitterInfo%nn       ! initialisierung
    do j=1,ne*GitterInfo%nn
       Ke_H%Ke(i,j)=0.0
    enddo
    Var_H%Ue(i)=0.0
    Var_H%Fe(i)=0.0   
  enddo
!  assume there is only one beam structure Meshbalk
  allocate(Ke_Be%Ke(1:((Meshbalk(1)%ibnn-1)*BEAM_ELMT_SIZE),1:((Meshbalk(1)%ibnn-1)*BEAM_ELMT_SIZE)), &
           Var_Be%Ue(1:((Meshbalk(1)%ibnn-1)*BEAM_ELMT_SIZE)), Var_Be%Fe(1:((Meshbalk(1)%ibnn-1)*BEAM_ELMT_SIZE)), &
             STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation Beam rigidity matrix !"
  endif  
  do i=1,(Meshbalk(1)%ibnn-1)*BEAM_ELMT_SIZE       ! initialisierung
    do j=1,(Meshbalk(1)%ibnn-1)*BEAM_ELMT_SIZE
       Ke_Be%Ke(i,j)=0.0
    enddo
    Var_Be%Ue(i)=0.0
    Var_Be%Fe(i)=0.0   
  enddo

! <-------------------- PARSE INPUT FILE FOR DEFINITION of points and imposed bc--------------->
  call Input_file(GitterInfo, MeshPunkte, MeshGen, Meshbalk, elmtabbdung, Dimponiert, Kraftimponiert)

  call InputSetUnknown(GitterInfo, MeshPunkte, elmtabbdung, Dunbekannt, & 
                       Dimponiert, Kraftbewgg, Kraftimponiert)

! ****************************** Definitions von Punkten und Elementen Typen ***********************************
  call createTrussGen(GitterInfo,MeshGen,MeshPunkte,elmtabbdung)

  OPEN(UNIT=27, FILE='Untitled3.log', ACTION='WRITE')
  write(27,701) "----inputTruss"

  call printCreateTrussGen(GitterInfo%nn,GitterInfo%nt,MeshGen,MeshPunkte,elmtabbdung)

!<---------------------- PARSE INPUT FILE FOR DEFINITION of section 2S mesh" ---------------------------->
  call allocMesh2DfromFile(Gitter2DInfo,Gitter2D)

  call createMesh2DfromFile(Gitter2DInfo,Gitter2D)

  call deallocMesh2DfromFile(Gitter2DInfo,Gitter2D)
 
  call createTrussGen(GitterInfo,MeshGen,MeshPunkte,elmtabbdung)

!************** First Part: Resolve the Laplacian for torsion stress  *****************************
! <-------------------- DO SOMETHING --------------------->
  write(27,701) "----- main2DSection"
! call Matrice_Ke_H_Torsion(Gitter2DInfo, Gitter2D, Ke_H_2D)
!************** Second (Main) Part: Build and Assembly Truss matrices                     *****************************
! <-------------------- DO SOMETHING --------------------->
  write(27,701) "----- mainTruss"
!      call Matrice_Ke_H_Truss(MeshGen(1), VarGen(1), Ke_H, 1,2, GitterInfo%nn, GitterInfo%ne, &
!                              elmtabbdung, &
!                              Dunbekannt, Dimponiert, &
!                              Kraftbewgg, Kraftimponiert)

!! THINK ABOUT TO GENERALIZE EXPRESSION BELOW: IT COULD HAVE ALSO A SPC condition

    !************* CALL ELEMENT BUILT MATRIZE **************
! km1=0
! do jj=1,nn
!  do kk=1,nn
!   if (kk.ne.km1) then
!      write(*,702) "----- mainTruss: elmt typ= ", MeshGen(km1)%typ
!      call Matrice_Ke_H_Truss(MeshGen(km1), VarGen(km1), Ke_H, jj,kk, n, ne, &
!                              elmtabbdung, &
!                              Dunbekannt, Dimponiert, &
!                              Kraftbewgg, Kraftimponiert)
!
!     call prntadjuststruct(9, 9, Ke_H, SubName)
!     km1=km1+1
!   endif
! enddo !! do kk
! enddo !! do jj

  ! copy a submatrix A=(K_PP,K_PL; K_LP, K_LL) from the full matrix Ke_H (L=imposed, S=null, P=unkowns)
  ! thus here for displacements: P={1:1, 1:2} L={} S={1:3, 2:1, 2:2, 2:3, 3:1, 3:2, 3:3}
  ! and forces: P={2:1, 2:2, 2:3 3:1, 3:2} L={1:2} S={1:1 1:3 3:3}


!****************** DEALLOCATE ******************
  deallocate(Ke_H%Ke, Var_H%Ue, Var_H%Fe, Exakt%Ux, Exakt%Fx, &
             elmtabbdung, STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error deallocation Matrizen !"
  endif  

  deallocate(Ke_Be%Ke, Var_Be%Ue, Var_Be%Fe, &
             STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error deallocation Matrizen !"
  endif 

  deallocate(Dunbekannt, Dimponiert, STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error deallocation Bewegungs definitionsmatrizen !"
  endif 
  deallocate(Kraftbewgg, Kraftimponiert, STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error deallocation Matrizen zur berechnenden Bewegungskraefte matrix !"
  endif 
!  deallocate(inull, M_rot, STAT = allocStat)
!  if (allocStat.ne.0) then
!     print *," Error deallocation Matrizen inull... !"
!  endif 
! not used any more:
!  do k=1,2
!    deallocate( MeshGen(k)%CoeffsH1, MeshGen(k)%CoeffsH2, &
!                MeshGen(k)%CoeffsH3, MeshGen(k)%CoeffsH4, STAT = allocStat)
!    if (allocStat.ne.0) then
!       write(*,702) " Error deallocation polynoms elmt k= !", k
!    endif
!  enddo
  deallocate(MeshGen, STAT=allocstat)
  if (allocstat.ne.0) then
    print *," Error allocate MeshGen null"
  endif

  deallocate(VarGen, STAT=allocstat)
  if (allocstat.ne.0) then
    print *," Error deallocate VarGen null"
  endif
  deallocate(MeshPunkte%x,MeshPunkte%y,MeshPunkte%z, STAT=allocstat)
  if (allocstat.ne.0) then
    print *," Error deallocate MeshPunkte null"
  endif
  deallocate(Meshbalk(1)%nodes, STAT=allocstat)
  if (allocstat.ne.0) then
    print *," Error deallocate Meshbalk(1)%nodes null"
  endif
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

705 format (5i5)

end program mainTruss
