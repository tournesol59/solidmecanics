program mainTruss
  use typesbalken
  use InputTruss_mod
  use RigidTrussContacts_mod
  use printabstractstruct_mod

  implicit none

  integer               :: ne=3 ! dof Freiheitsgraden Anzahl
  real                  :: anglescalar
  type(tMeshInfo)       :: MeshInfo
  type(tMeshCoord)      :: MeshPunkte
  type(tMeshElmt)       :: MeshT1, MeshT2
  type(tMeshElmt),dimension(:),allocatable  :: MeshGen(:)   ! Elementen 1-2 und 1-3
  type(tVarElmt),dimension(:),allocatable   :: VarGen(:)    ! Werten (U,F) in Elementen
  integer,pointer       :: elmtabbdung(:,:) 
  integer,pointer       :: meshverbindg(:,:) ! Definition von Verbindungen !
  real,pointer          :: M_rot(:)
  type(tRigidFullMat)   :: Ke_H   ! Whole global rigidity matrix
  type(tRigidFullMat)   :: A_H, B_H  ! Sub matrices from Ke_H (f=A.U,unknown + B.U,applied
  type(tVarFull)        :: Var_H  ! grouped solution in Ue,Fe vectors
  type(tExakt)          :: Exakt ! Struct fuer Exakte Loesung
  integer,pointer       :: Dunbekannt(:,:)  ! Unbekannte Bewegungen zu berechnen
  real,pointer          :: Dimponiert(:,:)  ! imponierte Bewegungen
  integer,pointer       :: Kraftbewgg(:,:)  ! Result of Calc: unbekannte Kraefte in bewegungsfrei Lenkung
  real,pointer          :: Kraftimponiert(:,:)  ! Bekannte imponiert Kraefte
  integer               :: jj,kk
  integer               :: i,j,k,km1,l,p,allocStat

  integer,pointer       :: inull(:)
  character(len=30)     :: SubName="mainTruss                     "

! ************** read only first line of input def file, the rest is parsed by Input_file procedure
  OPEN(UNIT=25, FILE='Untitled2.msh', ACTION='READ')
  read(25,705) MeshInfo%nn, MeshInfo%nt
  CLOSE(UNIT=25)

! **************************** allocate Geometry Strukturen ************************
  allocate(inull(1), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate null"
  endif
  allocate(M_rot(9), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate M_rot"
  endif
  allocate(MeshPunkte%x(1:MeshInfo%nn), MeshPukte%y(1:MeshInfo%nn), MeshPunkte%z(1:MeshInfo%nn), STAT=allocStat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate M_rot"
  endif

! ***************************** allocate Elementen Strukturen (allocatable) *****************************************
  allocate(Dunbekannt(1:MeshInfo%nt,1:7), Dimponiert(1:MeshInfo%nn,1:3), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation Bewegungs definitionsmatrizen !"
  endif 
  allocate(Kraftbewgg(1:MeshInfo%nt,1:7), Kraftimponiert(1:MeshInfo%nn,1:3), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation Matrizen zur berechnenden Bewegungskraefte matrix !"
  endif 

  allocate(MeshGen(1:2), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate MeshGen null"
  endif
  do k=1,2
      allocate( MeshGen(k)%CoeffsH1(1:4), MeshGen(k)%CoeffsH2(1:4), &
                 MeshGen(k)%CoeffsH3(1:4), MeshGen(k)%CoeffsH4(1:4), STAT=allocStat)
    if (allocStat.ne.0) then
       print *," Error allocation polynoms elmt 2 !"
    endif
  enddo
  allocate(VarGen(1:2), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate VarGen null"
  endif

! ****************************** Definitions von Punkten und zwei Elementen Typen ***********************************

  MeshT1%typ = 1           ! beam
  MeshT1%nraumx = 2        ! subdivision, default=2
  MeshT1%startx =  MeshPunkte%x(1)        ! node 1 x-koord  !
  MeshT1%endx =  MeshPunkte%x(3)          ! node 2 x-koord  !
  MeshT1%starty =  MeshPunkte%y(1)       ! node 1 y-koord  !
  MeshT1%endy =  MeshPunkte%y(3)          ! node 2 y-koord  !
  MeshT1%startz =  MeshPunkte%z(1)        ! node 2 z-koord  !
  MeshT1%endz =  MeshPunkte%z(3)          ! node 2 z-koord  !

  MeshT1%dlen = sqrt((MeshT1%endx - MeshT1%startx)**2+ &
                    (MeshT1%endy - MeshT1%starty)**2+ &
                    (MeshT1%endz - MeshT1%startz)**2)           ! Laengen des Rechengebiets (Balk) in x/y !
  MeshT1%SArea = 0.0004         ! Section area
  MeshT1%EY = 210000000            ! YOUNG Modulus
  MeshT1%vu = 0.3            ! Poisson Coeff
  MeshT1%CI= 0.005             ! kinetics inertia moment
  MeshT1%q1 = -0.0
  MeshT1%q2 = -0.0         ! verteilte  Lasten q(x)=x/xlen*q1+(1-x/xlen)*q2
!     real        :: q?
  MeshT1%node_1 = 1 ! shall be overwritten
  MeshT1%node_2 = 2
  allocate( MeshT1%CoeffsH1(1:4), MeshT1%CoeffsH2(1:4), MeshT1%CoeffsH3(1:4), &
             MeshT1%CoeffsH4(1:4), STAT=allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation polynoms elmt 1 !"
  endif  

  MeshT2%typ = 2           ! bar
  MeshT2%nraumx = 2        ! subdivision, default=2
  MeshT2%startx = MeshPunkte%x(1)        ! node 1 x-koord  !
  MeshT2%endx = MeshPunkte%x(2)          ! node 2 x-koord  !
  MeshT2%starty = MeshPunkte%y(1)       ! node 1 y-koord  !
  MeshT2%endy = MeshPunkte%y(2)          ! node 2 y-koord  !
  MeshT2%startz = MeshPunkte%z(1)        ! node 2 z-koord  !
  MeshT2%endz = MeshPunkte%z(2)          ! node 2 z-koord  !

 MeshT2%dlen = sqrt((MeshT2%endx - MeshT2%startx)**2+ &
                     (MeshT2%endy - MeshT2%starty)**2+ &
                     (MeshT2%endz - MeshT2%startz)**2)          ! Laengen des Rechengebiets (Balk) in x/y !
  MeshT2%SArea = 0.0004         ! Section area
  MeshT2%EY = 210000000            ! YOUNG Modulus
  MeshT2%vu = 0.3            ! Poisson Coeff
  MeshT2%CI= 0.0             ! kinetics inertia moment
  MeshT2%q1 = -0.25
  MeshT2%q2 = -0.25         ! verteilte  Lasten q(x)=x/xlen*q1+(1-x/xlen)*q2
!     real        :: q?
  MeshT2%node_1 = 1 ! shall be overwritten
  MeshT2%node_2 = 2
  allocate( MeshT2%CoeffsH1(1:4), MeshT2%CoeffsH2(1:4), &
            MeshT2%CoeffsH3(1:4), MeshT2%CoeffsH4(1:4), STAT=allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation polynoms elmt 2 !"
  endif

! ***************************** allocate definition arrays in the main program *********************************
  allocate(Ke_H%Ke(1:ne*MeshInfo%nn,1:ne*MeshInfo%nn), Var_H%Ue(1:ne*MeshInfo%nn), Var_H%Fe(1:ne*MeshInfo%nn), &
           Exakt%Ux(1:ne*MeshInfo%nn), Exakt%Fx(1:ne*MeshInfo%nn), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation global matrix !"
  endif  
  do i=1,ne*MeshInfo%nn       ! initialisierung
    do j=1,ne*MeshInfo%nn
       Ke_H%Ke(i,j)=0.0
    enddo
    Var_H%Ue(i)=0.0
    Var_H%Fe(i)=0.0   
  enddo

  allocate(elmtabbdung(1:MeshInfo%nn,1:MeshInfo%nn),  meshverbindg(1:MeshInfo%nn,1:MeshInfo%nn), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation connection matrix !"
  endif

! <-------------------- PARSE INPUT FILE FOR DEFINITION of points and imposed bc--------------->
  call Input_file(MeshInfo, MeshPunkte, elmtabbdung, meshverbindg, Dunbekannt, Dimponiert, Kraftbewgg, Kraftimponiert)

! <-------------------- DO SOMETHING --------------------->
  write(*,701) "----- mainTruss"
  k=1
  km1=1
  do jj=1,MeshInfo%nn
    do kk=jj,MeshInfo%nn
      if (elmtabbdung(jj,kk).eq.2) then
!! THINK ABOUT TO GENERALIZE EXPRESSION BELOW: IT COULD HAVE ALSO A SPC condition
     MeshGen(k)%node_1 = jj
     MeshGen(k)%node_2 = kk
     MeshGen(k)%typ = elmtabbdung(jj,kk)           ! bar or beam
     MeshGen(k)%nraumx = 2        ! subdivision, default=2
     MeshGen(k)%startx = MeshPunkte%x(jj)        ! node 1 x-koord  !
     MeshGen(k)%endx = MeshPunkte%x(kk)          ! node 2 x-koord  !
     MeshGen(k)%starty = MeshPunkte%y(jj)       ! node 1 y-koord  !
     MeshGen(k)%endy = MeshPunkte%y(kk)          ! node 2 y-koord  !
     MeshGen(k)%startz = MeshPunkte%z(jj)        ! node 2 z-koord  !
     MeshGen(k)%endz = MeshPunkte%z(kk)          ! node 2 z-koord  ! 
   
    MeshGen(k)%dlen = sqrt(( MeshGen(k)%endx -  MeshGen(k)%startx)**2+ &
                           ( MeshGen(k)%endy -  MeshGen(k)%starty)**2+ &
                           ( MeshGen(k)%endz -  MeshGen(k)%startz)**2)          ! Laengen des Rechengebiets (Balk) in x/y !*
  write(*,704) 2,MeshGen(k)%dlen
     MeshGen(k)%SArea = MeshT2%SArea        ! Section area
     MeshGen(k)%EY = MeshT2%EY             ! YOUNG Modulus
     MeshGen(k)%vu = MeshT2%vu             ! Poisson Coeff
     MeshGen(k)%CI=  MeshT2%CI              ! kinetics inertia moment
     MeshGen(k)%q1 = MeshT2%q1 
     MeshGen(k)%q2 = MeshT2%q2         ! verteilte  Lasten q(x)=x/xlen*q1+(1-x/xlen)*q2
!     real        :: q?

     k=k+1 ! important
   elseif (elmtabbdung(jj,kk).eq.1) then
     MeshGen(k)%node_1 = jj
     MeshGen(k)%node_2 = kk
     MeshGen(k)%typ = elmtabbdung(jj,kk)           ! bar or beam
     MeshGen(k)%nraumx = 2        ! subdivision, default=2
     MeshGen(k)%startx = MeshPunkte%x(jj)        ! node 1 x-koord  !
     MeshGen(k)%endx = MeshPunkte%x(kk)          ! node 2 x-koord  !
     MeshGen(k)%starty = MeshPunkte%y(jj)       ! node 1 y-koord  !
     MeshGen(k)%endy = MeshPunkte%y(kk)          ! node 2 y-koord  !
     MeshGen(k)%startz = MeshPunkte%z(jj)        ! node 2 z-koord  !
     MeshGen(k)%endz = MeshPunkte%z(kk)          ! node 2 z-koord  !
   
    MeshGen(k)%dlen = sqrt(( MeshGen(k)%endx -  MeshGen(k)%startx)**2+ &
                           ( MeshGen(k)%endy -  MeshGen(k)%starty)**2+ &
                           ( MeshGen(k)%endz -  MeshGen(k)%startz)**2)          ! Laengen des Rechengebiets (Balk) in x/y !*
  write(*,704) 2,MeshGen(k)%dlen
     MeshGen(k)%SArea = MeshT1%SArea        ! Section area
     MeshGen(k)%EY = MeshT1%EY             ! YOUNG Modulus
     MeshGen(k)%vu = MeshT1%vu             ! Poisson Coeff
     MeshGen(k)%CI=  MeshT1%CI              ! kinetics inertia moment
     MeshGen(k)%q1 = MeshT1%q1 
     MeshGen(k)%q2 = MeshT1%q2         ! verteilte  Lasten q(x)=x/xlen*q1+(1-x/xlen)*q2
!     real        :: q?

     k=k+1 ! important
   endif !! endif typ elmtabbdung
    !************* CALL ELEMENT BUILT MATRIZE **************
   if (k.ne.km1) then
      write(*,702) "----- mainTruss: elmt typ= ", MeshGen(km1)%typ
      call Matrice_Ke_H_Truss(MeshGen(km1), VarGen(km1), Ke_H, jj,kk, n, ne, &
                              elmtabbdung, &
                              Dunbekannt, Dimponiert, &
                              Kraftbewgg, Kraftimponiert)

     call prntadjuststruct(9, 9, Ke_H, SubName)
     km1=km1+1
   endif
 enddo !! do kk
 enddo !! do jj

  ! copy a submatrix A=(K_PP,K_PL; K_LP, K_LL) from the full matrix Ke_H (L=imposed, S=null, P=unkowns)
  ! thus here for displacements: P={1:1, 1:2} L={} S={1:3, 2:1, 2:2, 2:3, 3:1, 3:2, 3:3}
  ! and forces: P={2:1, 2:2, 2:3 3:1, 3:2} L={1:2} S={1:1 1:3 3:3}


!****************** DEALLOCATE ******************
  deallocate(Ke_H%Ke, Var_H%Ue, Var_H%Fe, Exakt%Ux, Exakt%Fx, &
             elmtabbdung, STAT = allocStat)
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
  deallocate(inull, M_rot, STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error deallocation Matrizen inull... !"
  endif 
  deallocate(  MeshT1%CoeffsH1, MeshT1%CoeffsH2, &
               MeshT1%CoeffsH3, MeshT1%CoeffsH4, STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error deallocation Element MeshT1 !"
  endif 
  deallocate(  MeshT2%CoeffsH1, MeshT2%CoeffsH2, &
               MeshT2%CoeffsH3,  MeshT2%CoeffsH4, STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error deallocation Element MeshT2 !"
  endif
  do k=1,2
    deallocate( MeshGen(k)%CoeffsH1, MeshGen(k)%CoeffsH2, &
                MeshGen(k)%CoeffsH3, MeshGen(k)%CoeffsH4, STAT = allocStat)
    if (allocStat.ne.0) then
       write(*,702) " Error deallocation polynoms elmt k= !", k
    endif
  enddo
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

  print *,"========== PROGRAM TERMINATED CORRECTLY =========="


701 format(a)
702 format(a,i10)
704 format(i10,f10.5)
705 format (i10,i10)

end program mainTruss
