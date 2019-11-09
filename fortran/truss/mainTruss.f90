program mainTruss
  use typesbalken
  use RigidTrussContacts_mod
  use printabstractstruct_mod

  implicit none

  integer               :: ne=3 ! Freiheitsgraden Anzahl
  integer               :: n=3  ! Nodes Anzahl
  real                  :: anglescalar
  type(tMeshCoord)      :: MeshPunkte
  type(tMeshElmt)       :: MeshT1, MeshT2
  type(tMeshElmt),dimension(:),allocatable  :: MeshGen(:)   ! Elementen 1-2 und 1-3
  type(tVarElmt),dimension(:),allocatable   :: VarGen(:)    ! Werten (U,F) in Elementen
  integer,pointer       :: elmtabbdung(:,:) 
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

  allocate(inull(1), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate null"
  endif
  allocate(M_rot(9), STAT=allocstat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate M_rot"
  endif
  allocate(MeshPunkte%x(1:n), MeshPunkte%y(1:n), MeshPunkte%z(1:n), STAT=allocStat)
  if (allocstat.ne.0) then
    print *,"mainTruss: error allocate M_rot"
  endif

! ***************************** allocate Elementen Strukturen (allocatable) *****************************************

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
  MeshPunkte%x(1)=1.0
  MeshPunkte%y(1)=0.0
  MeshPunkte%z(1)=0.0
  MeshPunkte%x(2)=0.0
  MeshPunkte%y(2)=0.0
  MeshPunkte%z(2)=0.0
  MeshPunkte%x(3)=0.0
  MeshPunkte%y(3)=0.5
  MeshPunkte%z(3)=0.0

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
  allocate(Ke_H%Ke(1:ne*n,1:ne*n), Var_H%Ue(1:ne*n), Var_H%Fe(1:ne*n), &
           Exakt%Ux(1:ne*n), Exakt%Fx(1:ne*n), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation global matrix !"
  endif  
  do i=1,ne*n       ! initialisierung
    do j=1,ne*n
       Ke_H%Ke(i,j)=0.0
    enddo
    Var_H%Ue(i)=0.0
    Var_H%Fe(i)=0.0   
  enddo

  allocate(elmtabbdung(1:n,1:n), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation element definition matrix !"
  endif 
  do i=1,n       ! initialisierung
    do j=1,n
       elmtabbdung(i,j)=0.0
    enddo
  enddo
  elmtabbdung(1,2)=2 ! 2=Balken
  elmtabbdung(2,1)=2 ! ..Symmetric
  elmtabbdung(1,3)=1 ! 1=Rod
  elmtabbdung(3,1)=1 ! ..Symmetric

  allocate(Dunbekannt(1:n,1:ne), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation unbekannt (frei rod or eingespannt) definition matrix !"
  endif 
  do i=1,n       ! initialisierung
    do j=1,ne
       Dunbekannt(i,j)=1 !  0=nicht unbekannt, 1=unbekannt fest eigespannnt,
                         ! 2=unbekannt translation aber freie rot
    enddo
  enddo
  Dunbekannt(1,1)=1
  Dunbekannt(1,2)=1 
  Dunbekannt(1,3)=1

  allocate(Dimponiert(1:n,1:ne), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation imposed freiheitsgrad matrix !"
  endif
  do i=1,n       ! initialisierung
    do j=1,ne
       Dimponiert(i,j)=1 !default
    enddo
  enddo
  Dimponiert(3,1)=0.0
  Dimponiert(3,2)=0.0
!  Dimponiert(3,3)=0.0
  Dimponiert(2,1)=0.0
  Dimponiert(2,2)=0.0
  Dimponiert(2,3)=0.0
! and 3 Unbekannte fuer Node 1 muessen vorher im Dunbekannt matrix definiert werden

  allocate(Kraftbewgg(1:n,1:ne), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation index zur berechnenden Bewegungskraefte matrix !"
  endif 
  do i=1,n       ! initialisierung
    do j=1,ne
       Kraftbewgg(i,j)=0    ! 1=zu bestimmen, 0=imponiert (siehe Kraftimponiert)
    enddo
  enddo
  Kraftbewgg(3,1)=1
  Kraftbewgg(3,2)=1
  Kraftbewgg(2,1)=1 
  Kraftbewgg(2,2)=1 
  Kraftbewgg(2,3)=1 
  Kraftbewgg(1,1)=1 

  allocate(Kraftimponiert(1:n,1:ne), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation imposed forces/moments in Bewegungspunkten matrix !"
  endif 
  do i=1,n       ! initialisierung
    do j=1,ne
       Kraftimponiert(i,j)=0.0    ! default
    enddo
  enddo
  Kraftimponiert(3,3)=0.0
  Kraftimponiert(1,2)=-10.0
  Kraftimponiert(1,3)=0.0
! es fehlt lediglich das verteiltes Last  (negative Y-Richtg) am Element (1-2)
! es wird in dem Struktur MeshT1 gespeichert



! <-------------------- DO SOMETHING --------------------->
  write(*,701) "----- mainTruss"
  k=1
  km1=1
  do jj=1,n
    do kk=jj,n
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

end program mainTruss
