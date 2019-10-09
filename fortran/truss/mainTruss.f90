program mainTruss
  use typesbalken
  use RigidTrussContacts
  implicit none

  integer               :: ne=3 ! Freiheitsgraden Anzahl
  integer               :: n=3  ! Nodes Anzahl
  type(tMeshElmt)       :: MeshT1, MeshT2 ! Elementen 1-2 und 1-3
  type(tVarElmt)        :: VarT1, VarT2   ! entsprechenden Werten (U,F) in Elementen 1-2, 1-3
  integer,pointer       :: elmtabbdung(:,:) 
  real,pointer          :: Ke_H(:,:)   ! Whole global rigidity matrix
  real,pointer          :: A_H(:,:), B_H(:,:)  ! Sub matrices from Ke_H (f=A.U,unknown + B.U,applied)
  integer,pointer       :: Dunbekannt(:,:)  ! Unbekannte Bewegungen zu berechnen
  real,pointer          :: Dimponiert(:,:)  ! imponierte Bewegungen
  integer,pointer       :: Kraftbewgg(:,:)  ! Result of Calc: unbekannte Kraefte in bewegungsfrei Lenkung
  real,pointer          :: Kraftimponiert(:,:)  ! Bekannte imponiert Kraefte

  integer               :: i,j,k,l,p,allocStat


! ****************************** Definitions der Elementen ***********************************
  MeshT1%typ = 2
  MeshT1%nraumx = 2        ! subdivsion, default=2
  MeshT1%startx = 1.0        ! node 1 x-koord  !
  MeshT1%endx = 0.0          ! node 2 x-koord  !
  MeshT1%starty = 0.0        ! node 1 y-koord  !
  MeshT1%endy = 0.5          ! node 2 y-koord  !
  MeshT1%startz = 0.0        ! node 2 z-koord  !
  MeshT1%endz = 0.0          ! node 2 z-koord  !

!  MeshT1%dlen = sqrt((MeshT1%endx - MeshT1%startx)**2+(MeshT1%endy - MeshT1%starty)**2+(MeshT1%endz - MeshT1%startz)**2)           ! Laengen des Rechengebiets (Balk) in x/y !
  MeshT1%SArea = 0.0004         ! Section area
  MeshT1%EY = 210000000            ! YOUNG Modulus
  MeshT1%vu = 0.3            ! Poisson Coeff
  MeshT1%CI= 0.005             ! kinetics inertia moment
  MeshT1%q1 = -0.25
  MeshT1%q2 = -0.25         ! verteilte  Lasten q(x)=x/xlen*q1+(1-x/xlen)*q2
!     real        :: q?
  MeshT1%node_1 = 1
  MeshT1%node_2 = 2
  allocate(  MeshT1%CoeffsH1(1:4),  MeshT1%CoeffsH2(1:4),  MeshT1%CoeffsH3(1:4),  MeshT1%CoeffsH4(1:4), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation polynoms elmt 1 !"
  endif  


  MeshT2%typ = 1
  MeshT2%nraumx = 2        ! subdivsion, default=2
  MeshT2%startx = 1.0        ! node 1 x-koord  !
  MeshT2%endx = 0.0          ! node 2 x-koord  !
  MeshT2%starty = 0.0        ! node 1 y-koord  !
  MeshT2%endy = 0.0          ! node 2 y-koord  !
  MeshT2%startz = 0.0        ! node 2 z-koord  !
  MeshT2%endz = 0.0          ! node 2 z-koord  !

!  MeshT2%dlen = sqrt((MeshT1%endx - MeshT1%startx)**2+(MeshT1%endy - MeshT1%starty)**2+(MeshT1%endz - MeshT1%startz)**2)           ! Laengen des Rechengebiets (Balk) in x/y !
  MeshT2%SArea = 0.0004         ! Section area
  MeshT2%EY = 210000000            ! YOUNG Modulus
  MeshT2%vu = 0.3            ! Poisson Coeff
  MeshT2%CI= 0.0             ! kinetics inertia moment
  MeshT2%q1 = 0.0
  MeshT2%q2 = 0.0         ! verteilte  Lasten q(x)=x/xlen*q1+(1-x/xlen)*q2
!     real        :: q?
  MeshT2%node_1 = 1
  MeshT2%node_2 = 2
  allocate(  MeshT2%CoeffsH1(1:4),  MeshT2%CoeffsH2(1:4),  MeshT2%CoeffsH3(1:4),  MeshT2%CoeffsH4(1:4), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation polynoms elmt 2 !"
  endif

! ***************************** allocate in the main program *********************************
  allocate(Ke_H(1:ne*n, 1:ne*n), STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation global matrix !"
  endif  
  do i=1,ne*n       ! initialisierung
    do j=1,ne*n
       Ke_H(i,j)=0.0
    enddo
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
       Dunbekannt(i,j)=1 ! 1=unbekannt, 0=nicht unbekannt
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
!  call Matrice_Ke_H_Truss(MeshT1, VarT1, Ke_H, 1,2, n, ne, elmtabbdung, Dunbekannt, Dimponiert, Kraftbewgg, Kraftimponiert)

!****************** DEALLOCATE ******************
  deallocate(Ke_H, elmtabbdung, STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error deallocation Matrizen !"
  endif  
  deallocate(Ke_H, STAT = allocStat)
  deallocate(Dunbekannt, Dimponiert, STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error deallocation Bewegungs definitionsmatrizen !"
  endif 
  deallocate(Kraftbewgg, Kraftimponiert, STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error deallocation Matrizen zur berechnenden Bewegungskraefte matrix !"
  endif 


  print *,"========== PROGRAM TERMINATED CORRECTLY =========="
end program mainTruss
