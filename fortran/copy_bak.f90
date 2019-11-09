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
  allocate(  MeshT1%CoeffsH1(1:4),  MeshT1%CoeffsH2(1:4),  MeshT1%CoeffsH3(1:4),  MeshT1%CoeffsH4(1:4) STAT = allocStat)
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
  MeshT2q1 = 0.0
  MeshT2%q2 = 0.0         ! verteilte  Lasten q(x)=x/xlen*q1+(1-x/xlen)*q2
!     real        :: q?
  MeshT2%node_1 = 1
  MeshT2%node_2 = 2
  allocate(  MeshT2%CoeffsH1(1:4),  MeshT2%CoeffsH2(1:4),  MeshT2%CoeffsH3(1:4),  MeshT2%CoeffsH4(1:4) STAT = allocStat)
  if (allocStat.ne.0) then
     print *," Error allocation polynoms elmt 2 !"
  endif
