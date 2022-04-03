MODULE KIRCHHOFFBELROUND2D_mod
 ! this module is associated with the other module LagrangeRound2d_mod
 ! to implement the method of shell triangle elements by Belyuschko or
 ! denoted C0 because it has C0 continuity for the angle of rotation
 ! it has 15 dofs: 
 ! [ ux1, uy1, wz1, thx1, thy1, ..., thy3 ] (it has three nodes)
 ! it uses linear interpolation between the points, hence the derivatives
 ! are constant and there are no integration point, however the vertical
 ! displacement wz due to flexure is partitionized into wb and wr, which are
 ! the disp. due to pure bending and the disp. due to pure rigid body
 ! In this view two other variables are introduced:
 !  w.h = wb + r v omega, 
 ! where wb is due to pure bending, omega is rotation vector (Ox, Oy) of the rigid body of the whole element
 ! The second pair of variable is theta.hat or { thx.h ; thy.h }
 ! th.h = th.b + R*omega
 ! where th.b is derived from a curvature interpolation of {kxx,kyy,kxy}
 ! and R is the 2-identity matrix grouped 3 times vertically
 ! wh and th.h are used to compute the shear displacement vector {gax, gay}
 ! {gax, gay}={thy + wz,x ; -thx + wz,y} (here the colon indic the derivative)
 ! {gax, gay}=Brd*{th} + Bsd*{wz}
 ! where Brd and Bsd are simple position and derivation interpolation matrices
 ! There is a very important assumption about the last expression
 ! that is (related to Kirchhoff theory): 
 !   {0, 0} = Brd*{th.h} + Bsd*{w.h}
 ! and (another approximation):
 !   Brd*th ~= Brd*th.h   leads finally to 
 ! {gax, gay}=Brd*{wz-w.h}
 ! Then two previous relations are used to express {th.h} and {w.h} wrt unknown ! quantities: 
 ! {th.h}=A*{th}+R*omega 
 ! and the minimalization of the norm of ||th-th.h|| leads to
 ! omega= -1/3*R*(id-A)*{th}
 ! {w.h}= -1/3*X*(id-A)*{th} where X=PCrossVect*R bec. w=0 + r v omega (nodes)
 ! The worflow for determining a finite element full stress with 15 displacem,
 ! based on Galerkin method is given below:
 ! ---------------------------|
 ! |                          |
 ! |> omega -(2)->----------- | --(omega,4)--> w.h --> w.s  <--- {wz}
 !               th.h -(3)-->th-th.h                    |
 !  {th}  -(1)-->  -----------^                         |
 !   |                                                  |
 !   --(Brd 5)------------------> {gax,gay} <---(Bsd 6)--
 !
 ! Finally, take
 ! {sigs.x, sigs.y}=G*diag[5/6, 5/6]*{gax, gay}  (7)
 ! and factorize by the dof displ and rot to obtain by integration 4 rigidity
 ! matrices:
 ! (8) Kb,ww = shear, wz to wz matrix
 ! (9) Kb,wO+Ks,wO = sum of coupling contributions of bending mom. and shear
 ! (10) Kb,Oz+Ks,Oz = should be symetric to the last
 ! (11) Kb,OO+Ks,OO = sum of contribution of normal bending mon. and shear

 ! (1) subroutine lgc0fillmatrixangang
 ! (2) subroutine lgc0fillmatrixangrot
 ! (3) subroutine lgc0fillmatrixangrbody
 ! (4) subroutine lgc0fillmatrixrotWS
 ! (5) subroutine lgc0fillmatrixintangsh
 ! (6) subroutine lgc0fillmatrixderWSzsh
 ! (8) subroutine lgc0integmatrixKww
 ! (9) subroutine lgc0integmatrixKwth
 ! (11) subroutine lgc0integmatrixKthth

 ! Finally a sparse matrix assemblier based on the two subroutines
 ! in LagrangeRound2D_mod module
 ! linlgassemblyPqrtAMem  and PartB but adapted to the bending and shear here
 !shall be used
 
 use TypesRound
 implicit none

 ! global variables
 real,pointer               :: matrixA(:,:) ! dimension(6,6)
 real,pointer               :: matrixB(:,:) ! dimension(3,6)
 real,pointer               :: matrixC(:,:) ! dimension(6,3)
 real,pointer               :: matrixR(:,:) ! dimension(2,6)
 real,pointer               :: matrixX(:,:) ! dimension(3,6)
 real,pointer               :: matrixBb(:,:) ! dimension(3,6)
 real,pointer               :: matrixBr(:,:) ! dimension(2,6)
 real,pointer               :: matrixBs(:,:) ! dimension(2,3)
 real,pointer               :: matrixX1A(:,:) ! dimension(3,6)

 real,pointer               :: matrixDb(:,:) ! dimension(3,3)
 real,pointer               :: matrixGs(:,:) ! dimension(2,2)

 real,pointer               :: matrixKwwindm1(:,:) ! dimension(2,3) 
 real,pointer               :: matrixKwwinc(:,:) ! dimension(3,3)

 real,pointer               :: matrixKwthindm1(:,:) ! dimension(3,6)
 real,pointer               :: matrixKwthindm2(:,:) ! dimension(2,2)
 real,pointer               :: matrixKwthinc(:,:) ! dimension(3,6)
 real,pointer               :: matrixKwthsum(:,:) ! dimension(3,6)
 
 real,pointer               :: matrixKththindm1(:,:) ! dimension(3,6)
 real,pointer               :: matrixKththindm2(:,:) ! dimension(2,6)
 real,pointer               :: matrixKththinc(:,:) ! dimension(6,6)
 real,pointer               :: matrixKththsum(:,:) ! dimension(6,6)

 real,pointer               :: matrixKall(:,:) ! dimension(9,9)
 
 public :: lgc0allocatematrices, lgc0fillmatrixangang, lgc0fillmatrixangrot, & 
         lgc0fillmatrixangrbody, lgc0fillmatrixrotWS, lgc0deallocatematrices
 
 contains

 !************************************ 
 subroutine lgc0allocatematrices
 use typesround

 implicit none
 integer                  :: allocStat

 allocate(matrixA(6,6), matrixB(3,6), matrixC(6,3), STAT=allocStat )

 ! dimension(6,3)
 real,pointer               :: matrixR(:,:) ! dimension(2,6)
 real,pointer               :: matrixX(:,:) ! dimension(3,6)
 real,pointer               :: matrixBb(:,:) ! dimension(3,6)
 real,pointer               :: matrixBr(:,:) ! dimension(2,6)
 real,pointer               :: matrixBs(:,:) ! dimension(2,3)
 real,pointer               :: matrixX1A(:,:) ! dimension(3,6)

 real,pointer               :: matrixDb(:,:) ! dimension(3,3)
 real,pointer               :: matrixGs(:,:) ! dimension(2,2)

 real,pointer               :: matrixKwwindm1(:,:) ! dimension(2,3) 
 real,pointer               :: matrixKwwinc(:,:) ! dimension(3,3)

 real,pointer               :: matrixKwthindm1(:,:) ! dimension(3,6)
 real,pointer               :: matrixKwthindm2(:,:) ! dimension(2,2)
 real,pointer               :: matrixKwthinc(:,:) ! dimension(3,6)
 real,pointer               :: matrixKwthsum(:,:) ! dimension(3,6)
 
 real,pointer               :: matrixKththindm1(:,:) ! dimension(3,6)
 real,pointer               :: matrixKththindm2(:,:) ! dimension(2,6)
 real,pointer               :: matrixKththinc(:,:) ! dimension(6,6)
 real,pointer               :: matrixKththsum(:,:) ! dimension(6,6)

 real,pointer               :: matrixKall(:,:) ! dimension(9,9)
 
 if (allocStat.ne.0) then
    print *,"---lgc0allocatematrices: error could not allocate memory"
 endif
 end subroutine lgc0allocatematrices


END MODULE KIRCHHOFFBELROUND2D_mod

