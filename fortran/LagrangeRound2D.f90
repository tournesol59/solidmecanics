! Interpolation means the calculation of displacements and rotations
! through the Element. It can be have a support of 6 nodes triangles
! (membrane effect with 3 integration points) or 3 nodes triangle 
! (D0 with Kirchhoff approx however if membrane is quadratic, 
! flexion is quadratic but limited to parabolic [ax^2+by, cy^2, dxy]).
!  
! The methods of this module shall be called in combination with a
! structure 'InterpolData' which contains two main things
!   a linked list of names of generic subroutines (differentiate
!     , sum, product of different matrix, evaluation results), that
!     have been already called or are to be called. Accessed via macros
!     this is not OOP but need to be organized.
!   a sparse matrix (see Module sparse CSR) that shall be filled 
!
! Use: interpolate with a0+a1*X+a2*Y+a3*X*Y
module LagrangeRound2D_mod
!---
implicit none
private
  interface interpolD0Curving
     module procedure interpolD0Curving
  end interface

public :: interpolD0Curving, interpolD0Membrane, interpolD0Shearing  
! , interpol3PMemEvalDerX, interpol3PMemEvalDerY, interpol3PProdMat,
! , interpol3PMemEvalPolXY, interpol3PCurvEvalPolXY, interpol3PShEvalKt
!, interpol3PShEvalKws

contains 
!********************************************************************!
!  add in functions /procedures used by interpolD0Curving !
!********************************************************************!
!********************************************************************!
!  procedure interpolD0Curving !
!********************************************************************!
  subroutine interpolD0Curving(llNodes, llFunctions, lokale)
    use TypesRound
    implicit none



  end subroutine interpolD0Curving
!********************************************************************!
!  procedure interpolD0Membrane !
!********************************************************************!
  
!********************************************************************!
!  procedure interpolD0Shearing !
!********************************************************************!

end module LagrangeRound2D_mod

