! Use: in Gesetz(1) f_ = Sum,k-i mu,k-i*Phi,i(x,y)
! where:   f_ ist die virtuelle Verschiebung, das durch einen finiten Vektor !
!          von Basisfunktionen (hier Lagrange fuer kontinuous Approx) representiert ist !
!          Nachdem die Loesung des Problems gefunden ist,
!          wird die Evaluation der Basis Verschiebungen zum realen 
!          Verschiebungen mit diesem Modul gemacht
module Lagrange2D_mod
!---
implicit none
private
!---
 ! interface
 !   module procedure linlgevalvector
 ! end interface


!********************************************************************!
!  Variable Tabelle der innere produkt von Komposition-Gradienten    !
!        je nach Indiz des Netzelements und der Position der         !
!        subelement  in der Platte                                   !
!                                                                    !
!********************************************************************!
  real,dimension(12)   :: GradientSurf = (/ (/ 1.0, 0.0 /), &
                             (/ 0.0, 1.0 /), &
                             (/-1.0, 1.0 /), &
                             (/-1.0, 0.0 /), &
                             (/ 0.0,-1.0 /), &
                             (/ 1.0,-1.0 /) /)
! Dann Tabelle der alle moeglichen Skalar Produkten zw. den surf. Gradienten Vektoren !
! nur ein Teil davon wird benutzt : !
  real,dimension(21)    :: GradientScalarProd = &
                          (/ 1.0,  0.0, -1.0, -1.0,  0.0,  1.0, &
                                   1.0,  1.0,  0.0, -1.0, -1.0, &
                                         2.0,  1.0, -1.0, -2.0, &
                                               1.0,  0.0, -1.0, &
                                                     1.0,  1.0, &
                                                           2.0 /)

  ! type(tMesh)       :: Gitter            ! schon jetzt im solidmain !
  ! type(tNumeric)    :: VarNum            ! schon jetzt im solidmain allokation !
  ! real, pointer :: VirtVeschiebg(:, :)  !  schon jetzt im CalcGradientSurf Modul !

public :: linlgevalvector 

contains 
  
!********************************************************************!
!  procedure linlgevalvector: nimmt an, dass der finite element probl.!
!     geloest ist und berechnet die approximierten Verschiebungen    !
!      in der Platte                                                 !
!********************************************************************!
 SUBROUTINE linlgevalvector(Gitter, VirtVerschiebg, VarNum)
 use types
 use interpol2d_mod
 implicit none
  type(tMesh)           :: Gitter
  real,pointer          :: VirtVerschiebg(:,:)   ! nraumx*nraumy
  type(tNumeric)        :: VarNum
  intent(in)            :: Gitter
  intent(inout)         :: VirtVerschiebg, VarNum
  ! lokal Variables
  real                  :: valround(1:6), facets(1:6)
  integer               :: i,j,k,l

  do j=0,Gitter%nraumy-1
    do i=0,Gitter%nraumx-1
      k=j*(Gitter%nraumx)+i+1  ! Indiz des Nodes in der Loesung VarNum%Verschiebung
       
!      if (((i==0).and.(j==0)).or.((i==Gitter%nraumx-1).and.(j==0)) &
!        .or.((i==Gitter%nraumx-1).and.(j==Gitter%nraumy-1)) &
!        .or.(i==0).and.(j=Gitter%nraumy-1)))
!      VarNum%Verschiebung(1,k)=VirtVesrchiebg(i,j)   ! vielleicht muss ich VarNum als Redundant einfach loeschen spaeter !

    enddo
  enddo
  END SUBROUTINE linlgevalvector

end module Lagrange2D_mod

