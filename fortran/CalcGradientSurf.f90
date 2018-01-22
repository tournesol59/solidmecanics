module CalcGradientSurf_mod

! Diese Module zeigt auf die Berechnung des Gradients eines tengentialen Tensors
! T auf die Ebene (x,y) definiert
! Gleichug: I1= Sum( grad(f_) && T__ds) auf die Flaeche S
! 	    I2= Sum( f_ & t0_dl) auf dem Kontur C2
!           I3=Sum( f_ & T__ & n_dl) auf dem Kontur C3
!
! Use: in Gesetz(1) : Sum(f_ & (p_ + div(T__(ks)) & div(M__(ks)) )ds
! where:   f_ ist die virtuelle Verschiebung (Vektor)
!          p_ ist die Schwerkraft
!          T__ ist den Tensor Tengentskraft auf die Ebene PI
!          M__ ist den Tensor Winkelbiegung??
!
!
use types
!use interpol2D
implicit none
!
private
 ! interface
!  module procedure allocateGradientSurf  ! dont know how to use it !
!          subroutine allocateGradientSurf       
!          end subroutine allocateGradientSurf        
!  end interface
 

!  interface
          !  module procedure deallocateGradientSurf
!          subroutine deallocateGradientSurf
!          end subroutine deallocateGradientSurf
!  end interface

  interface
          subroutine initGradientSurfReference
             ! nimmt das ReferenzNetz und differenziert alle Felder Biegung, Tension.. darauf
          end subroutine initGradientSurfReference
  end interface

  interface
          subroutine interpolateRealElement

          end subroutine interpolateRealElement
  end interface

  interface
          subroutine calcGradientInterpolation

          end subroutine calcGradientInterpolation
  end interface

  interface
          subroutine calcOverallGradient

          end subroutine calcOverallGradient
  end interface

  interface
          subroutine calcTensorMultGradient

          end subroutine calcTensorMultGradient
  end interface

! Alle diese Tablen sind 2-dimensionale, Skalarfelds fuer jeden Node
! except Mesh, der beinhaltet andere Werte (dy, dyq, dyy, dyyq)

    type(tMesh)   :: Gitter      ! u.a hat x(:), y(:) Werte der Geometrie   ! 
    real, pointer :: Biegung(:, :)  ! array pointer: Verschiebung in Z Richtung !
    real, pointer :: Biegmoment(:, :)  ! array pointer: BiegMoment auf X//Y Richtung !
    real, pointer :: Biegwinkel(:, :)  ! array pointer: Winkelaenderung in X//Y Achse !
    real, pointer :: VirtLoesung(:, :, :)  ! array pointer: Virtuelle Loesung in polynomial basis=sum(b,k*Ref,k(x,y)) !
    real, pointer    :: VirtVerschieb(:, :)  ! array pointer: virtuelle Verschiebung in Z Richtung !
    real, pointer :: Tension(:, :, :)   ! array pointer: Spannung tengential  !
    !

  public            :: allocateGradientSurf, deallocateGradientSurf
  public            :: initGradientSurfReference, interpolateRealElement

 contains
  !**************************************************************************!
  !                                                                          !
  ! subroutine allocateGradientSurf                                          !
  !                                                                          !
  !      Funktion: Allokieren der benötigten Felder                          !
  !                wird von XXX und den YYY                                  !
  !                Verfahren aufgerufen                                      !
  !                                                                          !
  !**************************************************************************!


  SUBROUTINE allocateGradientSurf( nraumx, nraumy )
  use types
  implicit none
    ! lokal Variable
    integer, INTENT(IN)                :: nraumx
    integer, INTENT(IN)                :: nraumy

    allocate( Biegung(nraumx, nraumy) ) 
    allocate( Tension(nraumx, nraumy, 2) )
    allocate( Biegmoment(nraumx, nraumy) ) 
    allocate( Biegwinkel(nraumx, nraumy) )
    allocate( VirtLoesung(nraumx, nraumy, 2) ) 
    allocate( VirtVerschieb(nraumx, nraumy) )
     
  END SUBROUTINE allocateGradientSurf

  !**************************************************************************!
  !                                                                          !
  ! subroutine deallocateGradientSurf                                        !
  !                                                                          !
  !      Funktion: deallokieren der benötigten Felder                        !
  !                wird von XXX und den YYY                                  !
  !                Verfahren aufgerufen                                      !
  !                                                                          !
  !**************************************************************************! 
  SUBROUTINE deallocateGradientSurf
  use types
  implicit none
  integer               :: IERR7

    deallocate( Biegung, STAT = IERR7) 
    deallocate( Tension, STAT = IERR7 )
    deallocate( Biegmoment, STAT = IERR7 ) 
    deallocate( Biegwinkel, STAT = IERR7 )
    deallocate( VirtLoesung, STAT = IERR7 ) 
    deallocate( VirtVerschieb, STAT = IERR7 ) 
     
  END SUBROUTINE deallocateGradientSurf


end module CalcGradientSurf_mod

