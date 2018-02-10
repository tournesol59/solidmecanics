module CalcGradientSurf_mod


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


! Alle diese Tablen sind 2-dimensionale, Skalarfelds fuer jeden Node
! except Mesh, der beinhaltet andere Werte (dy, dyq, dyy, dyyq)

    type(tMesh)   :: Gitter      ! u.a. hat x(:), y(:) Werte der Geometrie   ! 
    real, pointer :: VirtVeschiebg(:, :)  ! array pointer: Verschiebung Koeffizienten (nraumx,nraumy) !
			            !    die Koeffizienten mal die Basisfunktionen machen !
                                    !    die reale Verschiebung !

    real, pointer    :: VirtVerschiebMat(:, :)  ! array pointer: virtuelle Verschiebung-Funktionen !
                                                ! Scalar Produkten !
    !
  
  public            :: allocateGradientSurf, deallocateGradientSurf
  public            :: interpolateRealElement

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


  !**************************************************************************!
  !                                                                          !
  ! Subroutine lderivf                                                       !
  !                                                                          !
  !      procedure: berechnet die Integral des Verschiebungsgradients auf    !
  !      einer realen Flaeschenelement von Indiz k; Ausgang = Vektor(3)       ! 
  !                                                                          !
  !**************************************************************************!

  SUBROUTINE lgderivf(n,Gitter,ki,kj,dF,der_x,der_y,res)
     use types
     implicit none
     integer,intent(in)           :: n  ! Order der Interpolation2D 2 =default
     type(tMesh),intent(in)       :: Gitter ! Gitter definiert in "TypesDef"
     integer,intent(in)           :: ki
     integer,intent(in)           :: kj ! Indizen des Elements
     real,intent(in)              :: dF ! Flaeschen des Elements = (a1 x a2).n
     real,pointer,intent(in)      :: der_x(:) ! Eingang = Vektor(2) pointer von
                                     ! realen auf-x differenzierten Basisvektoren
     real,pointer,intent(in)      :: der_y(:) ! Eingang = Vektor(2) pointer von
                                     ! realen auf-y differenzierten Basisvektoren

     real,dimension(1:3),intent(out) :: res  ! Resultat der Berechnung
     ! lokal
     integer                      :: i,j,k
     integer                      :: nrows=5,ncols=5 ! Anzeige Parameters
     integer                      :: xs, ys    ! Anzeige step fuer den Matrix der_x
     CHARACTER(LEN=*), PARAMETER  :: FMT7I = "(4I5, 4I5, 4I5, 4I5, 4I5, 4I5, 4I5)" !, 4I5, 4I5, 4I5, 4I5)"
     CHARACTER(LEN=*), PARAMETER  :: FMT6E = "(4I5, 5E14.7, 5E14.7, 5E14.7, 5E14.7, 5E14.7, 5E14.7)" !, 5E14.7, 5E14.7)"
     xs = Gitter%nraumx / ncols
     ys = Gitter%nraumy / nrows
     WRITE(*,'(A12)')  "   lgderivf"    !  man zeigt die Vektor-Eingaenge der_x, 
     WRITE(*,'(A12)')  "vektor ax%dx" !  um die Werte mit den Ausgaenge zu pruefen
     WRITE(*,'(A12)')  " ==========="
     WRITE(*,FMT7I)  0, 1, 1+xs, 1+2*xs, 1+3*xs, 1+4*xs, 1+5*xs
     do i=1,nrows
       WRITE(*,FMT6E)  1+(i-1)*ys,der_x(1),der_x(1+xs),der_x(1+2*xs),der_x(1+3*xs),der_x(1+4*xs),der_x(1+5*xs)
     end do
     

     res(1)=0.0  ! stub procedure im Moment !
     res(2)=0.0
     res(3)=0.0
   END SUBROUTINE lgderivf

end module CalcGradientSurf_mod

