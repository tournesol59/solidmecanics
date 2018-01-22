program Testformat
  use types
  use Hermite2d_mod
  implicit none
  real          :: Coefficients(1:10, 50)  ! Der Input: 10 Koefficients x NNodes
  real          :: Values(81)  ! Der Resultat fuer ein Node

  call calchermitepol(2, Coefficients)

  call evalhermitepol(2, 1, Coefficients, Values)

end program Testformat
