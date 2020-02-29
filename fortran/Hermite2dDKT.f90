module Hermite2dDKT_mod
  ! Calculus expressions of first and second derivatives of Discrete Kirchhoff Theory
  use typesround

  implicit none

  private

  type tTrigElmtAngular
    real                :: cos1,cos2,cos3,sin1,sin2,sin3
    real                :: l1,l2,l3
  end type tTrigElmtAngular

  interface hermangular
     module procedure hermangular
  end interface hermangular

  interface hermderivativesX
     module procedure hermderivativesX
  end interface hermderivativesX

  interface hermderivativesY
     module procedure hermderivativesY
  end interface hermderivativesY

  interface hermsecondderivsX
     module procedure hermsecondderivsX
  end interface hermsecondderivsX

  interface hermsecondderivsY
     module procedure hermsecondderivsY
  end interface hermsecondderivsY
  
 public :: hermangular, hermderivativesX, hermderivativesY, hermsecondderivsX,hermsecondderivsY 

  contains
 !*********************************
 ! procedure hermangular, assumes uz is zero (plane coordinates transmitted)
  subroutine hermangular(ux,uy,uz,anglevalues)
  type(tNode)         :: ux,uy,uz
  type(tTrigElmtAngular) :: anglevalues
 ! use
  intent(in)          :: ux,uy,uz
  intent(inout)       :: anglevalues

  anglevalues%l1=0.
  anglevalues%l2=0.
  anglevalues%l3=0.

  anglevalues%l1=anglevalues%l1+(ux%vects(2)-ux%vects(1))**2
  anglevalues%l2=anglevalues%l2+(ux%vects(3)-ux%vects(1))**2
  anglevalues%l3=anglevalues%l3+(ux%vects(3)-ux%vects(2))**2

  anglevalues%l1=anglevalues%l1+(uy%vects(2)-uy%vects(1))**2
  anglevalues%l2=anglevalues%l2+(uy%vects(3)-uy%vects(1))**2
  anglevalues%l3=anglevalues%l3+(uy%vects(3)-uy%vects(2))**2

  anglevalues%l1=anglevalues%l1+(uz%vects(2)-uz%vects(1))**2
  anglevalues%l2=anglevalues%l2+(uz%vects(3)-uz%vects(1))**2
  anglevalues%l3=anglevalues%l3+(uz%vects(3)-uz%vects(2))**2

  anglevalues%l1=sqrt(anglevalues%l1)
  anglevalues%l2=sqrt(anglevalues%l2)
  anglevalues%l3=sqrt(anglevalues%l3)

  anglevalues%cos1=(ux%vects(3)-ux%vects(2))/anglevalues%l1
  anglevalues%cos2=(ux%vects(3)-ux%vects(1))/anglevalues%l2
  anglevalues%cos3=(ux%vects(2)-ux%vects(1))/anglevalues%l3

  anglevalues%sin1=(uy%vects(3)-uy%vects(2))/anglevalues%l1
  anglevalues%sin2=(uy%vects(3)-uy%vects(1))/anglevalues%l2
  anglevalues%sin3=(uy%vects(2)-uy%vects(1))/anglevalues%l3

  end subroutine hermangular


 !*********************************
 ! procedure hermderivativesX, assumes uz is zero (plane coordinates transmitted)
  subroutine hermderivativesX(anglevalues,ksi,eta,derivksi,deriveta)
  use typesround
  implicit none
  type(tTrigElmtAngular) :: anglevalues
  real                :: ksi, eta
  real,pointer        :: derivksi(:)
  real,pointer        :: deriveta(:)
  ! use
  intent(in)          :: ksi,eta,anglevalues
  intent(inout)       :: derivksi,deriveta
  ! lokale
  integer             :: i
  real                :: cos1,cos2,cos3,sin1,sin2,sin3
  real                :: l1,l2,l3
  
  l1=anglevalues%l1
  l2=anglevalues%l2
  l3=anglevalues%l3

  cos1=anglevalues%cos1
  cos2=anglevalues%cos2
  cos3=anglevalues%cos3
  sin1=anglevalues%sin1
  sin2=anglevalues%sin2
  sin3=anglevalues%sin3

  derivksi(1)=6.*(1-2*ksi-eta)*sin3/l3+6.*eta*sin2/l2
  deriveta(1)=-6.*ksi*sin3/l3+6.*(1-ksi-2*eta)*sin2/l2
  
  deriveta(2)=-2-3*(1-2*ksi-eta)*sin3**2+3*eta*sin2**2
  derivksi(2)=-2.+3*ksi*sin3**2-3*(1-ksi-2*eta)*sin2**2 

  derivksi(3)=3*(1-2*ksi-eta)*sin3*cos3-3*eta*sin2*cos2
  deriveta(3)=-3*ksi*sin3*cos3+3*(1-ksi-2*eta)*sin2*cos2

  derivksi(4)=6*(1-2*ksi-eta)*sin3/l3-6*eta*sin1/l1
  deriveta(4)=-6*ksi*(sin3/l3+sin1/l1)
  
  derivksi(5)=1.-3*(1-2*ksi-eta)*sin3**2-3*eta*sin1**2
  deriveta(5)=3*ksi*(sin3**2-sin1**2)

  derivksi(6)=3*(1-2*ksi-eta)*sin3*cos3+3*eta*sin1*cos1
  deriveta(6)=3*ksi*(sin1*cos1-sin3*cos3)

  derivksi(7)=6*eta*(sin2/l2-sin1/l1)
  deriveta(7)=-6*(1-2*eta*ksi)*sin2/l2-6*ksi*sin1/l1

  derivksi(8)=3*eta*(sin2**2-sin1**2)
  deriveta(8)=1.-3*(1-ksi-2*eta)*sin2**2-3*ksi*sin1**2

  derivksi(9)=3*eta*(sin1*cos1-sin2*cos2)
  deriveta(9)=3*(1-ksi-2*eta)*sin2*cos2-3*ksi*sin1*cos1

  end subroutine hermderivativesX


 !*********************************
 ! procedure hermseconderivs
  subroutine hermsecondderivsX(anglevalues,ksi,eta,deriv2ksi,deriv2ksieta,deriv2eta)
    use typesround
  implicit none
  type(tTrigElmtAngular) :: anglevalues
  real                :: ksi, eta
  real,pointer        :: deriv2ksi(:)
  real,pointer        :: deriv2ksieta(:)
  real,pointer        :: deriv2eta(:)
  ! use
  intent(in)          :: ksi,eta,anglevalues
  intent(inout)       :: deriv2ksi,deriv2ksieta,deriv2eta
  ! lokale
  integer             :: i
  real                :: cos1,cos2,cos3,sin1,sin2,sin3
  real                :: l1,l2,l3
  
  l1=anglevalues%l1
  l2=anglevalues%l2
  l3=anglevalues%l3

  cos1=anglevalues%cos1
  cos2=anglevalues%cos2
  cos3=anglevalues%cos3
  sin1=anglevalues%sin1
  sin2=anglevalues%sin2
  sin3=anglevalues%sin3


  deriv2ksi(1)=-12*sin3/l3
  deriv2ksieta(1)=6*(sin2/l2-sin3/l3)
  deriv2eta(1)=-12*sin2/l2

  deriv2ksi(2)=6*sin3**2
  deriv2ksieta(2)=3*(sin2**2+sin3**2)
  deriv2eta(2)=6*sin2**2

  deriv2ksi(3)=3*sin3**2
  deriv2ksieta(3)=3*sin2**2
  deriv2eta(3)=3*sin2**2

  deriv2ksi(4)=-12*sin3/l3
  deriv2ksieta(4)=-6*sin1/l1
  deriv2eta(4)=0.

  deriv2ksi(5)=6*sin3**2
  deriv2ksieta(5)=3*(sin3**2-sin1**2)
  deriv2eta(5)=0.

  deriv2ksi(6)=-6*sin3*cos3
  deriv2ksieta(6)=3*(sin1*cos1-sin3*cos3)
  deriv2eta(6)=0.

  deriv2ksi(7)=0.
  deriv2ksieta(7)=6*(sin2/l2-sin1/l1)
  deriv2eta(7)=12*sin2/l2

  deriv2ksi(8)=0.
  deriv2ksieta(8)=3*(sin2**2-sin1**2)
  deriv2eta(8)=6*sin2**2

  deriv2ksi(9)=0.
  deriv2ksieta(9)=3*(sin1*cos1-sin2*cos2)
  deriv2eta(8)=-6*sin2*cos2

  end subroutine hermsecondderivsX


 !*********************************
 ! procedure hermderivativesY, assumes uz is zero (plane coordinates transmitted)
  subroutine hermderivativesY(anglevalues,ksi,eta,derivksi,deriveta)
  use typesround
  implicit none
  type(tTrigElmtAngular) :: anglevalues
  real                :: ksi, eta
  real,pointer        :: derivksi(:)
  real,pointer        :: deriveta(:)
  ! use
  intent(in)          :: ksi,eta,anglevalues
  intent(inout)       :: derivksi,deriveta
  ! lokale
  integer             :: i
  real                :: cos1,cos2,cos3,sin1,sin2,sin3
  real                :: l1,l2,l3
  
  l1=anglevalues%l1
  l2=anglevalues%l2
  l3=anglevalues%l3

  cos1=anglevalues%cos1
  cos2=anglevalues%cos2
  cos3=anglevalues%cos3
  sin1=anglevalues%sin1
  sin2=anglevalues%sin2
  sin3=anglevalues%sin3

  derivksi(1)=6.*(1-2*ksi-eta)*cos3/l3-6.*eta*cos2/l2
  deriveta(1)=-6.*ksi*cos3/l3+6.*(1-ksi-2*eta)*cos2/l2
  
  deriveta(2)=3*(1-2*ksi-eta)*cos3*sin3-3*eta*sin2*cos2
  derivksi(2)=-3*ksi*sin3*cos3-3*(1-ksi-2*eta)*sin2*cos2

  derivksi(3)=-3*(1-2*ksi-eta)*cos3**2+3*eta*cos2**2
  deriveta(3)=-3*ksi*cos3**2-3*(1-ksi-2*eta)*cos2**2


  derivksi(4)=6*(1-2*ksi-eta)*cos3/l3-6*eta*cos1/l1
  deriveta(4)=-6*ksi*(cos3/l3+cos1/l1)
  
  derivksi(5)=3*(1-2*ksi-eta)*cos3*sin3+3*eta*cos1*sin1
  deriveta(5)=-3*ksi*cos3*sin3+3*ksi*cos1*sin1

  derivksi(6)=1.-3*(1-2*ksi-eta)*cos3**2-3*eta*cos1**2
  deriveta(6)=3*ksi*(cos3**2-cos1**2)

  derivksi(7)=-6*eta*(cos2/l2+cos1/l1)
  deriveta(7)=6*(1-2*eta*ksi)*cos2/l2-6*ksi*cos1/l1

  derivksi(8)=3*eta*(sin2*cos2-sin1*cos1)
  deriveta(8)=-3*(1-ksi-2*eta)*sin2*cos2-3*ksi*sin1*cos1

  derivksi(9)=3*eta*(cos2*2-cos1**2)
  deriveta(9)=1.-3*(1-ksi-2*eta)*cos2**2-3*ksi*cos1**2

  end subroutine hermderivativesY


 !*********************************
 ! procedure hermseconderivs
  subroutine hermsecondderivsY(anglevalues,ksi,eta,deriv2ksi,deriv2ksieta,deriv2eta)
    use typesround
  implicit none
  type(tTrigElmtAngular) :: anglevalues
  real                :: ksi, eta
  real,pointer        :: deriv2ksi(:)
  real,pointer        :: deriv2ksieta(:)
  real,pointer        :: deriv2eta(:)
  ! use
  intent(in)          :: ksi,eta,anglevalues
  intent(inout)       :: deriv2ksi,deriv2ksieta,deriv2eta
  ! lokale
  integer             :: i
  real                :: cos1,cos2,cos3,sin1,sin2,sin3
  real                :: l1,l2,l3
  
  l1=anglevalues%l1
  l2=anglevalues%l2
  l3=anglevalues%l3

  cos1=anglevalues%cos1
  cos2=anglevalues%cos2
  cos3=anglevalues%cos3
  sin1=anglevalues%sin1
  sin2=anglevalues%sin2
  sin3=anglevalues%sin3


  deriv2ksi(1)=-12*cos3/l3
  deriv2ksieta(1)=6*(cos2/l2-cos3/l3)
  deriv2eta(1)=-12*cos2/l2

  deriv2ksi(2)=-6*sin3*cos3-6*cos2*sin2  ! TO VERIFY !
  deriv2ksieta(2)=-3*sin3*cos3
  deriv2eta(2)=-6*sin2*cos2

  deriv2ksi(3)=6*cos3**2
  deriv2ksieta(3)=6*cos2**2
  deriv2eta(3)=3*cos2**2

  deriv2ksi(4)=12*cos3/l3
  deriv2ksieta(4)=0.
  deriv2eta(4)=6*(cos1/l1+cos3/l3)

  deriv2ksi(5)=-6*sin3*cos3
  deriv2ksieta(5)=0.
  deriv2eta(5)=3*(cos1*sin1-cos3*sin3)

  deriv2ksi(6)=6*cos3**2
  deriv2ksieta(6)=0.
  deriv2eta(6)=3*(cos3**2-cos1**2)

  deriv2ksi(7)=0.
  deriv2ksieta(7)=-6*(cos2/l2+cos1/l1)
  deriv2eta(7)=-12*cos2/l2

  deriv2ksi(8)=0.
  deriv2ksieta(8)=3*(cos2**2-cos1**2)
  deriv2eta(8)=-6*cos2*sin2

  deriv2ksi(9)=0.
  deriv2ksieta(9)=3*(cos2**2-cos1**2)
  deriv2eta(8)=-6*cos2**2

  end subroutine hermsecondderivsY


 end module Hermite2dDKT_mod
