! Use: in Gesetz(1) f_ = Sum,k-i mu,k-i*Phi,i(x,y)
! where:   f_ ist die virtuelle Verschiebung, das durch einen finiten Vektor !
!          von Basisfunktionen (hier Lagrange fuer kontinuous Approx) representiert ist !
!          Dieses Modul Sets the Problem based on 2D-Interpolation von Lagrange
module Lagrange2D_mod
!---
implicit none
private
!---
 ! interface
 !   module procedure linlg2dfillmatrixq
 ! end interface


!********************************************************************!
!  Variable Tabelle der innere produkt von Komposition-Gradienten    !
!        je nach Indiz des Netzelements und der Position der         !
!        subelement  in der Platte                                   !
!                                                                    !
!********************************************************************!

  real,dimension(8)   :: GaussPoints = (/ (/ 0.57735, 0.57735 /), &
                             (/ 0.57735, -0.57735 /), &
                             (/ -0.57735, 0.57735 /), &
                             (/ -0.57735, 0.57735 /) /)

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

public :: linlg2allocatesdns, computeGaussInt, linlg2dfillmatrixq

contains 

!********************************************************************!
!  procedure linlg2allocatesdns DENSE MATRICES: n=10--> 10000 memory !
!********************************************************************!
 SUBROUTINE linlg2allocatesdns(Gitter, SkalarProdMatrixQ, SkalarProdMatrixT)
 use types
 use interpol2D_mod
 use Lagrange1d_mod
 implicit none
  type(tMesh)         :: Gitter
  integer             :: allocStat  
  real,pointer        :: SkalarProdMatrixQ(:,:)
  real,pointer        :: SkalarProdMatrixT(:,:,:)
  intent(in)          :: Gitter
  intent(inout)       :: SkalarProdMatrixQ, SkalarProdMatrixT
  
 allocate(SkalarProdMatrixQ(1:Gitter%nraumx*Gitter%nraumy,1:Gitter%nraumx*Gitter%nraumy), & 
     SkalarProdMatrixT(1:Gitter%nraumx*Gitter%nraumy,1:Gitter%nraumx*Gitter%nraumy,1:4), & 
     STAT = allocStat)
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not allocate all variables!'
       STOP
    end if
 END SUBROUTINE linlg2allocatesdns

!********************************************************************!
!  procedure computeGaussInt   : Perform the filling of the SkalarProdMatrix   !
!              ! tt=0: Prop*Prop, 1:Prop*DerX, 2:Prop*DerY, 3: DerX*DerX, 4:DerX*DerY
!              ! ii: Art of the polynoms: position of one elemt to another
!********************************************************************!
 SUBROUTINE computeGaussInt(tt, ii, pol1, pol2, pol3, valend)
 use types
 use interpol2D_mod
 use Lagrange1d_mod
 implicit none
 integer,intent(in)	     :: tt, ii;

 real,intent(out)                    :: valend
 type(tPolynom),intent(in)           :: pol1, pol2, pol3
 real                                :: val1, val2, val3, val4
 real                                :: der1, der2, der3, der4
 integer                             :: j

 valend=0.0
 ! So polynoms sould be transmitted:
 ! pol1%d=2
 ! pol1%coeffs= (/ -1, 0, 1 /)
 ! pol2%d=2
 ! pol2%coeffs= (/ 1, 1, 0 /)
 ! pol3%d=2
 ! pol3%coeffs= (/ 1,-1, 0 /)
 select case(tt) ! Aim of the calculus: should we have Prop*Der, Prop*Prop, Der*Der..

 case(0)  !Prop*Prop
  do j=0,3 
   select case(ii)  ! Art of polynoms in the element
   case(0)   
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     !call evaluation(pol1, GaussPoints(2*j+1), val3)
     !call evaluation(pol1, GaussPoints(2*j+2), val4)
     valend=valend+val1*val2*val1*val2*1.0;
   case(-1)
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     call evaluation(pol2, GaussPoints(2*j+1), val3)
     call evaluation(pol2, GaussPoints(2*j+2), val4)
     valend=valend+val1*val2*val3*val4*1.0;
  case(+1)
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     call evaluation(pol3, GaussPoints(2*j+1), val3)
     call evaluation(pol3, GaussPoints(2*j+2), val4)
     valend=valend+val1*val2*val3*val4*1.0;
  end select
 enddo

 case(1)  ! Prop*DerX
  do j=0,3
   select case(ii)  ! Art of polynoms in the element
   case(0)      
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     call calcderivative(pol1, GaussPoints(2*j+1), der1)
     !call evaluation(pol1, GaussPoints(2*j+2), val3)
     valend=valend+val1*val2*der1*val2*1.0;
   case(-1)
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     call calcderivative(pol2, GaussPoints(2*j+1), der1)
     call evaluation(pol2, GaussPoints(2*j+2), val3)
     valend=valend+val1*val2*der1*val3*1.0;
  case(+1)
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     call calcderivative(pol3, GaussPoints(2*j+1), der3)
     call evaluation(pol3, GaussPoints(2*j+2), val3)
     valend=valend+val1*val2*der1*val3*1.0;
  end select
 enddo

 case(2)  ! Prop*DerY
  do j=0,3
   select case(ii)  ! Art of polynoms in the element
   case(0)      
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     !call evaluation(pol1, GaussPoints(2*j+1), val1)
     call calcderivative(pol1, GaussPoints(2*j+2), der1)
     valend=valend+val1*val2*val1*der1*1.0;
   case(-1)
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     call evaluation(pol2, GaussPoints(2*j+1), val3)
     call calcderivative(pol2, GaussPoints(2*j+2), der3)
     valend=valend+val1*val2*val3*der3*1.0;
  case(+1)
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     call evaluation(pol3, GaussPoints(2*j+1), val3)
     call calcderivative(pol3, GaussPoints(2*j+2), der3)
     valend=valend+val1*val2*val3*der3*1.0;
  end select
 enddo

 case(3)  ! DerX*DerX
  do j=0,3
   select case(ii)  ! Art of polynoms in the element
   case(0)      
     call calcderivative(pol1, GaussPoints(2*j+1), der1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     call calcderivative(pol1, GaussPoints(2*j+1), der2)
     !call evaluation(pol1, GaussPoints(2*j+2), val3)
     valend=valend+der1*val2*der2*val2*1.0;
   case(-1)
     call calcderivative(pol1, GaussPoints(2*j+1), der1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     call calcderivative(pol2, GaussPoints(2*j+1), der2)
     call evaluation(pol2, GaussPoints(2*j+2), val3)
     valend=valend+der1*val2*der2*val3*1.0;
  case(+1)
     call calcderivative(pol1, GaussPoints(2*j+1), der1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     call calcderivative(pol3, GaussPoints(2*j+1), der3)
     call evaluation(pol3, GaussPoints(2*j+2), val3)
     valend=valend+der1*val2*der3*val3*1.0;
  end select
 enddo

 case(4)  ! DerX*DerY
  do j=0,3
   select case(ii)  ! Art of polynoms in the element
   case(0)      
     call calcderivative(pol1, GaussPoints(2*j+1), der1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call calcderivative(pol1, GaussPoints(2*j+2), der2)
     valend=valend+der1*val2*val1*val2*1.0;
   case(-1)
     call calcderivative(pol1, GaussPoints(2*j+1), der1)
     call evaluation(pol1, GaussPoints(2*j+2), val1)
     call evaluation(pol2, GaussPoints(2*j+1), val2)
     call calcderivative(pol2, GaussPoints(2*j+2), der2)
     valend=valend+der1*val1*val2*der2*1.0;
  case(+1)
     call calcderivative(pol1, GaussPoints(2*j+1), der1)
     call evaluation(pol1, GaussPoints(2*j+2), val2)
     call evaluation(pol3, GaussPoints(2*j+1), val3)
     call calcderivative(pol3, GaussPoints(2*j+2), der3)
     valend=valend+der1*val2*val3*der3*1.0;
  end select
 enddo

 case(5)  ! DerY*DerY
  do j=0,3
   select case(ii)  ! Art of polynoms in the element
   case(0)      
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call calcderivative(pol1, GaussPoints(2*j+2), der1)
     call evaluation(pol1, GaussPoints(2*j+1), val2)
     call calcderivative(pol1, GaussPoints(2*j+2), der2)
     valend=valend+val1*der1*val2*der2*1.0;
   case(-1)
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call calcderivative(pol1, GaussPoints(2*j+2), der1)
     call evaluation(pol2, GaussPoints(2*j+1), val2)
     call calcderivative(pol2, GaussPoints(2*j+2), der2)
     valend=valend+val1*der1*val2*der2*1.0;
  case(+1)
     call evaluation(pol1, GaussPoints(2*j+1), val1)
     call calcderivative(pol1, GaussPoints(2*j+2), der1)
     call evaluation(pol3, GaussPoints(2*j+1), val3)
     call calcderivative(pol3, GaussPoints(2*j+2), der3)
     valend=valend+val1*der1*val3*der3*1.0;
  end select
 enddo

 end select
! return valend
 END SUBROUTINE computeGaussInt
 
!********************************************************************!
!  procedure linlg2dfillmatrixq:                                     !
!     
!********************************************************************!
 SUBROUTINE linlg2dfillmatrixq(Gitter, chicoeff, VarNum, SkalarProdMatrixQ)
 use types
 use interpol2D_mod
 use Lagrange1d_mod
 implicit none
 type(tMesh),intent(in)          :: Gitter
  real,intent(inout)             :: chicoeff(1:8,1:100)  ! 2nd Dimension ist
  type(tNumeric),intent(inout)   :: VarNum
  real,pointer,intent(inout)     :: SkalarProdMatrixQ(:,:)
  ! lokale Variablen
  integer                        :: i,j,k,r, l,m,n,p
  real                           :: K_x, K_y, costh, sinth, cosps, sinps, valuex, valuey
  type(tPolynom)                 :: pol1, pol2, pol3
   pol1%n=2
   pol1%coeffs= (/ -1.0, 0.0, 1.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
   pol2%n=2
   pol2%coeffs= (/ 1.0, 1.0, 0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
   pol3%n=2
   pol3%coeffs= (/ 1.0,-1.0, 0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

 ! Initialisierung
  do i=1,Gitter%nraumx*Gitter%nraumy
    do j=1,3
       VarNum%TensionTg(i,1)=0.0    ! sigma,xx; sigma,xy; sigma,yy
    enddo
  enddo  
  do i=1,Gitter%nraumx*Gitter%nraumy
     VarNum%Tension3(1,i)=0.0       ! sigma,xz
     VarNum%Tension3(2,i)=0.0       ! sigma,yz
  enddo 
  do i=1,Gitter%nraumx*Gitter%nraumy
    do j=1,4
       VarNum%BiegMoment(j,i)=0.0   ! z*sigma,xx, z*sigma,yy, z*sigma,xy=yx and letzte Teil fuer Vn=dMns/ds+Qn nur fuer RandBedingung
    enddo 
  enddo
  do i=1,(Gitter%nraumx*Gitter%nraumy)
    do j=1,(Gitter%nraumx*Gitter%nraumy)  ! big one
      SkalarProdMatrixQ(i,j)=0.0
    enddo
  enddo

 !  Lagrange2D polynom used for this interpolation has following coeff:
 !  Phi,i(X,Y)=(1+X)*(1-X)*(1-Y)*(1+Y)/(1)
 ! expanded:
 !  Phi,i(X,Y)= (1- X^2 - Y^2 + X^2*Y^2)
 !  Der%X = -2*X + 2*X*Y^2
 !  Der%Y = -2*Y + 2*Y*X^2


 do i=0,Gitter%nraumx-1
    do j=0,Gitter%nraumy-1
      k=j*Gitter%nraumx+i+1

      do m=0,Gitter%nraumx-1
        do n=0,Gitter%nraumy-1
          l=n*Gitter%nraumx+m+1
          p= l-k
          ! search the Kurvatures 
          K_x= chicoeff(1,k)
          K_y=chicoeff(8,k) 
          costh=cos((Gitter%x(m+2)-Gitter%x(m+1))*K_x)
          sinth=sin((Gitter%x(m+2)-Gitter%x(m+1))*K_x)
          cosps=cos((Gitter%y(n+2)-Gitter%y(n+1))*K_y)
          sinth=sin((Gitter%y(n+2)-Gitter%y(n+1))*K_y)

      ! first loop of the matrix with contribution of B_ x Tg(T_)= skalar
          select case(p) 
          case (0)  ! Elemts superposition, calc of Four pts Gauss Integral with Prop*DerX,Y
             call computeGaussInt(1, p, pol1, pol2, pol3, valuex) 
             call computeGaussInt(1, p, pol1, pol2, pol3, valuey) 
             SkalarProdMatrixQ(k,l)= (valuex*K_x + valuey*K_y)
          case (1)
             call computeGaussInt(1, p, pol1, pol2, pol3, valuex) 
             call computeGaussInt(1, p, pol1, pol2, pol3, valuey) 
             SkalarProdMatrixQ(k,l)= valuex*K_x + valuey*K_y
          case (-1)
             call computeGaussInt(1, p, pol1, pol2, pol3, valuex) 
             call computeGaussInt(1, p, pol1, pol2, pol3, valuey) 
              SkalarProdMatrixQ(k,l)= valuex*K_x + valuey*K_y  
          case (10)   !Gitter%nraumx ist nicht konstant aber fuer Praesent
             call computeGaussInt(1, p, pol1, pol2, pol3, valuex) 
             call computeGaussInt(1, p, pol1, pol2, pol3, valuey) 
              SkalarProdMatrixQ(k,l)= valuex*K_x + valuey*K_y
          case (-10)   !-Gitter%nraumx ist nicht konstant aber fuer Praesent
             call computeGaussInt(1, p, pol1, pol2, pol3, valuex) 
             call computeGaussInt(1, p, pol1, pol2, pol3, valuey) 
              SkalarProdMatrixQ(k,l)= valuex*K_x + valuey*K_y
          end select
        enddo
      enddo
    enddo
  enddo
! Display
    do j=0,Gitter%nraumy-1
      do i=0,Gitter%nraumx-1
         k=j*Gitter%nraumx+i+1
         write(*,302) '    Superline k   =' ,k,'  {'
         write(*,*)
         do l=0,Gitter%nraumy-1
             write(*,302) '    SkalarProd matrix column    =',j
             write(*,102) (SkalarProdMatrixQ(k,l*Gitter%nraumx+m+1), m=0,Gitter%nraumx-1)
        enddo
     enddo
   enddo
   write(*,*)

 102  format (f10.5)
 202  format (a,f10.5)
 302  format (a,i10)

 END SUBROUTINE linlg2dfillmatrixq

end module Lagrange2D_mod

