module accessories_mod

  use types

  !--------------------------------------------------------------------------!
  implicit none
  private
  !-------

  interface matmultiply
    module procedure matmultiply
  end interface

  interface sqmatinverse
    module procedure sqmatinverse
  end interface

  interface sqmattransform
    module procedure sqmattransform
  end interface

  !--------------------------------------------------------------------------!
  public :: matmultiply, sqmatinverse, sqmattransform
  !-------------------------------------------------------------------------!
 contains

  !--------------------------------------------------------------------------!
  ! matmultiply: multiply zwei Matrizen A(:,:) und B(:,:), angenommen dass dims
  !             (n,r) und (r,m) sind. Resultate ist in C(:,:) gespeichert
  !--------------------------------------------------------------------------!
  subroutine matmultiply(n, r, m, A, B, C)
  !--------------------------------------------------------------------------!
  ! In Variablen
  !   dims n, r, m
  !   2-arrays A,B und C
  integer, intent(in)         :: n,r,m
  real,pointer,intent(in)     :: A(:,:), B(:,:)
  real,pointer, intent(inout) :: C(:,:)
  !--------------------------------------------------------------------------!
  !  Lokal Variablen
  integer                     :: i,j,k

!  use accessories
   
   do i=1,n
     do j=1,m
       C(i,j)=0.0
       do k=1,r
         C(i,j)=C(i,j)+A(i,k)*B(k,j)
       enddo
     enddo
   enddo
   
  end subroutine matmultiply

  !--------------------------------------------------------------------------!
  !  sqmatinverse: invertiert eine square Matrix. Damit wird lapack sub
  !                 dgertf (LU dekomposition) und dgerti benutzt werden
  !--------------------------------------------------------------------------!
  subroutine sqmatinverse(n, A)
  !--------------------------------------------------------------------------!
  !  In Variablen
  integer,intent(in)          :: n
  real,pointer,intent(inout)  :: A
  !--------------------------------------------------------------------------!
  !  Lokal Variablen
  real,pointer                :: pivot(:), lapackworkspace(:)
  integer                     :: Info1,Info2,allocStat  

!  use accessories

  allocate(pivot(1:n), lapackworkspace(1:n*n), STAT=allocStat)
  call DGETRF(n, n, A, n, pivot, Info1)
  call DGETRI(n, A, n, pivot, n*n, lapackworkspace, Info2)

  end subroutine sqmatinverse

  !--------------------------------------------------------------------------!
  ! sqmattransform: eine square 3-Matrix A und eine 3*3 Basistransformation
  !                 T3 gegeben, diese berechnet die aequivalent Matrix 
  !                 inv(T3)*A*T=B
  !--------------------------------------------------------------------------!
  subroutine sqmattransform(A, T3, C)
  real,dimension(9),intent(in)     :: A, T3
  real,dimension(9),intent(inout)  :: C
  !--------------------------------------------------------------------------!
  !  Lokal Variablen
  integer                     :: i
  real,dimension(9)           :: B,T3inv

!  use accessories

  do i=1,9
    T3inv(i)=T3(i)
  enddo
  call sqmatinverse(3,T3inv)
  call matmultiply(3,3,3,A,T3,B)
  call matmultiply(3,3,3,T3inv,B,C)

  end subroutine sqmattransform


end module accessories_mod

