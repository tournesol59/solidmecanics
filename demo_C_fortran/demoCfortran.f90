! my rapid test for C-F interoperability !
module demoCfortran
use ISO_C_BINDING
implicit none

  type,bind(C) :: tInfo
    real(c_double),dimension(3) :: x
  end type tInfo

  private
    interface SUBTEST
      module procedure SUBTEST
    end interface

  public:: SUBTEST

  contains

! -------------------------------------------------------- !
!  SUBTEST tries to allocate a contingent bloc and
!    move a pointer got from type(c_ptr) to this bloc
!    (assumed arrays don't work)
!     a. declare pd
!     b. declare fixed array(size) as target and points C ptr, pc, to it
!     c. converts pc to F PTR pf
!     d. points pd to same target (=)
!     e. print test
!----------------------------------------------------------!

   subroutine SUBTEST(pd) bind(C, name="subtest")
   use ISO_C_BINDING
   implicit none

   !------------------------------------!
!   type(c_ptr),value  :: pd 
   type(c_ptr)        :: pd 
   intent(inout) :: pd   
   !------------------------------------!
   type(c_ptr) :: pc
   real(c_double), pointer  :: p1(:), pf1(:)  ! p1 would be allocatable but seems difficult
   double precision,target :: tf(10) = (/ 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. /)      ! let's try with this
   !------------------------------------!
   integer n, allocstat;
   n=10
   allocate(p1(1:n), stat=allocstat)
   if (allocstat.ne.0) then
     print *, 'Error not all fields allocated!'
   end if

   ! pf => p1
   ! pd = c_loc(p1)  ! marche pas
    pc = c_loc(tf)
    call c_f_pointer(pc,pf1,[10])
    write(*, 101) 'Field index 0 is: ', pf1(1)    
    
    pd = pc
   
!    write(*, 101) 'Field index 5 is: ', pd(5)
!    pf1 => p1               ! if pf not to array: Error: Different ranks in pointer assignment at (1)

101 format(a,f10.5)  

   end subroutine SUBTEST
   
end module demoCfortran
