! my rapid test for C-F interoperability !
module demoCfortran
use ISO_C_BINDING
implicit none

  type,bind(C) :: tInfo
    real(c_double),dimension(10) :: x
  end type tInfo

  private
    interface SUBTEST
      module procedure SUBTEST
    end interface

    interface SUBCOPY
      module procedure SUBCOPY
   end interface

  public:: SUBTEST, SUBCOPY

  contains

! -------------------------------------------------------- !
!  SUBTEST called from C tries to allocate a contingent bloc and
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
   real*8, pointer  :: p1(:), pf1(:)  ! p1 would be allocatable but seems difficult
   real*8,target :: tf(10) = (/ 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. /)      ! let's try with this
   !------------------------------------!
   integer n, allocstat;

!   n=10
!   allocate(p1(1:n), stat=allocstat)
!   if (allocstat.ne.0) then
!     print *, 'Error not all fields allocated!'
!   end if
!   p1(1)=1.41
   ! pf => p1
   ! pd = c_loc(p1)  ! marche pas
   
!    write(*, 101) 'Field index 5 is: ', pd(5)
!    pf1 => p1               ! if pf not to array: Error: Different ranks in pointer assignment at (1)
   
    pc = c_loc(tf)
    call c_f_pointer(pc,pf1,[10])  ! pc must be intermediate
    pd = pc

    write(*, 101) 'print from Fortran, Field index 1 is: ', pf1(2)    


101 format(a,f10.5)  

   end subroutine SUBTEST
   
!-------------------------------------------------
!  SUBCOPY is the contrary: it calls C subroutine
!  thus, SUBCOPY is part of a PROGRAM (see main.f90)
!  which fills in a (**double) to begin, then tInfo
!  and SUBCOPY copy it to a local array,dimension(10,2) 
!  please take care of Fortran colums!


   subroutine SUBCOPY(pd) bind(C, name="subcopy")
   use ISO_C_BINDING
   implicit none
 !-----------------------------------------!
!   real     :: pd            ! error: arg must be all passed by ref, thus ptr or character,kind
   type(c_ptr)        :: pd 
   intent(inout)      :: pd
 !-----------------------------------------!
   type(c_ptr)        :: pc
   real, dimension(10), allocatable  :: tfd(:)
   real,pointer                        :: pf1(:)
   integer n,l,i,j,allocstat
   
   call   c_f_pointer(pc,pf1,[10])
   pc = pd
   write(*,102) 'print from Fortran first value: ', pf1(1) 

   n=10
   l=2
   allocate(tfd(1:n), stat=allocstat)
   if (allocstat.ne.0) then
     print *, 'Error not all fields allocated!'
   end if
   
102 format(a,f10.5)  

   end subroutine SUBCOPY 

end module demoCfortran
