program demo_callC


use ISO_C_BINDING

implicit none

   type(c_ptr)        :: pc
  ! real, dimension(10), allocatable  :: tfd(:)
   real,pointer       :: tfd(:)
   real,pointer                        :: pf1(:)
   integer n,l,i,j,allocstat

   interface
      subroutine DEMO_CALLEE(pt3) bind(C, name="democallee")
         use ISO_C_BINDING, only : c_double
         implicit none
         real(c_double),intent(inout)  :: pt3
      end subroutine DEMO_CALLEE
   end interface

   n=1
   ! l=2
   allocate(tfd(1:n), stat=allocstat)
   if (allocstat.ne.0) then
     print *, 'Error not all fields allocated!'
   end if
   
   call DEMOCALLEE(tfd)
  
   write(*,103) "value transmitted from C: ", tfd(1)
103 format(a,f10.5)  

end program demo_callC
