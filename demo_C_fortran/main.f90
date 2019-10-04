program demo_callC


use ISO_C_BINDING

use demoCfortran
implicit none

   type(c_ptr)        :: pc
   real, dimension(10), allocatable  :: tfd(:)
   real,pointer                        :: pf1(:)
   integer n,l,i,j,allocstat

   interface DEMO_CALLEE
      subroutine DEMO_CALLEE(pt3)
      real,(c_double)  :: pt3
   end interface

   n=10
   ! l=2
   allocate(tfd(1:n), stat=allocstat)
   if (allocstat.ne.0) then
     print *, 'Error not all fields allocated!'
   end if
   
   call DEMO_CALLEE(  
 
103 format(a,f10.5)  

end program demo_callC
