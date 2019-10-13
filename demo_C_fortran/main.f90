program demo_callC


use ISO_C_BINDING

implicit none

   type(c_ptr)        :: pc
!!   real*8, dimension(10), allocatable  :: tfd(:)
   real*8,pointer       :: tfd(:)
  ! real,pointer                        :: pf1(:)
   integer n,l,i,j,allocstat
  integer ii,jj,kk
   common/ijk/ ii, jj, kk
   real*8 ff
   character*32 cc

   interface
      subroutine DEMO_CALLEE(d3) bind(C, name="democallee")
         use ISO_C_BINDING, only : c_double
         implicit none
         real(c_double),intent(inout)          :: d3
!         real(c_double),intent(inout),pointer  :: pt3(:)
      end subroutine DEMO_CALLEE
   end interface

   n=1
   ! l=2
   allocate(tfd(1:n), stat=allocstat)
   if (allocstat.ne.0) then
     print *, 'Error not all fields allocated!'
   end if
   do i=1,n
     tfd(i)=1.5
   enddo
   write(*,103) "value copied to C: ", tfd(1)

   call DEMO_CALLEE(tfd(1))


   ii=2
   jj=3
   kk=4
   ff=9.0567
   cc='Example of a character string'
   
   write(6,100) ii,ff
   call abc(ii)
   write(6,200) ii
   write(6,300) ii,jj,kk

   call doubleIJK(cc)
   write(6,300) ii,jj,kk
   write(6,400) cc

   stop 
  
103 format(a,f10.5)  
100 format('ii= ',i2,' ff= ', f10.4)
200 format('ii ', i2)
300 format('ii= ',i2, ' jj= ',i2,' k= ', i2)
400 format(a32)

end !program demo_callC

subroutine abc(jj)
  jj=jj*2
  return
end !subroutine abc
