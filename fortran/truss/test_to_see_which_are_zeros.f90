PROGRAM test_to_see_which_are_zeros

IMPLICIT NONE

 integer     :: a,b,c,d

 real        :: fa, fb, fc, fd

! read(*,"(4I5)") a,b,c,d

! write(*,*) a,b,c,d

  OPEN(UNIT=25, FILE='Untitled.in', ACTION='READ')
 
  read(25, "(3f8.7)") fa, fb, fc
! read(*, "(e10.2)" fa fb, fc, fd
  write(*,105) fa, fb, fc
  read(25, "(3e10.2)") fa, fb, fc
  write(*,105) fa, fb, fc
  fd=fa+fb+fc
  write(*,106) fd

  CLOSE(UNIT=25) 
105 format(4e10.2)
106 format(f10.5) 
END program test_to_see_which_are_zeros
