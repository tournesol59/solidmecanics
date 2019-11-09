!**************************************************************************!
!                                                                          !
! module verifymeshround                                                   !
!                                                                          !
!      Funktion: Verify that triangles are good triangles (not aligned elements)
!                                                                          !
!**************************************************************************!
!                                                                          !
!      trigelementcheck        Pruefen ein partikular Element              !
!      trigrotationxcheck      Pruefen dass eine Rotation je nach X axis   !
!                                nach einem anderen Element von gleichen Form fuehrt !
!       trigrotationycheck      Pruefen dass eine Rotation je nach Y axis ..
!**************************************************************************!

module verifymeshround_mod
  use TypesRound

  private
    real      :: TrigRadianToDegree = 57.297

  public :: verifyprntvects, verifyprntmatrix, trigcopyvects, trigcalcbasis, &
            trigcharacteristics, trigorientation, trigrotationxcheck

  contains
!********************************************
! subroutine verifyprntvects
  subroutine verifyprntvects(ux,uy,uz, SubName)
  use typesround
  implicit none

! uebergegebene Variabels
    type(tNode),intent(in)              :: ux, uy, uz  ! linien-Vektoren: 
    character(len=30),intent(in)        :: SubName
! lokale Var
   integer                   :: i,j    
   ! print basis

 write(*,209) "----- ", adjustl(SubName)
    do j=1,3
       select case (j)
         case (1)
           write(*,309) '     local basis input ux  j =',j
           write(*,109) (ux%vects(i), i=1,3)
         case (2)
           write(*,309) '     local basis input uy  j =',j
           write(*,109) (uy%vects(i), i=1,3)
         case (3)
           write(*,309) '     local basis input uz  j =',j
           write(*,109) (uz%vects(i), i=1,3)
        end select
    enddo

109 format(3f10.5)
209 format(a,a)
309 format(a,i10)
  end subroutine verifyprntvects

!********************************************
! subroutine verifyprntmatrix
  subroutine verifyprntmatrix(axout_xyz, SubName)
  use typesround
  implicit none

! uebergegebene Variabels
    type(tPassageMatrix),intent(in)     :: axout_xyz  ! linien-Vektoren: 
    character(len=30),intent(in)        :: SubName
! lokale Var
   integer                   :: i,j    
   ! print basis

 write(*,210) "----- ", adjustl(SubName)
    do i=1,3
      write(*,310) '     matrix line (row major)  i =',i
      write(*,110) (axout_xyz%mat(3*(j-1)+i), j=1,3)
    enddo
110 format(3f10.5)
210 format(a,a)
310 format(a,i10)
  end subroutine verifyprntmatrix

!********************************************
! subroutine trigcopyvects
!    simply copy the Passagematrix into three points,
!  whose xi-coordinates are into a vector 
!********************************************
  subroutine trigcopyvects(axin_xyz,ux,uy,uz)
  use typesround
  implicit none

! uebergegebene Variabels
    type(tNode)              :: ux, uy, uz  ! linien-Vektoren: ux hat nur x-achse Werten.
    type(tPassageMatrix)     :: axin_xyz
! Vorhabe
    intent(inout)            :: ux, uy, uz
    intent(in)               :: axin_xyz
! lokale Variable
    integer                  :: i,j

  do i=1,3
    do j=1,3
       select case (i)
         case (1)
           ux%vects(j)=axin_xyz%mat(3*(j-1)+1)  ! ux%vects(j) enthaelt x-achse Wert der j-te Column Vektor in mat
         case (2)
           uy%vects(j)=axin_xyz%mat(3*(j-1)+2)  ! da im Fortran der Column index 'oben' steht (j*dim), ist
         case (3)
           uz%vects(j)=axin_xyz%mat(3*(j-1)+3)  ! der Wert auch in mat(3*j+i) enthalten.
       end select
    enddo
  enddo

  end subroutine trigcopyvects

!********************************************
! subroutine trigcalcbasis
!    calculate the local basis linked to the 
!    triangle plane, which is defined by ux,uy,uz
!    and return this basis in Matrix format,
!    (coordinates of axis in canonical basis)
!*******************************************
subroutine trigcalcbasis(ux,uy,uz,axout_xyz)
  use typesround
  implicit none

! uebergegebene Variabels
    type(tNode)              :: ux, uy, uz  ! linien-Vektoren: ux hat nur x-achse Werten.
    type(tPassageMatrix)     :: axout_xyz
! Vorhabe
    intent(in)               :: ux, uy, uz
    intent(inout)            :: axout_xyz
! lokale Variable
    character(len=30)        :: Subname=" ----------trigcalcbasis       "
    integer                  :: i,j
    real                     :: unorm
    real,dimension(3)        :: cross3
    real,dimension(3)        :: cross2
    integer,dimension(1:3)     :: perm1=(/1, 1, 2 /)
    integer,dimension(1:3)     :: perm2=(/2, 3, 3 /)

  call verifyprntvects(ux, uy, uz, Subname)
  do i=1,3
    select case(i) 
      case(1)
        do j=1,2
           axout_xyz%mat(3*(j-1)+i) = ux%vects(perm2(j))-ux%vects(perm1(j)) ! j=1: edge(2-1), j=2: edge(3-1), j=3: no need
        enddo
      case(2)
        do j=1,2
          axout_xyz%mat(3*(j-1)+i) = uy%vects(perm2(j))-uy%vects(perm1(j)) ! j=1: edge(2-1), j=2: edge(3-1), j=3: no need
        enddo
      case(3)
        do j=1,2
          axout_xyz%mat(3*(j-1)+i) = uz%vects(perm2(j))-uz%vects(perm1(j)) ! j=1: edge(2-1), j=2: edge(3-1), j=3: no need
        enddo
    end select
  enddo
  axout_xyz%mat(3*2+1) = 0.0 ! just fill not hazardous value in third vector that will be re-calculated just below
  axout_xyz%mat(3*2+2) = 1.0
  axout_xyz%mat(3*2+3) = 1.0
  do j=1,3
    unorm = 0.0
    do i=1,3
      unorm = unorm + axout_xyz%mat(3*(j-1)+i)**2
    enddo
    unorm = sqrt(unorm)
    do i=1,3
      axout_xyz%mat(3*(j-1)+i) = axout_xyz%mat(3*(j-1)+i) / unorm
    enddo
  enddo

  ! mat(.,3)= mat(.,1) x mat(.,2) : (vector product)
  ! first calc: uj(k)-uj(l) is already computed and stored in axout_xyz%mat
  ! mat(1,3)= (u2(2)-u2(1))*(u3(3)-u3(1)) - (u3(2)-u3(1))*(u2(3)-u2(1))
  cross3(1) = (axout_xyz%mat(3*0+2)*axout_xyz%mat(3*1+3)) -  &
                        (axout_xyz%mat(3*1+2)*axout_xyz%mat(3*0+3))
  ! mat(2,3)= (u3(2)-u3(1))*(u1(3)-u1(1)) - (u1(2)-u1(1))*(u3(3)-u1(1))
  cross3(2) = (axout_xyz%mat(3*0+3)*axout_xyz%mat(3*1+1)) -  &
                        (axout_xyz%mat(3*0+1)*axout_xyz%mat(3*1+3))
  ! mat(3,3)= (u1(2)-u1(1))*(u2(3)-u2(1)) - (u2(2)-u2(1))*(u1(3)-u1(1))
  cross3(3) = (axout_xyz%mat(3*0+1)*axout_xyz%mat(3*1+2)) -  &
                        (axout_xyz%mat(3*0+2)*axout_xyz%mat(3*1+1))
 ! then recopy:
  axout_xyz%mat(3*2+1) = cross3(1)
  axout_xyz%mat(3*2+2) = cross3(2)
  axout_xyz%mat(3*2+3) = cross3(3)
  unorm = 0.0
  do i=1,3
    unorm = unorm + axout_xyz%mat(3*2+i)**2
  enddo
  unorm = sqrt(unorm)
  do i=1,3
      axout_xyz%mat(3*2+i) = axout_xyz%mat(3*2+i) / unorm
  enddo

! then re-calculate mat(2,.) to have orthonormal basis:
  ! first calc:
  ! mat(.,2)= mat(.,3) x mat(.,1) :

  ! mat(1,2)= mat(2,3)*mat(3,1) - mat(3,3)*mat(2,1)
  cross2(1) = (axout_xyz%mat(3*2+2)*axout_xyz%mat(3*0+3)) -  &
                        (axout_xyz%mat(3*2+3)*axout_xyz%mat(3*0+2))
  ! mat(2,2)= mat(3,3)*mat(1,1) - mat(1,3)*mat(3,1)
  cross2(2) = (axout_xyz%mat(3*2+3)*axout_xyz%mat(3*0+1)) -  &
                        (axout_xyz%mat(3*2+1)*axout_xyz%mat(3*0+3))
  ! mat(3,2)= mat(1,3)*mat(2,1) - mat(2,3)*mat(1,1)
  cross2(3) = (axout_xyz%mat(3*2+1)*axout_xyz%mat(3*0+2)) -  &
                        (axout_xyz%mat(3*2+2)*axout_xyz%mat(3*0+1))
 ! then recopy:
  axout_xyz%mat(3*1+1) = cross2(1)
  axout_xyz%mat(3*1+2) = cross2(2)
  axout_xyz%mat(3*1+3) = cross2(3)
  unorm = 0.0  !! end below not needed but for coherence
  do i=1,3
    unorm = unorm + axout_xyz%mat(3*1+i)**2
  enddo
  unorm = sqrt(unorm)
  do i=1,3
      axout_xyz%mat(3*1+i) = axout_xyz%mat(3*1+i) / unorm
  enddo
end subroutine trigcalcbasis

!********************************************
! subroutine trigcharacteristics
!    calculate the vertices and angles defining the triangle
!    and the triplet (x2, x3, y3) in the element-base system 
!    of axes (x,y), origin is at point 1
!*******************************************
subroutine trigcharacteristics(ux,uy,uz, dist, angle, triplet, area)
   use typesround
   implicit none

! uebergegebene Variabels
    type(tNode)              :: ux, uy, uz
    real,dimension(1:3)      :: dist
    real,dimension(1:3)      :: angle
    real,dimension(1:3)      :: angles_ref
    real,dimension(1:3)      :: triplet
    real,dimension(1:4)      :: area
! Vorhabe
    intent(in)               :: ux, uy, uz
    intent(inout)            :: dist, angle, triplet, area
! lokale Variable
    integer                  :: i,j
    character(len=30)        :: Subname=" ----------trigcharacteristics "

    call verifyprntvects(ux,uy,uz,Subname);
! calc L1, L2, L3:
    dist(1) = sqrt( (ux%vects(3) - ux%vects(2))*(ux%vects(3) - ux%vects(2)) + &
                    (uy%vects(3) - uy%vects(2))*(uy%vects(3) - uy%vects(2)) + &
                    (uz%vects(3) - uz%vects(2))*(uz%vects(3) - uz%vects(2)) )

    dist(2) = sqrt( (ux%vects(3) - ux%vects(1))*(ux%vects(3) - ux%vects(1)) + &
                    (uy%vects(3) - uy%vects(1))*(uy%vects(3) - uy%vects(1)) + &
                    (uz%vects(3) - uz%vects(1))*(uz%vects(3) - uz%vects(1)) ) 

    dist(3) = sqrt( (ux%vects(2) - ux%vects(1))*(ux%vects(2) - ux%vects(1)) + &
                    (uy%vects(2) - uy%vects(1))*(uy%vects(2) - uy%vects(1)) + &
                    (uz%vects(2) - uz%vects(1))*(uz%vects(2) - uz%vects(1)) )

    triplet(1) = dist(3) ! x2
    triplet(2) = ( (ux%vects(3) - ux%vects(1))*(ux%vects(2) - ux%vects(1)) + &
                   (uy%vects(3) - uy%vects(1))*(uy%vects(2) - uy%vects(1)) + &
                   (uz%vects(3) - uz%vects(1))*(uz%vects(2) - uz%vects(1)) ) / ( dist(3) )  ! x3= v31*v21/d3, d3=norm(v21)
    triplet(3) =  abs(( (ux%vects(2) - ux%vects(1))*(uy%vects(3) - uy%vects(1)) - &
                      (uy%vects(2) - uy%vects(1))*(ux%vects(3) - ux%vects(1))) / ( dist(3) ))  ! y3= (v31 X v21).nz/d3, d3=norm(v21)

 !!!   angle(1) = 57.297 * asin(triplet(3)/dist(2))  ! arcsinus(y3/d2) d2 length of segment opposite to pnt p2
    angle(3) = TrigRadianToDegree * ( atan(triplet(2)/triplet(3)) + atan((triplet(1)-triplet(2))/triplet(3)) )  
    angle(2) = TrigRadianToDegree * ( atan(triplet(3)/(triplet(1) - triplet(2))) )
    angle(1) = 180 - angle(3) - angle(2)
! AIRES A REVERIFER AVEC UNE FIGURE, area(i)=area of small triangle defining the center of gravity, and opposite to vertex i!
    area(4) =  triplet(3)*triplet(1)/2  ! x2*y3/2 
    area(3) = abs(triplet(3)/3*triplet(1)/2)
    area(2) = abs(triplet(3)*triplet(1)/4) - area(3)/2
    area(1)=area(4)-area(3)-area(2)
    !!!OTHER MEANS:
    !!!area(1) y3*(1/3*x3-2/3*(x3/2+x2/2)) - (x3-x2)*(1/3*y3-2/3*0) | 
    !!!area(3) y3*(1/3*x3-2/3*(x3/2+0)) - (x3-0)*(1/3*y3-2/3*0)  

   write(*, 103) "--------- Characteristics von Triangle Mesh: --------- "
   do i= 1,3 
     write(*, 503) "vertex: ", i, "an. di. ar.", angle(i), dist(i), area(i)
   enddo 
   i=3
   write(*, 403) "area tot: ", i, area(4)
   write(*, 503) "number ", i, "x2, x3, y3", triplet(1), triplet(2), triplet(3)


 103  format (a)
! 203  format (a,a)
! 303  format (f10.5)
 403  format (a,i10, f10.5)
 503  format(a, i10,a,f10.5,f10.5,f10.5)
end subroutine trigcharacteristics

!*************************************************************!
! subroutine trigorientation
!          soll vertices ueber Referenz (1,0,0;0,1,0;0,0,1) 
!          projecktzieren
!*************************************************************!
  subroutine trigorientation(ux,uy,uz,axout_xyz)  !, vect)
   use typesround
   use spacegeometry_mod
   implicit none

! uebergegebene Variabels
    type(tNode)              :: ux, uy, uz
!    real,dimension(1:3)      :: dist
    type(tPassageMatrix)     :: axout_xyz
  !  type(tNode)              :: vect
! Vorhabe
    intent(in)               :: ux, uy, uz  ! ,dist
    intent(inout)            :: axout_xyz    !vect
! lokale Variable
    type(tNode)              :: ux1, uy1, uz1
    integer                  :: i,j,allocStat
    real                     :: psi, theta, phi, errorint
    type(tPassageMatrix)     :: axinta_xyz, axintb_xyz, axend_xyz
    character(len=30)        :: Subname=" ---------trigorientation    "

    call trigcalcbasis(ux,uy,uz,axout_xyz)  ! calculate a local orthonormal basis

   ! print basis
    call verifyprntmatrix(axout_xyz, Subname)

! identify Euler angles:
    theta = acos(axout_xyz%mat(9))
    psi = asin(axout_xyz%mat(7)/sin(theta))
    phi = asin(axout_xyz%mat(7)/sin(theta))

    theta = TrigRadianToDegree * theta  ! nutation q
    psi = TrigRadianToDegree * psi      ! presession y
    phi = TrigRadianToDegree * phi      ! rotation j
! 
    print *, 'Terminated test trigorientation: OK'
! verify that the Euler rotations operated in inverse 
! order to the local orthonormal basis creates an identity matrix:
    allocate(ux1%vects(3), uy1%vects(3), uz1%vects(3), STAT=allocStat) 
    if (allocStat.NE.0) then
      print *, 'ERROR Allocate ux1 uy1 uz1: Could not allocate all variables!'
      STOP
    end if
    call trigcopyvects(axout_xyz,ux1,uy1,uz1) ! copy resulting matrix in 3-nodes points (vector)
    call trigoprotation(3,ux1,uy1,uz1,-psi,axinta_xyz) ! multiply by passage matrix

    call verifyprntmatrix(axinta_xyz, Subname)
    call trigcopyvects(axinta_xyz,ux1,uy1,uz1) ! copy resulting matrix in 3-nodes points (vector)
    call trigoprotation(1,ux1,uy1,uz1,-theta,axintb_xyz)  ! hence again

    call verifyprntmatrix(axintb_xyz, Subname)
    call trigcopyvects(axintb_xyz,ux1,uy1,uz1)  
    call trigoprotation(3,ux1,uy1,uz1,-phi,axend_xyz)
    ! print inverse operated basis

    call verifyprntmatrix(axend_xyz, Subname)

    deallocate(ux1%vects, uy1%vects, uz1%vects, STAT=allocStat) 
    if (allocStat.NE.0) then
      print *, 'ERROR Dellocate ux1 uy1 uz1: Could not deallocate all variables!'
      STOP
    end if
 102  format (3f10.5)
 202  format (a,f10.5)
 302  format (a,i10)
  end subroutine trigorientation

!***********************************************************!
! Subroutine CalcCrossProduct (endlich! muss Code oben/unten versetzten)
  subroutine CalcCrossProduct(ux, uy, uz, permute, crossvect)
  use typesround
  implicit none
! uebergegebene Variabels
  type(tNode)              :: ux, uy, uz, crossvect
  integer,dimension(3)     :: permute
! Vorhabe
  intent(in)               :: ux, uy, uz, permute
  intent(out)              :: crossvect
! lokale Variable
  integer                  :: i
 ! PAY ATTENTION: ux contains the x-coordinates of points 1,2 and 3 (indices in ux)
 ! and vice versa for uy,uz. The vectors have to be calculated by differenciation
 ! but the result crossvect indices 1,2,3 correspond to x,y,z values
  if ( (.not.(permute(1).eq.0)).and.(.not.(permute(2).eq.0)) ) then  ! Berechnet U(1-3)xU(2-3) Vektor
      crossvect%vects(1)= (uy%vects(1)-uy%vects(3))*(uz%vects(2)-uz%vects(3)) - &
                          (uz%vects(1)-uz%vects(3))*(uy%vects(2)-uy%vects(3))

      crossvect%vects(2)= (uz%vects(1)-uz%vects(3))*(ux%vects(2)-ux%vects(3)) - &
                          (ux%vects(1)-ux%vects(3))*(uz%vects(2)-uz%vects(3))

      crossvect%vects(3)= (ux%vects(1)-ux%vects(3))*(uy%vects(2)-uy%vects(3)) - &
                          (uy%vects(1)-uy%vects(3))*(ux%vects(2)-ux%vects(3))

   else if ( (.not.(permute(1).eq.0)).and.(.not.(permute(3).eq.0)) ) then  ! Berechnet U(1-2)xU(3-2) Vektor
      crossvect%vects(1)= (uy%vects(1)-uy%vects(2))*(uz%vects(3)-uz%vects(2)) - &
                          (uz%vects(1)-uz%vects(2))*(uy%vects(3)-uy%vects(3))

      crossvect%vects(2)= (uz%vects(1)-uz%vects(2))*(ux%vects(3)-ux%vects(2)) - &
                          (ux%vects(1)-ux%vects(2))*(uz%vects(3)-uz%vects(2))

      crossvect%vects(3)= (ux%vects(1)-ux%vects(2))*(uy%vects(3)-uy%vects(2)) - &
                          (uy%vects(1)-uy%vects(2))*(ux%vects(3)-ux%vects(2))

  else if ( (.not.(permute(2).eq.0)).and.(.not.(permute(3).eq.0)) ) then  ! Berechnet U(2-1)xU(3-1) Vektor
      crossvect%vects(1)= (uy%vects(2)-uy%vects(1))*(uz%vects(3)-uz%vects(1)) - &
                          (uz%vects(2)-uz%vects(1))*(uy%vects(3)-uy%vects(1))

      crossvect%vects(2)= (uz%vects(2)-uz%vects(1))*(ux%vects(3)-ux%vects(1)) - &
                          (ux%vects(2)-ux%vects(1))*(uz%vects(3)-uz%vects(1))

      crossvect%vects(3)= (ux%vects(2)-ux%vects(1))*(uy%vects(3)-uy%vects(1)) - &
                          (uy%vects(2)-uy%vects(1))*(ux%vects(3)-ux%vects(1))
  endif

end subroutine CalcCrossProduct

!***********************************************************!
! Subroutine CalcDotProduct (endlich! muss Code oben/unten versetzten)
! if not ambiguous, index i select the basis vector in axref_xyz to be multiplied by U1
  subroutine CalcDotProduct( U1, axref_xyz, i, val)
  use typesround
  implicit none
! uebergegebene Variable
  type(tNode)    :: U1
  type(tPassageMatrix) :: axref_xyz
  integer        :: i
  real           :: val
! Vorhabe
  intent(in)     :: U1, axref_xyz, i
  intent(out)    :: val
! locale

  val = U1%vects(1)*axref_xyz%mat(i) + &
        U1%vects(2)*axref_xyz%mat(3+i) + &
        U1%vects(3)*axref_xyz%mat(6+i)  


  end subroutine CalcDotProduct

!*************************************************************!
! subroutine trigrotationxcheck
!           soll 3x3 Matrizen ausnuetzen
!           computes the Passagematrix and two angles from the data nodes and ref axis
!           stores the results (local permutation indices, matrix and angles) in a structure
!           of type(tRoundlocal), which must be passed from/to the rigid finite Element 
!           computation routine.
!*************************************************************!
  subroutine trigrotationxcheck(ne, MeshR, axref_xyz, Varlocal)
    use typesround
    implicit none
! uebergegebene Variabels
    integer                  :: ne     ! Indiz to select
    type(tRoundMesh)         :: MeshR       ! Gitterwerte  
    type(tPassageMatrix)     :: axref_xyz   ! ref input, output
    type(tRoundlocal)        :: Varlocal
!    type(tNode),dimension(1:3) :: vects  ! rotated element to output
! Vorhabe
    intent(in)               :: ne, MeshR, axref_xyz
    intent(inout)            :: Varlocal
! lokale Variable
    integer                  :: i,j,k,l
    real                     :: coth, sith, coph, siph, crossval
    type(tPassageMatrix)     :: axint_xyz, axout_xyz
    type(tNode)              :: ux, uy, uz, ucopy, crossvect
    real,dimension(1:3)      :: dist
    real,dimension(1:3)      :: angle
    real,dimension(1:3)      :: triplet
    real,dimension(1:4)      :: area


   j=MeshR%elems(ne)%numeros(2)
   k=MeshR%elems(ne)%numeros(3)
   l=MeshR%elems(ne)%numeros(4) ! copy Nodes coordinates in three vectors ux,uy,uz, ith-component owns ith-node 
   do i=1,3
     select case (i)
       case (1)           !! There is an INCORRECT assumption below, that the nodes are ordered in 1,2,3 !!
         ux%vects(1)= MeshR%nodes(j)%vects(1)
         ux%vects(2)= MeshR%nodes(k)%vects(1)
         ux%vects(3)= MeshR%nodes(l)%vects(1)
       case (2)
         uy%vects(1)= MeshR%nodes(j)%vects(2)
         uy%vects(2)= MeshR%nodes(k)%vects(2)
         uy%vects(3)= MeshR%nodes(l)%vects(2)
       case (3)
         uz%vects(1)= MeshR%nodes(j)%vects(3)
         uz%vects(2)= MeshR%nodes(k)%vects(3)
         uz%vects(3)= MeshR%nodes(l)%vects(3)
    end select
   enddo
   ! Setzt die Ordnung von Punkten mit Kenntnis von referenz Z-Achsis !! Hum Hum (muss besser kunftig)
   ! Voraussetzt, dass Varlocal%permute schon partial erfuellt ist

!   if ( (.not.(Varlocal%permute(1).eq.0)).and.(.not.(Varlocal%permute(2).eq.0)).and.(.not.(Varlocal%permute(3).eq.0)) )
!      write (*,*)
!      print *, "At least one node of Varlocal%permute has to be equal to zero (uninitialized)" 
!   endif

   call CalcCrossProduct(ux, uy, uz, Varlocal%permute, crossvect)
 ! permute definiert, welche Vektor mit welchem anderen zum CrossProduct gemacht wird
   call CalcDotProduct( crossvect, axref_xyz, 3, crossval)  ! DotProduct mit Vektor 3 von axref_xyz

  ! ********************************************************************************************
  ! ******** Folgende wird kompliziert Code aber zum sehr einfachen Ziel: permutation der Punkten
   if ( (.not.(Varlocal%permute(1).eq.0)).and.(.not.(Varlocal%permute(2).eq.0)) ) then
      if (0.le.crossval) then ! do fast nothing, e.g. if permute(1)=b and permute(2)=c then permute(3)=a
                         ! that is do not change b,c
         if ( (Varlocal%permute(1).eq.1).and.(Varlocal%permute(2).eq.2) ) then
            Varlocal%permute(3)=3
         elseif ( (Varlocal%permute(1).eq.1).and.(Varlocal%permute(2).eq.3) ) then
            Varlocal%permute(3)=2
         else
            Varlocal%permute(3)=1
         endif
      else ! intervertice permute(1) and permute(2)
         if ( (Varlocal%permute(1).eq.1).and.(Varlocal%permute(2).eq.2) ) then
            Varlocal%permute(1)=2
            Varlocal%permute(2)=1
            Varlocal%permute(3)=3
         elseif ( (Varlocal%permute(1).eq.1).and.(Varlocal%permute(2).eq.3) ) then
            Varlocal%permute(1)=3
            Varlocal%permute(2)=1
            Varlocal%permute(3)=2
         else
            Varlocal%permute(1)=3
            Varlocal%permute(2)=2 !! found Varlocal%permute(3)=2 corrected:Varlocal%permute(2)=2 , ERROR?
            Varlocal%permute(3)=1
         endif
      endif
   endif
   if ( (.not.(Varlocal%permute(1).eq.0)).and.(.not.(Varlocal%permute(3).eq.0)) ) then
      if (0.le.crossval) then ! do fast nothing, e.g. if permute(1)=b and permute(2)=c then permute(3)=a
                         ! that is do not change b,c
         if ( (Varlocal%permute(1).eq.1).and.(Varlocal%permute(3).eq.2) ) then
            Varlocal%permute(2)=3
         elseif ( (Varlocal%permute(1).eq.1).and.(Varlocal%permute(3).eq.3) ) then
            Varlocal%permute(2)=2
         else
            Varlocal%permute(2)=1
         endif
      else ! intervertice permute(1) and permute(2)
         if ( (Varlocal%permute(1).eq.1).and.(Varlocal%permute(3).eq.2) ) then
            Varlocal%permute(1)=2
            Varlocal%permute(3)=1
            Varlocal%permute(2)=3
         elseif ( (Varlocal%permute(1).eq.1).and.(Varlocal%permute(3).eq.3) ) then
            Varlocal%permute(1)=3
            Varlocal%permute(3)=1
            Varlocal%permute(2)=2
         else
            Varlocal%permute(1)=3
            Varlocal%permute(3)=2
            Varlocal%permute(2)=1
         endif
       endif
   endif
   if ( (.not.(Varlocal%permute(2).eq.0)).and.(.not.(Varlocal%permute(3).eq.0)) ) then
      if (0.le.crossval) then ! do fast nothing, e.g. if permute(1)=b and permute(2)=c then permute(3)=a
                         ! that is do not change b,c
         if ( (Varlocal%permute(2).eq.1).and.(Varlocal%permute(3).eq.2) ) then
            Varlocal%permute(1)=3
         elseif ( (Varlocal%permute(1).eq.1).and.(Varlocal%permute(3).eq.3) ) then
            Varlocal%permute(1)=2
         else
            Varlocal%permute(1)=1
         endif
      else ! intervertice permute(1) and permute(2)
         if ( (Varlocal%permute(2).eq.1).and.(Varlocal%permute(3).eq.2) ) then
            Varlocal%permute(2)=2
            Varlocal%permute(3)=1
            Varlocal%permute(1)=3
         elseif ( (Varlocal%permute(2).eq.1).and.(Varlocal%permute(3).eq.3) ) then
            Varlocal%permute(2)=3
            Varlocal%permute(3)=1
            Varlocal%permute(1)=2
         else
            Varlocal%permute(2)=3
            Varlocal%permute(3)=2
            Varlocal%permute(1)=1
         endif
     endif
   endif
  ! Endlich
  ! ********************************************************************************************
  ! Noch mal komplizierten Kode zum einfachen Permutation (dieses Mal wirklich der Koordinaten):
  ! **************
   !if ((Varlocal%permute(1).eq.1).and.(Varlocal%permute(2).eq.2).and.(Varlocal%permute(3).eq.3)) then
     ! do nothing 
   !elseif
   if ((Varlocal%permute(1).eq.2).and.(Varlocal%permute(2).eq.1).and.(Varlocal%permute(3).eq.3)) then
     ! do something
       ucopy%vects(1)= ux%vects(1)
       ucopy%vects(2)= uy%vects(1)
       ucopy%vects(3)= uz%vects(1)

       ux%vects(1)= ux%vects(2)
       uy%vects(1)= uy%vects(2)
       uz%vects(1)= uz%vects(2)
 
       ux%vects(2)= ucopy%vects(1)
       uy%vects(2)= ucopy%vects(2)
       uz%vects(2)= ucopy%vects(3)

   elseif ((Varlocal%permute(1).eq.3).and.(Varlocal%permute(2).eq.2).and.(Varlocal%permute(3).eq.1)) then
     ! do something
       ucopy%vects(1)= ux%vects(1)
       ucopy%vects(2)= uy%vects(1)
       ucopy%vects(3)= uz%vects(1)

       ux%vects(1)= ux%vects(3)
       uy%vects(1)= uy%vects(3)
       uz%vects(1)= uz%vects(3)
 
       ux%vects(3)= ucopy%vects(1)
       uy%vects(3)= ucopy%vects(2)
       uz%vects(3)= ucopy%vects(3)

   elseif ((Varlocal%permute(1).eq.1).and.(Varlocal%permute(2).eq.3).and.(Varlocal%permute(3).eq.2)) then
     ! do something
       ucopy%vects(1)= ux%vects(2)
       ucopy%vects(2)= uy%vects(2)
       ucopy%vects(3)= uz%vects(2)

       ux%vects(2)= ux%vects(3)
       uy%vects(2)= uy%vects(3)
       uz%vects(2)= uz%vects(3)
 
       ux%vects(3)= ucopy%vects(1)
       uy%vects(3)= ucopy%vects(2)
       uz%vects(3)= ucopy%vects(3)

   elseif ((Varlocal%permute(1).eq.2).and.(Varlocal%permute(2).eq.3).and.(Varlocal%permute(3).eq.1)) then
     ! do something (1)
       ucopy%vects(1)= ux%vects(1)
       ucopy%vects(2)= uy%vects(1)
       ucopy%vects(3)= uz%vects(1)

       ux%vects(1)= ux%vects(2)
       uy%vects(1)= uy%vects(2)
       uz%vects(1)= uz%vects(2)
 
       ux%vects(2)= ucopy%vects(1)
       uy%vects(2)= ucopy%vects(2)
       uz%vects(2)= ucopy%vects(3)
     ! do something (2)
       ucopy%vects(1)= ux%vects(2) ! ie ux%vects(1)..
       ucopy%vects(2)= uy%vects(2)
       ucopy%vects(3)= uz%vects(2)

       ux%vects(2)= ux%vects(3)
       uy%vects(2)= uy%vects(3)
       uz%vects(2)= uz%vects(3)
 
       ux%vects(3)= ucopy%vects(1)
       uy%vects(3)= ucopy%vects(2)
       uz%vects(3)= ucopy%vects(3)

   elseif ((Varlocal%permute(1).eq.3).and.(Varlocal%permute(2).eq.1).and.(Varlocal%permute(3).eq.2)) then
     ! do something(1)
       ucopy%vects(1)= ux%vects(1)
       ucopy%vects(2)= uy%vects(1)
       ucopy%vects(3)= uz%vects(1)

       ux%vects(1)= ux%vects(3)
       uy%vects(1)= uy%vects(3)
       uz%vects(1)= uz%vects(3)
 
       ux%vects(3)= ucopy%vects(1)
       uy%vects(3)= ucopy%vects(2)
       uz%vects(3)= ucopy%vects(3)
     ! do something(2)
       ucopy%vects(1)= ux%vects(2) 
       ucopy%vects(2)= uy%vects(2)
       ucopy%vects(3)= uz%vects(2)

       ux%vects(2)= ux%vects(3)
       uy%vects(2)= uy%vects(3)
       uz%vects(2)= uz%vects(3)
 
       ux%vects(3)= ucopy%vects(1)
       uy%vects(3)= ucopy%vects(2)
       uz%vects(3)= ucopy%vects(3)
  endif
  ! Endlich
  ! ********************************************************************************************
  ! Folgendes macht der prinzipiell Arbeit dieser (..) Subroutine trigrotationxcheck:

   call trigcharacteristics(ux,uy,uz, dist, angle, triplet, area) ! calc Daten.. (siehe sub)

   call trigorientation(ux,uy,uz, axout_xyz) ! calc lokal Basis (einfach.. siehe sub)

  end subroutine trigrotationxcheck

end module verifymeshround_mod
