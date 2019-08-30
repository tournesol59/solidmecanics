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

  public :: trigcopyvects, trigcalcbasis, trigcharacteristics, trigrotationxcheck

  contains
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
! subroutine trigcopybasis
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
    integer                  :: i,j
    real                     :: unorm
    integer,dimension(3)     :: perm1=/(1, 1, 2)
    integer,dimension(3)     :: perm2=/(2, 3, 3)

  do i=1,3
   unorm = 0.0
    do j=1,3
       select case (i)
         case (1)
           norm= norm + ( ux%vects(perm2(j))-ux%vects(perm1(j)) )* &
                        ( ux%vects(perm2(j))-ux%vects(perm1(j)) )
         case (2)
           norm= norm + ( uy%vects(perm2(j))-uy%vects(perm1(j)) )* &
                        ( uy%vects(perm2(j))-uy%vects(perm1(j)) )
         case (3)
           norm= norm + ( uz%vects(perm2(j))-uz%vects(perm1(j)) )* &
                        ( uz%vects(perm2(j))-uz%vects(perm1(j)) )
       end select
    enddo
    unorm = sqrt(unorm)
    do j=1,3                  ! mat(i,j)=(ui(j+1)-ui(j))/normj (Fortran index und 1=x, 2=y, 3=z)
       select case (i)
         case (1)
           axout_xyz%mat(3*(j-1)+i) = ( ux%vects(perm2(j))-ux%vects(perm1(j)) )/unorm
         case (2)
           axout_xyz%mat(3*(j-1)+i) = ( uy%vects(perm2(j))-uy%vects(perm1(j)) )/unorm
!         case (3)
! nichts: drittes Vektor aus Vektor-Product aus ersten Vektoren
       end select
    enddo
  enddo
  ! mat(.,3)= mat(.,1) x mat(.,2) :

  ! mat(1,3)= (u2(2)-u2(1))*(u3(3)-u3(1)) - (u3(2)-u3(1))*(u2(3)-u2(1))
  axout_xyz%mat(3*2+1)= (axout_xyz%mat(3*1+2)*axout_xyz%mat(3*2+3) -  &
                        (axout_xyz%mat(3*2+2)*(axout_xyz%mat(3*1+3)
  ! mat(2,3)= (u3(2)-u3(1))*(u1(3)-u1(1)) - (u1(2)-u1(1))*(u3(3)-u1(1))
  axout_xyz%mat(3*2+2)= (axout_xyz%mat(3*2+2)*axout_xyz%mat(3*0+3) -  &
                        (axout_xyz%mat(3*0+2)*(axout_xyz%mat(3*2+3)
  ! mat(1,3)= (u1(2)-u1(1))*(u2(3)-u2(1)) - (u2(2)-u2(1))*(u1(3)-u1(1))
  axout_xyz%mat(3*2+3)= (axout_xyz%mat(3*0+2)*axout_xyz%mat(3*1+3) -  &
                        (axout_xyz%mat(3*1+2)*(axout_xyz%mat(3*0+3) 

! then re-calculate mat(2,.) to have orthonormal basis:
  ! mat(.,2)= mat(.,3) x mat(.,1) :

  ! mat(1,2)= mat(2,3)*mat(3,1) - mat(3,3)*mat(2,1)
  axout_xyz%mat(3*1+1)= (axout_xyz%mat(3*2+2)*axout_xyz%mat(3*0+3) -  &
                        (axout_xyz%mat(3*2+3)*(axout_xyz%mat(3*0+2)
  ! mat(2,2)= mat(3,3)*mat(1,1) - mat(1,3)*mat(3,1)
  axout_xyz%mat(3*1+2)= (axout_xyz%mat(3*2+3)*axout_xyz%mat(3*0+1) -  &
                        (axout_xyz%mat(3*2+1)*(axout_xyz%mat(3*0+3)
  ! mat(3,2)= mat(1,3)*mat(2,1) - mat(2,3)*mat(1,1)
  axout_xyz%mat(3*1+3)= (axout_xyz%mat(3*2+1)*axout_xyz%mat(3*0+2) -  &
                        (axout_xyz%mat(3*2+2)*(axout_xyz%mat(3*0+1)


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
    triplet(3) =  ( (uy%vects(2) - uy%vects(1))*(uz%vects(3) - uz%vects(1))*(ux%vects(2) - ux%vects(1)) + &
                    (uz%vects(2) - uz%vects(1))*(ux%vects(3) - ux%vects(1))*(uy%vects(2) -uy%vects(1)) + &
                    (ux%vects(2) - ux%vects(1))*(uy%vects(2) - uy%vects(1))*(uz%vects(2) - uz%vects(1))  ) / ( dist(3) )  ! y3= (v31 X v21)/d3, d3=norm(v21)

    angle(1) = TrigRadianToDegree * asin(triplet(3)/dist(2))  ! arcsinus(y3/d2) d2 length of segment opposite to pnt p2
    angle(2) = TrigRadianToDegree * asin(triplet(3)/(triplet(1) - triplet(3)) ) 
    angle(3) = 180 - angle(1) - angle(2)

    area(4) =  triplet(3)*triplet(2)/2  ! x2*y3/2 
    area(1) = abs( triplet(3)*(1/3*triplet(2)-2/3*(triplet(2)/2+triplet(1)/2)) - &
                  (triplet(2)-triplet(1))*1/3*triplet(3) )     ! | y3*(1/3*x3-2/3*(x3/2+x2/2)) - (x3-x2)*(1/3*y3-2/3*0) | 
    area(2) = abs( triplet(3)*(1/3*triplet(2)-2/3*triplet(2)/2) - &
                  (triplet(2))*1/3*triplet(3) )     ! | y3*(1/3*x3-2/3*(x3/2+0)) - (x3-0)*(1/3*y3-2/3*0) | 
    area(3) = area(4) - area(1) - area(2)
end subroutine trigcharacteristics

!*************************************************************!
! subroutine trigorientation
!          soll vertices ueber Referenz (1,0,0;0,1,0;0,0,1) 
!          projecktzieren
!*************************************************************!
  subroutine trigorientation(ux,uy,uz, dist, axout_xyz, vect)
   use typesround, spacegeometry_mod
   implicit none

! uebergegebene Variabels
    type(tNode)              :: ux, uy, uz
    real,dimension(1:3)      :: dist
    type(tPassageMatrix)     :: axout_xyz
    type(tNode)              :: vect
! Vorhabe
    intent(in)               :: ux, uy, uz, dist
    intent(inout)            :: axout_xyz, vects
! lokale Variable
    type(tNode)              :: ux1, uy1, uz1
    integer                  :: i,j
    real                     :: psi, theta, phi, errorint
    type(tPassageMatrix)     :: axinta_xyz, axintb_xyz, axend_xyz

    call trigcalcbasis(ux,uy,uz,axout_xyz)  ! calculate a local orthonormal basis
    ! print basis
    do j=1,3
         write(*,302) '     local basis matrix line j   =',j
         write(*,102) (axout_xyz%mat(j*3+i), i=1,3)
    enddo

! identify Euler angles:
    theta = acos(axout_xyz%mat(9))
    psi = asin(axout_xyz%mat(7)/sin(theta))
    pĥi = asin(axout_xyz%mat(7)/sin(theta))

    theta = TrigRadianToDegree * theta  ! nutation q
    psi = TrigRadianToDegree * psi      ! presession y
    phi = TrigRadianToDegree * phi      ! rotation j
! 

! verify that the Euler rotations operated in inverse 
! order to the local orthonormal basis creates an identity matrix:
    call trigcopyvects(axout_xyz,ux1,uy1,uz1) ! copy resulting matrix in 3-nodes points (vector)
    call trigoprotation(3,ux1,uy1,uz1,-psi,axinta_xyz) ! multiply by passage matrix
    call trigcopyvects(axinta_xyz,ux1,uy1,uz1) ! copy resulting matrix in 3-nodes points (vector)
    call trigoprotation(1,ux1,uy1,uz1,-theta,axintb_xyz)  ! hence again
    call trigcopyvects(axintb_xyz,ux1,uy1,uz1)  
    call trigoprotation(3,ux1,uy1,uz1,-phi,axend_xyz)
    ! print inverse operated basis
    do j=1,3
         write(*,302) '     inverse operated matrix line j   =',j
         write(*,102) (axend_xyz%mat(j*3+i), i=1,3)
    enddo
 102  format (f10.5)
 202  format (a,f10.5)
 302  format (a,i10)
  end subroutine trigorientation

!*************************************************************!
! subroutine trigrotationxcheck
!           soll 3x3 Matrizen ausnuetzen
!           computes the Passagematrix and two angles from the data nodes and ref axis
!*************************************************************!
  subroutine trigrotationxcheck(ne, MeshR, axref_xyz, axout_xyz, vects)
    use typesround
    implicit none
! uebergegebene Variabels
    integer                  :: ne     ! Indiz to select
    type(tRoundMesh)         :: MeshR       ! Gitterwerte  
    type(tPassageMatrix)     :: axref_xyz, axout_xyz ! ref input, output
    type(tNode),dimension(1:3) :: vects  ! rotated element to output
! Vorhabe
    intent(in)               :: ne, MeshR, axref_xyz
    intent(inout)            :: axout_xyz
! lokale Variable
    integer                  :: i,j,k,l
    real                     :: coth, sith, coph, siph
    type(tPassageMatrix)     :: axint_xyz
    type(tNode)              :: ux, uy, uz
    real,dimension(1:3)      :: dist
    real,dimension(1:3)      :: angle
    real,dimension(1:3)      :: triplet
    real,dimension(1:4)      :: area


   j=MeshR%elems(ne)%numeros(1)
   k=MeshR%elems(ne)%numeros(2)
   l=MeshR%elems(ne)%numeros(3)! copy Nodes coordinates in three vectors ux,uy,uz, ith-component owns ith-node 
   do i=1,3
     select case i
       case (1)
         ux%vects(j)= MeshR%nodes(j)%vects(1)
       case (2)
         uy%vects(j)= MeshR%nodes(k)%vects(2)
       case (3)
         uz%vects(j)= MeshR%nodes(l)%vects(3)
    end select
   enddo
   call trigcharacteristics(ux,uy,uz, dist, angle, triplet, area)
   call trigorientation(ux,uy,uz, dist, axout_xyz,vects)

  end subroutine trigrotationxcheck

end module verifymeshround_mod
