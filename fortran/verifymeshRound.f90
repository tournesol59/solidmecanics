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

  public :: trigcopyvects, trigcalcbasis, trigcharacteristics, &
             trigorientation, trigrotationxcheck

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
    integer,dimension(1:3)     :: perm1=(/1, 1, 2 /)
    integer,dimension(1:3)     :: perm2=(/2, 3, 3 /)

  do i=1,3
   unorm = 0.0
    do j=1,3
       select case (i)
         case (1)
           unorm= unorm + ( ux%vects(perm2(j))-ux%vects(perm1(j)) )* &
                        ( ux%vects(perm2(j))-ux%vects(perm1(j)) )
         case (2)
           unorm= unorm + ( uy%vects(perm2(j))-uy%vects(perm1(j)) )* &
                        ( uy%vects(perm2(j))-uy%vects(perm1(j)) )
         case (3)
           unorm= unorm + ( uz%vects(perm2(j))-uz%vects(perm1(j)) )* &
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
  axout_xyz%mat(3*2+1)= (axout_xyz%mat(3*1+2)*axout_xyz%mat(3*2+3)) -  &
                        (axout_xyz%mat(3*2+2)*axout_xyz%mat(3*1+3))
  ! mat(2,3)= (u3(2)-u3(1))*(u1(3)-u1(1)) - (u1(2)-u1(1))*(u3(3)-u1(1))
  axout_xyz%mat(3*2+2)= (axout_xyz%mat(3*2+2)*axout_xyz%mat(3*0+3)) -  &
                        (axout_xyz%mat(3*0+2)*axout_xyz%mat(3*2+3))
  ! mat(1,3)= (u1(2)-u1(1))*(u2(3)-u2(1)) - (u2(2)-u2(1))*(u1(3)-u1(1))
  axout_xyz%mat(3*2+3)= (axout_xyz%mat(3*0+2)*axout_xyz%mat(3*1+3)) -  &
                        (axout_xyz%mat(3*1+2)*axout_xyz%mat(3*0+3))

! then re-calculate mat(2,.) to have orthonormal basis:
  ! mat(.,2)= mat(.,3) x mat(.,1) :

  ! mat(1,2)= mat(2,3)*mat(3,1) - mat(3,3)*mat(2,1)
  axout_xyz%mat(3*1+1)= (axout_xyz%mat(3*2+2)*axout_xyz%mat(3*0+3)) -  &
                        (axout_xyz%mat(3*2+3)*axout_xyz%mat(3*0+2))
  ! mat(2,2)= mat(3,3)*mat(1,1) - mat(1,3)*mat(3,1)
  axout_xyz%mat(3*1+2)= (axout_xyz%mat(3*2+3)*axout_xyz%mat(3*0+1)) -  &
                        (axout_xyz%mat(3*2+1)*axout_xyz%mat(3*0+3))
  ! mat(3,2)= mat(1,3)*mat(2,1) - mat(2,3)*mat(1,1)
  axout_xyz%mat(3*1+3)= (axout_xyz%mat(3*2+1)*axout_xyz%mat(3*0+2)) -  &
                        (axout_xyz%mat(3*2+2)*axout_xyz%mat(3*0+1))


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
   print *,""
   write(*,503) "    Check input data first vector: ",1," : ",ux%vects(1),uy%vects(1),uz%vects(1)
   write(*,503) "    Check input data secnd vector: ",2," : ",ux%vects(2),uy%vects(2),uz%vects(2)
   write(*,503) "    Check input data third vector: ",3," : ",ux%vects(3),uy%vects(3),uz%vects(3)

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
    angle(3) = 57.297 * ( atan(triplet(2)/triplet(3)) + atan((triplet(1)-triplet(2))/triplet(3)) )  
    angle(2) = 57.297 * ( atan(triplet(3)/(triplet(1) - triplet(2))) )
    angle(1) = 180 - angle(3) - angle(2)
! AIRES A REVERIFER AVEC UNE FIGURE, area(i)=area of small triangle defining the center of gravity, and opposite to vertex i!
    area(4) =  triplet(3)*triplet(1)/2  ! x2*y3/2 
    area(3) = abs(triplet(3)/3*triplet(1)/2)
    area(2)=abs(triplet(3)*triplet(1)/4) - area(3)/2
    area(1)=area(4)-area(3)-area(2)
    !!!area(1) = abs( triplet(3)*(1/3*triplet(2)-2/3*(triplet(2)/2+triplet(1)/2)) - &
    !!!              (triplet(2)-triplet(1))*1/3*triplet(3) )     ! | y3*(1/3*x3-2/3*(x3/2+x2/2)) - (x3-x2)*(1/3*y3-2/3*0) | 
    !!!area(3) = abs( triplet(3)*(1/3*triplet(2)-2/3*triplet(2)/2) - &
    !!!              (triplet(2))*1/3*triplet(3) )     ! | y3*(1/3*x3-2/3*(x3/2+0)) - (x3-0)*(1/3*y3-2/3*0)  

   print *, 'Started test trigcharacteristics: OK'
   write(*, 103) "--------- Characteristis von Triangle Mesh: --------- "
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
    phi = asin(axout_xyz%mat(7)/sin(theta))

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

  if ( (.not.(permute(1).eq.0)).and.(.not.(permute(2).eq.0)) ) then  ! Berechnet U(1-3)xU(2-3) Vektor
      crossvect%vects(1)= (uy%vects(1)-uy%vects(3))*(uz%vects(2)-uz%vects(3)) - &
                          (uz%vects(1)-uz%vects(3))*(uy%vects(2)-uy%vects(3))

      crossvect%vects(2)= (uz%vects(1)-uz%vects(3))*(ux%vects(2)-ux%vects(3)) - &
                          (ux%vects(1)-ux%vects(3))*(uz%vects(2)-uz%vects(3))

      crossvect%vects(3)= (ux%vects(1)-ux%vects(3))*(uy%vects(2)-uy%vects(3)) - &
                          (uy%vects(1)-uy%vects(3))*(ux%vects(2)-ux%vects(3))

   else if ( (.not.(permute(1).eq.0)).and.(.not.(permute(3).eq.0)) ) then  ! Berechnet U(1-3)xU(2-3) Vektor
      crossvect%vects(1)= (uy%vects(1)-uy%vects(2))*(uz%vects(3)-uz%vects(2)) - &
                          (uz%vects(1)-uz%vects(2))*(uy%vects(3)-uy%vects(3))

      crossvect%vects(2)= (uz%vects(1)-uz%vects(2))*(ux%vects(3)-ux%vects(2)) - &
                          (ux%vects(1)-ux%vects(2))*(uz%vects(3)-uz%vects(2))

      crossvect%vects(3)= (ux%vects(1)-ux%vects(2))*(uy%vects(3)-uy%vects(2)) - &
                          (uy%vects(1)-uy%vects(2))*(ux%vects(3)-ux%vects(2))

  else if ( (.not.(permute(2).eq.0)).and.(.not.(permute(3).eq.0)) ) then  ! Berechnet U(1-3)xU(2-3) Vektor
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
            Varlocal%permute(3)=2
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
