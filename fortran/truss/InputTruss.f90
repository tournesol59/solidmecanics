module InputTruss_mod

 use TypesBalken
 implicit none

  private
!  interface Input_file
!    module procedure Input_file
!  end interface Input_file

!  interface InputSetUnknown
!    module procedure InputSetUnknown
!  end interface InputSetUnknown

  integer    :: DLIMIT=-1.00
  integer    :: FLIMIT=-1000.00

  public:: Input_file, InputSetUnknown

  contains
!******************************************************!
! read, parse Input file e.g. Untitled2.in
!******************************************************!
  subroutine Input_file(MeshInfo,MeshT,MeshGen,meshpattern,Dimposed,Fimposed)   !meshconnect

  use typesbalken

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tMeshInfo)            :: MeshInfo     ! Gitterwerte                   !
  type(tMeshCoord)           :: MeshT     ! Gitterwerte                      !
  type(tMeshElmt),pointer    :: MeshGen(:)
!  type(tMeshElmt),pointer   :: ElmtS(:)     ! Gitter Klassen fur jede Verbindung
  integer,pointer            :: meshpattern(:,:) ! Definition von Elementen !
!  integer,pointer           :: meshconnect(:,:) ! Definition von Verbindungen !
!  integer,pointer           :: Dunknowns(:,:)
  real,pointer               :: Dimposed(:,:)
!  integer,pointer           :: Fmovemt(:,:)
  real,pointer               :: Fimposed(:,:)
!  type(tExakt)              :: Exakt      ! Exakte Loesung                 !
!  type(tFileIO)             :: FileIO     ! Ausgabesteuerung               !
  !-------------------------------------------------------------------------!  
  character(LEN=13)          :: Mshfilename  ! Name of the input file (13 characters) in current dir !
  integer :: i,l,countmax !, nbound, cbound, bnddim, bndind ! counter integers
  character(LEN=6) :: bndname
  real             :: a, b, c, d
  OPEN(UNIT=25, FILE='Untitled.in', ACTION='READ')  !!!!!!!!!!! ONLY FOR DEBUG TEST
  i=1
  read(25, "(E6.2)") MeshT%x(i)
  write(*,"(E6.2)") MeshT%x(i)
  i=2
  read(25,"(E6.2)") a
  write(*,"(E6.2)") a
  i=3
  read(25, "(E6.2)") MeshT%x(i)
  write(*,"(E6.2)") MeshT%x(i)
  CLOSE(UNIT=25)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  OPEN(UNIT=25, FILE='Untitled3.in', ACTION='READ')
  read(25,*) ! #
  read(25,*) ! instead of
!!  read(25,307) MeshInfo%nn, MeshInfo%nt  ! already read
  read(25,205) MeshInfo%EY, MeshInfo%nu
  read(25,*) ! #coord nodes (nn)
  do i=1,MeshInfo%nn
   read(25,207) MeshT%x(i), MeshT%y(i), MeshT%z(i)   ! z Koordinate muss geschrieben werden, wird aber nicht benutzt
  enddo
  read(25,*) ! #
  do i=1,MeshInfo%nt
    read(25,310) meshpattern(i,1), meshpattern(i,2), meshpattern(i,3), MeshGen(i)%SArea, MeshGen(i)%CI ! Solely for SArea and CI a call to MeshGen(i) is done
  enddo
  read(25,*) ! # elmts consists in types of connection for element (nn,ne) see code TypesBalkDef.f90 to have a description of the code, 7 is not treated yet
  do i=1,MeshInfo%nt
    read(25,308) MeshT%elmts(i,1), MeshT%elmts(i,2), MeshT%elmts(i,3), MeshT%elmts(i,4), &
               MeshT%elmts(i,5), MeshT%elmts(i,6), MeshT%elmts(i,7)  ! 
  enddo
  read(25,*) ! #
  do i=1,MeshInfo%nn
    read(25,207) l, Dimposed(i,1), Dimposed(i,2), Dimposed(i,3)   ! l is not used
  enddo
  read(25,*) ! #
  do i=1,MeshInfo%nn
    read(25,207) l, Fimposed(i,1), Fimposed(i,2), Fimposed(i,3)  ! l is not used
  enddo

  do i=1,MeshInfo%nt
    read(25,206) l, MeshGen(i)%q1, MeshGen(i)%q2  ! l is not used, q1,q2 filled for all elemts also for non beam elmt
  enddo
  CLOSE(UNIT=25)

 105  format (a20)
 205  format (2f8.7)             
 206  format (i5,2f8.7)
 207  format (3f8.7)
 307  format (2i5)
 308  format (7i5)
 309  format (3i5)
 310  format (i5,i5,i5,f8.7,f8.7)

  end subroutine Input_file

!**********************************************************************************!
! InputSetUnknown: the idea behind is to have just a table (each line 
! an element rod or beam) with a single code number for whether displacement (resp. force)
! have to be determined or is imposed (resp. a force is provided at the node) and
! this procedure fills depending from that code number (integer) to booleans for each 
! triplet of dof/forces in a table per element
! 
!**********************************************************************************!
  subroutine InputSetUnknown(MeshInfo,MeshT,meshpattern,Dunknowns,Dimposed,Fmovemt,Fimposed)

  use typesbalken

  implicit none
  !                                                                          !
  type(tMeshInfo)           :: MeshInfo     ! Gitterwerte                    !
  type(tMeshCoord)           :: MeshT     ! Gitterwerte                      !
!  type(tMeshElmt),pointer    :: MeshElmts(:)     ! Gitter Klassen fur jede Verbindung
  integer,pointer            :: meshpattern(:,:) ! Definition von Elementen !
  integer,pointer            :: meshconnect(:,:) ! Definition von Verbindungen !
  integer,pointer            :: Dunknowns(:,:)
  integer,pointer            :: Fmovemt(:,:)
  real,pointer               :: Dimposed(:,:)
  real,pointer               :: Fimposed(:,:)
  ! lokale
  integer                    :: i,j,nj
  
! Remark: the program shall be completed
  do i=1, MeshInfo%nt
     do j=1,2
        nj=meshpattern(i,j+1)   ! index of Element
        select case(MeshT%elmts(i,2))
           case(1)                 ! + 1= >o--  bar articulated
               Dunknowns(i,3*(j-1)+1) = 1
               Dunknowns(i,3*(j-1)+2) = 1
               Dunknowns(i,3*(j-1)+3) = 0  ! angle free

               Fmovemt(i,3*(j-1)+1) = 0    ! force to determine
               Fmovemt(i,3*(j-1)+2) = 0
               Fmovemt(i,3*(j-1)+3) = 1    ! moment equal zero, thus imposed

           case(2)                 ! + 2= :>--  bar in a rail
               if (Dimposed(nj,1).le.(DLIMIT)) then ! artefact to show that no displacemt is imposed
                 Dunknowns(i,3*(j-1)+1) = 0
                 Dunknowns(i,3*(j-1)+2) = 1
               elseif (Dimposed(nj,2).le.(DLIMIT)) then
                 Dunknowns(i,3*(j-1)+1) = 1
                 Dunknowns(i,3*(j-1)+2) = 0
               endif
               Dunknowns(i,3*(j-1)+3) = 0  ! Moment free

               if (Fimposed(nj,1).le.(FLIMIT)) then
                 Fmovemt(i,3*(j-1)+1) = 0   ! artefact to show that no force is imposed
                 Fmovemt(i,3*(j-1)+2) = 1
               elseif (Fimposed(nj,2).le.(FLIMIT)) then
                 Fmovemt(i,3*(j-1)+1) = 1
                 Fmovemt(i,3*(j-1)+2) = 0    ! force to determine
               endif
               Fmovemt(i,3*(j-1)+3) = 1    ! moment equal zero, thus imposed

           case(3)                 ! + 3= //=   beam fixed
               Dunknowns(i,3*(j-1)+1) = 1
               Dunknowns(i,3*(j-1)+2) = 1
               Dunknowns(i,3*(j-1)+3) = 1

               Fmovemt(i,3*(j-1)+1) = 0    ! force to determine
               Fmovemt(i,3*(j-1)+2) = 0
               Fmovemt(i,3*(j-1)+3) = 0    ! moment to determine

           case(4)                 ! + 4= :/=   beam in rail
               if (Dimposed(nj,1).le.(DLIMIT)) then
                 Dunknowns(i,3*(j-1)+1) = 0
                 Dunknowns(i,3*(j-1)+2) = 1
               elseif (Dimposed(nj,2).le.(DLIMIT)) then
                 Dunknowns(i,3*(j-1)+1) = 1
                 Dunknowns(i,3*(j-1)+2) = 0
               endif
               Dunknowns(i,3*(j-1)+3) = 1  !

               if (Fimposed(nj,1).le.(FLIMIT)) then
                 Fmovemt(i,3*(j-1)+1) = 0   ! artefact to show that no force is imposed
                 Fmovemt(i,3*(j-1)+2) = 1
               elseif (Fimposed(nj,2).le.(FLIMIT)) then
                 Fmovemt(i,3*(j-1)+1) = 1
                 Fmovemt(i,3*(j-1)+2) = 0    ! force to determine
               endif
               Fmovemt(i,3*(j-1)+3) = 0    ! moment to determine 

!          case(5)                 ! + 5= |o|=  beam articulated

          case(6)                 ! + 6= ->o<- bar to bar
               Dunknowns(i,3*(j-1)+1) = 0  ! default, to determine, see after
               Dunknowns(i,3*(j-1)+2) = 0  ! default, to determine, see after
               if (Dimposed(nj,1).ge.(DLIMIT)) then
                 Dunknowns(i,3*(j-1)+1) = 1
               endif
               if (Dimposed(nj,2).ge.(DLIMIT)) then
                 Dunknowns(i,3*(j-1)+2) = 1
               endif
               Dunknowns(i,3*(j-1)+3) = 0  ! angle free

               Fmovemt(i,3*(j-1)+1) = 1    ! default, force provided, see after
               Fmovemt(i,3*(j-1)+2) = 1
               if (Fimposed(nj,1).le.(FLIMIT)) then
                 Fmovemt(i,3*(j-1)+1) = 0  
               endif
               if (Fimposed(nj,2).le.(FLIMIT)) then
                 Fmovemt(i,3*(j-1)+2) = 0  
               endif
               Fmovemt(i,3*(j-1)+3) = 1    ! moment equal zero, thus imposed

          case(7)                  ! + 7= =| |= beam to beam
               Dunknowns(i,3*(j-1)+1) = 1  ! default, imposed, see after
               Dunknowns(i,3*(j-1)+2) = 1
               if (Dimposed(nj,1).le.(DLIMIT)) then
                 Dunknowns(i,3*(j-1)+1) = 0
               endif
               if (Dimposed(nj,2).le.(DLIMIT)) then
                 Dunknowns(i,3*(j-1)+2) = 1
               endif
               Dunknowns(i,3*(j-1)+3) = 1  ! angle equal to adjacent beam element

               Fmovemt(i,3*(j-1)+1) = 0    ! default, force provided, see after
               Fmovemt(i,3*(j-1)+2) = 0
               if (Fimposed(nj,1).le.(FLIMIT)) then
                 Fmovemt(i,3*(j-1)+1) = 0  
               endif
               if (Fimposed(nj,2).le.(FLIMIT)) then
                 Fmovemt(i,3*(j-1)+2) = 0  
               endif
               Fmovemt(i,3*(j-1)+3) = 0  ! to determine with moment at the all end of the beams
        end select
    enddo
  enddo
  end subroutine InputSetUnknown


end module InputTruss_mod
