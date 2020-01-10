module InputTruss_mod

 use TypesBalken
 implicit none

  private
  interface Input_file
    module procedure Input_file
  end interface Input_file

  public:: Input_file

  contains

  subroutine Input_file(MeshInfo,MeshT,meshpattern,meshconnect,Dunknowns,Dimposed,Fmovemt,Fimposed)

  use typesbalken

  implicit none

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  ! Liste der übergebenen Argumente                                          !
  !                                                                          !
  type(tMeshInfo)           :: MeshInfo     ! Gitterwerte                    !
  type(tMeshCoord)           :: MeshT     ! Gitterwerte                      !
!  type(tMeshElmt),pointer    :: ElmtS(:)     ! Gitter Klassen fur jede Verbindung
  integer,pointer            :: meshpattern(:,:) ! Definition von Elementen !
  integer,pointer            :: meshconnect(:,:) ! Definition von Verbindungen !
  integer,pointer            :: Dunknowns(:,:)
  real,pointer               :: Dimposed(:,:)
  integer,pointer            :: Fmovemt(:,:)
  real,pointer               :: Fimposed(:,:)
!  type(tExakt)               :: Exakt      ! Exakte Loesung                 !
!  type(tFileIO)              :: FileIO     ! Ausgabesteuerung               !
  !--------------------------------------------------------------------------!  
  character(LEN=13)          :: Mshfilename  ! Name of the input file (13 characters) in current dir !
  integer :: i,countmax !, nbound, cbound, bnddim, bndind ! counter integers
  character(LEN=6) :: bndname

  OPEN(UNIT=25, FILE='Untitled1.msh', ACTION='READ')
  read(25,307) MeshInfo%nn, MeshInfo%ne

  do i=1,MeshInfo%nn
    read(25,207) MeshT%x(i), MeshT%y(i), MeshT%z(i)   ! z Koordinate muss geschrieben werden, wird aber nicht benutzt
  enddo

  do i=1,MeshInfo%ne
    read(25,309) meshpattern(i,1), meshpattern(i,2), meshpattern(i,3), meshpattern(i,4)
  enddo

  do i=1,MeshInfo%ne
    read(25,308) MeshT%elmts(i,1), MeshT%elmts(i,2), MeshT%elmts(i,3)  ! 
  enddo

  do i=1,MeshInfo%ne
    read(25,309) Dunknowns(i,1), Dunknowns(i,2), Dunknowns(i,3), Dunknowns(i,4)
  enddo

  do i=1,MeshInfo%ne
    read(25,207) Dimposed(i,1), Dimposed(i,2), Dimposed(i,3)
  enddo

  do i=1,MeshInfo%ne
    read(25,309) Fmovemt(i,1), Fmovemt(i,2), Fmovemt(i,3), Fmovemt(i,4)
  enddo

  do i=1,MeshInfo%ne
    read(25,207) Fimposed(i,1), Fimposed(i,2), Fimposed(i,3)
  enddo

  CLOSE(UNIT=25)

 105  format (a20)
 106  format (a,a)
 205  format (f18.16)
 206  format (a,f18.16)
 207  format (3f18.16)
 305  format (i10)
 306  format (a,i10)
 307  format (2i10)
 308  format (3i10)
 309  format (4i10)
 405  format (i10, i10, a6)
 406  format (i10, f18.16, f18.16, f18.16)
 507  format (i10, i10, i10, i10, i10, i10, i10)
 508  format (i10, i10, i10, i10, i10, i10, i10, i10)
 509 format (i10, i10, i10, i10, i10, i10, i10, i10, i10)

  end subroutine Input_file

end module InputTruss_mod
