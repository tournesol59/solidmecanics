program main

  use types
  use input_mod
  use allocatefields_mod
  use createmesh_mod
  use interpol2d_mod
!  use calcgradientsurf_mod

  implicit none
  ! Testen CalcGradientSurf, 

  ! Parameter des Gitter,netz: Ebene auf eine Fläche
 ! integer      :: NRAUMX=20
 ! integer      :: NRAUMY=10
 ! real         :: STARTX=0.0
 ! real         :: ENDX=30.0
 ! real         :: STARTY=0.0
 ! real         :: ENDY=15.0
  ! Referenz Flaeche ist ein Rechteck-quadratnetz mit NRAUMX*NRAUMY Elementen

  ! Parameter der Interpolation
  integer      :: ORDER=2

  !--------------------------------------------------------------------------!
  ! Variablendeklarationen                                                   !
  !--------------------------------------------------------------------------!
  type(tMesh)                :: Gitter     ! Gitterwerte                     !
  type(tRandbedingungen)     :: RB         ! Randbedingungen                 !
  type(tConstants)           :: Const      ! Konstanten                      !
  type(tFileIO)              :: FileIO     ! Ausgabesteuerung                !
  type(tExakt)               :: Exakt      ! Exakte Loesung                  !
  real,pointer               :: Uvar(:,:)  ! Feld mit der numerischen Loesung!!=Tr(spanng) !
  real,pointer               :: rhs(:,:)   ! Feld für die Rechte Seite       !
  type(tNumeric)             :: VarNum     ! Numerische Felder mit Loesungen   !
  !--------------------------------------------------------------------------!   

  real,pointer        :: chicoeff(:,:)
  real,pointer        :: der_a1xx_x(:),der_a1xx_y(:),der_a1xy_x(:),der_a1xy_y(:) &
                        ,der_a2yx_x(:),der_a2yx_y(:),der_a2yy_x(:),der_a2yy_y(:)

  integer             :: allocStat
  ! ------------------------------------------< Eingabegroessen einlesen >---!

  call input(Gitter,Exakt,Const,FileIO)       

  ! ------------------------------------------< Speicherplatz allokieren >---!

  call allocatefields(Gitter,RB,Uvar,rhs,Exakt,chicoeff,Const,VarNUm)  

  allocate(der_a1xx_x(100))
  allocate(der_a1xx_y(100))

  allocate(der_a1xy_x(100))
  allocate(der_a1xy_y(100))

  allocate(der_a2yx_x(100))
  allocate(der_a2yx_y(100))

  allocate(der_a2yy_x(100))
  allocate(der_a2yy_y(100))

  ! --------------------------------------------------< Gitter festlegen >---!
   if (Const%auto == 1) then
      call createMesh(Gitter)    
   else 
      call createMesh2(Gitter)
   endif
   call trychristoffei(Gitter, chicoeff)

  write(*,*) '================    Interpolation starts now    ================ '
  write(*,*)

  call basisdifferenz(Gitter,der_a1xx_x,der_a1xx_y,der_a1xy_x,der_a1xy_y &
                                 ,der_a2yy_x, der_a2yy_y, der_a2yx_x, der_a2yx_y)

  write(*,*) '     differenz der lokalen Basis kompletiert    '
  call trychristoffei(Gitter, chicoeff)


  call christoffei(Gitter, chicoeff, der_a1xx_x,der_a1xx_y,der_a1xy_x,der_a1xy_y &
                                 ,der_a2yy_x, der_a2yy_y, der_a2yx_x, der_a2yx_y)

  write(*,*) '     Christoffei Koeffizienten berechnet   '
  call trychristoffei(Gitter, chicoeff)


  ! ------------------------------------------< Speicherplatz allokieren >---!

  !call deallocateFields(Gitter,RB,Uvar,rhs,Exakt,chicoeff,Const) 
    deallocate(Gitter%x, Gitter%y, RB%randl, RB%randr, RB%rando, RB%randu, &
         Exakt%loesung, Uvar, rhs, chicoeff, STAT = allocStat )
    
    if (allocStat.NE.0) then
       print *, 'ERROR AllocateFields: Could not deallocate correctly!'
       STOP
    end if 
  
  deallocate(der_a1xx_x)
  deallocate(der_a1xx_y)

  deallocate(der_a1xy_x)
  deallocate(der_a1xy_y)

  deallocate(der_a2yx_x)
  deallocate(der_a2yx_y)

  deallocate(der_a2yy_x)
  deallocate(der_a2yy_y)

    write(*,*) 
    write(*,*) '============= Program terminated correctly ================ '

end program main
