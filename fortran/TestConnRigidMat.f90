! Module to predict the number of non-zero Elements of the Rigidity Elastic Matrix
! with the help of Linked List in C
! preceeded by a first AppelC modules, declaration of passed C functions

! RoadMap: complete C procedures based on functions: sorting, modifying k-th cell, inserting cell
!          start with creating list of connected elements in C
!          follow (1) with test of Kruskal algo in Fortran

!****************************************************

#ifndef _GFORTRAN03
#define _GFORTRAN03_
#endif

module TestConnAppelC
  use ISO_C_BINDING, only : C_CHAR, C_INT, C_DOUBLE, C_PTR
  implicit none

    type, bind(C) ::  edge
      integer(kind=C_INT)                         :: el_1, el_2
      real(kind=C_DOUBLE)                         :: weight
      type(C_PTR)                                 :: cel
    end type edge 

   interface
    type(C_PTR) function create_edge_f(n1, n2, r) bind (C,name='create_edge')
      import C_INT, C_DOUBLE, C_PTR
      real(kind=C_DOUBLE)                 :: r
      integer(kind=C_INT)                 :: n1,n2
    end function create_edge_f


    subroutine add_end_edge_f( cell, endcell) bind (C,name='add_end_edge')
      import C_INT, C_DOUBLE, C_PTR
      type(C_PTR)                         :: cell, endcell
    end subroutine add_end_edge_f


    integer function remove_end_edge_f(cell) bind (C,name='remove_end_edge')
      import C_INT, C_DOUBLE, C_PTR
      type(C_PTR)                         :: cell 
    end function remove_end_edge_f


    integer function removeheap_list_f(cell) bind (C,name='removeheap_list')
      import C_INT, C_DOUBLE, C_PTR
      type(C_PTR)                         :: cell 
    end function removeheap_list_f


    subroutine cp_edge_nodes_f(cell, copy) bind (C,name='remove_end_edge')
      import C_INT, C_DOUBLE, C_PTR
      type(C_PTR)                          :: cell, copy
    end subroutine cp_edge_nodes_f


    real(kind=C_DOUBLE) function get_weight_f(cell) bind (C,name='get_weight')
      import C_INT, C_DOUBLE, C_PTR
      type(C_PTR)                          :: cell
    end function get_weight_f


    type(C_PTR) function get_next_edge_f(cell) bind (C,name='get_next_edge')
      import C_INT, C_DOUBLE, C_PTR
      type(C_PTR)                          :: cell
    end function get_next_edge_f


!   type(C_INT) function length_edge_f(cell) bind (C,name='length_edge')
    integer(kind=C_INT) function length_edge_f(cell) bind (C,name='length_edge')
      import C_INT, C_DOUBLE, C_PTR
      type(C_PTR)                          :: cell
    end function length_edge_f


    type(C_PTR) function reverse_list_f(cell) bind (C,name='reverse_list')
      import C_INT, C_DOUBLE, C_PTR
      type(C_PTR)                          :: cell
    end function reverse_list_f


!    type(C_INT) function find_edge_f(cell,i,j) bind (C,name='find_edge')
    integer(kind=C_INT) function find_edge_f(cell,i,j) bind (C,name='find_edge')
      import C_INT, C_DOUBLE, C_PTR
      type(C_PTR)                          :: cell
      integer(kind=C_INT)                  :: i,j
    end function find_edge_f

   end interface

!   public ::  edge, create_edge_f, add_end_edge_f, remove_end_edge_f, removeheap_list_f, &
!             cp_edge_nodes_f, get_weight_f, get_next_edge_f, length_edge_f, reverse_list_f

end module TestConnAppelC

!****************************************************
module TestConnRigidMat
  use ISO_C_BINDING, only : C_CHAR, C_INT, C_DOUBLE, C_PTR
  use typesround
 ! use TestConnAppelC
  implicit none

  interface create_ll_connectivity
      module procedure create_ll_connectivity
  end interface

  public:: create_ll_connectivity, frontnodesConnLL, geti_ownby_el, geti_neighof_el, &
           copy_node
!          create_edge_f, add_end_edge_f, remove_end_edge_f, &
!           cp_edge_nodes_f, get_weight_f, get_next_edge_f, length_edge_f, &
!          reverse_list_f, find_edge_f  !,..


  contains
!****************************************
! access function 
integer function geti_ownby_el(MeshR, i,j)
  use typesround
  implicit none
  type(tRoundMesh), intent(in)                     :: MeshR
  integer,intent(in)                               :: i,j
  if ((j.le.4).and.(j.ge.1)) then
     geti_ownby_el = MeshR%elems(i)%numeros(j)
  else
     geti_ownby_el = -1
  endif
end function geti_ownby_el

!****************************************
! access function 
integer function geti_neighof_el(MeshR, i,j)
  use typesround
  implicit none
  type(tRoundMesh), intent(in)                     :: MeshR
  integer,intent(in)                               :: i,j

  if ((j.le.4).and.(j.ge.1)) then
     geti_neighof_el = MeshR%neighbours(i)%numeros(j)
  else
     geti_neighof_el = -1
  endif
end function geti_neighof_el

!****************************************
! access subroutine copy_node
subroutine copy_node(MeshR, i,k,node_1)
  use typesround
  implicit none
  type(tRoundMesh), intent(in)                     :: MeshR
  type(tNode), intent(inout)                       :: node_1
  integer, intent(in)                              :: i,k
  integer                                          :: l

  do l=1,3
    select case (i)
      case (1)
          node_1%vects(l)=MeshR%nodes(MeshR%elems(k)%numeros(2))%vects(l)
      case (2) 
          node_1%vects(l)=MeshR%nodes(MeshR%elems(k)%numeros(3))%vects(l)
      case (3) 
          node_1%vects(l)=MeshR%nodes(MeshR%elems(k)%numeros(4))%vects(l)
   end select
  enddo
end subroutine copy_node

!***************************************
! print list for testing
subroutine prntTestConnLL(Subname, list)
  use ISO_C_BINDING, only : C_CHAR, C_INT, C_DOUBLE, C_PTR
  use TestConnAppelC
  use typesround
  implicit none

  character(len=30),intent(in)  :: Subname
#ifndef _FORTRAN03_
  integer(kind=8),intent(in)    :: list
#else
  type(C_PTR),intent(in)        :: list
#endif
#ifndef _FORTRAN03_
  integer(kind=8)               :: cell, iterator, create_edge, get_next_edge ! remark C-functions address in f90
#else
  type(C_PTR)                   :: cell, iterator
#endif
  type(edge),pointer            :: elmt(:)
  integer(kind=C_INT)           :: leng, length_edge
  integer                       :: i, n1, n2
  real(kind=C_DOUBLE)           :: weight, weightex, get_weight

  print *, " ---------" , Subname
  leng = length_edge(list)
  if (leng.ge.1) then
    iterator =  list
    do i=1,leng
      weightex=dble(1.0)
      cell=create_edge(1,2,weightex)
      call cp_edge_nodes(iterator, cell)  ! yes with arg integer(kind=8)
      call c_f_pointer(cell,elmt)

      write(*,201) "  node = ", elmt(1)%el_1, elmt(1)%el_2
      weight=get_weight(iterator)
      write(*,301) "         ", (weight)
      iterator=get_next_edge(iterator)
    enddo
  endif

101 format(a,a)
201 format(a,2i10)
301 format(a,f10.6)

end subroutine prntTestConnLL

!***************************************
! real calc Weight
real function calcWeightConnLL(node_1, node_2, node_3, node_4)
  use typesround
  implicit none

  type(tNode)                 :: node_1, node_2, node_3, node_4 ! nodes 1 and 2 are common to 2 elmts (frontier)
  intent(in)                  :: node_1, node_2, node_3, node_4
  integer                     :: i
  real                        :: weight

  weight = 0.0
  do i=1,3
    weight = weight + (node_2%vects(i) - node_1%vects(i))**2
  enddo
  calcWeightConnLL = weight
end function calcWeightConnLL

!***************************************
! int identify frontier nodes
subroutine frontnodesConnLL(p,q,r,s,t,u, couplea, coupleb)
  implicit none

  integer                    :: p,q,r,s,t,u
  integer,dimension(2)       :: couplea,coupleb
  intent(in)                 :: p,q,r,s,t,u
  intent(inout)              :: couplea,coupleb

  if (((p.eq.s).and.(q.eq.t)).or.((p.eq.t).and.(q.eq.s))) then
      couplea(1)=p
      couplea(2)=q
      coupleb(1)=s
      coupleb(2)=t
  elseif (((p.eq.s).and.(q.eq.u)).or.((p.eq.u).and.(q.eq.s))) then
      couplea(1)=p
      couplea(2)=q
      coupleb(1)=s
      coupleb(2)=u
  elseif (((p.eq.t).and.(q.eq.u)).or.((p.eq.u).and.(q.eq.t))) then
      couplea(1)=p
      couplea(2)=q
      coupleb(1)=t
      coupleb(2)=u
  elseif (((r.eq.s).and.(p.eq.t)).or.((r.eq.t).and.(p.eq.s))) then
      couplea(1)=r
      couplea(2)=p
      coupleb(1)=t
      coupleb(2)=s
  elseif (((r.eq.s).and.(p.eq.u)).or.((r.eq.u).and.(p.eq.s))) then
      couplea(1)=r
      couplea(2)=p
      coupleb(1)=s
      coupleb(2)=u
  elseif (((r.eq.t).and.(p.eq.u)).or.((r.eq.u).and.(p.eq.t))) then
      couplea(1)=r
      couplea(2)=p
      coupleb(1)=t
      coupleb(2)=u
  elseif (((q.eq.s).and.(r.eq.t)).or.((q.eq.t).and.(r.eq.s))) then
      couplea(1)=q
      couplea(2)=r
      coupleb(1)=s
      coupleb(2)=t
  elseif (((q.eq.s).and.(r.eq.u)).or.((q.eq.u).and.(r.eq.s))) then
      couplea(1)=q
      couplea(2)=r
      coupleb(1)=s
      coupleb(2)=u
  elseif (((q.eq.t).and.(r.eq.u)).or.((q.eq.u).and.(r.eq.t))) then
      couplea(1)=q
      couplea(2)=r
      coupleb(1)=t
      coupleb(2)=u
  else
      couplea(1)=-1
      couplea(2)=-1
      coupleb(1)=-1
      coupleb(2)=-1
  endif


end subroutine frontnodesConnLL

!****************************************
! subroutine create_ll_connectivity
subroutine create_ll_connectivity(NTestMeshR, MeshR, infoconnect)
  use ISO_C_BINDING, only : C_CHAR, C_INT, C_DOUBLE, C_PTR 
  use typesround  ! in particular tMeshRound
!  use TestConnAppelC
  implicit none

! uebergegebene args
  type(tRoundMeshInfo)             :: NTestMeshR
  type(tRoundMesh)                 :: MeshR
#ifndef _FORTRAN03_
  integer(kind=8)                  :: infoconnect
#else
  type(C_PTR)                      :: infoconnect
#endif
  intent(in)                       :: NTestMeshR, MeshR
  intent(out)                      :: infoconnect
! lokale vars
#ifndef _FORTRAN03_
  integer(kind=8)                  :: cell , create_edge  ! remark C-functions address in f90
#else
  type(C_PTR)                      :: cell
#endif
  integer                          :: ne,countel, i,j,k,l, nb, p,q,r,s,t,u, allocinfo
  integer,dimension(:),allocatable :: visited_el
  real*8                           :: weight1, weightex
  type(tNode)                      :: node_1, node_2, node_3, node_4, node_5, node_6


  ne=NTestMeshR%nElem
!  allocate(node_1%vects(3),node_2%vects(3),node_3%vects(3),node_4%vects(3), &
!           node_5%vects(3), node_6%vects(3), visited_el(ne), STAT=allocinfo)
  allocate(visited_el(ne), STAT=allocinfo)
  if (allocinfo.ne.0) then
     STOP
     print *," --------- create_ll_connectivity: Error allocation copy node"
  endif
  do i=1,ne
    visited_el(i)=0
  enddo
  weightex=dble(1.0)
  infoconnect=create_edge(1,2,weightex)

  k=1 ! first element
  countel = 1

  do while (countel<ne)  ! parcourir les elements du tableau MeshR sachant que
! la variable k sera assigne "suivant les voisins selon le contenu de MeshR" a chaque valeur de l'intervalle 
! [1, ne] il est donc necessaire d'utiliser un tableau d'entier "visités"

    
   ! recuperer les numeros de chaque sommets
   ! ..copy_node(MeshR, k, 1, node_1)
   ! ..copy_node(MeshR, k, 2, node_2)
   ! ..copy_node(MeshR, k, 3, node_3)
   ! parcourir une 2eme boucle les 3 voisins et tester quel cote est communs pour chaque voisins a l element
    j=1
    nb=geti_neighof_el(MeshR, k,j)
    do while ( (j.le.3).and.(nb.ge.0) )
 ! if numeros(j).eq.(-1) it means no frontier with another elmt =+>  ajouter par appel a une C fonction l'edge forme par {elmt_k, elms_i}
      ! recuperer les numeros de chaque sommets
      ! ..copy_node(MeshR, (geti_neighof_el(MeshR, k,j), 1, node_4)
      ! ..copy_node(MeshR, (geti_neighof_el(MeshR, k,j), 2, node_5)
      ! ..copy_node(MeshR, (geti_neighof_el(MeshR, k,j), 3, node_6)

     if ( visited_el(nb).eq.0 ) then ! if neighbour not already visited
        weightex=dble(1.0)
        cell=create_edge(k, nb,weightex)
        call add_end_edge(cell, infoconnect)
        k=nb  ! next iteration

     elseif (((j.eq.3).and.((visited_el(nb)).eq.0)).or.((j.le.2).and.(nb.le.0))) then ! if not any more neighbour to visit
      ! explore all remaining index of elements
        l=1
        do while ( ((visited_el(l)).eq.0).and.(l.le.ne) )
           l=l+1
        enddo
        if (l.le.ne) then
           k=l
        else
          print *," --------- create_ll_connectivity: narrow pass: all neighbours visited"
        endif
      ! 2nd case: if a neighbour has been visited but not the connection (k,neighbour(j)): add it
     elseif ( (visited_el(nb)).eq.1 ) then !TODO  .and.(edge_not_present_f(infoconnect,j)
        weightex=dble(1.0)
        cell=create_edge(k, nb, weightex)
        call add_end_edge_f(cell, infoconnect)
     endif
      j=j+1
    enddo ! j

    ! out of if visited neighbour test: 
    ! mark visitation
    visited_el(k)=1
    countel=countel+1
     ! important choice of k for the next iteration has been done above
  enddo

! comment following preview a completion of this procedure by adding not consequent edges:
!   k=1
!   do while (k.le.n)
!    j=1
!    do while ( (j.le.3).and.((geti_neighof_el(MeshR, k,j)).ge.0) )
! imp not already implemented:
!        if (.not.(find_edge_f(infoconnect,k,j)) then 
!           cell=create_edge_f(k,(geti_neighof_el(MeshR, k,j)), 1.0)
!           call add_end_edge_f(cell, infoconnect)     
!           j=j+1
!        endif
!     enddo
!   enddo

! deallocate(node_1%vects, node_2%vects, node_3%vects, node_4%vects, &
!          node_5, node_6, visited_el, STAT=allocinfo)
  deallocate(visited_el, STAT=allocinfo)
  if (allocinfo.ne.0) then
     STOP
     print *," --------- create_ll_connectivity: Error deallocation copy node\n"
  endif
end subroutine  create_ll_connectivity

end module TestConnRigidMat
