program mainTestConn
 use ISO_C_BINDING, only : C_CHAR, C_INT, C_DOUBLE, C_PTR
 use typesround
 use TestConnAppelC

 implicit none

  character(len=30)             :: Subname

  type(C_PTR)                   :: list

  type(C_PTR)                   :: edge_cell, edge_elmts, edge_iterator, edgecopy ! create_edge, get_next_edge ! remark C-functions

 ! type(edge),pointer            :: edge_elmt(:)

 ! type(tree),pointer            :: tree_elmt(:)

  integer(kind=C_INT)           :: leng, length_edge
  integer                       :: i, n1, n2, info
  real(8)                       :: r
 
  edge_elmts = create_edge_f(1,2,dble(0.5))
  info=removeheap_list_f(edge_elmts)  ! trick to have a NULL list
  if (info.eq.0) then
    do i=1,5
      read(1001,*) n1,n2,r
       edge_cell = create_edge_f(n1,n2,r)
     call add_end_edge_f( edge_cell, edge_elmts) 
    enddo

    do i=1,5
      call cp_edge_nodes_f(edge_elmts, edgecopy)
      info=removeheap_list_f(edge_elmts)
      r=get_weight_f( edgecopy )   ! r=edgecopy%weight
 !     n1=edgecopy%el_1
 !     n2=edgecopy%el_2
      write(*,1002) i, r 
    enddo
  endif
1001 format(i10, i10, e10.4)
1002 format(i10, e10.4)
end program mainTestConn
