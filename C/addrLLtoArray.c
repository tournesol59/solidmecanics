#include "solidmecanicsC.h"

int prnt_edge_list(char Subname[30], ptr_EDGE list) {

  ptr_EDGE                 cell;

  int                      i, leng, n1, n2;
  double                   weight;

  printf(" ---------%s\n", Subname);
  leng = length_edge(list);
  cell = list;
  if (leng>1) {
    for (i=1; i<=leng; i++) {
      printf("  node = (%d, %d) \n", cell->node1, cell->node2);
      weight=get_weight(cell);
      printf("   %f \n",weight);
      
      cell=get_next_edge(cell);
    }
  }
}


int main() {
  int info, NEXT_ID, len_conn, len_camino, count;
  Index couplenode[2];
  Index testid;
  Number weight;
  ptr_EDGE conn_elmt, elmt, iter, copie;
  //ptr_TREE camino, feather;

  conn_elmt = create_edge(1, 2, 0.707);
  elmt = create_edge(2, 3, 1.0);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(1, 12, 1.0);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(12, 11, 0.707);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(11, 10, 1.0);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(10, 9, 0.707);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(9, 8, 1.0);
  add_end_edge(elmt, conn_elmt);

  elmt = create_edge(3, 4, 0.707);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(4, 5, 1.0);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(4, 6, 0.707);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(6, 7, 1.0);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(7, 8, 0.707);
  add_end_edge(elmt, conn_elmt);

  info = prnt_edge_list("main                          ", conn_elmt);

  sort_elmtedge_byweight(conn_elmt);

  info = prnt_edge_list("main                          ", conn_elmt);
  return 0;
}





