#include "solidmecanicsC.h"


int main() {
  int NEXT_ID, len_conn, len_camino, count;
  Index couplenode[2];
  Index testid;
  Number weight;
  ptr_EDGE conn_elmt, elmt, iter, copie;
  ptr_TREE camino, feather;

  conn_elmt = create_edge(1, 2, 0.707);
  elmt = create_edge(2, 3, 1.0);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(3, 4, 0.707);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(2, 5, 1.0);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(4, 7, 1.0);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(5, 6, 0.707);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(6, 7, 1.0);
  add_end_edge(elmt, conn_elmt);
  elmt = create_edge(7, 8, 0.707);
  add_end_edge(elmt, conn_elmt);

  elmt = create_edge(7, 8, 0.707);
  add_end_edge(elmt, conn_elmt);
  
  iter = conn_elmt;
  cp_edge_nodes(iter, copie);  
  couplenode[0]=copie->node1;
  couplenode[1]=copie->node2;
  printf("couple of nodes in position %d: (%d, %d)\n", 1, couplenode[0], couplenode[1]);
  len_conn = length_edge(conn_elmt);
  if (remove_end_edge(conn_elmt) != 0) {
     printf("error in removing last elmt edge in position: %d\n", len_conn);
  }
  else {
     len_conn = length_edge(conn_elmt);
     printf("success in removing last elmt edge, now edge has length: %d\n", len_conn);
  }

  camino =  create_tree(couplenode[0], couplenode[1], 1);
  count = 1;
  while (count < len_conn) {
     cp_edge_nodes(iter, copie);  
     couplenode[0]=copie->node1;
     couplenode[1]=copie->node2;
     feather = create_tree(couplenode[0], couplenode[1], 1);
     add_end_tree(feather, camino);
     iter = get_next_edge(iter);
     count=count+1;
  }
  len_camino = length_tree(camino);
  if (len_camino!=len_conn) {
     printf("error in building tree list(%d) from edge list(%d)", len_camino, len_conn);
  }
  else {
     printf("Sucess, program temminates\n");
  }

  return 0;
}


