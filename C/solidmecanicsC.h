#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/*
* Do not set any modifications to these structures
* as they are previewed for an Algorithm with list
* owned by different ClusterId.
*
*/
typedef int Index;
typedef double Number;

typedef struct node {
	
   Index ind;
   Number  vects[3];
} * ptr_NODE;
/* Definition of linked list <edge>
 */
typedef struct edge {

   Index node1, node2; 
   Number weight;
   struct edge *next;
   
} * ptr_EDGE;
// to implement are:
ptr_EDGE create_edge(Index nd1, Index nd2, Number we);
void add_end_edge(ptr_EDGE element, ptr_EDGE list);
int remove_end_edge(ptr_EDGE element);
int removeheap_list(ptr_EDGE cell);
void cp_edge_nodes(ptr_EDGE element, ptr_EDGE copie);
Number get_weight(ptr_EDGE element);
ptr_EDGE get_next_edge(ptr_EDGE element);
int length_edge(ptr_EDGE list);
ptr_EDGE reverse_list(ptr_EDGE list);
void sort_elmtedge_byweight(ptr_EDGE list);

/* Definition of linked list <tree>
*/ 
typedef struct tree {

   Index node1, node2;
   Index clusterId; // the final tree will be constructed by assembling tree of elements with different clusterIds
   struct tree *next;
} * ptr_TREE;

// to implement are:
ptr_TREE create_tree(Index nd1, Index nd2, Index id);
void add_end_tree(ptr_TREE element, ptr_TREE list);
int remove_end_tree(ptr_TREE element);
void cp_tree_nodes(ptr_TREE element, ptr_TREE copie);
Index get_clusterId(ptr_TREE element);
ptr_TREE get_next_tree(ptr_TREE element);
int length_tree(ptr_TREE list);

/*
 * Implementation of help sub routines for Kruskal
 * */

int sort_edge_weight(ptr_EDGE list);


