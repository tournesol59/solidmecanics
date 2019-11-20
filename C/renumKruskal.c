/******************************************************
 Implementation of several functions that will be called
 from Fortran, which realizes Kruskal algorithm: 
 the finding of a best (shortest edges sum of weight) 
 tree of connexion between Vertices of a CONNEX graph.
 
 This seems not to be unaffordable as set of C code (but
 a little much less readable in Fortran.) to get at the
 end a complete renumeration of the Finite Element given 
 by the mesh, which will be advantageous for sparse matrix
 inversion.

 Input:
   ++ from Fortran: two allocated array of elements (
       which will be the original graph) and nodes
       (for weight computation)

 Intern Data:
   ++ a linked list of the couple of connected elements
     with inverse of common frontier as weight
     (el1, el2, weight) and clusterId for tree of appartenance
     for Kruskal algorithm implementation.

   ++ a linked list of the final tree computed by the algorithm.
      
  Output:
  ++ from/for Fortran: an allocated array representing the final tree
     with usable methods:
     i: new Index of element
     tree[i] structure with Index of Elements in allocated arrays in Fortran
                       with 3 couples of (new Index, Fortran Index) of Nodes   
                       with Index i+1 or -1 if branch of the tree is ended
		       with Index nextb1, nextb2 (max two) of possible two other                         connections
 IPST CNAM		       
*******************************************************/

#include "solidmecanicsC.h"
#include <malloc.h>
/*
 * Definition of linked list edge
 */
ptr_EDGE create_edge(Index nd1, Index nd2, Number we) {

  struct edge* element=(ptr_EDGE) malloc(sizeof(struct edge));

  (*element).node1=nd1;
  (*element).node2=nd2;
  (*element).weight=we;
  (*element).next = NULL;

  return element;
}


int removeheap_list(ptr_EDGE cell) {  // cannot do else a
  int i, leng, res;
  ptr_EDGE iterator;

  leng = length_edge(cell);
  if (leng != 0) {
    for (i=0; i<leng; i++) {
       iterator=cell;
       res=remove_end_edge(iterator);
       if (res!=0) {
         printf("Error removeheap_list\n");
         exit(1);
       }
    }
//    free( cell );  // already done
  } //endif
  return res;
}

void add_end_edge(ptr_EDGE element, ptr_EDGE list) {
  
  ptr_EDGE iterator;

  iterator=list;
  if (iterator != NULL) {
    while ((*iterator).next != NULL) {
       iterator = (*iterator).next;
    }

    (*iterator).next = element;
  }
  else {
     list = element;
  }
}

int remove_end_edge(ptr_EDGE element) {
  
  ptr_EDGE iterator, el;
  int res;

  el=element;
  iterator=element;
  if (iterator != NULL) {
    while ((*iterator).next != NULL) {
       el = iterator;
       iterator = (*iterator).next;
    }

    free(iterator);
    (*el).next = NULL;
    res=0;
  }
  else {
    res=-1;
  }
  return res;
}

void cp_edge_nodes(ptr_EDGE element, ptr_EDGE copie) {

  if (element !=NULL) {
    (*copie).node1=(*element).node1;
    (*copie).node2=(*element).node2;
    //(*copie).weight=(*element).weight;
  }

}

Number get_weight(ptr_EDGE element) {
  
  Number res;

  res=(*element).weight;
  return res;
}

ptr_EDGE get_next_edge(ptr_EDGE element) {

  ptr_EDGE res;
	if (element != NULL) {
		res = (*element).next;
	}
	else {
		printf("Error get_next_edge argument 1 is NULL\n");
		res = NULL;
	}

  return res;
}

int length_edge(ptr_EDGE list) {
	  
  ptr_EDGE iterator;
  int count=0;

  iterator=list;
  if (iterator != NULL) {
    count=count+1;
    while ((*iterator).next != NULL) {
       iterator = (*iterator).next;
       count=count+1;
    }
  }
  return count;
}

ptr_EDGE reverse_list(ptr_EDGE list) {
  ptr_EDGE reverse=NULL;
  ptr_EDGE iterator, el=NULL;
  Index count, k;

  count = length_edge(list);
  iterator = list;

  for (k=0; k<count; k++) {
    el=create_edge(1, 2, 1.0);
    cp_edge_nodes(iterator, el);
    (*el).next=reverse;
    iterator=(*iterator).next;
  }
  free(iterator);
  free(el);
  return reverse; 
}

/* sorting */

void sort_elmtedge_byweight(ptr_EDGE list) {
// iterative sort O(N*(N-1))
  Index k,j,len;
  Number we, wemin;
  ptr_EDGE iterator, edgemin, copy;

  len = length_edge(list);
  edgemin = list;
  iterator = list;
  copy = create_edge(1,2,1.0);
//initial:
  for (k=1; k<(len+1); k++) {
    printf("i:%d, w:%f   ",k,get_weight(iterator) );
    iterator=get_next_edge(iterator);  
  }
  iterator = list;

  for (k=1; k<(len+1); k++) { // only len-1
     wemin=get_weight(edgemin);
     iterator = get_next_edge(edgemin);
     for (j=k+1; j<(length_edge(edgemin)+k); j++) { // up to len
        we=get_weight(iterator);
        if (we<wemin) {
           cp_edge_nodes(iterator, copy);  // affect only nodes records, nor weight nor next cell
           (*copy).weight = (*iterator).weight;//new wemin
           cp_edge_nodes(edgemin, iterator);
           (*iterator).weight = (*edgemin).weight;
           cp_edge_nodes(copy, edgemin);
           (*edgemin).weight = (*copy).weight; // new wemin
           wemin=we;
        }
        iterator=get_next_edge(iterator);
     }
     edgemin=get_next_edge(edgemin);
    // if (edgemin != NULL) { // not useful
     iterator=get_next_edge(edgemin);
    // }
  }
  // verify:
  iterator = list;
  for (k=1; k<(len+1); k++) {
    printf("i:%d, w:%f   ",k,get_weight(iterator) );
    iterator=get_next_edge(iterator);  
  }


  return;
}

/*
 * Definition of linked list tree
 */
ptr_TREE create_tree(Index nd1, Index nd2, Index id) {

  struct tree element;
  ptr_TREE res;

  element.node1=nd1;
  element.node2=nd2;
  element.clusterId=id;
  element.next = NULL;
  res = &element;
  return res;
}

void add_end_tree(ptr_TREE element, ptr_TREE list) {
  
  ptr_TREE iterator;

  iterator=list;
  if (iterator != NULL) {
    while ((*iterator).next != NULL) {
       iterator = (*iterator).next;
    }

    (*iterator).next = element;
  }
  else {
     list = element;
  }
}

int remove_end_tree(ptr_TREE element) {
  
  ptr_TREE iterator, el;
  Index res;

  el=element;
  iterator=element;
  if (iterator != NULL) {
    while ((*iterator).next != NULL) {
       el = iterator;
       iterator = (*iterator).next;
    }

    free(iterator);
    (*el).next = NULL;
    res = 0;
  }
  else {
    res = -1;
  }
  return res;
}

void cp_tree_nodes(ptr_TREE element, ptr_TREE copie) {

  if (element !=NULL) {
    (*copie).node1=(*element).node1;
    (*copie).node2=(*element).node2;
  }

}


Index get_clusterId(ptr_TREE element) {
  
  Index id;

  id=(*element).clusterId;
  return id;
}

ptr_TREE get_next_tree(ptr_TREE element) {

  ptr_TREE res;
	if (element != NULL) {
		res = (*element).next;
	}
	else {
		printf("Error get_next_tree argument 1 is NULL\n");
		res = NULL;
	}

  return res; 
}

int length_tree(ptr_TREE list) {
	  
  ptr_TREE iterator;
  int count=0;

  iterator=list;
  if (iterator != NULL) {
    count=count+1;
    while ((*iterator).next != NULL) {
       iterator = (*iterator).next;
       count=count+1;
    }
  }
  return count;
}

