/******** layout_qdp.c *********/
/* MIMD version 6 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */
/* This version uses QDP for layout */
/*
   setup_layout() does any initial setup.  When it is called the
     lattice dimensions nx,ny,nz and nt have been set.
     This routine sets the global variables "sites_on_node",
     "even_sites_on_node" and "odd_sites_on_node".
   node_number(x,y,z,t) returns the node number on which a site lives.
   node_index(x,y,z,t) returns the index of the site on the node - ie the
     site is lattice[node_index(x,y,z,t)].
   These routines will change as we change our minds about how to distribute
     sites among the nodes.  Hopefully the setup routines will work for any
     consistent choices. (ie node_index should return a different value for
     each site on the node.)
*/
#include "generic_includes.h"
#include <qdp.h>

void
setup_layout(void)
{
  int c[4];

  c[0] = nx;
  c[1] = ny;
  c[2] = nz;
  c[3] = nt;
  QDP_set_latsize(4, c);
  QDP_create_layout();
  sites_on_node = QDP_sites_on_node;
  even_sites_on_node = QDP_subset_len(QDP_even);
  odd_sites_on_node = QDP_subset_len(QDP_odd);
}

int
node_number(int x, int y, int z, int t)
{
  int c[4];

  c[0] = x;
  c[1] = y;
  c[2] = z;
  c[3] = t;
  return QDP_node_number(c);
}

int
node_index(int x, int y, int z, int t)
{
  int c[4];

  c[0] = x;
  c[1] = y;
  c[2] = z;
  c[3] = t;
  return QDP_index(c);
}
