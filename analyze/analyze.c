/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#define _FILE_OFFSET_BITS 64
#define _THREAD_SAFE
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#include "../compat.h"
#include "../graph500.h"
#include "../options.h"

static int64_t maxvtx, maxdeg, nIJ, mindeg, maxpart, totalpart;
static const struct packed_edge * restrict IJ;
static int64_t * restrict head, * restrict deg, * restrict next;

int 
create_graph_from_edgelist (struct packed_edge *IJ_in, int64_t nedge)
{
  int err = 0;

  IJ = IJ_in;
  nIJ = nedge;
  maxvtx = -1;
  maxdeg = -1;
  mindeg = INT64_MAX;
  maxpart = 0;
  totalpart = 0;

  int64_t k;
  for (k = 0; k < nedge; ++k) {
    if (get_v0_from_edge(&IJ[k]) > maxvtx)
      maxvtx = get_v0_from_edge(&IJ[k]);
    if (get_v1_from_edge(&IJ[k]) > maxvtx)
      maxvtx = get_v1_from_edge(&IJ[k]);
  }

  head = malloc ((2*(maxvtx+1) + 2*nIJ) * sizeof (int64_t));
  if (!head) return -1;
  deg = &head[maxvtx+1];
  next = &deg[maxvtx+1];

  for (k = 0; k <= maxvtx; ++k) {
    head[k] = -1;
    deg[k] = 0;
  }

  for (k = 0; k < nedge; ++k) {
    const int64_t i = get_v0_from_edge(&IJ[k]);
    const int64_t j = get_v1_from_edge(&IJ[k]);
    int64_t t_head, t;

    if (i >= 0 && j >= 0 && i != j) {
      next[2*k] = -1;
      next[1+2*k] = -1;
      t = 2*k+1; /* Point at the *other* end. */
      t_head = head[i];
      head[i] = t;
      assert (t_head < 2*nIJ);
      next[t] = t_head;
      ++deg[i];

      --t;
      t_head = head[j];
      head[j] = t;
      assert (t_head < 2*nIJ);
      next[t] = t_head;
      ++deg[j];
    }
  }

  assert (nodes > 0);
  if (nodes < 2) {
    printf("# v deg(v)\n");
  }
  for (int64_t kg = 0; kg <= maxvtx; ++kg) {
    if (nodes < 2) {
      printf("%lu %lu\n", kg, deg[kg]); 
    }
    if (deg[kg] > maxdeg)
      maxdeg = deg[kg];
    if (deg[kg] < mindeg)
      mindeg = deg[kg];
  }

  int64_t vpart = maxvtx/nodes;

  if (nodes >= 2) {
    printf("# partition deg\n");
  }
  for (int64_t p = 0; p < nodes; ++p) {
    int64_t edges = 0;
    for (int64_t v = 0; v < vpart; ++v) {
      edges += deg[(p*vpart)+v];
    }
    if (edges > maxpart)
      maxpart = edges;
    totalpart += edges;
    if (nodes >= 2) {
      printf("%lu %lu\n", p, edges); 
    }
  }
  printf("# maxvtx = %lu maxdeg = %lu mindeg = %lu nedges = %lu\n",
         maxvtx, maxdeg, mindeg, nedge);
  printf("# partitions = %lu vert/part = %lu\n", nodes, vpart);
  int64_t avgpart = totalpart/nodes;
  printf("# maxpart = %lu avgpart = %lu imbalance = %lf\n",
         maxpart, avgpart, (float)maxpart/avgpart);
  return err;
}

int
make_bfs_tree (int64_t *bfs_tree_out, int64_t *max_vtx_out,
	       int64_t srcvtx)
{
 return 0;
}

void
destroy_graph (void)
{
  free (head);
}
