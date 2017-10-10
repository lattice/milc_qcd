#ifndef _GENERIC_QUDA_H
#define _GENERIC_QUDA_H
/******************** generic_quda.h *********************************
*  MIMD version 7 	 				            *
*/

#include <quda_milc_interface.h>
#include "../include/openmp_defs.h"

#ifdef HAVE_QUDA
int initialize_quda(void);
#endif

static QudaMILCSiteArg_t newQudaMILCSiteArg() {
  QudaMILCSiteArg_t arg;
  arg.site = lattice;
  arg.link = NULL;
  arg.link_offset = (char*)lattice->link-(char*)lattice;
  arg.mom = NULL;
  arg.mom_offset = (char*)lattice->mom-(char*)lattice;
  arg.size = sizeof(*lattice);
  return arg;
}

#include <string.h>

static inline void fast_copy(void *dest, const void *src, size_t n) {
  memcpy(dest, src, n);
}

/*
  Allocate a pinned gauge-field array suitable for DMA transfer to the GPU
 */
static su3_matrix* create_G_quda(void) {
  return (su3_matrix*)qudaAllocatePinned(sites_on_node*4*sizeof(su3_matrix));
}

/*
  Extract the gauge field elements into a pinned array suitable for DMA transfer to the GPU
 */
static su3_matrix* create_G_from_site_quda(void) {
  su3_matrix *links = create_G_quda();
  int i;
  site *s;

  FORALLSITES_OMP(i,s,){
    fast_copy(links+4*i, s->link, 4*sizeof(su3_matrix));
  } END_LOOP_OMP

  return links;
}

/*
  Copy the momentum field elements into the site struct array
 */
static void copy_to_site_from_G_quda(su3_matrix *links) {
  int i;
  site *s;

  FORALLSITES_OMP(i,s,){
    fast_copy(s->link, links+4*i, 4*sizeof(su3_matrix));
  } END_LOOP_OMP
}

/*
  Free the pinned gauge-field array
 */
static void destroy_G_quda(su3_matrix *links) {
  qudaFreePinned(links);
}

/*
  Allocate a pinned momentum-field array suitable for DMA transfer to the GPU
 */
static anti_hermitmat* create_M_quda(void) {
  return (anti_hermitmat*)qudaAllocatePinned(sites_on_node*4*sizeof(anti_hermitmat));
}

/*
  Extract the momentum field elements into a pinned array suitable for DMA transfer to the GPU
 */
static anti_hermitmat* create_M_from_site_quda(void) {
  anti_hermitmat* momentum = create_M_quda();
  int i;
  site *s;

  FORALLSITES_OMP(i,s,){
    fast_copy(momentum+4*i, s->mom, 4*sizeof(anti_hermitmat));
  } END_LOOP_OMP

  return momentum;
}

/*
  Copy the momentum field elements into the site struct array
 */
static void copy_to_site_from_M_quda(anti_hermitmat *momentum) {
  int i;
  site *s;

  FORALLSITES_OMP(i,s,){
    fast_copy(s->mom, momentum+4*i, 4*sizeof(anti_hermitmat));
  } END_LOOP_OMP
}

/*
  Free the pinned gauge-field array
 */
static void destroy_M_quda(anti_hermitmat *momentum) {
  qudaFreePinned(momentum);
}

/*
  Return the most recent fermion link field passed to QUDA
  (defined in generic_ks/ks_multicg_offset_gpu.c)
*/
imp_ferm_links_t* get_fn_last();

/*
  Update the fermion link field passed to QUDA
  (defined in generic_ks/ks_multicg_offset_gpu.c)
*/
void set_fn_last(imp_ferm_links_t *fn_last_new);

#endif /* GENERIC_QUDA_H */
