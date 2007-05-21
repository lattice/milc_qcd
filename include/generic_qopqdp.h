#ifndef _GENERIC_QOPQDP_H
#define _GENERIC_QOPQDP_H
/******************** generic_qopqdp.h ******************************
*  MIMD version 7 	 				            *
*/

#include <qop.h>
#include "../include/su3.h"
#include <qla.h>

/* map_milc_to_qopqdp.c */

fsu3_matrix ** create_raw4_F_G (void);
fsu3_matrix ** create_raw4_F_F (void);
fsu3_vector * create_raw_F_V(void);
QLA_F3_DiracFermion * create_raw_F_D(void);

void destroy_raw4_F_G (fsu3_matrix *raw[]);
void destroy_raw4_F_F (fsu3_matrix *raw[]);
void destroy_raw_F_V (fsu3_vector *raw);
void destroy_raw_F_D (QLA_F3_DiracFermion *raw);

fsu3_matrix ** create_raw4_F_G_from_site(field_offset src, int milc_parity);
fsu3_matrix ** create_raw4_F_F_from_site(field_offset src, int milc_parity);
fsu3_vector * create_raw_F_V_from_site(field_offset src, int milc_parity);
QLA_F3_DiracFermion * create_raw_F_D_from_site(field_offset src, int milc_parity);

fsu3_matrix ** create_raw4_F_G_from_field(su3_matrix *src, int milc_parity);
fsu3_matrix ** create_raw4_F_F_from_field(anti_hermitmat *src, 
					  int milc_parity);
fsu3_vector * create_raw_F_V_from_field(su3_vector *src, int milc_parity);
QLA_F3_DiracFermion * create_raw_F_D_from_field(wilson_vector *src, int milc_parity);

void unload_raw4_F_G_to_site(field_offset dest, fsu3_matrix *raw[], 
			     int milc_parity);
void unload_raw4_F_F_to_site(field_offset dest, fsu3_matrix *raw[], 
			     int milc_parity);
void unload_raw_F_V_to_site(field_offset dest, fsu3_vector *raw, 
			    int milc_parity);
void unload_raw_F_D_to_site(field_offset dest, QLA_F3_DiracFermion *raw, 
			    int milc_parity);

void unload_raw4_F_G_to_field(su3_matrix *dest, fsu3_matrix *raw[], 
			      int milc_parity);
void unload_raw4_F_F_to_field(anti_hermitmat *dest, fsu3_matrix *raw[], 
			      int milc_parity);
void unload_raw_F_V_to_field(su3_vector *dest, fsu3_vector *raw, 
			     int milc_parity);
void unload_raw_F_D_to_field(wilson_vector *dest, QLA_F3_DiracFermion *raw, 
			     int milc_parity);

dsu3_matrix ** create_raw4_D_G (void);
void destroy_raw4_D_G (dsu3_matrix *raw[]);
dsu3_matrix ** create_raw4_D_F (void);
void destroy_raw4_D_F (dsu3_matrix *raw[]);
dsu3_vector * create_raw_D_V(void);
void destroy_raw_D_V (dsu3_vector *raw);
void destroy_raw_D_D (QLA_D3_DiracFermion *raw);
dsu3_matrix ** create_raw4_D_G_from_site(field_offset src, int milc_parity);
dsu3_matrix ** create_raw4_D_G_from_field(su3_matrix *src, int milc_parity);
dsu3_matrix ** create_raw4_D_F_from_site(field_offset src, int milc_parity);
dsu3_matrix ** create_raw4_D_F_from_field(anti_hermitmat *src, 
					  int milc_parity);
dsu3_vector * create_raw_D_V_from_site(field_offset src, int milc_parity);
dsu3_vector * create_raw_D_V_from_field(su3_vector *src, int milc_parity);
QLA_D3_DiracFermion * create_raw_D_D_from_site(field_offset src, int milc_parity);
QLA_D3_DiracFermion * create_raw_D_D_from_field(wilson_vector *src, int milc_parity);
void unload_raw4_D_G_to_site(field_offset dest, dsu3_matrix *raw[], 
			     int milc_parity);
void unload_raw4_D_G_to_field(su3_matrix *dest, dsu3_matrix *raw[], 
			      int milc_parity);
void unload_raw4_D_F_to_site(field_offset dest, dsu3_matrix *raw[], 
			     int milc_parity);
void unload_raw4_D_F_to_field(anti_hermitmat *dest, dsu3_matrix *raw[], 
			      int milc_parity);
void unload_raw_D_V_to_site(field_offset dest, dsu3_vector *raw, 
			    int milc_parity);
void unload_raw_D_V_to_field(su3_vector *dest, dsu3_vector *raw, 
			     int milc_parity);
void unload_raw_D_D_to_site(field_offset dest, QLA_D3_DiracFermion *raw, 
			    int milc_parity);
void unload_raw_D_D_to_field(wilson_vector *dest, QLA_D3_DiracFermion *raw, 
			     int milc_parity);

#endif /* GENERIC_QOPQDP_H */
