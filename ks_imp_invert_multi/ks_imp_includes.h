/****************** ks_imp_includes.h ******************************/
/*
*  Include files for Kogut-Susskind dynamical improved action application
*/

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "lattice.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic_ks.h"
#include "../include/generic.h"
#include "../include/dirs.h"
#include "../include/io_prop_ks.h"

#ifdef FN
#define dslash dslash_fn
#endif
#ifdef EO
#define dslash dslash_eo
#endif

/* prototypes for functions in high level code */
int setup();
int readin(int prompt);
int update();
void update_h( Real eps );
void update_u( Real eps );
double hmom_action( );
double fermion_action( );


void g_measure( void );
void gauge_field_copy(field_offset src,field_offset dest);
