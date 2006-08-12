/************************ grsource_rhmc.c *****************************/
/* MIMD version 7 */
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.  Assumes "dslash" has
   been defined to be the appropriate "dslash_fn" or "dslash_eo"
*/
#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include "../include/dslash_ks_redefine.h"

/* construct a gaussian random vector, g_rand, and phi=M(dagger)*M)^(nf/8)*g_rand  */
/* "parity" is EVEN, ODD, or EVENANDODD.  The parity is the parity at
    which phi is computed.  g_rand must always be computed at all sites. */

// residues, roots and order define rational function approximation for
 // x^(nf/8)
void grsource_imp_rhmc( field_offset dest, Real mass, Real *residues, Real *roots,
  int order, int parity, su3_vector **multi_x, su3_vector *sumvec ) {
register int i,j;
register site *s;
Real final_rsq;
/*TEMP*/ double sum;

    sum=0.0;
    FORSOMEPARITY(i,s,parity){
        for(j=0;j<3;j++){
#ifdef SITERAND
            s->g_rand.c[j].real = gaussian_rand_no(&(s->site_prn));
            s->g_rand.c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
            s->g_rand.c[j].real = gaussian_rand_no(&node_prn);
            s->g_rand.c[j].imag = gaussian_rand_no(&node_prn);
#endif
        }
	/*TEMP*/ sum += (double)magsq_su3vec( &(s->g_rand) );
    }
    /*TEMP*/g_doublesum( &sum);  node0_printf("GRSOURCE: sum = %.10e\n",sum);
#ifdef FN
    if( !valid_longlinks )load_longlinks();
    if( !valid_fatlinks )load_fatlinks();
#endif
    ks_ratinv( F_OFFSET(g_rand), multi_x, mass, roots, order, niter, rsqmin, parity, &final_rsq );
    ks_rateval( sumvec, F_OFFSET(g_rand), multi_x, residues, order, parity );
    FORSOMEPARITY(i,s,parity){ *(su3_vector *)F_PT(s,dest) = sumvec[i]; }
}/* grsource_rhmc */

// Old grsource, still used by measurement routines
void grsource_imp( field_offset dest, Real mass, int parity) {
register int i,j;
register site *s;
    FORALLSITES(i,s){
        for(j=0;j<3;j++){
#ifdef SITERAND
            s->g_rand.c[j].real = gaussian_rand_no(&(s->site_prn));
            s->g_rand.c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
            s->g_rand.c[j].real = gaussian_rand_no(&node_prn);
            s->g_rand.c[j].imag = gaussian_rand_no(&node_prn);
#endif
        }
    }
#ifdef FN
    if( !valid_longlinks )load_longlinks();
    if( !valid_fatlinks )load_fatlinks();
#endif
    dslash_site( F_OFFSET(g_rand), dest, parity);
    scalar_mult_latvec( dest, -1.0, dest, parity );
    scalar_mult_add_latvec( dest, F_OFFSET(g_rand), 2.0*mass,
	dest, parity );
}/* grsource */


/* construct a plain gaussian random vector in the site structure  */
/* "parity" is EVEN, ODD, or EVENANDODD. */

void grsource_plain( field_offset dest, int parity ) {
  int i,j;
  site *s;
  su3_vector *rand;
  FORSOMEPARITY(i,s,parity){
    for(j=0;j<3;j++){
#ifdef SITERAND
      rand = (su3_vector *)F_PT(s,dest);
            rand->c[j].real = gaussian_rand_no(&(s->site_prn));
            rand->c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
            rand->c[j].real = gaussian_rand_no(&node_prn);
            rand->c[j].imag = gaussian_rand_no(&node_prn);
#endif
    }
  }    
}/* grsource */


/* Check congrad by multiplying src by M, compare result to g_rand */
/* Before calling checkmul() you should call grsource(EVENANDODD) and
   congrad(...,EVENANDODD) */
void checkmul_imp( field_offset src, Real mass ) {
register int i,j;
register site *s;
    dslash_site( src, F_OFFSET(ttt), EVENANDODD);
    scalar_mult_add_latvec( F_OFFSET(ttt), src, 2.0*mass,
	F_OFFSET(ttt), EVENANDODD );
    FORALLSITES(i,s){
	printf("Site %d %d %d %d\n",s->x,s->y,s->z,s->t);
	for(j=0;j<3;j++){
	    printf("%d %d\t%e\t%e\t%e\n",i,j,(double)s->g_rand.c[j].real,
		(double)s->ttt.c[j].real,(double)s->g_rand.c[j].real - (double)s->ttt.c[j].real);
	    printf("%d %d\t%e\t%e\t%e\n",i,j,(double)s->g_rand.c[j].imag,
		(double)s->ttt.c[j].imag,(double)s->g_rand.c[j].imag - (double)s->ttt.c[j].imag);
	}
	printf("\n");
/**sleep(2);**/
    }
}/* checkmul */