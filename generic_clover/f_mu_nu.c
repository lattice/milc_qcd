/******************* f_mu_nu.c ****************************************/

/* MIMD version 6 */
/* This version uses gathers to get the neighbors */
/* version of 1/4/95 by UMH */
/* 1/24/00 combined with Schroedinger functional version - UMH */

/* Compute F_{mu,nu} used in the clover fermion action. */

/*  f_mn[site] = (G_mn - G_mn^dagger)/8  (antihermitian) where

    G_mn is given by the product of link matrices in the pattern shown.

          <-----^                  <-----^ 
          |     |		   |     | 
          |     |		   |     | 
  G_mn =  O----->   +  <-----O  +  v---->O  +  O<----^ 
		       |     | 		       |     | 
                       |     | 		       |     | 
     nu  ^             v-----> 		       v-----> 
         |       
         O---> mu

    The "site" is indicated by O.

    If we write U_mu = exp[i A_mu] then

    f_mn =approx= i F_{mu,nu}
  
 */

#include "generic_clover_includes.h"

void f_mu_nu(field_offset f_mn,int mu,int nu) {
register int i;
register site *s;
msg_tag *tag0,*tag1,*tag2,*tag3,*tag4;
su3_matrix tmat4;
int disp[4];	/* displacement vector for general gather */
int order_flag;

field_offset tempmat1,tempmat2,tempmat3;
tempmat1 = F_OFFSET(mp);
tempmat2 = F_OFFSET(tmp);
tempmat3 = F_OFFSET(htmp[0]);

    /* Want mu < nu, so that only nu can be TUP! */
    if(mu > nu){
	i = mu;
	mu = nu;
	nu = i;
	order_flag = 1;
    }
    else{
	order_flag = 0;
    }

    /* get link[nu] from direction +mu */
     tag0 = start_gather( F_OFFSET(link[nu]), sizeof(su3_matrix),
	 mu, EVENANDODD, gen_pt[0] );

     /* get link[mu] from direction +nu */
     tag1 = start_gather( F_OFFSET(link[mu]), sizeof(su3_matrix),
	 nu, EVENANDODD, gen_pt[1] );

     /* Make one corner with link[nu]^dagger link[mu] */
     FORALLSITES(i,s){
	 /* mult_su3_an( &(s->link[nu]), &(s->link[mu]), &(s->tempmat1) ); */
	 mult_su3_an( &(s->link[nu]), &(s->link[mu]),
	    ((su3_matrix *)F_PT(s,tempmat1)) );
    }

    /* Make one corner with link[nu](x+mu) link[mu](x+nu)^dagger
       and multiply the two corners together in the two different ways */
    /* Note f_mn is here used as a temporary! */
    wait_gather(tag0);
    wait_gather(tag1);
    FORALLSITES(i,s){
#ifdef SCHROED_FUN
	if(s->t==(nt-1) && nu==TUP){
	    mult_su3_na( (su3_matrix *)(gen_pt[0][i]),
		&(s->boundary[mu]), ((su3_matrix *)F_PT(s,f_mn)) );
	}
	else{
	    mult_su3_na( (su3_matrix *)(gen_pt[0][i]),
		(su3_matrix *)(gen_pt[1][i]), ((su3_matrix *)F_PT(s,f_mn)) );
	}
#else
        mult_su3_na( (su3_matrix *)(gen_pt[0][i]),
	    (su3_matrix *)(gen_pt[1][i]), ((su3_matrix *)F_PT(s,f_mn)) );
#endif
	mult_su3_nn( ((su3_matrix *)F_PT(s,tempmat1)),
	    ((su3_matrix *)F_PT(s,f_mn)),
	    ((su3_matrix *)F_PT(s,tempmat2)) );
	mult_su3_nn( ((su3_matrix *)F_PT(s,f_mn)),
	    ((su3_matrix *)F_PT(s,tempmat1)),
	    ((su3_matrix *)F_PT(s,tempmat3)) );
    }

    /* tempmat2 is the plaquette +mu -nu and must be gathered from -nu */
    tag2 = start_gather( tempmat2, sizeof(su3_matrix),
	OPP_DIR(nu), EVENANDODD, gen_pt[2] );

    /* tempmat3 is the plaquette -mu +nu and must be gather from -mu */
    tag3 = start_gather( tempmat3, sizeof(su3_matrix),
	OPP_DIR(mu), EVENANDODD, gen_pt[3] );

    /* Now make +mu +nu plaquette and put in f_mn */
    FORALLSITESDOMAIN(i,s){
        mult_su3_nn( &(s->link[mu]), ((su3_matrix *)F_PT(s,f_mn)), &tmat4 );
        mult_su3_na( &tmat4, &(s->link[nu]), ((su3_matrix *)F_PT(s,f_mn)) );
    }

    /* Now gather +mu -nu plaquette and add to f_mn */
    wait_gather(tag2);
    FORALLSITESDOMAIN(i,s){
	add_su3_matrix( ((su3_matrix *)F_PT(s,f_mn)),
	    (su3_matrix *)(gen_pt[2][i]), ((su3_matrix *)F_PT(s,f_mn)) );
    }

    /* For Schroedinger functional, tempmat2 is wrong for t=nt-1 if nu=TUP,
       but it will not be used after being moved to t=0! */
    FORALLSITES(i,s){
	mult_su3_an( (su3_matrix *)(gen_pt[1][i]),
	    ((su3_matrix *)F_PT(s,tempmat1)), &tmat4 );
	mult_su3_nn( &tmat4, (su3_matrix *)(gen_pt[0][i]),
	    ((su3_matrix *)F_PT(s,tempmat2)) );
    }

    /* tempmat2 is now plaquette -mu -nu and must be gathered with
       displacement -mu-nu */
    for(i=XUP;i<=TUP;i++)disp[i]=0;
    disp[mu] = -1;
    disp[nu] = -1;
    tag4 = start_general_gather( tempmat2, sizeof(su3_matrix),
	disp, EVENANDODD, gen_pt[4] );

    /* Now gather -mu +nu plaquette and add to f_mn */
    wait_gather(tag3);
    FORALLSITESDOMAIN(i,s){
	add_su3_matrix( ((su3_matrix *)F_PT(s,f_mn)),
	    (su3_matrix *)(gen_pt[3][i]), ((su3_matrix *)F_PT(s,f_mn)) );
    }
    cleanup_gather(tag0);
    cleanup_gather(tag1);
    cleanup_gather(tag2);
    cleanup_gather(tag3);

    /* Finally gather -mu -nu plaquette and add to f_mn */
    wait_general_gather(tag4);
    FORALLSITESDOMAIN(i,s){
	add_su3_matrix( ((su3_matrix *)F_PT(s,f_mn)),
	    (su3_matrix *)(gen_pt[4][i]), ((su3_matrix *)F_PT(s,f_mn)) );
    }
    cleanup_general_gather(tag4);

    /* F_{mu,nu} is now 1/8 of f_mn - f_mn^dagger */
    FORALLSITESDOMAIN(i,s){
	su3_adjoint( ((su3_matrix *)F_PT(s,f_mn)), &tmat4 );
	if(order_flag == 0){
	    sub_su3_matrix( ((su3_matrix *)F_PT(s,f_mn)), &tmat4, &tmat4 );
	}
	else{
	    sub_su3_matrix( &tmat4, ((su3_matrix *)F_PT(s,f_mn)), &tmat4 );
	}
	scalar_mult_su3_matrix( &tmat4, 0.125, ((su3_matrix *)F_PT(s,f_mn)) );
    }

} /* f_mu_nu */

