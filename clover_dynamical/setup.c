/******** setup.c *********/
/* MIMD version 6 */
/* Version for dynamical clover fermions with Symanzik/tadople
	improved gauge field	*/

/* Modifications ...
   Original 1995 by Matt Wingate.
   6/6/98 Version 5 port CD
*/

#include "cl_dyn_includes.h"
#include <string.h>
int initial_set();
#define IF_OK if(status==0)

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int setup()   {
  int prompt;

  /* print banner, get volume, nflavors, seed */
  prompt=initial_set();
  /* initialize the node random number generator */
  initialize_prn(&node_prn,iseed,volume+mynode());
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  /* set up neighbor pointers and comlink structures */
  make_nn_gathers();

  return(prompt);
}


/* SETUP ROUTINES */
int initial_set(){
  int prompt,status;
  /* On node zero, read lattice size, seed, nflavors and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 with Clover fermions\n");
    printf("Microcanonical simulation with refreshing\n");
    printf("MIMD version 6\n");
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
#ifdef HMC_ALGORITHM
    printf("Hybrid Monte Carlo algorithm\n");
#endif
#ifdef PHI_ALGORITHM
    printf("PHI algorithm\n");
#else
    printf("R algorithm\n");
#endif
#ifdef SPECTRUM
    printf("With spectrum measurements\n");
#endif
    time_stamp("start");

    status = get_prompt( &prompt );

    IF_OK status += get_i(prompt,"nflavors", &par_buf.nflavors );
#ifdef PHI_ALGORITHM
    IF_OK if(par_buf.nflavors != 2){
	    printf("Dummy! Use phi algorithm only for two flavors\n");
	    terminate(-1);
	}
#endif

    IF_OK status += get_i(prompt,"nx", &par_buf.nx );
    IF_OK status += get_i(prompt,"ny", &par_buf.ny );
    IF_OK status += get_i(prompt,"nz", &par_buf.nz );
    IF_OK status += get_i(prompt,"nt", &par_buf.nt );

    IF_OK status += get_i(prompt,"iseed", &par_buf.iseed );

    if(status>0) par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if(mynode()==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  
  if( par_buf.stopflag != 0 )
    normal_exit(0);


  nflavors=par_buf.nflavors;
  nx=par_buf.nx;
  ny=par_buf.ny;
  nz=par_buf.nz;
  nt=par_buf.nt;
  iseed=par_buf.iseed;

  this_node = mynode();
  number_of_nodes = numnodes();
  volume=nx*ny*nz*nt;
  total_iters=0;
  return(prompt);
}
    
/* read in parameters and coupling constants	*/
int readin(int prompt) {
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/

  int status;
  Real x;

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0){
    
    printf("\n\n");

    status = 0;
    
    /* warms, trajecs */
    IF_OK status += get_i(prompt,"warms", &par_buf.warms);
    IF_OK status += get_i(prompt,"trajecs", &par_buf.trajecs);
    
    /* trajectories between propagator measurements */
    IF_OK status += get_i(prompt,"traj_between_meas", 
			  &par_buf.propinterval);
    
    /* get couplings and broadcast to nodes	*/
    /* beta, kappa */
    IF_OK status += get_f(prompt,"beta", &par_buf.beta);
    IF_OK status += get_f(prompt,"kappa", &par_buf.kappa);
    
    /* Clover coefficient, u0 */
    IF_OK status += get_f(prompt,"clov_c", &par_buf.clov_c );
    IF_OK status += get_f(prompt,"u0", &par_buf.u0 );

    /* microcanonical time step */
    IF_OK status += get_f(prompt,"microcanonical_time_step", &par_buf.epsilon );
    
    /*microcanonical steps per trajectory */
    IF_OK status += get_i(prompt,"steps_per_trajectory", &par_buf.steps);

    /* maximum no. of conjugate gradient iterations */
    IF_OK status += get_i(prompt,"max_cg_iterations", &par_buf.niter);

    /* maximum no. of conjugate gradient restarts */
    IF_OK status += get_i(prompt,"max_cg_restarts", &par_buf.nrestart );
    
    /* error per site for conjugate gradient */
    IF_OK {
      status += get_f(prompt,"error_per_site", &x );
      par_buf.rsqmin = x*x;
    }

    /* error for propagator conjugate gradient */
    IF_OK {
      status += get_f(prompt,"error_for_propagator", &x );
      par_buf.rsqprop = x*x;
    }
    
    /* find out what kind of starting lattice to use */
    IF_OK status += ask_starting_lattice( prompt, &par_buf.startflag,
	par_buf.startfile );

    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_lattice( prompt, &(par_buf.saveflag),
			     par_buf.savefile );
    
    if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if(this_node==0) */
  
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  if( par_buf.stopflag != 0 )
    normal_exit(0);

  
  warms = par_buf.warms;
  trajecs = par_buf.trajecs;
  steps = par_buf.steps;
  propinterval = par_buf.propinterval;
  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  niter = par_buf.niter;
  nrestart = par_buf.nrestart;
  rsqmin = par_buf.rsqmin;
  rsqprop = par_buf.rsqprop;
  epsilon = par_buf.epsilon;
  beta = par_buf.beta;
  kappa = par_buf.kappa;
  clov_c = par_buf.clov_c;
  u0 = par_buf.u0;
  strcpy(startfile,par_buf.startfile);
  strcpy(savefile,par_buf.savefile);

  /* Load part of inversion control structure for generic inverters */
  qic.min = 0;
  qic.max = niter;
  qic.nrestart = nrestart;
  qic.resid = sqrt(rsqprop);  /* Note different convention for error */

  /* Load part of Dirac matrix parameters */
  dcp.Kappa = kappa;
  dcp.Clov_c = clov_c;
  dcp.U0 = u0;

  /* Do whatever is needed to get lattice */
  startlat_p = reload_lattice( startflag, startfile );
  
  return(0);
}

