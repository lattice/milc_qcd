/************************ setup.c ****************************/
/* MIMD version 7 */
/*			    -*- Mode: C -*-
// File: setup.c
// Created: Fri Aug  4 1995
// Authors: J. Hetrick & K. Rummukainen
// Modified for general improved action 5/24/97  DT
//
// Description: Setup routines for improved fermion lattices
//              Includes lattice structures for Naik imroved 
//              staggered Dirac operator
//         Ref: S. Naik, Nucl. Phys. B316 (1989) 238
//              Includes a parameter prompt for Lepage-Mackenzie 
//              tadpole improvement
//         Ref: Phys. Rev. D48 (1993) 2250
//  $Log: setup.c,v $
//  Revision 1.6  2006/09/01 03:58:08  detar
//  Updating error tolerances and fiducial samples.  Removing junk files.
//
//  Revision 1.5  2006/08/25 04:55:41  detar
//  Groom to remove declarations of unused variables and do strict initialization
//
//  Revision 1.4  2006/05/17 17:31:53  detar
//  Add links to Ludmila's chemical potential code.
//  Add npbp_reps to the su3_rmd_eos and su3_rmd_mu_eos sample inputs
//
//  Revision 1.3  2005/11/10 16:58:44  detar
//  Experimenting with tags and versions
//
//  Revision 1.2  2005/11/10 16:55:09  detar
//  Add CVS version and log tags to file
//
//
*/
/* MIMD version 7 */
#define IF_OK if(status==0)

#include "ks_imp_includes.h"	/* definitions files and prototypes */

EXTERN gauge_header start_lat_hdr;
gauge_file *gf;

gauge_file *r_parallel_i(char *);
void r_parallel(gauge_file *, field_offset);
void r_parallel_f(gauge_file *);

gauge_file *r_binary_i(char *);
void r_binary(gauge_file *);
void r_binary_f(gauge_file *);
void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int
setup()
{
  int initial_set();
  void make_3n_gathers();
#ifdef HAVE_QDP
  int i;
#endif
  int prompt;
  
  /* print banner, get volume, nflavors1,nflavors2, nflavors, seed */
  prompt = initial_set();
  /* initialize the node random number generator */
  initialize_prn( &node_prn, iseed, volume+mynode() );
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();

  /* Mark t_longlink and t_fatlink as unallocated */
  t_longlink = NULL;  
  t_fatlink = NULL;
  
#ifdef DBLSTORE_FN
  t_longbacklink = NULL;
  t_fatbacklink = NULL;
#endif

#ifdef DM_DU0
  /* Allocate space for u0 derivative of temp fatlinks */
  t_dfatlink_du0 = (su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));
  if(t_dfatlink_du0==NULL){
    printf("NODE %d: no room for t_dfatlink\n",this_node);
    terminate(1);
  }
#endif

  node0_printf("Made lattice\n"); fflush(stdout);
  /* set up neighbor pointers and comlink structures
     code for this routine is in com_machine.c  */
  make_nn_gathers();
  node0_printf("Made nn gathers\n"); fflush(stdout);
  /* set up 3rd nearest neighbor pointers and comlink structures
     code for this routine is below  */
  make_3n_gathers();
  node0_printf("Made 3nn gathers\n"); fflush(stdout);
  /* set up K-S phase vectors, boundary conditions */
  phaseset();
  
#ifdef HAVE_QDP
  for(i=0; i<8; i++) {
    implinks[i] = QDP_create_M();
  }
  fatlinks = implinks;
  longlinks = implinks + 4;
  for(i=0; i<4; ++i) {
    shiftdirs[i] = QDP_neighbor[i];
    shiftdirs[i+4] = neighbor3[i];
  }
  for(i=0; i<8; ++i) {
    shiftfwd[i] = QDP_forward;
    shiftbck[i] = QDP_backward;
  }
#endif

  node0_printf("Finished setup\n"); fflush(stdout);
  return( prompt );
}


/* SETUP ROUTINES */
int 
initial_set()
{
  int prompt,status;
  /* On node zero, read lattice size, seed, nflavors1, nflavors2,
     nflavors, and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 with improved KS action\n");
    printf("Microcanonical simulation with refreshing\n");
    printf("MIMD version 7 $Name:  $\n");
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
    status=get_prompt(&prompt);
#ifdef ONEMASS
    IF_OK status += get_i(prompt,"nflavors", &par_buf.nflavors );
#else
    IF_OK status += get_i(prompt,"nflavors1", &par_buf.nflavors1 );
    IF_OK status += get_i(prompt,"nflavors2", &par_buf.nflavors2 );
#endif
#ifdef PHI_ALGORITHM
#ifdef ONEMASS
    IF_OK if(par_buf.nflavors != 4){
      printf("Dummy! Use phi algorithm only for four flavors\n");
      status++;
    }
#else
    IF_OK if( par_buf.nflavors1 != 4 || par_buf.nflavors2 != 4 ){
      printf("Dummy! Use phi algorithm only for four flavors\n");
      status++;
    }
#endif
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
  
  nx=par_buf.nx;
  ny=par_buf.ny;
  nz=par_buf.nz;
  nt=par_buf.nt;
  iseed=par_buf.iseed;
#ifdef ONEMASS
  nflavors=par_buf.nflavors;
#else
  nflavors1=par_buf.nflavors1;
  nflavors2=par_buf.nflavors2;
#endif
  
  this_node = mynode();
  number_of_nodes = numnodes();
  volume=nx*ny*nz*nt;
  total_iters=0;
  return(prompt);
}

/* read in parameters and coupling constants	*/
int
readin(int prompt)
{
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/
  
  int status;
  Real x;
#ifdef SPECTRUM
  int i;
  char request_buf[MAX_SPECTRUM_REQUEST];
#endif
  
  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0) {
    
    printf("\n\n");
    status=0;
    
    /* warms, trajecs */
    IF_OK status += get_i(prompt,"warms", &par_buf.warms );
    IF_OK status += get_i(prompt,"trajecs", &par_buf.trajecs );
    
    /* trajectories between propagator measurements */
    IF_OK status += 
      get_i(prompt,"traj_between_meas", &par_buf.propinterval );
    
    /* get couplings and broadcast to nodes	*/
    /* beta, mass1, mass2 or mass */
    IF_OK status += get_f(prompt,"beta", &par_buf.beta );
#ifdef ONEMASS
    IF_OK status += get_f(prompt,"mass", &par_buf.mass );
#else
    IF_OK status += get_f(prompt,"mass1", &par_buf.mass1 );
    IF_OK status += get_f(prompt,"mass2", &par_buf.mass2 );
#endif
    IF_OK status += get_f(prompt,"u0", &par_buf.u0 );
    
    /* microcanonical time step */
    IF_OK status += 
      get_f(prompt,"microcanonical_time_step", &par_buf.epsilon );
    
    /*microcanonical steps per trajectory */
    IF_OK status += get_i(prompt,"steps_per_trajectory", &par_buf.steps );
    
    /* maximum no. of conjugate gradient iterations */
    IF_OK status += get_i(prompt,"max_cg_iterations", &par_buf.niter );
    
    /* error per site for conjugate gradient */
    IF_OK status += get_f(prompt,"error_per_site", &x );
    IF_OK par_buf.rsqmin = x*x;   /* rsqmin is r**2 in conjugate gradient */
    /* New conjugate gradient normalizes rsqmin by norm of source */
    
    /* error for propagator conjugate gradient */
    IF_OK status += get_f(prompt,"error_for_propagator", &x );
    IF_OK par_buf.rsqprop = x*x;
    
#ifdef NPBP_REPS
    /* number of random sources npbp_reps */
    IF_OK status += get_i(prompt,"npbp_reps", &par_buf.npbp_reps_in );
#endif
    
#ifdef SPECTRUM
    /* request list for spectral measurments */
    /* prepend and append a comma for ease in parsing */
    IF_OK status += get_s(prompt,"spectrum_request", request_buf );
    IF_OK strcpy(par_buf.spectrum_request,",");
    IF_OK strcat(par_buf.spectrum_request,request_buf);
    IF_OK strcat(par_buf.spectrum_request,",");
    
    /* source time slice and increment */
    IF_OK status += get_i(prompt,"source_start", &par_buf.source_start );
    IF_OK status += get_i(prompt,"source_inc", &par_buf.source_inc );
    IF_OK status += get_i(prompt,"n_sources", &par_buf.n_sources );
    
    /* Additional parameters for spectrum_multimom */
    if(strstr(par_buf.spectrum_request,",spectrum_multimom,") != NULL){
      IF_OK status += get_i(prompt,"spectrum_multimom_nmasses",
			    &par_buf.spectrum_multimom_nmasses );
      IF_OK status += get_f(prompt,"spectrum_multimom_low_mass",
			    &par_buf.spectrum_multimom_low_mass );
      IF_OK status += get_f(prompt,"spectrum_multimom_mass_step",
			    &par_buf.spectrum_multimom_mass_step );
    }
    /* Additional parameters for fpi */
    par_buf.fpi_nmasses = 0;
    if(strstr(par_buf.spectrum_request,",fpi,") != NULL){
      IF_OK status += get_i(prompt,"fpi_nmasses",
			    &par_buf.fpi_nmasses );
      if(par_buf.fpi_nmasses > MAX_FPI_NMASSES){
	printf("Maximum of %d exceeded.\n",MAX_FPI_NMASSES);
	terminate(1);
      }
      for(i = 0; i < par_buf.fpi_nmasses; i++){
	IF_OK status += get_f(prompt,"fpi_mass",
			      &par_buf.fpi_mass[i]);
      }
    }
    
#endif /*SPECTRUM*/
    
    /* find out what kind of starting lattice to use */
    IF_OK status += ask_starting_lattice( prompt, &(par_buf.startflag),
					  par_buf.startfile );
    
    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_lattice( prompt, &(par_buf.saveflag),
					par_buf.savefile );
    IF_OK status += ask_ildg_LFN( prompt, par_buf.saveflag,
				  par_buf.stringLFN );
    
    if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if(this_node==0) */
  
    /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  
  if( par_buf.stopflag != 0 )return par_buf.stopflag;
  
  warms = par_buf.warms;
  trajecs = par_buf.trajecs;
  steps = par_buf.steps;
  propinterval = par_buf.propinterval;
  niter = par_buf.niter;
  npbp_reps_in = par_buf.npbp_reps_in;
  rsqmin = par_buf.rsqmin;
  rsqprop = par_buf.rsqprop;
  epsilon = par_buf.epsilon;
  beta = par_buf.beta;
#ifdef ONEMASS
  mass = par_buf.mass;
#else
  mass1 = par_buf.mass1;
  mass2 = par_buf.mass2;
#endif
  u0 = par_buf.u0;
#ifdef SPECTRUM
  strcpy(spectrum_request,par_buf.spectrum_request);
  source_start = par_buf.source_start;
  source_inc = par_buf.source_inc;
  n_sources = par_buf.n_sources;
  spectrum_multimom_nmasses = par_buf.spectrum_multimom_nmasses;
  spectrum_multimom_low_mass = par_buf.spectrum_multimom_low_mass;
  spectrum_multimom_mass_step = par_buf.spectrum_multimom_mass_step;
  fpi_nmasses = par_buf.fpi_nmasses;
  for(i = 0; i < fpi_nmasses; i++){
    fpi_mass[i] = par_buf.fpi_mass[i];
  }
#endif /*SPECTRUM*/
  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile,par_buf.startfile);
  strcpy(savefile,par_buf.savefile);
  strcpy(stringLFN, par_buf.stringLFN);
  
  /* Do whatever is needed to get lattice */
  if( startflag == CONTINUE ){
    rephase( OFF );
  }
  startlat_p = reload_lattice( startflag, startfile );
  /* if a lattice was read in, put in KS phases and AP boundary condition */
  valid_fatlinks = valid_longlinks = 0;
  phases_in = OFF;
  rephase( ON );
  
  /* make table of coefficients and permutations of loops in gauge action */
  make_loop_table();
  /* make table of coefficients and permutations of paths in quark action */
  make_path_table();
  
  return(0);
}

/* Set up comlink structures for 3rd nearest gather pattern; 
   make_lattice() and  make_nn_gathers() must be called first, 
   preferably just before calling make_3n_gathers().
*/
void 
make_3n_gathers()
{
  int i;
#ifdef HAVE_QDP
  int disp[4]={0,0,0,0};
#endif
  
  for(i=XUP; i<=TUP; i++) {
    make_gather(third_neighbor, &i, WANT_INVERSE,
		ALLOW_EVEN_ODD, SWITCH_PARITY);
  }
  
  /* Sort into the order we want for nearest neighbor gathers,
     so you can use X3UP, X3DOWN, etc. as argument in calling them. */
  
  sort_eight_gathers(X3UP);

#ifdef HAVE_QDP
  for(i=0; i<4; i++) {
    disp[i] = 3;
    neighbor3[i] = QDP_create_shift(disp);
    disp[i] = 0;
  }
#endif
}

/* this routine uses only fundamental directions (XUP..TDOWN) as directions */
/* returning the coords of the 3rd nearest neighbor in that direction */

void 
third_neighbor(int x, int y, int z, int t, int *dirpt, int FB,
	       int *xp, int *yp, int *zp, int *tp)
     /* int x,y,z,t,*dirpt,FB;  coordinates of site, direction (eg XUP), and
	"forwards/backwards"  */
     /* int *xp,*yp,*zp,*tp;    pointers to coordinates of neighbor */
{
  int dir;
  dir = (FB==FORWARDS) ? *dirpt : OPP_DIR(*dirpt);
  *xp = x; *yp = y; *zp = z; *tp = t;
  switch(dir){
  case XUP: *xp = (x+3)%nx; break;
  case XDOWN: *xp = (x+4*nx-3)%nx; break;
  case YUP: *yp = (y+3)%ny; break;
  case YDOWN: *yp = (y+4*ny-3)%ny; break;
  case ZUP: *zp = (z+3)%nz; break;
  case ZDOWN: *zp = (z+4*nz-3)%nz; break;
  case TUP: *tp = (t+3)%nt; break;
  case TDOWN: *tp = (t+4*nt-3)%nt; break;
  default: printf("third_neighb: bad direction\n"); exit(1);
  }
}
