/* plaquette + 1x1x1 ("chair") loops, with two reps. */
#ifdef GAUGE_ACTION_PART1
/* defines NREPS NLOOP MAX_LENGTH MAX_NUM */
#define NREPS 2
#define NLOOP 2
#define MAX_LENGTH 6
#define MAX_NUM 16
#endif

#ifdef GAUGE_ACTION_PART2
    static int loop_ind[NLOOP][MAX_LENGTH] = {
    { XUP, YUP, XDOWN, YDOWN, NODIR, NODIR },
    { XUP, YUP, ZUP, XDOWN, YDOWN , ZDOWN},
    };
    static int loop_length_in[NLOOP] = {4,6};

    for(j=0;j<NLOOP;j++){
	loop_num[j] = 0;
	loop_length[j] = loop_length_in[j];
	for(i=0;i<NREPS;i++){
	    loop_coeff[j][i] = 0.0;
	}
    }

    /* set up the loop coefficients */
    loop_coeff[0][0]= 0.5;
    loop_coeff[0][1]= 0.1;
    loop_coeff[1][0]= 0.1;
    loop_coeff[1][1]= 0.04;
    strcpy(gauge_action_description,"\"1x1 plus chair, two reps\"");
    node0_printf("1x1 plus chair, two reps\n");
#endif
