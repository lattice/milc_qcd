/* 
 * Andreas Stathopoulos
 * binaryRecursiveColoring.cc
 * This file was originially written by Andreas Stathopoulos and adopted for use in milc_qcd by Heather Maria Switzer
 * 
 * Copyright (c) 2015, College of William & Mary All rights reserved.
 * 
 * Andreas Stathopoulos, Jesse Laeuchli, and Kostas Orginos. Hierarchical probing for estimating the matrix of the matrix inverse on toroidal lattices, 2013, SIAM J. Sci. Comput., 35(5), S299S322, arXiv:1302.4018 [hep-lat]
 * 
 */

// This is Andreas's file that has been modified slightly for chroma.
// Declarations have been put in the header file
// Heather Maria Switzer adapted from Arjun Singh Gambhir
#include "binaryRecursiveColoring.h"

/* The struct typedef is in the header
 ***********************************************************
 *typedef struct mashVars{
 *	unsigned int d;			// Number of dimensions
 *	unsigned int *ptsPerDim;	// Points per dimension
 *	unsigned int dmax;		// Which dimension has most points
 * 	unsigned int twoTod;		// 2^d
 *	unsigned int *logPts;		// log(ptsPerDim)
 *	char *RBorder;		//Red-black order of a d-dim torus with 2 pts per dim
 *	char *PTbits;		// Bits for d coordinates of a point (d x maxlog2(Pts))
 * } meshVars;

 ***********************************************************/

// Turn integer into binary string. log2(n) should be <= arraySize 
void int2bin(unsigned int n, char *array, unsigned int arraySize){
	// Keep dividing n by 2 until it becomes 0
	while(n){
		if(n & 0x01) //i-th bit is 1
			array[--arraySize] = '1';
		else
			array[--arraySize] = '0';
		n >>= 1;
	}
	// The rest
	while(arraySize) array[--arraySize]='0';
}

// Convert binary string to 10-base integer, Assumes array has at least arraySize
unsigned int bin2int(char *array, unsigned int arraySize){
	unsigned int result = 0;
	unsigned int power = 1;
	while(arraySize){
		if(array[--arraySize] == '1') result += power;
		power <<= 1;
	}
	return result;
}

/***********************************************************/
/* RBorder array: Permutation array of 2^d nodes.
 * Permutation nodes are kept in binary character strings
 * so it is a 1D array that can be viewed as:
 * 	  d columns
 * 2^d   <11000110>    -> Represents a # for 0 to 2^d-1 
 * rows  <01001010>
 *       <10111010>
 */
void findRBorder(struct meshVars *mesh){
	unsigned int i, j, reds, blacks, next, destination;
	unsigned int *Cols;

	Cols = (unsigned int *)malloc(sizeof(unsigned int)*mesh->twoTod);
	// Define RBorder
	Cols[0] = 0;
	next = 1;
	for(i = 1; i <= mesh->d; ++i){
		for(j = 0; j < next; ++j)
			Cols[next + j] = !Cols[j];
		next <<= 1;
	}

	reds = 0; blacks = mesh->twoTod>>1; // blacks is 2^{d-1}
	for(i = 0; i < mesh->twoTod; ++i){
		if(Cols[i] == 0)
			destinations = ++reds;
		else
			destination = ++blacks;
		int2bin(destination, &(mesh->RBorder[i*mesh->d]), mesh->d);
	}

	free(Cols);
}

/***********************************************************/
void hierOrderSetUp(struct meshVars *mesh){
	unsigned int i, j, k, n0;
	unsigned int *logPts;

	// Find 2^d
	mesh->twoTod = 1; for (i = 0; i < mesh->d; ++i) mesh->twoTod <<= 1;
	// Allocate auxiliary arrays
	// log2(ptsPerDim)
	mesh->logPts = (unsigned int *)malloc(sizeof(unsigned int)*mesh->d);
	// Binary RB order as a char array 2^d \times d
	mesh->RBorder = (char *)malloc(sizeof(char)*mesh->twoTod*mesh->d);

	// Basic parameters
	mesh->dmax = 0;
	for(i = 0; i < mesh->d; ++i){
		if(mesh->ptsPerDim[i] > mesh->ptsPerDim[mesh->dmax])
			mesh->dmax = i;		// Track the maximum sized dimension
		mesh->logPts[i] = 0;
		n0 = mesh->ptsPerDim[i];
		while(n0 >>= 1) ++mesh->logPts[i]; 	// Finds log2(ptsPerDim[i])
	}

	// Auxiliary bit array for all d coordinates of a point (d x maxlog2(Pts))
	mesh->PTbits =(char *)malloc(sizeof(char)*mesh->d*mesh->logPts[mesh->dmax]);

	// Find the red-black ordering in binary for a 2 point d-dimensional torus
	findRBorder(mesh);
}

/***********************************************************/
unsigned int hierOrderPoint(unsigned int *coord, struct meshVars *mesh){
	unsigned int logmax = mesh->logPts[mesh->dmax];
	unsigned int m, j, from, to, nbits, activeDims;
	char *location_str;		// Contains up to logPts(dmax) strings of d bits
	char *crossbits;		// The d j-th bits of all coordinates

	location_str = (char *)malloc(sizeof(char)*mesh->d*logmax);
	crossbits = (char *)malloc(sizeof(char)*mesh->d);
	
	for(j = 0; j < mesh->d; ++j){
		// Bin form of coordinate j into j-th row of PTbits
		int2bin(coord[j], &(mesh->PTbits[j*logmax]), logmax);
	}

	nbits = 0;			// Number of bits in location_str

	// Go over coordinate bits array from the least to most significant bit
	for(m = logmax; m > 0; --m){
		// Choose bits only from dimensions which can be subdivided further
		activeDims = 0;
		for(j = 0; j < mesh->d; ++j) crossbits[j] = '0'; // Zero out crossbits
		from = 0; 		// Starting at row 0 of PTbits
		for(j = 0; j < mesh->d; ++j){
			if(logmax - m < mesh->logPts[j])
				crossbits[++activeDims] = mesh->PTbits[from+m-1];// (j, m)
			from += logmax; // Go to next (j + 1) row in PTbits
		}
		// Take the first activeDims bits of the RBorder(crossbits)
		from = mesh->d * bin2int(crossbits, mesh->d);
		to = from + activeDims;
		for(j = from; j < to; ++j){
			location_str[++nbits] = mesh->RBorder[j];
		}
	}

	j = bin2int(location_str, nbits);
	free(location_str);
	free(crossbits);
	return j;
}

/***********************************************************/
// Deallocate mesh information
void freeMeshVars(struct meshVars *mesh){
	free(mesh->ptsPerDim);
	free(mesh->logPts);
	free(mesh->RBorder);
	free(mesh->PTbits);
}

/***********************************************************
 * Hadamard vector related functions 
 ***********************************************************/

/* hadaColPerm produces the index sequence of hadamard columns that
 * correspond to a red-black ordering of a power of 2 dimension
 *
 * **CAUTION** IT ONLY WORKS FOR N WHICH ARE A POWER OF 2 **
 * 
 * Input: 
 *	unsigned int N: Dimension of the vector space
 *	unsigned int Hpsize: Number of indices needed in the sequence
 * Output:
 * 	unsigned int *Hperm: the array containing the indices
 */
void hadaColPerm(unsigned int N, unsigned int* Hperm, unsigned int Hpsize){
	unsigned int n0, newIndices, step, i, k, icur, logN;
	logN = 0;

	if(N < 2) return;
	
	// logN = ceil(log2(N)), n0 = 2^logN
	n0 = N - 1;
	while(n0 >>= 1) ++logN;
	++logN;
	n0 = 1;
	for(i = 1; i < logN; ++i) n0 <<= 1;

	// Every step doubles the Hperm elements with step (1:2:newinds-1)*n0/newinds
	newindices = 2;
	step = n0/2;
	icur = 0;
	Hperm[++icur] = 0;
	for(k = 1; k <= logN; ++k){
		for(i = 1; i < newindices; i += 2){
			Hperm[++icur] = i * step;
			if(icur == Hpsize) return;
		}
		newindices = newindices <<= 1; // Times 2
		step = step >>= 1;	       // Divided by 2
	}
}

/***********************************************************/
unsigned int num1bits(unsigned int a){
	// Find out how many 1 bit are in the binary representation of a
	unsigned int bits = 0;
	unsigned int bit = 1;
	while(a >= bit){
		if(a & bit) ++bits;
		bit = bits << 1;
	}
	return bits;
}

int Hada_element(unsigned int i, unsigned int j){
	unsigned int bits = num1bits(i & j);
	if(bits % 2) return -1;
	return 1;
}

/***********************************************************
 * PROBLEM SPECIFIC FUNCTIONS
***********************************************************/

/***********************************************************/
void index2coord(unsigned int i, struct meshVars *mesh, unsigned int *coord){
	unsigned int j;
	coord[0] = i % mesh->ptPerDim[0];
	for(j = 1; j < mesh->d; ++j){
		i = (unsigned int)(i - coord[j-1])/mesh->ptPerDim[j-1];
		coord[j] = i % mesh->ptsPerDim[j];
	}
}

/***********************************************************/
// Create permutation array for all i rows local on this processor
// This is problem specific and SHOULD BE ADAPTED
void hierPerm(struct meshVars *mesh, unsigned int *perm, unsigned int N){

}

/***********************************************************/
// Not using for now, this will go in our inline measurements.
#if 0
int main()
{
  struct meshVars mesh;
  unsigned int i, N, Hpsize, sample;
  unsigned int *perm, *Hperm;
  double *RHS;
  
  /* Interactive input */
  printf("Give the mesh dimensions: ");
  scanf("%u", &(mesh.d));
  printf("Give number of points per dimension **powers of 2 only**\n");
  mesh.ptsPerDim = (unsigned int *)malloc(sizeof(unsigned int)*mesh.d);
  N=1;
  for (i=0;i<mesh.d; i++){
    printf("dim(%u): ",i);
    scanf("%u", &mesh.ptsPerDim[i]);
    N *= mesh.ptsPerDim[i];
  }  
  printf("\n");

  /* Set up the hierarchical data structs */
  hierOrderSetUp(&mesh);

  /* Find the row permutation once for each point in local mesh*/
  /* Here it's done for all nodes N. But can be done only for local ones */
  perm = (unsigned int *)malloc(sizeof(unsigned int)*N);
  hierPerm(&mesh, perm, N);
  
  /* Now with the perm obtained we do not need the mesh any more */
  freeMeshVars(&mesh);

  /* Create the column permutation of the Hadamard vectors. Do it once */
  /* Create only for as many columns as you need. For example 1024 or 4096 */
  /* CAUTION: N is the global size (total spatial dimension of the mesh) */
  /* CAUTION: N must be a power of 2 */
  Hpsize = 1024;
  Hperm = (unsigned int *)malloc(sizeof(unsigned int)*Hpsize);
  hadaColPerm(N, Hperm, Hpsize);

  /**********************************************************************/
  /* At this point we are ready to create the permuted Hadamard vectors */
  /* What is needed is perm and Hperm */
  /**********************************************************************/

  /* Example production of right hand sides */
  RHS = (double *)malloc(sizeof(double)*N);

  for (sample=0; sample<Hpsize; sample++) {
     /* I create all N here, but it can be done on local rows only */
     for (i=0;i<N;i++) 
        RHS[i] = (double) Hada_element(perm[i], Hperm[sample]);

     /* use RHS in trace computation */
  }

  /* When it is all done free */
  free(perm);
  free(Hperm);
  free(RHS);
}
#endif

