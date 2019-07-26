#ifndef __binaryRecursiveColoring_h_
#define __binaryRecursiveColoring_h_
#include <stdio.h>
#include <stdlib.h> /* mallocs, frees */
//#include <chroma.h> //XXX: Find milc_qcd equivalent

typedef struct mashVars{
	unsigned int d;			// Number of dimensions
	unsigned int *ptsPerDim;	// Points per dimension
	unsigned int dmax;		// Which dimension has most points
	unsigned int twoTod;		// 2^d
	unsigned int *logPts;		// log(ptsPerDim)
	char *RBorder;		//Red-black order of a d-dim torus with 2 pts per dim
	char *PTbits;		// Bits for d coordinates of a point (d x maxlog2(Pts))
} meshVars;

// Turns an integer into a binary string. log(n) should be <= arraySize
void int2bin(unsigned int n, char *array, unsigned int arraySize);

// Converts binary string into an integer. assumes array has at least arraySize
unsigned int bin2int(char *array, unsigned int arraySize);

/* RBorder array: Permutation array of 2^d nodes.
 * Permutation nodes are kept in binary character strings
 * so it is a 1D array that can be viewed as:
 * 	  d columns
 * 2^d   <11000110>    -> Represents a # for 0 to 2^d-1 
 * rows  <01001010>
 *       <10111010>
 */

void findRBorder(struct meshVars *mesh);

void hierOrderSetUp(struct meshVars *mesh);

unsigned int hierOrderPoint(unsigned int *coord, struct meshVars *mesh);

void freeMeshVars(struct meshVars *mesh);


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
void hadaColPerm(unsigned int N, unsigned int* Hperm, unsigned Hpsize);

unsigned int num1bits(unsigned int a);

int Hada_element(unsigned int i, undigned int j);

void index2coord(unsigned int i, struct meshVars *mesh, unsigned int *coord);

/****************************************************************************
 * Create the permutation array for all i rows local on this processor.
 * This is problem specific and SHOULD BE ADAPTED.
 */
void hierPerm(struct meshVars *mesh, unsigned int *perm, unsigned int N);

#endif


