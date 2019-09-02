/* Heather Switzer
 *
 * This code has been adopted for use in MILC_QCD from Arjun Gambhir's 
 * LaLiBe code found https://github.com/callat-qcd/lalibe/blob/cdde01a04b7418c1a6da0ccfcfcaa268af7a5545/lib/measurements/HP_prop_w.cc.
 *
 * This computes a Feyman-Hellmann fully diluted color Hierarchical Probing
 * (HP) propagator with an element-wise product with random noise
 *
 * INPUT:
 *      Type of noise for element-wise Hadamard product
 *      Random Noise Seed
 *      Start/End vector for HP
 * OUTPUT:
 * HP propagator
 *
 */

#include "generic_ks_includes.h"

void create_HP(unsigned int starting_vector, unsigned int ending_vector){

   /* Create the HP stuff, adopted from Andreas' main function */
   struct meshVars mesh;
   unsigned int i, N, Hpsize, sample;
   unsigned int *perm, *Hperm;

   mesh.d = 4; // Number of dimensions - need to change
   mesh.ptsPerDim = (unsigned int *)malloc(sizeof(unsigned int)*mesh.d);
   N = 1;
   
   unsigned int lattSize[4] = {nx, ny, nz, nt};
   for(i = 0; i < mesh.d; ++i){
        // Compute the next highest power of 2 of 32-bit v
        unsigned int v = lattSize[i];
        v--;
        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;
        v++;
        // mesh.ptsPerDim[i] = v;
        // Extend HP to beyond purely powers of 2
        mesh.ptsPerDim[i] = v;
        N *= mesh.ptsPerDim[i];
   }

   // Set up hierarchical data structs
   hierOrderSetUp(&mesh);

   // Find the row permuation once for each point in the local mesh
   // Here it's done for all nodes N. But can be done only for local ones
   perm = (unsigned int *)malloc(sizeof(unsigned int)*N);
   hierPerm(&mesh, perm, N);

   /* With the perm we do not need the mesh anymore */
   freeMeshVars(&mesh);

   /* Create the column permutation of the Hadamard vectors. Do it once.
   *  Create only for as many columns as you need.
   *  CAUTION: N is the global size (total spatial dimension of the mesh)
   *  CAUTION: N must be a power of 2 
   */
   Hpsize = ending_vector;
   Hperm = (unsigned int *)malloc(sizeof(unsigned int)*Hpsize);
   hadaColPerm(N, Hperm, Hpsize);
   
   /*********************************************************************
   * Now we need to create permuted Hadamard vectors.
   * We need perm and Hperm. 
   *********************************************************************/
   
   // *******************************************************************
   // TODO: Right now it is following along with Arjun's LaLiBe code
   // that creates all the vectors and inverts them at one time. This
   // needs to be changed to only store and comoute one HP vector at a
   // time.
   // *******************************************************************

   // Hold HP vectors we want to invert
   su3_vector* vectors; // Creates a vector of vectors of su3_vectors

   int HP_index, x, y, z, t;
   for(HP_index = starting_vector; HP_index <= ending_vector; HP_index++)
      for(x = 0; x < nx; x++)
         for(y = 0; y < ny; y++)
            for(z = 0; z < nz; z++)
               for(t = 0; t < nt; t++){
                  int chroma_coords[4] = {x, y, z, t};
               }

}
