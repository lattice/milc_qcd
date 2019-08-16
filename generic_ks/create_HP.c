/* Code adapted from https://github.com/callat-qcd/lalibe/blob/master/lib/measurements/HP_fh_prop_w.cc
 *
 * LaLiBe/lib/measurements/HP_fh_prop.cc by Arjun Singh Gambhir
 *
 * Code adapted by Heather Maria Switzer
 *
 */

#include "generic_ks_includes.h"

#define Nd 4

void initHP(int s){
        //Initialize HP stuff
        struct meshVars mesh;
        unsigned int i, N, v, Hpsize, sample;
        unsigned int *perm, *Hperm;

        mesh.d = Nd;
        mesh.ptsPerDim = (unsigned int *)malloc(sizeof(unsigned int)*mesh.d);
        N = 1;
        unsigned latticeDims = {nx, ny, nz, nt};

        for(i = 0; i < mesh.d; ++i){
                // Logic taken from Bit Twiddling Hacks by Sean Anderson, Sept. 5, 2010
                // Increases lattice to next power of 2 size
                v = latticeDims[i]; // 32-bit v
                v--;
                v |= v >> 1;
                v |= v >> 2;
                v |= v >> 4;
                v |= v >> 8;
                v |= v >> 16;
                v++;

                // Extend HP tp beyond purely powers of 2
                mesh.ptsPerDim[i] = v;
                N *= mesh.ptsPerDim[i];
        }

        // Set up hierarchical data structs
        hierOrderSetUp(&mesh);

        // Find row permutation for each point in local mesh
        // Done for all nodes N here, but can be done for only local ones
        perm = (unsigned int *)malloc(sizeof(unsigned int)*N);
        hierPerm(&mesh, perm, N);

        // Don't need mesh anymore now that we have perm
        freeMeshVars(&mesh);

        /* Create col perm of Hadamard vectors. Only need to do once.
         * Create only as many columns as needed
         * CAUTION: N is global size (total spatial dimension of mesh)
         * CAUTION: N must be a power 2
        */

        Hpsize = 1024; // Temporary
        Hperm = (unsigned int *)malloc(sizeof(unsigned int)*Hpsize);
        hadaColPerm(N, Hperm, Hpsize);
        /*****************************************************************
        * Here we create the permuted Hadamard vectors
        * We need Hperm and perm for this.
        *****************************************************************/
         
       //TODO
       //multi1d <LatticeInteger> vectors;
       //vector.resize(params.hpparam.ending_vector - params.hpparam.starting_vector + 1);
        unsigned int x, y, z, t;

        unsigned int HP_index;
        for(HP_index = 1; HP_index <= s; ++HP_index)
                for(x = 0; x < nx; ++x)
                        for(y = 0; y < ny; ++y)
                                for(z = 0; z < nz; ++z)
                                        for(t = 0; t < nt; ++t){
                                                su3_vector
                                        }

         
}


