#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */

#define SITERAND	/* Use site-based random number generators */

/* defines for 3rd nearest neighbor (NAIK) stuff */
#define X3UP 8
#define Y3UP 9
#define Z3UP 10
#define T3UP 11
#define T3DOWN 12
#define Z3DOWN 13
#define Y3DOWN 14
#define X3DOWN 15
 
#define OPP_3_DIR(dir) (23-(dir))
#define DIR3(dir) ((dir)+8)
#define FORALL3UPDIR(dir) for(dir=X3UP; dir<=T3UP; dir++)

#define INVERT_M 0
#define INVERT_MdaggerM 1

#endif /* _DEFINES_H */
