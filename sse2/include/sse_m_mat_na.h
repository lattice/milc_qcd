#define _inline_sse_mult_su3_na(aa,bb,cc) \
{ \
__asm__ __volatile__ ("movupd %0, %%xmm0 \n\t" \
                      "movupd %1, %%xmm1 \n\t" \
                      "movupd %2, %%xmm2 \n\t" \
                      "xorpd %3, %%xmm0 \n\t" \
                      "xorpd %4, %%xmm1 \n\t" \
                      "xorpd %5, %%xmm2" \
                      : \
                      : \
                      "m" ((bb)->e[0][0]), \
                      "m" ((bb)->e[0][1]), \
                      "m" ((bb)->e[0][2]), \
                      "m" (_sse_sgn4), \
                      "m" (_sse_sgn4), \
                      "m" (_sse_sgn4)); \
__asm__ __volatile__ ("movsd %0, %%xmm3 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm4 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "movsd %4, %%xmm5 \n\t" \
                      "unpcklpd %%xmm3, %%xmm3 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm4, %%xmm4 \n\t" \
                      "mulpd %%xmm0, %%xmm3 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "unpcklpd %%xmm5, %%xmm5 \n\t" \
                      "mulpd %%xmm0, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %5, %%xmm6" \
                      : \
                      : \
                      "m" ((aa)->e[0][0].real), \
                      "m" ((aa)->e[0][1].real), \
                      "m" ((aa)->e[1][0].real), \
                      "m" ((aa)->e[1][2].real), \
                      "m" ((aa)->e[2][0].real), \
                      "m" ((aa)->e[2][1].real)); \
__asm__ __volatile__ ("movsd %0, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm3 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm4 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movsd %3, %%xmm6 \n\t" \
                      "movsd %4, %%xmm7 \n\t" \
                      "xorpd %5, %%xmm0" \
                      : \
                      : \
                      "m" ((aa)->e[0][2].real), \
                      "m" ((aa)->e[1][1].real), \
                      "m" ((aa)->e[2][2].real), \
                      "m" ((aa)->e[0][0].imag), \
                      "m" ((aa)->e[1][1].imag), \
                      "m" (_sse_sgn4)); \
__asm__ __volatile__ ("xorpd %0, %%xmm1 \n\t" \
                      "xorpd %1, %%xmm2 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
                      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm6 \n\t" \
                      "mulpd %%xmm1, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %2, %%xmm6 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %4, %%xmm6 \n\t" \
                      "movsd %5, %%xmm7" \
                      : \
                      : \
                      "m" (_sse_sgn4), \
                      "m" (_sse_sgn4), \
                      "m" ((aa)->e[2][2].imag), \
                      "m" ((aa)->e[1][0].imag), \
                      "m" ((aa)->e[0][1].imag), \
                      "m" ((aa)->e[2][0].imag)); \
__asm__ __volatile__ ("unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movsd %0, %%xmm0 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm7 \n\t" \
                      "unpcklpd %%xmm0, %%xmm0 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm0 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm0, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      : \
                      : \
                      "m" ((aa)->e[0][2].imag), \
                      "m" ((aa)->e[2][1].imag), \
                      "m" ((aa)->e[1][2].imag)); \
__asm__ __volatile__ ("movupd %%xmm3, %0 \n\t" \
                      "movupd %%xmm4, %1 \n\t" \
                      "movupd %%xmm5, %2 \n\t" \
                      : \
                      "=m" ((cc)->e[0][0]), \
                      "=m" ((cc)->e[1][0]), \
                      "=m" ((cc)->e[2][0])); \
__asm__ __volatile__ ("movupd %0, %%xmm0 \n\t" \
                      "movupd %1, %%xmm1 \n\t" \
                      "movupd %2, %%xmm2 \n\t" \
                      "xorpd %3, %%xmm0 \n\t" \
                      "xorpd %4, %%xmm1 \n\t" \
                      "xorpd %5, %%xmm2" \
                      : \
                      : \
                      "m" ((bb)->e[1][0]), \
                      "m" ((bb)->e[1][1]), \
                      "m" ((bb)->e[1][2]), \
                      "m" (_sse_sgn4), \
                      "m" (_sse_sgn4), \
                      "m" (_sse_sgn4)); \
__asm__ __volatile__ ("movsd %0, %%xmm3 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm4 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "movsd %4, %%xmm5 \n\t" \
                      "unpcklpd %%xmm3, %%xmm3 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm4, %%xmm4 \n\t" \
                      "mulpd %%xmm0, %%xmm3 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "unpcklpd %%xmm5, %%xmm5 \n\t" \
                      "mulpd %%xmm0, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %5, %%xmm6" \
                      : \
                      : \
                      "m" ((aa)->e[0][0].real), \
                      "m" ((aa)->e[0][1].real), \
                      "m" ((aa)->e[1][0].real), \
                      "m" ((aa)->e[1][2].real), \
                      "m" ((aa)->e[2][0].real), \
                      "m" ((aa)->e[2][1].real)); \
__asm__ __volatile__ ("movsd %0, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm3 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm4 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movsd %3, %%xmm6 \n\t" \
                      "movsd %4, %%xmm7 \n\t" \
                      "xorpd %5, %%xmm0" \
                      : \
                      : \
                      "m" ((aa)->e[0][2].real), \
                      "m" ((aa)->e[1][1].real), \
                      "m" ((aa)->e[2][2].real), \
                      "m" ((aa)->e[0][0].imag), \
                      "m" ((aa)->e[1][1].imag), \
                      "m" (_sse_sgn4)); \
__asm__ __volatile__ ("xorpd %0, %%xmm1 \n\t" \
                      "xorpd %1, %%xmm2 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
                      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm6 \n\t" \
                      "mulpd %%xmm1, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %2, %%xmm6 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %4, %%xmm6 \n\t" \
                      "movsd %5, %%xmm7" \
                      : \
                      : \
                      "m" (_sse_sgn4), \
                      "m" (_sse_sgn4), \
                      "m" ((aa)->e[2][2].imag), \
                      "m" ((aa)->e[1][0].imag), \
                      "m" ((aa)->e[0][1].imag), \
                      "m" ((aa)->e[2][0].imag)); \
__asm__ __volatile__ ("unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movsd %0, %%xmm0 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm7 \n\t" \
                      "unpcklpd %%xmm0, %%xmm0 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm0 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm0, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      : \
                      : \
                      "m" ((aa)->e[0][2].imag), \
                      "m" ((aa)->e[2][1].imag), \
                      "m" ((aa)->e[1][2].imag)); \
__asm__ __volatile__ ("movupd %%xmm3, %0 \n\t" \
                      "movupd %%xmm4, %1 \n\t" \
                      "movupd %%xmm5, %2 \n\t" \
                      : \
                      "=m" ((cc)->e[0][1]), \
                      "=m" ((cc)->e[1][1]), \
                      "=m" ((cc)->e[2][1])); \
__asm__ __volatile__ ("movupd %0, %%xmm0 \n\t" \
                      "movupd %1, %%xmm1 \n\t" \
                      "movupd %2, %%xmm2 \n\t" \
                      "xorpd %3, %%xmm0 \n\t" \
                      "xorpd %4, %%xmm1 \n\t" \
                      "xorpd %5, %%xmm2" \
                      : \
                      : \
                      "m" ((bb)->e[2][0]), \
                      "m" ((bb)->e[2][1]), \
                      "m" ((bb)->e[2][2]), \
                      "m" (_sse_sgn4), \
                      "m" (_sse_sgn4), \
                      "m" (_sse_sgn4)); \
__asm__ __volatile__ ("movsd %0, %%xmm3 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm4 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "movsd %4, %%xmm5 \n\t" \
                      "unpcklpd %%xmm3, %%xmm3 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm4, %%xmm4 \n\t" \
                      "mulpd %%xmm0, %%xmm3 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "unpcklpd %%xmm5, %%xmm5 \n\t" \
                      "mulpd %%xmm0, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %5, %%xmm6" \
                      : \
                      : \
                      "m" ((aa)->e[0][0].real), \
                      "m" ((aa)->e[0][1].real), \
                      "m" ((aa)->e[1][0].real), \
                      "m" ((aa)->e[1][2].real), \
                      "m" ((aa)->e[2][0].real), \
                      "m" ((aa)->e[2][1].real)); \
__asm__ __volatile__ ("movsd %0, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm3 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm4 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movsd %3, %%xmm6 \n\t" \
                      "movsd %4, %%xmm7 \n\t" \
                      "xorpd %5, %%xmm0" \
                      : \
                      : \
                      "m" ((aa)->e[0][2].real), \
                      "m" ((aa)->e[1][1].real), \
                      "m" ((aa)->e[2][2].real), \
                      "m" ((aa)->e[0][0].imag), \
                      "m" ((aa)->e[1][1].imag), \
                      "m" (_sse_sgn4)); \
__asm__ __volatile__ ("xorpd %0, %%xmm1 \n\t" \
                      "xorpd %1, %%xmm2 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
                      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm6 \n\t" \
                      "mulpd %%xmm1, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %2, %%xmm6 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %4, %%xmm6 \n\t" \
                      "movsd %5, %%xmm7" \
                      : \
                      : \
                      "m" (_sse_sgn4), \
                      "m" (_sse_sgn4), \
                      "m" ((aa)->e[2][2].imag), \
                      "m" ((aa)->e[1][0].imag), \
                      "m" ((aa)->e[0][1].imag), \
                      "m" ((aa)->e[2][0].imag)); \
__asm__ __volatile__ ("unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movsd %0, %%xmm0 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm7 \n\t" \
                      "unpcklpd %%xmm0, %%xmm0 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm0 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm0, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      : \
                      : \
                      "m" ((aa)->e[0][2].imag), \
                      "m" ((aa)->e[2][1].imag), \
                      "m" ((aa)->e[1][2].imag)); \
__asm__ __volatile__ ("movupd %%xmm3, %0 \n\t" \
                      "movupd %%xmm4, %1 \n\t" \
                      "movupd %%xmm5, %2 \n\t" \
                      : \
                      "=m" ((cc)->e[0][2]), \
                      "=m" ((cc)->e[1][2]), \
                      "=m" ((cc)->e[2][2])); \
}
