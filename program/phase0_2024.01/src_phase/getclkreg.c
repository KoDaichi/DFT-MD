#include <stdio.h>
#include "getclkreg.h"

void getclkreg_(ticks *t)
{
#ifdef __sparcv9
     __asm__ __volatile__("rd %%tick, %0" : "=r" (*t));
#elif __x86_64__
     unsigned int eax, edx;
     __asm__ volatile("rdtsc" : "=a" (eax), "=d" (edx));
     *t = ((unsigned long)eax) | (((unsigned long)edx) << 32);
#endif
     return;
}
