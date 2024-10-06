#include <stdio.h>


void
get_delta_(void *a1, void *a2, int *i)
{
	long	b1 = (long)a1;
	long	b2 = (long)a2;
	long	delta;

	delta = b1-b2;
/*
	printf("%x-%x = %d-%d = %d = %x\n",b1,b2,b1,b2,delta,delta); fflush(stdout);
	printf("%x-%x = %d-%d = %d = %x\n",b1/8,b2/8,b1/8,b2/8,delta/8,delta/8); fflush(stdout);
*/
	*i = (long)delta;
}

void
print_address_(void *array)
{
	unsigned long	x = ((unsigned long) array);
	unsigned long 	y = ((unsigned long) array);
	unsigned long 	z = ((unsigned long) array);

	y=(y>>6)%((512*1024/8)>>6);
	z=(z>>6)%((8*1024/4)>>6);

	printf("0x%014lx  [%012lx] Cache(%03lx|%02lx)\n",
			x,(x>>3),
			y,z
			);
	fflush(stdout);
}

