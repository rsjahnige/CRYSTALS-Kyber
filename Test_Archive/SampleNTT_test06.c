#include "ml_kem.h"
#include <stdio.h>

int main() {
	union integer* a;
	union byte B[34];
	unsigned int it = 0;

	while (it < 7) {
		// Initialize B
		for (int i=0; i < 34; i++) {
			B[i].e = it*i + i;
		}
		
		a = SampleNTT(B);

		// Print a
		for (int i=0; i < 256; i++) {
			if (a[i].t != 0) 
				printf("%dx^%d + ", a[i].t, i);
		}
		
		free(a);
		it += 1;	
	}



	return 1;
}
