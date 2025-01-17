#include "ml_kem.h"
#include <stdio.h>

int main() {
	union integer* a;
	union byte B[64*3];

	// Initialize B
	for (int i=0; i < 64*3; i++) {
		B[i].e = i;
	}
	
	a = SamplePolyCBD(B, 3);

	// Print a
	for (int i=0; i < 256; i++) {
		if (a[i].t != 0) 
			printf("%dx^%d + ", a[i].t, i);
	}
	
	free(a);

	return 1;
}
