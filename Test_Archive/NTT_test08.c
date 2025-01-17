#include "ml_kem.h"
#include <stdio.h>

int main() {
	union byte B[34];
	union integer* f1;
	union integer* fh;
	union integer* f2;

	// Initialize byte array B
	for (int i=0; i < 34; i++) {
		B[i].e = i * 2;
	}

	f1 = SampleNTT(B);
	fh = NTT(f1);
	f2 = InverseNTT(fh);

	for (int i=0; i < 256; i++) {
		if (f1[i].t != f2[i].t) 
			printf("ERROR :: f1[%d] = %d :: f2[%d] = %d\n", i, f1[i].t, i, f2[i].t);
	}

	free(f1);
	free(fh);
	free(f2);

	printf("Test Complete!\n");

	return 1;
}
