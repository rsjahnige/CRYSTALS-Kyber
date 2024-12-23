#include "ml_kem.h"
#include <stdio.h>

void printBytes(union byte* B, unsigned int len) {
	for (int i=0; i < len; i++) {
		printf("B::Index %d: %d\n", i, B[i].e);
	}
}

void printIntegers(union integer* F) {
	for (int i=0; i < 256; i++) {
		printf("F::Index %d: %d\n", i, F[i].t);
	}
}

int main() {
	union integer* F;
	union byte* B;
	int failCount = 0;

	int multiple = 4096/256;
	int bitCount = 12;

	F = malloc(sizeof(union integer) * 256);

	for (int i=0; i < 256; i++) {
		F[i].t = i * multiple;
	}

	B = ByteEncode(F, bitCount);
	//printBytes(B, 32*bitCount);

	free(F); 	

	F = ByteDecode(B, bitCount);
	//printIntegers(F);

	for (int i=0; i < 256; i++) {
		if (F[i].t != (i * multiple)) {
			printf("Test Failed: %d\n", i);
			failCount += 1;
		}
	}

	if (failCount == 0) printf("Test Successful!\n");
	else printf("Fail Count: %d\n", failCount);

	free(F);
	free(B);
		
	return 0;
}
