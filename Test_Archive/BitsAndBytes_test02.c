#include "ml_kem.h"
#include <stdlib.h>
#include <stdio.h>

static void printBits(union bit* b, unsigned int l) {
	printf("Bits: ");
	for (int i=0; i < l; i++) {
		printf("%d", b[i].o);
	}
	printf("\n");
}

int main() {
	union bit* bitArray;
	union byte* byteArray;

	// Initialize bit array
	bitArray = malloc(sizeof(union bit) * 8);
	for (int i=0; i < 8; i++) {
		bitArray[i].o = i % 2;
	}

	// Print bitArray
	printBits(bitArray, 8);

	// Convert to bytes and print result
	byteArray = BitsToBytes(bitArray, 8);
	printf("Bytes: %d\n", byteArray[0].e);
	
	// Convert back to bits and print result
	bitArray = BytesToBits(byteArray, 1);
	printBits(bitArray, 8);

	free(byteArray);
	free(bitArray);
	return 0;
}
