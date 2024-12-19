#include "ml-kem.h"
#include "stdio.h"

void printBinary(union byte val) {
	union byte temp;

	printf("0b");
	for (int i=6; i >= 0; i--) {
		temp.s = (val.s >> i) & 0b0000001;
		printf("%d", temp.s);
	}
	printf("\n");
}

int main() {
	union byte val;

	for (int i = 0; i <= 127; i++) {
		val.s = i;
		printf("Initial value: ");
		printBinary(val);	
		printf("Reversed value: ");
		val = BitRev7(val);
		printBinary(val);
		printf("\n");
	}

	return 0;
}
