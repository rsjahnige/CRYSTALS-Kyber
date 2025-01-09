#include <stdio.h>
#include <stdlib.h>
#include "sha3.h"

int main() {
	
	FILE *fp;
	unsigned int byte;
	unsigned int low = 15;		// low four bits of a byte
	unsigned int it = 0;

	union hex *hex_array;
	union bit *bit_string;

	fp = fopen("hex_input.txt", "r");

	if (fp == NULL) {
		printf("ERROR :: could not open hex_input.txt\n");
		return 1;
	}

	fseek(fp, 0L, SEEK_END);
	hex_array = malloc(sizeof(union hex) * ftell(fp));
	rewind(fp);

	while (fscanf(fp, "%02x", &byte) != EOF) {
		hex_array[it].d = (byte >> 4) & low; 
		hex_array[it+1].d = byte & low;
		it += 2;
	}
	fclose(fp);

	bit_string = h2b(hex_array, it/2, 1600);
	free(hex_array);

	bit_string = Keccak_f(bit_string);
	hex_array = b2h(bit_string, 1600);

	printf("Hex String (After Permutation): \n");
	for (int i=0; i < it; i++) {
		if (((i % 2) == 0) && (i != 0)) printf(" ");
		if (((i % 32) == 0) && (i != 0)) printf("\n");
		printf("%x", hex_array[i].d);
	}
	printf("\n");

	free(bit_string);
	free(hex_array);
	return 0;
}
