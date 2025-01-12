#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sha3.h"

int main(int argc, char* argv[]) {
	
	FILE *in_fp;
	FILE *out_fp;

	unsigned long bin_len;
	unsigned char digit;
	unsigned int dgstLen;
	unsigned int capacity;

	union bit* bin_input;
	union bit* bin_output;
	union hex* hex_output;

	if (argc < 3) {
		printf("ERROR :: Invalid number of arguments provided: %d\n", argc);
		printf("Accepted format: ./test05 [Hash | XOF] [{224,256,384,512} | {128,256}]\n");
		exit(EXIT_FAILURE);
	}

	if (strcmp(argv[1], "XOF") == 0) {
		if (argc != 4) exit(EXIT_FAILURE);
		capacity = 2 * atoi(argv[2]);
		dgstLen = atoi(argv[3]);	
	} else {
		if (argc != 3) exit(EXIT_FAILURE);
		dgstLen = atoi(argv[2]);
	}

	// Open input file
	in_fp = fopen("bin_msg.txt", "r");
	if (in_fp == NULL) {
		printf("ERROR :: could not open bin_msg.txt\n");
		return 1;
	}

	// Determine the length of the input binary message
	fseek(in_fp, 0L, SEEK_END);
	bin_len = ftell(in_fp);
	rewind(in_fp);

	// Read binary digits from input file (bin_msg.txt)
	bin_input = malloc(sizeof(union bit) * bin_len);
	for (int i=0; i < bin_len; i++) {
		fscanf(in_fp, "%c", &digit);
		bin_input[i].b = digit & 0x01;
	}
	fclose(in_fp);

	// Perform desired Hash or XOF depending on user input
	if (strcmp(argv[1], "Hash") == 0) {
		bin_output = sha3_b(bin_input, bin_len, dgstLen, 2*dgstLen, (union bit[]){0,1,0,0});
	} else if (strcmp(argv[1], "XOF") == 0) {
		bin_output = sha3_b(bin_input, bin_len, dgstLen, capacity, (union bit[]){1,1,1,1});
	} else {
		printf("Invalid parameter provide as argv[1] :: Must be 'Hash' or 'XOF'\n");
		exit(EXIT_FAILURE);
	}

	// Convert resulting binary string to a hex array
	hex_output = b2h(bin_output, dgstLen); 

	// Open output file
	out_fp = fopen("act_output.txt", "w");
	if (out_fp == NULL) {
		printf("ERROR :: could not open act_output.txt.txt\n");
		return 1;
	}

	// Write resulting hex array to output file
	for (int i=0; i < dgstLen/4; i++) {
		fprintf(out_fp, "%x", hex_output[i].d);
	}
	fclose(out_fp);
		
	free(bin_input);
	free(bin_output);
	free(hex_output);

	return 0;
}
