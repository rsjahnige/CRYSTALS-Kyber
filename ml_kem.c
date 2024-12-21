/********************************************************
 * File: ml-kem.c
 * Author: Ryan Jahnige
 *
 * Description: Module-Lattice-Based Key-Encapsulation
 * 		Mechanism Standard (FIPS-203)
 * *****************************************************/

#include "ml_kem.h"

// How do you round without floating-point arithmetic?
//#define ROUND(x,y) ((x) < ((y) + (0.5)) ? (x) : ((y) + 1))

// Inline 7-bit reversal (see FIP-203:2.3)
inline union byte BitRev7(union byte r) {
	union byte tmp;

	// Swap bit 6 with bit 4 and bit 0 with bit 2
	tmp.s = (0b0101010 & r.s);	// preserve bits 1,3,5
	tmp.s |= (0b0010001 & (r.s >> 2)) | (0b1000100 & (r.s << 2));

	// Swap bits 0-2 with bits 5-7
	r.s = (0b0001000 & tmp.s);	// preserve bit 3
	r.s |= (0b0000111 & (tmp.s >> 4)) | (0b1110000 & (tmp.s << 4));

	return r;
}

/****************************************
 * Conversion and Compression Algorithms
 * *************************************/

// Converts a bit array (of length that is a multiple of 8) into an array of bytes
// 
// NOTE: The input bit array is in little-endian order
union byte* BitsToBytes(union bit* b, unsigned int l) {
	union byte* B;
	B = malloc(l/8);

	for (int i=0; i < l; i++) {
		if ((i % 8) == 0) B[i/8].e = 0b00000000;	// initialize byte array index
		B[i/8].e |= (b[i].o & 0b00000001) << (i % 8);
	}

	return B;
}

// Perfroms the inverse of BitToBytes, converting a byte array into a bit array
//
// NOTE: The output bit array is in little-endian order
union bit* BytesToBits(union byte* B, unsigned int L) {
	union bit* b;
	union byte* C;

	b = malloc(sizeof(union bit) * L * 8);
	C = malloc(L);

	for (int i=0; i < L; i++) {
		C[i].e = B[i].e;	// Copy B[i] into array C[i] in B^l
		for (int j=0; j < 8; j++) {
			b[(8 * i) + j].o = C[i].e & 0b00000001;
			C[i].e = C[i].e >> 1;
		}
	}

	free(C);
	return b;
}
