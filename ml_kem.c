/********************************************************
 * File: ml_kem.c
 * Author: Ryan Jahnige
 *
 * Description: Module-Lattice-Based Key-Encapsulation
 * 		Mechanism Standard (FIPS-203)
 * *****************************************************/

#include "ml_kem.h"

// How do you round without floating-point arithmetic?
//#define ROUND(x,y) ((x) < ((y) + (0.5)) ? (x) : ((y) + 1))

// 7-bit reversal (see FIP-203:2.3)
static union byte BitRev7(union byte r) {
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
static union byte* BitsToBytes(union bit* b, unsigned int l) {
	union byte* B;
	B = malloc((sizeof(union byte) * l) / 8);

	for (int i=0; i < l; i++) {
		if ((i % 8) == 0) B[i/8].e = 0b00000000;	// initialize byte array index
		B[i/8].e |= (b[i].o & 0b00000001) << (i % 8);	// B[i/8] = B[i/8] + b[i] * 2^(i%8)
	}

	return B;
}

// Perfroms the inverse of BitToBytes, converting a byte array into a bit array
//
// NOTE: The output bit array is in little-endian order
static union bit* BytesToBits(union byte* B, unsigned int L) {
	union bit* b;
	union byte* C;

	b = malloc(sizeof(union bit) * L * 8);
	C = malloc(sizeof(union byte) * L);

	for (int i=0; i < L; i++) {
		C[i].e = B[i].e;					// Copy B[i] into array C[i] in B^l
		for (int j=0; j < 8; j++) {
			b[(8 * i) + j].o = C[i].e & 0b00000001;		// b[8i+j] = C[i] % 2 
			C[i].e = C[i].e >> 1;				// C[i] = C[i] / 2
		}
	}

	free(C);
	return b;
}

// Encode an array of d-bit integers into a byte array for 1<=d<=12
//
// NOTE: F is in the finite field Zm^256; therefore, the lenght of F 
// 	shall always be 256 (see FIPS-203:2.3)
union byte* ByteEncode(union integer* F, unsigned int d) {
	union integer a;
	union bit* b;
	union byte* B;

	b = malloc(sizeof(union bit) * 256 * d);

	for (int i=0; i < 256; i++) {
		a.t = F[i].t;					// 'a' in Zm
		for (int j=0; j < d; j++) {	
			b[(i * d) + j].o = a.t & 0x001;		// b[i*d+j] = a % 2	
			a.t >>= 1;				// a = (a - b[i*d+j])/2
		}

	}
	
	B = BitsToBytes(b,256*d);
	free(b);

	return B;
}

// Decode a byte array into an array of d-bit integers for 1<=d<=12
//
// NOTE: For 1<=d<=11, the conversion is one-to-one; whereas, for d=12 
// 	it is no longer a one-to-one operation when the input byte array
// 	is NOT produced by ByteEncode()
union integer* ByteDecode(union byte* B, unsigned int d) {
	union bit* b;
	union integer* F;
	unsigned int m;

	b = BytesToBits(B,32*d);
	F = malloc(sizeof(union integer) * 256);

	// Set the value of 'm' appropriately
	if ((1 <= d) && (d < 12)) m = 0x001 << d;
	else if (d == 12) m = Q;
	else exit(EXIT_FAILURE);

	for (int i=0; i < 256; i++) {
		F[i].t = 0x000;		// initialize F[i]

		// F[i] = sum_j=0->(d-1){(b[i*d+j]*2^j) mod m}
		for (int j=0; j < d; j++){
			F[i].t |= (b[(i * d) + j].o * (0x001 << j)) % m;			}
	}

	free(b);

	return F;
}
