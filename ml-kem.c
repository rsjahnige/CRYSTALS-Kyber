/********************************************************
 * File: ml-kem.c
 * Author: Ryan Jahnige
 *
 * Description: Module-Lattice-Based Key-Encapsulation
 * 		Mechanism Standard (FIPS-203)
 * *****************************************************/

#include "ml-kem.h"

// Inline addition and multiplication functions defined for
// integers in the set Zm (see FIP-203:2.4.1)
inline int add(int x, int y, int m) {
	return (x + y) % m;
}

inline int mult(int x, int y, int m) {
	return (x * y) % m;
}

// Inline 7-bit reversal
union byte BitRev7(union byte r) {
	union byte tmp;

	tmp.s = (0b0101010 & r.s);
	tmp.s |= (0b0010001 & (r.s >> 2)) | (0b1000100 & (r.s << 2));

	r.s = (0b0001000 & tmp.s);
	r.s |= (0b0000111 & (tmp.s >> 4)) | (0b1110000 & (tmp.s << 4));

	return r;
}
