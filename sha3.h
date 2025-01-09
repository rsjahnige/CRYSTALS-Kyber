/*************************************************
 * File: sha3.h
 * Author: Ryan Jahnige
 *
 * Description: Header file for the SHA-3 family of
 * 		algorithms.
 * ***********************************************/

#ifndef SHA3_H
#define SHA3_H

#include <stdlib.h>		// Dynamic memory allocation

// SHA-3 Constants (see FIPS-202:3.4)
#define B	1600		// total length of the state array A
#define W	64		// length of z-axis/lane size (i.e., len(A)=B such that W=B/25)
#define Nr	24		// number of rounds
#define L	6		// binary logarithm of the lane size (i.e., log2(W)=L)

// Data type used to create bit arrays
union bit {
	unsigned int b : 1;	// binary digit in {0,1}	
};

union hex {
	unsigned int d : 4; 	// hexidecimal digit
};


union bit* h2b(const union hex* H, unsigned int m, unsigned int n);
union hex* b2h(const union bit* S, unsigned int n);

union bit* Keccak_f(union bit S[B]);

#endif
