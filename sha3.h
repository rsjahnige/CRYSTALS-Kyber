/*************************************************
 * File: sha3.h
 * Author: Ryan Jahnige
 *
 * Description: Header file for the SHA-3 family of
 * 		algorithms
 * ***********************************************/

#ifndef SHA3_H
#define SHA3_H

#include <stdlib.h>		// Dynamic memory allocation & exit()

// SHA-3 Constants (see FIPS-202:3.4)
#define B	1600		// total length of the state array A
#define W	64		// length of z-axis/lane size (i.e., len(A)=B such that W=B/25)
#define Nr	24		// number of rounds
#define L	6		// binary logarithm of the lane size (i.e., log2(W)=L)

// Data type used to create bit arrays
union bit {
	unsigned int b : 1;	// binary digit in {0,1}	
};

// Data type used to create hex arrays
union hex {
	unsigned int d : 4; 	// hexidecimal digit in {0-F}
};

// Conversion function from hexadecimal strings to the SHA-3 bit
// strings they represent (FIPS-202:A-B.1)
//
// NOTE: len(H) = 2m & len(S) = n <= 8m = len(T)
union bit* h2b(const union hex* H, unsigned int m, unsigned int n);

// Conversion function from SHA-3 bit strings to the hexadecimal
// strings they represent (FIPS-202:A-B.1)
//
// NOTE: len(S) = n & len(H) = 2 * ceiling(n/8) 
union hex* b2h(const union bit* S, unsigned int n);

// SHA-3 function that maps a bit string, bstr, of length n to a
// to a string of random bits of length d
//
// NOTE: The sfx argument is used to distinguish whether this function
// 	behaves as a Hash Function or an Extendable-Output Function (XOF); 
// 	for Hash functions c=2d and sfx={0,1,0,0} whereas for XOF c is 
// 	in {256,512} and sfx={1,1,1,1}
union bit* sha3_b(const union bit* bstr, unsigned int n, unsigned int d, 
			unsigned int c, union bit sfx[4]);

// SHA-3 function that maps a hex string, hstr, of length 2*m to a
// to a string of random hex values of length 2*(d/8)
//
// NOTE: The sfx argument is used to distinguish whether this function
// 	behaves as a Hash Function or an Extendable-Output Function (XOF); 
// 	for Hash functions c=2d and sfx={0,1,0,0} whereas for XOF c is 
// 	in {256,512} and sfx={1,1,1,1}
union hex* sha3_h(const union hex* hstr, unsigned int m, unsigned int d, 
			unsigned int c, union bit sfx[4]);

// SHA-3 function that maps a character string, cstr, of lenght m to 
// a string of random characters of length d/8
//
// NOTE: The sfx argument is used to distinguish whether this function
// 	behaves as a Hash Function or an Extendable-Output Function (XOF); 
// 	for Hash functions c=2d and sfx={0,1,0,0} whereas for XOF c is 
// 	in {256,512} and sfx={1,1,1,1}
unsigned char* sha3_s(const char* cstr, unsigned int m, unsigned int d, 
			unsigned int c, union bit sfx[4]);

#endif
