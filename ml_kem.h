/****************************************************
 * File: ml_kem.h
 * Author: Ryan Jahnige
 *
 * Description: Header file for the ML-KEM algorithm.
 * 		There is NO "passing by reference"
 * **************************************************/
#ifndef ML_KEM_H
#define ML_KEM_H

#include <stdlib.h>	// dynamic memory allocation

// Global integer constants - defined in FIPS-203:2.4
#define N 256
#define Q 3329

// Global integers - intialized when a parameter set 
// is selected (see FIP-203:8)
int k, n1, n2, du, dv;

// Data type used to create bit arrays
union bit {
	unsigned int o : 1;	// binary digit in {0,1}	
};

// "byte" is a loosely defined term, so for this program it will
// take on a couple different forms - be careful! 
// 
// NOTE: Unions take on the size of the largest data type; in this
// case, that would be 8 bits
union byte {
	unsigned int s : 7;	// seven bit integer (FIPS-203:2.3)
	unsigned int e : 8;	// eight bit integer
};

// Data type for integer arrays F in Zm^256, where m=2^d if 1<=d<12 
// and m=Q if d=12 (see FIPS-203:4.2.1)
//
// NOTE: It is likely more efficient if the default integer size
// 	for your system is used; however, this ensures that sufficient
// 	bits are allocated 
union integer {
	unsigned int t : 12;	// twelve bit integer 	
	unsigned int l : 24; 	// twenty-four bit integer - used for real number calculations
};

union integer Compress(union integer x, unsigned int d);
union integer Decompress(union integer y, unsigned int d);

#endif
