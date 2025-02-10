/****************************************************
 * File: ml_kem.h
 * Author: Ryan Jahnige
 *
 * Description: Header file for the ML-KEM algorithm.
 * 		There is NO "passing by reference"
 * **************************************************/
#ifndef ML_KEM_H
#define ML_KEM_H

#include "sha3.h"
#include <stdlib.h>	// dynamic memory allocation
#include <stdio.h>

// Global integer constants - defined in FIPS-203:2.4
#define N 256
#define Q 3329

// Global integers - intialized when a parameter set 
// is selected (see FIP-203:8)
//int k, n1, n2, du, dv;

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
union integer {
	unsigned int t : 12;	// twelve bit integer 	
	unsigned int l : 24; 	// twenty-four bit integer - used for real number calculations
};

struct ML_KEM {
	union byte k;			// dimensions of the matrix A and vectors s, e, y, and e1
	union byte n1;			// specifies the distribution for generating vectors s, e, and y
	union byte n2;			// specifies the distribution for generating vectors e1 and e2
	union byte du;			// parameter for funtion calls in Encrypt() and Decrypt() 
	union byte dv;			// parameter for funtion calls in Encrypt() and Decrypt()
};

struct PKE {
	union byte* ek;			// encryption key - public
	union byte* dk;			// decryption key - private
	unsigned int ek_len;		// encryption key length
	unsigned int dk_len;		// decryption key length
};

struct PKE KeyGen(const struct ML_KEM* params, const union byte* d);
union byte* Encrypt(const struct ML_KEM* params, union byte* ek, const union byte* m, const union byte* r);
union byte* Decrypt(const struct ML_KEM* params, union byte* dk, const union byte* c); 

#endif
