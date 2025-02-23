/****************************************************
 * File: ml_kem.h
 * Author: Ryan Jahnige
 *
 * Description: Header file for the ML-KEM algorithm
 * **************************************************/
#ifndef ML_KEM_H
#define ML_KEM_H

#include "sha3.h"
#include <stdlib.h>	// dynamic memory allocation
#include <stdio.h>	// printing of error messages

// Used to read from /dev/urandom
#ifdef __linux__
	#include <sys/stat.h>
	#include <unistd.h>
	#include <fcntl.h>
#endif

// Global integer constants - defined in FIPS-203:2.4
#define N 256
#define Q 3329

// ML-KEM error indicator
extern int ml_errno;

// "byte" is a loosely defined term, so for this program it will
// take on a couple different forms. Users do not need to concern 
// themselves with the 7-bit data type; it is only used for internal 
// computations. 
// 
// NOTE: Unions take on the size of the largest data type; in this
// case, that would be 8 bits
union byte {
	unsigned int s : 7;	// seven bit integer (FIPS-203:2.3)
	unsigned int e : 8;	// eight bit integer
};

// "Global integers" - intialized when a parameter set 
// is selected (see FIP-203:8)
struct PARAMS {
	union byte k;			// dimensions of the matrix A and vectors s, e, y, and e1
	union byte n1;			// specifies the distribution for generating vectors s, e, and y
	union byte n2;			// specifies the distribution for generating vectors e1 and e2
	union byte du, dv;		// determines ciphertext length
};
 
struct PKE {
	union byte* ek;			// encryption key - public 
	union byte* dk;			// decryption key - private (see SP 800-227)
	unsigned int ek_len;		// encryption key length - 384*k + 32
	unsigned int dk_len;		// decryption key length - 768*k + 96
};

struct KEM {
	union byte K[32];		// shared secret key - always 32-bytes
	union byte* c;			// ciphertext
	unsigned int c_len;		// ciphertext length - 32*(du*k + dv)
};

// Generate an encapsulation and decapsulation key pair based on the
// provided parameter set - parameter set shall be initalized using the
// init() function. Returns a PKE struct that contains the key pair
//
// NOTE: The calling application is responsible for deallocation PKE.ek
// 	and PKE.dk
struct PKE KEM_KeyGen(const struct PARAMS* params);

// Generate a shared secret key and ciphertext from the provided encapsulation key. The 
// ciphertext length is determined by the input parameter set - see the KEM struct above.
//
// NOTE: The calling application is responsible for deallocating the ciphertext KEM.c;
// 	however, the shared secret key KEM.k does not need to be deallocated since it
// 	is always of fixed length
struct KEM KEM_Encaps(const struct PARAMS* params, const union byte* ek, unsigned int ek_len);

// Generate a shared secret key of length 32 from the provided ciphertext using the
// provided decapsulation key.
//
// NOTE: The calling application is responsible for deallocating the resulting
// 	shared secret key.
union byte* KEM_Decaps(const struct PARAMS* params, const union byte* dk, unsigned int dk_len,
			const union byte* c, unsigned int c_len);

// Available parameter sets defined in FIPS-203
enum ML_KEM {
	ML_KEM_512 = 512, 
	ML_KEM_768 = 768, 
	ML_KEM_1024 = 1024
};

// Initialize a PARAMS struct based on the chosen parameter set
const struct PARAMS init(enum ML_KEM param_set);

#endif
