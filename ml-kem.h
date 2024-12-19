/****************************************************
 * File: ml-kem.h
 * Author: Ryan Jahnige
 *
 * Description: Header file for the ML-KEM algorithm.
 * 		There is NO "passing by reference"
 * **************************************************/

// Global integer constants - defined in FIPS-203:2.4
#define N 256
#define Q 3329

// Global integers - intialized when a parameter set 
// is selected (see FIP-203:8)
int k, n1, n2, du, dv;

// "byte" is a loosely defined term, so I'm making it 7 bits
// for this program
union byte {
	int s : 7;	// seven bit integer (FIPS-203:2.3)
};

union byte BitRev7(union byte r);
