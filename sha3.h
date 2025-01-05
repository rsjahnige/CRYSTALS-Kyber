/*************************************************
 * File: sha3.h
 * Author: Ryan Jahnige
 *
 * Description: Header file for the SHA-3 family of
 * 		algorithms.
 * ***********************************************/

#ifndef SHA3_H
#define SHA3_H

#include <stdlib.h>

// Data type used to create bit arrays
union bit {
	unsigned int b : 1;	// binary digit in {0,1}	
};

union hex {
	unsigned int d : 4; 	// hexidecimal digit
};

union bit* h2b(union hex* H, unsigned int m, unsigned int n);
union hex* b2h(union bit* S, unsigned int n);

void Theta(union bit* A, unsigned int w);

#endif
