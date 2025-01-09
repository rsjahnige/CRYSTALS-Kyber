#include "sha3.h"

// Conversion function from hexadecimal strings to the SHA-3 bit
// strings they represent (FIPS-202:A-B.1)
//
// NOTE: len(H) = 2m & len(S) = n <= 8m = len(T)
// WARNING: This function makes the assumption that (unsigned int)
// 		is a minimum of 8 bits
union bit* h2b(const union hex* H, unsigned int m, unsigned int n) {
	union bit* T;
	union bit* S;
	unsigned int h;

	T = malloc(sizeof(union bit) * 8 * m);

	// Steps 2-3: build T bit string from the hexadecimal string,
	// 		H, provided as input
	for (int i=0; i < m; i++) {
		h = 16 * H[2*i].d + H[2*i+1].d;
		for (int j=0; j < 8; j++) {
			T[8*i+j].b = h & 1;
			h = h >> 1;
		}	
	}

	// Step 4: build S bit string by truncating the T bit string 
	// 		to the size specified by n
	if (n < (8*m)) { 
		S = malloc(sizeof(union bit) * n);
		for (int i=0; i < n; i++) {
			S[i].b = T[i].b;
		}
		free(T);
	} else {
		S = T;
	}
	
	return S;
}	

// Conversion function from SHA-3 bit strings to the hexadecimal
// strings they represent (FIPS-202:A-B.1)
//
// NOTE: len(S) = n & len(H) = 2 * ceiling(n/8) 
// WARNING: This function makes the assumption that (unsigned int)
// 		is a minimum of 8 bits
union hex* b2h(const union bit* S, unsigned int n) {
	union bit* T;
	union hex* H;
	unsigned int h = n % 8;
	unsigned int m = (n + 8 - 1) / 8;	// ceiling calculation for n/8

	h = (h != 0) ? (8 - h) : h;
	T = malloc(sizeof(union bit) * (n + h));
	H = malloc(sizeof(union hex) * 2 * m);

	// Step 2: Copy SHA-3 bit string S into the padded SHA-3
	// 		bit string T; 'm' is set above
	for (int i=0; i < (n+h); i++) {
		T[i].b = (i < n) ? S[i].b : 0;
	}

	// Steps 3-4: Convert each index of the bit string T to its hexadecimal 
	// 		equivalent, H
	for (int i=0; i < m; i++) {
		h = 0;
		for (int j=0; j < 8; j++) {
			h += (1 << j) * T[8*i+j].b;	// (T[0] * 2^0) + ... + (T[7] * 2^7) 
		}
		// Swap higher 4 bits with lower 4 bits
		H[2*i].d = (h >> 4) & 15;
		H[2*i+1].d = h & 15;
	}

	return H;
}

// Theta function defined in FIPS-202:3.2.1 - the effect is to XOR each
// bit in the state array with the parities of two columns
static union bit* Theta(union bit* A) {
	union bit C[5][W];
	union bit D[5][W];

	// Step 1: build C array by taking the XOR sum of all
	// 		bits in each column
	for (int x=0; x < 5; x++) {
		for (int z=0; z < W; z++) {
			C[x][z].b = 0;
			for (int y=0; y < 5; y++) {
				C[x][z].b ^= A[W*(5*y+x)+z].b;
			}
		}
	}

	// Step 2: build D array by taking the XOR sum of 2 columns 
	// 		adjacent to each bit in the state array
	for (int x=0; x < 5; x++) {
		for (int z=0; z < W; z++) {
			D[x][z].b = C[(x+4)%5][z].b ^ C[(x+1)%5][(z+W-1)%W].b;
		}
	}

	// Step 3: Apply the theta step mapping to each bit of the 
	// 		state array, A
	for (int x=0; x < 5; x++) {
		for (int y=0; y < 5; y++) {
			for (int z=0; z < W; z++) {
				A[W*(5*y+x)+z].b ^= D[x][z].b;
			}
		}
	}
 
	return A;		// The state array, A, is modified in place
}

// Rho function defined in FIPS-202:3.2.2 - the effect is to rotate each
// bit in the lane by an offset dependent upon the fixed x and y coordinates
static union bit* Rho(union bit* A) {
	union bit Ap[5][5][W];
	unsigned int x, y, yt, offset;

	// Step 1 (modified): State array A is copied into state array
	// 			Ap so that A can be modified in place
	for (x=0; x < 5; x++) {
		for (y=0; y < 5; y++) {
			for (int z=0; z < W; z++) {
				Ap[x][y][z].b = A[W*(5*y+x)+z].b;
			}
		}
	}

	// Step 2: Set (x,y) = (1,0)
	x=1;
	y=0;

	// Step 3: Rotate the bits of each lane by the offset
	for (int t=0; t < 24; t++) {
		offset = (t+1)*(t+2)/2;
		for (int z=0; z < W; z++) {
			A[W*(5*y+x)+z].b = Ap[x][y][(z-offset)%W].b;
		}
		// Update (x,y) coordinates in accordance with step 3b
		yt = y;
		y = (2*x + 3*y) % 5;
		x = yt;
	}

	return A;		// The state array, A, has been updated appropriately 
}

// Pi function defined in FIPS-202:3.2.3 - the effect is to rearrange the
// positions of the lanes 
static union bit* Pi(union bit* A) {
	union bit Ap[5][5][W];

	// Copy state array A into state array Ap so that A
	// can be modified in place
	for (int x=0; x < 5; x++) {
		for (int y=0; y < 5; y++) {
			for (int z=0; z < W; z++) {
				Ap[x][y][z].b = A[W*(5*y+x)+z].b;
			}
		}
	}

	// Step 1: Rearrange the bits in each slice
	for (int x=0; x < 5; x++) {
		for (int y=0; y < 5; y++) {
			for (int z=0; z < W; z++) {
				A[W*(5*y+x)+z].b = Ap[(x+3*y)%5][x][z].b;
			}
		}
	}


	return A;		// The state array, A, has been updated appropriately
}

// Chi function defined int FIPS-202:3.2.4 - the effect is to XOR each bit
// with a non-linear function of two other bits in its row
static union bit* Chi(union bit* A) {
	union bit Ap[5][5][W];

	// Copy state array A into state array Ap so that A
	// can be modified in place
	for (int x=0; x < 5; x++) {
		for (int y=0; y < 5; y++) {
			for (int z=0; z < W; z++) {
				Ap[x][y][z].b = A[W*(5*y+x)+z].b;
			}
		}
	}
	
	// Step 1: XOR each bit with a function of two other bits in its row
	for (int x=0; x < 5; x++) {
		for (int y=0; y < 5; y++) {
			for (int z=0; z < W; z++) {
				A[W*(5*y+x)+z].b = 
					Ap[x][y][z].b ^ ((Ap[(x+1)%5][y][z].b ^ 1) * Ap[(x+2)%5][y][z].b);
			}
		}
	}

	return A;
}

// Round Constant function defined in FIPS-202:3.2.5 - Generate the round
// constant using a function that is based on a linear feedback shift
// register 
// 
// WARNING: This function makes the assumption that (unsigned int)
// 		is a minimum of 8 bits
static union bit rc(unsigned int t) {
	union bit R[9];

	// Step 1: Return 1
	if ((t % 255) == 0) {
		R[0].b = 1;
		return R[0];
	}

	// Step 2: set R = 0100000000
	for (int i=0; i < 9; i++) {
		R[i].b = 0;
	}
	R[1].b = 1;	

	// Step 3: Determine the round constant 
	for (int i=1; i <= (t % 255); i++) {
		R[0].b = 0;			// Equivalent to 0 || R
		R[0].b ^= R[8].b;
		R[4].b ^= R[8].b;
		R[5].b ^= R[8].b;
		R[6].b ^= R[8].b;
		
		// Equivalent to R=Trunc8[R] - R[0] is reset to 0 above
		for (int j=8; j > 0; j--) {
			R[j].b = R[j-1].b;
		}
	}

	return R[1];
}

// Iota function defined in FIP-202:3.2.5 - the effect is to modify some of the
// bits of lane (0,0) in a manner that depends on the round index ir
static union bit* Iota(union bit* A, unsigned int ir) {
	union bit RC[W];
	
	// Step 2: Initialize round constant to 0^W
	for (int i=0; i < W; i++) {
		RC[i].b = 0;
	}

	// Step 3: Populate RC array using the rc() function
	for (int j=0; j <= L; j++) {
		RC[(1<<j)-1] = rc(j+7*ir);
	}

	// Step 4: Modify lane (0,0) using the round constant RC
	for (int z=0; z < W; z++) {
		A[z].b ^= RC[z].b;
	}

	return A;
}

union bit* Keccak_f(union bit S[B]) {

	for (int i=(12+2*L-Nr); i < (12+2*L); i++) {
		S = Iota(Chi(Pi(Rho(Theta(S)))), i);
	}

	return S;
}
