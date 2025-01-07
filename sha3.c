#include "sha3.h"

// Iterator offset for x & y coordinates (see FIPS-202:3.1.4)
// -> Should only be needed when dealing with 'lane' operations
#define IT(i) ((i + 2) % 5)	

// Conversion function from hexadecimal strings to the SHA-3 bit
// strings they represent (FIPS-202:A-B.1)
//
// NOTE: len(H) = 2m & len(S) = n <= 8m = len(T)
// WARNING: This function makes the assumption that (unsigned int)
// 		is a minimum of 8 bits
union bit* h2b(union hex* H, unsigned int m, unsigned int n) {
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
union hex* b2h(union bit* S, unsigned int n) {
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
//
// NOTE: w = length of z-axis (i.e., len(A)=b such that w=b/25)
void Theta(union bit* A, unsigned int w) {
	union bit C[5][w];
	union bit D[5][w];

	// Step 1: build C array by taking the XOR sum of all
	// 		bits in each column
	for (int x=0; x < 5; x++) {
		for (int z=0; z < w; z++) {
			C[x][z].b = 0;
			for (int y=0; y < 5; y++) {
				C[x][z].b ^= A[w*(5*y+x)+z].b;
			}
		}
	}

	// Step 2: build D array by taking the XOR sum of 2 columns 
	// 		adjacent to each bit in the state array
	for (int x=0; x < 5; x++) {
		for (int z=0; z < w; z++) {
			D[x][z].b = C[(x+4)%5][z].b ^ C[(x+1)%5][(z+w-1)%w].b;
		}
	}

	// Step 3: Apply the theta step mapping to each bit of the 
	// 		state array, A
	for (int x=0; x < 5; x++) {
		for (int y=0; y < 5; y++) {
			for (int z=0; z < w; z++) {
				A[w*(5*y+x)+z].b ^= D[x][z].b;
			}
		}
	}
 
	return;		// The state array, A, is modified in place
}

// Rho function defined in FIPS-202:3.2.2 - 
//
// NOTE: w = length of z-axis (i.e., len(A)=b such that w=b/25)
void Rho(union bit* A, unsigned int w) {
	union bit Ap[5][5][w];
	unsigned int x, y, yt, offset;

	// Step 1 (modified): State array A is copied into state array
	// 			Ap so that A can be modified in place
	for (x=0; x < 5; x++) {
		for (y=0; y < 5; y++) {
			for (int z=0; z < w; z++) {
				Ap[x][y][z].b = A[w*(5*y+x)+z].b;
			}
		}
	}

	// Step 2: Set (x,y) = (1,0)
	x=1;
	y=0;

	// Step 3: Rotate the bits of each lane by the offset
	for (int t=0; t < 23; t++) {
		offset = (t+1)*(t+2)/2;
		for (int z=0; z < w; z++) {
			A[w*(5*y+x)+z].b = Ap[x][y][(z-offset)%w].b;
		}
		// Update (x,y) coordinates in accordance with step 3b
		yt = y;
		y = (2*x + 3*y) % 5;
		x = yt;
	}

	return;		// The state array, A, has bee updated appropriately 
}
