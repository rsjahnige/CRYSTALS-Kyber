#include "sha3.h"

// SHA-3 Constants (see FIPS-202:3.4)
#define B	1600		// total length of the state array A
#define W	64		// length of z-axis/lane size (i.e., len(A)=B such that W=B/25)
#define Nr	24		// number of rounds
#define L	6		// binary logarithm of the lane size (i.e., log2(W)=L)

/***************************************
 * Keccak_f Permutation Functions
 * ************************************/

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

// Iota function defined in FIPS-202:3.2.5 - the effect is to modify some of the
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

// Keccak_f permutation defined in FIPs-202:3.3/4 - the Keccak_f[1600] permutation
// underlies the six SHA-3 functions
//
// NOTE: Keccak_f[B] = Keccak_p[B,12+2L] 
static union bit* Keccak_f(union bit* S) {

	// Step 2: Apply the step mappings, in order, for the specified number 
	// 		of rounds (always 24 for SHA-3 family of functions)
	for (int i=(12+2*L-Nr); i < (12+2*L); i++) {
		S = Iota(Chi(Pi(Rho(Theta(S)))), i);
	}

	return S;
}

/*********************************
 * Sponge Construction Functions
 * ******************************/ 

// pad10*1 function defined in FIPS-202:5.1 - produce an output string of the
// form 10*1 such that the 0 bit is either omitted or repeaded as necessary 
//
// NOTE: m+len(pad(x,m)) is a positive multiple of x
static union bit* pad(unsigned int x, unsigned int m) {
	union bit* P;
	int j = x - ((m + 2) % x);	// Step 1: set j accordingly
	
	P = malloc(sizeof(union bit) * (j+2));

	// Step 2: Initialize P to 1||0^j||1
	P[0].b = 1;
	for (int i=1; i <= j; i++) {
		P[i].b = 0;
	}
	P[j+1].b = 1;

	return P;
}

// Sponge[f,pad,r] function defined in FIPS-202:4 - An arbitrary number of bits, N,
// of length m are absorbed into the function, then a fixed number of random output 
// bits, d, are squeezed out
//
// NOTE: f=Keccak_f[1600] => r=1600-c					(1)
// 	Futhermore, m+len(pad(x,m)) = q*x for some q in I
// 	And, j=(-m-2)%x such that len(pad(x,m)) = j+2
// 	=> len(pad(x,m)) <= x OR len(pad(x,m)) = x+1
// 	=> (m % x)+len(pad(x,m)) = q*x for q in {1,2}
// 	Since x=r for Sponge[f,pad,r]
// 	If q=1, then len(pad(x,m)) <= r					(2) 
// 		=> len(pad(x,m) = r - (m % r)
// 	If q=2, then len(pad(x,m)) = r+1				(3)
// 		=> -m-2 = -1 % r => -m % r = 1 % r
// 		=> m % r = r - 1
static union bit* Sponge(union bit* N, unsigned int m, 
			unsigned int d, unsigned int c) {
	union bit* P;					// The concatenation of N || pad(r, m)
	union bit* S;					// Stores pad() return array, then Keccak_f permutation
	union bit* Z;
	unsigned int n, it;				// n = len(P)/2; it = iterator for Z
	unsigned int r = B - c;				// See (1) above - r = rate
	unsigned int x = m % r;

	// Derive the length of pad(r,m)
	if (x == (r-1)) x = r + 1;			// See (3) above
	else x = r - x;					// See (2) above

	// Step 1: let P = N || pad(r, m)
	P = malloc(sizeof(union bit) * (m + x));
	S = pad(r, m);

	for (int i=0; i < m; i++) P[i].b = N[i].b;
	for (int i=m; i < (m+x); i++) P[i].b = S[i-m].b;

	free(S);

	// Steps 2-3: Let n=len(P)/r - 'c' is not used
	n = (m + x) / r;

	// Step 5: Let S = 0^B
	S = malloc(sizeof(union bit) * B);
	for (int i=0; i < B; i++) S[i].b = 0;

	// Step 6: Perform the Keccak_f permutation on S such that 
	// 		S = S[0-(r-1)] XOR P_i for each i in {0,...,(n-1)}
	for (int i=0; i < n; i++) {
		for (int j=0; j < r; j++) S[j].b ^= P[r*i+j].b;
		S = Keccak_f(S);
	}


	// Step 7: Initialize Z - empty string of bits
	Z = malloc(sizeof(union bit) * d);

	it = 0;
	while (it < d) {
		// Step 9: If d <= |Z|, then Trunc_d(Z)
		if ((it + r) >= d) {			
			for (int i=0; i < (d-it); i++)
				Z[it+i].b = S[i].b;
			it += (d - it);
		} else {
			// Step 8: Let Z = Z || Trunc_r(S)
			for (int i=0; i < r; i++)
				Z[it+i].b = S[i].b;
			it += r;
			S = Keccak_f(S);		// Step 10: Let S = f(S)
		}
	}

	free(P); 
	free(S);

	return Z;
}

/**********************************
 * Conversion Functions 
 *********************************/

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

/***************************************
 * SHA-3 Function Specifications
 **************************************/

// SHA-3 function that maps a bit string, bstr, of length n to a
// to a string of random bits of length d
//
// NOTE: The sfx argument is used to distinguish whether this function
// 	behaves as a Hash Function or an Extendable-Output Function (XOF); 
// 	for Hash functions c=2d whereas for XOF c is in {256,512}	
union bit* sha3_b(const union bit* bstr, unsigned int n, unsigned int d, 
			unsigned int c, union bit sfx[4]) {
	union bit* N;
	union bit* D;
	unsigned int s;			// Length of suffix - either 2 or 4

	if (sfx[2].b == 1) s = 4;	// Shake suffix length
	else s = 2;			// Hash and RawShake suffix length

	N = malloc(sizeof(union bit) * (n + s));

	// Set higher order bits of N equal to bstr
	for (int i=0; i < n; i++) N[i].b = bstr[i].b;

	// Append suffix to the end of N - N = bstr || sfx
	if (s == 2) {
		N[n].b = sfx[0].b;
		N[n+1].b = sfx[1].b;
	} else {
		for (int i=0; i < 4; i++)
			N[n+i].b = sfx[i].b;
	}

	// Map N to a random string of bits of length d
	D = Sponge(N, (n+s), d, c);
	free(N);

	return D;
}

// SHA-3 function that maps a hex string, hstr, of length 2*m to a
// to a string of random hex values of length 2*(d/8)
//
// Note: This is primarily an intermediary function the permits the 
// 	mapping between character strings and bit strings
union hex* sha3_h(const union hex* hstr, unsigned int m, unsigned int d, 
			unsigned int c, union bit sfx[4]) {
	union bit* N;
	union bit* M;
	union hex* D;

	N = h2b(hstr, m, 8*m); 			// Convert hex string to a string of bits	
	M = sha3_b(N, 8*m, d, c, sfx);		// Execute the SHA-3 Keccak sponge function
	free(N);

	D = b2h(M, d);				// Convert the output bit string to a hex string
     	free(M);	
	
	return D;
}

// SHA-3 function that maps a character string, cstr, of lenght m to 
// a string of random characters of length d/8
//
// NOTE: The sfx argument is used to distinguish whether this function
// 	behaves as a Hash Function or an Extendable-Output Function (XOF); 
// 	for Hash functions c=2d whereas for XOF c is in {256,512}	
unsigned char* sha3_s(const char* cstr, unsigned int m, unsigned int d, 
			unsigned int c, union bit sfx[4]) {
	union hex* H;
	union hex* Z;
	unsigned char* D;

	// Convert the input character string to a string of
	// hex values of length 2*m
	H = malloc(sizeof(union hex) * m * 2);
	for (int i=0; i < m; i++) {
		H[2*i].d = (cstr[i] >> 4) & 0x0F;
		H[2*i+1].d = cstr[i] & 0x0F;
	}

	// Execute the SSHA-3 Keccak sponge function
	Z = sha3_h(H, m, d, c, sfx);
	free(H);

	// Convert the output hex string to a string of
	// characters of length d/8
	D = malloc(sizeof(unsigned char) * d/8);
	for (int i=0; i < d/8; i++) {
		D[i] = (unsigned char) Z[2*i].d;
		D[i] <<= 4;
		D[i] ^= (unsigned char) Z[2*i+1].d;
	}
	free(Z);

	return D;	
}
