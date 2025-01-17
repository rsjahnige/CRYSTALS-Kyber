/********************************************************
 * File: ml_kem.c
 * Author: Ryan Jahnige
 *
 * Description: Module-Lattice-Based Key-Encapsulation
 * 		Mechanism Standard (FIPS-203)
 * *****************************************************/

#include "ml_kem.h"

// 7-bit reversal (see FIP-203:2.3)
static union byte BitRev7(union byte r) {
	union byte tmp;

	// Swap bit 6 with bit 4 and bit 0 with bit 2
	tmp.s = (0b0101010 & r.s);	// preserve bits 1,3,5
	tmp.s |= (0b0010001 & (r.s >> 2)) | (0b1000100 & (r.s << 2));

	// Swap bits 0-2 with bits 5-7
	r.s = (0b0001000 & tmp.s);	// preserve bit 3
	r.s |= (0b0000111 & (tmp.s >> 4)) | (0b1110000 & (tmp.s << 4));

	return r;
}

/****************************************
 * Conversion and Compression Algorithms
 * *************************************/

// Converts a bit array (of length that is a multiple of 8) into an array of bytes
// 
// NOTE: The input bit array is in little-endian order
static union byte* BitsToBytes(union bit* b, unsigned int l) {
	union byte* B;
	B = malloc((sizeof(union byte) * l) / 8);

	for (int i=0; i < l; i++) {
		if ((i % 8) == 0) B[i/8].e = 0b00000000;	// initialize byte array index
		B[i/8].e |= (b[i].b & 0b00000001) << (i % 8);	// B[i/8] = B[i/8] + b[i] * 2^(i%8)
	}

	return B;
}

// Perfroms the inverse of BitToBytes, converting a byte array into a bit array
//
// NOTE: The output bit array is in little-endian order
static union bit* BytesToBits(union byte* B, unsigned int L) {
	union bit* b;
	union byte* C;

	b = malloc(sizeof(union bit) * L * 8);
	C = malloc(sizeof(union byte) * L);

	for (int i=0; i < L; i++) {
		C[i].e = B[i].e;					// Copy B[i] into array C[i] in B^L
		for (int j=0; j < 8; j++) {
			b[(8 * i) + j].b = C[i].e & 0b00000001;		// b[8i+j] = C[i] % 2 
			C[i].e = C[i].e >> 1;				// C[i] = C[i] / 2
		}
	}

	free(C);
	return b;
}

// Compression for integers in the finite field Zq
//
// NOTE: The bit length of the input integer stays the same - the 
// 	value only changes for d < 12
static union integer Compress(union integer x, unsigned int d) {
	union integer quo, rem, div;

	if (d < 12) {
		div.l = (0x001 << d) * x.t;	// calculate dividend
		quo.t = div.l / Q;		// calculate quotient
		rem.t = div.l % Q;		// calculate remainder
		
		x.t = quo.t;
		if (rem.t > (Q / 2)) x.t += 1;	// round-up if remainder is > 0.5
		x.t %= (0x001 << d);
	}

	return x;
}

// Decompression for integers in the finite field Z(2^d) where 1<=d<12
//
// NOTE: Similar to the Compress() funtion, the Decompress() function has no
// 	affect on integers in Zq and bit lengths of the integers do not change,
// 	only the values do to mimic the desired behavior
static union integer Decompress(union integer y, unsigned int d) {
	union integer quo, rem, div, dsr;

	if (d < 12) {
		dsr.t = 0x001 << d;		// set divisor
		div.l = Q * y.t;		// calculate dividend

		quo.t = div.l / dsr.t;		// calculate quotient
		rem.t = div.l % dsr.t;		// calculate remainder

		y.t = quo.t;
		if (rem.t >= (dsr.t / 2)) y.t += 1;	// round-up if remainder is >= 0.5
	}

	return y;
}

// Encode an array of d-bit integers into a byte array for 1<=d<=12
//
// NOTE: F is in the finite field Zm^256; therefore, the lenght of F 
// 	shall always be 256 (see FIPS-203:2.3)
static union byte* ByteEncode(union integer* F, unsigned int d) {
	union integer a;
	union bit* b;
	union byte* B;

	b = malloc(sizeof(union bit) * 256 * d);

	for (int i=0; i < 256; i++) {
		a.t = F[i].t;					// 'a' in Zm
		for (int j=0; j < d; j++) {	
			b[(i * d) + j].b = a.t & 0x001;		// b[i*d+j] = a % 2	
			a.t >>= 1;				// a = (a - b[i*d+j])/2
		}

	}
	
	B = BitsToBytes(b,256*d);
	free(b);

	return B;
}

// Decode a byte array, B in B^32d, into an array of d-bit integers 
// for 1<=d<=12
//
// NOTE: For 1<=d<=11, the conversion is one-to-one; whereas, for d=12 
// 	it is no longer a one-to-one operation when the input byte array
// 	is NOT produced by ByteEncode()
static union integer* ByteDecode(union byte* B, unsigned int d) {
	union bit* b;
	union integer* F;
	unsigned int m;

	b = BytesToBits(B,32*d);
	F = malloc(sizeof(union integer) * 256);

	// Set the value of 'm' appropriately
	if ((1 <= d) && (d < 12)) m = 0x001 << d;
	else if (d == 12) m = Q;
	else exit(EXIT_FAILURE);

	for (int i=0; i < 256; i++) {
		F[i].t = 0x000;		// initialize F[i]

		// F[i] = sum_j=0->(d-1){(b[i*d+j]*2^j) mod m}
		for (int j=0; j < d; j++) {
			F[i].t |= (b[(i * d) + j].b * (0x001 << j)) % m;
		}			
	}

	free(b);

	return F;
}

/****************************************
 * Sampling Algorithms
 * *************************************/

// Uniform sampling of NTT representations - converts a seed together
// with two indexing bytes into a polynomial in the NTT domain (i.e., Tq)
//
// Note: The input byte array is of length 34 and the ouput is an array 
// 	in Zq^256 that contains the coefficients of the sampled element of 
// 	Tq (FIPS-203:4.2.2)
union integer* SampleNTT(union byte* B) {
	union integer* a;
	union bit* b;
	union bit* S;
	union byte* C;
	union integer I[3];
	union integer d1, d2;
	unsigned int j = 0; 
	int k = 0;

	a = malloc(sizeof(union integer) * 256);
	C = malloc(sizeof(union byte) * 3);

	b = BytesToBits(B, 34);						// Convert B to bit array to inject into SHA algo
	S = sha3_b(b, 34*8, 280*8*3, 256, (union bit[]){1,1,1,1});	// XOF.Init() & XOF.Absorb(B)	

	while (j < 256) {
		C = BitsToBytes(S+k, 8*3);				// XOF.Squeeze(3) 
		for (int i=0; i < 3; i++) I[i].t = C[i].e;		// Increase address space of C 

		// Note that 0 <= d1,d2 < 2^12
		d1.t = I[0].t + 256 * (I[1].t % 16);			// d1 <- C[0] + 256 * (C[1] % 16)
		d2.t = (I[1].t / 16) + (16 * I[2].t);			// d2 <- floor(C[1] / 16) + 16 * C[2]

		if (d1.t < Q) {
			a[j] = d1;
			j += 1;
		}

		if ((d2.t < Q) && (j < 256)) {
			a[j] = d2;
			j += 1;
		}

		// Sets a limit on the number of iterations -
		// see FIPS-203:A-B
		k += (8 * 3);
		if (k >= (280*8*3 - 8*3)) {
			k = -1;
			break;
		}
	}

	free(b);
	free(S);
	free(C);

	// If limit has been reached, then change the
	// last 2 indices of B and try again
	if (k == -1) {
		free(a);
		B[32].e += 1;
		B[33].e += 1;
		return SampleNTT(B);
	}

	return a;
}

// Sample from the Centered Binomial Distribution - given a stream of uniformly
// random bytes of length 64*n, derive the coefficient array of a polynomial f in 
// Rq according the distribution Dn(Rq) (see FIPS-203:4.2.2)
//
// Note: The distribution Dn(Rq) of polynomials in Rq are sometimes referred to
// 	as "errors" or "noise"
union integer* SamplePolyCBD(union byte* B, unsigned int n) {
	union integer* f;
	union bit* b;
	union integer x, y;

	b = BytesToBits(B, 64*n);					// Convert B to a bit array
	f = malloc(sizeof(union integer) * 256);			// f is in Zq^256

	for (int i=0; i < 256; i++) {
		x.t = 0;
		y.t = 0;
		for (int j=0; j < n; j++) x.t += b[2*i*n+j].b;		// 0 <= x <= n
		for (int j=0; j < n; j++) y.t += b[2*i*n+n+j].b;	// 0 <= y <= n

		// 0<=f[i]<=n OR Q-n<=f[i]<=Q-1
		if (x.t >= y.t) f[i].t = x.t - y.t;
		else f[i].t = Q - (y.t - x.t);
	}

	return f;
}

/****************************************
 * Number-Theoretic Transform
 ***************************************/

union integer* NTT(union integer* f) {
	union integer zeta, t;		// 24-bit integers
	union byte i, pow;		// 7-bits bytes 
	union integer* fh;		// fh = NTT of f in Zq^256

	// Copy input array into fh
	fh = malloc(sizeof(union integer) * 256);
	for (int j=0; j < 256; j++) fh[j] = f[j];

	i.s = 1;
	for (int len=128; len >= 2; len /= 2) {
		for (int start=0; start < 256; start += (2*len)) {

			// Derive the value of zeta using BitRev7(i) -
			// equivalent to zeta = (17^BitRev7(i)) % Q
			zeta.l = 1;
			pow = BitRev7(i);
			while (pow.s > 0) { 
				zeta.l = (zeta.l * 17) % Q;
				pow.s -= 1;
			}
			i.s += 1;

			// Update fh in place 
			for (int j=start; j < (start + len); j++) {
				// Calculate t - note that 3328*3328 < 2^24
				t.l = (zeta.l * fh[j+len].t) % Q;

				// Set fh[j+len] - workaround for negative mod arithmatic
				// with unsigned int
				if (fh[j].t >= t.l) fh[j+len].t = fh[j].t - t.l;	
				else fh[j+len].t = Q - (t.l - fh[j].t);

				// Set fh[j] - note that addition needs to be done
				// in the 24-bit address space
				t.l = (fh[j].t + t.l) % Q;
				fh[j].t = t.l;
			}
		}
	}

	return fh;
}

union integer* InverseNTT(union integer* fh) {
	union integer* f;
	union byte i, pow;		// 7-bit bytes
	union integer zeta, t, temp;	// 24-bit integers

	// Copy input array into f
	f = malloc(sizeof(union integer) * 256);
	for (int i=0; i < 256; i++) f[i] = fh[i];

	i.s = 127;
	for (int len=2; len <= 128; len *= 2) {
		for (int start=0; start < 256; start += (2*len)) {
			
			// Derive the value of zeta using BitRev7(i) -
			// equivalent to zeta = (17^BitRev7(i)) % Q
			zeta.l = 1;
			pow = BitRev7(i);
			while (pow.s > 0) { 
				zeta.l = (zeta.l * 17) % Q;
				pow.s -= 1;
			}
			i.s -= 1;

			for (int j=start; j < (start + len); j++) {
				t.l = f[j].t;		// Increase address space

				// Set fh[j] - note that addition needs to be done
				// in the 24-bit address space
				temp.l = (t.l + f[j+len].t) % Q;
				f[j].t = temp.l;

				// Set fh[j+len] - workaround for negative mod arithmatic
				// with unsigned int
				if (f[j+len].t >= t.l) temp.l = f[j+len].t - t.l;
				else temp.l = Q - (t.l - f[j+len].t);
				temp.l = (zeta.l * temp.l) % Q;
				f[j+len].t = temp.l;
			}
		}
	}

	for (int j=0; j < 256; j++) {
		temp.l = (f[j].t * 3303) % Q;
		f[j].t = temp.l;
	}

	return f;
}
