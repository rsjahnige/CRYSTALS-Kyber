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
static union byte* BitsToBytes(const union bit* b, unsigned int l) {
	union byte* B;
	B = malloc(sizeof(union byte) * l / 8);

	for (int i=0; i < l; i++) {
		if ((i % 8) == 0) B[i/8].e = 0b00000000;	// initialize byte array index
		B[i/8].e |= (b[i].b & 0b00000001) << (i % 8);	// B[i/8] = B[i/8] + b[i] * 2^(i%8)
	}

	return B;
}

// Perfroms the inverse of BitToBytes, converting a byte array into a bit array
//
// NOTE: The output bit array is in little-endian order
static union bit* BytesToBits(const union byte* B, unsigned int L) {
	union bit* b;
	union byte C[L];

	b = malloc(sizeof(union bit) * L * 8);

	for (int i=0; i < L; i++) {
		C[i].e = B[i].e;					// Copy B[i] into array C[i] in B^L
		for (int j=0; j < 8; j++) {
			b[(8 * i) + j].b = C[i].e & 0b00000001;		// b[8i+j] = C[i] % 2 
			C[i].e = C[i].e >> 1;				// C[i] = C[i] / 2
		}
	}

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
static union byte* ByteEncode(const union integer* F, unsigned int d) {
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
static union integer* ByteDecode(const union byte* B, unsigned int d) {
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
	b = BytesToBits(B, 34);						// Convert B to bit array
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
		
		free(C);
	}

	free(b);
	free(S);

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
union integer* SamplePolyCBD(const union byte* B, unsigned int n) {
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

	free(b);

	return f;
}

/****************************************
 * Number-Theoretic Transform Functions
 ***************************************/

// Transform an element f in Rq to an element fh in Tq - the ring
// Rq is isomorphic to the ring Tq 
//
// Note: Tq is a product of 128 rings that each consist of polynomials 
// of degree at most one (i.e., fh[i]=fh[2i]+fh[2i+1]X for i in 
// {0,...,127}) 
union integer* NTT(const union integer* f) {
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

// Transform an element fh in Tq to an element f in Rq - the isomorphism
// between the rings mentioned above implies f*g = InverseNTT(fh*gh)
//
// Note: Rq is the ring Zq[X]/(X^n+1) of polynomials of the form 
// f = f0 + f1*X + ... + f255*X^255 for fi in Zq for all i
union integer* InverseNTT(const union integer* fh) {
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

	// Multiply each enty by 3303 mod Q
	for (int j=0; j < 256; j++) {
		temp.l = (f[j].t * 3303) % Q;
		f[j].t = temp.l;
	}

	return f;
}

/**************************************
 * Multiplication in the NTT dommain
 * ***********************************/

// Computes the product of two degree-one polynomials with respect 
// to a quadratic modulus
//
// NOTE: The coefficients are a0+a1X and b0+b1X. The quadratic modulus is 
// 	X^2-gamma where gamma iss 17^(2*BitRev7(i)+1)
static union integer* BaseCaseMultiply(union integer a0, union integer a1, union integer b0,
					union integer b1, union integer gamma) {
	union integer* C;
	union integer temp;

	C = malloc(sizeof(union integer) * 2);

	temp.l = (a1.t * b1.t) % Q;
	temp.l = (temp.l * gamma.l) % Q;
	temp.l += (a0.t * b0.t) % Q;
	C[0].t = temp.l % Q;

	temp.l = (a0.t * b1.t) % Q;
	temp.l += (a1.t * b0.t) % Q;
	C[1].t = temp.l % Q;

	return C;
}

// Computes the product, in the ring Tq, of two NTT representations
static union integer* MultiplyNTTs(const union integer* fh, const union integer* gh) {
	union integer* C;
	union integer* h;
	union byte j, pow;		// j is 7 bit; pow is 8 bits
	union integer gamma;

	h = malloc(sizeof(union integer) * 256);

	for (int i=0; i < 128; i++) {
		// Derive the value of gamma using BitRev7(i) -
		// equivalent to gamma = (17^(2*BitRev7(j)+1)) % Q
		j.s = i;
		gamma.l = 1;
		j = BitRev7(j);
		pow.e = 2*j.s + 1;		// Note that BitRev7(j)<=127 so pow<=255
		while (pow.e > 0) { 
			gamma.l = (gamma.l * 17) % Q;
			pow.e -= 1;
		}

		C = BaseCaseMultiply(fh[2*i], fh[2*i+1], gh[2*i], gh[2*i+1], gamma);
		h[2*i].t = C[0].t;
		h[2*i+1].t = C[1].t;
		free(C);
	}

	return h;
}

/*************************************
 * K-PKE Component Scheme Algorithms
 ************************************/

// The G function takes a variable length input `seed` of length `len`
// and produces a 64-byte output that can be split into two 32-byte 
// outputs (FIPS-203:4.1)
static union byte* G(const union byte* seed, unsigned int len) {
	union bit* hash;
	union bit* input;
	union byte* result;

	input = BytesToBits(seed, len);
	hash = sha3_b(input, len*8, 512, 512*2, (union bit[]){0,1,0,0});
	result = BitsToBytes(hash, 512);

	free(input);
	free(hash);

	return result;
}

// The Pseudorandom Function (PRF) takes a parameter `n` in {2,3}, one
// 32-byte input, `s`, and one 1-byte input, `b`
//
// It produces a (64*n)-byte output using SHAKE256
static union byte* PRF(const union byte* s, union byte b, unsigned int n) {
	union bit* hash;
	union bit* input;
	union byte seed[33];
	union byte* result;

	for (int i=0; i < 32; i++) {
		seed[i].e = s[i].e;	
	}
	seed[32].e = b.e;

	input = BytesToBits(seed, 33);
	hash = sha3_b(input, 33*8, 8*64*n, 256, (union bit[]){1,1,1,1});
	result = BitsToBytes(hash, 8*64*n);

	free(input);
	free(hash);

	return result;
}

static union integer* PolyAddition(const union integer* u, const union integer* v) {
	union integer temp;
	union integer* z;

	z = malloc(sizeof(union integer) * 256);

	for (int i=0; i < 256; i++) {
		temp.l = u[i].t;			// increase address space for addition
		z[i].t = (temp.l + v[i].t) % Q;		// perform addition and cast back to 12-bit address space
	}

	return z;
}

static union integer* VectorMultiply(union integer** u, union integer** v,
					unsigned int k) {
	union integer* w;
	union integer* z;
	union integer* y;

	// Perform vector-vector multiplication in the NTT domain -
	// see FIPS-203:2.4.7
	w = MultiplyNTTs(u[0], v[0]);
	for (int i=1; i < k; i++) {
		z = MultiplyNTTs(u[i], v[i]);	
		y = PolyAddition(w, z);

		free(w);
		free(z);

		w = y;
	}

	return w;
}


struct PKE KeyGen(const struct ML_KEM* params, const union byte* d) {

	// Temporary arrays
	union byte* output;
	union integer* sample;
	union integer* temp;

	union byte rho[34];				// matrix `A` seed + 2 bytes
	union byte sigma[32];				// secret `s` and noise `e` seed
	union byte rand[33];				// randomness input + one byte (d || params.k)

	union byte n;
	union integer* A[params->k.e][params->k.e];	// pseudorandom coefficients
	union integer* s[params->k.e];			// secret 
	union integer* e[params->k.e];			// noise
	union integer* t[params->k.e];			// collection of "noisy" linear equations

	struct PKE keys;				// function output

	n.e = 0;	// initialize `n` to 0

	// Concatenate byte array d with params->k
	for (int i=0; i < 32; i++) rand[i] = d[i];
	rand[32] = params->k;

	output = G(rand, 33);

	// Split `output` into two pseudorandom 32-byte seeds
	for (int i=0; i < 32; i++) rho[i] = output[i];
	for (int i=0; i < 32; i++) sigma[i] = output[i+32];

	free(output);

	// Generate matrix A in (Zq^256)^(kxk)
	for (int i=0; i < params->k.e; i++) {
		for (int j=0; j < params->k.e; j++) {
			rho[32].e = j;
			rho[33].e = i;
			A[i][j] = SampleNTT(rho);
		}

	}
	
	// Generate vector s in the NTT domain
	for (int i=0; i < params->k.e; i++) {

		output = PRF(sigma, n, params->n1.e);
		sample = SamplePolyCBD(output, params->n1.e);

		s[i] = NTT(sample);
		n.e += 1;

		free(output);
		free(sample);
	}
	
	
	// Generate vector e in the NTT domain
	for (int i=0; i < params->k.e; i++) {

		output = PRF(sigma, n, params->n1.e);
		sample = SamplePolyCBD(output, params->n1.e);

		e[i] = NTT(sample);
		n.e += 1;

		free(output);
		free(sample);
	}

	// Generate vector t in the NTT domain	
	for (int i=0; i < params->k.e; i++) {
		
		temp = VectorMultiply(A[i], s, params->k.e);
		t[i] = PolyAddition(temp, e[i]);
		free(temp);
	}

	// Allocate memory for PKE struct
	keys.ek = malloc(sizeof(union byte*) * params->k.e);
	keys.dk = malloc(sizeof(union byte*) * params->k.e);
	keys.rho = malloc(sizeof(union byte) * 32);
	
	// Copy first 32 indices of local `rho` array into `keys.rho` struct
	for (int i=0; i < 32; i++) {
		keys.rho[i] = rho[i];
	}

	// Run ByteEncode k times to generate the encryption and
	// decryption keys
	for (int i=0; i < params->k.e; i++) {
		keys.ek[i] = ByteEncode(t[i], 12);
		keys.dk[i] = ByteEncode(s[i], 12);
	}

	// Free up dynamically allocated memory
	for (int i=0; i < params->k.e; i++) {
		for (int j=0; j < params->k.e; j++) {
			free(A[i][j]);
		}
		free(s[i]);
		free(e[i]);
		free(t[i]);
	}

	return keys;
}

union byte* Encrypt(const struct ML_KEM* params, const struct PKE* keys, 
			const union byte* m, const union byte* r) {

	// Ciphertext arrays
	union byte* c;			// c = c1 || c2
	union byte* c1[params->k.e];
	union byte* c2;
	unsigned int c1_len, c2_len;

	// Temporary arrays
	union byte* output;
	union integer* sample;
	union integer* temp1;
	union integer* temp2;
				
	union byte n;
	union byte rho[34];
	union integer* t[params->k.e];
	union integer* y[params->k.e];
	union integer* e1[params->k.e];
	union integer* e2;
	union integer* u[params->k.e];
	union integer* v;
	union integer mu[256];
	union integer* At[params->k.e][params->k.e];	// Transpose of matrix A

	n.e = 0;	// initialize n to 0

	// Run ByteDecode() k times to decode t in (Zq^256)^k
	for (int i=0; i < params->k.e; i++) {
		t[i] = ByteDecode(keys->ek[i], 12); 
	}

	// Copy 32-byte set to local array
	for (int i=0; i < 32; i++) {
		rho[i].e = keys->rho[i].e;
	}

	// Re-generate matrix A in (Zq^256)^(kxk) sampled in KenGen()
	for (int i=0; i < params->k.e; i++) {
		for (int j=0; j < params->k.e; j++) {
			rho[32].e = j;
			rho[33].e = i;
			At[j][i] = SampleNTT(rho);
		}
	}

	// Generate y i (Zq^256)^k sampled from the CBD in the NTT domain
	for (int i=0; i < params->k.e; i++) {
		
		output = PRF(r, n, params->n1.e);
		sample = SamplePolyCBD(output, params->n1.e);

		y[i] = NTT(sample);
		n.e += 1;

		free(output);
		free(sample);
	}

	// Generate e1 in (Zq^256)^k sampled from the CBD
	for (int i=0; i < params->k.e; i++) {

		output = PRF(r, n, params->n2.e);
		e1[i] = SamplePolyCBD(output, params->n2.e);
		n.e += 1;

		free(output);
	}

	// Sample e2 in Zq^256 from the CBD
	output = PRF(r, n, params->n2.e);
	e2 = SamplePolyCBD(output, params->n2.e);
	free(output);

	// Generate vector `u` in the ring Rq - NTT^-1(At * y) + e1
	for (int i=0; i < params->k.e; i++) {
		
		temp1 = VectorMultiply(At[i], y, params->k.e);	// Perfrom vector multiplication in the NTT domain
		temp2 = InverseNTT(temp1);			// Transform `temp1` from Tq to Rq

		// Add vector `e1[i]` to vector u[i] in the ring Rq
		u[i] = PolyAddition(temp2, e1[i]);

		free(temp1);
		free(temp2);
	}

	temp1 = ByteDecode(m, 1);
	for (int i=0; i < 256; i++) {
		mu[i] = Decompress(temp1[i], 1);
	}
	free(temp1);

	// Encode plaintext `m` into polynomial `v`
	temp1 = VectorMultiply(t, y, params->k.e);
	temp2 = InverseNTT(temp1);

	free(temp1);

	temp1 = PolyAddition(temp2, e2);
	v = PolyAddition(temp1, mu);

	free(temp1);
	free(temp2);

	// Generate ciphertext c1
	for (int i=0; i < params->k.e; i++) {

		temp1 = malloc(sizeof(union integer) * 256);

		for (int j=0; j < 256; j++) {
			temp1[j] = Compress(u[i][j], params->du.e);
		}

		c1[i] = ByteEncode(temp1, params->du.e);
		free(temp1);
	}

	// Generate ciphertext c2
	temp2 = malloc(sizeof(union integer) * 256);
	for (int i=0; i < 256; i++) {
		temp2[i] = Compress(v[i], params->dv.e);
	}
	c2 = ByteEncode(temp2, params->dv.e);
	free(temp2);

	// Concatenate ciphertext c1 and c2
	c1_len = 32 * params->du.e;
	c2_len = 32 * params->dv.e;
	c = malloc(sizeof(union byte) * (c1_len * params->k.e + c2_len));
	for (int i=0; i < params->k.e; i++) {
		for (int j=0; j < c1_len; j++) {
			c[i*c1_len+j] = c1[i][j];
		}
	}
	c1_len *= params->k.e;
	for (int i=0; i < c2_len; i++) {
		c[i+c1_len] = c2[i];
	}

	// Free up dynamically allocated memory
	for (int i=0; i < params->k.e; i++) {
		for (int j=0; j < params->k.e; j++) {
			free(At[i][j]);
		}
		free(y[i]);
		free(e1[i]);
		free(u[i]);
		free(c1[i]);
		free(t[i]);
	}
	free(e2);
	free(v);
	free(c2);

	return c;
}
