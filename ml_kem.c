/********************************************************
 * File: ml_kem.c
 * Author: Ryan Jahnige
 *
 * Description: Module-Lattice-Based Key-Encapsulation
 * 		Mechanism Standard (FIPS-203)
 * *****************************************************/

#include "ml_kem.h"

// Data type for integer arrays F in Zm^256, where m=2^d if 1<=d<12 
// and m=Q if d=12 (see FIPS-203:4.2.1)
union integer {
	unsigned int t : 12;	// twelve bit integer 	
	unsigned int l : 24; 	// twenty-four bit integer - used for real number calculations
};

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
// NOTE: F is in the finite field Zm^N; therefore, the lenght of F 
// 	shall always be N (see FIPS-203:2.3)
static union byte* ByteEncode(const union integer* F, unsigned int d) {
	union integer a;
	union bit* b;
	union byte* B;

	b = malloc(sizeof(union bit) * N * d);

	for (int i=0; i < N; i++) {
		a.t = F[i].t;					// 'a' in Zm
		for (int j=0; j < d; j++) {	
			b[(i * d) + j].b = a.t & 0x001;		// b[i*d+j] = a % 2	
			a.t >>= 1;				// a = (a - b[i*d+j])/2
		}

	}
	
	B = BitsToBytes(b,N*d);
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
	F = malloc(sizeof(union integer) * N);

	// Set the value of 'm' appropriately
	if ((1 <= d) && (d < 12)) m = 0x001 << d;
	else if (d == 12) m = Q;
	else exit(EXIT_FAILURE);

	for (int i=0; i < N; i++) {
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
// 	in Zq^N that contains the coefficients of the sampled element of 
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

	a = malloc(sizeof(union integer) * N);
	b = BytesToBits(B, 34);						// Convert B to bit array
	S = sha3_b(b, 34*8, 280*8*3, N, (union bit[]){1,1,1,1});	// XOF.Init() & XOF.Absorb(B)	

	while (j < N) {
		C = BitsToBytes(S+k, 8*3);				// XOF.Squeeze(3) 
		for (int i=0; i < 3; i++) I[i].t = C[i].e;		// Increase address space of C 

		// Note that 0 <= d1,d2 < 2^12
		d1.t = I[0].t + N * (I[1].t % 16);			// d1 <- C[0] + N * (C[1] % 16)
		d2.t = (I[1].t / 16) + (16 * I[2].t);			// d2 <- floor(C[1] / 16) + 16 * C[2]

		if (d1.t < Q) {
			a[j] = d1;
			j += 1;
		}

		if ((d2.t < Q) && (j < N)) {
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
	f = malloc(sizeof(union integer) * N);				// f is in Zq^N

	for (int i=0; i < N; i++) {
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
	union integer* fh;		// fh = NTT of f in Zq^N

	// Copy input array into fh
	fh = malloc(sizeof(union integer) * N);
	for (int j=0; j < N; j++) fh[j] = f[j];

	i.s = 1;
	for (int len=128; len >= 2; len /= 2) {
		for (int start=0; start < N; start += (2*len)) {

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
	f = malloc(sizeof(union integer) * N);
	for (int i=0; i < N; i++) f[i] = fh[i];

	i.s = 127;
	for (int len=2; len <= 128; len *= 2) {
		for (int start=0; start < N; start += (2*len)) {
			
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
	for (int j=0; j < N; j++) {
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
// 	X^2-gamma where gamma is 17^(2*BitRev7(i)+1)
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

	h = malloc(sizeof(union integer) * N);

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
 * Auxiliary Algorithms
 ************************************/

// Entropy source for random byte generation - the return value
// for each getRandomBytes() function is an array of bytes of
// length 32
//
// NOTE: The entropy sources SHALL BE approved RBG's as prescribed
// 	in SP 800-90A/B/C

#ifdef __linux__

// For linux machines, read entropy from /dev/urandom
static union byte* getRandomBytes() {
	int fd;
	unsigned int rnd[32];
	union byte* random;

	fd = open("/dev/urandom", O_RDONLY);
	if (fd < 0) return NULL;

	if (read(fd, rnd, 32*sizeof(unsigned int)) != 32*sizeof(unsigned int)) {
		close(fd);
		return NULL;
	}
	close(fd);

	random = malloc(sizeof(union byte) * 32);
	for (int i=0; i < 32; i++) {
		random[i].e = rnd[i] % N;
	}

	return random;
}	

#else

// Unknown entropy source for current operating system - returning NULL
// forces the program to terminate
static union byte* getRandomBytes() {
	printf("ml_kem.c:getRandomBytes() :: Operating system not supported\n");
	return NULL;
}

#endif

// The Pseudorandom Function (PRF) takes a parameter `n` in {2,3}, one
// 32-byte input, `s`, and one 1-byte input, `b` and produces a (64*n)-byte 
// output (FIPS-203:4.1)
//
// PRFn(s,b) := SHAKE256(s||b,8*64*n)
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
	hash = sha3_b(input, 33*8, 8*64*n, N, (union bit[]){1,1,1,1});
	result = BitsToBytes(hash, 8*64*n);

	free(input);
	free(hash);

	return result;
}

// The H function takes a variable length input `seed` of length `len`
// and produces a 32-byte output (FIPS-203:4.1)
//
// H(s) := SHA3-256(s)
static union byte* H(const union byte* seed, unsigned int len) {
	union bit* hash;
	union bit* input;
	union byte* result;

	input = BytesToBits(seed, len);
	hash = sha3_b(input, len*8, N, N*2, (union bit[]){0,1,0,0});
	result = BitsToBytes(hash, N);

	free(input);
	free(hash);

	return result;
}

// The J function takes a variable length input `seed` of length `len`
// and produces a 32-byte output (FIPS-203:4.1)
//
// J(s) := SHAKE256(s, 8*32)
static union byte* J(const union byte* seed, unsigned int len) {
	union bit* hash;
	union bit* input;
	union byte* result;

	input = BytesToBits(seed, len);
	hash = sha3_b(input, len*8, 8*32, N, (union bit[]){1,1,1,1});
	result = BitsToBytes(hash, 8*32);

	free(input);
	free(hash);

	return result;	
}

// The G function takes a variable length input `seed` of length `len`
// and produces a 64-byte output (FIPS-203:4.1)
//
// G(c) := SHA3-512(c)
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


// PolyAddition() takes as input two arrays in the field Zq^N and computes the additive
// result of each index modulo Q. The result is an element of Zq^N
//
// NOTE: Arithmatic with coefficent arrays in Rq or Tq are both performed coordinatewise, 
// so this function is applicable for either field (FIPS-203:2.4.5)
static union integer* PolyAddition(const union integer* u, const union integer* v) {
	union integer temp;
	union integer* z;

	z = malloc(sizeof(union integer) * N);

	for (int i=0; i < N; i++) {
		temp.l = u[i].t;			// increase address space for addition
		z[i].t = (temp.l + v[i].t) % Q;		// perform addition and assign to 12-bit address space
	}

	return z;
}

// PolySubtraction() takes as input two arrays in the field Zq^N and computes the
// difference between each index. The result is an element of Zq^N
//
// NOTE: Arithmatic with coefficent arrays in Rq or Tq are both performed coordinatewise, 
// so this function is applicable for either field (FIPS-203:2.4.5)
static union integer* PolySubtraction(const union integer* u, const union integer* v) {
	union integer* z;

	z = malloc(sizeof(union integer) * N);

	for (int i=0; i < N; i++) {
		if (u[i].t < v[i].t) {				// if u[i] - v[i] is negative
			z[i].t = Q - (v[i].t - u[i].t);		// negative int modulo Q logic
		} else {
			z[i].t = u[i].t - v[i].t;		// u[i] - v[i] >= 0 - find difference
		}	
	}

	return z;
}

// VectorMultiply() takes as input two arrays of length `k` of polynominals 
// in the NTT domain and produces a result in the base ring Tq, represented as 
// an element of Zq^N
static union integer* VectorMultiply(union integer** u, union integer** v,
					unsigned int k) {
	union integer* w;
	union integer* z;
	union integer* y;

	// Perform array-array multiplication in the NTT domain -
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

/******************************************
 * K-PKE Component Scheme Algorithms - NOT 
 * APPROVED for use in a standalone fashion
 *****************************************/

// Use randomness, d of length 32, to generate an encryption key and 
// corresponding decryption key. See getRandomBytes() functions above
// for generating randomness.
//
// NOTE: Any application that calls this function is responsible for
// 	dealocating the byte arrays PKE.ek an PKE.dk
static struct PKE PKE_KeyGen(const struct PARAMS* params, const union byte* d) {

	// Function output
	struct PKE keys;

	// Temporary arrays
	union byte* output;
	union integer* sample;
	union integer* temp;
	union byte rand[33];				// randomness input + one byte (d || params.k)

	// K-PKE.KeyGen() defined variables
	union byte n;
	union byte rho[34];				// matrix `A` seed + 2 bytes
	union byte sigma[32];				// secret `s` and noise `e` seed
	union integer* A[params->k.e][params->k.e];	// pseudorandom coefficients
	union integer* s[params->k.e];			// secret 
	union integer* e[params->k.e];			// noise
	union integer* t[params->k.e];			// collection of "noisy" linear equations

	n.e = 0;	// initialize n to 0

	// Concatenate byte array d with params->k
	for (int i=0; i < 32; i++) rand[i] = d[i];
	rand[32] = params->k;

	output = G(rand, 33);

	// Split `output` into two pseudorandom 32-byte seeds
	for (int i=0; i < 32; i++) rho[i] = output[i];
	for (int i=0; i < 32; i++) sigma[i] = output[i+32];

	free(output);

	// Generate matrix A in (Zq^N)^(kxk)
	for (int i=0; i < params->k.e; i++) {
		for (int j=0; j < params->k.e; j++) {
			rho[32].e = j;
			rho[33].e = i;
			A[i][j] = SampleNTT(rho);
		}

	}
	
	// Generate array s in the NTT domain
	for (int i=0; i < params->k.e; i++) {

		output = PRF(sigma, n, params->n1.e);
		sample = SamplePolyCBD(output, params->n1.e);

		s[i] = NTT(sample);
		n.e += 1;

		free(output);
		free(sample);
	}
	
	
	// Generate array e in the NTT domain
	for (int i=0; i < params->k.e; i++) {

		output = PRF(sigma, n, params->n1.e);
		sample = SamplePolyCBD(output, params->n1.e);

		e[i] = NTT(sample);
		n.e += 1;

		free(output);
		free(sample);
	}

	// Generate a noisy linear system, t, in the NTT domain	
	for (int i=0; i < params->k.e; i++) {
		
		temp = VectorMultiply(A[i], s, params->k.e);
		t[i] = PolyAddition(temp, e[i]);
		free(temp);
	}

	// Allocate memory for PKE struct
	keys.ek_len = 384 * params->k.e + 32;
	keys.dk_len = 384 * params->k.e;
	keys.ek = malloc(sizeof(union byte) * keys.ek_len);
	keys.dk = malloc(sizeof(union byte) * keys.dk_len);

	// Run ByteEncode k times to generate the encryption key
	for (int i=0; i < params->k.e; i++) {
		output = ByteEncode(t[i], 12);
		for (int j=0; j < 384; j++) {
			keys.ek[384*i+j] = output[j];
		}
		free(output);
	}

	// Append first 32 indices of local `rho` array to the end of `keys.ek`
	for (int i=0; i < 32; i++) {
		keys.ek[384*params->k.e+i] = rho[i];
	}

	// Run ByteEncode k times to generate the decryption key
	for (int i=0; i < params->k.e; i++) {
		output = ByteEncode(s[i], 12);
		for (int j=0; j < 384; j++) {
			keys.dk[384*i+j] = output[j];
		}
		free(output);
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

// Use the encryption key, ek, generated by PKE_KeyGen() to encrypt a plaintext
// message, m of length 32, using the randomnes r of length 32. See getRandomBytes()
// function description above.
//
// NOTE: The returned ciphertext, c, is of length 32*(du*k+dv)
static union byte* PKE_Encrypt(const struct PARAMS* params, const union byte* ek, 
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
	
	// K-PKE.Encrypt() defined variables	
	union byte n;
	union byte rho[34];
	union integer mu[N];				
	union integer* t[params->k.e];
	union integer* y[params->k.e];
	union integer* e1[params->k.e];
	union integer* e2;
	union integer* u[params->k.e];
	union integer* v;
	union integer* At[params->k.e][params->k.e];	// Transpose of matrix A

	n.e = 0;	// initialize n to 0

	// Run ByteDecode() k times to decode t in (Zq^N)^k
	for (int i=0; i < params->k.e; i++) {
		t[i] = ByteDecode(ek+(384*i), 12); 
	}

	// Copy 32-byte set to local array
	for (int i=0; i < 32; i++) {
		rho[i] = ek[384*params->k.e + i];
	}

	// Re-generate matrix A in (Zq^N)^(kxk) sampled in KenGen() - 
	// note that the A transpose is generated here
	for (int i=0; i < params->k.e; i++) {
		for (int j=0; j < params->k.e; j++) {
			rho[32].e = j;
			rho[33].e = i;
			At[j][i] = SampleNTT(rho);
		}
	}

	// Generate y i (Zq^N)^k sampled from the CBD in the NTT domain
	for (int i=0; i < params->k.e; i++) {
		
		output = PRF(r, n, params->n1.e);
		sample = SamplePolyCBD(output, params->n1.e);

		y[i] = NTT(sample);
		n.e += 1;

		free(output);
		free(sample);
	}

	// Generate e1 in (Zq^N)^k sampled from the CBD
	for (int i=0; i < params->k.e; i++) {

		output = PRF(r, n, params->n2.e);
		e1[i] = SamplePolyCBD(output, params->n2.e);
		n.e += 1;

		free(output);
	}

	// Sample e2 in Zq^N from the CBD
	output = PRF(r, n, params->n2.e);
	e2 = SamplePolyCBD(output, params->n2.e);
	free(output);

	// Generate array u in the ring Rq - NTT^-1(At * y) + e1
	for (int i=0; i < params->k.e; i++) {
		
		temp1 = VectorMultiply(At[i], y, params->k.e);	// NTT domain (i.e., Tq)
		temp2 = InverseNTT(temp1);			// Transform to Rq

		// Add noise array e1 to At*y - u[i] in the ring Rq
		u[i] = PolyAddition(temp2, e1[i]);

		free(temp1);
		free(temp2);
	}

	// Encode the input message m as mu
	temp1 = ByteDecode(m, 1);
	for (int i=0; i < N; i++) {
		mu[i] = Decompress(temp1[i], 1);
	}
	free(temp1);

	// Transform encoded message mu into polynomial v
	temp1 = VectorMultiply(t, y, params->k.e);		// NTT domain (i.e., Tq)
	temp2 = InverseNTT(temp1);				// Transform to Rq

	free(temp1);

	temp1 = PolyAddition(temp2, e2);
	v = PolyAddition(temp1, mu);

	free(temp1);
	free(temp2);

	// Generate ciphertext c1
	for (int i=0; i < params->k.e; i++) {

		temp1 = malloc(sizeof(union integer) * N);

		for (int j=0; j < N; j++) {
			temp1[j] = Compress(u[i][j], params->du.e);
		}

		c1[i] = ByteEncode(temp1, params->du.e);
		free(temp1);
	}

	// Generate ciphertext c2
	temp2 = malloc(sizeof(union integer) * N);
	for (int i=0; i < N; i++) {
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

// Use the decryption key, dk, generated by PKE_KeyGen() to decrypt a 
// ciphertext, c, generated by PKE_Encrypt()
//
// NOTE: The returned message, m, is of length 32
static union byte* PKE_Decrypt(const struct PARAMS* params, const union byte* dk, 
				const union byte* c) {

	// Function ouptut
	union byte* m;

	// Ciphertext variables
	union byte c1[params->k.e][32*params->du.e];
	union byte c2[32*params->dv.e];
	unsigned int c1_len, c2_len;

	// Temporary arrays
	union integer* temp1;
	union integer* temp2;

	// K-PKE.Decrypt() defined variables
	union integer* u[params->k.e];
	union integer* v;
	union integer* s[params->k.e];
	union integer* w;

	// Split ciphertext into 2 byte arrays
	c1_len = 32 * params->du.e;
	c2_len = 32 * params->dv.e;
	for (int i=0; i < params->k.e; i++) {
		for (int j=0; j < c1_len; j++) {
			c1[i][j] = c[i*c1_len + j];
		}
	}

	c1_len *= params->k.e;
	for (int i=0; i < c2_len; i++) {
		c2[i] = c[i+c1_len];
	}

	// Generate coefficient array u in the NTT domain - u in (R2^du)^k		
	for (int i=0; i < params->k.e; i++) {
		temp1 = ByteDecode(c1[i], params->du.e);
		
		for (int j=0; j < N; j++) {
			temp1[j] = Decompress(temp1[j], params->du.e);
		}

		u[i] = NTT(temp1);
		free(temp1);
	}

	// Generate the constant term v in the ring R2^dv
	v = ByteDecode(c2, params->dv.e);
	for (int i=0; i < N; i++) {
		v[i] = Decompress(v[i], params->dv.e);
	}

	// Decode the decryption key dk into the secrect variable s in Rq
	for (int i=0; i < params->k.e; i++) {
		s[i] = ByteDecode(dk+(384*i), 12);
	}

	// Calculate w = v - InvertseNTT(st * u) resulting in w in Rq
	temp1 = VectorMultiply(s, u, params->k.e);
	temp2 = InverseNTT(temp1);
	w = PolySubtraction(v, temp2);

	free(temp1);
	free(temp2);

	// Decode plaintext message m from w
	for (int i=0; i < N; i++) {
		w[i] = Compress(w[i], 1);
	}
	m = ByteEncode(w, 1);

	// Free up dynamically allocated memory
	for (int i=0; i < params->k.e; i++) {
		free(u[i]);
		free(s[i]);
	}
	free(v);
	free(w);

	return m;
}

/*************************************
 * ML-KEM Internal Algorithms
 ************************************/

// Use randomness, d and z both of length 32, to generate a KEM encapsulation key 
// of length 384*k+32 and a cooresponding decapsulation key of length 768*k+96. The 
// result is a PKE struct that contains both KEM keys and their respetive lengths.
//
// NOTE: The calling application is responsible for deallocating PKE.ek and PKE.dk
static struct PKE KeyGen_internal(const struct PARAMS* params, const union byte* d, 
					const union byte* z) {
	struct PKE pke_keys;	// output of PKE_KeyGen()
	struct PKE kem_keys;	// output KEM encaps and decaps keys

	union byte* H_ek;	// Hash of PKE encryption key
	unsigned int it;	// iterator used to append to KEM decaps key

	// Run the key generation algorithm for K-PKE
	pke_keys = PKE_KeyGen(params, d);

	// Set KEM encaps key equal to PKE encryption key
	kem_keys.ek = pke_keys.ek;
	kem_keys.ek_len = pke_keys.ek_len;

	// Build KEM decaps key
	kem_keys.dk_len = 768 * params->k.e + 96;
	kem_keys.dk = malloc(sizeof(union byte) * kem_keys.dk_len);
	
	// Append PKE decryption key to KEM decaps key
	for (int i=0; i < pke_keys.dk_len; i++) {
		kem_keys.dk[i] = pke_keys.dk[i];
	}

	// Append PKE encryption key to KEM decaps key
	it = pke_keys.dk_len;
	for (int i=0; i < kem_keys.ek_len; i++) {
		kem_keys.dk[i+it] = kem_keys.ek[i];
	}

	// Compute SHA3-256 hash of PKE encryption key
	H_ek = H(kem_keys.ek, kem_keys.ek_len);

	// Append hash to KEM decaps key
	it += kem_keys.ek_len;
	for (int i=0; i < 32; i++) {
		kem_keys.dk[i+it] = H_ek[i];
	}

	// Append randomness z to KEM decaps key
	it += 32;
	for (int i=0; i < 32; i++) {
		kem_keys.dk[i+it] = z[i];
	}

	// Free dynamically allocated memory
	free(pke_keys.dk);
	free(H_ek);

	return kem_keys;
}

// Use a KEM encapsulation key generated by KeyGen_internal() and randomness, m of
// length 32, to generate a shared secret key and an associated ciphertext. The 
// result is a KEM struct containing the generated values.
//
// NOTE: The calling application is responsible for deallocating the ciphertext 
// 	KEM.c; however, the shared secret key KEM.K does not need to be deallocated
// 	since it is always a fixed length of 32 bytes.
static struct KEM Encaps_internal(const struct PARAMS* params, const union byte* ek,
					const union byte* m) {
	struct KEM Kc;		// output shared secret key and ciphertext
	union byte r[32];	// randomness output from G function

	// Temporary arrays
	union byte* input;
	union byte* output;

	unsigned int ek_len = 384 * params->k.e + 32;

	// Set ciphertext length according to input parameter set
	Kc.c_len = 32 * (params->du.e * params->k.e + params->dv.e);

	// Compute hash of input encryption key
	output = H(ek, ek_len);
	
	// Concatenate m and H(ek)
	input = malloc(sizeof(union byte) * (32 + 32));
	for (int i=0; i < 32; i++) input[i] = m[i];
	for (int i=0; i < 32; i++) input[i+32] = output[i];
	free(output);

	// Derive shared secret key and randomness from hash
	// of m||H(ek)
	output = G(input, 64);
	free(input);

	// Split output of G function into Kc.K and r
	for (int i=0; i < 32; i++) Kc.K[i] = output[i];
	for (int i=0; i < 32; i++) r[i] = output[i+32];
	free(output);

	// Encrypt random message m with randomness r
	Kc.c = PKE_Encrypt(params, ek, m, r);

	return Kc;
}

// Use a KEM decapsulation key generated by KeyGen_internal() to produce a shared
// secret key from a ciphertext, c of length 32*(du*k+dv). 
//
// NOTE: The calling application is responsible for deallocating the return array
static union byte* Decaps_internal(const struct PARAMS* params, const union byte* dk,
					const union byte* c) {
	union byte* K1;		// output shared secret key
	union byte* K2;		// rejected shared secret key	
	struct PKE keys;	// PKE encryption/decryption keys extracted from dk

	// Temporary arrays
	union byte* input;
	union byte* output;

	union byte h[32];		// hash of PKE encryption key
	union byte z[32];		// implicit rejection value
	union byte* m;			// decrypted ciphertext
	union byte r[32];		// derived randomness from hash of m||h
	union byte* c2;			// re-encrypted ciphertext

	unsigned int c_len = 32 * (params->du.e * params->k.e + params->dv.e);

	// Extract PKE decryption key from KEM decaps key
	keys.dk_len = 384 * params->k.e;
	keys.dk = malloc(sizeof(union byte) * keys.dk_len);
	for (int i=0; i < keys.dk_len; i++) {
		keys.dk[i] = dk[i];
	}

	// Extract PKE encryption key from KEM decaps key
	keys.ek_len = (768 * params->k.e + 32) - keys.dk_len;
	keys.ek = malloc(sizeof(union byte) * keys.ek_len);
	for (int i=0; i < keys.ek_len; i++) {
		keys.ek[i] = dk[i + keys.dk_len];
	}

	// Extract hash of PKE encryption key from KEM decaps key
	for (int i=0; i < 32; i++) {
		h[i] = dk[i + keys.ek_len + keys.dk_len];
	}

	// Extract implicit rejection value from KEM decaps key
	for (int i=0; i < 32; i++) {
		z[i] = dk[i + keys.ek_len + keys.dk_len + 32];
	}

	// Decrypt ciphertext
	m = PKE_Decrypt(params, dk, c);

	// Concatenate m and h
	input = malloc(sizeof(union byte) * 64);
	for (int i=0; i < 32; i++) input[i] = m[i];
	for (int i=0; i < 32; i++) input[i+32] = h[i];

	// Compute hash of m || h
	output = G(input, 64);
	free(input);

	// Split output of G function into two 32-byte arrays
	K1 = malloc(sizeof(union byte) * 32);
	for (int i=0; i < 32; i++) K1[i] = output[i];
	for (int i=0; i < 32; i++) r[i] = output[i+32];
	free(output);

	// Concatenate z and c
	input = malloc(sizeof(union byte) * (32 + c_len));
	for (int i=0; i < 32; i++) input[i] = z[i];
	for (int i=0; i < c_len; i++) input[i+32] = c[i];
	
	// Compute hash of z || c
	K2 = J(input, 32+c_len);
	free(input);

	// Re-encrypt using the derived randomness r
	c2 = PKE_Encrypt(params, keys.ek, m, r);

	// If ciphertexts do not match, then implicitly reject
	for (int i=0; i < c_len; i++) {
		if (c[i].e != c2[i].e) {
			free(K1);
			K1 = K2; 
			break;
		}
	}

	// Free up dynamically allocated memory
	if (K1 != K2) free(K2);
	free(keys.ek);
	free(keys.dk);
	free(m);
	free(c2);

	return K1;
}

/********************************
 * ML-KEM External Algorithms
 * *****************************/

// Derives a KEM encapsulation and decapsulation key pair from
// random bytes using KeyGen_internal() 
struct PKE KEM_KeyGen(const struct PARAMS* params) {
	union byte* d;
	union byte* z;
	struct PKE result;

	d = getRandomBytes();
	z = getRandomBytes();

	if ((d == NULL) || (z == NULL)) {
		if (d != NULL) free(d);
		if (z != NULL) free(z);

		printf("ml_kem.c:KEM_KeyGen() :: Random bit generation failed\n");
		exit(EXIT_FAILURE);
	}

	result = KeyGen_internal(params, d, z);
	free(d);
	free(z);

	return result;
}

// Derives a shared secret key and associated ciphertext using the provided
// KEM encapsulation key. Primary purpose is to provide input checking for
// Encaps_internal().
struct KEM KEM_Encaps(const struct PARAMS* params, const union byte* ek, 
			unsigned int ek_len) {
	union byte* m;
	union byte* test;
	union integer* temp;
	struct KEM result;
	unsigned int len = 384 * params->k.e;
	
	// Type check
	if ((len + 32) != ek_len) {
		printf("ml_kem.c:KEM_Encaps() :: Type check failed\n");
		exit(EXIT_FAILURE);
	}

	// Modulus check
	test = malloc(sizeof(union byte) * len);
	for (int i=0; i < len; i++) test[i] = ek[i];


	temp = ByteDecode(test, 12);
	free(test);

	test = ByteEncode(temp, 12);
	free(temp);

	for (int i=0; i < len; i++) {
		if (test[i].e != ek[i].e) {
			printf("ml_kem.c:KEM_Encaps() :: Modulus check failed\n");
			free(test);
			exit(EXIT_FAILURE);
		}
	}
	free(test);

	// Generate an array of 32 random bytes
	m = getRandomBytes();
	if (m == NULL) {
		printf("ml_kem.c:KEM_Encaps() :: Random bit generation failed\n");
		exit(EXIT_FAILURE);
	}

	result = Encaps_internal(params, ek, m);
	free(m);

	return result;
}

// Derives a shared secret key from the proved ciphertext using the provided 
// KEM decapsulation key. Primary purpose is to provide input checking for 
// Decaps_internal()
union byte* KEM_Decaps(const struct PARAMS* params, const union byte* dk, unsigned int dk_len,
			const union byte* c, unsigned int c_len) {
	union byte* test;
	union byte* input;
	unsigned int len;
	unsigned int offset = 384 * params->k.e;

	// Ciphertext type check
	len = 32 * (params->du.e * params->k.e + params->dv.e);
	if (c_len != len) {
		printf("ml_kem.c:KEM_Decaps() :: Ciphertext type check failed\n");
		exit(EXIT_FAILURE);
	}

	// Decapsulation key type check
	len = 768 * params->k.e + 96;
	if (dk_len != len) {
		printf("ml_kem.c:KEM_Decaps() :: Decapsulation key type check failed\n");
		exit(EXIT_FAILURE);
	}

	// Hash check
	len = (768 * params->k.e + 32) - offset;
	input = malloc(sizeof(union byte) * len);

	for (int i=0; i < len; i++) input[i] = dk[i + offset];
	test = H(input, len);
	free(input);
	
	offset = len + offset;
	for (int i=0; i < 32; i++) {
		if (test[i].e != dk[i + offset].e) {
			printf("ml_kem.c:KEM_Decaps() :: Hash check failed\n");
			free(test);
			exit(EXIT_FAILURE);
		}
	}
	free(test);
	

	return Decaps_internal(params, dk, c);
}

// Initialize a PARAMS struct based on the chosen parameter set - see
// FIPS-203:8
const struct PARAMS init(enum ML_KEM param_set) {
	
	struct PARAMS params;

	switch (param_set) {
	case 512:
		params.k.e = 2;
		params.n1.e = 3;
		params.n2.e = 2;
		params.du.e = 10;
		params.dv.e = 4;
		break;
	case 768:
		params.k.e = 3;
		params.n1.e = 2;
		params.n2.e = 2;
		params.du.e = 10;
		params.dv.e = 4;
		break;
	case 1024:
		params.k.e = 4;
		params.n1.e = 2;
		params.n2.e = 2;
		params.du.e = 11;
		params.dv.e = 5;
		break;
	default:
		printf("ml_kem.c:init() :: Invalid paramater set provided\n");
		exit(EXIT_FAILURE);
	}

	return params;
}
