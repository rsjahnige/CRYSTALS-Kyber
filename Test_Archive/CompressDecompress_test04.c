#include "ml_kem.h"
#include <stdio.h>

int main() {
	union integer y, x;

	for (int i=0; i < 3329; i++) {
		y.t = i;

		if (y.t < 2) {
			x = Compress(Decompress(y, 1), 1);
			if (x.t != y.t) printf("ERROR-1:: y.t=%d - x.t=%d\n", y.t, x.t);
		} 
		
		if (y.t < 4) {
			x = Compress(Decompress(y, 2), 2);
			if (x.t != y.t) printf("ERROR-2:: y.t=%d - x.t=%d\n", y.t, x.t);
		} 
		
		if (y.t < 8) {
			x = Compress(Decompress(y, 3), 3);
			if (x.t != y.t) printf("ERROR-3:: y.t=%d - x.t=%d\n", y.t, x.t);
		} 
		
		if (y.t < 16) {			
			x = Compress(Decompress(y, 4), 4);
			if (x.t != y.t) printf("ERROR-4:: y.t=%d - x.t=%d\n", y.t, x.t);
		} 
		
		if (y.t < 32) {
			x = Compress(Decompress(y, 5), 5);
			if (x.t != y.t) printf("ERROR-5:: y.t=%d - x.t=%d\n", y.t, x.t);
		} 
		
		if (y.t < 64) {
			x = Compress(Decompress(y, 6), 6);
			if (x.t != y.t) printf("ERROR-6:: y.t=%d - x.t=%d\n", y.t, x.t);
		} 
		
		if (y.t < 128) {
			x = Compress(Decompress(y, 7), 7);
			if (x.t != y.t) printf("ERROR-7:: y.t=%d - x.t=%d\n", y.t, x.t);
		} 

		if (y.t < 256) {
			x = Compress(Decompress(y, 8), 8);
			if (x.t != y.t) printf("ERROR-8:: y.t=%d - x.t=%d\n", y.t, x.t);
		} 

		if (y.t < 512) {
			x = Compress(Decompress(y, 9), 9);
			if (x.t != y.t) printf("ERROR-9:: y.t=%d - x.t=%d\n", y.t, x.t);
		} 

		if (y.t < 1024) {
			x = Compress(Decompress(y, 10), 10);
			if (x.t != y.t) printf("ERROR-10:: y.t=%d - x.t=%d\n", y.t, x.t);
		} 

		if (y.t < 2048) {
			x = Compress(Decompress(y, 11), 11);
			if (x.t != y.t) printf("ERROR-11:: y.t=%d - x.t=%d\n", y.t, x.t);
		} 

		if (y.t < 3329) {
			x = Compress(Decompress(y, 12), 12);
			if (x.t != y.t) printf("ERROR-12:: y.t=%d - x.t=%d\n", y.t, x.t);
		} 
	}

	printf("Test Complete!\n");

	return 0;
}
