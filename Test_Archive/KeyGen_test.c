#include "ml_kem.h"
#include "stdio.h"

int main() {
	union byte input[32];
	struct ML_KEM params = {2, 3};	// ML-KEM-512
	struct PKE keys;

	for (int i=0; i < 32; i++) input[i].e = i;

	keys = KeyGen(&params, input);
 
	printf("Encryption key: ");
	for (int i=0; i < params.k.e; i++) {
		for (int j=0; j < 256; j++) {
			printf("%d ", keys.ek[i][j].e);
		}
	}
	printf("\n\n");

	printf("Decryption key: ");
	for (int i=0; i < params.k.e; i++) {
		for (int j=0; j < 256; j++) {
			printf("%d ", keys.dk[i][j].e);
		}
	}
	printf("\n\n");

	printf("Rho: ");
	for (int i=0; i < 32; i++) {
		printf("%d ", keys.rho[i].e);
	}
	printf("\n");

	// Free up memory
	for (int i=0; i < params.k.e; i++) {
		free(keys.ek[i]);
		free(keys.dk[i]);
	}
	free(keys.ek);
	free(keys.dk);
	free(keys.rho);
}
