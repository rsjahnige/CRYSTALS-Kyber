#include "ml_kem.h"
#include <stdio.h>

int main() {
	
	struct PARAMS params;
	struct PKE keys;
	union byte randomness[32];
	union byte message[32];
	union byte* ciphertext;
	union byte* dmessage;
	unsigned int len;

	params = init(ML_KEM_512);

	// Initialize randomness array
	for (int i=0; i < 32; i++) {
		randomness[i].e = i;
		message[i].e = i % 5;
	}

	printf("Plaintext: ");
	for (int i=0; i < 32; i++) {
		printf("%d ", message[i].e);
	}
	printf("\n\n");

	// Generate public-private key pair
	keys = PKE_KeyGen(&params, randomness);

	ciphertext = PKE_Encrypt(&params, keys.ek, message, randomness);	

	len = 32 * (params.du.e * params.k.e + params.dv.e);
	
	printf("Ciphertext: ");
	for (int i=0; i < len; i++) {
		printf("%d ", ciphertext[i].e);
	}	
	printf("\n\n");

	dmessage = PKE_Decrypt(&params, keys.dk, ciphertext);

	printf("Plaintext: ");
	for (int i=0; i < 32; i++) {
		printf("%d ", dmessage[i].e);
	}
	printf("\n\n");

	// Free up memory
	free(keys.ek);
	free(keys.dk);

	free(ciphertext);
	free(dmessage);

	return 0;
}
