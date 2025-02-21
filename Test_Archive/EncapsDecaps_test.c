#include "ml_kem.h"
#include <stdio.h>

int main() {

	struct PARAMS params;
	struct PKE keys;
	struct KEM kem;
	union byte* sym_key;

	params = init(ML_KEM_512);
	keys = KEM_KeyGen(&params);

	kem = KEM_Encaps(&params, keys.ek, keys.ek_len);
	sym_key = KEM_Decaps(&params, keys.dk, keys.dk_len, kem.c, kem.c_len);

	for (int i=0; i < 32; i++) {
		if (kem.K[i].e != sym_key[i].e) {
			printf("Test failed\n");
			goto free_mem;
		}
	}

	printf("Test successful\n");

	free_mem:
	free(keys.ek);
	free(keys.dk);
	free(kem.c);
	free(sym_key);

	return 0;
}
