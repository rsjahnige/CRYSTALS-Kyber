#include "ml_kem.h"
#include <stdio.h>

int main() {

	struct PARAMS params;
	struct PKE keys;

	params = init(ML_KEM_512);
	keys = KEM_KeyGen(&params);

	printf("Encapsulation key: ");
	for (int i=0; i < keys.ek_len; i++) {
		printf("%d ", keys.ek[i].e);
	}
	printf("\n\n");

	printf("Decapsulation key: ");
	for (int i=0; i < keys.dk_len; i++) {
		printf("%d ", keys.dk[i].e);
	}
	printf("\n\n");

	free(keys.ek);
	free(keys.dk);

	return 0;
}
