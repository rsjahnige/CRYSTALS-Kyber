CC = gcc
CFLAGS = -Wall -g

test10: sha3.o ml_kem.o PKE_EncryptDecrypt_test.c
	$(CC) $(CFLAGS) $^ -o $@

test09: sha3.o ml_kem.o KeyGen_test.c
	$(CC) $(CFLAGS) $^ -o $@

gammaTest: sha3.o ml_kem.o GammaLogic_test.c
	$(CC) $(CFLAGS) $^ -o $@

test08: sha3.o ml_kem.o NTT_test08.c
	$(CC) $(CFLAGS) $^ -o $@

zetaTest: sha3.o ml_kem.o ZetaLogic_test.c
	$(CC) $(CFLAGS) $^ -o $@

test07: sha3.o ml_kem.o SampleCBD_test07.c
	$(CC) $(CFLAGS) $^ -o $@

test06: sha3.o ml_kem.o SampleNTT_test06.c
	$(CC) $(CFLAGS) $^ -o $@

test05: sha3.o SHA_Testing.c
	$(CC) $(CFLAGS) $^ -o $@

test04: ml_kem.o CompressDecompress_test04.c
	$(CC) $(CFLAGS) $^ -o $@

test03: ml_kem.o EncodeDecode_test03.c
	$(CC) $(CFLAGS) $^ -o $@

test02: ml_kem.o BitsAndBytes_test02.c
	$(CC) $(CFLAGS) $^ -o $@

test01: ml_kem.o BitRev7_test01.c
	$(CC) $(CFLAGS) $^ -o $@

ml_kem.o: sha3.o ml_kem.h ml_kem.c
	$(CC) $(CFLAGS) -c $^

sha3.o: sha3.h sha3.c
	$(CC) $(CFLAGS) -c $^

clean:
	rm *.o *.gch
