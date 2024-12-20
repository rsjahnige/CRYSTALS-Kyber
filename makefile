CC = gcc
CFLAGS = -Wall

test01: ml_kem.o BitRev7_test01.c
	$(CC) $(CFLAGS) $^ -o $@

ml_kem.o: ml_kem.h ml_kem.c
	$(CC) $(CFLAGS) -c $^

clean:
	rm *.o *.gch
