CC = gcc
CFLAGS = -Wall

test01: ml-kem.o BitRev7_test01.c
	$(CC) $(CFLAGS) $^ -o $@

ml-kem.o: ml-kem.h ml-kem.c
	$(CC) $(CFLAGS) -c $^

clean:
	rm *.o *.gch
