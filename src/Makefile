
CC = gcc

all: test avx2 neon

test: test.c
	$(CC) $(CLFAGS) $< -o $@

avx2: avx2.c
	$(CC) $(CLFAGS) $< -o $@

neon: neon.c
	$(CC) $(CLFAGS) $< -o $@

.PHOYNY: clean

clean:
	rm -f test
	rm -f avx2
	rm -f neon
