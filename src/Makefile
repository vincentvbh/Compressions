
CC = gcc

all: test avx2_ref neon_ref

test: test.c
	$(CC) $(CLFAGS) $< -o $@

avx2_ref: avx2_ref.c
	$(CC) $(CLFAGS) $< -o $@

neon_ref: neon_ref.c
	$(CC) $(CLFAGS) $< -o $@

neon: neon.c
	$(CC) $(CLFAGS) $< -o $@

.PHOYNY: clean

clean:
	rm -f test
	rm -f avx2_ref
	rm -f neon_ref
	rm -f neon
