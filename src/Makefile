
CC = gcc

CFLAGS += -O3

AARCH64_FLAGS = -march=armv8-a+dotprod

all: test avx2_ref neon_ref

test: test.c
	$(CC) $(CLFAGS) $< -o $@

avx2_ref: avx2_ref.c
	$(CC) $(CLFAGS) $< -o $@

neon_ref: neon_ref.c
	$(CC) $(CLFAGS) $< -o $@

neon: neon.c
	$(CC) $(CLFAGS) $(AARCH64_FLAGS) $< -o $@

.PHOYNY: clean

clean:
	rm -f test
	rm -f avx2_ref
	rm -f neon_ref
	rm -f neon
