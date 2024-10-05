
CC = gcc

CFLAGS += -Wall -O3 -mtune=native -march=native

# AARCH64_FLAGS = -march=armv8-a+dotprod
AARCH64_FLAGS = -march=armv8-a
AVX2_FLAGS = -mavx2 -z noexecstack

JASMINC = jasminc

all: test avx2_ref neon_ref armv7e-m_ref

jazz_x86.s: jazz_x86.jazz
	$(JASMINC) $< -o $@

test: test.c
	$(CC) $(CLFAGS) $< -o $@

avx2_ref: avx2_ref.c
	$(CC) $(CLFAGS) $< -o $@

avx2: avx2.c
	$(CC) $(CLFAGS) $(AVX2_FLAGS) $< -o $@

jazz_x86: jazz_x86.c jazz_x86.jazz jazz_x86.s
	$(CC) $(CLFAGS) jazz_x86.s $< -o $@

neon_ref: neon_ref.c
	$(CC) $(CLFAGS) $< -o $@

neon: neon.c
	$(CC) $(CLFAGS) $(AARCH64_FLAGS) $< -o $@

armv7e-m_ref: armv7e-m_ref.c
	$(CC) $(CLFAGS) $< -o $@

.PHOYNY: clean

clean:
	rm -f test
	rm -f avx2_ref
	rm -f neon_ref
	rm -f neon
	rm -f armv7e-m_ref
	rm -f jazz_x86
