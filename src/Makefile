
CC = gcc

CFLAGS += -Wall -O3 -mtune=native

# AARCH64_FLAGS = -march=armv8-a+dotprod
AARCH64_FLAGS = -march=armv8-a
AVX2_FLAGS = -mavx2 -z noexecstack

JASMINC = jasminc

JAZZ_OBJ = jazz_x86.s jazz_avx2.s

all: test avx2_ref neon_ref armv7e-m_ref

jazz_x86.s: jazz_x86.jazz
	$(JASMINC) $< -o $@

jazz_avx2.s: jazz_avx2.jazz
	$(JASMINC) $< -o $@

test: test.c
	$(CC) $(CFLAGS) $< -o $@

avx2_ref: avx2_ref.c
	$(CC) $(CFLAGS) $< -o $@

avx2: avx2.c
	$(CC) $(CFLAGS) $(AVX2_FLAGS) $< -o $@

jazz_x86: jazz_x86.c jazz_x86.jazz jazz_x86.s
	$(CC) $(CFLAGS) jazz_x86.s $< -o $@

jazz_avx2: jazz_avx2.c jazz_avx2.jazz jazz_avx2.s
	$(CC) $(AVX2_FAGS) jazz_avx2.s $< -o $@

neon_ref: neon_ref.c
	$(CC) $(CFLAGS) $< -o $@

neon: neon.c
	$(CC) $(CFLAGS) $(AARCH64_FLAGS) $< -o $@

armv7e-m_ref: armv7e-m_ref.c
	$(CC) $(CFLAGS) $< -o $@

.PHOYNY: clean

clean:
	rm -f test
	rm -f avx2_ref
	rm -f avx2
	rm -f jazz_x86
	rm -f jazz_avx2
	rm -f $(JAZZ_OBJ)
	rm -f neon_ref
	rm -f neon
	rm -f armv7e-m_ref
