
CC = gcc

CFLAGS += -Wall -O3

# AARCH64_FLAGS = -march=armv8-a+dotprod
AARCH64_FLAGS = -march=armv8-a
AVX2_FLAGS =

all: test avx2_ref neon_ref armv7e-m_ref

test: test.c
	$(CC) $(CLFAGS) $< -o $@

avx2_ref: avx2_ref.c
	$(CC) $(CLFAGS) $< -o $@

avx2: avx2.c
	$(CC) $(CLFAGS) $(AVX2_FLAGS) $< -o $@

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
