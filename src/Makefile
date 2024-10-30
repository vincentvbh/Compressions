
CC ?= gcc

CFLAGS += -Wall -O3

AARCH32_CROSS_CC = arm-linux-gnueabi-gcc-12
AARCH32_CROSS_FLAGS = -march=armv7-a -z noexecstack -g -static
AARCH32_ASM = __asm_compress_armv7e-m.S

AARCH64_FLAGS = -march=armv8-a
AVX2_FLAGS = -mavx2 -z noexecstack

JASMINC = jasminc

JAZZ_X86_OBJ = jazz_x86.s

JAZZ_AVX2_OBJ = jazz_avx2.s

JAZZ_ARMv7E-M_OBJ = jazz_armv7e-m.s

all: test avx2_ref neon_ref armv7e-m_ref

jazz_x86.s: jazz_x86.jazz
	$(JASMINC) $< -o $@

jazz_avx2.s: jazz_avx2.jazz
	$(JASMINC) $< -o $@

jazz_armv7e-m.s: jazz_armv7e-m.jazz
	$(JASMINC) -arch=arm-m4 $< -o $@
	sh ./edit_jazz_obj.sh

test: test.c
	$(CC) $(CFLAGS) $< -o $@

avx2_ref: avx2_ref.c
	$(CC) $(CFLAGS) $< -o $@

avx2: avx2.c
	$(CC) $(CFLAGS) $(AVX2_FLAGS) $< -o $@

jazz_x86: jazz_x86.c jazz_x86.jazz jazz_x86.s
	$(CC) $(CFLAGS) jazz_x86.s $< -o $@

jazz_avx2: jazz_avx2.c jazz_avx2.jazz jazz_avx2.s
	$(CC) $(CFLAGS) $(AVX2_FAGS) jazz_avx2.s $< -o $@

jazz_armv7e-m: jazz_armv7e-m.c jazz_armv7e-m.jazz jazz_armv7e-m.s
	$(AARCH32_CROSS_CC) $(CFLAGS) $(AARCH32_CROSS_FLAGS) jazz_armv7e-m.s $< -o $@

neon_ref: neon_ref.c
	$(CC) $(CFLAGS) $< -o $@

neon: neon.c
	$(CC) $(CFLAGS) $(AARCH64_FLAGS) $< -o $@

armv7e-m_ref: armv7e-m_ref.c
	$(CC) $(CFLAGS) $< -o $@

armv7e-m: armv7e-m.c $(AARCH32_ASM)
	$(AARCH32_CROSS_CC) $(CFLAGS) $(AARCH32_CROSS_FLAGS) $(AARCH32_ASM) $< -o $@

.PHOYNY: clean

clean:
	rm -f test
	rm -f avx2_ref
	rm -f avx2
	rm -f jazz_x86
	rm -f jazz_avx2
	rm -f $(JAZZ_X86_OBJ)
	rm -f $(JAZZ_AVX2_OBJ)
	rm -f neon_ref
	rm -f neon
	rm -f armv7e-m_ref
	rm -f armv7e-m
	rm -f jazz_armv7e-m
	rm -f $(JAZZ_ARMv7E-M_OBJ)

