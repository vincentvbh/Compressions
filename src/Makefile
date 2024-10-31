
# platform independent flags
CC                     ?= gcc
CFLAGS                 += -Wall -O3

# cross-compilation targeting Armv7-A on linux
AARCH32_CROSS_CC        = arm-linux-gnueabi-gcc-12
AARCH32_CROSS_FLAGS     = -march=armv7-a -z noexecstack -g -static
AARCH32_ASM             = __asm_compress_armv7e-m.S

# cross-compilation targeting Armv7-M on macOS
ARMv7-M_ARCH_FLAGS      = -mthumb -mcpu=cortex-m4 -mfloat-abi=hard -mfpu=fpv4-sp-d16
ARMv7-M_DEFINES         = -DSTM32F4 -DSTM32F407VG

ARMv7-M_CROSS_PREFIX   ?= arm-none-eabi
ARMv7-M_CROSS_CC        = $(ARMv7-M_CROSS_PREFIX)-gcc
ARMv7-M_CROSS_LD        = $(ARMv7-M_CROSS_PREFIX)-gcc
ARMv7-M_CROSS_OBJCOPY   = $(ARMv7-M_CROSS_PREFIX)-objcopy
ARMv7-M_CROSS_OBJDUMP   = $(ARMv7-M_CROSS_PREFIX)-objdump
ARMv7-M_CROSS_GDB       = $(ARMv7-M_CROSS_PREFIX)-gdb
OPENCM3_DIR             = libopencm3
OPENCM3_OBJS            = hal-opencm3.o
OPENCM3_LIBNAME         = opencm3_stm32f4

ARMv7-M_LDSCRIPT        = stm32f4discovery.ld


ARMv7-M_CROSS_CFLAGS   += -O3 -Wall \
						  -I$(OPENCM3_DIR)/include \
						  -fno-common $(ARMv7-M_ARCH_FLAGS) -MD $(ARMv7-M_DEFINES)

ARMv7-M_CROSS_LDFLAGS  += --static -Wl,--start-group -lc -lgcc -lnosys -Wl,--end-group \
		 				  -T$(ARMv7-M_LDSCRIPT) -nostartfiles \
		 				  $(ARMv7-M_ARCH_FLAGS) \
		 				  -L$(OPENCM3_DIR)/lib

opencm3_lib:
	@if [ ! "`ls -A $(OPENCM3_DIR)`" ] ; then \
		printf "######## ERROR ########\n"; \
		printf "\tlibopencm3 is not initialized.\n"; \
		printf "\tPlease run (in the root directory):\n"; \
		printf "\t$$ git submodule init\n"; \
		printf "\t$$ git submodule update\n"; \
		printf "\tbefore running make.\n"; \
		printf "######## ERROR ########\n"; \
		exit 1; \
		fi
	make -C $(OPENCM3_DIR)

%.bin: %.elf
	$(ARMv7-M_CROSS_OBJCOPY) -Obinary $(*).elf $(*).bin

%.elf: %.o $(OPENCM3_OBJS) $(ARMv7-M_CROSS_LDSCRIPT)
	$(ARMv7-M_CROSS_LD) -o $(*).elf $(*).o $(OPENCM3_OBJS) $(ARMv7-M_CROSS_LDFLAGS) -l$(OPENCM3_LIBNAME)

%.o: %.c opencm3_lib
	$(ARMv7-M_CROSS_CC) $(ARMv7-M_CROSS_CFLAGS) -o $@ -c $<

%.o: %.S opencm3_lib
	$(ARMv7-M_CROSS_CC) $(ARMv7-M_CROSS_CFLAGS) -o $@ -c $<

ARMv7-M_ASM             = __asm_compress_armv7e-m.S

# aarch64 flags
AARCH64_FLAGS           = -march=armv8-a

# avx2 flags
AVX2_FLAGS              = -mavx2 -z noexecstack

# Jasmin compiler
JASMINC                 = jasminc

# Jasmin objects.
JAZZ_X86_OBJ            = jazz_x86.s
JAZZ_AVX2_OBJ           = jazz_avx2.s
JAZZ_ARMv7E-M_OBJ       = jazz_armv7e-m.s

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

jazz_armv7e-m_native.elf: jazz_armv7e-m_native.o jazz_armv7e-m.c jazz_armv7e-m.jazz jazz_armv7e-m.s $(OPENCM3_OBJS) $(ARMv7-M_CROSS_LDSCRIPT)
	$(ARMv7-M_CROSS_CC) $(OPENCM3_OBJS) jazz_armv7e-m.s $< -o $@ $(ARMv7-M_CROSS_LDFLAGS) -l$(OPENCM3_LIBNAME)

neon_ref: neon_ref.c
	$(CC) $(CFLAGS) $< -o $@

neon: neon.c
	$(CC) $(CFLAGS) $(AARCH64_FLAGS) $< -o $@

armv7e-m_ref: armv7e-m_ref.c
	$(CC) $(CFLAGS) $< -o $@

armv7e-m: armv7e-m.c $(AARCH32_ASM)
	$(AARCH32_CROSS_CC) $(CFLAGS) $(AARCH32_CROSS_FLAGS) $(AARCH32_ASM) $< -o $@

armv7e-m_native.elf: armv7e-m_native.o $(ARMv7-M_ASM) $(OPENCM3_OBJS) $(ARMv7-M_CROSS_LDSCRIPT)
	$(ARMv7-M_CROSS_CC) $(OPENCM3_OBJS) $(ARMv7-M_ASM) $< -o $@ $(ARMv7-M_CROSS_LDFLAGS) -l$(OPENCM3_LIBNAME)

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
	rm -f armv7e-m_native
	rm -f jazz_armv7e-m
	rm -f $(JAZZ_ARMv7E-M_OBJ)
	rm -f *.o *.d *.elf *.bin

