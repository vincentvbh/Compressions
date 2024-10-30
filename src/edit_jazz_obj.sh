
sed -i 's/	.global	___jazz_poly_compress1/	.global	___jazz_poly_compress1\n\t.type ___jazz_poly_compress1, %function/' jazz_armv7e-m.s
sed -i 's/	.global	__jazz_poly_compress1/	.global	__jazz_poly_compress1\n\t.type __jazz_poly_compress1, %function/' jazz_armv7e-m.s
