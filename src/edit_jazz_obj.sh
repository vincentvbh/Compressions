
if [[ "$OSTYPE" = "darwin"* ]]
then
    SED=gsed
else
    SED=sed
fi

$SED -i '1i .type ___jazz_poly_compress1, %function' jazz_armv7e-m.s
$SED -i '1i .type __jazz_poly_compress1, %function' jazz_armv7e-m.s
$SED -i '1i .type ___jazz_poly_compress4, %function' jazz_armv7e-m.s
$SED -i '1i .type __jazz_poly_compress4, %function' jazz_armv7e-m.s
$SED -i '1i .type ___jazz_poly_compress10, %function' jazz_armv7e-m.s
$SED -i '1i .type __jazz_poly_compress10, %function' jazz_armv7e-m.s

$SED -i '1i .type ___jazz_stack_poly_compress1, %function' jazz_armv7e-m.s
$SED -i '1i .type __jazz_stack_poly_compress1, %function' jazz_armv7e-m.s
$SED -i '1i .type ___jazz_stack_poly_compress4, %function' jazz_armv7e-m.s
$SED -i '1i .type __jazz_stack_poly_compress4, %function' jazz_armv7e-m.s
$SED -i '1i .type ___jazz_stack_poly_compress10, %function' jazz_armv7e-m.s
$SED -i '1i .type __jazz_stack_poly_compress10, %function' jazz_armv7e-m.s

